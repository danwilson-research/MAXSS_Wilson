# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 17:30:54 2026

@author: dw557
"""

import xarray as xr
import numpy as np
import os
import pandas as pd
from os import path
from glob import glob
from pyproj import Geod
from tqdm import tqdm

def calculate_ellipsoidal_areas(ds, res=0.25):
    geod = Geod(ellps="WGS84")
    storm_lat = ds.latitude.values
    lon_sample = ds.longitude.values[0]
    row_areas = []

    for lat in storm_lat:
        lats = [lat-res/2, lat-res/2, lat+res/2, lat+res/2, lat-res/2]
        lons = [lon_sample-res/2, lon_sample+res/2, lon_sample+res/2, lon_sample-res/2, lon_sample-res/2]
        poly_area, _ = geod.polygon_area_perimeter(lons, lats)
        row_areas.append(abs(poly_area))

    return xr.DataArray(row_areas, coords={'latitude': storm_lat}, dims='latitude')

def calc_hourly_flux(MAXSS_working_directory, output_base, netcdf_output_root, runs, MAXSS_regions, specified_storms, salinity_threshold=10.0):
    master_summary = {}
    
    for region in MAXSS_regions:
        reference_directory = path.join(MAXSS_working_directory, "output/MAXSS_RUN/maxss/storm-atlas/tropical/ibtracs", region)
        year_directory_list = glob(path.join(reference_directory, "*/"))

        for year_dir in year_directory_list:
            year_name = os.path.basename(os.path.normpath(year_dir))
            storm_directory_list = glob(path.join(year_dir, "*/"))

            for storm_path in storm_directory_list:
                storm_name = os.path.basename(os.path.normpath(storm_path))

                # --- Fragment skipping print statement ---
                if len(specified_storms) > 0:
                    if not any(fragment in storm_path for fragment in specified_storms):
                        print(f"Skipping storm (no fragment match): {storm_path}")
                        continue

                land_fraction_path = path.join(MAXSS_working_directory, "maxss/storm-atlas/tropical/ibtracs/", region, year_name, storm_name, 'Resampled_for_fluxengine_MAXSS_land_fraction.nc')

                # Load both land_proportion and sea_ice_proportion from the single mask file
                with xr.open_dataset(land_fraction_path) as lf_ds:
                    land_fraction = lf_ds.land_proportion.load()
                    if "sea_ice_proportion" in lf_ds:
                        sea_ice_proportion = lf_ds.sea_ice_proportion.load()
                    else:
                        print(f"   -> Warning: 'sea_ice_proportion' not found in {land_fraction_path}; assuming 0 ice.")
                        sea_ice_proportion = None
                
                land_fraction = land_fraction.where(land_fraction <= 1.0, 1.0)
                water_proportion = 1.0 - land_fraction

                grid_areas = None
                run_pbar = tqdm(runs, desc=f"Processing {storm_name}", leave=False)

                for run_name in run_pbar:
                    storm_run_path = path.join(MAXSS_working_directory, "output", run_name, "maxss/storm-atlas/tropical/ibtracs", region, year_name, storm_name)
                    fluxengine_files = sorted(glob(os.path.join(storm_run_path, "*.nc")))

                    # --- Missing run skip print statement ---
                    if not fluxengine_files:
                        print(f"   -> Run {run_name} not found. Skipping.")
                        continue

                    # Vectorized reading of ALL files for this storm run using xarray
                    # Note: If Dask isn't installed, swap xr.open_mfdataset for Solution 1's list comprehension
                    with xr.open_mfdataset(fluxengine_files, combine='by_coords', chunks=None, parallel=False) as ds:
                        if grid_areas is None:
                            grid_areas = calculate_ellipsoidal_areas(ds, res=0.25)
                                    
                        # --- VECTORIZED MATH BLOCK ---
                        current_water_map = water_proportion.reindex_like(ds, method="nearest")

                        # Align sea ice in (time, latitude, longitude) and scale linearly by (1 - f_ice)
                        if sea_ice_proportion is not None:
                            current_ice_map = sea_ice_proportion.reindex_like(ds, method="nearest").fillna(0.0).clip(min=0.0, max=1.0)
                            ice_free_proportion = 1.0 - current_ice_map
                        else:
                            ice_free_proportion = 1.0
                        
                        # Enforce FluxEngine's valid_min / valid_max (which xarray ignores by default)
                        valid_min = ds.OF.attrs["valid_min"]
                        valid_max = ds.OF.attrs["valid_max"]
                        clean_OF = ds.OF.where(np.isfinite(ds.OF) & (ds.OF >= valid_min) & (ds.OF <= valid_max))

                        # Exclude estuarine pixels if EITHER salinity input falls below salinity_threshold
                        if salinity_threshold is not None:
                            is_estuarine = xr.zeros_like(clean_OF, dtype=bool)
                            if "pco2_sss" in ds:
                                is_estuarine = is_estuarine | (ds.pco2_sss < salinity_threshold) | (~np.isfinite(ds.pco2_sss))
                            if "OKS1" in ds:
                                is_estuarine = is_estuarine | (ds.OKS1 < salinity_threshold) | (~np.isfinite(ds.OKS1))
                            clean_OF = clean_OF.where(~is_estuarine)
                          
                        mass_map = clean_OF * grid_areas * (1.0 / 24.0) * current_water_map * ice_free_proportion
                        hourly_fluxes_tg = mass_map.sum(dim=['latitude', 'longitude'], skipna=True) / 1e12
                        
                        hourly_fluxes_arr = hourly_fluxes_tg.values
                        time_vals = ds.time.values
                        storm_total_tg = float(hourly_fluxes_tg.sum(skipna=True))

                        # --- SAVE NETCDF INSIDE CONTEXT MANAGER ---
                        out_folder = path.join(netcdf_output_root, region, year_name)
                        os.makedirs(out_folder, exist_ok=True)
                        processedFilePath = os.path.join(out_folder, f"{storm_name}_{run_name}.nc")

                        #Add salinity threshold to the netcdf metadata
                        thresh_attr = f"{salinity_threshold} PSU" if salinity_threshold is not None else "None"
                        
                        ds_out = xr.Dataset(
                            data_vars={
                                "Hourly_flux": (["time"], hourly_fluxes_arr, {"units": "Tg C hr-1"}),
                                "Total_flux": ([], storm_total_tg, {"units": "Tg C", "salinity_threshold": thresh_attr})
                            },
                            coords={"time": time_vals}
                        )
                        ds_out['time'].encoding['units'] = "seconds since 1970-01-01 00:00:00"
                        ds_out['time'].encoding['dtype'] = 'i8'
                        ds_out.to_netcdf(processedFilePath)

                    # --- Run success summary print statement ---
                    print(f"      - {run_name}: {storm_total_tg:.6f} Tg")

                    # --- UPDATE MASTER SUMMARY ---
                    storm_id = (storm_name, year_name, region)
                    if storm_id not in master_summary:
                        master_summary[storm_id] = {}
                    master_summary[storm_id][run_name] = storm_total_tg

    # --- CSV Generation ---
    final_rows = []
    for (storm, year, region), run_data in master_summary.items():
        row = {'Storm': storm, 'Year': year, 'Region': region}
        row.update(run_data)
        final_rows.append(row)

    if final_rows:
        df_final = pd.DataFrame(final_rows)
        cols = ['Storm', 'Year', 'Region'] + [r for r in runs if r in df_final.columns]
        df_final = df_final[cols]
        
        # Extract unique years processed in this specific run to name the file safely
        unique_years = sorted(list(set(df_final['Year'].astype(str))))
        if len(unique_years) == 1:
            csv_name = f"storm_component_flux_summary_{unique_years[0]}.csv"
        else:
            csv_name = f"storm_component_flux_summary_{unique_years[0]}_to_{unique_years[-1]}.csv"

        output_path = path.join(output_base, csv_name)
        os.makedirs(output_base, exist_ok=True)
        df_final.to_csv(output_path, index=False)
        
        print(f"\n{'='*40}")
        print(f"SUCCESS: Summary saved to {output_path}")
        print(f"{'='*40}")
    else:
        print("\nNo data was processed. Check your directory paths.")
