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

def calc_hourly_flux(MAXSS_working_directory, output_base, netcdf_output_root, runs, MAXSS_regions, specified_storms):
    master_summary = {}
    
    for region in MAXSS_regions:
        reference_directory = path.join(MAXSS_working_directory, "output/MAXSS_RUN/maxss/storm-atlas/tropical/ibtracs", region)
        year_directory_list = glob(path.join(reference_directory, "*/"))

        for year_dir in year_directory_list:
            year_name = os.path.basename(os.path.normpath(year_dir))
            storm_directory_list = glob(path.join(year_dir, "*/"))

            for storm_path in storm_directory_list:
                storm_name = os.path.basename(os.path.normpath(storm_path))

                # --- RESTORED: Fragment skipping print statement ---
                if len(specified_storms) > 0:
                    if not any(fragment in storm_path for fragment in specified_storms):
                        print(f"Skipping storm (no fragment match): {storm_path}")
                        continue

                land_fraction_path = path.join(MAXSS_working_directory, "maxss/storm-atlas/tropical/ibtracs/", region, year_name, storm_name, 'Resampled_for_fluxengine_MAXSS_land_fraction.nc')

                with xr.open_dataset(land_fraction_path) as lf_ds:
                    land_fraction = lf_ds.land_proportion.load()
                
                land_fraction = land_fraction.where(land_fraction <= 1.0, 1.0)
                water_proportion = 1.0 - land_fraction

                grid_areas = None
                run_pbar = tqdm(runs, desc=f"Processing {storm_name}", leave=False)

                for run_name in run_pbar:
                    storm_run_path = path.join(MAXSS_working_directory, "output", run_name, "maxss/storm-atlas/tropical/ibtracs", region, year_name, storm_name)
                    fluxengine_files = sorted(glob(os.path.join(storm_run_path, "*.nc")))

                    # --- RESTORED: Missing run skip print statement ---
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
                        mass_map = ds.OF * grid_areas * (1.0 / 24.0) * current_water_map
                        hourly_fluxes_tg = mass_map.sum(dim=['latitude', 'longitude'], skipna=True) / 1e12
                        
                        hourly_fluxes_arr = hourly_fluxes_tg.values
                        time_vals = ds.time.values
                        storm_total_tg = float(hourly_fluxes_tg.sum(skipna=True))

                        # --- SAVE NETCDF INSIDE CONTEXT MANAGER ---
                        out_folder = path.join(netcdf_output_root, region, year_name)
                        os.makedirs(out_folder, exist_ok=True)
                        processedFilePath = os.path.join(out_folder, f"{storm_name}_{run_name}.nc")

                        ds_out = xr.Dataset(
                            data_vars={
                                "Hourly_flux": (["time"], hourly_fluxes_arr, {"units": "Tg C hr-1"}),
                                "Total_flux": ([], storm_total_tg, {"units": "Tg C"})
                            },
                            coords={"time": time_vals}
                        )
                        ds_out['time'].encoding['units'] = "seconds since 1970-01-01 00:00:00"
                        ds_out['time'].encoding['dtype'] = 'i8'
                        ds_out.to_netcdf(processedFilePath)

                    # --- RESTORED: Run success summary print statement ---
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
        
        output_path = path.join(output_base, "storm_component_flux_summary.csv")
        os.makedirs(output_base, exist_ok=True)
        df_final.to_csv(output_path, index=False)
        
        print(f"\n{'='*40}")
        print(f"SUCCESS: Summary saved to {output_path}")
        print(f"{'='*40}")
    else:
        print("\nNo data was processed. Check your directory paths.")