# -*- coding: utf-8 -*-

#Import required packages
import os
from os import path
from glob import glob
import numpy as np
import xarray as xr
from pyproj import Geod
from tqdm import tqdm
import sys

# ==============================================================================
# USER CONFIGURATION
# ==============================================================================
# Usage: python grid_fluxes.py Deike   (defaults to "Deike" if no argument is given)
PARAMETERIZATION = sys.argv[1] if len(sys.argv) > 1 else "Deike"

TEST_MODE = True # True = 3 storms on REF_RUN & MAXSS_RUN; False = all storms & all 8 runs

MAXSS_working_directory = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory"

PARAM_CONFIG = {
    "Deike": {
        "runs_root": path.join(MAXSS_working_directory, "Deike_output"),
        "output_dir": path.join(MAXSS_working_directory, "Deike_output", "Spatially_integrated_fluxes"),
    },
    "Nightingale": {
        "runs_root": path.join(MAXSS_working_directory, "Nightingale-output"),
        "output_dir": path.join(MAXSS_working_directory, "Nightingale-output", "Spatially_integrated_fluxes"),
    },
}

ALL_RUNS = [
    "REF_RUN",
    "SSS_RUN",
    "SST_WITH_GRADIENTS_RUN",
    "WIND_RUN",
    "MAXSS_RUN",
    "PRESSURE_RUN",
    "SST_NO_GRADIENTS_RUN",
    "V_GAS_RUN",
]
# ==============================================================================


def calculate_ellipsoidal_areas(lat_values, lon_sample, res=0.25):
    """Compute exact WGS84 ellipsoidal area (m^2) for each latitude band."""
    geod = Geod(ellps="WGS84")
    row_areas = []
    for lat in lat_values:
        lats = [lat - res / 2, lat - res / 2, lat + res / 2, lat + res / 2, lat - res / 2]
        lons = [lon_sample - res / 2, lon_sample + res / 2, lon_sample + res / 2, lon_sample - res / 2, lon_sample - res / 2]
        poly_area, _ = geod.polygon_area_perimeter(lons, lats)
        row_areas.append(abs(poly_area))
    return xr.DataArray(row_areas, coords={"latitude": lat_values}, dims="latitude")


def build_aligned_master_grid(sample_ds, lat_bounds=(-15.0, 85.0), lon_bounds=(-125.0, 30.0), res=0.25):
    """
    Build a stationary 0.25-deg North Atlantic master grid whose coordinate centers
    match the exact sub-pixel registration of the FluxEngine storm files.
    """
    lat_offset = float(sample_ds.latitude.values[0]) % res
    lon_offset = float(sample_ds.longitude.values[0]) % res

    lat_start = np.floor(lat_bounds[0] / res) * res + lat_offset
    lat_end = np.ceil(lat_bounds[1] / res) * res + lat_offset
    lon_start = np.floor(lon_bounds[0] / res) * res + lon_offset
    lon_end = np.ceil(lon_bounds[1] / res) * res + lon_offset

    master_lat = np.round(np.arange(lat_start, lat_end + res / 2, res), 6)
    master_lon = np.round(np.arange(lon_start, lon_end + res / 2, res), 6)
    return master_lat, master_lon


def generate_gridded_component_fluxes(
    MAXSS_working_directory,
    fluxengine_runs_root,
    output_nc_path,
    runs,
    parameterization_label="Deike",
    MAXSS_regions=("north-atlantic",),
    salinity_threshold=10.0,
    max_storms=None,
):
    # Collect all storm directories across years
    storm_entries = []
    for region in MAXSS_regions:
        ref_dir = path.join(fluxengine_runs_root, "MAXSS_RUN/maxss/storm-atlas/tropical/ibtracs", region)
        for year_dir in sorted(glob(path.join(ref_dir, "*/"))):
            year_name = os.path.basename(os.path.normpath(year_dir))
            for storm_path in sorted(glob(path.join(year_dir, "*/"))):
                storm_name = os.path.basename(os.path.normpath(storm_path))
                storm_entries.append((region, year_name, storm_name))

    if max_storms is not None:
        storm_entries = storm_entries[:max_storms]
        print(f"[{parameterization_label}] TEST MODE ACTIVE: Limiting run to {len(storm_entries)} storms and {len(runs)} runs: {runs}")
    else:
        print(f"[{parameterization_label}] FULL RUN: Processing {len(storm_entries)} storms across {len(runs)} runs.")

    if not storm_entries:
        raise RuntimeError(f"No storm directories found under {fluxengine_runs_root}.")

    # Initialize master grid from the first available FluxEngine file
    master_lat, master_lon = None, None
    for region, year_name, storm_name in storm_entries:
        for run_name in runs:
            sample_files = sorted(glob(path.join(
                fluxengine_runs_root, run_name,
                "maxss/storm-atlas/tropical/ibtracs", region, year_name, storm_name, "*.nc"
            )))
            if sample_files:
                with xr.open_dataset(sample_files[0]) as sample_ds:
                    master_lat, master_lon = build_aligned_master_grid(sample_ds, res=0.25)
                break
        if master_lat is not None:
            break

    if master_lat is None:
        raise RuntimeError("Could not find any FluxEngine NetCDF files to initialize master grid.")

    # Accumulators: cumulative flux (Tg C) and number of valid storm-hours per grid cell
    cumulative_flux = {
        run_name: np.zeros((len(master_lat), len(master_lon)), dtype=np.float64)
        for run_name in runs
    }
    valid_hour_counts = {
        run_name: np.zeros((len(master_lat), len(master_lon)), dtype=np.int32)
        for run_name in runs
    }

    # Main loop over storms
    for region, year_name, storm_name in tqdm(storm_entries, desc=f"Accumulating {parameterization_label} Grids"):
        land_fraction_path = path.join(
            MAXSS_working_directory,
            "maxss/storm-atlas/tropical/ibtracs",
            region,
            year_name,
            storm_name,
            "Resampled_for_fluxengine_MAXSS_land_fraction.nc",
        )

        with xr.open_dataset(land_fraction_path) as lf_ds:
            land_fraction = lf_ds.land_proportion.load()
            if "sea_ice_proportion" in lf_ds:
                sea_ice_proportion = lf_ds.sea_ice_proportion.load()
            else:
                print(f"   -> Warning: 'sea_ice_proportion' missing in {land_fraction_path}; assuming 0 ice.")
                sea_ice_proportion = None

        land_fraction = land_fraction.where(land_fraction <= 1.0, 1.0)
        water_proportion = 1.0 - land_fraction

        grid_areas = None

        for run_name in runs:
            storm_run_path = path.join(
                fluxengine_runs_root,
                run_name,
                "maxss/storm-atlas/tropical/ibtracs",
                region,
                year_name,
                storm_name,
            )
            fluxengine_files = sorted(glob(path.join(storm_run_path, "*.nc")))

            if not fluxengine_files:
                print(f"   -> Run {run_name} not found for {storm_name}. Skipping.")
                continue

            with xr.open_mfdataset(fluxengine_files, combine="by_coords", chunks=None, parallel=False) as ds:
                if grid_areas is None:
                    grid_areas = calculate_ellipsoidal_areas(ds.latitude.values, float(ds.longitude.values[0]), res=0.25)

                # 1. Align land and sea-ice masks onto the storm's FluxEngine grid
                current_water_map = water_proportion.reindex_like(ds, method="nearest")

                if sea_ice_proportion is not None:
                    current_ice_map = (
                        sea_ice_proportion.reindex_like(ds, method="nearest")
                        .fillna(0.0)
                        .clip(min=0.0, max=1.0)
                    )
                    ice_free_proportion = 1.0 - current_ice_map
                else:
                    ice_free_proportion = 1.0

                # 2. Enforce FluxEngine's valid_min / valid_max
                valid_min = ds.OF.attrs["valid_min"]
                valid_max = ds.OF.attrs["valid_max"]
                clean_OF = ds.OF.where(np.isfinite(ds.OF) & (ds.OF >= valid_min) & (ds.OF <= valid_max))

                # 3. Exclude estuarine pixels below salinity_threshold (10 PSU)
                if salinity_threshold is not None:
                    is_estuarine = xr.zeros_like(clean_OF, dtype=bool)
                    if "pco2_sss" in ds:
                        is_estuarine = is_estuarine | (ds.pco2_sss < salinity_threshold) | (~np.isfinite(ds.pco2_sss))
                    if "OKS1" in ds:
                        is_estuarine = is_estuarine | (ds.OKS1 < salinity_threshold) | (~np.isfinite(ds.OKS1))
                    clean_OF = clean_OF.where(~is_estuarine)

                # 4. Integrate hourly flux over the storm lifetime at each 0.25-deg pixel (Tg C)
                mass_map_tg = (clean_OF * grid_areas * (1.0 / 24.0) * current_water_map * ice_free_proportion) / 1e12
                storm_total_map = mass_map_tg.sum(dim="time", skipna=True)
                storm_valid_hours = np.isfinite(mass_map_tg).sum(dim="time")

                # 5. Verify the storm box fits completely inside the master grid
                if (
                    float(storm_total_map.latitude.min()) < master_lat[0]
                    or float(storm_total_map.latitude.max()) > master_lat[-1]
                    or float(storm_total_map.longitude.min()) < master_lon[0]
                    or float(storm_total_map.longitude.max()) > master_lon[-1]
                ):
                    print(f"   -> WARNING: {storm_name} ({run_name}) extends outside master grid bounds!")

                # 6. Snap onto stationary master grid (exact coordinate match via tolerance)
                aligned_flux = (
                    storm_total_map.reindex(
                        latitude=master_lat, longitude=master_lon, method="nearest", tolerance=1e-4, fill_value=0.0
                    )
                    .fillna(0.0)
                    .values
                )
                aligned_counts = (
                    storm_valid_hours.reindex(
                        latitude=master_lat, longitude=master_lon, method="nearest", tolerance=1e-4, fill_value=0
                    )
                    .fillna(0)
                    .values
                )

                cumulative_flux[run_name] += aligned_flux
                valid_hour_counts[run_name] += aligned_counts.astype(np.int32)

    # Build output NetCDF dataset
    data_vars = {}
    encoding = {}
    thresh_attr = f"{salinity_threshold} PSU" if salinity_threshold is not None else "None"

    for run_name in runs:
        # Mask cells that were never observed by any storm as NaN (land / outside domain)
        gridded_total = np.where(valid_hour_counts[run_name] > 0, cumulative_flux[run_name], np.nan)

        data_vars[run_name] = (
            ["latitude", "longitude"],
            gridded_total,
            {
                "units": "Tg C",
                "long_name": f"Cumulative TC-induced air-sea CO2 flux ({run_name}, {parameterization_label})",
                "parameterization": parameterization_label,
                "salinity_threshold": thresh_attr,
                "sea_ice_scaling": "linear (1 - sea_ice_proportion)",
            },
        )
        data_vars[f"{run_name}_valid_hours"] = (
            ["latitude", "longitude"],
            valid_hour_counts[run_name],
            {
                "units": "hours",
                "long_name": f"Number of valid storm-hours contributing to {run_name}",
            },
        )
        encoding[run_name] = {"zlib": True, "complevel": 4, "_FillValue": np.nan}
        encoding[f"{run_name}_valid_hours"] = {"zlib": True, "complevel": 4}

    ds_out = xr.Dataset(
        data_vars=data_vars,
        coords={
            "latitude": (["latitude"], master_lat, {"units": "degrees_north", "standard_name": "latitude"}),
            "longitude": (["longitude"], master_lon, {"units": "degrees_east", "standard_name": "longitude"}),
        },
        attrs={
            "title": f"Gridded Cumulative Tropical Cyclone CO2 Flux Components (0.25 deg) - {parameterization_label}",
            "parameterization": parameterization_label,
            "storms_processed": len(storm_entries),
            "salinity_threshold": thresh_attr,
        },
    )

    os.makedirs(path.dirname(output_nc_path), exist_ok=True)
    ds_out.to_netcdf(output_nc_path, encoding=encoding)

    print(f"\n{'='*50}")
    print(f"SUCCESS [{parameterization_label}]: Gridded NetCDF saved to:\n  {output_nc_path}")
    for run_name in runs:
        total_basin_tg = float(np.nansum(ds_out[run_name].values))
        print(f"  - {run_name}: {total_basin_tg:.6f} Tg C (basin sum)")
    print(f"{'='*50}")


if __name__ == "__main__":

    cfg = PARAM_CONFIG[PARAMETERIZATION]

    if TEST_MODE:
        runs_to_run = ["REF_RUN", "MAXSS_RUN"]
        max_storms = 3
        output_file = path.join(cfg["output_dir"], f"gridded_component_total_flux_{PARAMETERIZATION}_TEST_3storms.nc")
    else:
        runs_to_run = ALL_RUNS
        max_storms = None
        output_file = path.join(cfg["output_dir"], f"gridded_component_total_flux_{PARAMETERIZATION}_all_TCs.nc")

    generate_gridded_component_fluxes(
        MAXSS_working_directory=MAXSS_working_directory,
        fluxengine_runs_root=cfg["runs_root"],
        output_nc_path=output_file,
        runs=runs_to_run,
        parameterization_label=PARAMETERIZATION,
        MAXSS_regions=("north-atlantic",),
        salinity_threshold=10.0,
        max_storms=max_storms,
    )
