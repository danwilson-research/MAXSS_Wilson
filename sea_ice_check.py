import os
from glob import glob
import numpy as np
import netCDF4 as nc
from netCDF4 import num2date
from tqdm import tqdm

# ==============================================================================
# Configuration
# ==============================================================================
base_dir = "/nobackup/beegfs/workspace/dw557/MAXSS_project/MAXSS_working_directory/maxss/storm-atlas/tropical/ibtracs/north-atlantic"
VAR_NAME = "__eo_sea_ice_fraction"
ICE_TOLERANCE = 1e-5  # Avoids false positives from tiny floating-point noise

# Find all raw L4 wind files across all years and storms
raw_l4_files = sorted(glob(os.path.join(base_dir, "*", "*", "*_MAXSS_HIST_TC_L4.nc")))
print(f"Found {len(raw_l4_files)} raw L4 files to inspect.\n")

storms_with_ice = []
missing_var_files = []

for filepath in tqdm(raw_l4_files, desc="Scanning sea ice"):
    storm_folder = os.path.basename(os.path.dirname(filepath))
    year_folder = os.path.basename(os.path.dirname(os.path.dirname(filepath)))

    with nc.Dataset(filepath, "r") as ds:
        if VAR_NAME not in ds.variables:
            missing_var_files.append(filepath)
            continue

        ice_var = ds.variables[VAR_NAME]
        raw_data = ice_var[:]

        # Handle masked arrays and explicit _FillValue (1.0E20f)
        if np.ma.isMaskedArray(raw_data):
            ice_data = raw_data.filled(np.nan)
        else:
            ice_data = np.array(raw_data, dtype=np.float64)

        fill_val = getattr(ice_var, "_FillValue", 1.0e20)
        ice_data = np.where((ice_data >= fill_val * 0.9) | (~np.isfinite(ice_data)), np.nan, ice_data)

        # Find timesteps where any valid pixel has sea ice > ICE_TOLERANCE
        # Suppress All-NaN slice warnings for land-only timesteps if any exist
        with np.errstate(all="ignore"):
            max_ice_per_step = np.nanmax(ice_data, axis=(1, 2))

        affected_indices = np.where(max_ice_per_step > ICE_TOLERANCE)[0]

        if len(affected_indices) > 0:
            time_var = ds.variables["time"]
            dates = num2date(time_var[:], units=time_var.units)

            storms_with_ice.append((year_folder, storm_folder, len(affected_indices), float(np.nanmax(max_ice_per_step))))

            print(f"\n[SEA ICE DETECTED] {year_folder} / {storm_folder} ({len(affected_indices)} timesteps):")
            for idx in affected_indices:
                t_str = dates[idx].strftime("%Y-%m-%d %H:%M:%S")
                n_pixels = int(np.sum(ice_data[idx] > ICE_TOLERANCE))
                print(f"   - Step {idx:4d} | {t_str} | Max ice fraction: {max_ice_per_step[idx]:.4f} ({n_pixels} pixels)")

# ==============================================================================
# Final Summary Report
# ==============================================================================
print("\n" + "=" * 60)
print(f"SEA ICE CHECK COMPLETE ({len(raw_l4_files)} files scanned)")
print("=" * 60)

if not storms_with_ice:
    print("RESULT: Zero storms had sea ice > 0 at any timestep.")
else:
    print(f"RESULT: {len(storms_with_ice)} storm(s) had non-zero sea ice:")
    for yr, st, n_steps, max_val in storms_with_ice:
        print(f"  * {yr} - {st}: {n_steps} timesteps (peak fraction = {max_val:.4f})")

if missing_var_files:
    print(f"\nNote: {len(missing_var_files)} file(s) did not contain '{VAR_NAME}'.")
