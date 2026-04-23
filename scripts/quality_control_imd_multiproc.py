#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quality_control_imd_multiproc.py
=================================

Batch quality-control pipeline for IMD DWR Level-2 NetCDF files.

Data source
-----------
MOSDAC — Meteorological and Oceanographic Satellite Data Archival Centre
https://mosdac.gov.in/

Usage
-----
Place raw IMD NetCDF (.nc) files in ``<repo_root>/INPUT_DATA/`` and run::

    python scripts/quality_control_imd_multiproc.py

Output directories (VEL/, DBZ/) are created automatically relative to the
repository root.

Processing steps
----------------
For each NetCDF file:

1. Read and parse file metadata (scan mode, polarisation type, dimensions).
2. Apply clutter filtering to reflectivity (Z):
   - Single-pol: Gabella statistical filter.
   - Dual-pol:   Fuzzy-logic polarimetric classifier (ρ_hv, Φ_DP, Z_DR, V).
3. Apply custom bidirectional velocity dealiasing to radial velocity (V)
   for short-range scans.
4. Save corrected fields as space-separated text arrays per elevation sweep.

Tunable parameters
------------------
NUM_WORKERS : int
    Number of parallel worker processes. Default: 6.
ALPHA : float
    Shear threshold fraction of the Nyquist velocity used in dealiasing.
    Default: 0.8.
BETA : float
    Low-velocity threshold fraction of v_n for reference radial selection.
    Default: 0.3.
"""

import datetime as dt
import multiprocessing
import os
import warnings
from pathlib import Path

import numpy as np
import wradlib as wrl

from dwr_processing.io.netcdf_reader import (
    read_imd_netcdf,
    detect_polarisation_type,
    get_scan_mode,
)
from dwr_processing.clutter.filters import gabella_filter, polarimetric_filter
from dwr_processing.velocity.dealias import dealias_velocity
from dwr_processing.utils.sorting import sort_nicely

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parent.parent
RAWFILES_LOC = _REPO_ROOT / "INPUT_DATA"

OUT_VEL_LONG  = _REPO_ROOT / "VEL" / "long_range"
OUT_VEL_SHORT = _REPO_ROOT / "VEL" / "short_range"
OUT_DBZ_LONG  = _REPO_ROOT / "DBZ" / "long_range"
OUT_DBZ_SHORT = _REPO_ROOT / "DBZ" / "short_range"

for _d in (OUT_VEL_LONG, OUT_VEL_SHORT, OUT_DBZ_LONG, OUT_DBZ_SHORT):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Tunable parameters
# ---------------------------------------------------------------------------
NUM_WORKERS = 6
ALPHA = 0.8   # velocity shear threshold (fraction of Nyquist)
BETA  = 0.3   # low-velocity gate fraction for reference radial selection

# ---------------------------------------------------------------------------
# Collect and sort input files
# ---------------------------------------------------------------------------
nc_files = sort_nicely(
    [f for f in os.listdir(RAWFILES_LOC) if f.endswith(".nc")]
)

if not nc_files:
    raise FileNotFoundError(
        f"No NetCDF files found in {RAWFILES_LOC}. "
        "Please download IMD DWR data from https://mosdac.gov.in/ and place "
        "the .nc files in the INPUT_DATA/ directory."
    )

print(f"Found {len(nc_files)} NetCDF file(s) to process.")


# ---------------------------------------------------------------------------
# Per-file processing function
# ---------------------------------------------------------------------------
def process_raw_file(filename: str) -> None:
    """
    Quality-control a single IMD DWR NetCDF file.

    Applies clutter filtering (reflectivity) and velocity dealiasing
    (short-range scans only), then writes corrected fields to disk.

    Parameters
    ----------
    filename : str
        Basename of the NetCDF file located in ``INPUT_DATA/``.
    """
    filepath = RAWFILES_LOC / filename
    print(f"Processing: {filename}")

    raw = read_imd_netcdf(str(filepath))
    pol_type = detect_polarisation_type(raw)

    # --- Parse filename for output basename ---
    # Expected format: <PREFIX>_<SITE>_<YYYYMMDDHHMMSS>.nc  (adjust as needed)
    time_str  = filename[16:30]
    basename  = filename[11:16] + time_str

    # --- Scan geometry ---
    v = raw["variables"]
    d = raw["dimensions"]

    num_elev  = d["sweep"]["size"]
    num_bins  = d["bin"]["size"]
    num_azim  = d["radial"]["size"]
    range_res = float(v["gateSize"]["data"])
    sweep     = v["elevationNumber"]["data"]
    elev      = v["elevationAngle"]["data"]

    try:
        mode = get_scan_mode(num_elev)
    except ValueError as exc:
        print(f"  [SKIP] {exc}")
        return

    # --- Output directories for this mode ---
    vel_outdir = OUT_VEL_SHORT if mode == "short_range" else OUT_VEL_LONG
    dbz_outdir = OUT_DBZ_SHORT if mode == "short_range" else OUT_DBZ_LONG

    if mode != "short_range":
        # Long-range scans: apply only clutter filtering, no dealiasing
        dbz = v["Z"]["data"] if "Z" in v else None
        if dbz is not None and dbz.size != 0:
            rdata = gabella_filter(dbz) if pol_type == "single_pol" else gabella_filter(dbz)
            np.savetxt(
                dbz_outdir / f"{mode}_{basename}_corr_dbz_sweep_{sweep:02d}.txt",
                rdata,
            )
        return

    # --- Short-range: full QC pipeline ---
    vel    = v["V"]["data"]    if "V"    in v else None
    dbz    = v["Z"]["data"]    if "Z"    in v else None
    sigma  = v["W"]["data"]    if "W"    in v else None
    zdr    = v["ZDR"]["data"]  if "ZDR"  in v else None
    phidp  = v["PHIDP"]["data"] if "PHIDP" in v else None
    rhohv  = v["RHOHV"]["data"] if "RHOHV" in v else None
    v_n    = float(v["nyquist"]["data"])
    fill_value = v["V"].get("_FillValue") if "V" in v else None

    # --- Clutter filtering (reflectivity) ---
    if dbz is not None and dbz.size != 0 and num_azim == 360:
        if pol_type == "single_pol":
            rdata = gabella_filter(dbz)
        else:
            rdata = polarimetric_filter(dbz, rhohv, phidp, vel, zdr)
        np.savetxt(
            dbz_outdir / f"{mode}_{basename}_corr_dbz_sweep_{sweep:02d}.txt",
            rdata,
        )
    else:
        rdata = np.full((360, num_bins), np.nan)
        np.savetxt(
            dbz_outdir / f"{mode}_{basename}_corr_dbz_sweep_{sweep:02d}.txt",
            rdata,
        )

    # --- Velocity dealiasing ---
    if vel is not None and vel.size != 0 and num_azim == 360:
        # Mask clutter-flagged pixels in velocity
        if dbz is not None and dbz.size != 0:
            vel[np.where(np.isnan(rdata))] = np.nan
        if fill_value is not None:
            vel[vel == fill_value] = np.nan

        corr_vel = dealias_velocity(
            velocity=vel,
            v_n=v_n,
            range_resolution=range_res,
            alpha=ALPHA,
            beta=BETA,
        )
        np.savetxt(
            vel_outdir / f"{mode}_{basename}_corr_vel_sweep_{sweep:02d}.txt",
            corr_vel,
        )
        print(f"  Done: {filename} | mode={mode} | pol={pol_type} | elev={np.round(elev, 1)}")
    else:
        # Incomplete sweep — fill with NaN
        v_nan = np.full((360, num_bins), np.nan)
        np.savetxt(
            vel_outdir / f"{mode}_{basename}_corr_vel_sweep_{sweep:02d}.txt",
            v_nan,
        )


# ---------------------------------------------------------------------------
# Multiprocessing entry point
# ---------------------------------------------------------------------------
def mp_handler() -> None:
    """Distribute file processing across ``NUM_WORKERS`` parallel workers."""
    with multiprocessing.Pool(processes=NUM_WORKERS) as pool:
        pool.map(process_raw_file, nc_files)


if __name__ == "__main__":
    t_start = dt.datetime.now()
    mp_handler()
    elapsed = dt.datetime.now() - t_start
    print(f"\nQuality control complete. Total time: {elapsed}")
