# Data

Raw IMD DWR Level-2 NetCDF files are **not tracked** in this repository (they are excluded via `.gitignore`).

## Where to obtain data

Indian DWR data are distributed by **MOSDAC** (Meteorological and Oceanographic Satellite Data Archival Centre), operated by the Space Applications Centre (SAC), ISRO.

| Resource | Link |
|---|---|
| MOSDAC Portal | https://mosdac.gov.in/ |
| DWR Product Page | *(update with direct product link once available)* |

Registration on the MOSDAC portal is required for data access.

## Directory layout

Place your downloaded NetCDF files in a directory called `INPUT_DATA/` at the root of the repository (or update the path variable in `scripts/quality_control_imd_multiproc.py`):

```
dwr-qc/
└── INPUT_DATA/
    ├── IMD_RADAR_<SITE>_<YYYYMMDD>_<HHMMSS>.nc
    └── ...
```

## Expected NetCDF variable names

The pipeline expects the following standard IMD DWR variable names inside each NetCDF file:

| Variable | Description |
|---|---|
| `Z` | Reflectivity (dBZ) |
| `V` | Radial velocity (m/s) |
| `W` | Spectrum width (m/s) |
| `ZDR` | Differential reflectivity (dB) — dual-pol only |
| `PHIDP` | Differential phase (°) — dual-pol only |
| `RHOHV` | Co-polar correlation coefficient — dual-pol only |
| `siteLat` | Radar site latitude |
| `siteLon` | Radar site longitude |
| `siteAlt` | Radar site altitude (m) |
| `gateSize` | Range gate size (m) |
| `radialAzim` | Azimuth angles (°) |
| `radialElev` | Elevation angles (°) |
| `nyquist` | Nyquist velocity (m/s) |
| `elevationAngle` | Elevation angle for the sweep (°) |
| `elevationNumber` | Sweep index |
| `elevationList` | List of elevation angles |
| `firstGateRange` | Range to first gate (m) |
| `unambigRange` | Maximum unambiguous range (km) |
