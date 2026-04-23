# 🌩️ DWR Processing Toolkit

A Python toolkit for reading, quality-controlling, and processing **Doppler Weather Radar (DWR)** data from the Indian region (currently for datasets availaible from MOSDAC). The toolkit handles both **single-polarisation** and **dual-polarisation** radar products and is designed to work with the NetCDF format distributed by [MOSDAC (Meteorological and Oceanographic Satellite Data Archival Centre)](https://mosdac.gov.in/).

---

## ✨ Features

| Module | Description |
|---|---|
| **Clutter Filtering** | Gabella statistical filter (single-pol) and fuzzy-logic polarimetric classifier (dual-pol) |
| **Velocity Dealiasing** | Bidirectional (clockwise + anti-clockwise) radial-by-radial phase-unfolding algorithm |
| **Parallel Processing** | `multiprocessing`-based batch pipeline for large archives |
| **I/O** | NetCDF ingestion via `wradlib`; outputs corrected fields as text arrays |
| **Auto-mode detection** | Automatically distinguishes long-range (≤3 elevation sweeps) and short-range (10 sweeps) scan strategies |

---

## 📂 Repository Structure

```
dwr-qc/
├── dwr_processing/             # Core Python package
│   ├── __init__.py
│   ├── io/                     # Data ingestion and file utilities
│   │   ├── __init__.py
│   │   └── netcdf_reader.py
│   ├── clutter/                # Clutter identification and removal
│   │   ├── __init__.py
│   │   └── filters.py
│   ├── velocity/               # Velocity dealiasing
│   │   ├── __init__.py
│   │   └── dealias.py
│   └── utils/                  # Shared utilities
│       ├── __init__.py
│       └── sorting.py
├── scripts/
│   └── quality_control_imd_multiproc.py   # Batch QC pipeline (entry point)
├── notebooks/
│   └── 01_quickstart_demo.ipynb            # Quickstart walkthrough
├── data/                       # Place raw NetCDF files here (not tracked)
│   └── README.md
├── docs/
│   └── algorithm_notes.md      # Notes on dealiasing and clutter algorithms
├── requirements.txt
├── .gitignore
├── LICENSE
└── README.md
```

---

## 🛰️ Data

Indian DWR Level-2 NetCDF data are distributed by **MOSDAC** (Meteorological and Oceanographic Satellite Data Archival Centre), operated by the Space Applications Centre (SAC), ISRO.

- **Portal:** [https://mosdac.gov.in/](https://mosdac.gov.in/)
- **Product:** IMD Doppler Weather Radar – Level-2 (Base/Derived moments in NetCDF)
- Access requires registration on the MOSDAC portal. Dataset links can be updated in `data/README.md`.

> **Citation:** If you use IMD DWR data in your research, please acknowledge MOSDAC:
> *"DWR Level-2 data were obtained from the MOSDAC data portal (https://mosdac.gov.in/), operated by SAC/ISRO, India."*

---

## ⚙️ Installation

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/dwr-qc.git
cd dwr-qc
```

### 2. Create a virtual environment (recommended)

```bash
python3 -m venv .venv
source .venv/bin/activate        # Linux / macOS
# .venv\Scripts\activate         # Windows
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

> **GDAL / osgeo note:** `gdal` often requires system-level libraries. On Ubuntu/Debian:
> ```bash
> sudo apt-get install libgdal-dev
> pip install gdal==$(gdal-config --version)
> ```
> On macOS with Homebrew: `brew install gdal && pip install gdal`.

---

## 🚀 Quick Start

### Batch quality control (recommended)

Place your raw IMD NetCDF files in `INPUT_DATA/` (or update the path in the script), then:

```bash
python scripts/quality_control_imd_multiproc.py
```

Output directories are created automatically:

```
VEL/
  long_range/    ← corrected radial velocity (long-range scans)
  short_range/   ← corrected radial velocity (short-range scans)
DBZ/
  long_range/    ← clutter-filtered reflectivity (long-range scans)
  short_range/   ← clutter-filtered reflectivity (short-range scans)
```

Each output file is a space-separated text array of shape `(num_azimuth × num_gates)` for a single elevation sweep.

### Using the package directly

```python
import wradlib as wrl
from dwr_processing.clutter.filters import gabella_filter, polarimetric_filter
from dwr_processing.velocity.dealias import dealias_velocity
from dwr_processing.io.netcdf_reader import read_imd_netcdf

# Load a single file
raw = read_imd_netcdf("INPUT_DATA/IMD_RADAR_20200317_000000.nc")

# Apply clutter filter (single-pol example)
dbz = raw["variables"]["Z"]["data"]
dbz_clean = gabella_filter(dbz)

# Dealias velocity
vel = raw["variables"]["V"]["data"]
nyquist = raw["variables"]["nyquist"]["data"]
vel_corrected = dealias_velocity(vel, nyquist, range_resolution=raw["variables"]["gateSize"]["data"])
```

---

## 🧠 Algorithm Details

### Clutter Filtering

| Radar Type | Method |
|---|---|
| Single-polarisation | **Gabella filter** (`wradlib.clutter.filter_gabella`) — texture-based statistical clutter identifier |
| Dual-polarisation | **Fuzzy-logic echo classifier** (`wradlib.clutter.classify_echo_fuzzy`) using ρ_hv, Φ_DP, Z_DR, Doppler velocity, and the Gabella clutter map as features |

After classification, isolated speckle is removed with `wradlib.dp.linear_despeckle`.

### Velocity Dealiasing

A custom two-pass bidirectional dealiasing algorithm:

1. **Initial radial selection** — identifies a reference radial near the zero-velocity isodop using sign-change detection in azimuthally-averaged mean velocities.
2. **First pass** — gate-by-gate phase correction propagated simultaneously clockwise and anti-clockwise from the reference radial, in both azimuthal and radial directions. The Nyquist interval correction integer *k* is estimated as `round((v_ref − v_obs) / (2 × v_n))`.
3. **Second pass** — a second sweep over unresolved gates (flag = 0) using an expanded search window (*M* = 3 neighbouring radials) for residual aliasing.

Key parameters (adjustable at the top of `scripts/quality_control_imd_multiproc.py`):

| Parameter | Default | Description |
|---|---|---|
| `alpha` | 0.8 | Shear threshold as a fraction of Nyquist velocity |
| `beta` | 0.3 | Low-velocity gate fraction for reference radial selection |

---

## 📋 Requirements

See [`requirements.txt`](requirements.txt) for the full pinned list. Core dependencies:

- **[wradlib](https://docs.wradlib.org/)** ≥ 1.19 — radar I/O, clutter filtering, polarimetric processing
- **[Py-ART](https://arm-doe.github.io/pyart/)** ≥ 1.12 — ARM/DOE radar toolkit (used in extended analysis notebooks)
- **NumPy**, **SciPy**, **Matplotlib** — numerical and plotting stack
- **GDAL / osgeo** — geospatial coordinate handling
- **pathlib** — standard library (Python ≥ 3.4)

---

## 📓 Notebooks

| Notebook | Description |
|---|---|
| [`01_quickstart_demo.ipynb`](notebooks/01_quickstart_demo.ipynb) | Load a single IMD NetCDF file, apply QC, and visualise PPI plots |

---

## 🤝 Contributing

Contributions, bug reports, and feature requests are welcome. Please open an issue or submit a pull request. When contributing code, please follow PEP 8 and include docstrings for any new functions.

---

## 📄 License

This project is licensed under the MIT License. See [`LICENSE`](LICENSE) for details.

---

## 🙏 Acknowledgements

- **MOSDAC / SAC / ISRO** for providing IMD DWR Level-2 data via the [MOSDAC portal](https://mosdac.gov.in/).
- The **[wradlib](https://docs.wradlib.org/)** development team for the open-source radar processing library.
- The **[Py-ART](https://arm-doe.github.io/pyart/)** team at ARM/DOE.
