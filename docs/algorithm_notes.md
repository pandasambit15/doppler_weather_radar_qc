# Algorithm Notes

## 1. Clutter Filtering

### 1a. Gabella Filter (Single-Polarisation)

Applied when only base moments (Z, V, W) are available.

The **Gabella filter** (`wradlib.clutter.filter_gabella`) is a texture-based statistical method that identifies non-meteorological echoes by computing a local texture metric over a sliding window:

- **Window size (`wsize`):** 5 gates × 5 radials
- **No-rain threshold (`thrsnorain`):** 0.0 dBZ
- **Texture ratio threshold (`tr1`):** 10.0
- **Minimum fraction of non-rain neighbours (`n_p`):** 8
- **Second texture threshold (`tr2`):** 1.3

Gates flagged as clutter are set to NaN.

**Reference:**  
Gabella, M., & Notarpietro, R. (2002). Ground clutter characterization and elimination in mountainous terrain. *Proceedings of ERAD*, 305–311.

---

### 1b. Fuzzy-Logic Polarimetric Filter (Dual-Polarisation)

Applied when ZDR, ΦDKP, and ρ_hv are all present and non-empty.

Uses `wradlib.clutter.classify_echo_fuzzy` with the following feature weights:

| Feature | Weight |
|---|---|
| ρ_hv (`rho`) | 0.5 |
| Φ_DP (`phi`) | 0.3 |
| Reflectivity (`ref`, post-Gabella) | — |
| Radial velocity (`dop`) | 0.3 |
| Z_DR (`zdr`) | 0.5 |
| ρ_hv² (`rho2`) | 0.3 |
| Gabella clutter map (`map`) | 0.5 |

Classification threshold: **0.5** (gates with fuzzy score ≥ 0.5 are clutter).

After classification, `wradlib.dp.linear_despeckle` with a 5-gate window removes isolated non-meteorological pixels.

---

## 2. Velocity Dealiasing

The algorithm corrects **range-velocity ambiguity** (aliasing) that occurs when the true radial velocity exceeds the **Nyquist velocity** (v_n). Aliased velocities are folded by multiples of 2v_n. The correction integer *k* satisfies:

```
v_corrected = v_observed + k × 2 × v_n
```

where *k* is chosen to minimise the difference between the observed gate velocity and a locally-derived reference velocity *v_ref*.

### Phase 0 — Initial Radial Selection

The algorithm requires a **reference radial** near the zero-velocity isodop (the radial where mean Doppler velocity transitions from positive to negative, or vice versa). This radial is unlikely to contain aliased velocities, making it a safe starting point.

**Search procedure:**

1. For each azimuth, radials with large velocity shears (`> α × v_n`) are flagged as unreliable.
2. Among reliable radials, gates with small absolute velocities (`< β × v_n`) are collected.
3. The zero-crossing transition (sign flip in the mean velocity v_m) over 4 consecutive reliable azimuths identifies the reference radial.
4. If no zero-crossing is found, a relaxed search is performed with progressively looser gate-count requirements (starting from 40, stepping down by 5).

Key parameters:

| Parameter | Symbol | Default | Role |
|---|---|---|---|
| `alpha` | α | 0.8 | Shear threshold fraction of v_n |
| `beta` | β | 0.3 | Low-velocity gate fraction |

### Phase 1 — First-Pass Dealiasing

Starting from the reference radial, velocity correction propagates in **two directions simultaneously**:

- **Clockwise sweep:** azimuths θ = (loc₂ + 1) to (loc₂ + 1 + 180)
- **Anti-clockwise sweep:** azimuths θ = (loc₁ − 1) to (loc₁ − 1 − 180)

Within each azimuth, correction proceeds in both **azimuthal** and **radial** directions:

- *Azimuthal step:* For each gate, a reference velocity v_r is computed from up to 3 previously-corrected neighbouring radials. If the neighbouring shears are consistent (`< α × v_n`), the correction integer is applied.
- *Radial step:* An initial gate (anchor) is identified where 3 consecutive azimuthal neighbours are consistent. Corrections then propagate outward and inward from that anchor gate in 20-km segments.

A `gate_flag` array tracks which gates have been successfully corrected (flag = 1).

### Phase 2 — Second-Pass Dealiasing

A second sweep is performed over all gates still flagged as uncorrected (flag = 0). The search window is extended to *M* = 3 additional radials, retrying with a broader azimuthal neighbourhood.

### Output

Corrected velocity fields `corr_vel` of shape `(num_azimuth, num_bins)` are saved as space-separated text files per elevation sweep.

---

## 3. Scan Mode Detection

The number of elevation sweeps (`num_elev`) determines the processing mode:

| Condition | Mode | Description |
|---|---|---|
| `num_elev ≤ 3` | `long_range` | Volume scan optimised for maximum range |
| `num_elev == 10` | `short_range` | Dense volumetric scan for convective analysis |

Velocity dealiasing is currently only applied to `short_range` scans (where high-resolution Doppler data are most critical).
