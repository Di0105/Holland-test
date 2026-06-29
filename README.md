# Holland-test: Parametric TC Wind/Pressure Forcing for AtmoSurge DL Storm Surge Model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.6+](https://img.shields.io/badge/python-3.6%2B-blue.svg)](https://www.python.org/)

**Holland parametric tropical cyclone wind field reconstruction → AtmoSurge deep-learning storm surge hindcast for Hong Kong.**

---

## Overview

This repository implements the **STORM-lite Holland (1980) parametric wind model** — upgraded with Coriolis force, translation asymmetry, and Willoughby-Rahn RMW estimation — to generate station-level wind speed, wind direction, and local pressure forcing for **6 Hong Kong tide-gauge station pairs** during **8 historical tropical cyclones** (2001–2023).

The output is designed as meteorological forcing input for the **AtmoSurge** deep-learning storm surge prediction model.

### Scientific Context

- **8 TCs**: Utor (2001), Hagupit (2008), Vicente (2012), Usagi (2013), Hato (2017), Mangkhut (2018), Kompasu (2021), Saola (2023)
- **6 Hong Kong stations**: Kai Tak (SE), Cheung Chau (CCH), Lau Fau Shan (LFS), Sai Kung (SKG), Tai Mei Tuk (PLC), Waglan Island (WGL)
- **Comparison**: Holland-driven AtmoSurge vs observed-meteorology-driven AtmoSurge vs SLOSH historical report error
- **Key finding**: Raw Holland forcing systematically underestimates station wind at $r \gg RMW$; calibration tools included

---

## Repository Structure

```
Holland-test/
├── README.md                          ← English README (this file)
├── README_zh.md                       ← Chinese README (中文说明)
├── .gitignore
├── scripts/
│   ├── holland_atmosurge_pipeline_calculation.py   ← Main pipeline
│   └── calibrate_holland_wind.py                   ← Wind calibration tool
├── outputs/
│   ├── README.md                      ← Output file descriptions
│   └── wind_unit_audit.csv            ← MaxWind unit verification
├── docs/
│   └── formula_explanation.md         ← Full physics equations
└── data/
    └── README.md                      ← Input data specification
```

---

## Quick Start

### Prerequisites

- Python 3.6+
- No external dependencies (pure stdlib)

### Run the Pipeline

```bash
python scripts/holland_atmosurge_pipeline_calculation.py \
    --track /path/to/AtmoSurge_test1_track.csv \
    --judy-table /path/to/JUDYTABLE3.2.csv \
    --tc-nontc /path/to/TC_NonTC_Comparison_Table.csv \
    --output-dir ./outputs
```

### Print Formula Explanation

```bash
python scripts/holland_atmosurge_pipeline_calculation.py --explain-formula
```

### Calibrate Holland Wind Against Observations

```bash
python scripts/calibrate_holland_wind.py
```

---

## Key Design Principles

### 1. Rmax ≠ R34 (Strict Separation)

| Field | Meaning | Use in Holland |
|-------|---------|---------------|
| `Rmax_average_km` | Radius of Maximum Wind (RMW) | **Holland RMW input** |
| `R34_*` | Gale-force wind radius (烈風圈) | **QC only, never used as RMW** |

### 2. MaxWind Unit: kt, NOT km/h

The track CSV `MaxWind` is in **nautical miles per hour (knots)**. Conversion: `Vmax_ms = MaxWind × 0.514444`. Verified by unit audit (`outputs/wind_unit_audit.csv`).

### 3. AtmoSurge Pressure Input = Station Local Pressure P(r), NOT Pc

$$P(r) = P_c + \Delta P \cdot \exp(-(RMW/r)^B)$$

### 4. AtmoSurge Wind Input = Station Wind, NOT TC Center Vmax

Station wind = symmetric Holland gradient wind + translation background flow, reduced to 10 m.

---

## RMW Gap-Fill Strategy

1. **Observed**: `Rmax_average_km` from track CSV
2. **Interpolated**: Linear interpolation within the same TC event (no cross-TC interpolation)
3. **Formula**: Willoughby & Rahn (2004) — $R_{max} = 51.6 \cdot \exp(-0.0223 \cdot V_{max} + 0.0281 \cdot |lat|)$

RMW source is recorded in `holland_inputs.csv` column `rmw_source`.

---

## Known Limitations & Improvement Roadmap

| Issue | Status | Solution |
|-------|--------|----------|
| Wind underestimation at $r \gg RMW$ | 🔴 Active | Unclamp B, use R34 constraint, distance correction, ERA5 blend |
| Fixed Penv = 1013.25 hPa | 🟡 Planned | ERA5 outer-band MSLP per event |
| B clamp [0.8, 2.5] too restrictive | 🔴 Active | Allow B ∈ [0.3, 3.0] with R34 validation |
| No time-series adapter for AtmoSurge | 🟡 Planned | 10-min to 1-hour forcing sequences |
| Single Holland profile limitation | 🟡 Planned | Dual-B or GAHM approach |

---

## Output Files

| File | Rows × Columns | Description |
|------|---------------|-------------|
| `holland_inputs.csv` | 358 × 35 | Per-timestep B, RMW, Pc, Vmax, dP, Coriolis f, translation |
| `holland_station_forcing.csv` | 2148 × 22 | Station-level wind speed/direction, pressure P(r), distance to center |
| `event_pattern_coverage.csv` | 8 × 12 | RMW coverage fraction, Willoughby-Rahn usage per TC |
| `largest_surge_reference_comparison.csv` | — | Max observed surge vs max Holland wind per station |

---

## Related Repositories

| Repository | Purpose |
|-----------|---------|
| `Di0105/AtmoSurge_8TC_consolidated_package` | Final consolidated CSV + SLOSH + figures |
| `Di0105/CiteAgent-Copilot` | AI citation assistant for VS Code Copilot |

---

## References

- Jelesnianski, C.P. et al. (1992). *SLOSH: Sea, Lake, and Overland Surges from Hurricanes*. NOAA TR NWS 48.
- Holland, G.J. (1980). An analytic model of the wind and pressure profiles in hurricanes. *MWR*, 108(8).
- Holland, G.J. et al. (2010). A revised model for radial profiles of hurricane winds. *MWR*, 138(12).
- Willoughby, H.E. & Rahn, M.E. (2004). Parametric representation of the primary hurricane vortex. *MWR*, 132(12).
- Shashank, V.G. et al. (2022). Improvements in wind field hindcast for storm surge predictions. *Applied Ocean Research*.
- Liu, F. & Sasaki, J. (2019). Hybrid methods combining reanalysis and parametric typhoon model. *Scientific Reports*, 9.

---

## License

MIT License — see LICENSE file.

## Author

Judy Zhu (Di0105) — Hong Kong TC storm surge hindcast project, 2025–2026.