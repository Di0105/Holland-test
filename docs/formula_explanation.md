# Formula: STORM-lite Holland with Coriolis + Translation

> Version: `storm_lite_coriolis_translation_wr_rmw_v1`

## Constants

| Symbol | Value | Unit | Description |
|--------|-------|------|-------------|
| $P_{env}$ | 1013.25 | hPa | Environmental pressure (⚠️ first-pass; use ERA5 for production) |
| $\rho$ | 1.15 | kg/m³ | Air density |
| $\omega$ | $7.292\times10^{-5}$ | rad/s | Earth angular velocity |
| $C_{sfc}$ | 0.85 | — | Surface wind reduction factor (gradient → 10 m) |
| $\alpha$ | 0.55 | — | Background flow fraction of translation speed |
| $B_{min}$, $B_{max}$ | 0.8, 2.5 | — | Holland-B clamping range |
| $r_{min}$ | 1000 | m | Minimum radius floor |

## Input from Track CSV

```
Pc   = Min_Pressure                                [hPa]
Vmax = MaxWind [kt] × 0.514444                     [kt → m/s]
RMW  = Rmax_average_km × 1000                       [km → m, preferred]
     = Willoughby-Rahn if missing: 51.6 × exp(-0.0223·Vmax + 0.0281·|lat|)
R34  = R34_* fields                                 [QC only, NOT RMW]
Penv = 1013.25 (first-pass) or ERA5 outer-band MSLP [hPa]
dP   = Penv − Pc                                    [hPa]
```

## Coriolis & Translation

$$f = |2\omega\sin(\text{lat})|$$

$$V_{trans} = \frac{\text{Haversine distance between consecutive track points}}{\Delta t}$$

## Holland B (STORM-lite)

$$V_{surf} = V_{max} / C_{sfc}$$
$$V_{sym,max} = V_{surf} - \alpha \cdot V_{trans}$$
$$vv = \left(V_{sym,max} + \frac{f \cdot RMW}{2}\right)^2 - \frac{f^2 \cdot RMW^2}{4}$$
$$B_{raw} = \frac{vv \cdot e \cdot \rho}{\Delta P}$$
$$B = \text{clamp}(B_{raw}, 0.8, 2.5)$$

## Pressure Profile P(r)

$$P(r) = P_c + \Delta P \cdot \exp\left(-\left(\frac{RMW}{r}\right)^B\right)$$

## Gradient Wind V_g(r) with Coriolis

$$V_g(r) = \sqrt{\left(\frac{RMW}{r}\right)^B \cdot B \cdot \frac{\Delta P}{\rho} \cdot \exp\left(-\left(\frac{RMW}{r}\right)^B\right) + \frac{r^2 f^2}{4}} - \frac{f \cdot r}{2}$$

## 10 m Wind (Surface Wind Reduction)

$$V_{10,sym} = V_g(r) \times C_{sfc}$$

## Station-Level Total Wind

$$\vec{V}_{total} = \vec{V}_{sym} + \alpha \cdot \vec{V}_{trans}$$

Wind direction: cyclonic tangential (NH) + translation background flow → meteorological "wind from" direction.

## Known Limitation

When $r \gg RMW$ (station far from TC core), the single Holland profile decays too fast, systematically underestimating station wind speed. Recommended mitigations:
1. Unclamp B to allow $B < 0.8$ for outer wind fitting
2. Use R34 (gale radius) to constrain outer wind speed
3. Apply distance-dependent correction: $V_{corr} = V_{Holland} \times f(r/RMW)$
4. Blend with ERA5/NWP outer wind field (Liu & Sasaki, 2019)

## References

- Holland, G.J. (1980). *Monthly Weather Review*, 108(8), 1212–1218.
- Holland, G.J., Belanger, J.I., & Fritz, A. (2010). *Monthly Weather Review*, 138(12), 4393–4401.
- Willoughby, H.E. & Rahn, M.E. (2004). *Monthly Weather Review*, 132(12).
- Liu, F. & Sasaki, J. (2019). *Scientific Reports*, 9, 12209.
