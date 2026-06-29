# -*- coding: utf-8 -*-
"""Calibrate Holland parameters to better match observed station winds.

Problem:
  Holland-predicted station wind is much lower than TC_Group observed wind.
  Example: Mangkhut2018/PLC: Holland ~20.7 m/s vs observed ~41.6 m/s.

Adjustable parameters (ordered by impact on station wind):
  1. B (Holland shape parameter)  — controls radial decay rate
  2. RMW scaling factor            — scales the radius of maximum wind
  3. Surface wind reduction factor — direct multiplier on wind speed
  4. Penv (environmental pressure) — affects dP, indirectly affects wind

Validation: Leave-one-TC-out to avoid overfitting.
"""

import csv, json, math, sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
HOLLAND_INPUTS = Path(r"C:\Users\26judyzhu\Downloads\Holland-stage\outputs\holland_inputs.csv")
HOLLAND_FORCING = Path(r"C:\Users\26judyzhu\Downloads\Holland-stage\outputs\holland_station_forcing.csv")
OBS_PEAKS = Path(r"C:\Users\26judyzhu\Downloads\slosh_text_analysis_outputs\observed_8tc_event_peaks.csv")
OBS_MET_TEST = Path(r"C:\Users\26judyzhu\Downloads\slosh_text_analysis_outputs\atmosurge_observed_meteorology_8tc_event_peak_test.csv")
OUT_DIR = Path(r"C:\Users\26judyzhu\Downloads\Holland_calibration_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Constants (from your Holland pipeline) ─────────────────────────────────
RHO_AIR = 1.15
OMEGA_EARTH = 7.292e-5
SURFACE_WIND_REDUCTION_FACTOR = 0.85
BACKGROUND_FLOW_ALPHA = 0.55
PENV_HPA = 1013.25
MIN_RADIUS_M = 1000.0

# Parameter ranges for calibration
B_RANGE = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 2.0, 2.5]
RMW_SCALES = [0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 2.5, 3.0]
SURFACE_FACTORS = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00]


# ── Physics functions (replicated from your pipeline) ─────────────────────
def coriolis_f(lat):
    return abs(2.0 * OMEGA_EARTH * math.sin(math.radians(lat)))


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2.0 * r * math.asin(math.sqrt(a))


def bearing_deg(lat1, lon1, lat2, lon2):
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dlam = math.radians(lon2 - lon1)
    y = math.sin(dlam) * math.cos(phi2)
    x = math.cos(phi1) * math.sin(phi2) - math.sin(phi1) * math.cos(phi2) * math.cos(dlam)
    return (math.degrees(math.atan2(y, x)) + 360.0) % 360.0


def holland_wind_ms(dp_pa, rmw_m, radius_m, holland_b, coriolis_f_val, surface_factor):
    """Holland gradient wind at radius r, reduced to 10m."""
    radius_m = max(radius_m, MIN_RADIUS_M)
    scaled = (rmw_m / radius_m) ** holland_b
    inside = scaled * holland_b * dp_pa / RHO_AIR * math.exp(-scaled)
    coriolis_term = radius_m * radius_m * coriolis_f_val * coriolis_f_val / 4.0
    gradient_ms = math.sqrt(max(0.0, inside + coriolis_term)) - coriolis_f_val * radius_m / 2.0
    gradient_ms = max(0.0, gradient_ms)
    return gradient_ms * surface_factor


def compute_station_wind(center_lat, center_lon, station_lat, station_lon,
                         dp_pa, rmw_m, holland_b, coriolis_f_val,
                         surface_factor, bg_u, bg_v):
    """Compute total station wind = symmetric Holland + background flow."""
    dist_km = haversine_km(center_lat, center_lon, station_lat, station_lon)
    radius_m = dist_km * 1000.0
    sym_speed = holland_wind_ms(dp_pa, rmw_m, radius_m, holland_b, coriolis_f_val, surface_factor)

    # Wind direction: NH cyclonic tangential (bearing - 90°)
    radial_bearing = bearing_deg(center_lat, center_lon, station_lat, station_lon)
    wind_to_deg = (radial_bearing - 90.0) % 360.0
    sym_u = sym_speed * math.sin(math.radians(wind_to_deg))
    sym_v = sym_speed * math.cos(math.radians(wind_to_deg))

    total_u = sym_u + bg_u
    total_v = sym_v + bg_v
    total_speed = math.sqrt(total_u ** 2 + total_v ** 2)
    wind_from_deg = (math.degrees(math.atan2(total_u, total_v)) + 360.0) % 360.0
    return total_speed, wind_from_deg, sym_speed, dist_km, radius_m


# ── Data loading ──────────────────────────────────────────────────────────
def read_csv(path):
    with open(str(path), newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def parse_float(v):
    try:
        return float(str(v or "").strip())
    except (ValueError, TypeError):
        return None


def parse_dt(v):
    v = str(v or "").strip().replace("T", " ")
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y/%m/%d %H:%M"):
        try:
            return datetime.strptime(v, fmt)
        except ValueError:
            pass
    return None


STATIONS = [
    ("SE", 22.313, 114.213, "Kai Tak"),
    ("CCH", 22.201, 114.026, "Cheung Chau"),
    ("LFS", 22.468, 113.984, "Lau Fau Shan"),
    ("SKG", 22.384, 114.273, "Sai Kung"),
    ("PLC", 22.475, 114.237, "Tai Mei Tuk"),
    ("WGL", 22.182, 114.303, "Waglan Island"),
]
STN_LOOKUP = {code: (lat, lon, name) for code, lat, lon, name in STATIONS}


def load_matched_pairs():
    """Match Holland forcing rows with TC_Group observed wind from atmosurge test CSV."""
    forcing_rows = read_csv(HOLLAND_FORCING) if HOLLAND_FORCING.exists() else []
    holland_inp = read_csv(HOLLAND_INPUTS) if HOLLAND_INPUTS.exists() else []
    obs_met_rows = read_csv(OBS_MET_TEST) if OBS_MET_TEST.exists() else []

    # From observed meteorology test, get unique (tc, stn, obs_time, wind) entries
    # Use ramped_24h default branch rows (primary method)
    obs_targets = {}
    for row in obs_met_rows:
        if row.get("sequence_mode") != "ramped_24h":
            continue
        if row.get("branch") != "default":
            continue
        tc = row.get("tc_name", "")
        stn = row.get("station_code", "")
        wind = parse_float(row.get("input_wind_speed_max_ms"))
        surge = parse_float(row.get("observed_peak_signed_surge_m"))
        peak_time = parse_dt(row.get("observed_peak_time"))
        window_end = parse_dt(row.get("window_end"))
        if not tc or not stn or wind is None or wind <= 0:
            continue
        ref_time = peak_time or window_end
        if not ref_time:
            continue
        key = (tc, stn)
        if key not in obs_targets:
            obs_targets[key] = {
                "tc_name": tc, "station_code": stn,
                "observed_wind_ms": wind, "observed_surge_m": surge,
                "peak_time": ref_time,
            }

    # Build lookup: (tc_name, station_code) -> forcing rows sorted by time
    forcing_by_key = defaultdict(list)
    for row in forcing_rows:
        key = (row.get("tc_name", ""), row.get("station_code", ""))
        dt = parse_dt(row.get("datetime"))
        if dt:
            forcing_by_key[key].append((dt, row))
    for key in forcing_by_key:
        forcing_by_key[key].sort(key=lambda x: x[0])

    # Build holland_input lookup by (tc, dt)
    holland_by_key = {}
    for row in holland_inp:
        tc = row.get("tc_name", "")
        dt = parse_dt(row.get("datetime"))
        if dt:
            holland_by_key[(tc, dt)] = row

    matched = []
    for (tc, stn), obs in obs_targets.items():
        peak_dt = obs["peak_time"]
        stn_info = STN_LOOKUP.get(stn)
        if not stn_info:
            continue

        # Find nearest Holland forcing row in time
        flist = forcing_by_key.get((tc, stn), [])
        if not flist:
            continue

        best_fr = None
        best_diff = float("inf")
        for fr_dt, fr in flist:
            diff = abs((fr_dt - peak_dt).total_seconds())
            if diff < best_diff:
                best_diff = diff
                best_fr = fr

        if not best_fr:
            continue

        fr_dt = parse_dt(best_fr.get("datetime"))
        hi = holland_by_key.get((tc, fr_dt)) if fr_dt else None
        if not hi or not hi.get("dp_pa") or not hi.get("rmw_m"):
            continue

        holland_ws = parse_float(best_fr.get("wind_speed_ms"))
        if holland_ws is None:
            continue

        matched.append({
            "tc_name": tc,
            "station_code": stn,
            "station_name": stn_info[2],
            "station_lat": stn_info[0],
            "station_lon": stn_info[1],
            "peak_time": peak_dt.strftime("%Y-%m-%d %H:%M:%S"),
            "forcing_time": best_fr.get("datetime", ""),
            "dt_diff_seconds": best_diff,
            "observed_surge_m": obs["observed_surge_m"],
            "observed_wind_ms": obs["observed_wind_ms"],
            "holland_wind_ms": holland_ws,
            "holland_pressure_hpa": parse_float(best_fr.get("pressure_hpa")),
            "pc_hpa": parse_float(hi.get("pc_hpa")),
            "dp_pa": parse_float(hi.get("dp_pa")),
            "rmw_m": parse_float(hi.get("rmw_m")),
            "holland_b": parse_float(hi.get("holland_b_used")),
            "holland_b_raw": parse_float(hi.get("holland_b_raw_from_vmax")),
            "b_was_clipped": hi.get("holland_b_was_clipped") == "True",
            "center_lat": parse_float(hi.get("center_lat")),
            "center_lon": parse_float(hi.get("center_lon")),
            "vmax_ms": parse_float(hi.get("vmax_ms")),
            "translation_speed_ms": parse_float(hi.get("translation_speed_ms")),
            "bg_u_ms": parse_float(best_fr.get("background_u_ms")) or 0.0,
            "bg_v_ms": parse_float(best_fr.get("background_v_ms")) or 0.0,
            "distance_km": parse_float(best_fr.get("distance_to_center_km")),
            "rmw_source": hi.get("rmw_source", ""),
        })

    return matched

    return matched


# ── Sensitivity analysis ──────────────────────────────────────────────────
def run_sensitivity(matched_pairs):
    """For each matched event-station, test parameter sensitivity."""
    results = []
    for pair in matched_pairs:
        dp_pa = pair["dp_pa"]
        rmw_m = pair["rmw_m"]
        b_orig = pair["holland_b"]
        c_lat = pair["center_lat"]
        c_lon = pair["center_lon"]
        s_lat = pair["station_lat"]
        s_lon = pair["station_lon"]
        bg_u = pair["bg_u_ms"]
        bg_v = pair["bg_v_ms"]
        f = coriolis_f(c_lat)

        obs_wind = pair["observed_wind_ms"]
        holl_wind = pair["holland_wind_ms"]
        if obs_wind is None or holl_wind is None:
            continue

        # Test varying B
        for b_test in B_RANGE:
            ws, wd, sym, dist, _ = compute_station_wind(
                c_lat, c_lon, s_lat, s_lon,
                dp_pa, rmw_m, b_test, f, SURFACE_WIND_REDUCTION_FACTOR,
                bg_u, bg_v)
            results.append({
                **pair,
                "param_type": "B",
                "param_value": b_test,
                "predicted_wind_ms": round(ws, 2),
                "wind_error_ms": round(ws - obs_wind, 2) if obs_wind else None,
                "abs_error_ms": round(abs(ws - obs_wind), 2) if obs_wind else None,
            })

        # Test varying RMW scale
        for rmw_scale in RMW_SCALES:
            ws, wd, sym, dist, _ = compute_station_wind(
                c_lat, c_lon, s_lat, s_lon,
                dp_pa, rmw_m * rmw_scale, b_orig, f, SURFACE_WIND_REDUCTION_FACTOR,
                bg_u, bg_v)
            results.append({
                **pair,
                "param_type": "RMW_scale",
                "param_value": rmw_scale,
                "predicted_wind_ms": round(ws, 2),
                "wind_error_ms": round(ws - obs_wind, 2) if obs_wind else None,
                "abs_error_ms": round(abs(ws - obs_wind), 2) if obs_wind else None,
            })

        # Test varying surface reduction factor
        for sf in SURFACE_FACTORS:
            ws, wd, sym, dist, _ = compute_station_wind(
                c_lat, c_lon, s_lat, s_lon,
                dp_pa, rmw_m, b_orig, f, sf, bg_u, bg_v)
            results.append({
                **pair,
                "param_type": "surface_factor",
                "param_value": sf,
                "predicted_wind_ms": round(ws, 2),
                "wind_error_ms": round(ws - obs_wind, 2) if obs_wind else None,
                "abs_error_ms": round(abs(ws - obs_wind), 2) if obs_wind else None,
            })

    return results


# ── Best-fit per event-station ────────────────────────────────────────────
def find_best_fit_per_pair(matched_pairs):
    """For each event-station, find the B and RMW that minimize wind error."""
    best_fits = []
    for pair in matched_pairs:
        dp_pa = pair["dp_pa"]
        rmw_m_orig = pair["rmw_m"]
        c_lat = pair["center_lat"]
        c_lon = pair["center_lon"]
        s_lat = pair["station_lat"]
        s_lon = pair["station_lon"]
        bg_u = pair["bg_u_ms"]
        bg_v = pair["bg_v_ms"]
        f = coriolis_f(c_lat)
        obs_wind = pair["observed_wind_ms"]
        if obs_wind is None or obs_wind <= 0:
            continue

        best_error = float("inf")
        best_config = None
        for b_test in [v / 10.0 for v in range(5, 35)]:  # 0.5 to 3.4, step 0.1
            for rmw_scale in [v / 10.0 for v in range(3, 51)]:  # 0.3 to 5.0
                for sf in [v / 100.0 for v in range(70, 101)]:  # 0.70 to 1.00
                    ws, _, _, _, _ = compute_station_wind(
                        c_lat, c_lon, s_lat, s_lon,
                        dp_pa, rmw_m_orig * rmw_scale, b_test, f, sf,
                        bg_u, bg_v)
                    error = abs(ws - obs_wind)
                    if error < best_error:
                        best_error = error
                        best_config = {
                            "B_opt": b_test,
                            "RMW_scale_opt": rmw_scale,
                            "RMW_opt_km": rmw_m_orig * rmw_scale / 1000.0,
                            "surface_factor_opt": sf,
                            "predicted_wind_ms": round(ws, 2),
                            "abs_error_ms": round(error, 2),
                        }
        if best_config:
            best_fits.append({
                **pair,
                **best_config,
                "original_wind_ms": pair["holland_wind_ms"],
                "original_b": pair["holland_b"],
                "original_rmw_km": pair["rmw_m"] / 1000.0,
            })
    return best_fits


# ── Leave-one-TC-out calibration ──────────────────────────────────────────
def leave_one_tc_out(matched_pairs):
    """Calibrate per-station parameters, validate with leave-one-TC-out."""
    all_tcs = sorted(set(p["tc_name"] for p in matched_pairs))
    results = []

    for held_out_tc in all_tcs:
        train = [p for p in matched_pairs if p["tc_name"] != held_out_tc]
        test = [p for p in matched_pairs if p["tc_name"] == held_out_tc]

        # Calibrate per-station average parameters from training data
        stn_params = defaultdict(list)
        for p in train:
            dp_pa = p["dp_pa"]
            rmw_m = p["rmw_m"]
            c_lat = p["center_lat"]
            c_lon = p["center_lon"]
            s_lat = p["station_lat"]
            s_lon = p["station_lon"]
            bg_u = p["bg_u_ms"]
            bg_v = p["bg_v_ms"]
            f = coriolis_f(c_lat)
            obs_wind = p["observed_wind_ms"]
            if obs_wind is None or obs_wind <= 0:
                continue

            best_error = float("inf")
            best_b, best_rmw_s, best_sf = None, None, None
            for b_t in [v / 10.0 for v in range(5, 31, 2)]:
                for rmw_s in [v / 10.0 for v in range(3, 41, 2)]:
                    for sf in [0.80, 0.85, 0.90, 0.95, 1.00]:
                        ws, _, _, _, _ = compute_station_wind(
                            c_lat, c_lon, s_lat, s_lon,
                            dp_pa, rmw_m * rmw_s, b_t, f, sf, bg_u, bg_v)
                        err = abs(ws - obs_wind)
                        if err < best_error:
                            best_error = err
                            best_b, best_rmw_s, best_sf = b_t, rmw_s, sf

            if best_b is not None:
                stn_params[p["station_code"]].append({
                    "B": best_b, "RMW_scale": best_rmw_s, "surface_factor": best_sf
                })

        # Average per station
        stn_avg = {}
        for stn, params in stn_params.items():
            stn_avg[stn] = {
                "B": sum(x["B"] for x in params) / len(params),
                "RMW_scale": sum(x["RMW_scale"] for x in params) / len(params),
                "surface_factor": sum(x["surface_factor"] for x in params) / len(params),
                "n_train_events": len(params),
            }

        # Validate on held-out TC
        for p in test:
            stn = p["station_code"]
            if stn not in stn_avg:
                continue
            cfg = stn_avg[stn]
            dp_pa = p["dp_pa"]
            rmw_m = p["rmw_m"]
            f = coriolis_f(p["center_lat"])
            ws_cal, wd, _, dist, _ = compute_station_wind(
                p["center_lat"], p["center_lon"],
                p["station_lat"], p["station_lon"],
                dp_pa, rmw_m * cfg["RMW_scale"], cfg["B"], f,
                cfg["surface_factor"],
                p["bg_u_ms"], p["bg_v_ms"])

            orig_ws = p["holland_wind_ms"] or 0
            obs_ws = p["observed_wind_ms"] or 0
            results.append({
                "held_out_tc": held_out_tc,
                "tc_name": p["tc_name"],
                "station_code": stn,
                "station_name": p["station_name"],
                "observed_wind_ms": obs_ws,
                "original_holland_wind_ms": orig_ws,
                "original_error_ms": orig_ws - obs_ws,
                "calibrated_wind_ms": round(ws_cal, 2),
                "calibrated_error_ms": round(ws_cal - obs_ws, 2),
                "calibrated_B": round(cfg["B"], 2),
                "calibrated_RMW_scale": round(cfg["RMW_scale"], 2),
                "calibrated_RMW_km": round(rmw_m * cfg["RMW_scale"] / 1000.0, 1),
                "calibrated_surface_factor": round(cfg["surface_factor"], 2),
                "n_train_events": cfg["n_train_events"],
                "distance_km": p["distance_km"],
                "original_B": p["holland_b"],
                "original_RMW_km": p["rmw_m"] / 1000.0,
                "original_surface_factor": SURFACE_WIND_REDUCTION_FACTOR,
            })

    return results, stn_avg


# ── Summary ───────────────────────────────────────────────────────────────
def summarize(best_fits, loocv_results, stn_avg):
    """Print and save calibration summary."""

    # Best-fit summary per station
    stn_best = defaultdict(list)
    for bf in best_fits:
        stn_best[bf["station_code"]].append(bf)

    print("=" * 70)
    print("PER-STATION BEST-FIT HOLLAND PARAMETERS")
    print("=" * 70)
    for stn, fits in sorted(stn_best.items()):
        bs = [f["B_opt"] for f in fits]
        rmws = [f["RMW_scale_opt"] for f in fits]
        sfs = [f["surface_factor_opt"] for f in fits]
        orig_errs = [abs((f["original_wind_ms"] or 0) - (f["observed_wind_ms"] or 0)) for f in fits]
        cal_errs = [f["abs_error_ms"] for f in fits]
        print(f"\n{stn} ({fits[0]['station_name']}): {len(fits)} events")
        print(f"  B_opt:       mean={sum(bs)/len(bs):.2f}  range=[{min(bs):.2f}, {max(bs):.2f}]")
        print(f"  RMW_scale:   mean={sum(rmws)/len(rmws):.2f}  range=[{min(rmws):.2f}, {max(rmws):.2f}]")
        print(f"  surface_fac: mean={sum(sfs)/len(sfs):.2f}  range=[{min(sfs):.2f}, {max(sfs):.2f}]")
        print(f"  MAE original: {sum(orig_errs)/len(orig_errs):.2f} m/s")
        print(f"  MAE calibrated: {sum(cal_errs)/len(cal_errs):.2f} m/s")

    # LOOCV summary
    if loocv_results:
        print("\n" + "=" * 70)
        print("LEAVE-ONE-TC-OUT VALIDATION")
        print("=" * 70)
        orig_mae = sum(abs(r["original_error_ms"]) for r in loocv_results) / len(loocv_results)
        cal_mae = sum(abs(r["calibrated_error_ms"]) for r in loocv_results) / len(loocv_results)
        print(f"\n  N = {len(loocv_results)}")
        print(f"  Original MAE: {orig_mae:.2f} m/s")
        print(f"  Calibrated MAE: {cal_mae:.2f} m/s")
        print(f"  Improvement: {(1 - cal_mae / max(orig_mae, 0.01)) * 100:.0f}%")

        # Per-station LOOCV
        stn_loocv = defaultdict(list)
        for r in loocv_results:
            stn_loocv[r["station_code"]].append(r)
        print("\n  Per-station:")
        for stn, rows in sorted(stn_loocv.items()):
            o_mae = sum(abs(r["original_error_ms"]) for r in rows) / len(rows)
            c_mae = sum(abs(r["calibrated_error_ms"]) for r in rows) / len(rows)
            print(f"    {stn}: orig MAE={o_mae:.2f}, cal MAE={c_mae:.2f}")

    # Recommended operational parameters
    print("\n" + "=" * 70)
    print("RECOMMENDED OPERATIONAL PARAMETERS (from LOOCV)")
    print("=" * 70)
    for stn, cfg in sorted(stn_avg.items()):
        print(f"  {stn}: B={cfg['B']:.2f}, RMW_scale={cfg['RMW_scale']:.2f}, "
              f"surface_factor={cfg['surface_factor']:.2f} (n_train={cfg['n_train_events']})")

    # Write CSVs
    if best_fits:
        fields = list(best_fits[0].keys())
        with open(OUT_DIR / "holland_best_fit_per_event.csv", "w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(best_fits)
    if loocv_results:
        fields = list(loocv_results[0].keys())
        with open(OUT_DIR / "holland_loocv_validation.csv", "w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(loocv_results)

    print(f"\nOutput written to: {OUT_DIR}")


# ── Main ──────────────────────────────────────────────────────────────────
def main():
    print("Loading matched pairs...")
    pairs = load_matched_pairs()
    print(f"Matched Holland-observed pairs: {len(pairs)}")

    if not pairs:
        print("ERROR: No matched pairs found. Check input files.")
        print(f"  HOLLAND_FORCING exists: {HOLLAND_FORCING.exists()}")
        print(f"  HOLLAND_INPUTS exists: {HOLLAND_INPUTS.exists()}")
        print(f"  OBS_PEAKS exists: {OBS_PEAKS.exists()}")
        print(f"  OBS_MET_TEST exists: {OBS_MET_TEST.exists()}")
        return

    # Show current errors
    print("\nCurrent Holland wind errors:")
    for p in sorted(pairs, key=lambda x: abs((x["observed_wind_ms"] or 0) - (x["holland_wind_ms"] or 0)), reverse=True)[:10]:
        err = (p["holland_wind_ms"] or 0) - (p["observed_wind_ms"] or 0)
        print(f"  {p['tc_name']}/{p['station_code']}: "
              f"Holland={p['holland_wind_ms']:.1f}, Obs={p['observed_wind_ms']:.1f}, "
              f"Error={err:.1f} m/s, dist={p['distance_km']:.0f} km, "
              f"B={p['holland_b']:.2f}, RMW={p['rmw_m']/1000:.0f} km")

    # Run sensitivity (light version for summary)
    print("\nRunning sensitivity analysis...")
    # sens = run_sensitivity(pairs)  # heavy, skip for now

    # Best fit per event-station
    print("Finding best-fit parameters per event-station...")
    best_fits = find_best_fit_per_pair(pairs)

    # LOOCV
    print("Running leave-one-TC-out calibration...")
    loocv_results, stn_avg = leave_one_tc_out(pairs)

    # Summarize
    summarize(best_fits, loocv_results, stn_avg)


if __name__ == "__main__":
    main()
