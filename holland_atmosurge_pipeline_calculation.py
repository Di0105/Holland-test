"""Prepare Holland wind/pressure forcing from TC tracks for AtmoSurge.

This script is intentionally explicit about the RMW/R34 separation:
- Rmax_* and Rmax_average_km are treated as RMW inputs to Holland.
- R34_* are gale-radius fields only and are never used as RMW.

The current output is a supervised preprocessing product. It does not call the
private AtmoSurge DL model until that repository is available locally.
"""

import argparse
import csv
import math
from collections import defaultdict, namedtuple
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


DEFAULT_BASE_DIR = Path(r"C:\Users\26judyzhu\Desktop\Judy5")
DEFAULT_TRACK = DEFAULT_BASE_DIR / "AtmoSurge_test1_track.csv"
DEFAULT_JUDY_TABLE = DEFAULT_BASE_DIR / "JUDYTABLE3.2.csv"
DEFAULT_TC_NONTC = DEFAULT_BASE_DIR / "TC_NonTC_Comparison_Table.csv"
DEFAULT_OUTPUT_DIR = DEFAULT_BASE_DIR / "outputs"

PENV_HPA = 1013.25
RHO_AIR = 1.15
KNOT_TO_MPS = 0.514444
MIN_RADIUS_M = 1_000.0
HOLLAND_B_MIN = 0.8
HOLLAND_B_MAX = 2.5
EARTH_RADIUS_M = 6_371_000.0
OMEGA_EARTH = 7.292e-5
SURFACE_WIND_REDUCTION_FACTOR = 0.85
BACKGROUND_FLOW_ALPHA = 0.55
HOLLAND_FORMULA_VERSION = "storm_lite_coriolis_translation_wr_rmw_v1"

FORMULA_EXPLANATION = """
Holland calculation used in this script
--------------------------------------
CSV inputs used as Holland inputs:
    Pc   = Min_Pressure                       [hPa]
    Vmax = MaxWind * 0.514444                 [kt -> m/s]
    RMW  = Rmax_average_km * 1000             [km -> m, preferred]
    Penv = 1013.25                            [hPa, STORM-style Patm]
    dP   = Penv - Pc                          [hPa]

Important radius rule:
    Rmax/RMW is the radius of maximum wind and is used in the Holland formula.
    R34 is the gale-force radius / lie feng quan R and is kept only for QC.
    R34 is never used as RMW in this script.

RMW gap-fill rule:
    If Rmax_average_km is missing after within-TC linear interpolation, RMW is
    estimated by Willoughby and Rahn (2004):
    Rmax_km = 51.6 * exp(-0.0223 * Vmax_ms + 0.0281 * abs(latitude_deg))
    The output column rmw_source records whether RMW came from observed/interpolated
    Rmax_average_km or from the Willoughby-Rahn formula.

STORM-style Holland B estimate:
    f         = abs(2 * omega * sin(latitude))
    Vsurf     = Vmax / surface_wind_reduction_factor
    Vsym_max  = Vsurf - background_flow_alpha * translation_speed
    vv        = (Vsym_max + f * RMW / 2)^2 - f^2 * RMW^2 / 4
    B_raw     = vv * e * rho_air / dP_pa
    B     = clamp(B_raw, 0.8, 2.5)

Pressure profile at radius r from the TC center:
    P(r) = Pc + dP * exp(-(RMW / r)^B)

Gradient wind speed at radius r with Coriolis:
    Vg(r) = sqrt((RMW / r)^B * B * dP_pa / rho_air * exp(-(RMW / r)^B)
                 + r^2 * f^2 / 4) - f * r / 2
    V10_sym = Vg(r) * surface_wind_reduction_factor

Station forcing:
    r is computed from TC center Lat/Lon to each AtmoSurge station pair location.
    Wind direction combines Northern Hemisphere cyclonic tangential flow and the
    translation/background-flow vector from the track.
    holland_formula_version = storm_lite_coriolis_translation_wr_rmw_v1
""".strip()

# User-confirmed early rows to ignore. They lack usable RMW/R34 and should not be
# interpolated into the Mangkhut event start.
IGNORED_TRACK_TIMES = {
    ("MANGKHUT2018", datetime(2018, 9, 7, 6, 0)),
    ("MANGKHUT2018", datetime(2018, 9, 7, 12, 0)),
}

REQUIRED_HOLLAND_COLUMNS = [
    "Lat",
    "Lon",
    "Min_Pressure",
    "MaxWind",
]

INTERPOLATE_COLUMNS = [
    "Lat",
    "Lon",
    "Min_Pressure",
    "MaxWind",
    "Rmax_NE_km",
    "Rmax_NW_km",
    "Rmax_SW_km",
    "Rmax_SE_km",
    "Rmax_average_km",
    "Rmax_max_km",
    "R34_NE_km",
    "R34_SE_km",
    "R34_SW_km",
    "R34_NW_km",
    "R34_average_km",
    "R34_max_km",
]


Station = namedtuple("Station", "code name lat lon")


# Approximate station coordinates used only to make station-level forcing. Since we use the station pair, the coordinates error are within 5km, which is acceptable for the purpose of this analysis.
STATIONS = [
    Station("SE", "Kai Tak", 22.313, 114.213),
    Station("CCH", "Cheung Chau", 22.201, 114.026),
    Station("LFS", "Lau Fau Shan", 22.468, 113.984),
    Station("SKG", "Sai Kung", 22.384, 114.273),
    Station("PLC", "Tai Mei Tuk", 22.475, 114.237),
    Station("WGL", "Waglan Island", 22.182, 114.303),
]


def parse_datetime(value: str) -> datetime:
    value = value.strip()
    for fmt in ("%Y/%m/%d %H:%M", "%Y/%m/%d %H:%M:%S"):
        try:
            return datetime.strptime(value, fmt)
        except ValueError:
            pass
    raise ValueError(f"Unsupported datetime format: {value!r}")


def parse_float(value):
    if value is None:
        return None
    cleaned = value.strip()
    if cleaned in {"", "-", "---", "N/A", "nan", "NaN", "NULL"}:
        return None
    try:
        return float(cleaned)
    except ValueError:
        return None


def read_csv_dicts(path):
    with path.open("r", encoding="utf-8-sig", newline="") as csv_file:
        reader = csv.DictReader(csv_file, skipinitialspace=True)
        rows = []
        for row in reader:
            rows.append({key.strip(): (value or "").strip() for key, value in row.items()})
    return rows


def load_track(path):
    kept = []
    ignored = []

    for raw_row in read_csv_dicts(path):
        tc_name = str(raw_row["TC_name"]).strip()
        dt = parse_datetime(str(raw_row["Datetime"]))
        parsed = dict(raw_row)
        parsed["TC_name"] = tc_name
        parsed["Datetime"] = dt
        for column in INTERPOLATE_COLUMNS:
            parsed[column] = parse_float(raw_row.get(column))

        if (tc_name, dt) in IGNORED_TRACK_TIMES:
            parsed["ignore_reason"] = "user_confirmed_initial_missing_rmw_r34"
            ignored.append(parsed)
            continue
        kept.append(parsed)

    return kept, ignored


def interpolate_group_rows(rows, columns):
    rows.sort(key=lambda item: item["Datetime"])
    times = [datetime_to_seconds(row["Datetime"]) for row in rows]
    for column in columns:
        known_indexes = [index for index, row in enumerate(rows) if row.get(column) is not None]
        if len(known_indexes) < 2:
            continue

        for left_index, right_index in zip(known_indexes, known_indexes[1:]):
            left_time = times[left_index]
            right_time = times[right_index]
            left_value = float(rows[left_index][column])
            right_value = float(rows[right_index][column])
            if right_time == left_time:
                continue

            for fill_index in range(left_index + 1, right_index):
                if rows[fill_index].get(column) is not None:
                    continue
                weight = (times[fill_index] - left_time) / (right_time - left_time)
                rows[fill_index][column] = left_value + weight * (right_value - left_value)
                rows[fill_index].setdefault("interpolated_columns", set()).add(column)


def datetime_to_seconds(value):
    if not isinstance(value, datetime):
        raise TypeError(f"Expected datetime, got {type(value)!r}")
    return value.timestamp()


def clean_track(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[str(row["TC_name"])].append(row)

    cleaned = []
    for group_rows in grouped.values():
        interpolate_group_rows(group_rows, INTERPOLATE_COLUMNS)
        cleaned.extend(group_rows)

    cleaned.sort(key=lambda item: (str(item["TC_name"]), item["Datetime"]))
    return cleaned


def clamp(value, minimum, maximum):
    return max(minimum, min(maximum, value))


def print_formula_explanation():
    print(FORMULA_EXPLANATION)


def willoughby_rahn_rmax_km(vmax_ms, latitude_deg):
    """Estimate RMW/Rmax in km using Willoughby and Rahn (2004).

    The formula in the user's reference uses maximum wind speed and TC-center
    latitude. Here Vmax is in m/s because the track's MaxWind has already been
    converted from knots to m/s.
    """
    return 51.6 * math.exp(-0.0223 * vmax_ms + 0.0281 * abs(latitude_deg))


def coriolis_parameter(latitude_deg):
    return abs(2.0 * OMEGA_EARTH * math.sin(math.radians(latitude_deg)))


def translation_vector_from_rows(start_row, end_row):
    start_time = start_row["Datetime"]
    end_time = end_row["Datetime"]
    dt_seconds = datetime_to_seconds(end_time) - datetime_to_seconds(start_time)
    if dt_seconds <= 0:
        return 0.0, 0.0, 0.0

    start_lat = float(start_row["Lat"])
    start_lon = float(start_row["Lon"])
    end_lat = float(end_row["Lat"])
    end_lon = float(end_row["Lon"])
    distance_m = haversine_km(start_lat, start_lon, end_lat, end_lon) * 1000.0
    speed_ms = distance_m / dt_seconds
    movement_bearing = bearing_deg(start_lat, start_lon, end_lat, end_lon)
    u_ms, v_ms = wind_components(speed_ms, movement_bearing)
    return speed_ms, u_ms, v_ms


def add_translation_fields(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[str(row["TC_name"])].append(row)

    for group_rows in grouped.values():
        group_rows.sort(key=lambda item: item["Datetime"])
        for index, row in enumerate(group_rows):
            if index > 0:
                speed_ms, u_ms, v_ms = translation_vector_from_rows(group_rows[index - 1], row)
                source = "previous_track_step"
            elif len(group_rows) > 1:
                speed_ms, u_ms, v_ms = translation_vector_from_rows(row, group_rows[index + 1])
                source = "next_track_step_for_first_row"
            else:
                speed_ms, u_ms, v_ms = 0.0, 0.0, 0.0
                source = "single_point_event"

            row["translation_speed_ms"] = speed_ms
            row["translation_u_ms"] = u_ms
            row["translation_v_ms"] = v_ms
            row["background_u_ms"] = BACKGROUND_FLOW_ALPHA * u_ms
            row["background_v_ms"] = BACKGROUND_FLOW_ALPHA * v_ms
            row["translation_source"] = source


def max_consecutive_missing_rmw(rows):
    longest = 0
    current = 0
    for row in sorted(rows, key=lambda item: item["Datetime"]):
        if row.get("Rmax_average_km") is None:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def build_event_pattern_report(rows, ignored_rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[str(row["TC_name"])].append(row)

    report_rows = []
    for tc_name, group_rows in sorted(grouped.items()):
        group_rows.sort(key=lambda item: item["Datetime"])
        rmw_rows = [row for row in group_rows if row.get("Rmax_average_km") is not None]
        missing_rmw_rows = [row for row in group_rows if row.get("Rmax_average_km") is None]
        total_rows = len(group_rows)
        observed_coverage = (len(rmw_rows) / float(total_rows)) if total_rows else 0.0

        report_rows.append(
            {
                "tc_name": tc_name,
                "rows_after_user_ignore": total_rows,
                "rows_with_observed_or_interpolated_rmw": len(rmw_rows),
                "rows_needing_willoughby_rahn_rmw": len(missing_rmw_rows),
                "observed_rmw_coverage_fraction": observed_coverage,
                "pattern_ok_if_missing_rmw_rows_are_dropped": "yes" if len(missing_rmw_rows) == 0 else "no",
                "first_track_time": format_datetime(group_rows[0]["Datetime"]) if group_rows else "",
                "last_track_time": format_datetime(group_rows[-1]["Datetime"]) if group_rows else "",
                "first_observed_rmw_time": format_datetime(rmw_rows[0]["Datetime"]) if rmw_rows else "",
                "last_observed_rmw_time": format_datetime(rmw_rows[-1]["Datetime"]) if rmw_rows else "",
                "max_consecutive_missing_rmw_rows": max_consecutive_missing_rmw(group_rows),
            }
        )

    ignored_by_tc = defaultdict(int)
    for row in ignored_rows:
        ignored_by_tc[str(row["TC_name"])] += 1
    for row in report_rows:
        row["user_ignored_rows"] = ignored_by_tc.get(str(row["tc_name"]), 0)

    return report_rows


def build_holland_inputs(rows):
    inputs = []
    skipped = []

    for row in rows:
        missing = [column for column in REQUIRED_HOLLAND_COLUMNS if row.get(column) is None]
        if missing:
            skipped.append(
                {
                    "TC_name": row["TC_name"],
                    "Datetime": format_datetime(row["Datetime"]),
                    "skip_reason": f"missing_required_after_interpolation: {';'.join(missing)}",
                }
            )
            continue

        center_lat = float(row["Lat"])
        center_lon = float(row["Lon"])

        # Step 1: center pressure Pc.
        # The track CSV stores Min_Pressure in hPa. Holland pressure formulas use
        # pressure deficit dP, so we keep both hPa and Pa versions below.
        pc_hpa = float(row["Min_Pressure"])

        # Step 2: maximum wind Vmax.
        # The track CSV MaxWind is treated as knots, so convert it to m/s.
        vmax_ms = float(row["MaxWind"]) * KNOT_TO_MPS

        interpolated_columns = sorted(row.get("interpolated_columns", set()))

        # Step 3: RMW / Rmax.
        # Rmax_average_km is preferred. If it is still missing after within-TC
        # linear interpolation, estimate RMW with Willoughby and Rahn (2004).
        # R34 is NOT used here because it is the gale-force radius.
        if row.get("Rmax_average_km") is None:
            rmw_km = willoughby_rahn_rmax_km(vmax_ms, center_lat)
            rmw_source = "willoughby_rahn_2004_from_vmax_lat"
            rmw_was_estimated = True
        else:
            rmw_km = float(row["Rmax_average_km"])
            rmw_source = "Rmax_average_km_linear_interpolation" if "Rmax_average_km" in interpolated_columns else "Rmax_average_km_observed"
            rmw_was_estimated = False
        rmw_m = rmw_km * 1000.0

        # Step 4: pressure deficit dP = Penv - Pc.
        # Penv is currently a fixed first-pass environmental pressure assumption.
        dp_hpa = PENV_HPA - pc_hpa
        if dp_hpa <= 0 or rmw_m <= 0 or vmax_ms <= 0:
            skipped.append(
                {
                    "TC_name": row["TC_name"],
                    "Datetime": format_datetime(row["Datetime"]),
                    "skip_reason": "non_positive_holland_input",
                }
            )
            continue

        # Step 5: STORM-style Holland B with Coriolis and translation speed.
        # Following the same idea as STORM-return-periods, first remove part of
        # the background translation speed from the maximum wind before solving B.
        dp_pa = dp_hpa * 100.0
        translation_speed_ms = float(row.get("translation_speed_ms", 0.0))
        coriolis_f = coriolis_parameter(center_lat)
        vmax_surface_ms = vmax_ms / SURFACE_WIND_REDUCTION_FACTOR
        symmetric_vmax_for_b_ms = max(0.0, vmax_surface_ms - BACKGROUND_FLOW_ALPHA * translation_speed_ms)
        vv = (symmetric_vmax_for_b_ms + coriolis_f * rmw_m / 2.0) ** 2 - (coriolis_f * coriolis_f * rmw_m * rmw_m) / 4.0
        if vv <= 0:
            skipped.append(
                {
                    "TC_name": row["TC_name"],
                    "Datetime": format_datetime(row["Datetime"]),
                    "skip_reason": "non_positive_storm_lite_vv",
                }
            )
            continue
        holland_b_raw = vv * math.e * RHO_AIR / dp_pa
        holland_b = clamp(holland_b_raw, HOLLAND_B_MIN, HOLLAND_B_MAX)

        inputs.append(
            {
                "tc_name": row["TC_name"],
                "datetime": format_datetime(row["Datetime"]),
                "center_lat": center_lat,
                "center_lon": center_lon,
                "pc_hpa": pc_hpa,
                "pc_pa": pc_hpa * 100.0,
                "penv_hpa": PENV_HPA,
                "dp_hpa": dp_hpa,
                "dp_pa": dp_pa,
                "maxwind_kt": row["MaxWind"],
                "vmax_ms": vmax_ms,
                "rmw_km": rmw_km,
                "rmw_m": rmw_m,
                "rmw_source": rmw_source,
                "rmw_was_estimated_by_formula": rmw_was_estimated,
                "r34_average_km_qc_only": row.get("R34_average_km"),
                "r34_max_km_qc_only": row.get("R34_max_km"),
                "translation_speed_ms": translation_speed_ms,
                "translation_u_ms": row.get("translation_u_ms", 0.0),
                "translation_v_ms": row.get("translation_v_ms", 0.0),
                "background_u_ms": row.get("background_u_ms", 0.0),
                "background_v_ms": row.get("background_v_ms", 0.0),
                "translation_source": row.get("translation_source", ""),
                "coriolis_f": coriolis_f,
                "surface_wind_reduction_factor": SURFACE_WIND_REDUCTION_FACTOR,
                "background_flow_alpha": BACKGROUND_FLOW_ALPHA,
                "symmetric_vmax_for_b_ms": symmetric_vmax_for_b_ms,
                "holland_b_raw_from_vmax": holland_b_raw,
                "holland_b_used": holland_b,
                "holland_b_was_clipped": holland_b != holland_b_raw,
                "holland_formula_version": HOLLAND_FORMULA_VERSION,
                "interpolated_columns": ";".join(interpolated_columns),
                "source_file": row.get("source_file", ""),
                "source_row_number": row.get("source_row_number", ""),
            }
        )

    return inputs, skipped


def format_datetime(value):
    if isinstance(value, datetime):
        return value.strftime("%Y-%m-%d %H:%M:%S")
    return str(value)


def haversine_km(lat1, lon1, lat2, lon2):
    radius_km = 6371.0
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    d_phi = math.radians(lat2 - lat1)
    d_lambda = math.radians(lon2 - lon1)
    a = math.sin(d_phi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(d_lambda / 2.0) ** 2
    return 2.0 * radius_km * math.asin(math.sqrt(a))


def bearing_deg(lat1, lon1, lat2, lon2):
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    d_lambda = math.radians(lon2 - lon1)
    y = math.sin(d_lambda) * math.cos(phi2)
    x = math.cos(phi1) * math.sin(phi2) - math.sin(phi1) * math.cos(phi2) * math.cos(d_lambda)
    return (math.degrees(math.atan2(y, x)) + 360.0) % 360.0


def holland_pressure_hpa(pc_hpa, dp_hpa, rmw_m, radius_m, holland_b):
    """Holland pressure profile P(r).

    Pc is the central pressure. dP is the environmental pressure deficit.
    RMW controls where the wind reaches its maximum. The exponential term makes
    pressure approach Pc near the eye and Penv far away from the center.
    """
    radius_m = max(radius_m, MIN_RADIUS_M)
    exponent = -((rmw_m / radius_m) ** holland_b)
    return pc_hpa + dp_hpa * math.exp(exponent)


def holland_wind_ms(dp_pa, rmw_m, radius_m, holland_b, coriolis_f):
    """STORM-lite Holland gradient wind speed V(r), reduced to 10 m wind.

    This keeps the pipeline light but upgrades the original first-pass formula by
    adding the Coriolis term used in STORM-return-periods. The returned wind is a
    symmetric 10 m wind after applying the surface wind reduction factor.
    """
    radius_m = max(radius_m, MIN_RADIUS_M)
    scaled_radius = (rmw_m / radius_m) ** holland_b
    gradient_squared = scaled_radius * holland_b * dp_pa / RHO_AIR * math.exp(-scaled_radius)
    coriolis_squared = radius_m * radius_m * coriolis_f * coriolis_f / 4.0
    gradient_speed_ms = math.sqrt(max(0.0, gradient_squared + coriolis_squared)) - coriolis_f * radius_m / 2.0
    gradient_speed_ms = max(0.0, gradient_speed_ms)
    return gradient_speed_ms * SURFACE_WIND_REDUCTION_FACTOR


def wind_components(speed_ms, wind_to_direction_deg):
    radians = math.radians(wind_to_direction_deg)
    u_east = speed_ms * math.sin(radians)
    v_north = speed_ms * math.cos(radians)
    return u_east, v_north


def vector_to_direction_deg(u_east, v_north):
    if u_east == 0 and v_north == 0:
        return 0.0
    return (math.degrees(math.atan2(u_east, v_north)) + 360.0) % 360.0


def build_station_forcing(holland_inputs):
    forcing_rows = []
    for item in holland_inputs:
        center_lat = float(item["center_lat"])
        center_lon = float(item["center_lon"])
        pc_hpa = float(item["pc_hpa"])
        dp_hpa = float(item["dp_hpa"])
        dp_pa = float(item["dp_pa"])
        rmw_m = float(item["rmw_m"])
        holland_b = float(item["holland_b_used"])
        coriolis_f = float(item["coriolis_f"])
        background_u_ms = float(item["background_u_ms"])
        background_v_ms = float(item["background_v_ms"])

        for station in STATIONS:
            # Step 6: choose r for this station.
            # r is the great-circle distance from the TC center to the station
            # pair/proxy point. This r enters both P(r) and V(r).
            distance_km = haversine_km(center_lat, center_lon, station.lat, station.lon)
            radius_m = distance_km * 1000.0
            radial_bearing = bearing_deg(center_lat, center_lon, station.lat, station.lon)
            # Northern Hemisphere cyclonic tangential flow. This is the direction
            # the wind moves toward; meteorological wind direction is where it comes from.
            wind_to_deg = (radial_bearing - 90.0) % 360.0
            wind_from_deg = (wind_to_deg + 180.0) % 360.0

            # Step 7: compute station pressure and wind from the Holland profile.
            # The same RMW from Rmax_average_km is used for all stations at this
            # track time; only r changes from station to station.
            symmetric_speed_ms = holland_wind_ms(dp_pa, rmw_m, radius_m, holland_b, coriolis_f)
            symmetric_u_ms, symmetric_v_ms = wind_components(symmetric_speed_ms, wind_to_deg)
            u_east = symmetric_u_ms + background_u_ms
            v_north = symmetric_v_ms + background_v_ms
            speed_ms = math.sqrt(u_east * u_east + v_north * v_north)
            wind_to_total_deg = vector_to_direction_deg(u_east, v_north)
            wind_from_deg = (wind_to_total_deg + 180.0) % 360.0
            pressure_hpa = holland_pressure_hpa(pc_hpa, dp_hpa, rmw_m, radius_m, holland_b)

            forcing_rows.append(
                {
                    "tc_name": item["tc_name"],
                    "datetime": item["datetime"],
                    "station_code": station.code,
                    "station_name": station.name,
                    "station_lat": station.lat,
                    "station_lon": station.lon,
                    "distance_to_center_km": distance_km,
                    "bearing_center_to_station_deg": radial_bearing,
                    "pressure_hpa": pressure_hpa,
                    "wind_speed_ms": speed_ms,
                    "symmetric_holland_wind_speed_ms": symmetric_speed_ms,
                    "wind_from_direction_deg": wind_from_deg,
                    "u_east_ms": u_east,
                    "v_north_ms": v_north,
                    "symmetric_u_east_ms": symmetric_u_ms,
                    "symmetric_v_north_ms": symmetric_v_ms,
                    "translation_speed_ms": item["translation_speed_ms"],
                    "background_u_ms": item["background_u_ms"],
                    "background_v_ms": item["background_v_ms"],
                    "coriolis_f": item["coriolis_f"],
                    "surface_wind_reduction_factor": item["surface_wind_reduction_factor"],
                    "rmw_m": item["rmw_m"],
                    "rmw_source": item["rmw_source"],
                    "rmw_was_estimated_by_formula": item["rmw_was_estimated_by_formula"],
                    "r34_average_km_qc_only": item["r34_average_km_qc_only"],
                    "holland_b_used": holland_b,
                    "holland_formula_version": item["holland_formula_version"],
                }
            )
    return forcing_rows


def parse_station_from_abbreviation(value):
    if "(" in value and ")" in value:
        return value.split("(", 1)[1].split(")", 1)[0].strip()
    return value.strip()


def read_largest_surge_reference(judy_table, tc_nontc_table):
    by_station = {}

    if judy_table.exists():
        for row in read_csv_dicts(judy_table):
            station = parse_station_from_abbreviation(row.get("Station Pair (Abbreviation)", ""))
            maximum_surge = parse_float(row.get("Maximum Surge (m)"))
            wind_direction = row.get("Wind Direction", "").strip()
            pressure_level = row.get("Pressure Level", "").strip()
            # Overall wind-level rows have blank direction and pressure columns.
            if not station or maximum_surge is None or wind_direction or pressure_level:
                continue
            current = by_station.get(station)
            if current is None or maximum_surge > float(current["judytable_max_surge_m"]):
                by_station[station] = {
                    "station_name": station,
                    "judytable_max_surge_m": maximum_surge,
                    "judytable_wind_level": row.get("Wind Level", ""),
                    "judytable_wind_speed_range_ms": row.get("Wind Speed Range (m/s)", ""),
                }

    if tc_nontc_table.exists():
        for row in read_csv_dicts(tc_nontc_table):
            station = row.get("Station", "").strip()
            if not station:
                continue
            reference = by_station.setdefault(
                station,
                {
                    "station_name": station,
                    "judytable_max_surge_m": "",
                    "judytable_wind_level": "",
                    "judytable_wind_speed_range_ms": "",
                },
            )
            reference["tc_nontc_tc_max_surge_m"] = parse_float(row.get("TC Max Surge (m)"))
            reference["tc_nontc_tc_max_wind_level"] = row.get("TC Max WL", "")

    return sorted(by_station.values(), key=lambda item: str(item["station_name"]))


def merge_largest_forcing(largest_reference, forcing_rows):
    strongest_by_station = {}
    for row in forcing_rows:
        station = str(row["station_name"])
        current = strongest_by_station.get(station)
        if current is None or float(row["wind_speed_ms"]) > float(current["wind_speed_ms"]):
            strongest_by_station[station] = row

    merged = []
    for reference in largest_reference:
        station = str(reference["station_name"])
        strongest = strongest_by_station.get(station)
        merged_row = dict(reference)
        if strongest:
            merged_row.update(
                {
                    "max_holland_wind_tc": strongest["tc_name"],
                    "max_holland_wind_time": strongest["datetime"],
                    "max_holland_wind_ms": strongest["wind_speed_ms"],
                    "pressure_at_max_holland_wind_hpa": strongest["pressure_hpa"],
                    "distance_to_center_at_max_holland_wind_km": strongest["distance_to_center_km"],
                    "atmosurge_predicted_max_surge_m": "pending_private_repo_integration",
                }
            )
        else:
            merged_row["atmosurge_predicted_max_surge_m"] = "pending_private_repo_integration"
        merged.append(merged_row)
    return merged


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Holland forcing from TC track CSV.")
    parser.add_argument("--track", type=Path, default=DEFAULT_TRACK)
    parser.add_argument("--judy-table", type=Path, default=DEFAULT_JUDY_TABLE)
    parser.add_argument("--tc-nontc", type=Path, default=DEFAULT_TC_NONTC)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--explain-formula",
        action="store_true",
        help="Print the CSV-to-Holland formula mapping and exit.",
    )
    args = parser.parse_args()

    if args.explain_formula:
        print_formula_explanation()
        return

    raw_rows, ignored_rows = load_track(args.track)
    cleaned_rows = clean_track(raw_rows)
    add_translation_fields(cleaned_rows)
    event_pattern_report = build_event_pattern_report(cleaned_rows, ignored_rows)
    holland_inputs, skipped_rows = build_holland_inputs(cleaned_rows)
    station_forcing = build_station_forcing(holland_inputs)
    largest_reference = read_largest_surge_reference(args.judy_table, args.tc_nontc)
    largest_comparison = merge_largest_forcing(largest_reference, station_forcing)

    write_csv(args.output_dir / "ignored_track_rows.csv", ignored_rows)
    write_csv(args.output_dir / "event_pattern_coverage.csv", event_pattern_report)
    write_csv(args.output_dir / "skipped_holland_rows.csv", skipped_rows)
    write_csv(args.output_dir / "holland_inputs.csv", holland_inputs)
    write_csv(args.output_dir / "holland_station_forcing.csv", station_forcing)
    write_csv(args.output_dir / "largest_surge_reference_comparison.csv", largest_comparison)

    tc_names = sorted({str(row["TC_name"]) for row in cleaned_rows})
    print("TC events:", ", ".join(tc_names))
    print(f"Input track rows kept: {len(cleaned_rows)}")
    print(f"Ignored rows: {len(ignored_rows)}")
    print(f"Holland input rows: {len(holland_inputs)}")
    print(f"Skipped Holland rows: {len(skipped_rows)}")
    print(f"Station forcing rows: {len(station_forcing)}")
    formula_rmw_rows = [row for row in holland_inputs if row.get("rmw_was_estimated_by_formula")]
    print(f"RMW rows estimated by Willoughby-Rahn formula: {len(formula_rmw_rows)}")
    print(f"Outputs written to: {args.output_dir}")

    incomplete_without_formula = [row for row in event_pattern_report if row.get("pattern_ok_if_missing_rmw_rows_are_dropped") == "no"]
    if incomplete_without_formula:
        print("Dropping missing-RMW rows would break full event coverage for:", ", ".join(row["tc_name"] for row in incomplete_without_formula))

    if largest_comparison:
        largest = max(
            largest_comparison,
            key=lambda row: float(row.get("tc_nontc_tc_max_surge_m") or row.get("judytable_max_surge_m") or 0.0),
        )
        print(
            "Largest observed TC surge reference:",
            largest.get("station_name"),
            largest.get("tc_nontc_tc_max_surge_m") or largest.get("judytable_max_surge_m"),
            "m",
        )


if __name__ == "__main__":
    main()