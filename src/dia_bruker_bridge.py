#!/usr/bin/env python3

from __future__ import annotations

import argparse
import contextlib
import json
import logging
import sys
import warnings


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract Bruker TIMS slice summaries through alphaTims."
    )
    parser.add_argument("--bruker", required=True, help="Path to a Bruker .d folder")
    parser.add_argument("--bruker-so", required=True, help="Path to timsdata.so")
    parser.add_argument("--mz", required=True, type=float, help="Target m/z center")
    parser.add_argument("--mz-ppm", required=True, type=float, help="m/z tolerance in ppm")
    parser.add_argument("--mz-da", default=None, type=float, help="absolute m/z tolerance")
    parser.add_argument("--mz-bins", default=160, type=int, help="m/z profile bin count")
    parser.add_argument("--rt-min", default=None, type=float, help="RT lower bound in minutes")
    parser.add_argument("--rt-max", default=None, type=float, help="RT upper bound in minutes")
    parser.add_argument("--im-min", default=None, type=float, help="mobility lower bound")
    parser.add_argument("--im-max", default=None, type=float, help="mobility upper bound")
    parser.add_argument("--quad-min", default=None, type=float, help="quad lower bound")
    parser.add_argument("--quad-max", default=None, type=float, help="quad upper bound")
    return parser.parse_args()


def require_bounds(lo: float | None, hi: float | None, label: str) -> None:
    if (lo is None) != (hi is None):
        raise SystemExit(f"{label} bounds must specify both min and max")
    if lo is not None and hi is not None and lo >= hi:
        raise SystemExit(f"{label} min must be smaller than max")


def mz_bounds(args: argparse.Namespace) -> tuple[float, float]:
    if args.mz_da is not None:
        delta = args.mz_da
    else:
        delta = args.mz * args.mz_ppm / 1_000_000.0
    return args.mz - delta, args.mz + delta


def patch_bruker_runtime(bruker_so: str) -> None:
    import alphatims.bruker as bruker

    original_open = bruker.open_bruker_d_folder
    bruker.BRUKER_DLL_FILE_NAME = bruker_so

    @contextlib.contextmanager
    def open_with_override(bruker_d_folder_name: str, bruker_dll_file_name=bruker_so):
        with original_open(
            bruker_d_folder_name,
            bruker_dll_file_name=bruker_so,
        ) as state:
            yield state

    bruker.open_bruker_d_folder = open_with_override


def float_slice(lo: float | None, hi: float | None):
    if lo is None:
        return slice(None)
    return slice(lo, hi)


def seconds_slice_from_minutes(lo: float | None, hi: float | None):
    if lo is None:
        return slice(None)
    return slice(lo * 60.0, hi * 60.0)


def main() -> int:
    args = parse_args()
    require_bounds(args.rt_min, args.rt_max, "RT")
    require_bounds(args.im_min, args.im_max, "mobility")
    require_bounds(args.quad_min, args.quad_max, "quadrupole")
    if args.mz_bins <= 0:
        raise SystemExit("--mz-bins must be at least 1")

    warnings.simplefilter("ignore", FutureWarning)

    try:
        import numpy as np
        import alphatims.utils
        from alphatims.bruker import TimsTOF
    except Exception as exc:
        raise SystemExit(
            "The Bruker backend requires alphaTims with numpy in the selected Python environment"
        ) from exc

    patch_bruker_runtime(args.bruker_so)
    alphatims.utils.set_progress_callback(None)
    logging.getLogger().handlers.clear()
    logging.getLogger().setLevel(logging.CRITICAL)

    tims = TimsTOF(
        args.bruker,
        slice_as_dataframe=True,
        use_calibrated_mz_values_as_default=2,
        use_hdf_if_available=True,
        mmap_detector_events=True,
    )
    mz_min, mz_max = mz_bounds(args)

    selection = {
        "rt_values": seconds_slice_from_minutes(args.rt_min, args.rt_max),
        "mobility_values": float_slice(args.im_min, args.im_max),
        "quad_mz_values": float_slice(args.quad_min, args.quad_max),
        "mz_values": slice(mz_min, mz_max),
    }
    df = tims[selection, "df"]
    intensity_col = (
        "corrected_intensity_values"
        if "corrected_intensity_values" in df.columns
        else "intensity_values"
    )

    frames = tims.frames
    if "MsMsType" in frames.columns:
        frames = frames[frames["MsMsType"] != 0]
    if args.rt_min is not None:
        frames = frames[
            (frames["Time"] >= args.rt_min * 60.0)
            & (frames["Time"] <= args.rt_max * 60.0)
        ]
    valid_frame_indices = set(int(value) for value in frames.index.tolist())
    frames_considered = int(len(valid_frame_indices))
    if not df.empty:
        df = df[df["frame_indices"].isin(valid_frame_indices)].reset_index(drop=True)

    if df.empty:
        result = {
            "acquisition_mode": tims.acquisition_mode,
            "intensity_column": intensity_col,
            "mz_min": float(mz_min),
            "mz_max": float(mz_max),
            "frames_considered": frames_considered,
            "frames_with_signal": 0,
            "matched_events": 0,
            "rt_profile": [],
            "mz_profile": [],
            "im_profile": [],
        }
        json.dump(result, fp=sys.stdout)
        sys.stdout.write("\n")
        return 0

    rt_profile = (
        df.groupby(["frame_indices", "rt_values_min"], as_index=False)
        .agg(
            summed_intensity=(intensity_col, "sum"),
            matched_events=(intensity_col, "size"),
        )
        .sort_values(["frame_indices", "rt_values_min"])
        .reset_index(drop=True)
    )
    im_profile = (
        df.groupby("mobility_values", as_index=False)
        .agg(summed_intensity=(intensity_col, "sum"))
        .sort_values("mobility_values")
        .reset_index(drop=True)
    )

    edges = np.linspace(mz_min, mz_max, args.mz_bins + 1)
    weights = df[intensity_col].to_numpy(dtype=float)
    mz_hist, _ = np.histogram(df["mz_values"].to_numpy(dtype=float), bins=edges, weights=weights)
    mz_centers = (edges[:-1] + edges[1:]) / 2.0

    result = {
        "acquisition_mode": tims.acquisition_mode,
        "intensity_column": intensity_col,
        "mz_min": float(mz_min),
        "mz_max": float(mz_max),
        "frames_considered": frames_considered,
        "frames_with_signal": int(df["frame_indices"].nunique()),
        "matched_events": int(len(df)),
        "rt_profile": [
            {
                "frame_index": int(row.frame_indices),
                "rt_minutes": float(row.rt_values_min),
                "summed_intensity": float(row.summed_intensity),
                "matched_events": int(row.matched_events),
            }
            for row in rt_profile.itertuples(index=False)
        ],
        "mz_profile": [
            {
                "mz_center": float(center),
                "summed_intensity": float(intensity),
            }
            for center, intensity in zip(mz_centers.tolist(), mz_hist.tolist())
        ],
        "im_profile": [
            {
                "mobility": float(row.mobility_values),
                "summed_intensity": float(row.summed_intensity),
            }
            for row in im_profile.itertuples(index=False)
        ],
    }
    json.dump(result, fp=sys.stdout)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
