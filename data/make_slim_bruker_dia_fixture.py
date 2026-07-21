#!/usr/bin/env python3
"""Build a small Bruker diaPASEF .d fixture from a frame range.

This keeps the selected frame blobs from analysis.tdf_bin, rewrites Frames.TimsId
to the new compact binary offsets, renumbers Frames.Id to 1..N, and updates
frame-referenced SQLite tables. It is intended for local fixture generation, not
general Bruker data conversion.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sqlite3
import struct
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-d", required=True, type=Path, help="Source Bruker .d folder")
    parser.add_argument("--out-d", required=True, type=Path, help="Output slim .d folder")
    parser.add_argument("--frame-start", required=True, type=int, help="First source frame Id")
    parser.add_argument("--frame-end", required=True, type=int, help="Last source frame Id")
    parser.add_argument(
        "--target-label",
        default="",
        help="Optional human-readable target label written to manifest JSON",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing output folder",
    )
    return parser.parse_args()


def read_blob(source_bin, offset: int) -> bytes:
    source_bin.seek(offset)
    header = source_bin.read(4)
    if len(header) != 4:
        raise ValueError(f"could not read byte count at source offset {offset}")
    byte_count = struct.unpack("<I", header)[0]
    source_bin.seek(offset)
    blob = source_bin.read(byte_count)
    if len(blob) != byte_count:
        raise ValueError(
            f"short read for source offset {offset}: expected {byte_count}, got {len(blob)}"
        )
    return blob


def table_columns(conn: sqlite3.Connection, table: str) -> set[str]:
    return {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}


def sqlite_tables(conn: sqlite3.Connection) -> list[str]:
    return [
        row[0]
        for row in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"
        )
    ]


def copy_data_file(source: Path, destination: Path) -> None:
    shutil.copy2(source, destination)
    destination.chmod(0o644)


def copy_optional_sidecars(source_d: Path, out_d: Path) -> None:
    for name in ("SampleInfo.xml", "HyStarMetadata.xml"):
        source = source_d / name
        if source.exists():
            copy_data_file(source, out_d / name)


def main() -> None:
    args = parse_args()
    source_d = args.source_d
    out_d = args.out_d
    source_tdf = source_d / "analysis.tdf"
    source_bin = source_d / "analysis.tdf_bin"
    if not source_tdf.exists() or not source_bin.exists():
        raise SystemExit(f"source .d must contain analysis.tdf and analysis.tdf_bin: {source_d}")
    if args.frame_end < args.frame_start:
        raise SystemExit("--frame-end must be >= --frame-start")
    if out_d.exists():
        if not args.overwrite:
            raise SystemExit(f"output exists; pass --overwrite to replace: {out_d}")
        shutil.rmtree(out_d)
    out_d.mkdir(parents=True)
    copy_optional_sidecars(source_d, out_d)
    copy_data_file(source_tdf, out_d / "analysis.tdf")

    conn = sqlite3.connect(out_d / "analysis.tdf")
    conn.row_factory = sqlite3.Row
    frames = conn.execute(
        """
        SELECT Id, TimsId, Time, MsMsType
        FROM Frames
        WHERE Id BETWEEN ? AND ?
        ORDER BY Id
        """,
        (args.frame_start, args.frame_end),
    ).fetchall()
    if not frames:
        raise SystemExit("no source frames selected")

    id_map = {int(row["Id"]): idx + 1 for idx, row in enumerate(frames)}
    offset_map: dict[int, int] = {}
    with source_bin.open("rb") as src, (out_d / "analysis.tdf_bin").open("wb") as dst:
        for row in frames:
            old_id = int(row["Id"])
            old_offset = int(row["TimsId"])
            new_offset = dst.tell()
            blob = read_blob(src, old_offset)
            dst.write(blob)
            offset_map[old_id] = new_offset

    conn.execute("CREATE TEMP TABLE slim_frame_map(old_id INTEGER PRIMARY KEY, new_id INTEGER, new_tims_id INTEGER)")
    conn.executemany(
        "INSERT INTO slim_frame_map(old_id, new_id, new_tims_id) VALUES (?, ?, ?)",
        [(old_id, id_map[old_id], offset_map[old_id]) for old_id in id_map],
    )

    conn.execute(
        "DELETE FROM Frames WHERE Id NOT IN (SELECT old_id FROM slim_frame_map)"
    )
    for old_id in sorted(id_map):
        conn.execute(
            "UPDATE Frames SET Id = ?, TimsId = ? WHERE Id = ?",
            (id_map[old_id], offset_map[old_id], old_id),
        )

    for table in sqlite_tables(conn):
        columns = table_columns(conn, table)
        if "Frame" in columns:
            conn.execute(
                f"DELETE FROM {table} WHERE Frame NOT IN (SELECT old_id FROM slim_frame_map)"
            )
            conn.execute(
                f"""
                UPDATE {table}
                SET Frame = (
                    SELECT new_id FROM slim_frame_map
                    WHERE old_id = {table}.Frame
                )
                WHERE Frame IN (SELECT old_id FROM slim_frame_map)
                """
            )

    if table_columns(conn, "Segments") >= {"FirstFrame", "LastFrame"}:
        conn.execute("UPDATE Segments SET FirstFrame = 1, LastFrame = ?", (len(frames),))

    conn.commit()
    conn.execute("VACUUM")
    conn.close()

    manifest = {
        "source_d": str(source_d),
        "out_d": str(out_d),
        "source_frame_start": args.frame_start,
        "source_frame_end": args.frame_end,
        "frame_count": len(frames),
        "source_rt_start_min": float(frames[0]["Time"]) / 60.0,
        "source_rt_end_min": float(frames[-1]["Time"]) / 60.0,
        "target_label": args.target_label,
        "renumbered_frames": "Frames.Id rewritten to 1..N; source frame ids retained only in this manifest",
    }
    (out_d.parent / f"{out_d.name}.manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
