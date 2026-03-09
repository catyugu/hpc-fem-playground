#!/usr/bin/env python3
"""Compare MPFEM text results with COMSOL reference results."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ResultRow:
    x: float
    y: float
    z: float
    v: float
    t: float
    disp: float


def parse_result_file(file_path: Path) -> list[ResultRow]:
    rows: list[ResultRow] = []
    with file_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("%"):
                continue
            if line.startswith("x") or line.startswith("X"):
                continue

            parts = line.split()
            if len(parts) < 6:
                continue

            rows.append(
                ResultRow(
                    x=float(parts[0]),
                    y=float(parts[1]),
                    z=float(parts[2]),
                    v=float(parts[3]),
                    t=float(parts[4]),
                    disp=float(parts[5]),
                )
            )
    return rows


def compute_metrics(reference: list[float], current: list[float]) -> tuple[float, float, float, float]:
    if len(reference) != len(current):
        raise ValueError("Result vectors have different lengths")

    n = len(reference)
    if n == 0:
        raise ValueError("No rows available for comparison")

    sq_sum = 0.0
    ref_sq_sum = 0.0
    max_abs = 0.0
    max_rel = 0.0
    for i in range(n):
        diff = current[i] - reference[i]
        abs_diff = abs(diff)
        sq_sum += diff * diff
        ref_sq_sum += reference[i] * reference[i]
        if abs_diff > max_abs:
            max_abs = abs_diff

        denom = max(abs(reference[i]), 1e-16)
        rel = abs_diff / denom
        if rel > max_rel:
            max_rel = rel

    l2 = math.sqrt(sq_sum / n)
    if ref_sq_sum > 0.0:
        l2_relative = math.sqrt(sq_sum / ref_sq_sum)
    else:
        l2_relative = 0.0
    return l2, max_abs, max_rel, l2_relative


def coordinate_key(x: float, y: float, z: float) -> tuple[int, int, int]:
    scale = 10**12
    return (int(round(x * scale)), int(round(y * scale)), int(round(z * scale)))


def align_by_coordinates(reference_rows: list[ResultRow],
                         current_rows: list[ResultRow]) -> tuple[list[ResultRow], list[ResultRow]]:
    reference_map: dict[tuple[int, int, int], ResultRow] = {}
    for row in reference_rows:
        key = coordinate_key(row.x, row.y, row.z)
        if key in reference_map:
            raise ValueError(f"Duplicate coordinate in reference file: ({row.x}, {row.y}, {row.z})")
        reference_map[key] = row

    current_map: dict[tuple[int, int, int], ResultRow] = {}
    for row in current_rows:
        key = coordinate_key(row.x, row.y, row.z)
        if key in current_map:
            raise ValueError(f"Duplicate coordinate in current file: ({row.x}, {row.y}, {row.z})")
        current_map[key] = row

    if len(reference_map) != len(current_map):
        raise ValueError(
            f"Point count mismatch after coordinate mapping: reference={len(reference_map)} current={len(current_map)}"
        )

    missing_in_current = [key for key in reference_map if key not in current_map]
    if missing_in_current:
        raise ValueError(f"Current result is missing {len(missing_in_current)} reference points")

    missing_in_reference = [key for key in current_map if key not in reference_map]
    if missing_in_reference:
        raise ValueError(f"Reference result is missing {len(missing_in_reference)} current points")

    ordered_keys = sorted(reference_map.keys())
    aligned_reference = [reference_map[key] for key in ordered_keys]
    aligned_current = [current_map[key] for key in ordered_keys]
    return aligned_reference, aligned_current


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare MPFEM and COMSOL results")
    parser.add_argument("reference", type=Path, help="COMSOL reference result.txt")
    parser.add_argument("current", type=Path, help="MPFEM output result file")
    args = parser.parse_args()

    reference_rows = parse_result_file(args.reference)
    current_rows = parse_result_file(args.current)

    reference_rows, current_rows = align_by_coordinates(reference_rows, current_rows)

    ref_v = [row.v for row in reference_rows]
    cur_v = [row.v for row in current_rows]
    ref_t = [row.t for row in reference_rows]
    cur_t = [row.t for row in current_rows]
    ref_d = [row.disp for row in reference_rows]
    cur_d = [row.disp for row in current_rows]

    v_metrics = compute_metrics(ref_v, cur_v)
    t_metrics = compute_metrics(ref_t, cur_t)
    d_metrics = compute_metrics(ref_d, cur_d)

    print("field\tL2\tLinf\tmax_relative\tL2_relative")
    print(f"V\t{v_metrics[0]:.6e}\t{v_metrics[1]:.6e}\t{v_metrics[2]:.6e}\t{v_metrics[3]:.6e}")
    print(f"T\t{t_metrics[0]:.6e}\t{t_metrics[1]:.6e}\t{t_metrics[2]:.6e}\t{t_metrics[3]:.6e}")
    print(f"disp\t{d_metrics[0]:.6e}\t{d_metrics[1]:.6e}\t{d_metrics[2]:.6e}\t{d_metrics[3]:.6e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
