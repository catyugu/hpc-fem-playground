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


def compute_metrics(reference: list[float], current: list[float]) -> tuple[float, float, float]:
    if len(reference) != len(current):
        raise ValueError("Result vectors have different lengths")

    n = len(reference)
    if n == 0:
        raise ValueError("No rows available for comparison")

    sq_sum = 0.0
    max_abs = 0.0
    max_rel = 0.0
    for i in range(n):
        diff = current[i] - reference[i]
        abs_diff = abs(diff)
        sq_sum += diff * diff
        if abs_diff > max_abs:
            max_abs = abs_diff

        denom = max(abs(reference[i]), 1e-16)
        rel = abs_diff / denom
        if rel > max_rel:
            max_rel = rel

    l2 = math.sqrt(sq_sum / n)
    return l2, max_abs, max_rel


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare MPFEM and COMSOL results")
    parser.add_argument("reference", type=Path, help="COMSOL reference result.txt")
    parser.add_argument("current", type=Path, help="MPFEM output result file")
    args = parser.parse_args()

    reference_rows = parse_result_file(args.reference)
    current_rows = parse_result_file(args.current)

    if len(reference_rows) != len(current_rows):
        raise ValueError(
            f"Row count mismatch: reference={len(reference_rows)} current={len(current_rows)}"
        )

    ref_v = [row.v for row in reference_rows]
    cur_v = [row.v for row in current_rows]
    ref_t = [row.t for row in reference_rows]
    cur_t = [row.t for row in current_rows]
    ref_d = [row.disp for row in reference_rows]
    cur_d = [row.disp for row in current_rows]

    v_metrics = compute_metrics(ref_v, cur_v)
    t_metrics = compute_metrics(ref_t, cur_t)
    d_metrics = compute_metrics(ref_d, cur_d)

    print("field\tL2\tLinf\tmax_relative")
    print(f"V\t{v_metrics[0]:.6e}\t{v_metrics[1]:.6e}\t{v_metrics[2]:.6e}")
    print(f"T\t{t_metrics[0]:.6e}\t{t_metrics[1]:.6e}\t{t_metrics[2]:.6e}")
    print(f"disp\t{d_metrics[0]:.6e}\t{d_metrics[1]:.6e}\t{d_metrics[2]:.6e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
