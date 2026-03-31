#!/usr/bin/env python3
"""Search simple affine-toric covering blocks for equal-ball packings in B^4.

A 2x2 integer matrix A with det(A)=n>0 defines a degree-n torus covering.
Using the standard cotangent-lift construction on action-angle coordinates,
this yields an explicit toric block packing of n unit balls into the toric
domain associated with a translated copy of A Delta(1), where
    Delta(1) = {(x1, x2) : x1 >= 0, x2 >= 0, x1 + x2 <= 1}.

This module searches the simplest such one-block family:
    T = t + A Delta(1) subset Delta(mu),
where t is chosen to place the triangle in the positive quadrant. The minimal
mu for a fixed A can be computed in closed form.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass


@dataclass
class OneBlockCandidate:
    determinant: int
    matrix: list[list[int]]
    translation: list[float]
    required_capacity_mu: float
    achieved_ratio: float
    triangle_vertices: list[list[float]]


def triangle_vertices_from_matrix(a11: int, a12: int, a21: int, a22: int) -> list[tuple[int, int]]:
    return [(0, 0), (a11, a21), (a12, a22)]


def minimal_translation_and_capacity(vertices: list[tuple[int, int]]) -> tuple[tuple[float, float], float]:
    min_x = min(v[0] for v in vertices)
    min_y = min(v[1] for v in vertices)
    tx = float(-min_x) if min_x < 0 else 0.0
    ty = float(-min_y) if min_y < 0 else 0.0
    mu = max((tx + ty + vx + vy) for vx, vy in vertices)
    return (tx, ty), float(mu)


def candidate_from_matrix(a11: int, a12: int, a21: int, a22: int) -> OneBlockCandidate | None:
    det = a11 * a22 - a12 * a21
    if det <= 0:
        return None
    vertices = triangle_vertices_from_matrix(a11, a12, a21, a22)
    translation, mu = minimal_translation_and_capacity(vertices)
    shifted_vertices = [[translation[0] + vx, translation[1] + vy] for vx, vy in vertices]
    return OneBlockCandidate(
        determinant=det,
        matrix=[[a11, a12], [a21, a22]],
        translation=[translation[0], translation[1]],
        required_capacity_mu=mu,
        achieved_ratio=det / (mu * mu),
        triangle_vertices=shifted_vertices,
    )


def best_one_block_candidate(det_target: int, entry_bound: int) -> OneBlockCandidate:
    best: OneBlockCandidate | None = None
    for a11 in range(-entry_bound, entry_bound + 1):
        for a12 in range(-entry_bound, entry_bound + 1):
            for a21 in range(-entry_bound, entry_bound + 1):
                for a22 in range(-entry_bound, entry_bound + 1):
                    cand = candidate_from_matrix(a11, a12, a21, a22)
                    if cand is None or cand.determinant != det_target:
                        continue
                    if best is None:
                        best = cand
                        continue
                    if cand.required_capacity_mu < best.required_capacity_mu - 1e-12:
                        best = cand
                        continue
                    if abs(cand.required_capacity_mu - best.required_capacity_mu) <= 1e-12 and cand.achieved_ratio > best.achieved_ratio:
                        best = cand
    if best is None:
        raise RuntimeError(f"no determinant-{det_target} matrix found with entry bound {entry_bound}")
    return best


def summary_row(det_target: int, entry_bound: int) -> dict:
    best = best_one_block_candidate(det_target, entry_bound)
    return {
        "k": det_target,
        "sqrt_k": math.sqrt(det_target),
        "best_one_block": asdict(best),
        "mu_gap_vs_sqrt_k": best.required_capacity_mu - math.sqrt(det_target),
        "ratio_gap_vs_full": 1.0 - best.achieved_ratio,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=None)
    parser.add_argument("--entry-bound", type=int, default=8)
    parser.add_argument("--table-k-max", type=int, default=None)
    args = parser.parse_args()

    if args.table_k_max is not None:
        rows = [summary_row(k, args.entry_bound) for k in range(1, args.table_k_max + 1)]
        print(json.dumps({"rows": rows}, indent=2, sort_keys=True))
        return

    if args.k is None:
        parser.error("pass either --k or --table-k-max")
    print(json.dumps(summary_row(args.k, args.entry_bound), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
