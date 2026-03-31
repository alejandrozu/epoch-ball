#!/usr/bin/env python3
"""Attempt-016 search over multistage symplectic transvection networks.

This is a constructive follow-up to attempt 015.

Attempt 015 ruled out the single-fold family

    x -> L(F_{lambda,d}(x))

where a single one-sided transvection fold is followed by one common linear
symplectic factor. The obstruction relied on a preserved symplectic 2-plane
slice. The next live escape route is to destroy that slice by composing several
noncommuting folds.

This script searches a new family:

    Phi = T_m o ... o T_1,

where each stage is the exact time-1 map of the piecewise Hamiltonian

    H_j(x) = -(lambda_j / 2) * max(omega(x, a_j) - d_j, 0)^2

with stage map

    T_j(x) = x + lambda_j * max(omega(x, a_j) - d_j, 0) * a_j.

Unlike attempt 015, different stages can use different isotropic directions
`a_j`, so the preserved-slice mechanism is no longer automatic.

The script performs a beam search over low-complexity stage dictionaries,
audits the best multistage candidates on sampled boundary points of the union
of k source balls, and reports the best sampled radius after allowing one final
global translation (found by Ritter's enclosing-ball heuristic).
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


J = np.array(
    [
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [-1.0, 0.0, 0.0, 0.0],
        [0.0, -1.0, 0.0, 0.0],
    ],
    dtype=np.float64,
)


@dataclass
class StageSpec:
    direction_name: str
    direction: list[float]
    quantile: float
    d_value: float
    lam: float


@dataclass
class CandidateSummary:
    stage_count: int
    coarse_radius: float
    coarse_ratio: float
    center: list[float]
    stages: list[StageSpec]


def omega_values(points: np.ndarray, a: np.ndarray) -> np.ndarray:
    return points @ (J @ a)


def apply_stage(points: np.ndarray, a: np.ndarray, lam: float, d_value: float) -> np.ndarray:
    levels = np.maximum(omega_values(points, a) - d_value, 0.0)
    return points + lam * levels[:, None] * a[None, :]


def sample_boundary_points(k: int, n_pts_per_ball: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    clouds = []
    for i in range(k):
        pts = rng.normal(size=(n_pts_per_ball, 4))
        pts /= np.linalg.norm(pts, axis=1, keepdims=True)
        pts[:, 0] += 3.0 * i
        clouds.append(pts)
    return np.concatenate(clouds, axis=0)


def normalize(v: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(v))
    if norm <= 1e-12:
        raise RuntimeError("zero direction")
    return v / norm


def direction_library() -> list[tuple[str, np.ndarray]]:
    e_q1 = np.array([1.0, 0.0, 0.0, 0.0])
    e_q2 = np.array([0.0, 1.0, 0.0, 0.0])
    e_p1 = np.array([0.0, 0.0, 1.0, 0.0])
    e_p2 = np.array([0.0, 0.0, 0.0, 1.0])
    raw = {
        "e_q1": e_q1,
        "e_q2": e_q2,
        "e_p1": e_p1,
        "e_p2": e_p2,
        "e_q1+e_p1": e_q1 + e_p1,
        "e_q1-e_p1": e_q1 - e_p1,
        "e_q2+e_p2": e_q2 + e_p2,
        "e_q2-e_p2": e_q2 - e_p2,
        "e_p1+e_p2": e_p1 + e_p2,
        "e_p1-e_p2": e_p1 - e_p2,
        "e_q1+e_q2": e_q1 + e_q2,
        "e_q1-e_q2": e_q1 - e_q2,
        "e_p1+e_q2": e_p1 + e_q2,
        "e_p2+e_q1": e_p2 + e_q1,
        "e_p1+e_q1+e_p2": e_p1 + e_q1 + e_p2,
        "e_p1+e_q1+e_q2": e_p1 + e_q1 + e_q2,
    }
    directions: list[tuple[str, np.ndarray]] = []
    for name, vec in raw.items():
        unit = normalize(vec.astype(np.float64))
        directions.append((name, unit))
        directions.append((f"-{name}", -unit))
    return directions


def ritter_enclosing_ball(points: np.ndarray) -> tuple[np.ndarray, float]:
    idx0 = 0
    p0 = points[idx0]
    idx1 = int(np.argmax(np.sum((points - p0) ** 2, axis=1)))
    p1 = points[idx1]
    idx2 = int(np.argmax(np.sum((points - p1) ** 2, axis=1)))
    p2 = points[idx2]
    center = 0.5 * (p1 + p2)
    radius = float(np.linalg.norm(p1 - center))

    for point in points:
        diff = point - center
        dist = float(np.linalg.norm(diff))
        if dist <= radius + 1e-12:
            continue
        new_radius = 0.5 * (radius + dist)
        center = center + ((new_radius - radius) / dist) * diff
        radius = new_radius

    return center, radius


def evaluate_points(points: np.ndarray) -> tuple[np.ndarray, float]:
    center, radius = ritter_enclosing_ball(points)
    return center, radius


def apply_network(points: np.ndarray, stages: list[StageSpec]) -> np.ndarray:
    out = points.copy()
    for stage in stages:
        out = apply_stage(out, np.array(stage.direction), stage.lam, stage.d_value)
    return out


def beam_search(
    points: np.ndarray,
    k: int,
    max_stages: int,
    beam_width: int,
    quantiles: list[float],
    lambdas: list[float],
) -> dict:
    directions = direction_library()
    beam: list[tuple[list[StageSpec], np.ndarray, np.ndarray, float]] = [([], points, np.zeros(4), float("inf"))]
    stage_rows = []

    for stage_index in range(max_stages):
        expanded: list[tuple[list[StageSpec], np.ndarray, np.ndarray, float]] = []
        for prefix, prefix_points, _, _ in beam:
            for name, direction in directions:
                values = omega_values(prefix_points, direction)
                for quantile in quantiles:
                    d_value = float(np.quantile(values, quantile))
                    for lam in lambdas:
                        mapped = apply_stage(prefix_points, direction, lam, d_value)
                        center, radius = evaluate_points(mapped)
                        stage = StageSpec(
                            direction_name=name,
                            direction=direction.tolist(),
                            quantile=quantile,
                            d_value=d_value,
                            lam=lam,
                        )
                        expanded.append((prefix + [stage], mapped, center, radius))

        expanded.sort(key=lambda row: row[3])
        beam = expanded[:beam_width]
        stage_rows.append(
            {
                "stage_index": stage_index + 1,
                "beam": [
                    {
                        "stage_count": len(stages),
                        "coarse_radius": radius,
                        "coarse_ratio": None if radius <= 0 else k / (radius**4),
                        "center": center.tolist(),
                        "stages": [asdict(stage) for stage in stages],
                    }
                    for stages, _, center, radius in beam
                ],
            }
        )

    final_rows = []
    for stages, _, center, radius in beam:
        final_rows.append(
            CandidateSummary(
                stage_count=len(stages),
                coarse_radius=radius,
                coarse_ratio=0.0 if radius <= 0 else k / (radius**4),
                center=center.tolist(),
                stages=stages,
            )
        )
    final_rows.sort(key=lambda row: row.coarse_radius)
    return {
        "search_settings": {
            "max_stages": max_stages,
            "beam_width": beam_width,
            "quantiles": quantiles,
            "lambdas": lambdas,
            "direction_count": len(directions),
        },
        "stage_rows": stage_rows,
        "top_candidates": [asdict(row) for row in final_rows],
    }


def refine_candidates(
    candidates: list[CandidateSummary],
    k: int,
    n_pts_per_ball: int,
    seeds: list[int],
) -> list[dict]:
    refined = []
    for candidate in candidates:
        seed_rows = []
        radii = []
        for seed in seeds:
            points = sample_boundary_points(k=k, n_pts_per_ball=n_pts_per_ball, seed=seed)
            mapped = apply_network(points, candidate.stages)
            center, radius = evaluate_points(mapped)
            radii.append(radius)
            seed_rows.append(
                {
                    "seed": seed,
                    "radius": radius,
                    "ratio": k / (radius**4),
                    "center": center.tolist(),
                }
            )
        refined.append(
            {
                "stage_count": candidate.stage_count,
                "coarse_radius": candidate.coarse_radius,
                "coarse_ratio": candidate.coarse_ratio,
                "stages": [asdict(stage) for stage in candidate.stages],
                "refine_seed_rows": seed_rows,
                "worst_radius": max(radii),
                "worst_ratio": k / (max(radii) ** 4),
                "best_radius": min(radii),
                "best_ratio": k / (min(radii) ** 4),
                "mean_radius": float(np.mean(radii)),
                "mean_ratio": float(np.mean([k / (radius**4) for radius in radii])),
            }
        )
    refined.sort(key=lambda row: row["worst_radius"])
    return refined


def random_network_search(
    points: np.ndarray,
    k: int,
    n_candidates: int,
    stage_count: int,
    seed: int,
    quantiles: list[float],
    lambdas: list[float],
    keep_top_n: int,
) -> dict:
    rng = np.random.default_rng(seed)
    directions = direction_library()
    rows: list[CandidateSummary] = []

    for _ in range(n_candidates):
        mapped = points.copy()
        stages: list[StageSpec] = []
        for _stage_idx in range(stage_count):
            name, direction = directions[int(rng.integers(len(directions)))]
            quantile = float(quantiles[int(rng.integers(len(quantiles)))])
            lam = float(lambdas[int(rng.integers(len(lambdas)))])
            values = omega_values(mapped, direction)
            d_value = float(np.quantile(values, quantile))
            stages.append(
                StageSpec(
                    direction_name=name,
                    direction=direction.tolist(),
                    quantile=quantile,
                    d_value=d_value,
                    lam=lam,
                )
            )
            mapped = apply_stage(mapped, direction, lam, d_value)
        center, radius = evaluate_points(mapped)
        rows.append(
            CandidateSummary(
                stage_count=stage_count,
                coarse_radius=radius,
                coarse_ratio=0.0 if radius <= 0 else k / (radius**4),
                center=center.tolist(),
                stages=stages,
            )
        )

    rows.sort(key=lambda row: row.coarse_radius)
    return {
        "search_settings": {
            "n_candidates": n_candidates,
            "stage_count": stage_count,
            "seed": seed,
            "direction_count": len(directions),
            "quantiles": quantiles,
            "lambdas": lambdas,
        },
        "top_candidates": [asdict(row) for row in rows[:keep_top_n]],
    }


def build_report(args: argparse.Namespace) -> dict:
    coarse_points = sample_boundary_points(k=args.k, n_pts_per_ball=args.coarse_n_pts_per_ball, seed=args.coarse_seed)
    quantiles = [0.2, 0.35, 0.5, 0.65, 0.8]
    lambdas = [-2.0, -1.0, -0.5, 0.5, 1.0, 2.0]
    search = beam_search(
        points=coarse_points,
        k=args.k,
        max_stages=args.max_stages,
        beam_width=args.beam_width,
        quantiles=quantiles,
        lambdas=lambdas,
    )
    random_search = random_network_search(
        points=coarse_points,
        k=args.k,
        n_candidates=args.random_candidates,
        stage_count=args.random_stage_count,
        seed=args.random_seed,
        quantiles=quantiles,
        lambdas=lambdas,
        keep_top_n=args.refine_top_n,
    )
    top_candidates = []
    seen = set()
    combined_rows = search["top_candidates"][: args.refine_top_n] + random_search["top_candidates"][: args.refine_top_n]
    for row in combined_rows:
        key = tuple(
            (
                stage_row["direction_name"],
                round(stage_row["quantile"], 6),
                round(stage_row["d_value"], 6),
                round(stage_row["lam"], 6),
            )
            for stage_row in row["stages"]
        )
        if key in seen:
            continue
        seen.add(key)
        top_candidates.append(
            CandidateSummary(
                stage_count=row["stage_count"],
                coarse_radius=row["coarse_radius"],
                coarse_ratio=row["coarse_ratio"],
                center=row["center"],
                stages=[StageSpec(**stage_row) for stage_row in row["stages"]],
            )
        )
    refined = refine_candidates(
        candidates=top_candidates,
        k=args.k,
        n_pts_per_ball=args.refined_n_pts_per_ball,
        seeds=args.refine_seeds,
    )
    best = refined[0]
    return {
        "family": "multistage_one_sided_symplectic_transvection_network",
        "stage_formula": "T(x) = x + lambda * max(omega(x,a) - d, 0) * a",
        "coarse_search": search,
        "random_search": random_search,
        "refined_candidates": refined,
        "best_candidate": best,
        "conclusion": (
            "This new multistage transvection family destroys the fixed-slice mechanism of "
            "attempt 015 and gives a concrete constructive search space. The saved best "
            "candidate is only boundary-sampled evidence, not a verified packing."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--coarse-n-pts-per-ball", type=int, default=96)
    parser.add_argument("--refined-n-pts-per-ball", type=int, default=384)
    parser.add_argument("--coarse-seed", type=int, default=123)
    parser.add_argument("--refine-seeds", type=int, nargs="+", default=[11, 29, 47])
    parser.add_argument("--max-stages", type=int, default=3)
    parser.add_argument("--beam-width", type=int, default=8)
    parser.add_argument("--refine-top-n", type=int, default=8)
    parser.add_argument("--random-candidates", type=int, default=3000)
    parser.add_argument("--random-stage-count", type=int, default=3)
    parser.add_argument("--random-seed", type=int, default=2026)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/research_loop/attempt_016_multistage_transvection_search.json"),
    )
    args = parser.parse_args()

    report = build_report(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True))
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
