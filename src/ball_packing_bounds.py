#!/usr/bin/env python3
"""Bounds and feasibility checks for equal-ball packings in B^4.

This module works in normalized capacity coordinates:
 - each source ball has capacity 1
 - the target ball has capacity mu

In these units, the filled volume fraction is k / mu^2.

For k >= 3, feasibility is checked using the standard Cremona reduction
algorithm for vectors (mu; 1, ..., 1). This is the algebraic-geometric
criterion used in the 4-dimensional equal-ball packing problem.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass


TOL = 1e-12


@dataclass
class ReductionSummary:
    input_mu: float
    final_mu: float
    final_weights: list[float]
    is_reduced: bool
    positive_entries: bool
    positive_square: bool
    iterations: int
    defect: float


def sort_weights(weights: list[float]) -> list[float]:
    return sorted(weights, reverse=True)


def is_reduced(mu: float, weights: list[float], tol: float = TOL) -> bool:
    weights = sort_weights(weights)
    if len(weights) < 3:
        return True
    return weights[0] + weights[1] + weights[2] <= mu + tol


def positive_square(mu: float, weights: list[float], tol: float = TOL) -> bool:
    return mu * mu - sum(w * w for w in weights) >= -tol


def standard_cremona_step(mu: float, weights: list[float]) -> tuple[float, list[float]]:
    weights = sort_weights(weights)
    defect = weights[0] + weights[1] + weights[2] - mu
    if defect <= 0.0:
        return mu, weights
    new_mu = mu - defect
    new_weights = weights[:]
    new_weights[0] -= defect
    new_weights[1] -= defect
    new_weights[2] -= defect
    return new_mu, sort_weights(new_weights)


def cremona_reduce(mu: float, weights: list[float], tol: float = TOL, max_steps: int = 10000) -> ReductionSummary:
    cur_mu = float(mu)
    cur_weights = sort_weights([float(w) for w in weights])
    iterations = 0
    while iterations < max_steps and len(cur_weights) >= 3 and not is_reduced(cur_mu, cur_weights, tol=tol):
        cur_mu, cur_weights = standard_cremona_step(cur_mu, cur_weights)
        iterations += 1
    defect = 0.0
    if len(cur_weights) >= 3:
        defect = cur_weights[0] + cur_weights[1] + cur_weights[2] - cur_mu
    return ReductionSummary(
        input_mu=float(mu),
        final_mu=cur_mu,
        final_weights=cur_weights,
        is_reduced=is_reduced(cur_mu, cur_weights, tol=tol),
        positive_entries=min([cur_mu] + cur_weights) >= -tol,
        positive_square=positive_square(cur_mu, cur_weights, tol=tol),
        iterations=iterations,
        defect=defect,
    )


def equal_ball_feasible(mu: float, k: int, tol: float = TOL) -> bool:
    if k <= 0:
        raise ValueError("k must be positive")
    if k == 1:
        return mu >= 1.0 - tol
    if k == 2:
        return mu >= 2.0 - tol
    summary = cremona_reduce(mu, [1.0] * k, tol=tol)
    return summary.positive_entries and summary.is_reduced and summary.positive_square


def equal_ball_optimal_capacity(k: int, tol: float = 1e-10) -> float:
    if k == 1:
        return 1.0
    if k == 2:
        return 2.0
    lo = math.sqrt(k)
    hi = max(4.0, lo + 1.0)
    while not equal_ball_feasible(hi, k):
        hi *= 2.0
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if equal_ball_feasible(mid, k, tol=tol):
            hi = mid
        else:
            lo = mid
    return hi


def equal_ball_volume_fraction(k: int) -> float:
    mu = equal_ball_optimal_capacity(k)
    return k / (mu * mu)


def task_capacity(k: int, eps: float) -> float:
    if not (0.0 < eps < 1.0):
        raise ValueError("eps must lie in (0,1)")
    return math.sqrt(k / (1.0 - eps))


def task_summary(k: int, eps: float, achieved_ratio: float | None = None) -> dict:
    mu_task = task_capacity(k, eps)
    full_mu = equal_ball_optimal_capacity(k)
    full_fraction = k / (full_mu * full_mu)
    reduction = cremona_reduce(mu_task, [1.0] * k)
    summary = {
        "k": k,
        "eps": eps,
        "task_target_ratio": 1.0 - eps,
        "task_target_capacity_mu": mu_task,
        "task_vector_reduction": asdict(reduction),
        "full_packing_infimum_capacity_mu": full_mu,
        "full_packing_volume_fraction": full_fraction,
        "abstract_minimum_eps": max(0.0, 1.0 - full_fraction),
        "volume_surplus_over_full": mu_task * mu_task - k,
    }
    if achieved_ratio is not None:
        achieved_mu = math.sqrt(k / achieved_ratio)
        summary["achieved_ratio"] = achieved_ratio
        summary["achieved_capacity_mu"] = achieved_mu
        summary["mu_gap_vs_task"] = achieved_mu - mu_task
        summary["mu_gap_factor_vs_full"] = achieved_mu / full_mu
    return summary


def reduction_table(k_max: int) -> list[dict]:
    rows = []
    for k in range(1, k_max + 1):
        mu = equal_ball_optimal_capacity(k)
        fraction = k / (mu * mu)
        rows.append(
            {
                "k": k,
                "optimal_capacity_mu": mu,
                "volume_fraction": fraction,
                "abstract_minimum_eps": max(0.0, 1.0 - fraction),
                "full_packing": abs(fraction - 1.0) <= 1e-9,
                "reduction_at_infimum": asdict(cremona_reduce(mu + 1e-9, [1.0] * k)) if k >= 3 else None,
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table-k-max", type=int, default=None, help="Emit a table of optimal equal-ball packings up to this k.")
    parser.add_argument("--k", type=int, default=None)
    parser.add_argument("--eps", type=float, default=None)
    parser.add_argument("--achieved-ratio", type=float, default=None)
    args = parser.parse_args()

    if args.table_k_max is not None:
        print(json.dumps({"table": reduction_table(args.table_k_max)}, indent=2, sort_keys=True))
        return

    if args.k is None or args.eps is None:
        parser.error("either pass --table-k-max or pass both --k and --eps")
    print(json.dumps(task_summary(args.k, args.eps, achieved_ratio=args.achieved_ratio), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
