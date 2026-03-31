#!/usr/bin/env python3
"""Affine symplectic block packings for equal 4-balls.

This is a new constructive family, distinct from the toric/action-angle
packings of attempt_003.

Each block of size a*b is realized by:
  1. translating each source ball in the block to the origin,
  2. applying the affine symplectic map
       (q1, q2, p1, p2) -> (a q1, b q2, p1 / a, p2 / b),
  3. translating the resulting ellipsoids to an a-by-b grid in the thin
     (p1, p2)-directions.

The whole block then occupies the q-rectangle
    [c1-a, c1+a] x [c2-b, c2+b]
and the p-box
    [-1, 1] x [-1, 1].

Different blocks are made disjoint by placing their q-rectangles disjointly.
If all q-rectangles fit in the q-plane disk of radius rho, the whole 4D image
fits in the Euclidean 4-ball of radius sqrt(rho^2 + 2).
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import asdict, dataclass

import torch


torch.set_default_dtype(torch.float64)


@dataclass
class AffineBlockPlacement:
    n_balls: int
    side_a: int
    side_b: int
    q_center: list[float]


@dataclass
class AffinePackingSummary:
    k: int
    rho: float
    target_radius: float
    exact_radius: float
    predicted_ratio: float
    exact_ratio: float
    sampled_radius: float
    sampled_ratio: float
    block_limit: int
    grid_step: float
    partition: list[int]
    shape_choice: list[list[int]]
    placements: list[AffineBlockPlacement]
    seed: int
    n_pts_per_ball: int


def sample_boundary_points(k: int, n_pts_per_ball: int, seed: int) -> torch.Tensor:
    torch.manual_seed(seed)
    pts = []
    for i in range(k):
        x = torch.randn(n_pts_per_ball, 4)
        x = x / x.norm(dim=1, keepdim=True)
        x[:, 0] += 3.0 * i
        pts.append(x)
    return torch.cat(pts, dim=0)


def integer_partitions(total: int, minimum: int = 1) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []

    def recurse(remaining: int, start: int, prefix: list[int]) -> None:
        if remaining == 0:
            out.append(tuple(prefix))
            return
        for value in range(start, remaining + 1):
            prefix.append(value)
            recurse(remaining - value, value, prefix)
            prefix.pop()

    recurse(total, minimum, [])
    return out


def factor_shapes(n: int) -> list[tuple[int, int]]:
    shapes = set()
    for a in range(1, n + 1):
        if n % a != 0:
            continue
        b = n // a
        if a >= b:
            shapes.add((a, b))
        else:
            shapes.add((b, a))
    return sorted(shapes, reverse=True)


def rectangle_corners(cx: float, cy: float, a: int, b: int) -> list[tuple[float, float]]:
    return [
        (cx - a, cy - b),
        (cx - a, cy + b),
        (cx + a, cy - b),
        (cx + a, cy + b),
    ]


def rectangle_inside_q_disk(cx: float, cy: float, a: int, b: int, rho: float, tol: float = 1e-12) -> bool:
    return all(x * x + y * y <= rho * rho + tol for x, y in rectangle_corners(cx, cy, a, b))


def rectangles_overlap(
    cx1: float,
    cy1: float,
    a1: int,
    b1: int,
    cx2: float,
    cy2: float,
    a2: int,
    b2: int,
    tol: float = 1e-12,
) -> bool:
    if abs(cx1 - cx2) >= a1 + a2 - tol:
        return False
    if abs(cy1 - cy2) >= b1 + b2 - tol:
        return False
    return True


def candidate_centers(a: int, b: int, rho: float, grid_step: float) -> list[tuple[float, float]]:
    if rho < math.sqrt(a * a + b * b):
        return []
    max_x = int(math.floor((rho - a) / grid_step + 1e-12))
    max_y = int(math.floor((rho - b) / grid_step + 1e-12))
    out = []
    for ix in range(-max_x, max_x + 1):
        cx = round(ix * grid_step, 10)
        for iy in range(-max_y, max_y + 1):
            cy = round(iy * grid_step, 10)
            if rectangle_inside_q_disk(cx, cy, a, b, rho=rho):
                out.append((cx, cy))
    return out


def q_corner_radius_sq(placements: list[AffineBlockPlacement]) -> float:
    max_val = 0.0
    for placement in placements:
        cx, cy = placement.q_center
        for x, y in rectangle_corners(cx, cy, placement.side_a, placement.side_b):
            max_val = max(max_val, x * x + y * y)
    return max_val


def exact_ball_radius(center: list[float], scales: list[float], tol: float = 1e-12) -> float:
    eigvals = [scale * scale for scale in scales]
    drive = [center_i * scale for center_i, scale in zip(center, scales)]
    center_sq = sum(value * value for value in center)
    lam_max = max(eigvals)

    if max(abs(value) for value in drive) <= tol:
        return math.sqrt(center_sq + lam_max)

    def constraint(mu: float) -> float:
        total = 0.0
        for eigval, drive_i in zip(eigvals, drive):
            total += (drive_i / (mu - eigval)) ** 2
        return total

    lo = lam_max + 1e-12
    hi = max(lo + 1.0, lam_max + 1.0)
    while constraint(hi) > 1.0:
        hi *= 2.0

    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if constraint(mid) > 1.0:
            lo = mid
        else:
            hi = mid

    mu = hi
    correction = 0.0
    for eigval, drive_i in zip(eigvals, drive):
        correction += (drive_i * drive_i) / (mu - eigval)
    return math.sqrt(center_sq + mu + correction)


def parse_rho_values(rho_min: float, rho_max: float, rho_step: float) -> list[float]:
    count = int(math.floor((rho_max - rho_min) / rho_step + 1e-12))
    return [round(rho_min + i * rho_step, 10) for i in range(count + 1)]


def search_affine_blocks(k: int, rho_values: list[float], grid_step: float, block_limit: int) -> dict:
    partitions = [p for p in integer_partitions(k) if len(p) <= block_limit]
    factor_cache = {n: factor_shapes(n) for n in range(1, k + 1)}
    best_result: dict | None = None

    for rho in rho_values:
        for partition in partitions:
            shape_lists = [factor_cache[n] for n in partition]
            for chosen_shapes in itertools.product(*shape_lists):
                blocks = sorted(
                    zip(partition, chosen_shapes),
                    key=lambda item: (item[1][0] * item[1][1], item[1][0] * item[1][0] + item[1][1] * item[1][1]),
                    reverse=True,
                )
                candidates = []
                feasible = True
                for n_balls, (a, b) in blocks:
                    centers = candidate_centers(a, b, rho=rho, grid_step=grid_step)
                    if not centers:
                        feasible = False
                        break
                    candidates.append((n_balls, a, b, centers))
                if not feasible:
                    continue

                placed: list[AffineBlockPlacement] = []

                def backtrack(index: int) -> bool:
                    if index == len(candidates):
                        return True
                    n_balls, a, b, centers = candidates[index]
                    for cx, cy in centers:
                        if any(
                            rectangles_overlap(cx, cy, a, b, other.q_center[0], other.q_center[1], other.side_a, other.side_b)
                            for other in placed
                        ):
                            continue
                        placed.append(AffineBlockPlacement(n_balls=n_balls, side_a=a, side_b=b, q_center=[cx, cy]))
                        if backtrack(index + 1):
                            return True
                        placed.pop()
                    return False

                if not backtrack(0):
                    continue

                corner_radius_sq = q_corner_radius_sq(placed)
                target_radius = math.sqrt(corner_radius_sq + 2.0)
                exact_radius = exact_packing_radius(placed, k)
                result = {
                    "found": True,
                    "k": k,
                    "rho": rho,
                    "target_radius": target_radius,
                    "exact_radius": exact_radius,
                    "predicted_ratio": k / (target_radius**4),
                    "exact_ratio": k / (exact_radius**4),
                    "grid_step": grid_step,
                    "block_limit": block_limit,
                    "partition": list(partition),
                    "shape_choice": [[a, b] for a, b in chosen_shapes],
                    "placements": [asdict(block) for block in placed],
                }
                if best_result is None or result["exact_ratio"] > best_result["exact_ratio"]:
                    best_result = result
                break
            if best_result is not None and abs(best_result["rho"] - rho) <= 1e-12:
                break
        if best_result is not None and abs(best_result["rho"] - rho) <= 1e-12:
            break

    if best_result is not None:
        return best_result
    return {
        "found": False,
        "k": k,
        "rho_values": rho_values,
        "grid_step": grid_step,
        "block_limit": block_limit,
    }


def build_assignment_table(placements: list[AffineBlockPlacement], k: int) -> list[dict]:
    table: list[dict] = []
    ball_index = 0
    for placement in placements:
        a = placement.side_a
        b = placement.side_b
        cx, cy = placement.q_center
        for u in range(a):
            for v in range(b):
                if ball_index >= k:
                    break
                table.append(
                    {
                        "ball_index": ball_index,
                        "a": a,
                        "b": b,
                        "q_center": [cx, cy],
                        "p_center": [(-a + 1 + 2 * u) / float(a), (-b + 1 + 2 * v) / float(b)],
                    }
                )
                ball_index += 1
    if ball_index != k:
        raise RuntimeError(f"assignment table expected {k} balls, got {ball_index}")
    return table


def exact_packing_radius(placements: list[AffineBlockPlacement], k: int) -> float:
    assignment_table = build_assignment_table(placements, k)
    best = 0.0
    for entry in assignment_table:
        a = float(entry["a"])
        b = float(entry["b"])
        radius = exact_ball_radius(
            center=[
                float(entry["q_center"][0]),
                float(entry["q_center"][1]),
                float(entry["p_center"][0]),
                float(entry["p_center"][1]),
            ],
            scales=[a, b, 1.0 / a, 1.0 / b],
        )
        best = max(best, radius)
    return best


def direct_affine_map(points: torch.Tensor, result: dict) -> torch.Tensor:
    if not result.get("found", False):
        raise ValueError("packing result does not contain a found placement")
    placements = [AffineBlockPlacement(**entry) for entry in result["placements"]]
    assignment_table = build_assignment_table(placements, int(result["k"]))

    outputs = torch.empty_like(points)
    for entry in assignment_table:
        i = entry["ball_index"]
        mask = (points[:, 0] >= 3.0 * i - 1.5) & (points[:, 0] <= 3.0 * i + 1.5)
        local = points[mask].clone()
        local[:, 0] -= 3.0 * i
        a = float(entry["a"])
        b = float(entry["b"])
        local_q1 = a * local[:, 0]
        local_q2 = b * local[:, 1]
        local_p1 = local[:, 2] / a + float(entry["p_center"][0])
        local_p2 = local[:, 3] / b + float(entry["p_center"][1])
        outputs[mask, 0] = local_q1 + float(entry["q_center"][0])
        outputs[mask, 1] = local_q2 + float(entry["q_center"][1])
        outputs[mask, 2] = local_p1
        outputs[mask, 3] = local_p2
    return outputs


def verify_affine_packing(result: dict, n_pts_per_ball: int, seed: int) -> AffinePackingSummary:
    points = sample_boundary_points(int(result["k"]), n_pts_per_ball=n_pts_per_ball, seed=seed)
    mapped = direct_affine_map(points, result)
    sampled_radius = float(mapped.norm(dim=1).max())
    sampled_ratio = int(result["k"]) / (sampled_radius**4)
    placements = [AffineBlockPlacement(**entry) for entry in result["placements"]]
    return AffinePackingSummary(
        k=int(result["k"]),
        rho=float(result["rho"]),
        target_radius=float(result["target_radius"]),
        exact_radius=float(result["exact_radius"]),
        predicted_ratio=float(result["predicted_ratio"]),
        exact_ratio=float(result["exact_ratio"]),
        sampled_radius=sampled_radius,
        sampled_ratio=sampled_ratio,
        block_limit=int(result["block_limit"]),
        grid_step=float(result["grid_step"]),
        partition=list(result["partition"]),
        shape_choice=[list(x) for x in result["shape_choice"]],
        placements=placements,
        seed=seed,
        n_pts_per_ball=n_pts_per_ball,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--rho-min", type=float, required=True)
    parser.add_argument("--rho-max", type=float, required=True)
    parser.add_argument("--rho-step", type=float, default=0.1)
    parser.add_argument("--grid-step", type=float, default=0.1)
    parser.add_argument("--block-limit", type=int, default=4)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--n-pts", type=int, default=2000)
    args = parser.parse_args()

    result = search_affine_blocks(
        k=args.k,
        rho_values=parse_rho_values(args.rho_min, args.rho_max, args.rho_step),
        grid_step=args.grid_step,
        block_limit=args.block_limit,
    )
    if not result.get("found", False):
        print(json.dumps(result, indent=2, sort_keys=True))
        return

    summary = verify_affine_packing(result, n_pts_per_ball=args.n_pts, seed=args.seed)
    payload = {
        "search_result": result,
        "verification": asdict(summary),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
