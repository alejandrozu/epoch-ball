#!/usr/bin/env python3
"""Explicit toric square-packing baseline for equal 4-balls.

This module does not produce a Hamiltonian. It provides a direct symplectic
embedding candidate on the disjoint source union by translating each ball to
the origin, converting to action-angle coordinates, applying the standard
l^2-slot toric packing, and converting back to Cartesian coordinates.

In the normalized coordinates used here,
    x_j = q_j^2 + p_j^2,
    y_j = angle(q_j, p_j) / (2*pi),
so the source unit ball satisfies x_1 + x_2 <= 1 and the target ball of
capacity mu has Euclidean radius sqrt(mu). The square-packing map with
l = ceil(sqrt(k)) sends each source ball to one of the l^2 angle slots inside
the target ball of capacity l, hence achieves the exact theoretical ratio
    k / l^2.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass

import torch


torch.set_default_dtype(torch.float64)


TWO_PI = 2.0 * math.pi


@dataclass
class PackingSummary:
    k: int
    l: int
    target_capacity_mu: float
    target_radius: float
    predicted_ratio: float
    sampled_radius: float
    sampled_ratio: float
    max_excess_over_predicted_radius: float
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


def infer_ball_indices(q1: torch.Tensor, k: int) -> torch.Tensor:
    indices = torch.round(q1 / 3.0).to(torch.int64)
    return indices.clamp_(0, k - 1)


def cartesian_to_action_angle(q: torch.Tensor, p: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    x1 = q[:, 0].square() + p[:, 0].square()
    x2 = q[:, 1].square() + p[:, 1].square()
    y1 = torch.atan2(p[:, 0], q[:, 0]) / TWO_PI
    y2 = torch.atan2(p[:, 1], q[:, 1]) / TWO_PI
    y = torch.stack((torch.remainder(y1, 1.0), torch.remainder(y2, 1.0)), dim=1)
    x = torch.stack((x1, x2), dim=1)
    return x, y


def action_angle_to_cartesian(x: torch.Tensor, y: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    r = torch.sqrt(torch.clamp_min(x, 0.0))
    angles = TWO_PI * y
    cosines = torch.cos(angles)
    sines = torch.sin(angles)
    q = r * cosines
    p = r * sines
    return q, p


def square_slot_offsets(indices: torch.Tensor, l: int) -> torch.Tensor:
    row = torch.div(indices, l, rounding_mode="floor")
    col = torch.remainder(indices, l)
    return torch.stack((row, col), dim=1).to(torch.get_default_dtype())


def translate_union_to_origin(points: torch.Tensor, k: int) -> tuple[torch.Tensor, torch.Tensor]:
    q = points[:, :2].clone()
    p = points[:, 2:].clone()
    indices = infer_ball_indices(q[:, 0], k)
    q[:, 0] -= 3.0 * indices.to(q.dtype)
    return torch.cat((q, p), dim=1), indices


def square_packing_map(points: torch.Tensor, k: int) -> tuple[torch.Tensor, dict]:
    if points.ndim != 2 or points.shape[1] != 4:
        raise ValueError(f"expected points of shape (N, 4), got {tuple(points.shape)}")
    if k < 1:
        raise ValueError("k must be positive")

    centered, indices = translate_union_to_origin(points, k)
    q = centered[:, :2]
    p = centered[:, 2:]
    x, y = cartesian_to_action_angle(q, p)

    l = math.ceil(math.sqrt(k))
    offsets = square_slot_offsets(indices, l)
    x_out = float(l) * x
    y_out = (y + offsets) / float(l)
    q_out, p_out = action_angle_to_cartesian(x_out, y_out)
    return torch.cat((q_out, p_out), dim=1), {
        "k": k,
        "l": l,
        "target_capacity_mu": float(l),
        "target_radius": math.sqrt(l),
        "predicted_ratio": k / float(l * l),
    }


def verify_square_packing(k: int, n_pts_per_ball: int, seed: int) -> PackingSummary:
    points = sample_boundary_points(k, n_pts_per_ball, seed)
    mapped, meta = square_packing_map(points, k)
    radii = mapped.norm(dim=1)
    sampled_radius = float(radii.max())
    sampled_ratio = k / (sampled_radius**4)
    predicted_radius = meta["target_radius"]
    return PackingSummary(
        k=k,
        l=meta["l"],
        target_capacity_mu=meta["target_capacity_mu"],
        target_radius=predicted_radius,
        predicted_ratio=meta["predicted_ratio"],
        sampled_radius=sampled_radius,
        sampled_ratio=sampled_ratio,
        max_excess_over_predicted_radius=sampled_radius - predicted_radius,
        seed=seed,
        n_pts_per_ball=n_pts_per_ball,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--n-pts", type=int, default=2000, help="Boundary points sampled per source ball.")
    args = parser.parse_args()

    summary = verify_square_packing(k=args.k, n_pts_per_ball=args.n_pts, seed=args.seed)
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
