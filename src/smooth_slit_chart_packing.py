#!/usr/bin/env python3
"""Smooth nonlinear slit-chart packing for equal 4-balls.

This attempt uses a genuinely nonlinear smooth planar canonical chart instead
of the singular action-angle map from attempt_003 and instead of the affine
local-piece family from attempt_004.

For a fixed exterior point a > 1 on the q-axis, each unit disk in a symplectic
plane is mapped to a rectangle [0,1] x [0,pi] by:
  1. smooth "exterior log-polar" coordinates around the point (a,0),
  2. an exact triangular flattening of the resulting curvilinear band.

The product chart on the two symplectic planes is smooth and symplectic on the
whole 4-ball. We then apply the standard l^2 slot map in the rectangle
coordinates and map back.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass

import torch


torch.set_default_dtype(torch.float64)


PI = math.pi


@dataclass
class PackingSummary:
    k: int
    l: int
    a: float
    sampled_radius: float
    sampled_ratio: float
    seed: int
    n_pts_per_ball: int
    chart_roundtrip_max_error: float
    chart_area_total: float


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


def linear_interp(query: torch.Tensor, x_grid: torch.Tensor, y_grid: torch.Tensor) -> torch.Tensor:
    idx = torch.bucketize(query, x_grid)
    idx = idx.clamp(1, x_grid.numel() - 1)
    x0 = x_grid[idx - 1]
    x1 = x_grid[idx]
    y0 = y_grid[idx - 1]
    y1 = y_grid[idx]
    t = (query - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)


class SmoothSlitChart:
    def __init__(self, a: float, n_grid: int = 20001):
        if a <= 1.0:
            raise ValueError("a must be > 1")
        self.a = float(a)
        self.beta = 2.0 * math.asin(1.0 / self.a)
        self.v0 = PI - 0.5 * self.beta
        self.v1 = self.v0 + self.beta
        self.eps = 1e-10

        y_grid = torch.linspace(0.0, 1.0, n_grid)
        phi = self.v0 + self.beta * y_grid
        t = -self.a * torch.cos(phi)
        disc = torch.sqrt(torch.clamp_min(1.0 - (self.a * torch.sin(phi)).square(), 0.0))
        rho_minus = t - disc
        rho_plus = t + disc
        width = 0.5 * self.beta * (rho_plus.square() - rho_minus.square())

        dy = 1.0 / (n_grid - 1)
        w_mid = 0.5 * (width[1:] + width[:-1])
        cumulative = torch.empty_like(y_grid)
        cumulative[0] = 0.0
        cumulative[1:] = torch.cumsum(w_mid, dim=0) * dy

        self.y_grid = y_grid
        self.width_grid = width
        self.W_grid = cumulative
        self.area_total = float(cumulative[-1])

    def _ray_data_from_phi(self, phi: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        t = -self.a * torch.cos(phi)
        disc = torch.sqrt(torch.clamp_min(1.0 - (self.a * torch.sin(phi)).square(), 0.0))
        rho_minus = t - disc
        rho_plus = t + disc
        width = 0.5 * self.beta * (rho_plus.square() - rho_minus.square())
        return rho_minus, rho_plus, width

    def forward(self, q: torch.Tensor, p: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        dx = q - self.a
        phi = torch.atan2(p, dx)
        phi = torch.where(phi < self.v0, phi + 2.0 * PI, phi)
        y = (phi - self.v0) / self.beta
        y = y.clamp(self.eps, 1.0 - self.eps)
        phi = self.v0 + self.beta * y
        rho = torch.sqrt(dx.square() + p.square())
        rho_minus, _, width_raw = self._ray_data_from_phi(phi)
        width = linear_interp(y, self.y_grid, self.width_grid)
        x_band = 0.5 * self.beta * (rho.square() - rho_minus.square())
        X = x_band / width
        Y = linear_interp(y, self.y_grid, self.W_grid)
        return X, Y

    def inverse(self, X: torch.Tensor, Y: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        Y = Y.clamp(self.eps, self.area_total - self.eps)
        y = linear_interp(Y, self.W_grid, self.y_grid)
        y = y.clamp(self.eps, 1.0 - self.eps)
        phi = self.v0 + self.beta * y
        rho_minus, _, _ = self._ray_data_from_phi(phi)
        width = linear_interp(y, self.y_grid, self.width_grid)
        x_band = X * width
        rho_sq = rho_minus.square() + (2.0 * x_band / self.beta)
        rho = torch.sqrt(torch.clamp_min(rho_sq, 0.0))
        q = self.a + rho * torch.cos(phi)
        p = rho * torch.sin(phi)
        return q, p

    def roundtrip_error(self, n_r: int = 121, n_theta: int = 241) -> float:
        rs = torch.linspace(0.0, 1.0, n_r)
        thetas = torch.linspace(-PI, PI, n_theta)
        rr, tt = torch.meshgrid(rs, thetas, indexing="ij")
        q = rr * torch.cos(tt)
        p = rr * torch.sin(tt)
        X, Y = self.forward(q.reshape(-1), p.reshape(-1))
        q2, p2 = self.inverse(X, Y)
        err = torch.sqrt((q2 - q.reshape(-1)).square() + (p2 - p.reshape(-1)).square())
        return float(err.max())


def translate_union_to_origin(points: torch.Tensor, k: int) -> tuple[torch.Tensor, torch.Tensor]:
    q = points[:, :2].clone()
    p = points[:, 2:].clone()
    indices = infer_ball_indices(q[:, 0], k)
    q[:, 0] -= 3.0 * indices.to(q.dtype)
    return torch.cat((q, p), dim=1), indices


def smooth_square_slot_map(points: torch.Tensor, k: int, chart: SmoothSlitChart) -> tuple[torch.Tensor, dict]:
    centered, indices = translate_union_to_origin(points, k)
    q = centered[:, :2]
    p = centered[:, 2:]

    X1, Y1 = chart.forward(q[:, 0], p[:, 0])
    X2, Y2 = chart.forward(q[:, 1], p[:, 1])

    l = math.ceil(math.sqrt(k))
    row = torch.div(indices, l, rounding_mode="floor").to(q.dtype)
    col = torch.remainder(indices, l).to(q.dtype)

    X1_out = float(l) * X1
    X2_out = float(l) * X2
    Y1_out = Y1 / float(l) + row * (PI / float(l))
    Y2_out = Y2 / float(l) + col * (PI / float(l))

    q1_out, p1_out = chart.inverse(X1_out, Y1_out)
    q2_out, p2_out = chart.inverse(X2_out, Y2_out)
    mapped = torch.stack((q1_out, q2_out, p1_out, p2_out), dim=1)
    return mapped, {
        "k": k,
        "l": l,
        "a": chart.a,
    }


def verify_candidate(k: int, a: float, n_pts_per_ball: int, seed: int) -> PackingSummary:
    chart = SmoothSlitChart(a=a)
    points = sample_boundary_points(k, n_pts_per_ball=n_pts_per_ball, seed=seed)
    mapped, meta = smooth_square_slot_map(points, k=k, chart=chart)
    sampled_radius = float(mapped.norm(dim=1).max())
    sampled_ratio = k / (sampled_radius**4)
    return PackingSummary(
        k=k,
        l=meta["l"],
        a=a,
        sampled_radius=sampled_radius,
        sampled_ratio=sampled_ratio,
        seed=seed,
        n_pts_per_ball=n_pts_per_ball,
        chart_roundtrip_max_error=chart.roundtrip_error(),
        chart_area_total=chart.area_total,
    )


def search_a_values(k: int, a_values: list[float], n_pts_per_ball: int, seed: int) -> dict:
    rows = []
    best = None
    for a in a_values:
        summary = verify_candidate(k=k, a=a, n_pts_per_ball=n_pts_per_ball, seed=seed)
        row = asdict(summary)
        rows.append(row)
        if best is None or row["sampled_ratio"] > best["sampled_ratio"]:
            best = row
    return {
        "k": k,
        "seed": seed,
        "n_pts_per_ball": n_pts_per_ball,
        "rows": rows,
        "best": best,
    }


def parse_a_values(a_min: float, a_max: float, num: int) -> list[float]:
    if num == 1:
        return [a_min]
    return [a_min + (a_max - a_min) * i / (num - 1) for i in range(num)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--a", type=float, default=None)
    parser.add_argument("--a-min", type=float, default=None)
    parser.add_argument("--a-max", type=float, default=None)
    parser.add_argument("--num-a", type=int, default=5)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--n-pts", type=int, default=2000)
    args = parser.parse_args()

    if args.a is not None:
        summary = verify_candidate(k=args.k, a=args.a, n_pts_per_ball=args.n_pts, seed=args.seed)
        print(json.dumps(asdict(summary), indent=2, sort_keys=True))
        return

    if args.a_min is None or args.a_max is None:
        parser.error("pass either --a or both --a-min and --a-max")
    result = search_a_values(
        k=args.k,
        a_values=parse_a_values(args.a_min, args.a_max, args.num_a),
        n_pts_per_ball=args.n_pts,
        seed=args.seed,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
