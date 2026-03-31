#!/usr/bin/env python3
"""Piecewise toric origami tilings in moment coordinates.

This attempt uses a new topological/origami family. Each source moment
triangle Delta(1) is subdivided into an m-by-m lattice triangulation with
m^2 tiny triangles. Those tiny triangles are reassigned to tiny target
triangles inside Delta(n / m), where n = ceil(m * sqrt(k)).

On each tiny triangle the map is an affine symplectic toric chart
    x' = A x + b,   y' = A^{-T} y,
with A determined by the source and target tiny triangles.

The family achieves a near-full combinatorial ratio k m^2 / n^2, but the
piecewise map is generally discontinuous across fold lines. This script
quantifies both the density and the continuity failure.
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
class TinyTriangle:
    piece_id: int
    orientation: str
    lattice_vertices: list[list[int]]
    real_vertices: list[list[float]]
    level: int


@dataclass
class OrigamiSummary:
    k: int
    m: int
    n: int
    target_mu: float
    combinatorial_ratio: float
    exact_radius: float
    exact_ratio: float
    sampled_radius: float
    sampled_ratio: float
    seed: int
    n_pts_per_ball: int
    pieces_per_ball: int
    used_target_pieces: int
    source_adjacencies: int
    continuous_adjacencies_ball0: int
    max_edge_midpoint_jump_ball0: float
    mean_edge_midpoint_jump_ball0: float


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
    y = torch.stack((y1, y2), dim=1)
    x = torch.stack((x1, x2), dim=1)
    return x, y


def action_angle_to_cartesian(x: torch.Tensor, y: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    r = torch.sqrt(torch.clamp_min(x, 0.0))
    angles = TWO_PI * y
    q = r * torch.cos(angles)
    p = r * torch.sin(angles)
    return q, p


def subdivide_triangle(side: int, denominator: int) -> list[TinyTriangle]:
    pieces: list[TinyTriangle] = []
    piece_id = 0
    scale = float(denominator)

    for i in range(side):
        for j in range(side - i):
            verts = [[i, j], [i + 1, j], [i, j + 1]]
            pieces.append(
                TinyTriangle(
                    piece_id=piece_id,
                    orientation="up",
                    lattice_vertices=verts,
                    real_vertices=[[a / scale, b / scale] for a, b in verts],
                    level=i + j + 1,
                )
            )
            piece_id += 1

    for i in range(side - 1):
        for j in range(side - 1 - i):
            verts = [[i + 1, j], [i + 1, j + 1], [i, j + 1]]
            pieces.append(
                TinyTriangle(
                    piece_id=piece_id,
                    orientation="down",
                    lattice_vertices=verts,
                    real_vertices=[[a / scale, b / scale] for a, b in verts],
                    level=i + j + 2,
                )
            )
            piece_id += 1

    return pieces


def classify_source_piece_ids(x: torch.Tensor, m: int, up_grid: torch.Tensor, down_grid: torch.Tensor) -> torch.Tensor:
    u = x * float(m)
    i = torch.floor(u[:, 0]).to(torch.int64).clamp(0, m - 1)
    j = torch.floor(u[:, 1]).to(torch.int64).clamp(0, m - 1)
    du = u[:, 0] - i.to(u.dtype)
    dv = u[:, 1] - j.to(u.dtype)
    is_up = (i + j >= m - 1) | (du + dv <= 1.0)
    return torch.where(is_up, up_grid[i, j], down_grid[i, j])


def piece_id_grids(m: int, pieces: list[TinyTriangle]) -> tuple[torch.Tensor, torch.Tensor]:
    up_grid = torch.full((m, m), -1, dtype=torch.int64)
    down_grid = torch.full((m, m), -1, dtype=torch.int64)
    for piece in pieces:
        a, b = piece.lattice_vertices[0]
        if piece.orientation == "up":
            up_grid[a, b] = piece.piece_id
        else:
            # Down triangles are identified by the lower-left square cell.
            down_grid[a - 1, b] = piece.piece_id
    return up_grid, down_grid


def affine_map_from_triangles(source: TinyTriangle, target: TinyTriangle) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    s = torch.tensor(source.real_vertices, dtype=torch.get_default_dtype())
    t = torch.tensor(target.real_vertices, dtype=torch.get_default_dtype())
    s0 = s[0]
    t0 = t[0]
    s_mat = torch.stack((s[1] - s0, s[2] - s0), dim=1)
    t_mat = torch.stack((t[1] - t0, t[2] - t0), dim=1)
    A = t_mat @ torch.linalg.inv(s_mat)
    b = t0 - A @ s0
    A_inv_t = torch.linalg.inv(A).transpose(0, 1)
    return A, b, A_inv_t


def source_adjacencies(pieces: list[TinyTriangle]) -> list[tuple[int, int, list[list[float]]]]:
    edge_to_piece: dict[tuple[tuple[int, int], tuple[int, int]], list[int]] = {}
    for piece in pieces:
        verts = piece.lattice_vertices
        for idx in range(3):
            u = tuple(verts[idx])
            v = tuple(verts[(idx + 1) % 3])
            key = tuple(sorted((u, v)))
            edge_to_piece.setdefault(key, []).append(piece.piece_id)

    out = []
    for (u, v), ids in edge_to_piece.items():
        if len(ids) == 2:
            out.append((ids[0], ids[1], [[u[0], u[1]], [v[0], v[1]]]))
    return out


def min_m_for_ratio(k: int, target_ratio: float, m_cap: int = 10000) -> int:
    sqrt_k = math.sqrt(k)
    for m in range(1, m_cap + 1):
        n = math.ceil(m * sqrt_k)
        ratio = (k * m * m) / float(n * n)
        if ratio > target_ratio:
            return m
    raise RuntimeError(f"failed to find m <= {m_cap} for target ratio {target_ratio}")


def build_assignment(k: int, m: int) -> dict:
    source_pieces = subdivide_triangle(side=m, denominator=m)
    up_grid, down_grid = piece_id_grids(m, source_pieces)

    n = math.ceil(m * math.sqrt(k))
    target_pieces = subdivide_triangle(side=n, denominator=m)
    target_pieces_sorted = sorted(target_pieces, key=lambda piece: (piece.level, piece.orientation, piece.piece_id))
    used_targets = target_pieces_sorted[: k * m * m]

    A_list = []
    b_list = []
    A_inv_t_list = []
    max_level = 0

    for ball_idx in range(k):
        for local_piece_id, source_piece in enumerate(source_pieces):
            target_piece = used_targets[ball_idx * (m * m) + local_piece_id]
            A, b, A_inv_t = affine_map_from_triangles(source_piece, target_piece)
            A_list.append(A)
            b_list.append(b)
            A_inv_t_list.append(A_inv_t)
            max_level = max(max_level, target_piece.level)

    return {
        "k": k,
        "m": m,
        "n": n,
        "source_pieces": source_pieces,
        "up_grid": up_grid,
        "down_grid": down_grid,
        "A": torch.stack(A_list, dim=0),
        "b": torch.stack(b_list, dim=0),
        "A_inv_t": torch.stack(A_inv_t_list, dim=0),
        "used_targets": used_targets,
        "source_adjacencies": source_adjacencies(source_pieces),
        "max_level": max_level,
    }


def translate_union_to_origin(points: torch.Tensor, k: int) -> tuple[torch.Tensor, torch.Tensor]:
    q = points[:, :2].clone()
    p = points[:, 2:].clone()
    indices = infer_ball_indices(q[:, 0], k)
    q[:, 0] -= 3.0 * indices.to(q.dtype)
    return torch.cat((q, p), dim=1), indices


def origami_map(points: torch.Tensor, assignment: dict) -> torch.Tensor:
    k = assignment["k"]
    m = assignment["m"]
    centered, ball_indices = translate_union_to_origin(points, k)
    q = centered[:, :2]
    p = centered[:, 2:]
    x, y = cartesian_to_action_angle(q, p)

    local_piece_ids = classify_source_piece_ids(x, m, assignment["up_grid"], assignment["down_grid"])
    global_piece_ids = ball_indices * (m * m) + local_piece_ids

    A = assignment["A"][global_piece_ids]
    b = assignment["b"][global_piece_ids]
    A_inv_t = assignment["A_inv_t"][global_piece_ids]

    x_out = torch.einsum("nij,nj->ni", A, x) + b
    y_out = torch.einsum("nij,nj->ni", A_inv_t, y)
    q_out, p_out = action_angle_to_cartesian(x_out, y_out)
    return torch.cat((q_out, p_out), dim=1)


def continuity_audit_ball0(assignment: dict) -> tuple[int, int, float, float]:
    m = assignment["m"]
    jumps = []
    for u_id, v_id, edge in assignment["source_adjacencies"]:
        s0 = torch.tensor(edge[0], dtype=torch.get_default_dtype()) / float(m)
        s1 = torch.tensor(edge[1], dtype=torch.get_default_dtype()) / float(m)
        midpoint = 0.5 * (s0 + s1)
        Au = assignment["A"][u_id]
        bu = assignment["b"][u_id]
        Av = assignment["A"][v_id]
        bv = assignment["b"][v_id]
        jump = float(torch.norm(Au @ midpoint + bu - (Av @ midpoint + bv)))
        jumps.append(jump)

    continuous = sum(jump <= 1e-12 for jump in jumps)
    max_jump = max(jumps) if jumps else 0.0
    mean_jump = sum(jumps) / len(jumps) if jumps else 0.0
    return len(jumps), continuous, max_jump, mean_jump


def verify_origami_candidate(k: int, m: int, n_pts_per_ball: int, seed: int) -> OrigamiSummary:
    assignment = build_assignment(k=k, m=m)
    target_mu = assignment["n"] / float(m)
    combinatorial_ratio = (k * m * m) / float(assignment["n"] * assignment["n"])

    exact_radius = math.sqrt(assignment["max_level"] / float(m))
    exact_ratio = k / (exact_radius**4)

    points = sample_boundary_points(k, n_pts_per_ball=n_pts_per_ball, seed=seed)
    mapped = origami_map(points, assignment)
    sampled_radius = float(mapped.norm(dim=1).max())
    sampled_ratio = k / (sampled_radius**4)

    adjacency_count, continuous_count, max_jump, mean_jump = continuity_audit_ball0(assignment)

    return OrigamiSummary(
        k=k,
        m=m,
        n=assignment["n"],
        target_mu=target_mu,
        combinatorial_ratio=combinatorial_ratio,
        exact_radius=exact_radius,
        exact_ratio=exact_ratio,
        sampled_radius=sampled_radius,
        sampled_ratio=sampled_ratio,
        seed=seed,
        n_pts_per_ball=n_pts_per_ball,
        pieces_per_ball=m * m,
        used_target_pieces=k * m * m,
        source_adjacencies=adjacency_count,
        continuous_adjacencies_ball0=continuous_count,
        max_edge_midpoint_jump_ball0=max_jump,
        mean_edge_midpoint_jump_ball0=mean_jump,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--eps", type=float, default=None)
    parser.add_argument("--m", type=int, default=None)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--n-pts", type=int, default=500)
    args = parser.parse_args()

    if args.m is None:
        if args.eps is None:
            parser.error("pass either --m or --eps")
        m = min_m_for_ratio(args.k, 1.0 - args.eps)
    else:
        m = args.m

    summary = verify_origami_candidate(k=args.k, m=m, n_pts_per_ball=args.n_pts, seed=args.seed)
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
