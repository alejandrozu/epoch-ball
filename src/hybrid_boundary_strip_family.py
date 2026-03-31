#!/usr/bin/env python3
"""Analyze a hybrid near-full target decomposition for k=10.

This family is qualitatively different from the previous single-shape search.
Instead of asking one connected target region to do all the work, it splits the
near-full target set into:

1. an inner side-(3m) triangle, which can be tiled exactly by 9 rigid side-m
   triangles; and
2. a remaining boundary region of exactly m^2 tiny triangles, intended for the
   tenth ball via a non-simplicial strip complex.

The script makes that decomposition explicit on the triangular grid and audits
the boundary region combinatorially: size, connectivity, degree profile, path
status, and simple Hamiltonian-path obstructions coming from degree-1 vertices.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict, deque
from dataclasses import asdict, dataclass


@dataclass
class TinyTriangle:
    piece_id: int
    orientation: str
    lattice_vertices: list[list[int]]
    level: int


@dataclass
class HybridSummary:
    k: int
    m: int
    n: int
    target_mu: float
    exact_ratio: float
    used_target_pieces: int
    discarded_target_pieces: int
    inner_piece_count: int
    inner_max_level: int
    inner_macro_triangle_count: int
    boundary_piece_count: int
    boundary_connected: bool
    boundary_degree_histogram: dict[str, int]
    boundary_degree1_count: int
    boundary_degree3_count: int
    boundary_is_path_graph: bool
    boundary_hamiltonian_path_degree_obstruction: bool
    boundary_level_histogram: dict[str, int]
    source_up_triangle_count: int
    source_down_triangle_count: int
    source_rhombus_tiling_possible_by_parity: bool
    boundary_path_head: list[int]
    boundary_path_tail: list[int]
    inner_macro_triangles: list[dict]


def subdivide_triangle(side: int) -> list[TinyTriangle]:
    pieces: list[TinyTriangle] = []
    piece_id = 0
    for i in range(side):
        for j in range(side - i):
            pieces.append(
                TinyTriangle(
                    piece_id=piece_id,
                    orientation="up",
                    lattice_vertices=[[i, j], [i + 1, j], [i, j + 1]],
                    level=i + j + 1,
                )
            )
            piece_id += 1
    for i in range(side - 1):
        for j in range(side - 1 - i):
            pieces.append(
                TinyTriangle(
                    piece_id=piece_id,
                    orientation="down",
                    lattice_vertices=[[i + 1, j], [i + 1, j + 1], [i, j + 1]],
                    level=i + j + 2,
                )
            )
            piece_id += 1
    return pieces


def build_adjacency(pieces: list[TinyTriangle]) -> dict[int, set[int]]:
    edge_to_piece: dict[tuple[tuple[int, int], tuple[int, int]], list[int]] = defaultdict(list)
    for piece in pieces:
        verts = [tuple(v) for v in piece.lattice_vertices]
        for idx in range(3):
            u = verts[idx]
            v = verts[(idx + 1) % 3]
            edge_to_piece[tuple(sorted((u, v)))].append(piece.piece_id)

    adj = {piece.piece_id: set() for piece in pieces}
    for ids in edge_to_piece.values():
        if len(ids) == 2:
            u, v = ids
            adj[u].add(v)
            adj[v].add(u)
    return adj


def connected_component_size(adj: dict[int, set[int]], start: int) -> int:
    seen = {start}
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                queue.append(v)
    return len(seen)


def path_order_if_path(adj: dict[int, set[int]]) -> list[int]:
    endpoints = [u for u, nbrs in adj.items() if len(nbrs) == 1]
    if len(endpoints) != 2:
        return []
    order = [endpoints[0]]
    prev = None
    cur = endpoints[0]
    while True:
        nxts = [v for v in adj[cur] if v != prev]
        if not nxts:
            break
        nxt = nxts[0]
        order.append(nxt)
        prev, cur = cur, nxt
    if len(order) != len(adj):
        return []
    return order


def coarse_inner_macro_triangles(m: int) -> list[dict]:
    coarse = subdivide_triangle(side=3)
    out = []
    for piece in coarse:
        out.append(
            {
                "orientation": piece.orientation,
                "lattice_vertices": [[m * a, m * b] for a, b in piece.lattice_vertices],
                "piece_id": piece.piece_id,
            }
        )
    return out


def analyze_hybrid_family(m: int, k: int = 10) -> HybridSummary:
    if k != 10:
        raise ValueError("this hybrid family analyzer currently implements only the k=10 split as 9+1")

    n = math.ceil(m * math.sqrt(k))
    all_target = subdivide_triangle(side=n)
    target_sorted = sorted(all_target, key=lambda piece: (piece.level, piece.orientation, piece.piece_id))
    used_target = target_sorted[: k * m * m]
    inner_target = target_sorted[: 9 * m * m]
    boundary_target = target_sorted[9 * m * m : 10 * m * m]

    boundary_adj = build_adjacency(boundary_target)
    boundary_ids = [piece.piece_id for piece in boundary_target]
    connected = connected_component_size(boundary_adj, boundary_ids[0]) == len(boundary_ids) if boundary_ids else False
    boundary_deg_hist = Counter(len(boundary_adj[u]) for u in boundary_adj)
    boundary_levels = Counter(piece.level for piece in boundary_target)
    boundary_is_path = connected and boundary_deg_hist.get(1, 0) == 2 and boundary_deg_hist.get(3, 0) == 0 and boundary_deg_hist.get(0, 0) == 0
    boundary_order = path_order_if_path(boundary_adj) if boundary_is_path else []

    source = subdivide_triangle(side=m)
    source_orient = Counter(piece.orientation for piece in source)

    return HybridSummary(
        k=k,
        m=m,
        n=n,
        target_mu=n / float(m),
        exact_ratio=(k * m * m) / float(n * n),
        used_target_pieces=len(used_target),
        discarded_target_pieces=len(all_target) - len(used_target),
        inner_piece_count=len(inner_target),
        inner_max_level=max(piece.level for piece in inner_target),
        inner_macro_triangle_count=len(coarse_inner_macro_triangles(m)),
        boundary_piece_count=len(boundary_target),
        boundary_connected=connected,
        boundary_degree_histogram={str(key): value for key, value in sorted(boundary_deg_hist.items())},
        boundary_degree1_count=boundary_deg_hist.get(1, 0),
        boundary_degree3_count=boundary_deg_hist.get(3, 0),
        boundary_is_path_graph=boundary_is_path,
        boundary_hamiltonian_path_degree_obstruction=boundary_deg_hist.get(1, 0) > 2,
        boundary_level_histogram={str(key): value for key, value in sorted(boundary_levels.items())},
        source_up_triangle_count=source_orient["up"],
        source_down_triangle_count=source_orient["down"],
        source_rhombus_tiling_possible_by_parity=(source_orient["up"] == source_orient["down"]),
        boundary_path_head=boundary_order[:12],
        boundary_path_tail=boundary_order[-12:] if boundary_order else [],
        inner_macro_triangles=coarse_inner_macro_triangles(m),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--m", type=int, nargs="+", required=True)
    parser.add_argument("--k", type=int, default=10)
    args = parser.parse_args()

    summaries = [asdict(analyze_hybrid_family(m=value, k=args.k)) for value in args.m]
    if len(summaries) == 1:
        print(json.dumps(summaries[0], indent=2, sort_keys=True))
    else:
        print(json.dumps({"rows": summaries}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
