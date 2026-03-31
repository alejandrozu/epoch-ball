#!/usr/bin/env python3
"""Search source-side strip realizations for the attempt-010 hybrid family.

This tool studies the concrete `m=6` proof-of-concept that came out of
`hybrid_boundary_strip_family.py`. The target side-19 boundary region is a
36-triangle strip. The question here is whether the source triangle Delta(1)
admits a continuous piecewise-affine equal-area realization of the same
boundary-strip triangulation.

The current model is intentionally narrow:

* keep the exact target strip combinatorics;
* place all source vertices on the boundary of Delta(1);
* search over the only corner placements that avoid immediate collinearity;
* optimize the boundary spacings with PyTorch.

This does not prove correctness of the hybrid construction, but it gives
reproducible evidence about whether the `m=6` source-side strip is plausible
or immediately impossible.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import asdict, dataclass

import torch

from hybrid_boundary_strip_family import analyze_hybrid_family, build_adjacency, subdivide_triangle


torch.set_default_dtype(torch.float64)


@dataclass
class RealizationResult:
    reverse: bool
    corners: list[int]
    counts: list[int]
    side_of_arc: list[int]
    mse: float
    min_area: float
    max_area: float
    max_err: float
    mean_err: float
    all_positive: bool
    coords: list[list[float]]


@dataclass
class SearchSummary:
    m: int
    k: int
    target_face_area: float
    boundary_vertex_count: int
    face_count: int
    open_feasible_corner_triple_count: int
    closed_feasible_corner_triple_count: int
    tested_configs: int
    best_overall: RealizationResult
    best_positive: RealizationResult | None
    closed_feasible_triples: list[list[int]]
    target_boundary_cycle: list[list[int]]
    target_face_indices: list[list[int]]


def extract_strip_complex(m: int = 6, k: int = 10, reverse: bool = False) -> tuple[list[tuple[int, int]], list[tuple[int, int, int]]]:
    """Return the boundary cycle and face indices for the attempt-010 strip."""
    summary = analyze_hybrid_family(m=m, k=k)
    all_target = subdivide_triangle(summary.n)
    target_sorted = sorted(all_target, key=lambda piece: (piece.level, piece.orientation, piece.piece_id))
    boundary = target_sorted[9 * m * m : 10 * m * m]
    adj = build_adjacency(boundary)

    endpoints = [u for u, nbrs in adj.items() if len(nbrs) == 1]
    if len(endpoints) != 2:
        raise RuntimeError("expected a path strip with two endpoints")

    cur = endpoints[0]
    prev = None
    path = [cur]
    while True:
        nxts = [v for v in adj[cur] if v != prev]
        if not nxts:
            break
        nxt = nxts[0]
        path.append(nxt)
        prev, cur = cur, nxt

    pieces = {piece.piece_id: [tuple(v) for v in piece.lattice_vertices] for piece in boundary}

    edge_counts: dict[tuple[tuple[int, int], tuple[int, int]], int] = defaultdict(int)
    for pid in path:
        tri = pieces[pid]
        for idx in range(3):
            edge = tuple(sorted((tri[idx], tri[(idx + 1) % 3])))
            edge_counts[edge] += 1

    boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
    nb: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    for u, v in boundary_edges:
        nb[u].append(v)
        nb[v].append(u)

    start = min(nb)
    cycle = [start]
    prev = None
    cur = start
    while True:
        nxts = [v for v in nb[cur] if v != prev]
        if not nxts:
            break
        nxt = nxts[0] if len(nxts) == 1 or nxts[0] != start else nxts[1]
        if nxt == start:
            break
        cycle.append(nxt)
        prev, cur = cur, nxt

    index_of = {vertex: idx for idx, vertex in enumerate(cycle)}
    faces = [tuple(sorted(index_of[v] for v in pieces[pid])) for pid in path]
    if reverse:
        n_vertices = len(cycle)
        faces = [tuple(sorted((n_vertices - idx) % n_vertices for idx in face)) for face in faces]
    return cycle, faces


def corner_triple_counts(faces: list[tuple[int, int, int]], n_vertices: int) -> tuple[list[tuple[int, int, int]], list[tuple[int, int, int]]]:
    """Return corner placements under open-side and closed-side tests."""
    open_ok: list[tuple[int, int, int]] = []
    closed_ok: list[tuple[int, int, int]] = []

    for a in range(n_vertices):
        for b in range(a + 1, n_vertices):
            for c in range(b + 1, n_vertices):
                def side_open(idx: int) -> int:
                    if a <= idx < b:
                        return 0
                    if b <= idx < c:
                        return 1
                    return 2

                if all(len({side_open(idx) for idx in face}) > 1 for face in faces):
                    open_ok.append((a, b, c))

                side0 = set(range(a, b + 1))
                side1 = set(range(b, c + 1))
                side2 = set(list(range(c, n_vertices)) + list(range(0, a + 1)))
                if all(not (set(face) <= side0 or set(face) <= side1 or set(face) <= side2) for face in faces):
                    closed_ok.append((a, b, c))

    return open_ok, closed_ok


def boundary_arcs(corners: tuple[int, int, int], n_vertices: int) -> list[list[int]]:
    a, b, c = corners
    return [
        list(range(a, b + 1)),
        list(range(b, c + 1)),
        list(range(c, n_vertices)) + list(range(0, a + 1)),
    ]


def candidate_side_orders(corners: tuple[int, int, int], n_vertices: int) -> list[tuple[int, int, int]]:
    counts = [len(arc) - 1 for arc in boundary_arcs(corners, n_vertices)]
    long_arc = max(range(3), key=lambda idx: counts[idx])
    short_arcs = [idx for idx in range(3) if idx != long_arc]
    out = []
    for short_order in ((0, 2), (2, 0)):
        side_of_arc = [None, None, None]
        side_of_arc[long_arc] = 1
        side_of_arc[short_arcs[0]] = short_order[0]
        side_of_arc[short_arcs[1]] = short_order[1]
        out.append(tuple(side_of_arc))
    return out


def optimize_realization(
    faces: list[tuple[int, int, int]],
    n_vertices: int,
    corners: tuple[int, int, int],
    side_of_arc: tuple[int, int, int],
    steps: int,
    lr: float,
    target_area: float,
) -> RealizationResult:
    arcs = boundary_arcs(corners, n_vertices)
    face_tensor = torch.tensor(faces, dtype=torch.long)
    params = [torch.nn.Parameter(torch.zeros(len(arc) - 2)) for arc in arcs]
    opt = torch.optim.Adam(params, lr=lr)

    def build_coords() -> torch.Tensor:
        coords = torch.zeros(n_vertices, 2, dtype=torch.get_default_dtype())
        for arc_id, arc in enumerate(arcs):
            logits = torch.cat([params[arc_id], torch.zeros(1)])
            weights = torch.softmax(logits, dim=0)
            cumulative = torch.cat([torch.zeros(1), torch.cumsum(weights, dim=0)])
            side_type = side_of_arc[arc_id]
            for local, vid in enumerate(arc):
                t = cumulative[local]
                if side_type == 0:
                    coords[vid] = torch.stack([t, torch.tensor(0.0)])
                elif side_type == 1:
                    coords[vid] = torch.stack([1.0 - t, t])
                else:
                    coords[vid] = torch.stack([torch.tensor(0.0), 1.0 - t])
        return coords

    def face_areas(coords: torch.Tensor) -> torch.Tensor:
        pts = coords[face_tensor]
        v1 = pts[:, 1] - pts[:, 0]
        v2 = pts[:, 2] - pts[:, 0]
        return 0.5 * (v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0])

    best_snapshot: dict | None = None
    for step in range(steps):
        opt.zero_grad()
        coords = build_coords()
        areas = face_areas(coords)
        mse = ((areas - target_area) ** 2).mean()
        barrier = 50.0 * torch.relu(torch.tensor(1e-3) - areas).square().mean()
        entropy = torch.tensor(0.0)
        for param in params:
            weights = torch.softmax(torch.cat([param, torch.zeros(1)]), dim=0)
            entropy = entropy + 1e-4 * (1.0 / weights).mean()
        loss = mse + barrier + entropy
        loss.backward()
        opt.step()

        if step in (steps // 3, (2 * steps) // 3, steps - 1):
            with torch.no_grad():
                coords = build_coords()
                areas = face_areas(coords)
                errors = (areas - target_area).abs()
                snapshot = {
                    "mse": float(((areas - target_area) ** 2).mean()),
                    "min_area": float(areas.min()),
                    "max_area": float(areas.max()),
                    "max_err": float(errors.max()),
                    "mean_err": float(errors.mean()),
                    "coords": coords.tolist(),
                }
                if best_snapshot is None or snapshot["mse"] < best_snapshot["mse"]:
                    best_snapshot = snapshot

    if best_snapshot is None:
        raise RuntimeError("optimizer produced no snapshot")

    return RealizationResult(
        reverse=False,  # overwritten by caller when needed
        corners=list(corners),
        counts=[len(arc) - 1 for arc in arcs],
        side_of_arc=list(side_of_arc),
        mse=best_snapshot["mse"],
        min_area=best_snapshot["min_area"],
        max_area=best_snapshot["max_area"],
        max_err=best_snapshot["max_err"],
        mean_err=best_snapshot["mean_err"],
        all_positive=best_snapshot["min_area"] > 0.0,
        coords=best_snapshot["coords"],
    )


def search_realizations(m: int, k: int, config_limit: int, steps: int, lr: float) -> SearchSummary:
    boundary_cycle, faces_fwd = extract_strip_complex(m=m, k=k, reverse=False)
    _, faces_rev = extract_strip_complex(m=m, k=k, reverse=True)
    n_vertices = len(boundary_cycle)
    target_area = 1.0 / (2.0 * m * m)

    open_ok, closed_ok = corner_triple_counts(faces_fwd, n_vertices)
    closed_ok = sorted(closed_ok, key=lambda triple: (max(abs((triple[1] - triple[0]) - 19), abs((triple[2] - triple[1]) - 19)), max([
        triple[1] - triple[0],
        triple[2] - triple[1],
        n_vertices - (triple[2] - triple[0]),
    ])))

    results: list[RealizationResult] = []
    for reverse, faces in ((False, faces_fwd), (True, faces_rev)):
        for corners in closed_ok[:config_limit]:
            for side_of_arc in candidate_side_orders(corners, n_vertices):
                result = optimize_realization(
                    faces=faces,
                    n_vertices=n_vertices,
                    corners=corners,
                    side_of_arc=side_of_arc,
                    steps=steps,
                    lr=lr,
                    target_area=target_area,
                )
                result.reverse = reverse
                results.append(result)

    results.sort(key=lambda item: (item.mse, -item.min_area, item.max_err))
    best_overall = results[0]
    best_positive = next((item for item in results if item.all_positive), None)
    return SearchSummary(
        m=m,
        k=k,
        target_face_area=target_area,
        boundary_vertex_count=n_vertices,
        face_count=len(faces_fwd),
        open_feasible_corner_triple_count=len(open_ok),
        closed_feasible_corner_triple_count=len(closed_ok),
        tested_configs=len(results),
        best_overall=best_overall,
        best_positive=best_positive,
        closed_feasible_triples=[list(triple) for triple in closed_ok],
        target_boundary_cycle=[list(vertex) for vertex in boundary_cycle],
        target_face_indices=[list(face) for face in faces_fwd],
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--m", type=int, default=6)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--config-limit", type=int, default=6)
    parser.add_argument("--steps", type=int, default=3000)
    parser.add_argument("--lr", type=float, default=0.02)
    args = parser.parse_args()

    if args.m != 6 or args.k != 10:
        raise ValueError("this analyzer currently targets the saved attempt-010 m=6, k=10 proof-of-concept")

    summary = search_realizations(
        m=args.m,
        k=args.k,
        config_limit=args.config_limit,
        steps=args.steps,
        lr=args.lr,
    )
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
