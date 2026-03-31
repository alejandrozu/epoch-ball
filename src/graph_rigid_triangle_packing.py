#!/usr/bin/env python3
"""Graph search for continuous simplicial rigid-triangle packings.

This attempt treats the continuity-preserving simplicial origami family as a
graph problem on the triangular lattice.

Key local fact:
    once a source unit triangle is fixed, continuity and unit-triangle
    compatibility force every later row of the side-m triangulation.

So each continuous simplicial embedding of the m-subdivided source triangle is
just a rigid lattice triangle of side m. The remaining combinatorial question is
whether k such rigid side-m triangles can be packed into the target side-n
triangle, where n = ceil(m * sqrt(k)).
"""

from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass

import networkx as nx

Point = tuple[int, int]
TriangleId = tuple[Point, Point, Point]

DIRS: tuple[Point, ...] = (
    (1, 0),
    (0, 1),
    (-1, 1),
    (-1, 0),
    (0, -1),
    (1, -1),
)


@dataclass
class RigiditySummary:
    m_max: int
    admissible_seed_count: int
    all_rows_affine: bool
    all_rows_unique: bool
    sample_last_row: list[list[int]]


@dataclass
class PackingSummary:
    k: int
    m: int
    n: int
    exact_ratio_if_perfect: float
    placement_count: int
    best_found_count: int
    best_found_ratio: float
    best_found_unit_triangles: int
    exact_max_count: int | None
    exact_max_ratio: float | None
    exact_target_feasible: bool | None
    exact_target_ratio: float | None
    exact_search_used: bool
    random_trials: int
    seed: int


def add(u: Point, v: Point) -> Point:
    return (u[0] + v[0], u[1] + v[1])


def sub(u: Point, v: Point) -> Point:
    return (u[0] - v[0], u[1] - v[1])


def neighbors(p: Point) -> set[Point]:
    return {add(p, d) for d in DIRS}


def common_neighbors(u: Point, v: Point) -> set[Point]:
    return neighbors(u) & neighbors(v)


def is_unit_edge(u: Point, v: Point) -> bool:
    return sub(v, u) in DIRS or sub(u, v) in DIRS


def in_target_triangle(p: Point, n: int) -> bool:
    return p[0] >= 0 and p[1] >= 0 and p[0] + p[1] <= n


def canonical_triangle(vertices: tuple[Point, Point, Point]) -> TriangleId:
    return tuple(sorted(vertices))


def all_seed_rows() -> list[list[Point]]:
    origin = (0, 0)
    out = []
    neigh = list(neighbors(origin))
    for u in neigh:
        for v in neigh:
            if u == v:
                continue
            if is_unit_edge(u, v):
                out.append([u, v])
    return out


def other_common_neighbor(u: Point, v: Point, exclude: Point) -> Point:
    candidates = common_neighbors(u, v) - {exclude}
    if len(candidates) != 1:
        raise ValueError(f"expected one common neighbor, got {candidates} for edge {(u, v)}")
    return next(iter(candidates))


def propagate_rows(seed_row1: list[Point], m: int, origin: Point = (0, 0)) -> list[list[Point]]:
    if m < 1:
        raise ValueError("m must be positive")
    if len(seed_row1) != 2:
        raise ValueError("seed_row1 must have length 2")
    if not is_unit_edge(origin, seed_row1[0]) or not is_unit_edge(origin, seed_row1[1]) or not is_unit_edge(seed_row1[0], seed_row1[1]):
        raise ValueError("seed_row1 must form a unit triangle with origin")

    rows = [[origin], list(seed_row1)]
    prev = rows[0]
    curr = rows[1]

    for r in range(1, m):
        nxt: list[Point | None] = [None] * (r + 2)
        for i in range(1, r + 1):
            nxt[i] = other_common_neighbor(curr[i - 1], curr[i], prev[i - 1])

        nxt[0] = other_common_neighbor(curr[0], nxt[1], curr[1]) if r >= 1 else None
        nxt[r + 1] = other_common_neighbor(nxt[r], curr[r], curr[r - 1]) if r >= 1 else None
        rows.append([p for p in nxt if p is not None])
        prev = curr
        curr = rows[-1]

    return rows[: m + 1]


def rows_are_affine(rows: list[list[Point]]) -> bool:
    origin = rows[0][0]
    b = sub(rows[1][0], origin)
    a = sub(rows[1][1], origin)
    for r, row in enumerate(rows):
        if len(row) != r + 1:
            return False
        for i, p in enumerate(row):
            expected = add(origin, add((a[0] * i, a[1] * i), (b[0] * (r - i), b[1] * (r - i))))
            if p != expected:
                return False
    return True


def rows_have_unique_vertices(rows: list[list[Point]]) -> bool:
    seen: set[Point] = set()
    for row in rows:
        for p in row:
            if p in seen:
                return False
            seen.add(p)
    return True


def placement_triangles(rows: list[list[Point]]) -> list[TriangleId]:
    m = len(rows) - 1
    out: list[TriangleId] = []
    for r in range(m):
        for i in range(r + 1):
            up = canonical_triangle((rows[r][i], rows[r + 1][i + 1], rows[r + 1][i]))
            out.append(up)
        if r >= 1:
            for i in range(r):
                down = canonical_triangle((rows[r][i + 1], rows[r + 1][i + 1], rows[r][i]))
                out.append(down)
    return out


def placement_from_seed(seed0: Point, seed_row1: list[Point], m: int, n: int) -> dict | None:
    translated_row1 = [add(seed0, p) for p in seed_row1]
    rows = propagate_rows(translated_row1, m=m, origin=seed0)
    if any(not in_target_triangle(p, n) for row in rows for p in row):
        return None
    tris = placement_triangles(rows)
    if len(tris) != m * m or len(set(tris)) != m * m:
        return None
    return {
        "rows": rows,
        "triangles": frozenset(tris),
    }


def enumerate_placements(m: int, n: int) -> list[dict]:
    seeds = all_seed_rows()
    vertices = [(i, j) for i in range(n + 1) for j in range(n + 1 - i)]
    placements_by_triangles: dict[frozenset[TriangleId], dict] = {}

    for seed0 in vertices:
        for seed_row1 in seeds:
            placement = placement_from_seed(seed0=seed0, seed_row1=seed_row1, m=m, n=n)
            if placement is None:
                continue
            placements_by_triangles.setdefault(placement["triangles"], placement)

    placements = list(placements_by_triangles.values())
    placements.sort(key=lambda item: sorted(item["triangles"]))
    return placements


def greedy_random_packing(placements: list[dict], k: int, n_trials: int, seed: int) -> tuple[int, list[int]]:
    rng = random.Random(seed)
    triangle_sets = [set(item["triangles"]) for item in placements]
    usage = {}
    for idx, tris in enumerate(triangle_sets):
        for tri in tris:
            usage.setdefault(tri, []).append(idx)

    best_count = 0
    best_indices: list[int] = []

    scores = [sum(len(usage[tri]) for tri in tris) for tris in triangle_sets]

    for _ in range(n_trials):
        order = list(range(len(placements)))
        rng.shuffle(order)
        order.sort(key=lambda idx: (scores[idx], rng.random()))

        chosen: list[int] = []
        occupied: set[TriangleId] = set()
        for idx in order:
            if triangle_sets[idx] & occupied:
                continue
            chosen.append(idx)
            occupied.update(triangle_sets[idx])
            if len(chosen) == k:
                return len(chosen), chosen

        if len(chosen) > best_count:
            best_count = len(chosen)
            best_indices = list(chosen)

    return best_count, best_indices


def compatibility_graph(placements: list[dict]) -> nx.Graph:
    g = nx.Graph()
    g.add_nodes_from(range(len(placements)))
    triangle_sets = [set(item["triangles"]) for item in placements]
    for i in range(len(placements)):
        ti = triangle_sets[i]
        for j in range(i + 1, len(placements)):
            if not (ti & triangle_sets[j]):
                g.add_edge(i, j)
    return g


def exact_max_packing_count(placements: list[dict], target_k: int | None = None) -> int:
    g = compatibility_graph(placements)
    best = 0
    for clique in nx.find_cliques(g):
        if len(clique) > best:
            best = len(clique)
        if target_k is not None and best >= target_k:
            return best
    return best


def rigid_triangle_ratio(k: int, m: int, count: int) -> tuple[int, float]:
    n = math.ceil(m * math.sqrt(k))
    covered = count * m * m
    ratio = covered / float(n * n)
    return n, ratio


def verify_rigidity(m_max: int) -> RigiditySummary:
    seeds = all_seed_rows()
    last_rows = propagate_rows(seeds[0], m=m_max)
    all_affine = True
    all_unique = True
    for seed in seeds:
        for m in range(1, m_max + 1):
            rows = propagate_rows(seed, m=m)
            all_affine &= rows_are_affine(rows)
            all_unique &= rows_have_unique_vertices(rows)

    return RigiditySummary(
        m_max=m_max,
        admissible_seed_count=len(seeds),
        all_rows_affine=all_affine,
        all_rows_unique=all_unique,
        sample_last_row=[list(p) for p in last_rows[-1]],
    )


def search_rigid_packing(k: int, m: int, n_trials: int, seed: int, exact_limit: int) -> PackingSummary:
    n = math.ceil(m * math.sqrt(k))
    placements = enumerate_placements(m=m, n=n)
    best_count, _ = greedy_random_packing(placements=placements, k=k, n_trials=n_trials, seed=seed)
    exact_used = len(placements) <= exact_limit
    exact_max = exact_max_packing_count(placements, target_k=k) if exact_used else None
    return PackingSummary(
        k=k,
        m=m,
        n=n,
        exact_ratio_if_perfect=(k * m * m) / float(n * n),
        placement_count=len(placements),
        best_found_count=best_count,
        best_found_ratio=(best_count * m * m) / float(n * n),
        best_found_unit_triangles=best_count * m * m,
        exact_max_count=exact_max,
        exact_max_ratio=((exact_max * m * m) / float(n * n)) if exact_max is not None else None,
        exact_target_feasible=(exact_max >= k) if exact_max is not None else None,
        exact_target_ratio=((min(exact_max, k) * m * m) / float(n * n)) if exact_max is not None else None,
        exact_search_used=exact_used,
        random_trials=n_trials,
        seed=seed,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--m", type=int, default=None)
    parser.add_argument("--m-max", type=int, default=None)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--trials", type=int, default=2000)
    parser.add_argument("--exact-limit", type=int, default=300)
    parser.add_argument("--mode", choices=("rigidity", "packing", "sweep"), default="sweep")
    args = parser.parse_args()

    if args.mode == "rigidity":
        if args.m_max is None:
            parser.error("--m-max is required in rigidity mode")
        print(json.dumps(asdict(verify_rigidity(args.m_max)), indent=2, sort_keys=True))
        return

    if args.mode == "packing":
        if args.m is None:
            parser.error("--m is required in packing mode")
        print(
            json.dumps(
                asdict(search_rigid_packing(k=args.k, m=args.m, n_trials=args.trials, seed=args.seed, exact_limit=args.exact_limit)),
                indent=2,
                sort_keys=True,
            )
        )
        return

    if args.m_max is None:
        parser.error("--m-max is required in sweep mode")

    payload = {
        "rigidity": asdict(verify_rigidity(args.m_max)),
        "packing": [
            asdict(search_rigid_packing(k=args.k, m=m, n_trials=args.trials, seed=args.seed, exact_limit=args.exact_limit))
            for m in range(1, args.m_max + 1)
        ],
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
