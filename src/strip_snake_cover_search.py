#!/usr/bin/env python3
"""Probabilistic search over non-simplicial strip-snake toric cover shapes.

This family is intentionally different from the previous simplicial origami
attempts. A single ball is modeled by a connected snake polyomino of area 1/2
in action space; the shape is then pushed through an integer linear toric cover
matrix A with det(A) >= k. Distinct source balls occupy different angle slots
of the cover, so the achieved packing ratio is

    ratio = k / mu^2,

where mu is the minimal target capacity for the transformed shape.

The search is only a surrogate for the eventual smooth Hamiltonian, but any
simple connected polygonal snake shape can in principle be realized by a
piecewise-affine strip-fold homeomorphism of Delta(1), so the family is a
concrete non-simplicial continuation of the origami line.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass


Cell = tuple[int, int]
Point = tuple[float, float]


@dataclass
class MatrixCandidate:
    determinant: int
    matrix: list[list[int]]


@dataclass
class SnakeCandidate:
    n_cells: int
    determinant: int
    matrix: list[list[int]]
    required_capacity_mu: float
    achieved_ratio: float
    scale: float
    cells: list[list[int]]
    transformed_corner_count: int


@dataclass
class SearchSummary:
    k: int
    det_min: int
    det_max: int
    entry_bound: int
    shape_trials: int
    seed: int
    exact_enumeration: bool
    best_overall: SnakeCandidate
    best_by_n_cells: list[SnakeCandidate]
    matrix_count: int
    unique_shape_count_by_n_cells: list[dict[str, int]]


def normalize_cells(cells: set[Cell]) -> tuple[Cell, ...]:
    min_x = min(x for x, _ in cells)
    min_y = min(y for _, y in cells)
    shifted = sorted((x - min_x, y - min_y) for x, y in cells)
    return tuple(shifted)


def cell_neighbors(cell: Cell) -> list[Cell]:
    x, y = cell
    return [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]


def cell_adjacency_count(cells: tuple[Cell, ...]) -> tuple[int, dict[Cell, int]]:
    cell_set = set(cells)
    degrees = {cell: 0 for cell in cells}
    adjacency_edges = 0
    for cell in cells:
        x, y = cell
        for nb in ((x + 1, y), (x, y + 1)):
            if nb in cell_set:
                degrees[cell] += 1
                degrees[nb] += 1
                adjacency_edges += 1
    return adjacency_edges, degrees


def is_path_polyomino(cells: tuple[Cell, ...]) -> bool:
    if len(cells) <= 2:
        return True
    adjacency_edges, degrees = cell_adjacency_count(cells)
    if adjacency_edges != len(cells) - 1:
        return False
    degree_values = sorted(degrees.values())
    return degree_values.count(1) == 2 and all(value <= 2 for value in degree_values)


def random_snake_polyomino(n_cells: int, rng: random.Random, max_restarts: int = 200) -> tuple[Cell, ...]:
    if n_cells < 1:
        raise ValueError("n_cells must be positive")
    if n_cells == 1:
        return ((0, 0),)

    for _ in range(max_restarts):
        path = [(0, 0)]
        used = {path[0]}
        while len(path) < n_cells:
            frontier = [nb for nb in cell_neighbors(path[-1]) if nb not in used]
            if not frontier:
                break
            nxt = rng.choice(frontier)
            path.append(nxt)
            used.add(nxt)
        if len(path) == n_cells:
            return normalize_cells(set(path))

    # Fallback to a straight snake when random growth fails.
    return tuple((i, 0) for i in range(n_cells))


def canonical_seed_shapes(n_cells: int) -> list[tuple[Cell, ...]]:
    out: list[tuple[Cell, ...]] = []
    out.append(tuple((i, 0) for i in range(n_cells)))

    staircase: list[Cell] = [(0, 0)]
    x = 0
    y = 0
    for step in range(1, n_cells):
        if step % 2 == 1:
            x += 1
        else:
            y += 1
        staircase.append((x, y))
    out.append(normalize_cells(set(staircase)))
    return list(dict.fromkeys(out))


def enumerate_snake_polyominoes(n_cells: int) -> list[tuple[Cell, ...]]:
    """Enumerate all translation-normalized path-polyomino shapes with n cells."""
    if n_cells < 1:
        raise ValueError("n_cells must be positive")

    shapes: set[tuple[Cell, ...]] = set()
    path = [(0, 0)]
    used = {path[0]}

    def dfs() -> None:
        if len(path) == n_cells:
            shapes.add(normalize_cells(set(path)))
            return
        for nxt in cell_neighbors(path[-1]):
            if nxt in used:
                continue
            used.add(nxt)
            path.append(nxt)
            dfs()
            path.pop()
            used.remove(nxt)

    dfs()
    return sorted(shapes)


def shape_corners(cells: tuple[Cell, ...]) -> list[Point]:
    corners: set[tuple[int, int]] = set()
    for x, y in cells:
        corners.add((x, y))
        corners.add((x + 1, y))
        corners.add((x, y + 1))
        corners.add((x + 1, y + 1))
    return [(float(x), float(y)) for x, y in sorted(corners)]


def enumerate_cover_matrices(det_min: int, det_max: int, entry_bound: int) -> list[MatrixCandidate]:
    out: list[MatrixCandidate] = []
    seen: set[tuple[int, int, int, int]] = set()
    for a11 in range(-entry_bound, entry_bound + 1):
        for a12 in range(-entry_bound, entry_bound + 1):
            for a21 in range(-entry_bound, entry_bound + 1):
                for a22 in range(-entry_bound, entry_bound + 1):
                    det = a11 * a22 - a12 * a21
                    if det < det_min or det > det_max:
                        continue
                    key = (a11, a12, a21, a22)
                    if key in seen:
                        continue
                    seen.add(key)
                    out.append(MatrixCandidate(determinant=det, matrix=[[a11, a12], [a21, a22]]))
    out.sort(key=lambda cand: (cand.determinant, max(abs(x) for row in cand.matrix for x in row), cand.matrix))
    return out


def required_capacity_for_shape(corners: list[Point], scale: float, matrix: MatrixCandidate) -> float:
    a11, a12 = matrix.matrix[0]
    a21, a22 = matrix.matrix[1]
    xs: list[float] = []
    ys: list[float] = []
    sums: list[float] = []
    for u, v in corners:
        x = scale * (a11 * u + a12 * v)
        y = scale * (a21 * u + a22 * v)
        xs.append(x)
        ys.append(y)
        sums.append(x + y)
    tx = max(0.0, -min(xs))
    ty = max(0.0, -min(ys))
    return max(s + tx + ty for s in sums)


def evaluate_shape(cells: tuple[Cell, ...], matrices: list[MatrixCandidate], k: int) -> SnakeCandidate:
    n_cells = len(cells)
    scale = 1.0 / math.sqrt(2.0 * n_cells)
    corners = shape_corners(cells)
    best: SnakeCandidate | None = None
    for matrix in matrices:
        mu = required_capacity_for_shape(corners, scale, matrix)
        ratio = k / (mu * mu)
        candidate = SnakeCandidate(
            n_cells=n_cells,
            determinant=matrix.determinant,
            matrix=matrix.matrix,
            required_capacity_mu=mu,
            achieved_ratio=ratio,
            scale=scale,
            cells=[list(cell) for cell in cells],
            transformed_corner_count=len(corners),
        )
        if best is None or candidate.achieved_ratio > best.achieved_ratio + 1e-12:
            best = candidate
    if best is None:
        raise RuntimeError("no matrix candidates available")
    return best


def search_family(
    k: int,
    det_min: int,
    det_max: int,
    entry_bound: int,
    n_cells_values: list[int],
    shape_trials: int,
    seed: int,
    exact_enumeration: bool,
    path_only: bool,
) -> SearchSummary:
    matrices = enumerate_cover_matrices(det_min=det_min, det_max=det_max, entry_bound=entry_bound)
    rng = random.Random(seed)

    best_overall: SnakeCandidate | None = None
    best_by_n_cells: list[SnakeCandidate] = []
    unique_shape_count_by_n_cells: list[dict[str, int]] = []

    for n_cells in n_cells_values:
        best_for_n: SnakeCandidate | None = None
        seen_shapes: set[tuple[Cell, ...]] = set()

        for shape in canonical_seed_shapes(n_cells):
            if path_only and not is_path_polyomino(shape):
                continue
            seen_shapes.add(shape)
            cand = evaluate_shape(shape, matrices, k=k)
            if best_for_n is None or cand.achieved_ratio > best_for_n.achieved_ratio + 1e-12:
                best_for_n = cand

        if exact_enumeration:
            for shape in enumerate_snake_polyominoes(n_cells):
                if path_only and not is_path_polyomino(shape):
                    continue
                if shape in seen_shapes:
                    continue
                seen_shapes.add(shape)
                cand = evaluate_shape(shape, matrices, k=k)
                if best_for_n is None or cand.achieved_ratio > best_for_n.achieved_ratio + 1e-12:
                    best_for_n = cand
        else:
            for _ in range(shape_trials):
                shape = random_snake_polyomino(n_cells=n_cells, rng=rng)
                if path_only and not is_path_polyomino(shape):
                    continue
                if shape in seen_shapes:
                    continue
                seen_shapes.add(shape)
                cand = evaluate_shape(shape, matrices, k=k)
                if best_for_n is None or cand.achieved_ratio > best_for_n.achieved_ratio + 1e-12:
                    best_for_n = cand

        if best_for_n is None:
            raise RuntimeError(f"failed to find any candidate for n_cells={n_cells}")
        best_by_n_cells.append(best_for_n)
        unique_shape_count_by_n_cells.append({"n_cells": n_cells, "count": len(seen_shapes)})
        if best_overall is None or best_for_n.achieved_ratio > best_overall.achieved_ratio + 1e-12:
            best_overall = best_for_n

    if best_overall is None:
        raise RuntimeError("search produced no candidates")

    return SearchSummary(
        k=k,
        det_min=det_min,
        det_max=det_max,
        entry_bound=entry_bound,
        shape_trials=shape_trials,
        seed=seed,
        exact_enumeration=exact_enumeration,
        best_overall=best_overall,
        best_by_n_cells=best_by_n_cells,
        matrix_count=len(matrices),
        unique_shape_count_by_n_cells=unique_shape_count_by_n_cells,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--det-min", type=int, default=None)
    parser.add_argument("--det-max", type=int, default=None)
    parser.add_argument("--entry-bound", type=int, default=6)
    parser.add_argument("--shape-trials", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--exact-enumeration", action="store_true")
    parser.add_argument("--path-only", action="store_true")
    parser.add_argument("--n-cells", type=int, nargs="+", default=[3, 4, 5, 6, 7, 8, 10, 12, 16])
    args = parser.parse_args()

    det_min = args.k if args.det_min is None else args.det_min
    det_max = args.k + 4 if args.det_max is None else args.det_max

    summary = search_family(
        k=args.k,
        det_min=det_min,
        det_max=det_max,
        entry_bound=args.entry_bound,
        n_cells_values=args.n_cells,
        shape_trials=args.shape_trials,
        seed=args.seed,
        exact_enumeration=args.exact_enumeration,
        path_only=args.path_only,
    )
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
