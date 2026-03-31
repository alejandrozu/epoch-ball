#!/usr/bin/env python3
"""Heuristic search for multi-block toric packings using axis-aligned triangles.

Each block of size n = a*b corresponds to the standard axis-aligned toric block
with moment triangle
    conv((0,0), (a,0), (0,b)),
which packs n unit balls exactly.

This script searches for disjoint translated copies of such blocks inside the
target simplex Delta(mu) = {(x,y) : x >= 0, y >= 0, x + y <= mu}. The search
is constructive but heuristic because translations are restricted to a fixed
grid.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import asdict, dataclass


@dataclass
class BlockPlacement:
    n_balls: int
    side_a: int
    side_b: int
    translation: list[float]
    vertices: list[list[float]]


def translate_triangle(a: int, b: int, tx: float, ty: float) -> list[tuple[float, float]]:
    return [(tx, ty), (tx + a, ty), (tx, ty + b)]


def inside_simplex(triangle: list[tuple[float, float]], mu: float, tol: float = 1e-12) -> bool:
    for x, y in triangle:
        if x < -tol or y < -tol or x + y > mu + tol:
            return False
    return True


def projection_interval(axis: tuple[float, float], triangle: list[tuple[float, float]]) -> tuple[float, float]:
    values = [axis[0] * x + axis[1] * y for x, y in triangle]
    return min(values), max(values)


def triangles_overlap(
    triangle_a: list[tuple[float, float]],
    triangle_b: list[tuple[float, float]],
    tol: float = 1e-12,
) -> bool:
    axes: list[tuple[float, float]] = []
    for triangle in (triangle_a, triangle_b):
        for i in range(3):
            x1, y1 = triangle[i]
            x2, y2 = triangle[(i + 1) % 3]
            dx, dy = x2 - x1, y2 - y1
            axes.append((-dy, dx))
    for axis in axes:
        lo_a, hi_a = projection_interval(axis, triangle_a)
        lo_b, hi_b = projection_interval(axis, triangle_b)
        if hi_a <= lo_b + tol or hi_b <= lo_a + tol:
            return False
    return True


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


def candidate_translations(a: int, b: int, mu: float, grid_step: float) -> list[tuple[float, float, list[tuple[float, float]]]]:
    max_side = max(a, b)
    if mu < max_side:
        return []
    slack = mu - max_side
    max_index = int(math.floor(slack / grid_step + 1e-12))
    out = []
    for i in range(max_index + 1):
        tx = round(i * grid_step, 10)
        for j in range(max_index + 1):
            ty = round(j * grid_step, 10)
            if tx + ty > slack + 1e-12:
                continue
            triangle = translate_triangle(a, b, tx, ty)
            if inside_simplex(triangle, mu):
                out.append((tx, ty, triangle))
    return out


def search_axis_blocks(k: int, mu_values: list[float], grid_step: float, block_limit: int) -> dict:
    partitions = [p for p in integer_partitions(k) if len(p) <= block_limit]
    factor_cache = {n: factor_shapes(n) for n in range(1, k + 1)}
    best_result: dict | None = None

    for mu in mu_values:
        for partition in partitions:
            shape_lists = [factor_cache[n] for n in partition]
            for chosen_shapes in itertools.product(*shape_lists):
                blocks = sorted(
                    zip(partition, chosen_shapes),
                    key=lambda item: (item[1][0] * item[1][1], max(item[1])),
                    reverse=True,
                )
                grids = []
                feasible = True
                for n_balls, (a, b) in blocks:
                    translations = candidate_translations(a, b, mu=mu, grid_step=grid_step)
                    if not translations:
                        feasible = False
                        break
                    grids.append((n_balls, a, b, translations))
                if not feasible:
                    continue

                placed: list[BlockPlacement] = []

                def backtrack(index: int) -> bool:
                    if index == len(grids):
                        return True
                    n_balls, a, b, translations = grids[index]
                    for tx, ty, triangle in translations:
                        if any(triangles_overlap(triangle, existing.vertices) for existing in placed):
                            continue
                        placed.append(
                            BlockPlacement(
                                n_balls=n_balls,
                                side_a=a,
                                side_b=b,
                                translation=[tx, ty],
                                vertices=[[x, y] for x, y in triangle],
                            )
                        )
                        if backtrack(index + 1):
                            return True
                        placed.pop()
                    return False

                if not backtrack(0):
                    continue

                result = {
                    "k": k,
                    "mu": mu,
                    "ratio": k / (mu * mu),
                    "grid_step": grid_step,
                    "block_limit": block_limit,
                    "partition": list(partition),
                    "shape_choice": [[a, b] for a, b in chosen_shapes],
                    "placements": [asdict(block) for block in placed],
                }
                if best_result is None or result["ratio"] > best_result["ratio"]:
                    best_result = result
                break
            if best_result is not None and abs(best_result["mu"] - mu) <= 1e-12:
                break
        if best_result is not None and abs(best_result["mu"] - mu) <= 1e-12:
            break

    if best_result is None:
        return {
            "k": k,
            "mu_values": mu_values,
            "grid_step": grid_step,
            "block_limit": block_limit,
            "found": False,
        }
    best_result["found"] = True
    return best_result


def parse_mu_values(mu_min: float, mu_max: float, mu_step: float) -> list[float]:
    out = []
    count = int(math.floor((mu_max - mu_min) / mu_step + 1e-12))
    for i in range(count + 1):
        out.append(round(mu_min + i * mu_step, 10))
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--mu-min", type=float, required=True)
    parser.add_argument("--mu-max", type=float, required=True)
    parser.add_argument("--mu-step", type=float, default=0.1)
    parser.add_argument("--grid-step", type=float, default=0.1)
    parser.add_argument("--block-limit", type=int, default=4)
    args = parser.parse_args()

    result = search_axis_blocks(
        k=args.k,
        mu_values=parse_mu_values(args.mu_min, args.mu_max, args.mu_step),
        grid_step=args.grid_step,
        block_limit=args.block_limit,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
