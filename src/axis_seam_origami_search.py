#!/usr/bin/env python3
"""Attempt-013 axis-seam origami search.

This attempt follows the obstruction from attempt 012:

* in `B^4`, derivative jumps can only be hidden on seams whose source or target
  image lies on a true coordinate axis; the diagonal wall does not collapse a
  torus factor;
* therefore the next plausible piecewise-toric family is an "axis-seam
  origami", where different rigid affine pieces glue only along seams mapped to
  `x1 = 0` or `x2 = 0`.

This script tests the simplest genuinely new member of that family:

* split one source side-`m` triangulated moment triangle into two connected
  rigid pieces by an internal straight lattice cut;
* map each piece by a rigid triangular-lattice automorphism plus translation
  into the target side-`n` simplex, where `n = ceil(m * sqrt(k))`;
* force the common cut to map to a coordinate axis; and
* require the two target images to be interior-disjoint.

Why only straight cuts?
If two different affine maps agree on a bent seam containing two nonparallel
segments, their difference vanishes on a 2D affine span, so the maps are equal.
Hence a nontrivial two-piece rigid origami can only meet along a straight
lattice segment. Those are exactly the internal lines `i = c`, `j = c`, and
`i + j = c`.

The search below is exact and finite. It enumerates all rigid placements of the
two source pieces, groups them by the induced axis image of the common seam, and
checks whether any pair of different rigid maps gives disjoint target images.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from functools import lru_cache


Point = tuple[int, int]
Triangle = tuple[Point, Point, Point]
Matrix = tuple[tuple[int, int], tuple[int, int]]


def standard_unit_triangle(triangle: Triangle) -> bool:
    a, b, c = triangle
    twice_area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    if abs(twice_area) != 1:
        return False
    directions = {
        (1, 0),
        (0, 1),
        (-1, 1),
        (-1, 0),
        (0, -1),
        (1, -1),
    }
    edges = (
        (b[0] - a[0], b[1] - a[1]),
        (c[0] - b[0], c[1] - b[1]),
        (a[0] - c[0], a[1] - c[1]),
    )
    return all(edge in directions for edge in edges)


def lattice_automorphisms() -> list[Matrix]:
    out: list[Matrix] = []
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            for c in (-1, 0, 1):
                for d in (-1, 0, 1):
                    det = a * d - b * c
                    if abs(det) != 1:
                        continue
                    triangle = tuple(sorted(((0, 0), (a, c), (b, d))))
                    if standard_unit_triangle(triangle):
                        out.append(((a, b), (c, d)))
    out = sorted(set(out))
    if len(out) != 12:
        raise RuntimeError(f"expected 12 triangular-lattice automorphisms, got {len(out)}")
    return out


AUTOMORPHISMS = lattice_automorphisms()


def matvec(matrix: Matrix, point: Point) -> Point:
    return (
        matrix[0][0] * point[0] + matrix[0][1] * point[1],
        matrix[1][0] * point[0] + matrix[1][1] * point[1],
    )


def add_points(left: Point, right: Point) -> Point:
    return (left[0] + right[0], left[1] + right[1])


def triangle_barycenter_numerators(triangle: Triangle) -> tuple[int, int]:
    return (
        triangle[0][0] + triangle[1][0] + triangle[2][0],
        triangle[0][1] + triangle[1][1] + triangle[2][1],
    )


def subdivide_triangle(side: int) -> tuple[Triangle, ...]:
    out: list[Triangle] = []
    for i in range(side):
        for j in range(side - i):
            out.append(((i, j), (i + 1, j), (i, j + 1)))
    for i in range(side - 1):
        for j in range(side - 1 - i):
            out.append(((i + 1, j), (i + 1, j + 1), (i, j + 1)))
    return tuple(out)


def cut_partition(side: int, family: str, c: int) -> tuple[tuple[Triangle, ...], tuple[Triangle, ...]]:
    left: list[Triangle] = []
    right: list[Triangle] = []
    for triangle in subdivide_triangle(side):
        sx, sy = triangle_barycenter_numerators(triangle)
        if family == "i":
            value = sx - 3 * c
        elif family == "j":
            value = sy - 3 * c
        elif family == "sum":
            value = sx + sy - 3 * c
        else:
            raise ValueError(f"unknown cut family {family}")
        if value < 0:
            left.append(triangle)
        elif value > 0:
            right.append(triangle)
        else:
            raise RuntimeError(f"unexpected triangle barycenter on cut family={family} c={c}")
    return tuple(left), tuple(right)


def cut_vertices(side: int, family: str, c: int) -> tuple[Point, ...]:
    if family == "i":
        return tuple((c, j) for j in range(side - c + 1))
    if family == "j":
        return tuple((i, c) for i in range(side - c + 1))
    if family == "sum":
        return tuple((i, c - i) for i in range(c + 1))
    raise ValueError(f"unknown cut family {family}")


def map_triangle(matrix: Matrix, shift: Point, triangle: Triangle) -> Triangle:
    return tuple(sorted(add_points(matvec(matrix, vertex), shift) for vertex in triangle))


def inside_target(triangle: Triangle, side: int) -> bool:
    return all(vertex[0] >= 0 and vertex[1] >= 0 and vertex[0] + vertex[1] <= side for vertex in triangle)


def axis_line_key(image_vertices: tuple[Point, ...]) -> tuple[str, tuple[Point, ...]] | None:
    if all(vertex[0] == 0 for vertex in image_vertices):
        return ("x1=0", image_vertices)
    if all(vertex[1] == 0 for vertex in image_vertices):
        return ("x2=0", image_vertices)
    return None


def line_image(matrix: Matrix, shift: Point, vertices: tuple[Point, ...]) -> tuple[Point, ...]:
    return tuple(add_points(matvec(matrix, vertex), shift) for vertex in vertices)


def triangle_to_json(triangle: Triangle) -> list[list[int]]:
    return [[point[0], point[1]] for point in triangle]


def matrix_to_json(matrix: Matrix) -> list[list[int]]:
    return [[matrix[0][0], matrix[0][1]], [matrix[1][0], matrix[1][1]]]


def placement_to_json(matrix: Matrix, shift: Point) -> dict:
    return {
        "A": matrix_to_json(matrix),
        "b": [shift[0], shift[1]],
    }


@lru_cache(maxsize=None)
def rigid_region_placements(region: tuple[Triangle, ...], target_side: int) -> tuple[tuple[Matrix, Point, frozenset[Triangle]], ...]:
    placements: list[tuple[Matrix, Point, frozenset[Triangle]]] = []
    for matrix in AUTOMORPHISMS:
        mapped_without_shift = [map_triangle(matrix, (0, 0), triangle) for triangle in region]
        vertices = {vertex for triangle in mapped_without_shift for vertex in triangle}
        min_i = min(vertex[0] for vertex in vertices)
        max_i = max(vertex[0] for vertex in vertices)
        min_j = min(vertex[1] for vertex in vertices)
        max_j = max(vertex[1] for vertex in vertices)
        for bi in range(-min_i, target_side - max_i + 1):
            for bj in range(-min_j, target_side - max_j + 1):
                shift = (bi, bj)
                mapped = [map_triangle(matrix, shift, triangle) for triangle in region]
                if not all(inside_target(triangle, target_side) for triangle in mapped):
                    continue
                mapped_set = frozenset(mapped)
                if len(mapped_set) != len(region):
                    continue
                placements.append((matrix, shift, mapped_set))
    return tuple(placements)


def cut_summary(k: int, m: int, family: str, c: int) -> dict:
    n = math.ceil(m * math.sqrt(k))
    left_region, right_region = cut_partition(m, family, c)
    seam_vertices = cut_vertices(m, family, c)

    grouped: list[dict[tuple[str, tuple[Point, ...]], list[tuple[Matrix, Point, frozenset[Triangle]]]]] = [
        defaultdict(list),
        defaultdict(list),
    ]

    side_axis_counts = [0, 0]
    for side_idx, region in enumerate((left_region, right_region)):
        for matrix, shift, mapped_set in rigid_region_placements(region, n):
            key = axis_line_key(line_image(matrix, shift, seam_vertices))
            if key is None:
                continue
            grouped[side_idx][key].append((matrix, shift, mapped_set))
            side_axis_counts[side_idx] += 1

    matching_pair_count = 0
    valid_pair_count = 0
    min_overlap_count: int | None = None
    best_overlap_witness: dict | None = None
    first_valid_pair: dict | None = None

    for key, left_items in grouped[0].items():
        right_items = grouped[1].get(key, [])
        if not right_items:
            continue
        for matrix_left, shift_left, tris_left in left_items:
            left_set = set(tris_left)
            for matrix_right, shift_right, tris_right in right_items:
                if matrix_left == matrix_right and shift_left == shift_right:
                    continue
                matching_pair_count += 1
                overlap = left_set & set(tris_right)
                overlap_count = len(overlap)
                if min_overlap_count is None or overlap_count < min_overlap_count:
                    min_overlap_count = overlap_count
                    best_overlap_witness = {
                        "axis_family": key[0],
                        "axis_line_vertices": [[point[0], point[1]] for point in key[1]],
                        "left": placement_to_json(matrix_left, shift_left),
                        "right": placement_to_json(matrix_right, shift_right),
                        "overlap_count": overlap_count,
                        "overlap_triangles_sample": [triangle_to_json(triangle) for triangle in sorted(overlap)[:5]],
                    }
                if overlap_count == 0:
                    valid_pair_count += 1
                    if first_valid_pair is None:
                        first_valid_pair = {
                            "axis_family": key[0],
                            "axis_line_vertices": [[point[0], point[1]] for point in key[1]],
                            "left": placement_to_json(matrix_left, shift_left),
                            "right": placement_to_json(matrix_right, shift_right),
                        }

    return {
        "cut_family": family,
        "cut_parameter": c,
        "cut_length_edges": len(seam_vertices) - 1,
        "left_piece_triangle_count": len(left_region),
        "right_piece_triangle_count": len(right_region),
        "left_axis_aligned_placement_count": side_axis_counts[0],
        "right_axis_aligned_placement_count": side_axis_counts[1],
        "matching_axis_line_pair_count": matching_pair_count,
        "valid_nonoverlapping_pair_count": valid_pair_count,
        "min_overlap_count": min_overlap_count,
        "best_overlap_witness": best_overlap_witness,
        "first_valid_pair": first_valid_pair,
    }


def run_search(k: int, m_values: list[int]) -> dict:
    rows = []
    total_valid_pairs = 0
    for m in m_values:
        n = math.ceil(m * math.sqrt(k))
        cut_rows = []
        global_min_overlap: int | None = None
        best_cut_witness: dict | None = None
        for family in ("i", "j", "sum"):
            for c in range(1, m):
                row = cut_summary(k=k, m=m, family=family, c=c)
                cut_rows.append(row)
                total_valid_pairs += row["valid_nonoverlapping_pair_count"]
                overlap = row["min_overlap_count"]
                if overlap is not None and (global_min_overlap is None or overlap < global_min_overlap):
                    global_min_overlap = overlap
                    best_cut_witness = {
                        "cut_family": row["cut_family"],
                        "cut_parameter": row["cut_parameter"],
                        "cut_length_edges": row["cut_length_edges"],
                        "left_piece_triangle_count": row["left_piece_triangle_count"],
                        "right_piece_triangle_count": row["right_piece_triangle_count"],
                        "best_overlap_witness": row["best_overlap_witness"],
                    }
        rows.append(
            {
                "m": m,
                "n": n,
                "straight_cut_count": len(cut_rows),
                "axis_aligned_left_placements_total": sum(row["left_axis_aligned_placement_count"] for row in cut_rows),
                "axis_aligned_right_placements_total": sum(row["right_axis_aligned_placement_count"] for row in cut_rows),
                "matching_axis_line_pair_count_total": sum(row["matching_axis_line_pair_count"] for row in cut_rows),
                "valid_nonoverlapping_pair_count_total": sum(row["valid_nonoverlapping_pair_count"] for row in cut_rows),
                "global_min_overlap_count": global_min_overlap,
                "best_cut_witness": best_cut_witness,
                "cuts": cut_rows,
            }
        )

    return {
        "family": "axis_seam_two_piece_rigid_origami",
        "description": (
            "Two rigid triangular-lattice pieces glued along an internal straight cut whose "
            "common image is forced onto x1=0 or x2=0."
        ),
        "k": k,
        "tested_m_values": m_values,
        "triangular_lattice_automorphism_count": len(AUTOMORPHISMS),
        "triangular_lattice_automorphisms": [matrix_to_json(matrix) for matrix in AUTOMORPHISMS],
        "m_rows": rows,
        "total_valid_nonoverlapping_pairs": total_valid_pairs,
        "exists_axis_seam_two_piece_example": total_valid_pairs > 0,
        "conclusion": (
            "No two-piece rigid axis-seam origami was found in the tested near-full k=10 regime."
            if total_valid_pairs == 0
            else "At least one two-piece rigid axis-seam origami exists in the tested regime."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--m-values", type=int, nargs="+", default=[3, 4, 5, 6])
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()

    report = run_search(k=args.k, m_values=args.m_values)
    if args.output is None:
        print(json.dumps(report, indent=2, sort_keys=True))
        return

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
