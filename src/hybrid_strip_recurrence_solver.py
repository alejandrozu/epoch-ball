#!/usr/bin/env python3
"""Exact recurrence solver and toric gluing audit for the attempt-011 strip.

This exactifies the `m=6` hybrid proof-of-concept from attempt 010.

The script does three things:

1. searches for a rational seed for the boundary-only recurrence;
2. certifies the resulting source strip realization exactly in `Fraction`
   arithmetic, including the face areas; and
3. audits the induced piecewise toric charts across internal edges.

The key distinction is between:

* moment-space continuity, where the piecewise-affine action map glues exactly
  on the shared edges; and
* full 4D toric continuity, which also requires the angle linear maps
  `A^{-T}` to match on any internal edge carrying a full `T^2` fiber.

For the exact `m=6` strip, the first passes and the second fails on every
internal edge, so this exactification still does not produce a valid 4D
symplectic embedding of the tenth ball.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction

try:
    from src.hybrid_strip_realization import extract_strip_complex
except ModuleNotFoundError:
    from hybrid_strip_realization import extract_strip_complex


M = 6
K = 10
TARGET_DOUBLE_AREA = Fraction(1, M * M)
TARGET_AREA = Fraction(1, 2 * M * M)


def frac_to_str(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


def point_to_json(point: tuple[Fraction, Fraction]) -> list[str]:
    return [frac_to_str(point[0]), frac_to_str(point[1])]


def matrix_to_json(matrix: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]) -> list[list[str]]:
    return [[frac_to_str(entry) for entry in row] for row in matrix]


def matrix_key(matrix: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    return matrix


def add_vectors(u: tuple[Fraction, Fraction], v: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return (u[0] + v[0], u[1] + v[1])


def sub_vectors(u: tuple[Fraction, Fraction], v: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return (u[0] - v[0], u[1] - v[1])


def scale_vector(scale: Fraction, u: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return (scale * u[0], scale * u[1])


def midpoint(u: tuple[Fraction, Fraction], v: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return scale_vector(Fraction(1, 2), add_vectors(u, v))


def determinant2(
    u: tuple[Fraction, Fraction],
    v: tuple[Fraction, Fraction],
    w: tuple[Fraction, Fraction],
) -> Fraction:
    return (v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0])


def invert2(
    matrix: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]
) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    (a, b), (c, d) = matrix
    det = a * d - b * c
    if det == 0:
        raise ZeroDivisionError("singular matrix")
    return ((d / det, -b / det), (-c / det, a / det))


def matmul2(
    left: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    right: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def matvec2(
    matrix: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    vector: tuple[Fraction, Fraction],
) -> tuple[Fraction, Fraction]:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )


def transpose2(
    matrix: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]
) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    return ((matrix[0][0], matrix[1][0]), (matrix[0][1], matrix[1][1]))


def recurrence_from_seed(seed: Fraction) -> tuple[dict[int, Fraction], dict[int, Fraction], dict[int, Fraction]] | None:
    """Evaluate the exact boundary recurrence for the `m=6` strip."""
    c = Fraction(1, M * M)
    x = {0: Fraction(0, 1)}
    y = {0: Fraction(0, 1), 9: seed}

    try:
        for i in range(1, 10):
            previous_y = y[10 - i]
            if not (0 < previous_y < 1):
                return None
            x[i] = x[i - 1] + c / previous_y
            if i < 9 and not (0 < x[i] < 1):
                return None
            if i <= 8:
                y[9 - i] = previous_y - c / (1 - x[i])

        y[10] = y[9] + c
        z = {0: Fraction(0, 1)}
        for j in range(1, 10):
            previous_y = y[9 + j]
            if not (0 < previous_y < 1):
                return None
            z[j] = z[j - 1] + c / (1 - previous_y)
            if j < 9 and not (0 < z[j] < 1):
                return None
            y[10 + j] = previous_y + c / (1 - z[j])
    except ZeroDivisionError:
        return None

    z[10] = Fraction(1, 1)
    return x, y, z


def find_exact_roots(max_denominator: int) -> list[Fraction]:
    roots: set[Fraction] = set()
    for denominator in range(1, max_denominator + 1):
        for numerator in range(1, denominator):
            candidate = Fraction(numerator, denominator)
            data = recurrence_from_seed(candidate)
            if data is None:
                continue
            x, y, _ = data
            if x[9] == 1 and y[19] == 1:
                roots.add(candidate)
    return sorted(roots)


def build_source_boundary_vertices(
    x: dict[int, Fraction], y: dict[int, Fraction], z: dict[int, Fraction]
) -> list[tuple[Fraction, Fraction]]:
    vertices: list[tuple[Fraction, Fraction]] = []
    for idx in range(38):
        if idx <= 8:
            t = z[9 - idx]
            vertices.append((Fraction(0, 1), t))
        elif idx <= 18:
            vertices.append((x[idx - 9], Fraction(0, 1)))
        else:
            t = y[idx - 18]
            vertices.append((1 - t, t))
    return vertices


def target_boundary_vertices(boundary_cycle: list[tuple[int, int]]) -> list[tuple[Fraction, Fraction]]:
    return [(Fraction(a, M), Fraction(b, M)) for a, b in boundary_cycle]


def face_order(
    face: tuple[int, int, int],
    vertices: list[tuple[Fraction, Fraction]],
) -> tuple[int, int, int]:
    a, b, c = face
    det = determinant2(vertices[a], vertices[b], vertices[c])
    if det > 0:
        return (a, b, c)
    if det < 0:
        return (a, c, b)
    raise ValueError(f"degenerate face {face}")


def affine_data(
    face: tuple[int, int, int],
    source_vertices: list[tuple[Fraction, Fraction]],
    target_vertices: list[tuple[Fraction, Fraction]],
) -> dict:
    order = face_order(face, source_vertices)
    s0, s1, s2 = [source_vertices[idx] for idx in order]
    t0, t1, t2 = [target_vertices[idx] for idx in order]

    source_matrix = (
        (s1[0] - s0[0], s2[0] - s0[0]),
        (s1[1] - s0[1], s2[1] - s0[1]),
    )
    target_matrix = (
        (t1[0] - t0[0], t2[0] - t0[0]),
        (t1[1] - t0[1], t2[1] - t0[1]),
    )
    action_matrix = matmul2(target_matrix, invert2(source_matrix))
    offset = sub_vectors(t0, matvec2(action_matrix, s0))
    angle_matrix = transpose2(invert2(action_matrix))
    return {
        "order": order,
        "A": action_matrix,
        "b": offset,
        "A_inv_t": angle_matrix,
    }


def internal_edge_adjacencies(faces: list[tuple[int, int, int]]) -> list[tuple[tuple[int, int], int, int]]:
    edge_to_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
    for face_id, face in enumerate(faces):
        a, b, c = face
        for edge in ((a, b), (a, c), (b, c)):
            edge_to_faces[tuple(sorted(edge))].append(face_id)

    out: list[tuple[tuple[int, int], int, int]] = []
    for edge, owners in sorted(edge_to_faces.items()):
        if len(owners) == 2:
            out.append((edge, owners[0], owners[1]))
    return out


def build_report(max_seed_denominator: int) -> dict:
    boundary_cycle, faces = extract_strip_complex(m=M, k=K, reverse=False)
    boundary_cycle = [tuple(vertex) for vertex in boundary_cycle]
    faces = [tuple(face) for face in faces]

    roots = find_exact_roots(max_denominator=max_seed_denominator)
    if roots != [Fraction(17, 36)]:
        raise RuntimeError(f"unexpected exact seed set: {roots}")

    seed = roots[0]
    recurrence_data = recurrence_from_seed(seed)
    if recurrence_data is None:
        raise RuntimeError("exact seed did not evaluate")
    x, y, z = recurrence_data
    if x[9] != 1 or y[19] != 1:
        raise RuntimeError("end constraints failed")

    source_vertices = build_source_boundary_vertices(x, y, z)
    target_vertices = target_boundary_vertices(boundary_cycle)

    source_double_areas = [determinant2(source_vertices[a], source_vertices[b], source_vertices[c]) for a, b, c in faces]
    target_double_areas = [determinant2(target_vertices[a], target_vertices[b], target_vertices[c]) for a, b, c in faces]
    if any(area != TARGET_DOUBLE_AREA for area in source_double_areas):
        raise RuntimeError("source face area check failed")
    if any(area != TARGET_DOUBLE_AREA for area in target_double_areas):
        raise RuntimeError("target face area check failed")

    affine_pieces = [affine_data(face, source_vertices, target_vertices) for face in faces]
    distinct_action_matrices = {
        matrix_key(piece["A"])
        for piece in affine_pieces
    }
    distinct_angle_matrices = {
        matrix_key(piece["A_inv_t"])
        for piece in affine_pieces
    }

    adjacency_rows = []
    first_obstruction = None
    x_midpoint_zero_count = 0
    full_torus_midpoint_count = 0
    same_action_matrix_count = 0
    same_angle_matrix_count = 0

    for shared_edge, left_face, right_face in internal_edge_adjacencies(faces):
        u, v = shared_edge
        edge_midpoint = midpoint(source_vertices[u], source_vertices[v])
        left_piece = affine_pieces[left_face]
        right_piece = affine_pieces[right_face]

        left_action_mid = add_vectors(matvec2(left_piece["A"], edge_midpoint), left_piece["b"])
        right_action_mid = add_vectors(matvec2(right_piece["A"], edge_midpoint), right_piece["b"])
        midpoint_jump = sub_vectors(left_action_mid, right_action_mid)
        midpoint_jump_zero = midpoint_jump == (Fraction(0, 1), Fraction(0, 1))
        full_torus_here = edge_midpoint[0] > 0 and edge_midpoint[1] > 0 and edge_midpoint[0] + edge_midpoint[1] < 1
        same_action_matrix = left_piece["A"] == right_piece["A"]
        same_angle_matrix = left_piece["A_inv_t"] == right_piece["A_inv_t"]

        if midpoint_jump_zero:
            x_midpoint_zero_count += 1
        if full_torus_here:
            full_torus_midpoint_count += 1
        if same_action_matrix:
            same_action_matrix_count += 1
        if same_angle_matrix:
            same_angle_matrix_count += 1

        row = {
            "shared_edge": [u, v],
            "faces": [left_face, right_face],
            "midpoint": point_to_json(edge_midpoint),
            "midpoint_has_full_torus_fiber": full_torus_here,
            "action_midpoint_jump": point_to_json(midpoint_jump),
            "same_action_matrix": same_action_matrix,
            "same_angle_matrix": same_angle_matrix,
        }
        adjacency_rows.append(row)

        if first_obstruction is None and full_torus_here and not same_angle_matrix:
            first_obstruction = {
                "shared_edge": [u, v],
                "faces": [left_face, right_face],
                "midpoint": point_to_json(edge_midpoint),
                "left_action_matrix": matrix_to_json(left_piece["A"]),
                "right_action_matrix": matrix_to_json(right_piece["A"]),
                "left_angle_matrix": matrix_to_json(left_piece["A_inv_t"]),
                "right_angle_matrix": matrix_to_json(right_piece["A_inv_t"]),
            }

    if first_obstruction is None:
        raise RuntimeError("expected a toric gluing obstruction")

    return {
        "m": M,
        "k": K,
        "target_face_area": frac_to_str(TARGET_AREA),
        "target_face_double_area": frac_to_str(TARGET_DOUBLE_AREA),
        "seed_search": {
            "max_denominator": max_seed_denominator,
            "exact_roots": [frac_to_str(root) for root in roots],
            "unique_root_count": len(roots),
        },
        "exact_seed": frac_to_str(seed),
        "end_constraints": {
            "x9": frac_to_str(x[9]),
            "y19": frac_to_str(y[19]),
        },
        "recurrence_coordinates": {
            "x_0_to_9": [frac_to_str(x[idx]) for idx in range(10)],
            "y_0_to_19": [frac_to_str(y[idx]) for idx in range(20)],
            "z_0_to_10": [frac_to_str(z[idx]) for idx in range(11)],
        },
        "source_boundary_vertices": [
            {
                "index": idx,
                "point": point_to_json(point),
            }
            for idx, point in enumerate(source_vertices)
        ],
        "target_boundary_vertices": [
            {
                "index": idx,
                "point": point_to_json(point),
            }
            for idx, point in enumerate(target_vertices)
        ],
        "face_area_check": {
            "face_count": len(faces),
            "all_source_double_areas_equal_target": True,
            "source_double_area_values": sorted({frac_to_str(area) for area in source_double_areas}),
            "target_double_area_values": sorted({frac_to_str(area) for area in target_double_areas}),
        },
        "toric_gluing_audit": {
            "internal_edge_count": len(adjacency_rows),
            "action_midpoint_zero_jump_count": x_midpoint_zero_count,
            "full_torus_internal_edge_count": full_torus_midpoint_count,
            "same_action_matrix_count": same_action_matrix_count,
            "same_angle_matrix_count": same_angle_matrix_count,
            "distinct_action_matrix_count": len(distinct_action_matrices),
            "distinct_angle_matrix_count": len(distinct_angle_matrices),
            "all_action_midpoint_jumps_zero": x_midpoint_zero_count == len(adjacency_rows),
            "all_full_torus_edges_fail_angle_match": (
                full_torus_midpoint_count > 0 and same_angle_matrix_count == 0
            ),
            "first_obstruction": first_obstruction,
            "adjacencies": adjacency_rows,
        },
        "conclusion": {
            "exact_moment_realization_exists": True,
            "piecewise_action_map_is_continuous": True,
            "piecewise_toric_chart_is_continuous": False,
            "reason": (
                "The action maps glue on all 35 internal edges, but every internal edge has a full "
                "T^2 fiber at its midpoint and no adjacent pair has matching A^{-T}. Constant angle "
                "translations cannot repair a nonzero linear mismatch, so the boundary-only strip is "
                "still not a valid 4D toric piecewise-affine symplectic embedding."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--max-seed-denominator",
        type=int,
        default=72,
        help="Search rational recurrence seeds up to this denominator.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="results/research_loop/attempt_011_exact_strip_gluing.json",
        help="Optional JSON output path.",
    )
    args = parser.parse_args()

    report = build_report(max_seed_denominator=args.max_seed_denominator)
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.write("\n")
    print(text)


if __name__ == "__main__":
    main()
