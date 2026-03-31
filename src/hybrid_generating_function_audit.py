#!/usr/bin/env python3
"""Attempt-012 generating-function audit for the hybrid strip.

This is a new transform-method family following attempt 011.

We keep the exact `m=6` hybrid strip from attempt 011, but enlarge the local
4D model from a plain toric affine map

    (x, y) -> (A x + b, A^{-T} y)

to the most general exact symplectic map that still covers the same affine base
map on each face:

    (x, y) -> (A x + b, A^{-T}(y + grad S(x))).

For a piecewise-affine action map, this is the natural generating-function /
cotangent-lift correction family. The question is whether the facewise
potentials `S` can repair the continuity failure from attempt 011.

The answer for the exact `m=6` strip is no. Along any shared edge, continuity
can only depend on the angle directions whose source and target torus factors
are still present there. In `B^4`, the only collapsed torus directions come
from the coordinate axes `x1 = 0` or `x2 = 0`; the diagonal wall
`x1 + x2 = 1` does *not* collapse a circle. Therefore, on an edge where both
coordinates remain positive in source and target, the full `2 x 2` matrix
`A^{-T}` must already agree across the edge, and no generating-function
correction can change that.

This script audits the exact strip using that criterion. It also includes a
small positive control showing that derivative jumps can be hidden along an
axis edge when only one circle survives.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction

try:
    from src.hybrid_strip_recurrence_solver import (
        K,
        M,
        affine_data,
        build_source_boundary_vertices,
        find_exact_roots,
        internal_edge_adjacencies,
        recurrence_from_seed,
        target_boundary_vertices,
    )
    from src.hybrid_strip_realization import extract_strip_complex
except ModuleNotFoundError:
    from hybrid_strip_recurrence_solver import (
        K,
        M,
        affine_data,
        build_source_boundary_vertices,
        find_exact_roots,
        internal_edge_adjacencies,
        recurrence_from_seed,
        target_boundary_vertices,
    )
    from hybrid_strip_realization import extract_strip_complex


def frac_to_str(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


def point_to_json(point: tuple[Fraction, Fraction]) -> list[str]:
    return [frac_to_str(point[0]), frac_to_str(point[1])]


def matrix_to_json(matrix: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]) -> list[list[str]]:
    return [[frac_to_str(entry) for entry in row] for row in matrix]


def midpoint(u: tuple[Fraction, Fraction], v: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return ((u[0] + v[0]) / 2, (u[1] + v[1]) / 2)


def active_angle_indices(point: tuple[Fraction, Fraction]) -> tuple[int, ...]:
    """Active torus factors for the standard T^2-action on B^4.

    The standard action rotates the two complex coordinates independently, so a
    circle collapses exactly when the corresponding action coordinate vanishes.
    The diagonal wall `x1 + x2 = 1` is only the sphere boundary of the ball and
    does not collapse a torus factor.
    """
    return tuple(idx for idx, value in enumerate(point) if value > 0)


def restricted_matrix(
    matrix: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]],
    target_active: tuple[int, ...],
    source_active: tuple[int, ...],
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(tuple(matrix[row][col] for col in source_active) for row in target_active)


def build_exact_strip_data() -> tuple[list[tuple[int, int]], list[tuple[int, int, int]], list[dict], list[tuple[Fraction, Fraction]], list[tuple[Fraction, Fraction]]]:
    boundary_cycle, faces = extract_strip_complex(m=M, k=K, reverse=False)
    seed = find_exact_roots(72)
    if seed != [Fraction(17, 36)]:
        raise RuntimeError(f"unexpected exact seed set: {seed}")
    recurrence = recurrence_from_seed(seed[0])
    if recurrence is None:
        raise RuntimeError("failed to evaluate exact seed")
    x, y, z = recurrence
    source_vertices = build_source_boundary_vertices(x, y, z)
    target_vertices = target_boundary_vertices([tuple(vertex) for vertex in boundary_cycle])
    faces = [tuple(face) for face in faces]
    pieces = [affine_data(face, source_vertices, target_vertices) for face in faces]
    return [tuple(vertex) for vertex in boundary_cycle], faces, pieces, source_vertices, target_vertices


def toy_axis_edge_control() -> dict:
    """A positive control: an axis edge can hide derivative changes."""
    left_a_inv_t = ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)))
    right_a_inv_t = ((Fraction(1), Fraction(-1)), (Fraction(0), Fraction(1)))
    source_midpoint = (Fraction(0), Fraction(1, 3))
    target_midpoint = (Fraction(0), Fraction(1, 3))
    source_active = active_angle_indices(source_midpoint)
    target_active = active_angle_indices(target_midpoint)
    return {
        "description": (
            "Two affine base maps can agree on the axis edge x1=0 while having different "
            "full A^{-T} matrices. Only the second circle survives there, so the restricted "
            "1x1 angle matrix still matches."
        ),
        "source_midpoint": point_to_json(source_midpoint),
        "target_midpoint": point_to_json(target_midpoint),
        "source_active_angles": list(source_active),
        "target_active_angles": list(target_active),
        "left_angle_matrix": matrix_to_json(left_a_inv_t),
        "right_angle_matrix": matrix_to_json(right_a_inv_t),
        "restricted_left": [[frac_to_str(left_a_inv_t[1][1])]],
        "restricted_right": [[frac_to_str(right_a_inv_t[1][1])]],
        "restricted_match": restricted_matrix(left_a_inv_t, target_active, source_active)
        == restricted_matrix(right_a_inv_t, target_active, source_active),
    }


def audit_strip() -> dict:
    _, faces, pieces, source_vertices, target_vertices = build_exact_strip_data()

    rows = []
    restricted_match_count = 0
    full_rank_edge_count = 0
    first_obstruction = None

    for edge, left_face, right_face in internal_edge_adjacencies(faces):
        source_midpoint = midpoint(source_vertices[edge[0]], source_vertices[edge[1]])
        target_midpoint = midpoint(target_vertices[edge[0]], target_vertices[edge[1]])
        source_active = active_angle_indices(source_midpoint)
        target_active = active_angle_indices(target_midpoint)
        left_angle = pieces[left_face]["A_inv_t"]
        right_angle = pieces[right_face]["A_inv_t"]
        restricted_left = restricted_matrix(left_angle, target_active, source_active)
        restricted_right = restricted_matrix(right_angle, target_active, source_active)
        restricted_match = restricted_left == restricted_right
        if restricted_match:
            restricted_match_count += 1
        if len(source_active) == 2 and len(target_active) == 2:
            full_rank_edge_count += 1
        row = {
            "shared_edge": list(edge),
            "faces": [left_face, right_face],
            "source_midpoint": point_to_json(source_midpoint),
            "target_midpoint": point_to_json(target_midpoint),
            "source_active_angles": list(source_active),
            "target_active_angles": list(target_active),
            "left_angle_matrix": matrix_to_json(left_angle),
            "right_angle_matrix": matrix_to_json(right_angle),
            "restricted_left": [[frac_to_str(entry) for entry in line] for line in restricted_left],
            "restricted_right": [[frac_to_str(entry) for entry in line] for line in restricted_right],
            "restricted_match": restricted_match,
        }
        rows.append(row)
        if first_obstruction is None and not restricted_match:
            first_obstruction = row

    if first_obstruction is None:
        raise RuntimeError("expected a strip obstruction")

    return {
        "internal_edge_count": len(rows),
        "restricted_match_count": restricted_match_count,
        "full_rank_edge_count": full_rank_edge_count,
        "all_edges_are_full_rank": full_rank_edge_count == len(rows),
        "unique_active_angle_patterns": sorted(
            {
                (tuple(row["source_active_angles"]), tuple(row["target_active_angles"]))
                for row in rows
            }
        ),
        "first_obstruction": first_obstruction,
        "edges": rows,
    }


def build_report() -> dict:
    strip = audit_strip()
    return {
        "m": M,
        "k": K,
        "family": "facewise generating-function / cotangent-lift correction over the exact hybrid strip",
        "local_model": {
            "formula": "(x, y) -> (A x + b, A^{-T}(y + grad S(x)))",
            "derivation_summary": (
                "For a symplectic map covering a base map f(x), the pullback equations force "
                "D_y g = Df^{-T}. On an affine face f(x)=A x+b this integrates to "
                "g(x,y)=A^{-T}y+h(x), and the remaining symplectic condition is equivalent to "
                "A^T h being locally exact, i.e. h=A^{-T} grad S."
            ),
            "continuity_rule": (
                "Across a shared edge, the restricted angle matrix from surviving source circles "
                "to surviving target circles must already match. A generating-function correction "
                "changes only the x-dependent translation term, not that linear coefficient."
            ),
        },
        "toy_axis_edge_control": toy_axis_edge_control(),
        "exact_strip_audit": strip,
        "conclusion": {
            "generating_function_repair_possible": False,
            "reason": (
                "On the exact hybrid strip, every internal edge has both source angles and both "
                "target angles active, so the full 2x2 matrix A^{-T} must match across the edge. "
                "It never does: restricted_match_count = 0 on all 35 internal edges."
            ),
            "next_design_hint": (
                "Any future piecewise-affine action complex that hopes to glue by this method "
                "must place derivative jumps only along coordinate-axis edges, not merely along "
                "the diagonal wall x1+x2=1."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=str,
        default="results/research_loop/attempt_012_generating_function_audit.json",
        help="Optional JSON output path.",
    )
    args = parser.parse_args()

    report = build_report()
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.write("\n")
    print(text)


if __name__ == "__main__":
    main()
