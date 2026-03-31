#!/usr/bin/env python3
"""Attempt-014 local obstruction for rigid axis-seam origami.

This attempt is a proof-style follow-up to attempt 013 rather than another
placement search.

Attempt 013 exact-searched the simplest rigid axis-seam family and found:

* many candidate placements whose common seam lands on `x1 = 0` or `x2 = 0`;
* zero interior-disjoint realizations.

The local reason is simple and stronger than the search:

1. In the target simplex, every unit boundary edge on `x1 = 0` or `x2 = 0`
   has exactly one incident unit triangle.
2. If two different rigid pieces share a seam edge and both maps agree on that
   seam edge, then the two source unit triangles adjacent to the seam edge must
   map to unit triangles incident to the same target edge.
3. On an axis boundary edge there is only one such target triangle, so the two
   images coincide. This forces overlap immediately.

Therefore any rigid internal seam whose image lies on a true coordinate axis
forces overlap edge-by-edge. In particular, the whole rigid axis-seam strategy
is impossible, not just the two-piece straight-cut cases from attempt 013.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

try:
    from src.axis_seam_origami_search import cut_vertices, subdivide_triangle
except ModuleNotFoundError:
    from axis_seam_origami_search import cut_vertices, subdivide_triangle


Point = tuple[int, int]
Edge = tuple[Point, Point]
Triangle = tuple[Point, Point, Point]


def canonical_edge(left: Point, right: Point) -> Edge:
    return tuple(sorted((left, right)))


def triangle_edges(triangle: Triangle) -> tuple[Edge, Edge, Edge]:
    a, b, c = triangle
    return (
        canonical_edge(a, b),
        canonical_edge(b, c),
        canonical_edge(a, c),
    )


def edge_incidence(side: int) -> dict[Edge, list[Triangle]]:
    incidence: dict[Edge, list[Triangle]] = defaultdict(list)
    for triangle in subdivide_triangle(side):
        for edge in triangle_edges(triangle):
            incidence[edge].append(triangle)
    return dict(incidence)


def edge_class(edge: Edge, side: int) -> str:
    (a1, a2), (b1, b2) = edge
    if a1 == 0 and b1 == 0:
        return "x1=0"
    if a2 == 0 and b2 == 0:
        return "x2=0"
    if a1 + a2 == side and b1 + b2 == side:
        return "diag"
    return "interior"


def point_to_json(point: Point) -> list[int]:
    return [point[0], point[1]]


def edge_to_json(edge: Edge) -> list[list[int]]:
    return [point_to_json(edge[0]), point_to_json(edge[1])]


def triangle_to_json(triangle: Triangle) -> list[list[int]]:
    return [point_to_json(vertex) for vertex in triangle]


def local_boundary_controls(sample_side: int) -> dict:
    incidence = edge_incidence(sample_side)
    axis_x1_edge = canonical_edge((0, sample_side - 1), (0, sample_side))
    axis_x2_edge = canonical_edge((sample_side - 1, 0), (sample_side, 0))
    diagonal_edge = canonical_edge((0, sample_side), (1, sample_side - 1))
    interior_edge = canonical_edge((1, 1), (1, 2))
    controls = {}
    for name, edge in (
        ("x1_axis_edge", axis_x1_edge),
        ("x2_axis_edge", axis_x2_edge),
        ("diagonal_boundary_edge", diagonal_edge),
        ("interior_edge", interior_edge),
    ):
        triangles = incidence[edge]
        controls[name] = {
            "edge": edge_to_json(edge),
            "edge_class": edge_class(edge, sample_side),
            "incident_triangle_count": len(triangles),
            "incident_triangles": [triangle_to_json(triangle) for triangle in triangles],
        }
    return controls


def verify_boundary_uniqueness(max_side: int) -> dict:
    violations = []
    counts = {
        "x1=0": 0,
        "x2=0": 0,
        "diag": 0,
        "interior": 0,
    }
    for side in range(2, max_side + 1):
        incidence = edge_incidence(side)
        for edge, triangles in incidence.items():
            cls = edge_class(edge, side)
            counts[cls] += 1
            expected = 2 if cls == "interior" else 1
            if len(triangles) != expected:
                violations.append(
                    {
                        "side": side,
                        "edge": edge_to_json(edge),
                        "edge_class": cls,
                        "incident_triangle_count": len(triangles),
                    }
                )
    return {
        "max_side_checked": max_side,
        "edge_class_counts": counts,
        "all_counts_match_expected": not violations,
        "violations": violations,
    }


def seam_edges_from_vertices(vertices: tuple[Point, ...]) -> list[Edge]:
    return [canonical_edge(vertices[idx], vertices[idx + 1]) for idx in range(len(vertices) - 1)]


def canonical_axis_edges(length: int, axis: str) -> list[Edge]:
    if axis == "x1=0":
        return [canonical_edge((0, j), (0, j + 1)) for j in range(length)]
    if axis == "x2=0":
        return [canonical_edge((j, 0), (j + 1, 0)) for j in range(length)]
    raise ValueError(f"unknown axis {axis}")


def distinct_boundary_triangles_for_axis_path(length: int, axis: str) -> dict:
    side = length + 2
    incidence = edge_incidence(side)
    edges = canonical_axis_edges(length, axis)
    triangles = [incidence[edge][0] for edge in edges]
    return {
        "axis": axis,
        "path_length_edges": length,
        "edge_count": len(edges),
        "distinct_incident_triangle_count": len(set(triangles)),
        "edges": [edge_to_json(edge) for edge in edges],
        "incident_triangles": [triangle_to_json(triangle) for triangle in triangles],
    }


def verify_axis_path_overlap_bound(max_length: int) -> dict:
    rows = []
    violations = []
    for axis in ("x1=0", "x2=0"):
        for length in range(1, max_length + 1):
            row = distinct_boundary_triangles_for_axis_path(length, axis)
            row["overlap_lower_bound_equals_length"] = row["distinct_incident_triangle_count"] == length
            rows.append(row)
            if not row["overlap_lower_bound_equals_length"]:
                violations.append(row)
    return {
        "max_length_checked": max_length,
        "all_axis_paths_force_at_least_length_many_distinct_overlaps": not violations,
        "violations": violations,
        "rows": rows,
    }


def audit_attempt_013_artifact(path: Path) -> dict:
    data = json.loads(path.read_text())
    rows = []
    failures = []
    total_valid_pairs = 0
    for m_row in data["m_rows"]:
        m = m_row["m"]
        for cut in m_row["cuts"]:
            seam_vertices = cut_vertices(m, cut["cut_family"], cut["cut_parameter"])
            seam_edges = seam_edges_from_vertices(seam_vertices)
            lower_bound = len(seam_edges)
            result_row = {
                "m": m,
                "cut_family": cut["cut_family"],
                "cut_parameter": cut["cut_parameter"],
                "seam_length_edges": lower_bound,
                "artifact_min_overlap_count": cut["min_overlap_count"],
                "artifact_valid_pair_count": cut["valid_nonoverlapping_pair_count"],
                "passes_overlap_bound": cut["min_overlap_count"] >= lower_bound,
                "passes_zero_valid_pair_check": cut["valid_nonoverlapping_pair_count"] == 0,
            }
            rows.append(result_row)
            total_valid_pairs += cut["valid_nonoverlapping_pair_count"]
            if not result_row["passes_overlap_bound"] or not result_row["passes_zero_valid_pair_check"]:
                failures.append(result_row)
    return {
        "artifact_path": str(path),
        "tested_m_values": data["tested_m_values"],
        "audited_cut_count": len(rows),
        "total_valid_pair_count_from_artifact": total_valid_pairs,
        "all_cuts_satisfy_local_overlap_bound": not failures,
        "failures": failures,
        "rows": rows,
    }


def theorem_text() -> str:
    return (
        "If two different rigid unimodular-affine pieces share an internal seam edge "
        "whose common image is a unit edge on x1=0 or x2=0 in the target simplex, then "
        "the two source unit triangles adjacent to that seam edge must map to the same "
        "unique target unit triangle incident to that boundary edge. Hence each seam "
        "edge forces one overlapping target triangle; a seam of length L forces at least "
        "L overlaps."
    )


def build_report(max_side: int, max_length: int, artifact_path: Path) -> dict:
    boundary_check = verify_boundary_uniqueness(max_side=max_side)
    path_check = verify_axis_path_overlap_bound(max_length=max_length)
    artifact_audit = audit_attempt_013_artifact(artifact_path)
    return {
        "family": "local_boundary_uniqueness_obstruction_for_rigid_axis_seams",
        "theorem": theorem_text(),
        "sample_local_controls": local_boundary_controls(sample_side=max(4, min(8, max_side))),
        "boundary_uniqueness_check": boundary_check,
        "axis_path_overlap_check": path_check,
        "attempt_013_artifact_audit": artifact_audit,
        "conclusion": (
            "Every rigid internal seam mapped to x1=0 or x2=0 forces overlap, so the "
            "rigid axis-seam origami strategy is impossible."
            if boundary_check["all_counts_match_expected"]
            and path_check["all_axis_paths_force_at_least_length_many_distinct_overlaps"]
            and artifact_audit["all_cuts_satisfy_local_overlap_bound"]
            else "The local obstruction checks did not all pass."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-side", type=int, default=32)
    parser.add_argument("--max-length", type=int, default=16)
    parser.add_argument(
        "--artifact",
        type=Path,
        default=Path("results/research_loop/attempt_013_axis_seam_two_piece_search.json"),
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    report = build_report(max_side=args.max_side, max_length=args.max_length, artifact_path=args.artifact)
    if args.output is None:
        print(json.dumps(report, indent=2, sort_keys=True))
        return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
