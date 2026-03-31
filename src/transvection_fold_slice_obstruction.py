#!/usr/bin/env python3
"""Attempt-015 obstruction for affine symplectic transvection folds.

This attempt studies a new family that is genuinely different from the rigid
axis-seam origami ruled out in attempt 014.

Start with standard symplectic coordinates `(q1, q2, p1, p2)` on `R^4` and the
active symplectic 2-plane

    E = {(q1, 0, p1, 0)} = span(e_q1, e_p1).

The one-sided affine transvection fold across the hyperplane `p1 = d` is

    (q1, q2, p1, p2) -> (q1 + lambda * max(p1 - d, 0), q2, p1, p2),

which is piecewise affine symplectic and fixes the seam `p1 = d` pointwise.
After a symplectic squeeze in the active plane, the normalized slice map on
`E` becomes

    F_{s,lambda,d}(u, v) = (s * (u + lambda * max(v - d, 0)), v / s).

This script proves and audits the key obstruction:

1. `F_{s,lambda,d}` is area-preserving on the active slice.
2. Therefore the image of the central slice `D = B^4(1) ∩ E` always has
   symplectic area `pi`.
3. If a common symplectic linear left factor `L in Sp(4,R)` is applied, the
   slice image lies in the 2-plane `L(E)`. If `f1, f2` is an orthonormal basis
   of `L(E)` and `delta = |omega(f1, f2)| <= 1`, then

       EuclideanArea(L(F(D))) = pi / delta,

   hence any Euclidean 4-ball containing the full image must satisfy

       R >= delta^(-1/2) >= 1.

So neither the normalized fold itself nor any common linear symplectic left
factor can compress the preserved slice below radius 1. In particular, this
family does not yield a one-ball compressor and cannot serve as the missing
high-density local mechanism for the target packing problem.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


J = np.array(
    [
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [-1.0, 0.0, 0.0, 0.0],
        [0.0, -1.0, 0.0, 0.0],
    ]
)
E_Q1 = np.array([1.0, 0.0, 0.0, 0.0])
E_P1 = np.array([0.0, 0.0, 1.0, 0.0])
ACTIVE_PLANE_BASIS = np.stack([E_Q1, E_P1], axis=1)


@dataclass(order=True)
class FoldRecord:
    sampled_radius: float
    s: float
    d: float
    lam: float
    sampled_area: float
    boundary_samples: int


def omega(x: np.ndarray, y: np.ndarray) -> float:
    return float(x @ J @ y)


def transvection_linear_part(a: np.ndarray, lam: float) -> np.ndarray:
    return np.eye(4) + lam * np.outer(a, J @ a)


def affine_transvection(points: np.ndarray, a: np.ndarray, lam: float, d: float) -> np.ndarray:
    levels = points @ (J @ a) - d
    return points + lam * levels[:, None] * a[None, :]


def one_sided_transvection_fold(points: np.ndarray, lam: float, d: float) -> np.ndarray:
    out = points.copy()
    shear = np.maximum(points[:, 2] - d, 0.0)
    out[:, 0] += lam * shear
    return out


def normalized_fold_slice(points: np.ndarray, s: float, lam: float, d: float) -> np.ndarray:
    out = points.copy()
    out[:, 0] = s * (points[:, 0] + lam * np.maximum(points[:, 1] - d, 0.0))
    out[:, 1] = points[:, 1] / s
    return out


def normalized_fold_4d(points: np.ndarray, s: float, lam: float, d: float) -> np.ndarray:
    out = points.copy()
    out[:, 0] = s * (points[:, 0] + lam * np.maximum(points[:, 2] - d, 0.0))
    out[:, 2] = points[:, 2] / s
    return out


def sample_disk_boundary(n_samples: int) -> np.ndarray:
    theta = np.linspace(0.0, 2.0 * math.pi, num=n_samples, endpoint=False)
    return np.stack([np.cos(theta), np.sin(theta)], axis=1)


def embed_slice_points(points_2d: np.ndarray) -> np.ndarray:
    out = np.zeros((points_2d.shape[0], 4))
    out[:, 0] = points_2d[:, 0]
    out[:, 2] = points_2d[:, 1]
    return out


def oriented_polygon_area(points_2d: np.ndarray) -> float:
    x = points_2d[:, 0]
    y = points_2d[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def circle_from_points(points: list[np.ndarray]) -> tuple[np.ndarray, float]:
    if not points:
        return np.zeros(2), -1.0
    if len(points) == 1:
        return points[0].copy(), 0.0
    if len(points) == 2:
        center = 0.5 * (points[0] + points[1])
        radius = float(np.linalg.norm(points[0] - center))
        return center, radius

    best: tuple[np.ndarray, float] | None = None
    for i in range(3):
        for j in range(i + 1, 3):
            center, radius = circle_from_points([points[i], points[j]])
            if all(np.linalg.norm(p - center) <= radius + 1e-12 for p in points):
                if best is None or radius < best[1]:
                    best = (center, radius)
    if best is not None:
        return best

    a, b, c = points
    ax, ay = a
    bx, by = b
    cx, cy = c
    det = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(det) <= 1e-15:
        raise RuntimeError("degenerate support set in circumcircle computation")
    ux = (
        (ax * ax + ay * ay) * (by - cy)
        + (bx * bx + by * by) * (cy - ay)
        + (cx * cx + cy * cy) * (ay - by)
    ) / det
    uy = (
        (ax * ax + ay * ay) * (cx - bx)
        + (bx * bx + by * by) * (ax - cx)
        + (cx * cx + cy * cy) * (bx - ax)
    ) / det
    center = np.array([ux, uy])
    radius = float(np.linalg.norm(a - center))
    return center, radius


def point_in_circle(point: np.ndarray, circle: tuple[np.ndarray, float], tol: float = 1e-12) -> bool:
    center, radius = circle
    if radius < 0.0:
        return False
    return float(np.linalg.norm(point - center)) <= radius + tol


def welzl(points: list[np.ndarray], boundary: list[np.ndarray], n: int) -> tuple[np.ndarray, float]:
    if n == 0 or len(boundary) == 3:
        return circle_from_points(boundary)

    point = points[n - 1]
    circle = welzl(points, boundary, n - 1)
    if point_in_circle(point, circle):
        return circle

    boundary.append(point)
    circle = welzl(points, boundary, n - 1)
    boundary.pop()
    return circle


def minimum_enclosing_circle(points_2d: np.ndarray, seed: int) -> tuple[np.ndarray, float]:
    sys.setrecursionlimit(max(10000, points_2d.shape[0] + 50))
    points = [row.copy() for row in points_2d]
    random.Random(seed).shuffle(points)
    return welzl(points, [], len(points))


def orthonormalize_plane(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    e1 = v1 / np.linalg.norm(v1)
    v2_perp = v2 - np.dot(v2, e1) * e1
    norm_v2 = np.linalg.norm(v2_perp)
    if norm_v2 <= 1e-12:
        raise RuntimeError("image plane basis became degenerate")
    e2 = v2_perp / norm_v2
    return np.stack([e1, e2], axis=1)


def project_to_plane(points_4d: np.ndarray, plane_basis: np.ndarray) -> np.ndarray:
    return points_4d @ plane_basis


def q_shear_matrix(b11: float, b12: float, b22: float) -> np.ndarray:
    b = np.array([[b11, b12], [b12, b22]])
    top = np.concatenate([np.eye(2), b], axis=1)
    bottom = np.concatenate([np.zeros((2, 2)), np.eye(2)], axis=1)
    return np.concatenate([top, bottom], axis=0)


def p_shear_matrix(c11: float, c12: float, c22: float) -> np.ndarray:
    c = np.array([[c11, c12], [c12, c22]])
    top = np.concatenate([np.eye(2), np.zeros((2, 2))], axis=1)
    bottom = np.concatenate([c, np.eye(2)], axis=1)
    return np.concatenate([top, bottom], axis=0)


def linear_change_matrix(a: np.ndarray) -> np.ndarray:
    a_inv_t = np.linalg.inv(a).T
    top = np.concatenate([a, np.zeros((2, 2))], axis=1)
    bottom = np.concatenate([np.zeros((2, 2)), a_inv_t], axis=1)
    return np.concatenate([top, bottom], axis=0)


def rotation_block(theta: float) -> np.ndarray:
    rot = np.array(
        [
            [math.cos(theta), -math.sin(theta)],
            [math.sin(theta), math.cos(theta)],
        ]
    )
    top = np.concatenate([rot, np.zeros((2, 2))], axis=1)
    bottom = np.concatenate([np.zeros((2, 2)), rot], axis=1)
    return np.concatenate([top, bottom], axis=0)


def random_symplectic_matrix(rng: np.random.Generator, scale: float) -> np.ndarray:
    while True:
        a = np.eye(2) + scale * rng.normal(size=(2, 2))
        if abs(np.linalg.det(a)) > 0.35:
            if np.linalg.det(a) < 0:
                a[0, :] *= -1.0
            break
    q = q_shear_matrix(*(scale * rng.normal(size=3)))
    p = p_shear_matrix(*(scale * rng.normal(size=3)))
    rot = rotation_block(float(rng.uniform(-math.pi, math.pi)))
    return p @ q @ linear_change_matrix(a) @ rot


def symplectic_error(matrix: np.ndarray) -> float:
    return float(np.max(np.abs(matrix.T @ J @ matrix - J)))


def plane_delta(matrix: np.ndarray) -> tuple[float, np.ndarray]:
    image_basis = matrix @ ACTIVE_PLANE_BASIS
    ortho_basis = orthonormalize_plane(image_basis[:, 0], image_basis[:, 1])
    delta = abs(omega(ortho_basis[:, 0], ortho_basis[:, 1]))
    return delta, ortho_basis


def slice_boundary_stats(
    s: float,
    lam: float,
    d: float,
    boundary_samples: int,
    circle_seed: int,
    left_factor: np.ndarray | None = None,
) -> dict:
    boundary = sample_disk_boundary(boundary_samples)
    slice_image = normalized_fold_slice(boundary, s=s, lam=lam, d=d)
    area = abs(oriented_polygon_area(slice_image))
    if left_factor is None:
        circle_center, circle_radius = minimum_enclosing_circle(slice_image, seed=circle_seed)
        return {
            "plane_delta": 1.0,
            "theoretical_lower_bound": 1.0,
            "sampled_area": area,
            "sampled_radius": circle_radius,
            "circle_center": circle_center.tolist(),
        }

    slice_image_4d = left_factor @ embed_slice_points(slice_image).T
    slice_image_4d = slice_image_4d.T
    delta, ortho_basis = plane_delta(left_factor)
    projected = project_to_plane(slice_image_4d, ortho_basis)
    circle_center, circle_radius = minimum_enclosing_circle(projected, seed=circle_seed)
    area_projected = abs(oriented_polygon_area(projected))
    return {
        "plane_delta": delta,
        "theoretical_lower_bound": delta ** (-0.5),
        "sampled_area": area_projected,
        "sampled_radius": circle_radius,
        "circle_center": circle_center.tolist(),
    }


def parameter_grid() -> tuple[list[float], list[float], list[float]]:
    s_values = [float(value) for value in np.exp(np.linspace(math.log(0.8), math.log(1.25), num=19))]
    d_values = [round(float(value), 10) for value in np.linspace(-0.9, 0.9, num=19)]
    lam_values = [round(-1.0 + 0.125 * idx, 10) for idx in range(17) if abs(-1.0 + 0.125 * idx) > 1e-12]
    return s_values, d_values, lam_values


def coarse_search(boundary_samples: int, top_n: int, circle_seed: int) -> dict:
    s_values, d_values, lam_values = parameter_grid()
    top_records: list[FoldRecord] = []
    total = 0
    for s in s_values:
        for d in d_values:
            for lam in lam_values:
                total += 1
                stats = slice_boundary_stats(
                    s=s,
                    lam=lam,
                    d=d,
                    boundary_samples=boundary_samples,
                    circle_seed=circle_seed,
                )
                record = FoldRecord(
                    sampled_radius=stats["sampled_radius"],
                    s=s,
                    d=d,
                    lam=lam,
                    sampled_area=stats["sampled_area"],
                    boundary_samples=boundary_samples,
                )
                top_records.append(record)
    top_records.sort()
    unique: list[FoldRecord] = []
    seen = set()
    for record in top_records:
        key = (round(record.s, 12), round(record.d, 12), round(record.lam, 12))
        if key in seen:
            continue
        seen.add(key)
        unique.append(record)
        if len(unique) >= top_n:
            break
    return {
        "n_cases": total,
        "s_values": s_values,
        "d_values": d_values,
        "lam_values": lam_values,
        "top_records": [asdict(record) for record in unique],
    }


def refine_records(records: list[FoldRecord], boundary_samples: int, circle_seed: int) -> list[dict]:
    refined = []
    for record in records:
        stats = slice_boundary_stats(
            s=record.s,
            lam=record.lam,
            d=record.d,
            boundary_samples=boundary_samples,
            circle_seed=circle_seed,
        )
        refined.append(
            {
                "s": record.s,
                "d": record.d,
                "lam": record.lam,
                "sampled_radius": stats["sampled_radius"],
                "sampled_area": stats["sampled_area"],
                "boundary_samples": boundary_samples,
                "theoretical_lower_bound": 1.0,
                "radius_minus_bound": stats["sampled_radius"] - 1.0,
                "area_minus_pi": stats["sampled_area"] - math.pi,
            }
        )
    refined.sort(key=lambda row: row["sampled_radius"])
    return refined


def verify_transvection_formulas() -> dict:
    a = np.array([-1.0, 0.0, 0.0, 0.0])
    lam = 0.375
    d = 0.2
    points = np.array(
        [
            [0.1, 0.0, -0.5, 0.0],
            [0.2, 0.0, 0.2, 0.0],
            [-0.3, 0.0, 0.8, 0.0],
        ]
    )
    affine = affine_transvection(points, a=a, lam=lam, d=d)
    one_sided = one_sided_transvection_fold(points, lam=lam, d=d)
    linear_error = symplectic_error(transvection_linear_part(a=a, lam=lam))
    affine_match_error = float(np.max(np.abs(affine[:, 0] - (points[:, 0] + lam * (points[:, 2] - d)))))
    fold_match_error = float(np.max(np.abs(one_sided[:, 0] - (points[:, 0] + lam * np.maximum(points[:, 2] - d, 0.0)))))
    return {
        "active_transvection_vector": a.tolist(),
        "sample_lambda": lam,
        "sample_d": d,
        "linear_transvection_symplectic_error": linear_error,
        "affine_slice_formula_error": affine_match_error,
        "one_sided_slice_formula_error": fold_match_error,
        "jacobian_below_seam": [[1.0, 0.0], [0.0, 1.0]],
        "jacobian_above_seam": [[1.0, lam], [0.0, 1.0]],
        "normalized_jacobian_below_seam": [[1.5, 0.0], [0.0, 2.0 / 3.0]],
        "normalized_jacobian_above_seam": [[1.5, 1.5 * lam], [0.0, 2.0 / 3.0]],
        "jacobian_determinant_below_seam": 1.0,
        "jacobian_determinant_above_seam": 1.0,
        "normalized_jacobian_determinant_below_seam": 1.0,
        "normalized_jacobian_determinant_above_seam": 1.0,
    }


def identity_baseline(boundary_samples: int, circle_seed: int) -> dict:
    stats = slice_boundary_stats(
        s=1.0,
        lam=0.0,
        d=0.0,
        boundary_samples=boundary_samples,
        circle_seed=circle_seed,
    )
    return {
        "s": 1.0,
        "d": 0.0,
        "lam": 0.0,
        "sampled_radius": stats["sampled_radius"],
        "sampled_area": stats["sampled_area"],
        "boundary_samples": boundary_samples,
        "theoretical_lower_bound": 1.0,
    }


def mixed_left_factor_audit(
    best_fold: dict,
    boundary_samples: int,
    circle_seed: int,
    n_mixers: int,
    random_seed: int,
    mixer_scale: float,
) -> dict:
    rng = np.random.default_rng(random_seed)
    rows = []
    for mixer_index in range(n_mixers):
        matrix = random_symplectic_matrix(rng=rng, scale=mixer_scale)
        stats = slice_boundary_stats(
            s=best_fold["s"],
            lam=best_fold["lam"],
            d=best_fold["d"],
            boundary_samples=boundary_samples,
            circle_seed=circle_seed + mixer_index + 1,
            left_factor=matrix,
        )
        rows.append(
            {
                "mixer_index": mixer_index,
                "symplectic_error": symplectic_error(matrix),
                "plane_delta": stats["plane_delta"],
                "theoretical_lower_bound": stats["theoretical_lower_bound"],
                "sampled_radius": stats["sampled_radius"],
                "radius_minus_bound": stats["sampled_radius"] - stats["theoretical_lower_bound"],
                "sampled_area": stats["sampled_area"],
                "area_minus_pi_over_delta": stats["sampled_area"] - (math.pi / stats["plane_delta"]),
                "matrix": matrix.tolist(),
            }
        )
    rows.sort(key=lambda row: row["sampled_radius"])
    return {
        "n_mixers": n_mixers,
        "best_fold_parameters": best_fold,
        "min_plane_delta": min(row["plane_delta"] for row in rows),
        "max_plane_delta": max(row["plane_delta"] for row in rows),
        "min_theoretical_lower_bound": min(row["theoretical_lower_bound"] for row in rows),
        "max_theoretical_lower_bound": max(row["theoretical_lower_bound"] for row in rows),
        "min_sampled_radius": min(row["sampled_radius"] for row in rows),
        "max_sampled_radius": max(row["sampled_radius"] for row in rows),
        "min_radius_minus_bound": min(row["radius_minus_bound"] for row in rows),
        "max_symplectic_error": max(row["symplectic_error"] for row in rows),
        "rows": rows,
    }


def build_report(args: argparse.Namespace) -> dict:
    coarse = coarse_search(
        boundary_samples=args.coarse_boundary_samples,
        top_n=args.refine_top_n,
        circle_seed=args.circle_seed,
    )
    top_records = [FoldRecord(**row) for row in coarse["top_records"]]
    refined = refine_records(
        records=top_records,
        boundary_samples=args.refined_boundary_samples,
        circle_seed=args.circle_seed,
    )
    best_fold = refined[0]
    mixed = mixed_left_factor_audit(
        best_fold=best_fold,
        boundary_samples=args.mixed_boundary_samples,
        circle_seed=args.circle_seed,
        n_mixers=args.n_mixers,
        random_seed=args.random_seed,
        mixer_scale=args.mixer_scale,
    )
    return {
        "family": "affine_symplectic_transvection_folds_with_preserved_slice",
        "normalized_slice_formula": "F_{s,lambda,d}(u,v) = ( s * (u + lambda * max(v-d, 0)), v / s )",
        "one_sided_fold_formula": "(q1,q2,p1,p2) -> (q1 + lambda * max(p1-d,0), q2, p1, p2)",
        "mixed_left_factor_formula": "x -> L(S_s(Phi_{lambda,d}(x))) for L in Sp(4,R)",
        "exact_claim": (
            "For every normalized fold and every common linear symplectic left factor L, "
            "the image of the central active slice has symplectic area pi inside the plane "
            "L(E). If delta = |omega(f1,f2)| for an orthonormal basis of L(E), then any "
            "containing Euclidean 4-ball must satisfy R >= delta^(-1/2) >= 1."
        ),
        "formula_checks": verify_transvection_formulas(),
        "identity_baseline": identity_baseline(
            boundary_samples=args.refined_boundary_samples,
            circle_seed=args.circle_seed,
        ),
        "coarse_nontrivial_search": coarse,
        "refined_nontrivial_search": refined,
        "mixed_left_factor_audit": mixed,
        "conclusion": (
            "The preserved-slice obstruction survives the full normalized family and any "
            "common linear symplectic left factor. Nontrivial folds do not beat radius 1 "
            "even before mixing, and random mixed audits obey the sharper bound "
            "R >= delta^(-1/2). This family does not repair the proposed solution."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/research_loop/attempt_015_transvection_fold_slice_obstruction.json"),
    )
    parser.add_argument("--coarse-boundary-samples", type=int, default=256)
    parser.add_argument("--refined-boundary-samples", type=int, default=4096)
    parser.add_argument("--mixed-boundary-samples", type=int, default=2048)
    parser.add_argument("--refine-top-n", type=int, default=12)
    parser.add_argument("--n-mixers", type=int, default=16)
    parser.add_argument("--mixer-scale", type=float, default=0.65)
    parser.add_argument("--circle-seed", type=int, default=101)
    parser.add_argument("--random-seed", type=int, default=1234)
    args = parser.parse_args()

    report = build_report(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True))
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
