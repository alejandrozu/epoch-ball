#!/usr/bin/env python3
"""Attempt-017 search over cotangent-lift folding Hamiltonians.

This family uses smooth Hamiltonians of the form

    H(q, p, t) = <p, X(q)>

with a smooth base-plane vector field X on R^2. The time-1 Hamiltonian flow is
the cotangent lift of the time-1 flow of X on the q-plane:

    q_dot = X(q),
    p_dot = -(DX(q))^T p.

The concrete attempt-017 ansatz bends the long source strip in the q-plane
toward a spiral ribbon. This is inspired by the symplectic-folding line from
the literature, but kept in a directly Hamiltonian and easily auditable form.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class SpiralSpec:
    long_scale: float
    twist: float
    radial_offset: float
    radial_growth: float
    normal_scale: float
    strength: float
    support_x: float
    support_y: float
    x_start: float
    x_center: float


def sample_boundary_points(k: int, n_pts_per_ball: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    clouds = []
    for i in range(k):
        pts = rng.normal(size=(n_pts_per_ball, 4))
        pts /= np.linalg.norm(pts, axis=1, keepdims=True)
        pts[:, 0] += 3.0 * i
        clouds.append(pts)
    return np.concatenate(clouds, axis=0)


def ritter_enclosing_ball(points: np.ndarray) -> tuple[np.ndarray, float]:
    idx0 = 0
    p0 = points[idx0]
    idx1 = int(np.argmax(np.sum((points - p0) ** 2, axis=1)))
    p1 = points[idx1]
    idx2 = int(np.argmax(np.sum((points - p1) ** 2, axis=1)))
    p2 = points[idx2]
    center = 0.5 * (p1 + p2)
    radius = float(np.linalg.norm(p1 - center))

    for point in points:
        diff = point - center
        dist = float(np.linalg.norm(diff))
        if dist <= radius + 1e-12:
            continue
        new_radius = 0.5 * (radius + dist)
        center = center + ((new_radius - radius) / dist) * diff
        radius = new_radius

    return center, radius


def spiral_target_and_derivatives(q: np.ndarray, spec: SpiralSpec) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = q[:, 0]
    y = q[:, 1]
    u = spec.long_scale * (x - spec.x_start)
    theta = spec.twist * u
    radius = spec.radial_offset + spec.radial_growth * u

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    centerline = np.stack([radius * cos_theta, radius * sin_theta], axis=1)

    v = np.stack(
        [
            spec.radial_growth * cos_theta - radius * spec.twist * sin_theta,
            spec.radial_growth * sin_theta + radius * spec.twist * cos_theta,
        ],
        axis=1,
    )
    speed = np.linalg.norm(v, axis=1, keepdims=True)
    speed = np.maximum(speed, 1e-12)
    tangent = v / speed
    normal = np.stack([-tangent[:, 1], tangent[:, 0]], axis=1)

    accel = np.stack(
        [
            -2.0 * spec.radial_growth * spec.twist * sin_theta - radius * (spec.twist**2) * cos_theta,
            2.0 * spec.radial_growth * spec.twist * cos_theta - radius * (spec.twist**2) * sin_theta,
        ],
        axis=1,
    )
    tangent_accel = np.sum(tangent * accel, axis=1, keepdims=True)
    tangent_u_deriv = (accel - tangent_accel * tangent) / speed
    normal_u_deriv = np.stack([-tangent_u_deriv[:, 1], tangent_u_deriv[:, 0]], axis=1)

    target = centerline + spec.normal_scale * y[:, None] * normal
    dtarget_dx = spec.long_scale * (v + spec.normal_scale * y[:, None] * normal_u_deriv)
    dtarget_dy = spec.normal_scale * normal
    return target, dtarget_dx, dtarget_dy


def cutoff_and_gradients(q: np.ndarray, spec: SpiralSpec) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = q[:, 0]
    y = q[:, 1]
    sx = (x - spec.x_center) / spec.support_x
    sy = y / spec.support_y
    exponent = -(sx**8) - (sy**8)
    cutoff = np.exp(exponent)
    dcutoff_dx = cutoff * (-8.0 * (sx**7) / spec.support_x)
    dcutoff_dy = cutoff * (-8.0 * (sy**7) / spec.support_y)
    return cutoff, dcutoff_dx, dcutoff_dy


def vector_field_and_jacobian(q: np.ndarray, spec: SpiralSpec) -> tuple[np.ndarray, np.ndarray]:
    target, dtarget_dx, dtarget_dy = spiral_target_and_derivatives(q, spec)
    cutoff, dcutoff_dx, dcutoff_dy = cutoff_and_gradients(q, spec)
    delta = target - q
    field = spec.strength * cutoff[:, None] * delta

    jac = np.zeros((q.shape[0], 2, 2), dtype=np.float64)
    ex = np.array([1.0, 0.0], dtype=np.float64)
    ey = np.array([0.0, 1.0], dtype=np.float64)
    jac[:, :, 0] = spec.strength * (
        dcutoff_dx[:, None] * delta + cutoff[:, None] * (dtarget_dx - ex[None, :])
    )
    jac[:, :, 1] = spec.strength * (
        dcutoff_dy[:, None] * delta + cutoff[:, None] * (dtarget_dy - ey[None, :])
    )
    return field, jac


def rhs(q: np.ndarray, p: np.ndarray, spec: SpiralSpec) -> tuple[np.ndarray, np.ndarray]:
    field, jac = vector_field_and_jacobian(q, spec)
    p_dot = -np.einsum("nji,nj->ni", jac, p)
    return field, p_dot


def rk4_cotangent_flow(q: np.ndarray, p: np.ndarray, spec: SpiralSpec, n_steps: int) -> tuple[np.ndarray, np.ndarray]:
    dt = 1.0 / float(n_steps)
    q_cur = q.copy()
    p_cur = p.copy()
    for _ in range(n_steps):
        k1q, k1p = rhs(q_cur, p_cur, spec)
        k2q, k2p = rhs(q_cur + 0.5 * dt * k1q, p_cur + 0.5 * dt * k1p, spec)
        k3q, k3p = rhs(q_cur + 0.5 * dt * k2q, p_cur + 0.5 * dt * k2p, spec)
        k4q, k4p = rhs(q_cur + dt * k3q, p_cur + dt * k3p, spec)
        q_cur = q_cur + (dt / 6.0) * (k1q + 2.0 * k2q + 2.0 * k3q + k4q)
        p_cur = p_cur + (dt / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p)
    return q_cur, p_cur


def evaluate_points(points: np.ndarray, spec: SpiralSpec, n_steps: int) -> dict:
    q0 = points[:, :2]
    p0 = points[:, 2:]
    q1, p1 = rk4_cotangent_flow(q0, p0, spec, n_steps=n_steps)
    mapped = np.concatenate([q1, p1], axis=1)
    center, radius = ritter_enclosing_ball(mapped)
    return {
        "center": center.tolist(),
        "radius": radius,
        "ratio": 0.0 if radius <= 0.0 else None,
        "mapped_q_mean": q1.mean(axis=0).tolist(),
        "mapped_p_mean": p1.mean(axis=0).tolist(),
    }


def build_specs(k: int) -> list[SpiralSpec]:
    x_start = -1.0
    x_center = 1.5 * 3.0 * (k - 1)
    specs: list[SpiralSpec] = []
    for long_scale in [0.2, 0.3, 0.4, 0.5]:
        for twist in [0.8, 1.2, 1.6, 2.0]:
            for radial_offset in [0.8, 1.2, 1.6]:
                for radial_growth in [0.0, 0.08, 0.16]:
                    for normal_scale in [0.5, 0.75, 1.0]:
                        for strength in [0.8, 1.2]:
                            specs.append(
                                SpiralSpec(
                                    long_scale=long_scale,
                                    twist=twist,
                                    radial_offset=radial_offset,
                                    radial_growth=radial_growth,
                                    normal_scale=normal_scale,
                                    strength=strength,
                                    support_x=18.0,
                                    support_y=3.5,
                                    x_start=x_start,
                                    x_center=x_center,
                                )
                            )
    return specs


def coarse_search(k: int, n_pts_per_ball: int, seed: int, n_steps: int, keep_top_n: int) -> dict:
    points = sample_boundary_points(k=k, n_pts_per_ball=n_pts_per_ball, seed=seed)
    rows = []
    for spec in build_specs(k):
        eval_row = evaluate_points(points=points, spec=spec, n_steps=n_steps)
        radius = eval_row["radius"]
        ratio = 0.0 if radius <= 0.0 else k / (radius**4)
        rows.append(
            {
                "spec": asdict(spec),
                "coarse_center": eval_row["center"],
                "coarse_radius": radius,
                "coarse_ratio": ratio,
                "mapped_q_mean": eval_row["mapped_q_mean"],
                "mapped_p_mean": eval_row["mapped_p_mean"],
            }
        )
    rows.sort(key=lambda row: row["coarse_radius"])
    return {
        "search_settings": {
            "candidate_count": len(rows),
            "n_pts_per_ball": n_pts_per_ball,
            "seed": seed,
            "n_steps": n_steps,
            "keep_top_n": keep_top_n,
        },
        "top_candidates": rows[:keep_top_n],
    }


def refine_candidates(top_candidates: list[dict], k: int, n_pts_per_ball: int, seeds: list[int], n_steps: int) -> list[dict]:
    refined = []
    for row in top_candidates:
        spec = SpiralSpec(**row["spec"])
        seed_rows = []
        radii = []
        for seed in seeds:
            points = sample_boundary_points(k=k, n_pts_per_ball=n_pts_per_ball, seed=seed)
            eval_row = evaluate_points(points=points, spec=spec, n_steps=n_steps)
            radius = eval_row["radius"]
            ratio = k / (radius**4)
            radii.append(radius)
            seed_rows.append(
                {
                    "seed": seed,
                    "radius": radius,
                    "ratio": ratio,
                    "center": eval_row["center"],
                    "mapped_q_mean": eval_row["mapped_q_mean"],
                    "mapped_p_mean": eval_row["mapped_p_mean"],
                }
            )
        refined.append(
            {
                "spec": row["spec"],
                "coarse_radius": row["coarse_radius"],
                "coarse_ratio": row["coarse_ratio"],
                "refine_seed_rows": seed_rows,
                "best_radius": min(radii),
                "best_ratio": k / (min(radii) ** 4),
                "worst_radius": max(radii),
                "worst_ratio": k / (max(radii) ** 4),
                "mean_radius": float(np.mean(radii)),
                "mean_ratio": float(np.mean([k / (radius**4) for radius in radii])),
            }
        )
    refined.sort(key=lambda row: row["worst_radius"])
    return refined


def render_hamiltonian(spec: SpiralSpec) -> str:
    return f"""import torch

torch.set_default_dtype(torch.float64)

_LONG_SCALE = {spec.long_scale!r}
_TWIST = {spec.twist!r}
_RADIAL_OFFSET = {spec.radial_offset!r}
_RADIAL_GROWTH = {spec.radial_growth!r}
_NORMAL_SCALE = {spec.normal_scale!r}
_STRENGTH = {spec.strength!r}
_SUPPORT_X = {spec.support_x!r}
_SUPPORT_Y = {spec.support_y!r}
_X_START = {spec.x_start!r}
_X_CENTER = {spec.x_center!r}


def _vector_field(Q):
    x = Q[:, 0]
    y = Q[:, 1]
    u = _LONG_SCALE * (x - _X_START)
    theta = _TWIST * u
    radius = _RADIAL_OFFSET + _RADIAL_GROWTH * u

    cos_theta = torch.cos(theta)
    sin_theta = torch.sin(theta)

    centerline = torch.stack([radius * cos_theta, radius * sin_theta], dim=1)

    v = torch.stack(
        [
            _RADIAL_GROWTH * cos_theta - radius * _TWIST * sin_theta,
            _RADIAL_GROWTH * sin_theta + radius * _TWIST * cos_theta,
        ],
        dim=1,
    )
    speed = torch.linalg.norm(v, dim=1, keepdim=True).clamp_min(1e-12)
    tangent = v / speed
    normal = torch.stack([-tangent[:, 1], tangent[:, 0]], dim=1)

    accel = torch.stack(
        [
            -2.0 * _RADIAL_GROWTH * _TWIST * sin_theta - radius * (_TWIST ** 2) * cos_theta,
            2.0 * _RADIAL_GROWTH * _TWIST * cos_theta - radius * (_TWIST ** 2) * sin_theta,
        ],
        dim=1,
    )
    tangent_accel = (tangent * accel).sum(dim=1, keepdim=True)
    tangent_u_deriv = (accel - tangent_accel * tangent) / speed
    normal_u_deriv = torch.stack([-tangent_u_deriv[:, 1], tangent_u_deriv[:, 0]], dim=1)

    target = centerline + _NORMAL_SCALE * y[:, None] * normal
    sx = (x - _X_CENTER) / _SUPPORT_X
    sy = y / _SUPPORT_Y
    cutoff = torch.exp(-(sx.pow(8)) - (sy.pow(8)))
    return _STRENGTH * cutoff[:, None] * (target - Q)


def Hamiltonian(Q, P, t):
    return (P * _vector_field(Q)).sum(dim=1)
"""


def build_report(args: argparse.Namespace) -> dict:
    coarse = coarse_search(
        k=args.k,
        n_pts_per_ball=args.coarse_n_pts_per_ball,
        seed=args.coarse_seed,
        n_steps=args.coarse_n_steps,
        keep_top_n=args.refine_top_n,
    )
    refined = refine_candidates(
        top_candidates=coarse["top_candidates"],
        k=args.k,
        n_pts_per_ball=args.refined_n_pts_per_ball,
        seeds=args.refine_seeds,
        n_steps=args.refined_n_steps,
    )
    best = refined[0]
    best_spec = SpiralSpec(**best["spec"])
    dense = refine_candidates(
        top_candidates=[{"spec": best["spec"], "coarse_radius": best["coarse_radius"], "coarse_ratio": best["coarse_ratio"]}],
        k=args.k,
        n_pts_per_ball=args.dense_n_pts_per_ball,
        seeds=args.dense_seeds,
        n_steps=args.dense_n_steps,
    )[0]
    return {
        "family": "cotangent_lift_spiral_fold_hamiltonian",
        "hamiltonian_form": "H(q,p,t) = <p, X(q)>",
        "vector_field_form": "X(q) = strength * cutoff(q) * (target_spiral_ribbon(q) - q)",
        "coarse_search": coarse,
        "refined_candidates": refined,
        "best_candidate": best,
        "dense_validation": dense,
        "generated_hamiltonian": {
            "path": str(args.out_hamiltonian),
            "spec": asdict(best_spec),
        },
        "conclusion": (
            "This attempt tests a smooth Hamiltonian family whose flow is a cotangent lift of a "
            "folding vector field on the q-plane. The saved candidate is exact and smooth, but "
            "still requires numerical audit to judge packing quality."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--coarse-n-pts-per-ball", type=int, default=64)
    parser.add_argument("--refined-n-pts-per-ball", type=int, default=256)
    parser.add_argument("--dense-n-pts-per-ball", type=int, default=1024)
    parser.add_argument("--coarse-seed", type=int, default=123)
    parser.add_argument("--refine-seeds", type=int, nargs="+", default=[11, 29, 47])
    parser.add_argument("--dense-seeds", type=int, nargs="+", default=[11, 29, 47, 59, 83])
    parser.add_argument("--coarse-n-steps", type=int, default=64)
    parser.add_argument("--refined-n-steps", type=int, default=128)
    parser.add_argument("--dense-n-steps", type=int, default=256)
    parser.add_argument("--refine-top-n", type=int, default=8)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/research_loop/attempt_017_cotangent_lift_folding_search.json"),
    )
    parser.add_argument("--out-hamiltonian", type=Path, default=Path("Hamiltonian.py"))
    args = parser.parse_args()

    report = build_report(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True))

    best_spec = SpiralSpec(**report["best_candidate"]["spec"])
    args.out_hamiltonian.write_text(render_hamiltonian(best_spec))
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
