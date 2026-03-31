#!/usr/bin/env python3
"""Attempt-019 search over localized quadratic squeeze/rotation gate trees.

This family replaces attempt 018's control-translation gates by compactly
supported quadratic Hamiltonians whose linear part is a metaplectic-style local
gate:

    H = chi * [ lam * (z * p1 + eta * q2 * p2) + theta * (z * p2 - q2 * p1) ]

where z = q1 - q_mid for a pair midpoint. The `lam` term is a local squeeze and
the `theta` term rotates the active and passive planes together. The supports
are arranged in a recursive pair tree.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


WINDOW_NORM = 0.007029858406609658


@dataclass(frozen=True)
class GateSpec:
    lam0: float
    lam_growth: float
    eta: float
    theta0: float
    theta_growth: float
    support_scale: float
    support_r: float
    support_p: float


@dataclass(frozen=True)
class PairNode:
    level: int
    q_mid: float
    support_q: float


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


def bump_unit_interval(u: float) -> float:
    if not (0.0 < u < 1.0):
        return 0.0
    return math.exp(-1.0 / (u * (1.0 - u))) / WINDOW_NORM


def stage_weight(t: float, n_stages: int) -> tuple[int, float]:
    if not (0.0 < t < 1.0):
        return -1, 0.0
    s = t * n_stages
    stage = int(math.floor(s))
    if stage >= n_stages:
        stage = n_stages - 1
    u = s - stage
    return stage, n_stages * bump_unit_interval(u)


def sample_specs() -> list[GateSpec]:
    specs: list[GateSpec] = []
    for lam0 in [-0.35, -0.6, -0.9]:
        for lam_growth in [1.0, 1.25]:
            for eta in [0.0, 1.0]:
                for theta0 in [0.0, 0.35, 0.7]:
                    for theta_growth in [1.0, 1.3]:
                        for support_scale in [1.1, 1.4]:
                            for support_r in [1.5]:
                                for support_p in [1.8, 2.6]:
                                    specs.append(
                                        GateSpec(
                                            lam0=lam0,
                                            lam_growth=lam_growth,
                                            eta=eta,
                                            theta0=theta0,
                                            theta_growth=theta_growth,
                                            support_scale=support_scale,
                                            support_r=support_r,
                                            support_p=support_p,
                                        )
                                    )
    return specs


def build_level_groups(k: int, spec: GateSpec) -> list[list[PairNode]]:
    clusters = [{"center_x": 3.0 * i, "size": 1} for i in range(k)]
    level_groups: list[list[PairNode]] = []
    level = 0
    while len(clusters) > 1:
        group: list[PairNode] = []
        next_clusters = []
        i = 0
        while i < len(clusters):
            if i + 1 >= len(clusters):
                next_clusters.append(clusters[i])
                i += 1
                continue
            left = clusters[i]
            right = clusters[i + 1]
            half_gap = 0.5 * (right["center_x"] - left["center_x"])
            group.append(
                PairNode(
                    level=level,
                    q_mid=0.5 * (left["center_x"] + right["center_x"]),
                    support_q=spec.support_scale * half_gap + 1.25,
                )
            )
            next_clusters.append(
                {
                    "center_x": 0.5 * (left["center_x"] + right["center_x"]),
                    "size": left["size"] + right["size"],
                }
            )
            i += 2
        level_groups.append(group)
        clusters = next_clusters
        level += 1
    return level_groups


def level_params(level: int, spec: GateSpec) -> tuple[float, float]:
    lam = spec.lam0 * (spec.lam_growth**level)
    theta = spec.theta0 * ((-1.0) ** level) * (spec.theta_growth**level)
    return lam, theta


def rhs_level(
    q1: np.ndarray,
    q2: np.ndarray,
    p1: np.ndarray,
    p2: np.ndarray,
    nodes: list[PairNode],
    spec: GateSpec,
    weight: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    q1_dot = np.zeros_like(q1)
    q2_dot = np.zeros_like(q2)
    p1_dot = np.zeros_like(p1)
    p2_dot = np.zeros_like(p2)

    if not nodes:
        return q1_dot, q2_dot, p1_dot, p2_dot

    lam, theta = level_params(nodes[0].level, spec)
    for node in nodes:
        z = q1 - node.q_mid
        u = z / node.support_q
        v = q2 / spec.support_r
        w1 = p1 / spec.support_p
        w2 = p2 / spec.support_p
        cutoff = np.exp(-(u**4) - (v**4) - (w1**4) - (w2**4))

        G = lam * (z * p1 + spec.eta * q2 * p2) + theta * (z * p2 - q2 * p1)

        dcut_q1 = cutoff * (-4.0 * (u**3) / node.support_q)
        dcut_q2 = cutoff * (-4.0 * (v**3) / spec.support_r)
        dcut_p1 = cutoff * (-4.0 * (w1**3) / spec.support_p)
        dcut_p2 = cutoff * (-4.0 * (w2**3) / spec.support_p)

        dG_p1 = lam * z - theta * q2
        dG_p2 = lam * spec.eta * q2 + theta * z
        dG_q1 = lam * p1 + theta * p2
        dG_q2 = lam * spec.eta * p2 - theta * p1

        q1_dot += weight * (dcut_p1 * G + cutoff * dG_p1)
        q2_dot += weight * (dcut_p2 * G + cutoff * dG_p2)
        p1_dot += -weight * (dcut_q1 * G + cutoff * dG_q1)
        p2_dot += -weight * (dcut_q2 * G + cutoff * dG_q2)

    return q1_dot, q2_dot, p1_dot, p2_dot


def rhs(q: np.ndarray, p: np.ndarray, t: float, level_groups: list[list[PairNode]], spec: GateSpec) -> tuple[np.ndarray, np.ndarray]:
    n_main_stages = len(level_groups)
    stage, weight = stage_weight(t, n_main_stages)
    if weight == 0.0:
        return np.zeros_like(q), np.zeros_like(p)

    q1 = q[:, 0]
    q2 = q[:, 1]
    p1 = p[:, 0]
    p2 = p[:, 1]
    q1_dot, q2_dot, p1_dot, p2_dot = rhs_level(q1, q2, p1, p2, level_groups[stage], spec, weight)
    q_dot = np.stack([q1_dot, q2_dot], axis=1)
    p_dot = np.stack([p1_dot, p2_dot], axis=1)
    return q_dot, p_dot


def rk4_flow(q: np.ndarray, p: np.ndarray, level_groups: list[list[PairNode]], spec: GateSpec, n_steps: int) -> tuple[np.ndarray, np.ndarray]:
    dt = 1.0 / float(n_steps)
    q_cur = q.copy()
    p_cur = p.copy()
    for step in range(n_steps):
        t0 = step * dt
        k1q, k1p = rhs(q_cur, p_cur, t0, level_groups, spec)
        k2q, k2p = rhs(q_cur + 0.5 * dt * k1q, p_cur + 0.5 * dt * k1p, t0 + 0.5 * dt, level_groups, spec)
        k3q, k3p = rhs(q_cur + 0.5 * dt * k2q, p_cur + 0.5 * dt * k2p, t0 + 0.5 * dt, level_groups, spec)
        k4q, k4p = rhs(q_cur + dt * k3q, p_cur + dt * k3p, t0 + dt, level_groups, spec)
        q_cur = q_cur + (dt / 6.0) * (k1q + 2.0 * k2q + 2.0 * k3q + k4q)
        p_cur = p_cur + (dt / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p)
    return q_cur, p_cur


def evaluate_points(points: np.ndarray, level_groups: list[list[PairNode]], spec: GateSpec, n_steps: int) -> dict:
    q0 = points[:, :2]
    p0 = points[:, 2:]
    q1, p1 = rk4_flow(q0, p0, level_groups, spec, n_steps=n_steps)
    mapped = np.concatenate([q1, p1], axis=1)
    center, radius = ritter_enclosing_ball(mapped)
    return {
        "center": center.tolist(),
        "radius": radius,
        "mapped_q_mean": q1.mean(axis=0).tolist(),
        "mapped_p_mean": p1.mean(axis=0).tolist(),
    }


def level_payload(level_groups: list[list[PairNode]], spec: GateSpec) -> tuple[list[list[dict]], list[dict]]:
    groups = [[asdict(node) for node in group] for group in level_groups]
    params = []
    for level in range(len(level_groups)):
        lam, theta = level_params(level, spec)
        params.append({"lam": lam, "theta": theta})
    return groups, params


def render_hamiltonian(level_groups: list[list[PairNode]], spec: GateSpec, translation_shift: list[float]) -> str:
    groups, params = level_payload(level_groups, spec)
    return f"""import math
import torch

torch.set_default_dtype(torch.float64)

_WINDOW_NORM = {WINDOW_NORM!r}
_MAIN_LEVELS = {groups!r}
_LEVEL_PARAMS = {params!r}
_ETA = {spec.eta!r}
_SUPPORT_R = {spec.support_r!r}
_SUPPORT_P = {spec.support_p!r}
_TRANSLATE = torch.tensor({translation_shift!r}, dtype=torch.float64)


def _bump_unit_interval(u):
    if not (0.0 < u < 1.0):
        return 0.0
    return math.exp(-1.0 / (u * (1.0 - u))) / _WINDOW_NORM


def _stage_weight(t, n_stages):
    if not (0.0 < t < 1.0):
        return -1, 0.0
    s = t * n_stages
    stage = int(math.floor(s))
    if stage >= n_stages:
        stage = n_stages - 1
    u = s - stage
    return stage, n_stages * _bump_unit_interval(u)


def _level_hamiltonian(Q, P, nodes, params):
    q1 = Q[:, 0]
    q2 = Q[:, 1]
    p1 = P[:, 0]
    p2 = P[:, 1]
    lam = params["lam"]
    theta = params["theta"]
    H = 0.0 * (q1 + q2 + p1 + p2)
    for node in nodes:
        z = q1 - node["q_mid"]
        u = z / node["support_q"]
        v = q2 / _SUPPORT_R
        w1 = p1 / _SUPPORT_P
        w2 = p2 / _SUPPORT_P
        cutoff = torch.exp(-(u.pow(4)) - (v.pow(4)) - (w1.pow(4)) - (w2.pow(4)))
        gate = lam * (z * p1 + _ETA * q2 * p2) + theta * (z * p2 - q2 * p1)
        H = H + cutoff * gate
    return H


def Hamiltonian(Q, P, t):
    n_main_stages = len(_MAIN_LEVELS)
    total_stages = n_main_stages + 1
    stage, weight = _stage_weight(t, total_stages)
    H = 0.0 * (Q[:, 0] + P[:, 0])
    if weight == 0.0:
        return H
    if stage < n_main_stages:
        return weight * _level_hamiltonian(Q, P, _MAIN_LEVELS[stage], _LEVEL_PARAMS[stage])
    dq = _TRANSLATE[:2].to(Q.device)
    dp = _TRANSLATE[2:].to(P.device)
    return weight * ((P * dq).sum(dim=1) - (Q * dp).sum(dim=1))
"""


def coarse_search(k: int, n_pts_per_ball: int, seed: int, n_steps: int, keep_top_n: int) -> dict:
    points = sample_boundary_points(k=k, n_pts_per_ball=n_pts_per_ball, seed=seed)
    rows = []
    for spec in sample_specs():
        level_groups = build_level_groups(k, spec)
        eval_row = evaluate_points(points, level_groups, spec, n_steps=n_steps)
        radius = eval_row["radius"]
        rows.append(
            {
                "spec": asdict(spec),
                "node_counts_by_level": [len(group) for group in level_groups],
                "coarse_center": eval_row["center"],
                "coarse_radius": radius,
                "coarse_ratio": k / (radius**4),
                "mapped_q_mean": eval_row["mapped_q_mean"],
                "mapped_p_mean": eval_row["mapped_p_mean"],
            }
        )
    rows.sort(key=lambda row: row["coarse_radius"])
    best_no_rotation = next(row for row in rows if abs(row["spec"]["theta0"]) < 1e-12)
    return {
        "search_settings": {
            "candidate_count": len(rows),
            "n_pts_per_ball": n_pts_per_ball,
            "seed": seed,
            "n_steps": n_steps,
            "keep_top_n": keep_top_n,
        },
        "top_candidates": rows[:keep_top_n],
        "best_no_rotation_candidate": best_no_rotation,
    }


def refine_candidates(candidates: list[dict], k: int, n_pts_per_ball: int, seeds: list[int], n_steps: int) -> list[dict]:
    refined = []
    for row in candidates:
        spec = GateSpec(**row["spec"])
        level_groups = build_level_groups(k, spec)
        seed_rows = []
        radii = []
        for seed in seeds:
            points = sample_boundary_points(k=k, n_pts_per_ball=n_pts_per_ball, seed=seed)
            eval_row = evaluate_points(points, level_groups, spec, n_steps=n_steps)
            radius = eval_row["radius"]
            radii.append(radius)
            seed_rows.append(
                {
                    "seed": seed,
                    "radius": radius,
                    "ratio": k / (radius**4),
                    "center": eval_row["center"],
                    "mapped_q_mean": eval_row["mapped_q_mean"],
                    "mapped_p_mean": eval_row["mapped_p_mean"],
                }
            )
        refined.append(
            {
                "spec": row["spec"],
                "node_counts_by_level": row["node_counts_by_level"],
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


def dense_evaluation(spec: GateSpec, k: int, n_pts_per_ball: int, seeds: list[int], n_steps: int) -> dict:
    level_groups = build_level_groups(k, spec)
    centers = []
    rows = []
    radii = []
    for seed in seeds:
        points = sample_boundary_points(k=k, n_pts_per_ball=n_pts_per_ball, seed=seed)
        eval_row = evaluate_points(points, level_groups, spec, n_steps=n_steps)
        radius = eval_row["radius"]
        radii.append(radius)
        centers.append(np.array(eval_row["center"], dtype=np.float64))
        rows.append(
            {
                "seed": seed,
                "radius": radius,
                "ratio": k / (radius**4),
                "center": eval_row["center"],
                "mapped_q_mean": eval_row["mapped_q_mean"],
                "mapped_p_mean": eval_row["mapped_p_mean"],
            }
        )
    translation_shift = (-np.mean(np.stack(centers, axis=0), axis=0)).tolist()
    return {
        "spec": asdict(spec),
        "node_counts_by_level": [len(group) for group in level_groups],
        "seed_rows": rows,
        "best_radius": min(radii),
        "best_ratio": k / (min(radii) ** 4),
        "worst_radius": max(radii),
        "worst_ratio": k / (max(radii) ** 4),
        "mean_radius": float(np.mean(radii)),
        "mean_ratio": float(np.mean([k / (radius**4) for radius in radii])),
        "translation_shift": translation_shift,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--eps", type=float, default=0.1)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/research_loop/attempt_019_quadratic_gate_search.json"),
    )
    parser.add_argument("--emit", type=Path, default=Path("Hamiltonian.py"))
    args = parser.parse_args()

    coarse = coarse_search(k=args.k, n_pts_per_ball=64, seed=11, n_steps=96, keep_top_n=10)
    refined = refine_candidates(
        coarse["top_candidates"],
        k=args.k,
        n_pts_per_ball=128,
        seeds=[11, 29, 47],
        n_steps=160,
    )
    best_spec = GateSpec(**refined[0]["spec"])
    best_dense = dense_evaluation(best_spec, k=args.k, n_pts_per_ball=192, seeds=[11, 29, 47, 59, 83], n_steps=224)

    no_rotation_spec = GateSpec(**coarse["best_no_rotation_candidate"]["spec"])
    no_rotation_dense = dense_evaluation(
        no_rotation_spec,
        k=args.k,
        n_pts_per_ball=192,
        seeds=[11, 29, 47, 59, 83],
        n_steps=224,
    )

    best_groups = build_level_groups(args.k, best_spec)
    args.emit.write_text(render_hamiltonian(best_groups, best_spec, best_dense["translation_shift"]))

    summary = {
        "attempt": "attempt_019",
        "family": "localized_quadratic_squeeze_rotation_tree",
        "k": args.k,
        "eps": args.eps,
        "target_ratio": 1.0 - args.eps,
        "coarse_search": coarse,
        "refined_candidates": refined,
        "best_dense_evaluation": best_dense,
        "best_no_rotation_dense_evaluation": no_rotation_dense,
    }
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
