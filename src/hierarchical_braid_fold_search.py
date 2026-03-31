#!/usr/bin/env python3
"""Attempt-018 hierarchical coupled braid / control-translation search.

This family is designed to be qualitatively different from both:

1. attempt 016's half-space transvection networks, and
2. attempt 017's cotangent-lift folds `H(q, p, t) = <p, X_t(q)>`.

The new primitive is a genuinely coupled four-dimensional Hamiltonian:

    H_A = alpha(q1, p1) * (dq2 * p2 - dp2 * q2)
    H_B = beta(q1, q2, p2) * collapse * p1

Stage A translates the second symplectic plane by an amount controlled by the
first plane. Stage B uses the resulting second-plane branch signal to translate
the first plane back toward the pair midpoint. The supports are organized in a
binary pairing tree, so the construction behaves like a recursive braid / fold.
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
class BraidSpec:
    shift0: float
    shift_growth: float
    collapse_frac: float
    return_frac: float
    angle_step: float
    support_scale: float
    sign_scale: float
    target_sign_scale: float
    p1_scale: float
    angle0: float = 0.0


@dataclass(frozen=True)
class PairNode:
    level: int
    q_mid: float
    support_q: float
    sign_width: float
    target_sign_width: float
    shift_q2: float
    shift_p2: float
    collapse_q1: float


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


def build_level_groups(k: int, spec: BraidSpec) -> list[list[PairNode]]:
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
            gap = right["center_x"] - left["center_x"]
            half_gap = 0.5 * gap
            angle = spec.angle0 + spec.angle_step * level
            shift_amp = spec.shift0 * (spec.shift_growth**level)
            node = PairNode(
                level=level,
                q_mid=0.5 * (left["center_x"] + right["center_x"]),
                support_q=spec.support_scale * half_gap + 1.25,
                sign_width=spec.sign_scale * half_gap + 0.35,
                target_sign_width=spec.target_sign_scale * shift_amp + 0.2,
                shift_q2=shift_amp * math.cos(angle),
                shift_p2=shift_amp * math.sin(angle),
                collapse_q1=spec.collapse_frac * half_gap,
            )
            group.append(node)
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


def rhs_stage_a(
    q1: np.ndarray,
    q2: np.ndarray,
    p1: np.ndarray,
    p2: np.ndarray,
    nodes: list[PairNode],
    weight: float,
    p1_scale: float,
    shift_scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    q1_dot = np.zeros_like(q1)
    q2_dot = np.zeros_like(q2)
    p1_dot = np.zeros_like(p1)
    p2_dot = np.zeros_like(p2)

    for node in nodes:
        shift_q2 = shift_scale * node.shift_q2
        shift_p2 = shift_scale * node.shift_p2
        z = q1 - node.q_mid
        a = z / node.support_q
        b = p1 / p1_scale
        window = np.exp(-(a**4) - (b**4))
        side_arg = z / node.sign_width
        side = np.tanh(side_arg)
        coeff = window * side
        lform = shift_q2 * p2 - shift_p2 * q2

        dwindow_dq1 = window * (-4.0 * (a**3) / node.support_q)
        dwindow_dp1 = window * (-4.0 * (b**3) / p1_scale)
        sech2 = 1.0 - side**2
        dside_dq1 = sech2 / node.sign_width
        dcoeff_dq1 = dwindow_dq1 * side + window * dside_dq1
        dcoeff_dp1 = dwindow_dp1 * side

        q1_dot += weight * dcoeff_dp1 * lform
        p1_dot += -weight * dcoeff_dq1 * lform
        q2_dot += weight * coeff * shift_q2
        p2_dot += weight * coeff * shift_p2

    return q1_dot, q2_dot, p1_dot, p2_dot


def rhs_stage_b(
    q1: np.ndarray,
    q2: np.ndarray,
    p1: np.ndarray,
    p2: np.ndarray,
    nodes: list[PairNode],
    weight: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    q1_dot = np.zeros_like(q1)
    q2_dot = np.zeros_like(q2)
    p1_dot = np.zeros_like(p1)
    p2_dot = np.zeros_like(p2)

    for node in nodes:
        z = q1 - node.q_mid
        a = z / node.support_q
        window = np.exp(-(a**4))

        shift_norm = math.hypot(node.shift_q2, node.shift_p2)
        uq = node.shift_q2 / shift_norm
        up = node.shift_p2 / shift_norm
        proj = uq * q2 + up * p2
        side_arg = proj / node.target_sign_width
        side = np.tanh(side_arg)
        coeff = window * side

        dwindow_dq1 = window * (-4.0 * (a**3) / node.support_q)
        sech2 = 1.0 - side**2
        dside_dproj = sech2 / node.target_sign_width

        dcoeff_dq1 = dwindow_dq1 * side
        dcoeff_dq2 = window * dside_dproj * uq
        dcoeff_dp2 = window * dside_dproj * up

        q1_dot += -weight * node.collapse_q1 * coeff
        p1_dot += weight * node.collapse_q1 * dcoeff_dq1 * p1
        q2_dot += -weight * node.collapse_q1 * dcoeff_dp2 * p1
        p2_dot += weight * node.collapse_q1 * dcoeff_dq2 * p1

    return q1_dot, q2_dot, p1_dot, p2_dot


def rhs(
    q: np.ndarray,
    p: np.ndarray,
    t: float,
    level_groups: list[list[PairNode]],
    spec: BraidSpec,
) -> tuple[np.ndarray, np.ndarray]:
    n_main_stages = 3 * len(level_groups)
    stage, weight = stage_weight(t, n_main_stages)
    if weight == 0.0:
        return np.zeros_like(q), np.zeros_like(p)

    q1 = q[:, 0]
    q2 = q[:, 1]
    p1 = p[:, 0]
    p2 = p[:, 1]

    level = stage // 3
    nodes = level_groups[level]
    stage_type = stage % 3
    if stage_type == 0:
        q1_dot, q2_dot, p1_dot, p2_dot = rhs_stage_a(q1, q2, p1, p2, nodes, weight, spec.p1_scale, shift_scale=1.0)
    elif stage_type == 1:
        q1_dot, q2_dot, p1_dot, p2_dot = rhs_stage_b(q1, q2, p1, p2, nodes, weight)
    else:
        q1_dot, q2_dot, p1_dot, p2_dot = rhs_stage_a(
            q1,
            q2,
            p1,
            p2,
            nodes,
            weight,
            spec.p1_scale,
            shift_scale=-spec.return_frac,
        )

    q_dot = np.stack([q1_dot, q2_dot], axis=1)
    p_dot = np.stack([p1_dot, p2_dot], axis=1)
    return q_dot, p_dot


def rk4_flow(q: np.ndarray, p: np.ndarray, level_groups: list[list[PairNode]], spec: BraidSpec, n_steps: int) -> tuple[np.ndarray, np.ndarray]:
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


def evaluate_points(points: np.ndarray, level_groups: list[list[PairNode]], spec: BraidSpec, n_steps: int) -> dict:
    q0 = points[:, :2]
    p0 = points[:, 2:]
    q1, p1 = rk4_flow(q0, p0, level_groups, spec, n_steps=n_steps)
    mapped = np.concatenate([q1, p1], axis=1)
    center, radius = ritter_enclosing_ball(mapped)
    return {
        "center": center.tolist(),
        "radius": radius,
        "ratio": 0.0 if radius <= 0.0 else None,
        "mapped_q_mean": q1.mean(axis=0).tolist(),
        "mapped_p_mean": p1.mean(axis=0).tolist(),
    }


def build_specs() -> list[BraidSpec]:
    specs: list[BraidSpec] = []
    for shift0 in [0.6, 0.9, 1.2]:
        for shift_growth in [1.0, 1.2, 1.4]:
            for collapse_frac in [0.8, 0.95]:
                for angle_step in [0.6, 0.8, 1.0]:
                    for support_scale in [1.1, 1.3]:
                        for target_sign_scale in [0.15, 0.3, 0.55]:
                            for p1_scale in [3.0, 6.0]:
                                specs.append(
                                    BraidSpec(
                                        shift0=shift0,
                                        shift_growth=shift_growth,
                                        collapse_frac=collapse_frac,
                                        return_frac=0.0,
                                        angle_step=angle_step,
                                        support_scale=support_scale,
                                        sign_scale=0.35,
                                        target_sign_scale=target_sign_scale,
                                        p1_scale=p1_scale,
                                    )
                                )
    return specs


def coarse_search(k: int, n_pts_per_ball: int, seed: int, n_steps: int, keep_top_n: int) -> dict:
    points = sample_boundary_points(k=k, n_pts_per_ball=n_pts_per_ball, seed=seed)
    rows = []
    for spec in build_specs():
        level_groups = build_level_groups(k, spec)
        eval_row = evaluate_points(points=points, level_groups=level_groups, spec=spec, n_steps=n_steps)
        radius = eval_row["radius"]
        ratio = 0.0 if radius <= 0.0 else k / (radius**4)
        rows.append(
            {
                "spec": asdict(spec),
                "node_counts_by_level": [len(group) for group in level_groups],
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
        spec = BraidSpec(**row["spec"])
        level_groups = build_level_groups(k, spec)
        seed_rows = []
        radii = []
        for seed in seeds:
            points = sample_boundary_points(k=k, n_pts_per_ball=n_pts_per_ball, seed=seed)
            eval_row = evaluate_points(points=points, level_groups=level_groups, spec=spec, n_steps=n_steps)
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


def level_group_payload(level_groups: list[list[PairNode]]) -> list[list[dict]]:
    payload = []
    for group in level_groups:
        payload.append([asdict(node) for node in group])
    return payload


def render_hamiltonian(
    level_groups: list[list[PairNode]],
    translation_shift: list[float],
    p1_scale: float,
    return_frac: float,
) -> str:
    groups = level_group_payload(level_groups)
    return f"""import math
import torch

torch.set_default_dtype(torch.float64)

_WINDOW_NORM = {WINDOW_NORM!r}
_MAIN_LEVELS = {groups!r}
_P1_SCALE = {p1_scale!r}
_RETURN_FRAC = {return_frac!r}
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


def _stage_a(Q, P, nodes, shift_scale):
    q1 = Q[:, 0]
    q2 = Q[:, 1]
    p1 = P[:, 0]
    p2 = P[:, 1]
    H = 0.0 * (q1 + q2 + p1 + p2)
    for node in nodes:
        shift_q2 = shift_scale * node["shift_q2"]
        shift_p2 = shift_scale * node["shift_p2"]
        z = q1 - node["q_mid"]
        a = z / node["support_q"]
        b = p1 / _P1_SCALE
        window = torch.exp(-(a.pow(4)) - (b.pow(4)))
        side = torch.tanh(z / node["sign_width"])
        coeff = window * side
        H = H + coeff * (shift_q2 * p2 - shift_p2 * q2)
    return H


def _stage_b(Q, P, nodes):
    q1 = Q[:, 0]
    q2 = Q[:, 1]
    p1 = P[:, 0]
    p2 = P[:, 1]
    H = 0.0 * (q1 + q2 + p1 + p2)
    for node in nodes:
        z = q1 - node["q_mid"]
        a = z / node["support_q"]
        window = torch.exp(-(a.pow(4)))
        shift_norm = math.hypot(node["shift_q2"], node["shift_p2"])
        uq = node["shift_q2"] / shift_norm
        up = node["shift_p2"] / shift_norm
        proj = uq * q2 + up * p2
        side = torch.tanh(proj / node["target_sign_width"])
        H = H - node["collapse_q1"] * window * side * p1
    return H


def Hamiltonian(Q, P, t):
    n_main_stages = 3 * len(_MAIN_LEVELS)
    total_stages = n_main_stages + 1
    stage, weight = _stage_weight(t, total_stages)
    H = 0.0 * (Q[:, 0] + P[:, 0])
    if weight == 0.0:
        return H
    if stage < n_main_stages:
        level = stage // 3
        stage_type = stage % 3
        if stage_type == 0:
            return weight * _stage_a(Q, P, _MAIN_LEVELS[level], 1.0)
        if stage_type == 1:
            return weight * _stage_b(Q, P, _MAIN_LEVELS[level])
        return weight * _stage_a(Q, P, _MAIN_LEVELS[level], -_RETURN_FRAC)
    dq = _TRANSLATE[:2].to(Q.device)
    dp = _TRANSLATE[2:].to(P.device)
    return weight * ((P * dq).sum(dim=1) - (Q * dp).sum(dim=1))
"""


def dense_evaluation_for_best(spec: BraidSpec, k: int, n_pts_per_ball: int, seeds: list[int], n_steps: int) -> dict:
    level_groups = build_level_groups(k, spec)
    centers = []
    rows = []
    radii = []
    for seed in seeds:
        points = sample_boundary_points(k=k, n_pts_per_ball=n_pts_per_ball, seed=seed)
        eval_row = evaluate_points(points=points, level_groups=level_groups, spec=spec, n_steps=n_steps)
        radius = eval_row["radius"]
        ratio = k / (radius**4)
        radii.append(radius)
        center = np.array(eval_row["center"], dtype=np.float64)
        centers.append(center)
        rows.append(
            {
                "seed": seed,
                "radius": radius,
                "ratio": ratio,
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
        "worst_radius": max(radii),
        "mean_radius": float(np.mean(radii)),
        "best_ratio": k / (min(radii) ** 4),
        "worst_ratio": k / (max(radii) ** 4),
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
        default=Path("results/research_loop/attempt_018_hierarchical_braid_search.json"),
    )
    parser.add_argument(
        "--emit",
        type=Path,
        default=Path("Hamiltonian.py"),
    )
    args = parser.parse_args()

    coarse = coarse_search(k=args.k, n_pts_per_ball=80, seed=11, n_steps=128, keep_top_n=8)
    refined = refine_candidates(
        coarse["top_candidates"],
        k=args.k,
        n_pts_per_ball=160,
        seeds=[11, 29, 47],
        n_steps=192,
    )
    best_spec = BraidSpec(**refined[0]["spec"])
    dense = dense_evaluation_for_best(
        best_spec,
        k=args.k,
        n_pts_per_ball=256,
        seeds=[11, 29, 47, 59, 83],
        n_steps=256,
    )

    level_groups = build_level_groups(args.k, best_spec)
    args.emit.write_text(
        render_hamiltonian(level_groups, dense["translation_shift"], best_spec.p1_scale, best_spec.return_frac)
    )

    summary = {
        "attempt": "attempt_018",
        "family": "hierarchical_coupled_braid_control_translation",
        "k": args.k,
        "eps": args.eps,
        "target_ratio": 1.0 - args.eps,
        "coarse_search": coarse,
        "refined_candidates": refined,
        "best_dense_evaluation": dense,
    }
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
