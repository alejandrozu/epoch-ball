#!/usr/bin/env python3
"""Audit a generated Hamiltonian.py against the packing target."""

from __future__ import annotations

import argparse
import importlib.util
import json
import statistics
import time
from pathlib import Path

import torch


torch.set_default_dtype(torch.float64)


def load_hamiltonian(module_path: Path):
    spec = importlib.util.spec_from_file_location("audited_hamiltonian", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "Hamiltonian"):
        raise AttributeError(f"{module_path} does not define Hamiltonian(Q, P, t)")
    return module.Hamiltonian


def sample_boundary_points(k: int, n_pts: int, seed: int) -> torch.Tensor:
    torch.manual_seed(seed)
    pts = []
    for i in range(k):
        x = torch.randn(n_pts, 4)
        x = x / x.norm(dim=1, keepdim=True)
        x[:, 0] += 3.0 * i
        pts.append(x)
    return torch.cat(pts, dim=0)


def hamiltonian_rhs(hamiltonian, q: torch.Tensor, p: torch.Tensor, t: float):
    q_var = q.detach().requires_grad_(True)
    p_var = p.detach().requires_grad_(True)
    h = hamiltonian(q_var, p_var, t)
    if h.shape != (q.shape[0],):
        raise ValueError(f"Hamiltonian returned shape {tuple(h.shape)}, expected {(q.shape[0],)}")
    gq, gp = torch.autograd.grad(h.sum(), [q_var, p_var], allow_unused=True)
    if gq is None:
        gq = torch.zeros_like(q_var)
    if gp is None:
        gp = torch.zeros_like(p_var)
    return gp, -gq


def rk4_flow(hamiltonian, q: torch.Tensor, p: torch.Tensor, n_steps: int):
    dt = 1.0 / n_steps
    for step in range(n_steps):
        t0 = step * dt
        k1q, k1p = hamiltonian_rhs(hamiltonian, q, p, t0)
        k2q, k2p = hamiltonian_rhs(hamiltonian, q + 0.5 * dt * k1q, p + 0.5 * dt * k1p, t0 + 0.5 * dt)
        k3q, k3p = hamiltonian_rhs(hamiltonian, q + 0.5 * dt * k2q, p + 0.5 * dt * k2p, t0 + 0.5 * dt)
        k4q, k4p = hamiltonian_rhs(hamiltonian, q + dt * k3q, p + dt * k3p, t0 + dt)
        q = q + (dt / 6.0) * (k1q + 2.0 * k2q + 2.0 * k3q + k4q)
        p = p + (dt / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p)
    return q, p


def audit_seed(hamiltonian, k: int, seed: int, n_pts: int, n_steps: int):
    x0 = sample_boundary_points(k, n_pts, seed)
    q0 = x0[:, :2].clone()
    p0 = x0[:, 2:].clone()
    q1, p1 = rk4_flow(hamiltonian, q0, p0, n_steps=n_steps)
    rsq = (q1.square()).sum(dim=1) + (p1.square()).sum(dim=1)
    radius = float(rsq.max().sqrt())
    ratio = k / (radius**4)
    return {
        "seed": seed,
        "n_pts": n_pts,
        "n_steps": n_steps,
        "radius": radius,
        "ratio": ratio,
    }


def benchmark_call(hamiltonian, batch_size: int):
    q = torch.randn(batch_size, 2)
    p = torch.randn(batch_size, 2)
    t0 = time.perf_counter()
    h = hamiltonian(q, p, 0.37)
    dt = time.perf_counter() - t0
    return {
        "batch_size": batch_size,
        "seconds": dt,
        "mean": float(h.mean()),
        "std": float(h.std()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", type=Path, required=True, help="Path to Hamiltonian.py")
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--eps", type=float, required=True)
    parser.add_argument("--n-pts", type=int, default=500)
    parser.add_argument("--n-steps", type=int, default=500)
    parser.add_argument("--seeds", type=int, nargs="+", default=[123])
    parser.add_argument("--batch-size", type=int, default=5000)
    args = parser.parse_args()

    hamiltonian = load_hamiltonian(args.module.resolve())
    bench = benchmark_call(hamiltonian, args.batch_size)
    results = [audit_seed(hamiltonian, args.k, seed, args.n_pts, args.n_steps) for seed in args.seeds]

    ratios = [entry["ratio"] for entry in results]
    summary = {
        "module": str(args.module.resolve()),
        "k": args.k,
        "eps": args.eps,
        "target_ratio": 1.0 - args.eps,
        "benchmark": bench,
        "runs": results,
        "min_ratio": min(ratios),
        "max_ratio": max(ratios),
        "mean_ratio": statistics.fmean(ratios),
        "passes_all_runs": min(ratios) >= 1.0 - args.eps,
        "implied_max_eps_from_min_ratio": 1.0 - min(ratios),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
