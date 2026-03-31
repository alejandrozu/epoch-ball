#!/usr/bin/env python3
"""Witness-guided smooth kick-drift Hamiltonian fitting.

This attempt uses a new family distinct from the repository's original global
polynomial optimization. The flow is constrained to a composition of exact
Hamiltonian kicks and drifts:

  - kick segment: H(t, q, p) = sigma(t) * V(q)
  - drift segment: H(t, q, p) = sigma(t) * T(p)

Each kick or drift is integrated exactly over its time window because the
window integrates to 1. We train the segment parameters against the explicit
square-slot witness from attempt_003 on sampled boundary points.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import torch

from explicit_square_packing import square_packing_map


torch.set_default_dtype(torch.float64)


SQRT_PI_OVER_2 = math.sqrt(math.pi / 2.0)
SQRT_2 = math.sqrt(2.0)


@dataclass
class FitSummary:
    k: int
    n_blocks: int
    n_pts_per_ball: int
    n_iters: int
    seed: int
    train_radius: float
    train_ratio: float
    witness_radius: float
    witness_ratio: float
    mse_to_witness: float
    sigma_q1: float
    scale_log: float


def sample_boundary_points(k: int, n_pts_per_ball: int, seed: int) -> torch.Tensor:
    torch.manual_seed(seed)
    pts = []
    for i in range(k):
        x = torch.randn(n_pts_per_ball, 4)
        x = x / x.norm(dim=1, keepdim=True)
        x[:, 0] += 3.0 * i
        pts.append(x)
    return torch.cat(pts, dim=0)


def centered_centers(k: int, device: torch.device | None = None) -> torch.Tensor:
    ctr = 1.5 * (k - 1)
    values = [3.0 * i - ctr for i in range(k)]
    return torch.tensor(values, dtype=torch.get_default_dtype(), device=device)


def target_square_witness(points: torch.Tensor, k: int) -> torch.Tensor:
    mapped, _ = square_packing_map(points, k)
    return mapped


class KickDriftModel(torch.nn.Module):
    def __init__(self, k: int, n_blocks: int, sigma_q1: float = 0.7):
        super().__init__()
        self.k = k
        self.n_blocks = n_blocks
        self.sigma_q1 = float(sigma_q1)
        self.center_shift = -1.5 * (k - 1)
        self.scale_log = torch.nn.Parameter(torch.tensor(-1.15, dtype=torch.get_default_dtype()))

        # Kick coefficients per block and per source ball.
        self.kick_a = torch.nn.Parameter(torch.zeros(n_blocks, k))
        self.kick_b = torch.nn.Parameter(torch.zeros(n_blocks, k))
        self.kick_c = torch.nn.Parameter(torch.zeros(n_blocks, k))
        self.kick_d = torch.nn.Parameter(torch.zeros(n_blocks, k))

        # Global drift coefficients per block.
        self.drift_lin = torch.nn.Parameter(torch.zeros(n_blocks, 2))
        self.drift_quad = torch.nn.Parameter(torch.zeros(n_blocks, 3))
        self.drift_cubic = torch.nn.Parameter(torch.zeros(n_blocks, 4))

        self._initialize()

    def _initialize(self) -> None:
        with torch.no_grad():
            # Use a weak row-like initialization in q2 so the model does not
            # start from the degenerate translate-only seed.
            l = math.ceil(math.sqrt(self.k))
            for i in range(self.k):
                row = i // l
                col = i % l
                self.kick_b[0, i] = 0.35 * (row - 0.5 * (l - 1))
                self.kick_a[0, i] = -0.15 * (col - 0.5 * (l - 1))

            # First drift recenters the union in q1.
            self.drift_lin[0, 0] = self.center_shift

    def _kick_grad(self, q: torch.Tensor, block: int) -> tuple[torch.Tensor, torch.Tensor]:
        centers = centered_centers(self.k, device=q.device)
        q1 = q[:, 0:1]
        q2 = q[:, 1:2]
        s = (q1 - centers.view(1, -1)) / self.sigma_q1
        bump = torch.exp(-0.5 * s.square())

        a = self.kick_a[block].view(1, -1)
        b = self.kick_b[block].view(1, -1)
        c = self.kick_c[block].view(1, -1)
        d = self.kick_d[block].view(1, -1)

        poly = b * q2 + 0.5 * c * q2.square() + (d / 3.0) * q2.pow(3)
        dVdq1 = (bump * a).sum(dim=1) + ((-s / self.sigma_q1) * bump * poly).sum(dim=1)
        dVdq2 = (bump * (b + c * q2 + d * q2.square())).sum(dim=1)
        return dVdq1, dVdq2

    def _drift_grad(self, p: torch.Tensor, block: int) -> tuple[torch.Tensor, torch.Tensor]:
        p1 = p[:, 0]
        p2 = p[:, 1]

        lin1, lin2 = self.drift_lin[block]
        q11, q22, q12 = self.drift_quad[block]
        c1, c2, c12, c21 = self.drift_cubic[block]

        dTdp1 = lin1 + q11 * p1 + q12 * p2 + c1 * p1.square() + 2.0 * c12 * p1 * p2 + c21 * p2.square()
        dTdp2 = lin2 + q22 * p2 + q12 * p1 + c2 * p2.square() + c12 * p1.square() + 2.0 * c21 * p1 * p2
        return dTdp1, dTdp2

    def discrete_map(self, points: torch.Tensor) -> torch.Tensor:
        q = points[:, :2].clone()
        p = points[:, 2:].clone()

        q = torch.stack((q[:, 0] + self.center_shift, q[:, 1]), dim=1)
        scale = torch.exp(self.scale_log)
        q = scale * q
        p = p / scale

        for block in range(self.n_blocks):
            dVdq1, dVdq2 = self._kick_grad(q, block)
            p = torch.stack((p[:, 0] - dVdq1, p[:, 1] - dVdq2), dim=1)

            dTdp1, dTdp2 = self._drift_grad(p, block)
            q = torch.stack((q[:, 0] + dTdp1, q[:, 1] + dTdp2), dim=1)

        return torch.cat((q, p), dim=1)

    def segment_dict(self) -> dict:
        return {
            "k": self.k,
            "n_blocks": self.n_blocks,
            "sigma_q1": self.sigma_q1,
            "center_shift": self.center_shift,
            "scale_log": float(self.scale_log.detach().cpu()),
            "kick_a": self.kick_a.detach().cpu().tolist(),
            "kick_b": self.kick_b.detach().cpu().tolist(),
            "kick_c": self.kick_c.detach().cpu().tolist(),
            "kick_d": self.kick_d.detach().cpu().tolist(),
            "drift_lin": self.drift_lin.detach().cpu().tolist(),
            "drift_quad": self.drift_quad.detach().cpu().tolist(),
            "drift_cubic": self.drift_cubic.detach().cpu().tolist(),
        }


def surrogate_radius_loss(mapped: torch.Tensor, beta: float) -> torch.Tensor:
    rsq = mapped.square().sum(dim=1)
    return torch.logsumexp(beta * rsq, dim=0) / beta


def fit_model(
    k: int,
    n_blocks: int,
    n_pts_per_ball: int,
    seed: int,
    n_iters: int,
    lr: float,
    sigma_q1: float,
    mse_weight: float,
) -> tuple[KickDriftModel, FitSummary, dict]:
    points = sample_boundary_points(k, n_pts_per_ball=n_pts_per_ball, seed=seed)
    witness = target_square_witness(points, k)
    witness_radius = float(witness.norm(dim=1).max())
    witness_ratio = k / (witness_radius**4)

    model = KickDriftModel(k=k, n_blocks=n_blocks, sigma_q1=sigma_q1)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    best = {
        "radius": float("inf"),
        "ratio": 0.0,
        "mse": float("inf"),
        "state": None,
    }

    for it in range(n_iters):
        beta = min(0.2 + 0.003 * it, 6.0)
        mapped = model.discrete_map(points)
        rsq_loss = surrogate_radius_loss(mapped, beta=beta)
        mse = (mapped - witness).square().mean()
        reg = 1e-5 * sum(param.square().mean() for param in model.parameters())
        loss = rsq_loss + mse_weight * mse + reg

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(list(model.parameters()), max_norm=10.0)
        optimizer.step()

        with torch.no_grad():
            radius = float(mapped.norm(dim=1).max())
            ratio = k / (radius**4)
            mse_value = float(mse)
            if radius < best["radius"]:
                best["radius"] = radius
                best["ratio"] = ratio
                best["mse"] = mse_value
                best["state"] = {name: tensor.detach().clone() for name, tensor in model.state_dict().items()}

        if it % 200 == 0 or it == n_iters - 1:
            print(
                f"iter={it:4d} radius={radius:.6f} ratio={ratio:.6f} "
                f"mse={mse_value:.6f} beta={beta:.3f}",
                flush=True,
            )

    if best["state"] is not None:
        model.load_state_dict(best["state"])

    with torch.no_grad():
        mapped = model.discrete_map(points)
        train_radius = float(mapped.norm(dim=1).max())
        train_ratio = k / (train_radius**4)
        mse_to_witness = float((mapped - witness).square().mean())

    summary = FitSummary(
        k=k,
        n_blocks=n_blocks,
        n_pts_per_ball=n_pts_per_ball,
        n_iters=n_iters,
        seed=seed,
        train_radius=train_radius,
        train_ratio=train_ratio,
        witness_radius=witness_radius,
        witness_ratio=witness_ratio,
        mse_to_witness=mse_to_witness,
        sigma_q1=sigma_q1,
        scale_log=float(model.scale_log.detach().cpu()),
    )
    return model, summary, model.segment_dict()


def c_infty_normalizer(n_grid: int = 200001) -> float:
    x = torch.linspace(0.0, 1.0, n_grid)
    y = torch.zeros_like(x)
    mask = (x > 0.0) & (x < 1.0)
    xm = x[mask]
    y[mask] = torch.exp(-1.0 / (xm * (1.0 - xm)))
    dx = 1.0 / (n_grid - 1)
    return float(torch.trapezoid(y, dx=dx))


def emit_hamiltonian(module_path: Path, payload: dict) -> None:
    normalizer = c_infty_normalizer()
    code = f"""import math
import torch

torch.set_default_dtype(torch.float64)

_K = {payload['k']}
_N_BLOCKS = {payload['n_blocks']}
_SIGMA_Q1 = {payload['sigma_q1']!r}
_CENTER_SHIFT = {payload['center_shift']!r}
_SCALE_LOG = {payload['scale_log']!r}
_KICK_A = {json.dumps(payload['kick_a'])}
_KICK_B = {json.dumps(payload['kick_b'])}
_KICK_C = {json.dumps(payload['kick_c'])}
_KICK_D = {json.dumps(payload['kick_d'])}
_DRIFT_LIN = {json.dumps(payload['drift_lin'])}
_DRIFT_QUAD = {json.dumps(payload['drift_quad'])}
_DRIFT_CUBIC = {json.dumps(payload['drift_cubic'])}
_WINDOW_NORM = {normalizer!r}
_SQRT_PI_OVER_2 = {SQRT_PI_OVER_2!r}
_SQRT_2 = {SQRT_2!r}
_CENTERS = torch.tensor([3.0 * i + _CENTER_SHIFT for i in range(_K)], dtype=torch.float64)


def _bump_window(t, start, end):
    if not (start < t < end):
        return 0.0
    u = (t - start) / (end - start)
    if not (0.0 < u < 1.0):
        return 0.0
    return math.exp(-1.0 / (u * (1.0 - u))) / ((end - start) * _WINDOW_NORM)


def _kick_potential(Q, block):
    q1 = Q[:, 0:1]
    q2 = Q[:, 1:2]
    s = (q1 - _CENTERS.view(1, -1).to(Q.device)) / _SIGMA_Q1
    bump = torch.exp(-0.5 * s.square())
    primitive = _SIGMA_Q1 * _SQRT_PI_OVER_2 * torch.erf(s / _SQRT_2)
    a = torch.tensor(_KICK_A[block], dtype=Q.dtype, device=Q.device).view(1, -1)
    b = torch.tensor(_KICK_B[block], dtype=Q.dtype, device=Q.device).view(1, -1)
    c = torch.tensor(_KICK_C[block], dtype=Q.dtype, device=Q.device).view(1, -1)
    d = torch.tensor(_KICK_D[block], dtype=Q.dtype, device=Q.device).view(1, -1)
    poly = b * q2 + 0.5 * c * q2.square() + (d / 3.0) * q2.pow(3)
    return (primitive * a + bump * poly).sum(dim=1)


def _drift_potential(P, block):
    p1 = P[:, 0]
    p2 = P[:, 1]
    lin1, lin2 = _DRIFT_LIN[block]
    q11, q22, q12 = _DRIFT_QUAD[block]
    c1, c2, c12, c21 = _DRIFT_CUBIC[block]
    return (
        lin1 * p1
        + lin2 * p2
        + 0.5 * q11 * p1.square()
        + 0.5 * q22 * p2.square()
        + q12 * p1 * p2
        + (c1 / 3.0) * p1.pow(3)
        + (c2 / 3.0) * p2.pow(3)
        + c12 * p1.square() * p2
        + c21 * p1 * p2.square()
    )


def Hamiltonian(Q, P, t):
    pre_size = 1.0 / (2 * _N_BLOCKS + 2)
    block_size = pre_size
    H = 0.0 * (Q[:, 0] + P[:, 0])
    w_center = _bump_window(t, 0.0, pre_size)
    w_scale = _bump_window(t, pre_size, 2.0 * pre_size)
    if w_center != 0.0:
        H = H + w_center * (_CENTER_SHIFT * P[:, 0])
    if w_scale != 0.0:
        H = H + w_scale * _SCALE_LOG * (
            Q[:, 0] * P[:, 0] + Q[:, 1] * P[:, 1]
        )
    for block in range(_N_BLOCKS):
        kick_start = (2 + 2 * block) * block_size
        kick_end = kick_start + block_size
        drift_start = kick_end
        drift_end = drift_start + block_size
        wk = _bump_window(t, kick_start, kick_end)
        wd = _bump_window(t, drift_start, drift_end)
        if wk != 0.0:
            H = H + wk * _kick_potential(Q, block)
        if wd != 0.0:
            H = H + wd * _drift_potential(P, block)
    return H
"""
    module_path.write_text(code)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--n-blocks", type=int, default=4)
    parser.add_argument("--n-pts", type=int, default=160)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--n-iters", type=int, default=1200)
    parser.add_argument("--lr", type=float, default=3e-2)
    parser.add_argument("--sigma-q1", type=float, default=0.7)
    parser.add_argument("--mse-weight", type=float, default=0.05)
    parser.add_argument("--out", type=Path, default=Path("Hamiltonian.py"))
    parser.add_argument("--summary-out", type=Path, default=None)
    args = parser.parse_args()

    model, summary, payload = fit_model(
        k=args.k,
        n_blocks=args.n_blocks,
        n_pts_per_ball=args.n_pts,
        seed=args.seed,
        n_iters=args.n_iters,
        lr=args.lr,
        sigma_q1=args.sigma_q1,
        mse_weight=args.mse_weight,
    )

    emit_hamiltonian(args.out, payload)
    result = {
        "fit_summary": asdict(summary),
        "payload": payload,
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    print(text)
    if args.summary_out is not None:
        args.summary_out.write_text(text)


if __name__ == "__main__":
    main()
