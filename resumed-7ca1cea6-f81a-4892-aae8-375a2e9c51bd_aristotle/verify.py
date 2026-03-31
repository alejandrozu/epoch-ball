#!/usr/bin/env python3
"""Verify the generated Hamiltonian by integrating the flow."""
import torch, math, sys
torch.set_default_dtype(torch.float64)

from Hamiltonian import Hamiltonian

def verify(k, eps, n_pts=500, n_steps=500):
    R_tgt = (k/(1-eps))**0.25
    print(f"Verifying k={k}, eps={eps}, R_target={R_tgt:.4f}")

    # Sample boundary of each ball
    torch.manual_seed(123)
    pts = []
    for i in range(k):
        x = torch.randn(n_pts, 4)
        x = x / x.norm(dim=1, keepdim=True)
        x[:,0] += 3.0*i
        pts.append(x)
    X = torch.cat(pts)
    Q = X[:,:2].clone(); P = X[:,2:].clone()

    # RK4 integration
    dt = 1.0 / n_steps
    for step in range(n_steps):
        t = step * dt

        def dqp(Q, P, t_val):
            Q_ = Q.detach().requires_grad_(True)
            P_ = P.detach().requires_grad_(True)
            H = Hamiltonian(Q_, P_, t_val)
            H_sum = H.sum()
            gQ, gP = torch.autograd.grad(H_sum, [Q_, P_])
            return gP, -gQ  # dQ/dt = dH/dP, dP/dt = -dH/dQ

        k1Q, k1P = dqp(Q, P, t)
        k2Q, k2P = dqp(Q+0.5*dt*k1Q, P+0.5*dt*k1P, t+0.5*dt)
        k3Q, k3P = dqp(Q+0.5*dt*k2Q, P+0.5*dt*k2P, t+0.5*dt)
        k4Q, k4P = dqp(Q+dt*k3Q, P+dt*k3P, t+dt)

        Q = Q + dt/6*(k1Q+2*k2Q+2*k3Q+k4Q)
        P = P + dt/6*(k1P+2*k2P+2*k3P+k4P)

    rsq = (Q**2).sum(1) + (P**2).sum(1)
    R = rsq.max().sqrt().item()
    ratio = k / R**4
    print(f"R = {R:.6f}, R^4 = {R**4:.4f}, ratio = {ratio:.6f}")
    print(f"Target ratio = {1-eps:.6f}")
    print(f"{'PASS' if ratio >= 1-eps else 'FAIL'}")
    return ratio >= 1-eps

if __name__ == '__main__':
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    eps = float(sys.argv[2]) if len(sys.argv) > 2 else 0.99
    verify(k, eps)
