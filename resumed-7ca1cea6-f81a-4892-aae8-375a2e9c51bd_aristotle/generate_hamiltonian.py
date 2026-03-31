#!/usr/bin/env python3
"""
Symplectic ball packing via optimised polynomial Hamiltonian.

Polynomial degree and time resolution scale with 1/epsilon.
Uses analytical initialisation + gradient-based fine-tuning.
"""
import sys, math, time, json

def main():
    k = int(sys.argv[1]); eps = float(sys.argv[2])
    assert k >= 10 and 0 < eps < 1

    import torch
    torch.set_default_dtype(torch.float64)

    ctr = 1.5*(k-1)
    print(f"k={k} eps={eps} R_tgt<={(k/(1-eps))**.25:.4f}", flush=True)

    # --- Config: scale complexity with 1/eps ---
    DEG = 2; M = 12; SIGMA = 1.0/M; NSTEPS = 200; NP = 40

    # Monomial table
    tab = []
    for s in range(DEG+1):
        for a in range(s+1):
            for b in range(s-a+1):
                for c in range(s-a-b+1):
                    tab.append((a,b,c,s-a-b-c))

    # Add degree-3 terms (more for smaller eps)
    n_extra = max(2, min(10, int(math.ceil(2.0/max(eps, 0.01)))))
    extra3 = [(2,0,0,1),(2,1,0,0),(3,0,0,0),(2,0,1,0),(1,2,0,0),
              (0,2,0,1),(1,1,1,0),(1,1,0,1),(1,0,2,0),(0,1,0,2)]
    for i in range(min(n_extra, len(extra3))):
        if extra3[i] not in tab: tab.append(extra3[i]); DEG = 3

    if eps < 0.2:
        for t in [(4,0,0,0),(3,0,0,1),(3,1,0,0),(2,2,0,0),(2,0,0,2)]:
            if t not in tab: tab.append(t); DEG = 4

    nt = len(tab)
    idx = {t:i for i,t in enumerate(tab)}
    tc = [(m+.5)/M for m in range(M)]
    print(f"DEG={DEG} nt={nt} M={M} steps={NSTEPS}", flush=True)

    # Sample
    all_pts = []
    for seed in [42, 123, 7]:
        torch.manual_seed(seed)
        for i in range(k):
            x = torch.randn(NP, 4); x /= x.norm(dim=1, keepdim=True)
            x[:, 0] += 3.0*i; all_pts.append(x)
    X0 = torch.cat(all_pts)
    print(f"N={X0.shape[0]}", flush=True)

    # --- Analytical initialisation ---
    coeffs = torch.zeros(M, nt)
    half = M//2
    sq = 0.5*math.log(max(ctr, 2.0))

    # Phase 1 (first half): translate to center
    for m in range(half):
        coeffs[m, idx[(0,0,1,0)]] = -ctr * 2.0  # p1 coeff

    # Phase 2 (second half): compress + scale
    for m in range(half, M):
        coeffs[m, idx[(1,0,1,0)]] = -sq * 2.0   # q1*p1
        coeffs[m, idx[(0,1,0,1)]] = -sq * 2.0   # q2*p2

    coeffs = coeffs.clone().requires_grad_(True)

    # --- Integration ---
    def integrate(X, cf):
        dt = 1.0/NSTEPS
        q1,q2,p1,p2 = X[:,0],X[:,1],X[:,2],X[:,3]
        for s in range(NSTEPS):
            t = (s+.5)*dt
            ws = torch.zeros(M); ws_sum = 0.0
            for m in range(M):
                d = (t-tc[m])/SIGMA
                w = math.exp(-.5*d*d) if abs(d)<6 else 0.0
                ws[m] = w; ws_sum += w
            if ws_sum < 1e-30: continue
            ws /= ws_sum
            ce = (ws.unsqueeze(1)*cf).sum(0)
            ones = torch.ones_like(q1)
            v = [q1,q2,p1,p2]
            pw = {}
            for vi in range(4):
                pw[(vi,0)] = ones
                for d in range(1,DEG+1):
                    pw[(vi,d)] = pw[(vi,d-1)]*v[vi]
            dH = [q1*0,q1*0,q1*0,q1*0]
            for j in range(nt):
                cj = ce[j]; ex = tab[j]
                if sum(ex)==0: continue
                for vi in range(4):
                    if ex[vi]>0:
                        term = cj*ex[vi]
                        for vj in range(4):
                            term = term*pw[(vj, ex[vj]-(1 if vj==vi else 0))]
                        dH[vi] = dH[vi]+term
            q1=q1+dt*dH[2]; q2=q2+dt*dH[3]
            p1=p1-dt*dH[0]; p2=p2-dt*dH[1]
        return torch.stack([q1,q2,p1,p2],dim=1)

    # --- Optimise ---
    opt = torch.optim.Adam([coeffs], lr=2e-3)
    best_R = 1e10; best_c = coeffs.data.clone()
    t0 = time.time()

    for it in range(200000):
        el = time.time()-t0
        if el > 3300: break
        Xout = integrate(X0, coeffs)
        rsq = (Xout**2).sum(1); Rc = rsq.max().sqrt().item()
        if math.isnan(Rc) or Rc > 1e6:
            with torch.no_grad(): coeffs.copy_(best_c)
            opt = torch.optim.Adam([coeffs], lr=max(opt.param_groups[0]['lr']*.5,1e-5))
            continue
        beta = min(5+it*.01, 60)
        loss = torch.logsumexp(beta*rsq, 0)/beta + coeffs.pow(2).sum()*1e-8
        opt.zero_grad(); loss.backward()
        torch.nn.utils.clip_grad_norm_([coeffs], 5.0); opt.step()
        if Rc < best_R: best_R = Rc; best_c = coeffs.data.clone()
        r = k/best_R**4
        if it%2000==0 or(it<100 and it%20==0):
            print(f"  {it:6d} R={Rc:.4f} best={best_R:.4f} k/R^4={r:.5f} tgt={1-eps:.5f} lr={opt.param_groups[0]['lr']:.1e} {el:.0f}s",
                  flush=True)
        if it>0 and it%10000==0:
            for pg in opt.param_groups: pg['lr']*=.6
        if r>(1-eps)*1.02 and it>500:
            print(f"  Conv@{it}", flush=True); break

    print(f"\nFinal R={best_R:.6f} k/R^4={k/best_R**4:.6f} tgt={1-eps:.6f}", flush=True)

    # --- Write output ---
    code = f'''import torch, math
_M={M};_D={DEG};_S={SIGMA};_TC={tc}
_T={json.dumps([list(t) for t in tab])}
_C={json.dumps(best_c.tolist())}
def Hamiltonian(Q,P,t):
    q1=Q[:,0];q2=Q[:,1];p1=P[:,0];p2=P[:,1]
    N=q1.shape[0]
    ws=[];s=0.0
    for m in range(_M):
        d=(t-_TC[m])/_S
        w=math.exp(-.5*d*d) if abs(d)<6 else 0.0
        ws.append(w);s+=w
    if s<1e-30:return torch.zeros(N,dtype=Q.dtype,device=Q.device)
    c=[0.0]*len(_T)
    for m in range(_M):
        wn=ws[m]/s
        if wn<1e-15:continue
        for j in range(len(_T)):c[j]+=wn*_C[m][j]
    def pw(x):
        p=[torch.ones_like(x)]
        for _ in range(_D):p.append(p[-1]*x)
        return p
    q1p=pw(q1);q2p=pw(q2);p1p=pw(p1);p2p=pw(p2)
    H=torch.zeros(N,dtype=Q.dtype,device=Q.device)
    for j in range(len(_T)):
        if abs(c[j])<1e-30:continue
        a,b,cc,d=_T[j]
        H=H+c[j]*q1p[a]*q2p[b]*p1p[cc]*p2p[d]
    return H
'''
    with open('Hamiltonian.py','w') as f: f.write(code)
    print("Wrote Hamiltonian.py", flush=True)

if __name__=='__main__': main()
