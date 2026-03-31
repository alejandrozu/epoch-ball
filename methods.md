# Tried Method Families

## attempt_001

- Family: global optimization of a time-dependent low-degree polynomial Hamiltonian
- Mathematical modes: optimization, dynamical systems
- Core idea: train polynomial coefficients against sampled boundary points so the time-1 Hamiltonian flow compresses the ball union into a smaller target ball
- Outcome: failed; the optimized Hamiltonian underperformed its own linear seed and only handled a near-trivial epsilon regime

## attempt_002

- Family: geometric/algebraic obstruction analysis via equal-ball packing theory and Cremona reduction
- Mathematical modes: geometric, algebraic, inequalities/bounding
- Core idea: replace ad hoc optimization with the standard packing benchmark for `B^4`, compute the abstract packing threshold, and compare the repository candidate to that benchmark
- Outcome: succeeded as analysis; reproduced the classical packing table and showed the task regime `k>=10` has full abstract packing, but did not yet yield an explicit Hamiltonian construction

## attempt_003

- Family: explicit constructive toric/action-angle packing with square-slot embeddings and multi-block moment-triangle search
- Mathematical modes: discrete <-> continuous mapping, geometric, algorithmic
- Core idea: realize actual symplectic embeddings on the source union using toric coordinates, first via the standard `l^2` square packing and then via combinatorial searches over toric blocks inside the target simplex
- Outcome: partial success; built a validated explicit baseline with exact ratio `k / ceil(sqrt(k))^2`, found a better `k=10` multi-block arrangement with ratio about `0.6925`, but still did not reach arbitrary `1-epsilon` density and did not yet produce the required smooth global Hamiltonian

## Next family should avoid

- Another global polynomial or neural-style optimization ansatz
- Another purely obstruction-based analysis pass with no explicit constructive component
- Another purely square-grid or one-block toric pass that leaves the polar singularity unresolved

## Most promising next family

- Smooth constructive packing from toric/simplex geometry using slit-disk or Traynor-style regularized charts, plus richer non-axis multi-block layouts or localized Hamiltonian pieces
