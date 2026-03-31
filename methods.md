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

## Next family should avoid

- Another global polynomial or neural-style optimization ansatz
- Another purely obstruction-based analysis pass with no explicit constructive component

## Most promising next family

- Explicit constructive packing from toric/simplex geometry or from compositions of compactly supported local Hamiltonian moves
