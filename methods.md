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

## attempt_004

- Family: compactly supported local affine symplectic pieces with rectangle-packed laminate blocks
- Mathematical modes: linear algebra, combinatorial search, inequalities/bounding
- Core idea: replace the nonlinear toric chart by explicit affine symplectic scalings and translations on each ball, group balls into affine blocks, and pack those blocks by q-plane rectangle placements
- Outcome: failed as a viable solution path; even with exact-radius evaluation the best `k=10` ratio found was only about `0.0364`, and the square case `k=16` still only reached about `0.0232`, showing the affine family is far too rigid

## attempt_005

- Family: smooth nonlinear slit-chart packing via exterior log-polar canonical coordinates and rectangle slotting
- Mathematical modes: calculus, discrete <-> continuous mapping, geometric
- Core idea: replace the singular toric polar chart by a smooth planar symplectic chart from an exterior point, flatten each unit disk to a rectangle by an exact area-preserving triangular map, and then reuse the standard `l^2` square-slot packing in those smooth coordinates
- Outcome: failed as a viable packing family; the planar chart itself is numerically correct with area `≈ pi` and tiny roundtrip error near the best parameters, but the best sampled ratios were only about `0.0769` for `k=10` and `0.1230` for `k=16`, showing that generic area coordinates destroy the toric packing mechanism

## Next family should avoid

- Another global polynomial or neural-style optimization ansatz
- Another purely obstruction-based analysis pass with no explicit constructive component
- Another purely square-grid or one-block toric pass that leaves the polar singularity unresolved
- Another purely affine local-piece or rectangular lamination ansatz
- Another generic slit-chart or exterior-log-polar rectangle-packing ansatz that replaces the toric action variables by arbitrary area coordinates and then applies the same square-slot recipe

## Most promising next family

- A smooth constructive family that preserves a simplex-type action relation while regularizing the singular toric chart, or a fundamentally different non-affine folding construction that retains the square-root compression mechanism
