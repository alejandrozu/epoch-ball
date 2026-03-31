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

## attempt_006

- Family: witness-guided smooth kick-drift Hamiltonian factorization with source-indexed Gaussian kicks
- Mathematical modes: recursion / recurrence, dynamical systems
- Core idea: replace generic global polynomials by a composition of exact kick and drift Hamiltonian segments, fit the segment parameters against the explicit square-slot witness on sampled boundary points, and emit a real smooth `Hamiltonian.py`
- Outcome: failed as a viable Hamiltonian family; although the emitted Hamiltonian was smooth and fast enough to evaluate, the best discrete surrogate ratio was only about `3.45e-4` and the audited RK4 ratio dropped further to about `6.60e-5`, far below even the repository's original poor polynomial baseline

## attempt_007

- Family: topological origami moment-tiling by tiny-triangle reassignment in toric coordinates
- Mathematical modes: topological, discrete <-> continuous mapping, combinatorial
- Core idea: subdivide each source moment triangle into `m^2` tiny lattice triangles and map them piecewise-affinely into tiny target triangles inside `Delta(ceil(m sqrt(k)) / m)`, using local affine symplectic toric charts on each tiny triangle
- Outcome: partial but invalid as a Hamiltonian solution; the exact density formula was validated up to ratio `0.9999269539810081` for `k=10`, `eps=1e-4`, but every nontrivial `k=10` test had zero continuous internal source adjacencies and `O(1)` fold jumps, so the current assignment is not even continuous

## attempt_008

- Family: graph/symmetry reduction of continuity-preserving simplicial origami to rigid lattice-triangle packing
- Mathematical modes: graph / network, symmetry / group
- Core idea: impose continuity at the simplicial unit-triangle level, prove by row propagation that the source subdivision can only map to a rigid side-m lattice triangle, and then solve the reduced k-ball problem as an exact compatibility-clique search over rigid triangle placements inside the target triangular grid
- Outcome: strong negative evidence; the local rigidity reduction held exactly for all tested seeds and all `m <= 20`, and exact searches showed that for `k=10` the reduced family fits at most 9 rigid side-m triangles for every tested `m = 3, ..., 11`

## attempt_009

- Family: exact / probabilistic search over connected non-simplicial polyomino action-shapes with integer toric covers
- Mathematical modes: probabilistic, combinatorial, discrete <-> continuous mapping
- Core idea: replace the rigid simplicial image by a connected equal-area polyomino in moment space, search exact and sampled self-avoiding-walk polyominoes, and combine them with determinant-`>= k` integer cover matrices to minimize the target simplex capacity `mu`
- Outcome: useful near-match but negative overall; the exact unrestricted `n_cells=10` optimum stabilized at ratio `0.6920415224913494` for both `entry_bound=6` and `entry_bound=8`, still below attempt 003, while the stricter path-only strip-snake subfamily dropped to ratio `0.5540166204986149`, showing that the best surrogate witness depends on branching polyomino geometry and still lacks an explicit continuous/smooth realization

## Next family should avoid

- Another global polynomial or neural-style optimization ansatz
- Another purely obstruction-based analysis pass with no explicit constructive component
- Another purely square-grid or one-block toric pass that leaves the polar singularity unresolved
- Another purely affine local-piece or rectangular lamination ansatz
- Another generic slit-chart or exterior-log-polar rectangle-packing ansatz that replaces the toric action variables by arbitrary area coordinates and then applies the same square-slot recipe
- Another separable kick-drift boundary-fitting ansatz that tries to learn the branch structure of the toric witness through source-indexed Gaussian kicks and global drifts
- Another arbitrary tiny-triangle origami permutation that ignores global vertex and edge compatibility across fold lines
- Another continuity-preserving unit-triangle simplicial origami ansatz whose single-ball images are forced to rigid lattice triangles and then reduced to rigid triangle set packing
- Another single connected polyomino / strip-snake action-shape search with integer toric covers whose exact `n_cells=10` optimum already stabilized below attempt_003 and whose path-only subfamily collapsed to ratio `0.5540`

## Most promising next family

- A non-simplicial construction with explicit continuous gluing from the start, such as a smooth multi-chart shear / generating-function method or a multi-block non-simplicial complex that is not reduced to one connected action-shape
