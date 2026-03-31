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

## attempt_010

- Family: hybrid set/logic `9 rigid + 1 strip` decomposition on the triangular grid
- Mathematical modes: set / logic, algorithmic
- Core idea: for `k=10`, split the near-full target set into an inner side-`3m` triangle that tiles exactly into 9 rigid side-`m` triangles plus one remaining boundary strip/tree complex for the 10th ball, then test whether the `m=6` boundary path admits a continuous equal-area source-side strip realization
- Outcome: promising but still incomplete; at `m=6` the boundary region is exactly a 36-face path strip and the source-side boundary-realization search found 36 genuinely non-collinear corner placements plus a strictly positive best candidate with counts `[9,19,10]`, but the face areas still stalled at MSE about `4.65e-05` and max error about `0.0111`, while the naive boundary-prefix region also branches immediately for larger `m`

## attempt_011

- Family: algebraic / number-theoretic exact recurrence solver for the hybrid boundary strip
- Mathematical modes: algebraic, number theory, dynamical systems
- Core idea: exactify the saved `m=6` boundary-only strip by solving its equal-area boundary recurrence in rational arithmetic, then compute the induced local toric affine maps and audit whether the angle matrices glue across internal edges
- Outcome: strong exact negative evidence; the unique rational seed `17/36` certifies an exact 36-face equal-area moment realization, but every one of the 35 internal full-torus edges has mismatched `A^{-T}` matrices and the 36 faces all carry distinct affine action matrices, so the boundary-only strip is not a continuous 4D toric piecewise-affine symplectic embedding

## attempt_012

- Family: transform-method / generating-function cotangent-lift audit for piecewise-affine hybrid strips
- Mathematical modes: transform methods, calculus, topological
- Core idea: enlarge each affine face map to the most general fiber-preserving exact symplectomorphism over the same base map, derive the restricted-angle continuity rule on shared edges, and test whether facewise generating functions can repair the exact `m=6` hybrid strip
- Outcome: stronger negative evidence; the toy control shows that only coordinate-axis edges can hide derivative jumps because those are the only loci where a torus factor collapses in `B^4`, while the exact hybrid strip has 35 internal edges all with source and target active angle sets `[0,1]` and zero restricted `A^{-T}` matches, so no facewise generating-function correction can repair it

## attempt_013

- Family: axis-seam two-piece rigid origami search in the standard simplex
- Mathematical modes: game-theoretic, information-theoretic, set / logic
- Core idea: accept attempt 012's coordinate-axis restriction and test the simplest nontrivial family it leaves open: split one source side-`m` triangle into two rigid affine pieces along an internal straight cut, force the common seam to land on `x1=0` or `x2=0`, and enumerate every rigid placement of both pieces inside the near-full target simplex
- Outcome: stronger negative evidence on the axis-seam line; for `k=10` and every tested `m = 3,4,5,6,7,8`, the exact search found many matching axis-line candidates but zero interior-disjoint pairs, so even the simplest two-piece rigid coordinate-axis seam does not produce a valid origami fold inside the simplex

## attempt_014

- Family: local boundary-uniqueness obstruction for rigid axis-seam origami
- Mathematical modes: topological, inequalities / bounding, set / logic
- Core idea: replace the brute-force search by a local combinatorial theorem on the target simplex: every unit boundary edge has a unique incident unit triangle, so if an internal rigid seam edge lands on `x1=0` or `x2=0`, the two source unit triangles adjacent to that seam edge must map to the same target unit triangle and overlap; audit the saved attempt-013 artifact against the resulting lower bound `overlap >= seam length`
- Outcome: stronger general negative evidence; the local theorem rules out the entire rigid axis-seam strategy, not just the searched two-piece straight-cut cases, and all 81 cuts in the saved attempt-013 artifact satisfy the theoremic lower bound

## attempt_015

- Family: slice-preserving affine transvection-fold audit with symplectic left-factor mixing
- Mathematical modes: linear algebra, geometric, dynamical systems
- Core idea: replace rigid seam origami by a one-sided affine symplectic transvection fold across a hyperplane, normalize it to the active-slice map `F_{s,lambda,d}`, and test whether a common linear symplectic left factor that mixes the active and passive planes can evade the preserved-slice obstruction
- Outcome: stronger exact negative evidence on this new non-rigid family; the central active slice keeps symplectic area `pi`, so any containing ball radius obeys `R >= delta^(-1/2) >= 1`, where `delta = |omega(f1,f2)|` for an orthonormal basis of the mixed slice plane. The new script exact-checks the transvection formulas, runs a 5776-case nontrivial grid search with best sampled radius `1.001083448194763`, and audits 16 random symplectic mixers with theoretical lower bounds between `1.0083706647215258` and `2.1715166455758785` and minimum sampled radius `1.4699275946793848`. A naive attempt to extend the same area-`pi` argument to arbitrary linear pre-factors failed numerically, so pre-mixing remains a live escape route rather than something already ruled out

## attempt_016

- Family: multistage one-sided symplectic transvection network search
- Mathematical modes: algorithmic, dynamical systems, geometric
- Core idea: leave the slice-preserving/common-left-factor family by composing several exact one-sided transvection folds with stage-dependent isotropic directions, then search beam/random stage dictionaries to see whether noncommuting folds can genuinely compress the full `k=10` ball chain
- Outcome: negative constructive evidence; the saved 3-stage winner has refined worst radius `5.940129794736341`, a bounded 4-stage extension improves this only to `5.844588636564803`, and dense 2048-point validation over seven seeds leaves the best 4-stage candidate at worst radius `5.86484829393902` and worst ratio `0.008452260414214097`. The best networks still overwhelmingly reuse `q1/p1` folds and only add a tiny late `q2` correction, while more strongly cross-plane candidates are worse

## attempt_017

- Family: cotangent-lift spiral folding Hamiltonian search
- Mathematical modes: calculus, dynamical systems, geometric
- Core idea: use smooth Hamiltonians of the explicit form `H(q,p,t) = <p, X_t(q)>`, so the phase flow is the cotangent lift of a smooth folding isotopy on the `q`-plane; search spiral-ribbon vector fields inspired by symplectic folding, emit a real `Hamiltonian.py`, and audit that module directly
- Outcome: negative constructive evidence; after fixing an off-center support bug and adding an explicit post-fold translation stage, the best dense search candidate still has worst radius `5.884403405260652`, and the generated `Hamiltonian.py` audits at radii `5.608910263565111` to `5.67220344010714` with min ratio `0.009660348728953829`. The family is smooth and explicit but pure cotangent-lift folding appears limited by the q-compression versus p-fiber expansion tradeoff

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
- Another naive boundary-prefix hybrid `9 rigid + 1 strip` decomposition that only gives a path at `m=6` and does not yet solve the exact source equal-area strip fit or the larger-`m` scaling problem
- Another boundary-only exact recurrence solve inside the same piecewise toric affine chart class; attempt 011 showed that exact area matching is not the bottleneck because angle gluing already fails on every full-torus internal edge
- Another facewise fiber-preserving generating-function repair over the same affine action map when the complex still has full-rank internal edges; attempt 012 showed this entire cotangent-lift family cannot fix the saved exact strip
- Another two-piece rigid axis-seam origami search with a single straight cut; attempt 013 already exact-searched all such coordinate-axis seams for `m = 3,4,5,6,7,8` in the near-full `k=10` regime and found only overlapping pairs
- Another rigid axis-seam origami construction of any tree/leaf type whose internal seams land on target boundary edges; attempt 014 showed that any such seam forces overlap locally because each boundary edge has a unique incident target unit triangle
- Another slice-preserving affine transvection fold with only a common postcomposition linear symplectic mixer; attempt 015 showed that the preserved symplectic 2-plane slice forces `R >= delta^(-1/2) >= 1`, and generic mixing only worsens the lower bound
- Another multistage one-sided transvection network built from repeated half-space folds; attempt 016 pushed this family to a bounded 4-stage search and dense validation, yet the best worst-case ratio stayed near `0.00845`
- Another pure cotangent-lift folding family that only folds the q-plane and accepts the induced inverse-transpose action on p; attempt 017 made this family smooth and explicit but still stalled near ratio `0.01`

## Most promising next family

- A genuinely four-dimensional explicit folding construction, such as direct ellipsoid/polydisk folding or another coupled q-p Hamiltonian architecture that is not constrained to be a cotangent lift of a base-plane diffeomorphism
