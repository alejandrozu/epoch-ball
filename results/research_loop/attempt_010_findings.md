# Attempt 010 Findings

## Method family

- Hybrid set/logic strip-complex decomposition on the triangular grid.
- Mathematical modes: set / logic, algorithmic.
- Core idea: for `k=10`, split the near-full target set from attempt 007 into an inner rigid side-`3m` triangle, which tiles exactly into 9 rigid side-`m` triangles, plus one remaining `m^2`-cell boundary region for the 10th ball. This avoids asking all 10 balls to live in one connected target shape.

## What was implemented

- `src/hybrid_boundary_strip_family.py`
  - Builds the target triangular grid at `n = ceil(m sqrt(10))`.
  - Extracts the first `9 m^2` cells as the inner rigid side-`3m` triangle.
  - Extracts the next `m^2` cells as the boundary candidate for the 10th ball.
  - Audits the boundary candidate for connectivity, degree histogram, path status, and simple Hamiltonian-path obstructions.
  - Emits the explicit 9 macro side-`m` triangles inside the inner rigid region.
- `src/hybrid_strip_realization.py`
  - Extracts the exact simplicial strip complex of the saved `m=6` boundary path.
  - Distinguishes naive open-side corner placements from genuinely non-collinear closed-side corner placements on the source triangle boundary.
  - Uses a PyTorch boundary-spacing optimizer to test whether the source triangle can realize the same 36-face strip triangulation with positive face areas near the target area `1/72`.

## Main results

- Exact near-full target decomposition at `m=6`, `n=19`:
  - Ratio stays `0.997229916897507`, the same near-full density from attempt 007.
  - The first `324 = 9 * 36` target cells form the exact inner side-18 triangle, which tiles into 9 rigid side-6 triangles.
  - The remaining `36` target cells form a connected boundary region with degree histogram `{1: 2, 2: 34}`.
  - So the 10th-ball boundary region is literally a path graph of 36 tiny triangles, not just connected.
- Scaling sweep for `m = 6, 12, 18, 24, 30`:
  - The same `9 + 1` split persists exactly.
  - The boundary region stays connected in every tested case.
  - But it ceases to be a path immediately after `m=6`: degree-1 counts become `4, 9, 16, 25` and degree-3 counts become `68, 209, 422, 707` for `m = 12, 18, 24, 30`.
  - So the naive boundary-prefix choice gives an exact strip path only in the base `m=6` case and becomes strongly branching as `m` grows.
- Source-side strip realization audit at `m=6`:
  - The 36-cell target path is a boundary-only simplicial strip with `38` boundary vertices and `36` faces; in boundary-cycle indices the faces follow the deterministic pattern
    `(17,18,19), (17,19,20), (16,17,20), ... , (0,36,37)`.
  - A naive open-side corner test counts `140` candidate source corner triples, but most are fake because they still force collinear triangles when a corner is counted on its adjacent sides.
  - After the stricter closed-side collinearity test, only `36` source corner triples survive, all of the form `(a,18,37)` or `(18,b,37)`.
  - The saved optimizer artifact `results/research_loop/attempt_010_strip_realization_m6.json` tested the strongest closed-feasible branch first and found a strictly positive source realization with corners `[9,18,37]`, side counts `[9,19,10]`, and side assignment `[0,1,2]`.
  - In that best positive realization all `36` source face areas stay positive, with minimum `0.0027991314071759907`, maximum `0.024940329133336117`, target `1/72 = 0.013888888888888888`, mean absolute error `0.005915893299600265`, and MSE `4.652374300853456e-05`.
  - So the `m=6` source strip is not combinatorially impossible, but the current boundary-only optimizer does not yet achieve exact equal-area matching.

## Consequence for the task

- This is the first saved family in the repo that realizes the near-full `m=6` target arithmetic with 9 balls already handled by explicit rigid side-6 triangles and only one ball left to solve by a non-simplicial path-strip map.
- However, the family is not yet a correct solution:
  - the source-side non-simplicial map for that 10th ball has not yet been exactified into an equal-area continuous piecewise-affine map;
  - the naive source-side rhombus pairing route is parity-obstructed, since the side-`m` source triangle has unequal up/down tiny-triangle counts;
  - the best boundary-only source strip fit still has face-area error on the order of `1e-2`, so it is not yet a symplectic/area-preserving witness;
  - and the simple boundary-prefix region branches for larger `m`, so this exact path-strip target does not yet scale with `epsilon`.

## Strongest next move

- Keep the hybrid `9 rigid + 1 non-simplicial` idea, but first solve the concrete `m=6` source strip exactly, either by an analytic equal-area boundary recurrence or by allowing one extra interior vertex chain instead of a boundary-only realization.
- In parallel, replace the naive boundary-prefix region by a searched annular subset whose dual graph stays a path or another explicitly gluable strip/tree complex for larger `m`.
