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

## Consequence for the task

- This is the first saved family in the repo that realizes the near-full `m=6` target arithmetic with 9 balls already handled by explicit rigid side-6 triangles and only one ball left to solve by a non-simplicial path-strip map.
- However, the family is not yet a correct solution:
  - the source-side non-simplicial map for that 10th ball has not been constructed;
  - the naive source-side rhombus pairing route is parity-obstructed, since the side-`m` source triangle has unequal up/down tiny-triangle counts;
  - and the simple boundary-prefix region branches for larger `m`, so this exact path-strip target does not yet scale with `epsilon`.

## Strongest next move

- Keep the hybrid `9 rigid + 1 non-simplicial` idea, but replace the naive boundary-prefix region by a searched boundary subset whose dual graph is a path or another explicitly gluable strip/tree complex for larger `m`.
- In parallel, build a source-side non-simplicial strip complex for the side-6 source triangle and test an actual continuous piecewise-affine map onto the saved 36-cell boundary path at `m=6`.
