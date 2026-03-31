# Attempt 003 Findings

## Method family

- Explicit constructive toric/action-angle packing.
- Mathematical modes: discrete <-> continuous mapping, geometric, algorithmic.
- Core idea: build real symplectic embeddings on the source union using action-angle coordinates, first through the standard `l^2` square-slot construction and then through a search over multi-block toric moment-triangle layouts.

## What was implemented

- `src/explicit_square_packing.py`
  - Implements the standard toric square packing map on the disjoint source union.
  - Verifies the exact predicted ratio `k / ceil(sqrt(k))^2` on sampled boundary points.
- `src/toric_linear_family.py`
  - Searches one-block affine toric covering constructions coming from integer `2 x 2` matrices of determinant `k`.
  - Reports the best capacity `mu` found for each determinant.
- `src/axis_block_packing_search.py`
  - Searches for disjoint translated axis-aligned toric blocks inside the target simplex.
  - The search is constructive but heuristic because translations are restricted to a finite grid.

## Main results

- Square-slot baseline:
  - `k=10`: exact constructive ratio `0.625` at `mu=4`, confirmed numerically to machine precision.
  - `k=16`: exact full packing ratio `1.0` at `mu=4`, also confirmed numerically.
- One-block affine toric covering:
  - Exact on perfect squares in the searched range.
  - For `k=10`, the best one-block construction found requires `mu=5`, so the ratio is only `0.4`.
- Multi-block axis-aligned toric search:
  - For `k=10`, a four-block partition `1 + 2 + 3 + 4` was found inside `Delta(3.8)`.
  - This gives ratio `10 / 3.8^2 = 0.6925207756`, improving the square-slot baseline but still far from the task target for small `epsilon`.

## Consequence for the task

- This attempt produced the first explicit constructive embeddings in the repo that materially improve on the failed polynomial Hamiltonian family.
- However, the construction still does not solve the stated problem:
  - the best explicit `k=10` ratio found here is about `0.6925`, not arbitrarily close to `1`;
  - the current maps are written in singular action-angle coordinates and have not been converted into the required explicit smooth global Hamiltonian.

## Strongest next move

- Keep the toric constructive direction.
- Replace the singular polar/action-angle chart with a smooth slit-disk or Traynor-style area-preserving chart so the embedding can be promoted to a genuine smooth Hamiltonian-piece construction.
- Broaden the search from axis-aligned blocks to non-axis Delzant triangles and more flexible multi-block layouts, starting again with `k=10`.
