# Attempt 009 Findings

## Method family

- Probabilistic / exact search over non-simplicial connected polyomino action-shapes with integer toric covers.
- Mathematical modes: probabilistic, combinatorial, discrete <-> continuous mapping.
- Core idea: replace the rigid simplicial-triangle image of attempt 008 by a connected equal-area polyomino in moment space, then apply an integer matrix cover of determinant at least `k` and measure the smallest target simplex `Delta(mu)` that contains the translated image.

## What was implemented

- `src/strip_snake_cover_search.py`
  - Searches connected polyomino shapes generated from self-avoiding walks.
  - Enumerates integer cover matrices in a determinant window and evaluates the exact simplex capacity `mu` from transformed cell-corner data.
  - Supports exact enumeration for fixed `n_cells`.
  - Supports a `--path-only` mode for the stricter strip-snake subfamily whose cell-adjacency graph is a path.
  - Fixes a seed bug discovered in-session: the old "staircase" seed was disconnected and has been replaced by a connected alternating right/up path.

## Main results

- Exact unrestricted `n_cells=10`, `entry_bound=6`:
  - Best ratio `0.6920415224913494`.
  - Best capacity `mu = 3.8013155617496426`.
  - Best determinant `10`.
  - This is slightly below the constructive baseline from attempt 003, which had ratio `0.6925207756232687` at `mu = 3.8`.
- Exact unrestricted `n_cells=10`, `entry_bound=8`:
  - Best ratio stayed exactly `0.6920415224913494`.
  - The best witness was unchanged from the `entry_bound=6` exact run.
  - So enlarging the integer-matrix search range from `2280` to `4584` matrices did not improve the family.
- Exact path-only `n_cells=10`, `entry_bound=6`:
  - Best ratio `0.5540166204986149`.
  - Best capacity `mu = 4.2485291572496005`.
  - So the cleaner strip-snake subfamily is much worse than the unrestricted connected-polyomino surrogate.

## Consequence for the task

- This family does not improve the repository's best constructive ratio.
- The best unrestricted witness found here is also not a true strip snake: its cell-adjacency graph branches, so the strongest current result in this family depends on a more general connected polyomino geometry than the original strip-fold intuition.
- The session still does not produce the required explicit smooth Hamiltonian on `R^4 x [0,1]`; it only provides a surrogate search over candidate moment-space shapes.

## Strongest next move

- Treat this family as a useful near-match / negative result and move on to a different non-simplicial construction that is not based on a single connected polyomino action-shape.
- In particular, the next family should not rely on a single connected action-shape whose best exact witness already plateaus below attempt 003 and whose realizable path-like subfamily collapses to ratio `0.5540`.
