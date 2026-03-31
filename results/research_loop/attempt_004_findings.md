# Attempt 004 Findings

## Method family

- Compactly supported local affine symplectic pieces.
- Mathematical modes: linear algebra, combinatorial search, inequalities/bounding.
- Core idea: replace the toric/nonlinear chart mechanism with explicit affine symplectic scalings and translations on each ball, grouped into blocks and packed via axis-aligned q-plane rectangles.

## What was implemented

- `src/affine_rectangular_packing.py`
  - Defines affine symplectic block maps
    - `(q1, q2, p1, p2) -> (a q1, b q2, p1 / a, p2 / b)`
    - followed by thin-direction grid translations in `(p1, p2)`.
  - Searches for disjoint q-plane rectangle placements of these blocks inside a q-disk.
  - Computes both
    - a conservative corner-based radius bound,
    - and an exact radius for each affine-image ball via a trust-region calculation.
  - Verifies the resulting direct map on sampled boundary points.

## Main results

- `k=10`
  - Best found affine layout in the searched family:
    - partition `1 + 1 + 1 + 3 + 4`
    - exact radius `4.0699761526`
    - exact ratio `0.0364445458`
    - sampled ratio `0.0365906280`
  - Artifact: `results/research_loop/affine_rectangular_packing_k10.json`
- `k=16`
  - Best found affine layout in the searched family:
    - partition `1 + 1 + 2 + 12`
    - exact radius `5.1231056256`
    - exact ratio `0.0232266589`
    - sampled ratio `0.0233886477`
  - Artifact: `results/research_loop/affine_rectangular_packing_k16.json`

## Consequence for the task

- This family is not competitive.
- For `k=10`, its exact ratio `0.03644` is about
  - `2.94x` better than the repository's failed polynomial Hamiltonian audit ratio `0.012385`,
  - but only about `5.26%` of the constructive toric ratio `0.69252` from attempt 003.
- The failure persists in the square case:
  - the toric family packs `k=16` exactly with ratio `1`,
  - while this affine family only reaches about `0.0232`.
- So the missing mechanism is structural, not a search bug:
  local affine scalings grow the large axes linearly, while the successful toric family gets the crucial square-root behavior in Euclidean radius.

## Strongest next move

- Do not spend more time on purely affine local Hamiltonian pieces.
- Return to genuinely nonlinear smooth constructions:
  - regularized action-angle / slit-disk charts,
  - or non-affine folding constructions that reproduce the square-root compression mechanism.
