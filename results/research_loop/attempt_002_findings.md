# Attempt 002 Findings

## Method family

- Geometric / algebraic obstruction analysis via standard Cremona reduction for equal-ball packings in the 4-ball.
- This is intentionally different from attempt 001, which used optimization of a global polynomial Hamiltonian.

## What was implemented

- `src/ball_packing_bounds.py`
  - Computes equal-ball packing feasibility in normalized capacity coordinates.
  - Uses the standard Cremona reduction algorithm on vectors `(mu; 1, ..., 1)`.
  - Computes the infimum target capacity `mu` and corresponding volume fraction `k / mu^2`.

## What the code reproduced

- The classical small-`k` equal-ball packing fractions:
  - `k=3`: `3/4`
  - `k=5`: `20/25`
  - `k=6`: `24/25`
  - `k=7`: `63/64`
  - `k=8`: `288/289`
  - `k>=9`: full packing fraction `1`
- Artifact: `results/research_loop/reduction_table.json`

## Consequence for the task

- For `k=10`, the infimum target capacity is `mu = sqrt(10)`, i.e. full volume packing is abstractly allowed.
- For every `eps in (0,1)`, the task target `mu = sqrt(10/(1-eps))` is already reduced and has positive square.
- Therefore the target statement is consistent with known 4-dimensional packing theory.

## Consequence for the repository's proposed Hamiltonian

- The current explicit Hamiltonian is still incorrect as a general solution.
- At `k=10`, its audited effective capacity is about `28.415`, whereas the theoretical full-packing threshold is only `sqrt(10) ~= 3.1623`.
- Artifact: `results/research_loop/task_feasibility_k10.json`

## Strongest next move

- Keep the reduction tool as the abstract benchmark.
- The next construction attempt should be explicit and localized: toric/simplex or compactly supported Hamiltonian pieces, not another global low-degree polynomial ansatz.
