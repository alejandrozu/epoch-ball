# Attempt 005 Findings

- Family tested: smooth nonlinear slit-chart packing via exterior log-polar canonical coordinates and rectangle slotting.
- Main implementation: `src/smooth_slit_chart_packing.py`.
- Chart sanity checks:
  - `python -m py_compile src/smooth_slit_chart_packing.py` passed.
  - For `a` near `1`, the planar chart had `area_total ≈ pi` and roundtrip max error below about `2e-11`.
- Best measured results:
  - `k=10`, near-`a=1` sweep: best at `a=1.01` with sampled radius `3.376970763778154` and sampled ratio `0.07689370674740484`.
  - `k=16`, near-`a=1` sweep: best at `a=1.01` with sampled radius `3.376970763778154` and sampled ratio `0.12302993079584774`.
- Interpretation:
  - These sampled ratios are optimistic upper bounds on the true ratio because they come from sampled boundary maxima.
  - The perfect-square case `k=16` should be the easiest benchmark, but this family still reaches only about `12.3%` of full packing, while the toric square-slot baseline from `attempt_003` reaches `1.0`.
  - The chart is smooth and symplectic on each plane, so the failure is structural rather than a remaining coding bug.
  - The likely obstruction is that the slit-chart coordinates do not preserve the toric moment-simplex relation analogous to `x1 + x2 <= 1`, so the square-slot recipe no longer produces the `sqrt(l)` radius control that made the toric family competitive.
- Conclusion:
  - Keep the artifacts as a negative ablation.
  - Do not revisit generic slit-chart rectangle packings.
  - The next candidate should either regularize toric action variables while preserving a simplex-type relation, or switch to a non-affine folding family.
