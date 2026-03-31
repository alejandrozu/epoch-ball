# Attempt 018 Ablation Notes

## 1. Compression-sign bug (discarded implementation state)

- Family state: first implementation of the hierarchical braid/control-translation search before fixing the sign in stage B.
- What happened:
  - the second-plane branch signal translated each sibling **away** from its pair midpoint instead of toward it.
  - this made the construction anti-fold in the control plane.
- Observed outcome from the shell run:
  - best dense radius about `14.627626842692163`
  - worst dense radius about `14.661673224981904`
  - best dense ratio about `2.1842608865790334e-4`
- Interpretation:
  - this was a real implementation bug, not evidence about the intended family.

## 2. Partial-return stage ablation

- Family state: after fixing the sign bug, add a third per-level stage that tries to return part of the temporary `(q2,p2)` displacement after reading it out.
- Search setting:
  - `return_frac in {0.75, 1.0}`
- Best saved shell outcome:
  - best dense radius `12.354363330928901`
  - worst dense radius `12.454715489557286`
  - best dense ratio `4.2925812315411713e-4`
  - worst dense ratio `4.15589680864139e-4`
- Interpretation:
  - partially undoing the temporary workspace shift made the family worse, not better.
  - the best candidate in this ablation used `return_frac = 0.75`, but it still underperformed the no-return family.

## 3. Sharp-readout local sweep

- Family state: no-return family, but broaden the search around smaller shifts and smaller `target_sign_scale`.
- Saved artifact:
  - `results/research_loop/attempt_018_hierarchical_braid_search.json`
- Best dense result in that saved search:
  - best dense radius `11.917169580518717`
  - worst dense radius `12.029543054273557`
  - best dense ratio `4.958011417233913e-4`
  - worst dense ratio `4.775330963760927e-4`
- Important pattern:
  - even in this sharper-readout search, the best candidate still selected the largest tested `target_sign_scale = 0.55` instead of the sharper `0.15` or `0.3`.
  - this indicates that simply sharpening the branch detector does not remove the family's need for large temporary auxiliary-plane motion.
