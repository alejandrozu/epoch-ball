# Attempt 006 Findings

- Family tested: witness-guided smooth kick-drift Hamiltonian factorization.
- Main implementation: `src/witness_guided_kick_drift.py`.
- Core construction:
  - represent the candidate flow as exact alternating kick segments `H(t,q,p)=sigma(t)V(q)` and drift segments `H(t,q,p)=sigma(t)T(p)`,
  - train the segment parameters against the explicit square-slot witness from `attempt_003`,
  - emit a smooth `Hamiltonian.py` using `C^\infty` time windows.
- Results:
  - Witness-guided smoke fit (`mse_weight=0.05`): discrete surrogate ratio `2.4880203276104803e-05`.
  - Radius-only smoke fit (`mse_weight=0`): discrete surrogate ratio `3.4465526940376647e-04`.
  - End-to-end RK4 audit of the emitted `Hamiltonian.py`: radius `19.72808660462524`, ratio `6.601766155665936e-05`, benchmark call time `0.3401749539998491` seconds for batch size `5000`.
- Interpretation:
  - The family is smooth and computationally admissible, but its geometry is catastrophically wrong.
  - The audited ratio is only about `0.53%` of the repository's already-bad polynomial baseline from `attempt_001`, and about `0.0106%` of the toric square witness ratio `0.625`.
  - So the failure is not a fine-tuning issue. The separable kick-drift architecture with source-indexed Gaussian kicks does not reproduce the branch-cut or square-root compression structure of the toric witness.
- Conclusion:
  - Keep the code and artifacts as a negative ablation.
  - Do not revisit this separable kick-drift fitting family.
  - The next viable family should encode the branch geometry explicitly, likely through a folding/origami or multi-chart construction rather than boundary imitation alone.
