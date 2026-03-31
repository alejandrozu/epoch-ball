# Strongest Current Path

- Score: 0.57
- Ready for paper: no
- Hypothesis: An algebraic exact-recurrence solver can remove the remaining uncertainty in the `m=6` hybrid strip: if the boundary-only source strip admits a rational equal-area realization and the induced local toric charts glue acro

## Why it is currently best

Implemented `src/hybrid_strip_recurrence_solver.py`, a new algebraic exactification tool for the hybrid `9 rigid + 1 strip` family. The script scans rational seeds up to denominator 72, finds a unique exact recurrence seed `17/36`, reconstructs the 38 source boundary vertices in exact `Fraction` arithmetic, and verifies that all 36 source faces have exact area `1/72`, matching the target strip exa

## Known limits

The exactification does not make the proposed solution correct. The toric gluing audit shows that all 35 internal shared edges have midpoints strictly inside `Delta(1)`, so each carries a full `T^2` fiber, but `same_angle_matrix_count = 0`: no adjacent pair has matching `A^{-T}`. In fact the 36 faces carry 36 distinct action matrices. Because constant angle translations cannot repair a nonzero lin

## Next move

Leave boundary-only exact strip fitting inside the same toric chart class. The next hybrid attempt should add a genuinely different angle/gluing mechanism, such as a generating-function layer, a non-toric local model, or a formal proof that any connected non-triangular boundary-only strip in this piecewise toric affine class forces global affineness and therefore cannot work.
