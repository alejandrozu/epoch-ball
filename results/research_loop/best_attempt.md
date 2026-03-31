# Strongest Current Path

- Score: 0.58
- Ready for paper: no
- Hypothesis: The exact `m=6` hybrid strip is ruled out even after enlarging to the full facewise generating-function / cotangent-lift family over the same affine action map. Any viable near-full hybrid witness must move derivative jumps onto coordinate-axis seams or leave the fiber-preserving toric framework.

## Why it is currently best

Implemented `src/hybrid_generating_function_audit.py`, which takes the exact strip from attempt 011 and enlarges each face map to the most general fiber-preserving exact symplectomorphism over the same affine action map. The script derives the restricted-angle continuity rule, includes a toy axis-edge positive control showing that derivative jumps can survive only when a circle really collapses, and proves on the saved exact strip that every one of the 35 internal edges still fails the required restricted `A^{-T}` match. This upgrades the earlier obstruction from “plain toric affine charts fail” to “the whole facewise generating-function repair family fails on this strip.”

## Known limits

This attempt is still a stronger negative result rather than a constructive embedding. No Hamiltonian has been produced, and the audit only rules out fiber-preserving generating-function corrections of the saved affine action map. A future successful family may still need a different action complex, a non-fiber-preserving local model, or a fully different smooth construction.

## Next move

Stop trying to repair the saved strip inside the same fiber-preserving toric framework. The next live direction is to redesign the action complex so derivative jumps occur only along coordinate-axis seams where a torus factor really collapses, or else move to a non-toric local model where the facewise map does not preserve the action-base projection.
