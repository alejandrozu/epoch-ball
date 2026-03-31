# Strongest Current Path

- Score: 0.62
- Ready for paper: no
- Hypothesis: A local boundary-uniqueness obstruction may upgrade attempt 013 from a bounded search failure to a general rigid-family impossibility: if an internal seam edge of a rigid piecewise-affine origami lands on `x1=0` or `x2=0

## Why it is currently best

Implemented `src/axis_seam_local_obstruction.py`, a proof-style local obstruction script for rigid axis seams. The script enumerates the unit-triangle incidence structure of the target simplex, verifies exactly up to side 40 that every unit boundary edge on `x1=0`, `x2=0`, or the diagonal has one incident unit triangle while interior edges have two, and verifies up to axis-path length 20 that a bo

## Known limits

This broader rigid family still does not make the proposed solution correct. Attempt 014 shows that any rigid internal seam mapped to a true coordinate axis in the target simplex forces overlap immediately, edge-by-edge. So the entire rigid axis-seam origami strategy is impossible, not merely the specific two-piece straight-cut instances from attempt 013. Coordinate-axis seams were the only algebr

## Next move

Leave rigid axis-seam origami entirely. The next live family must either use non-rigid within-piece maps, a non-fiber-preserving/base-mixing local model, or a construction whose discontinuity set does not map to boundary edges of the target simplex.
