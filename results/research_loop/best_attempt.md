# Strongest Current Path

- Score: 0.6
- Ready for paper: no
- Hypothesis: A genuinely new axis-seam origami family may survive the attempt-012 obstruction: partition one source moment triangle into rigid affine pieces and allow derivative jumps only along seams whose common image lies on a tru

## Why it is currently best

Implemented `src/axis_seam_origami_search.py`, an exact search for the simplest nontrivial axis-seam origami family. The script enumerates the 12 rigid triangular-lattice automorphisms, splits a source side-`m` moment triangle into two connected rigid pieces along every internal straight cut `i=c`, `j=c`, or `i+j=c`, and tests all rigid placements of the two pieces inside the target side-`n` simpl

## Known limits

The new axis-seam family still does not make the proposed solution correct. In the bounded two-piece rigid setting, coordinate-axis seams never produce a genuine fold inside the target simplex. Every matched pair of rigid placements that sends the common cut to `x1=0` or `x2=0` overlaps in at least one target unit triangle, so no interior-disjoint two-piece axis-seam origami exists in the tested n

## Next move

Stay on the axis-seam line and move to a stronger structural step: either prove a general overlap lemma for any rigid leaf piece attached along a coordinate-axis seam, or extend the search to multi-piece axis-seam trees to see whether a nontrivial leaf attachment can ever avoid the same overlap mechanism.
