# Strongest Current Path

- Score: 0.59
- Ready for paper: no
- Hypothesis: A genuinely new axis-seam origami family may survive the attempt-012 obstruction: partition one source moment triangle into rigid affine pieces and allow derivative jumps only along seams whose common image lies o

## Why it is currently best

Implemented `src/axis_seam_origami_search.py`, a new exact search for the simplest nontrivial post-attempt-012 axis-seam family. The script enumerates the 12 rigid triangular-lattice automorphisms, splits a source side-`m` triangle into two connected rigid pieces along every internal straight cut `i=c`, `j=c`, or `i+j=c`, and checks all rigid placements of both pieces in the target side-`n` simplex with `n = ceil(m * sqrt(10))`. The saved artifact `results/research_loop/attempt_013_axis_seam_two_piece_search.json` shows zero valid nonoverlapping pairs for every tested `m = 3,4,5,6,7`, even though the number of matching axis-line candidates grows substantially with `m`.

## Known limits

This broader family still does not make the proposed solution correct. The search only covers the two-piece rigid axis-seam case, not full multi-piece axis-seam trees or non-rigid / base-mixing local models. What it does show is that coordinate-axis seams alone are not enough to produce even the simplest rigid fold inside the standard simplex: every matched pair of axis-aligned rigid placements still overlaps in at least one target unit triangle.

## Next move

Stay on the axis-seam line and either prove a general overlap lemma for any rigid leaf piece attached along a coordinate-axis seam, or extend the exact search to multi-piece axis-seam trees to see whether any nontrivial leaf attachment can ever avoid the same overlap mechanism.
