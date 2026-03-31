# Attempt 007 Findings

- Family: topological/origami moment-tiling via tiny-triangle reassignment in toric action-angle coordinates.
- Positive result: the exact density formula works. For `k=10`, `eps=0.01`, the script chooses `m=6`, `n=19`, and achieves exact ratio `0.9972299168975068`. For `eps=1e-4`, it chooses `m=117`, `n=370`, and reaches exact ratio `0.9999269539810081`.
- Negative result: the map is not continuous. Every nontrivial `k=10` run tested had `continuous_adjacencies_ball0 = 0`, including `0/45` at `eps=0.01` and `0/20358` at `eps=1e-4`.
- Fold magnitude does not decay under refinement. The mean edge-midpoint jump stays near `0.45`, and the max jump stays around `1.0`, so this is not a small smoothing defect.
- Sampling alone is misleading for this family. On finer runs the sampled radius misses the worst fold points, so sampled ratios can exceed `1`; the exact combinatorial radius is the relevant audit.
- Interpretation: arbitrary permutation of tiny target triangles captures the desired area arithmetic but destroys the simplicial compatibility needed for a continuous map. Any viable continuation of the origami idea must enforce a consistent vertex map and connected target subcomplex for each source ball.
