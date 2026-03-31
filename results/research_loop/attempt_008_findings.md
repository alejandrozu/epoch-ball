# Attempt 008 Findings

- Family: graph/symmetry reduction of continuity-preserving simplicial origami to rigid lattice-triangle packing.
- Local rigidity result: once the first unit triangle is fixed, every later row is forced. The saved rigidity audit confirms affine, vertex-injective propagation for all 12 admissible seeds and all `m <= 20`.
- Consequence: a continuous simplicial unit-triangle embedding of the subdivided source triangle is not a folded polyiamond at all; it is just a rigid side-`m` lattice triangle.
- Exact packing result for `k=10`: rigid side-`m` triangles fit exactly for `m=1` and `m=2`, but the exact maximum is only `9` for `m=3,4,5,6,7,8,9,10,11`.
- Key example: at `m=6`, the perfect area ratio would be `0.997229916897507`, but the exact graph search caps the rigid family at ratio `0.8975069252077562`.
- Heuristic extension: the saved mixed sweep still finds only `9` placements at `m=12`, with ratio `0.8975069252077562`.
- Interpretation: attempt_007 failed because arbitrary permutations destroyed continuity; attempt_008 shows that restoring continuity inside the unit-triangle simplicial class overcorrects and makes the single-ball image rigid. The next viable family has to be non-simplicial.
