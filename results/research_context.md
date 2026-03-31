# Research Context

- Stage: research_loop
- Note: session 12, starting attempt_012
- Agent: ARC1
- Model: gpt-5.4
- Reasoning effort: xhigh
- Semantic Scholar access: runtime-configured proxy + Codex MCP tools
- Semantic Scholar env vars: `SEMANTIC_SCHOLAR_API_BASE_URL`, `SEMANTIC_SCHOLAR_GRAPH_API_URL`, `SEMANTIC_SCHOLAR_RECOMMENDATIONS_API_URL`
- Semantic Scholar auth token available: yes
- Best score: 0.58
- Ready for paper: no
- Current hypothesis: The exact `m=6` hybrid strip cannot be repaired by any facewise fiber-preserving generating-function correction over the same affine action map. In `B^4`, only coordinate-axis edges collapse torus directions; the diagonal wall does not. Since all 35 internal edges of the exact strip keep both source and target circles active, every edge would require full `A^{-T}` agreement, which never occurs. The next viable family must either force derivative jumps onto coordinate-axis seams or leave the fiber-preserving toric framework.

## Latest Summary

- Worked: Implemented `src/hybrid_generating_function_audit.py`, a new transform-method audit for the exact `m=6`, `k=10` hybrid strip. The script enlarges each face map to the most general fiber-preserving exact symplectomorphism over the saved affine action map, derives the restricted-angle continuity rule, includes a toy axis-edge positive control, and saves the exact audit to `results/research_loop/attempt_012_generating_function_audit.json`.
- Did not work: The generating-function enlargement still cannot repair the exact strip. For `B^4`, only the coordinate axes collapse circles; the diagonal wall does not. All 35 internal edges of the strip have source and target active angle sets `[0,1]`, so the full `2 x 2` restricted matrix `A^{-T}` must match across each edge. It never does.
- Next step: Stop trying to repair the saved strip inside the same fiber-preserving toric framework. The next live direction is to redesign the action complex so derivative jumps occur only along coordinate-axis seams, or else move to a non-toric local model where the facewise map does not preserve the action-base projection.

## Strongest Current Path

- Hypothesis: The exact `m=6` hybrid strip is now ruled out even after enlarging to the full facewise generating-function / cotangent-lift family over the same affine action map; any viable near-full hybrid witness must move derivative jumps onto coordinate-axis seams or leave the fiber-preserving toric framework.
- Worked: Implemented `src/hybrid_generating_function_audit.py`, which derives the restricted-angle continuity rule for the cotangent-lift local model, includes an axis-edge positive control, and proves on the saved exact strip that every one of the 35 internal edges still fails the required restricted `A^{-T}` match.
- Limits: This is still a negative result. It rules out the current strip more strongly but does not yet produce a working embedding or Hamiltonian, and it does not address how to construct a large near-full complex whose derivative jumps lie only on genuine coordinate-axis collapse edges.
