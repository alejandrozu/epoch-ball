# Research Context

- Stage: research_loop
- Note: session 15, starting attempt_015
- Agent: ARC1
- Model: gpt-5.4
- Reasoning effort: xhigh
- Semantic Scholar access: runtime-configured proxy + Codex MCP tools
- Semantic Scholar env vars: `SEMANTIC_SCHOLAR_API_BASE_URL`, `SEMANTIC_SCHOLAR_GRAPH_API_URL`, `SEMANTIC_SCHOLAR_RECOMMENDATIONS_API_URL`
- Semantic Scholar auth token available: yes
- Best score: 0.62
- Ready for paper: no
- Current hypothesis: A local boundary-uniqueness obstruction may upgrade attempt 013 from a bounded search failure to a general rigid-family impossibility: if an internal seam edge of a rigid piecewise-affine origami lands on `x1=0` or `x2=0

## Latest Summary

- Worked: Implemented `src/axis_seam_local_obstruction.py`, a proof-style local obstruction script for rigid axis seams. The script enumerates the unit-triangle incidence structure of the target simplex, verifies exactly up to sid
- Did not work: This broader rigid family still does not make the proposed solution correct. Attempt 014 shows that any rigid internal seam mapped to a true coordinate axis in the target simplex forces overlap immediately, edge-by-edge.
- Next step: Leave rigid axis-seam origami entirely. The next live family must either use non-rigid within-piece maps, a non-fiber-preserving/base-mixing local model, or a construction whose discontinuity set does not map to boundary

## Strongest Current Path

- Hypothesis: A local boundary-uniqueness obstruction may upgrade attempt 013 from a bounded search failure to a general rigid-family impossibility: if an internal seam edge of a rigid piecewise-affine origami lands on `x1=0` or `x2=0
- Worked: Implemented `src/axis_seam_local_obstruction.py`, a proof-style local obstruction script for rigid axis seams. The script enumerates the unit-triangle incidence structure of the target simplex, verifies exactly up to sid
- Limits: This broader rigid family still does not make the proposed solution correct. Attempt 014 shows that any rigid internal seam mapped to a true coordinate axis in the target simplex forces overlap immediately, edge-by-edge.
