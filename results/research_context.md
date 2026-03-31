# Research Context

- Stage: research_loop
- Note: session 13, starting attempt_013
- Agent: ARC1
- Model: gpt-5.4
- Reasoning effort: xhigh
- Semantic Scholar access: runtime-configured proxy + Codex MCP tools
- Semantic Scholar env vars: `SEMANTIC_SCHOLAR_API_BASE_URL`, `SEMANTIC_SCHOLAR_GRAPH_API_URL`, `SEMANTIC_SCHOLAR_RECOMMENDATIONS_API_URL`
- Semantic Scholar auth token available: yes
- Best score: 0.59
- Ready for paper: no
- Current hypothesis: A genuinely new axis-seam origami family may survive the attempt-012 obstruction: partition one source moment triangle into rigid affine pieces and allow derivative jumps only along seams whose common image lies o

## Latest Summary

- Worked: Implemented `src/axis_seam_origami_search.py`, a new exact search for the simplest nontrivial axis-seam origami family. The script enumerates the 12 rigid triangular-lattice automorphisms, splits a source side-`m` trian
- Did not work: The new family still does not make the proposed solution correct. Even when derivative jumps are restricted to true coordinate-axis seams, the standard simplex appears too one-sided to host a genuine rigid fold. For eac
- Next step: Stay on the axis-seam line and either prove a general overlap lemma for any rigid leaf piece attached along a coordinate-axis seam, or extend the exact search to multi-piece axis-seam trees to see whether any nontrivial

## Strongest Current Path

- Hypothesis: A genuinely new axis-seam origami family may survive the attempt-012 obstruction: partition one source moment triangle into rigid affine pieces and allow derivative jumps only along seams whose common image lies o
- Worked: Implemented `src/axis_seam_origami_search.py`, a new exact search for the simplest nontrivial axis-seam origami family. The script enumerates the 12 rigid triangular-lattice automorphisms, splits a source side-`m` trian
- Limits: The search only covers the two-piece rigid axis-seam case, not full multi-piece axis-seam trees or non-rigid/base-mixing local models. It nevertheless shows that coordinate-axis seams alone still do not produce a genuine rig
