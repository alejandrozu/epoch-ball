# Attempt 002 Literature Notes

Semantic Scholar MCP search was attempted first for this session but returned HTTP 429 from the runtime proxy. I therefore used fallback web retrieval for the key references below.

## 1. Equal ball packings in the 4-ball

- Source: McDuff, "Symplectic embeddings of 4-dimensional ellipsoids", Mathematische Annalen (2025), theorem statement / table in the article preview.
- Link: https://link.springer.com/article/10.1007/s00208-025-03221-7
- Why it matters: the preview states the equal-ball packing numbers for the 4-ball, including
  - `v(B^4,8) = 288/289`
  - `v(B^4,k) = 1` for all `k >= 9`
- Relevance here: this means the target problem is abstractly feasible for every `k >= 10`; the failure of the repository's polynomial Hamiltonian is not explained by a packing obstruction.

## 2. Cremona reduction and reduced vectors

- Source: Karshon and Kessler, Journal of Symplectic Geometry 15(4) (2017), theorem 1.9 / theorem 6.3 in the preview PDF.
- Link: https://intlpress.com/site/pub/files/_fulltext/journals/jsg/2017/0015/0004/JSG-2017-0015-0004-a005.pdf
- Why it matters: the preview states that a reduced vector with positive square encodes a blowup form. This is the reduction framework behind the ball-packing feasibility test.
- Relevance here: for equal weights `(mu; 1^k)` with `k >= 10`, any task target `mu = sqrt(k/(1-eps))` satisfies `mu > sqrt(k) > 3`, so the vector is already reduced and has positive square. This matches the full-packing statement above.

## 3. Packing stability threshold

- Source: Obuse PDF "Packings of equal balls" (contains statement `N_stab(CP^2)=9`).
- Link: https://math.indianapolis.iu.edu/~obuse/Papers/007-packings.pdf
- Why it matters: the note summarizes the packing stability threshold for `CP^2`, consistent with the `k >= 9` full-packing regime.
- Relevance here: provides an independent secondary pointer that the qualitative threshold near 9 is the right one.
