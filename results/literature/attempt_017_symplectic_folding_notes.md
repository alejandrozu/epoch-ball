# Attempt 017 Literature Note: Symplectic Folding

## Query context

- Goal: pick a post-attempt-016 family that is explicitly constructive and qualitatively different from toric origami and one-sided transvection folds.
- Retrieval path: Semantic Scholar MCP tools.

## Useful hits

### 1. Felix Schlenk, "On symplectic folding" (1999)

- Semantic Scholar paper id: `9a78bc87da22be3b22108a122c72a76a7c1439e4`
- External ids: arXiv `math/9903086`
- Why it matters:
  - The abstract explicitly says symplectic folding is used to construct nearly optimal embeddings of ellipsoids and polydisks into balls and cubes.
  - It also says connected finite-volume symplectic manifolds may be asymptotically filled with skinny ellipsoids or polydisks.
- Attempt-017 use:
  - Treat folding as a legitimate constructive design pattern, not just an existence theorem.
  - Motivate trying a smooth global folding Hamiltonian family rather than another piecewise toric or half-space fold ansatz.

### 2. Christianson-Nelson, "Symplectic embeddings of four-dimensional polydisks into balls" (2016)

- Semantic Scholar paper id: `d662945a38abf4a5fc19e8e5bd7b4e672d7482a3`
- DOI: `10.2140/agt.2018.18.2151`
- Why it matters:
  - The abstract states that Schlenk's folding construction makes their bound on `P(a,1) -> B(c)` optimal on the range they study.
  - This gives a concrete modern signal that folding is not merely asymptotic folklore; it can be sharp for explicit four-dimensional embedding problems.
- Attempt-017 use:
  - Supports choosing a folding-inspired family as the next constructive experiment.
  - Suggests using intermediate strip / product-like geometries as a computational guide, even if the final family is implemented directly in Euclidean coordinates.

## Interpretation for the repo

- These papers do **not** by themselves provide an explicit Hamiltonian for the repository task.
- They do justify a new research direction:
  - use smooth folding-style maps,
  - prefer explicit symplectic constructions over pure obstruction work,
  - and test whether a global folding mechanism can outperform the failed transvection family.

## Practical design decision

- Attempt 017 will not try to reimplement Schlenk's full published folding formulas directly.
- Instead it will test a nearby explicit family with the same qualitative flavor:
  - smooth folding of a long strip-like support in configuration space,
  - lifted to a genuine symplectic map by cotangent lift,
  - which gives an explicit Hamiltonian of the form `H(q,p,t) = p · X_t(q)`.
