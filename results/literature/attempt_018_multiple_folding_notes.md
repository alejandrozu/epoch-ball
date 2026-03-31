# Attempt 018 Literature Note: Multiple Folding

## Query context

- Goal: motivate a recursive four-dimensional folding family that is not a cotangent lift and not a half-space transvection network.
- Retrieval path: Semantic Scholar MCP tools first, with runtime-proxy fallback available if needed.

## Useful evidence recovered

### 1. Felix Schlenk, "Multiple symplectic folding in four dimensions" (2005)

- Semantic Scholar paper id: `35a3ceb997813c5b4174758c5956f08e0df6368d`
- DOI: `10.1515/9783110199697.52`
- Why it mattered for attempt 018:
  - The title itself is enough to justify a move from single-pass folding to a recursive or multistage fold tree.
  - It supports testing a construction whose main structural idea is repeated local folds rather than one global bend.

### 2. Felix Schlenk, "Symplectic folding in higher dimensions" (2005)

- Semantic Scholar paper id: `8286415a6bb7a1b09868576e8413bd4c31cd7c5a`
- Why it mattered for attempt 018:
  - The title suggests that genuine higher-dimensional folding, not only base-plane folding, is a legitimate design target.
  - That directly motivated leaving attempt 017's cotangent-lift restriction and trying a coupled `q,p` construction.

## Retrieval limitations

- MCP calls during attempt 018 hit intermittent `429 Too Many Requests` limits.
- The first paper lookup above succeeded; the second lookup and snippet search were rate-limited during the run.
- Because of that, this note is intentionally modest: it records the structural guidance recovered from the literature rather than claiming detailed theorem transfer into the repo.

## Design consequence for the repo

- Attempt 018 adopted a recursive pairing tree with repeated local fold operations.
- The concrete repo family was:
  - stage A: use the first symplectic plane to create a temporary branch signal in the second plane,
  - stage B: use that second-plane signal to collapse the first plane toward each pair midpoint,
  - optional stage C ablation: partially return the temporary second-plane displacement.
- This is not Schlenk's published formula. It is a nearby algorithmic interpretation of the multiple-folding idea.
