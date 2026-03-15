# Benchmark methodology

## Scope
This note records the design considerations for headless structural and cold-plate thermal FEA evaluation.

## Current project context
The repository currently implements STEP meshing, CalculiX structural grading, notch stress, buckling, cold-plate pressure drop and heat transfer. The main code paths are `solvers/ccx_structural_grader.py` and `solvers/coldplate_cht_eval.py`. The principal engineering concerns are mesh convergence, element order, boundary conditions, solver exit status, thermal resistance, O-ring gland.

## Design questions
- What inputs, units, coordinate conventions, and preconditions does this part of the system require?
- Which state is authoritative, and which values are derived?
- How are invalid, stale, incomplete, or physically implausible inputs represented?
- Which assumptions are supported by source code or measured evidence, and which remain provisional?
- What observable evidence allows a reviewer to reproduce the behavior?

## Recommended design record
Document the selected approach, alternatives considered, rationale, expected failure modes, and any measurable trade-offs. Keep numerical constants linked to source or calibration evidence; do not copy unexplained values into a second source of truth.

## Verification
Pair this design note with the relevant scenarios under `docs/verification/scenarios/`. A design decision is not considered verified merely because the code compiles.

## Status
Working reference. Confirm implementation details against the source before using this note as a hardware, clinical, or production specification.
