# Collision/contact or physical interaction edge cases

## Project context
headless structural and cold-plate thermal FEA evaluation. This scenario exercises the project's core behavior: STEP meshing, CalculiX structural grading, notch stress, buckling, cold-plate pressure drop and heat transfer.

## Objective
Define a repeatable check for **collision/contact or physical interaction edge cases** without assuming unverified hardware behavior or undocumented numerical limits.

## Relevant implementation
- Primary implementation: `solvers/ccx_structural_grader.py`
- Supporting implementation: `solvers/coldplate_cht_eval.py`
- Main engineering concerns: mesh convergence, element order, boundary conditions, solver exit status, thermal resistance, O-ring gland.

## Procedure
1. Record the source revision, toolchain/runtime versions, and configuration.
2. Prepare a controlled input fixture and preserve its original bytes.
3. Execute the relevant function or command using the documented interface.
4. Capture raw input, output, diagnostics, and run metadata.
5. Repeat with the boundary or fault condition relevant to this scenario.
6. Compare the result with an independently reviewed expectation.

## Expected evidence
- Reproducible command or test harness invocation.
- Input fixture identifier and checksum.
- Observed output and diagnostic log.
- Pass/fail rationale tied to explicit acceptance criteria.

## Acceptance criteria
The behavior is deterministic for identical inputs, invalid data is handled according to the documented contract, and no safety-relevant failure is silently converted into a valid result.

## Open items
- Add measured acceptance limits when bench evidence is available.
- Link a checked-in fixture and automated test when implemented.
- Record unresolved assumptions rather than treating them as verified facts.
