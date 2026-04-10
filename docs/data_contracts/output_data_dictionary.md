# Output data dictionary

## Scope
Data handling reference for headless structural and cold-plate thermal FEA evaluation.

## Relevant sources
- Input/measurement source: inspect `solvers/coldplate_cht_eval.py` and the repository's existing data files.
- Consumer/algorithm: `solvers/ccx_structural_grader.py`.
- Domain concerns: mesh convergence, element order, boundary conditions, solver exit status, thermal resistance, O-ring gland.

## Contract requirements
- Preserve source values and original units.
- State column/field names, types, units, valid ranges, and whether a field is optional.
- Define timestamp origin, sampling assumptions, and ordering requirements where applicable.
- Distinguish missing, invalid, saturated, and intentionally masked values.
- Keep calibration identifiers and acquisition metadata alongside measurements.
- Reject malformed records with an actionable diagnostic; do not silently invent defaults.

## Provenance
For every derived dataset, record its source, transformation script, parameters, timestamp, and checksum. Mark simulated or synthetic data clearly and never present it as a measured capture.

## Review status
Template for project-specific completion. Confirm exact field names and constraints against the actual parser and checked-in fixtures before publishing a formal schema.
