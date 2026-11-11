# physics-fea-bench

Headless structural and thermal FEA for hardware that has to survive shipping.
No GUI, no license server, no engineer babysitting a solve.

## What is here

`solvers/ccx_structural_grader.py`: STEP in, quadratic C3D10 mesh out, CalculiX
`ccx` solve, `.frd` nodal stress parse. Reports von Mises peak notch stress,
stress concentration factor Kt, and the first N linear buckling multipliers from
`*BUCKLE`. Mesh sizing is deterministic, so two runs of one part are diffable.

`solvers/coldplate_cht_eval.py`: 1,200 W direct-to-chip microchannel cold plate.
Conjugate heat transfer resistance stack-up, Darcy-Weisbach pressure drop with
Dean-vortex bend correction, and an ISO 3601 EPDM O-ring gland check.

`benchmarks/sample_bracket_metrics.json`: ground truth against a failure trap.
The trap is the point. Linear C3D4 tets under-predict the notch peak by 39.5%,
Kt reads 1.74 instead of 2.87, and every gate passes. That bracket cracks in
fatigue. A grader that agrees with the trap is visibly wrong.

## Run

```
pip install gmsh numpy
python solvers/ccx_structural_grader.py part.step --nominal-mpa 41.2
python solvers/coldplate_cht_eval.py
```

`ccx` must be on PATH for the structural path. The cold plate evaluator is pure
Python.

## Why the traps matter

Most FEA regressions are not solver bugs. They are element order, mesh density at
a fillet, or a boundary condition that quietly removes the load path. All three
pass a convergence check. The benchmark encodes the wrong answer next to the
right one so the grader is scored on catching it, not on matching it.
