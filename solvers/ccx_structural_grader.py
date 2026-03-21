"""Headless STEP -> Gmsh C3D10 -> CalculiX (ccx) -> structural metrics.

Import STEP, heal, mesh quadratic tets with fillet refinement, emit a ccx deck
with a *STATIC step and a *BUCKLE step, run ccx headless, parse the .frd nodal
stress block. Reports von Mises peak notch stress, stress concentration factor
Kt, and the first N linear buckling eigenvalues. Deterministic, no GUI.
"""

from __future__ import annotations

import math
import re
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import gmsh

MIN_ELEM_MM = 0.20
CHORD_TOL_MM = 0.02
REFINE_TAGS = ("FILLET", "NOTCH", "HOLE")


@dataclass(frozen=True)
class MeshSpec:
    global_size_mm: float = 1.60
    order: int = 2  # C3D10

    def size_for_radius(self, radius_mm: float) -> float:
        """Chord-error sizing on a circular arc."""
        chord = math.sqrt(8.0 * radius_mm * CHORD_TOL_MM)
        return max(min(chord, self.global_size_mm), MIN_ELEM_MM)


def mesh_step(step: Path, work: Path, spec: MeshSpec = MeshSpec()) -> tuple[Path, int]:
    """Mesh a STEP solid to quadratic tets; returns (msh path, node count)."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.model.add("part")
        gmsh.model.occ.importShapes(str(step))
        gmsh.model.occ.synchronize()
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMax", spec.global_size_mm)
        gmsh.option.setNumber("Mesh.MeshSizeMin", MIN_ELEM_MM)
        gmsh.option.setNumber("Mesh.ElementOrder", spec.order)
        gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT
        for dim, tag in gmsh.model.getEntities(2):
            if not any(k in gmsh.model.getEntityName(dim, tag).upper() for k in REFINE_TAGS):
                continue
            dist = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", [tag])
            thr = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(thr, "InField", dist)
            gmsh.model.mesh.field.setNumber(thr, "SizeMin", MIN_ELEM_MM)
            gmsh.model.mesh.field.setNumber(thr, "SizeMax", spec.global_size_mm)
            gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.15)
            gmsh.model.mesh.field.setNumber(thr, "DistMax", 1.20)
            gmsh.model.mesh.field.setAsBackgroundMesh(thr)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Netgen")
        out = work / "part.msh"
        gmsh.write(str(out))
        return out, len(gmsh.model.mesh.getNodes()[0])
    finally:
        gmsh.finalize()


DECK = """*HEADING
structural grader, quadratic tets
*INCLUDE, INPUT={mesh}
*ELEMENT, TYPE=C3D10, ELSET=EALL
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.30
*DENSITY
7.85e-09
*BOUNDARY
{clamped}, 1, 3, 0.0
*STEP, NAME=STATIC
*STATIC
1., 1.
*CLOAD
{load_node}, 3, {load_n}
*NODE FILE
S
*END STEP
*STEP, NAME=BUCKLE
*BUCKLE
{n_modes}, 1e-4
*END STEP
"""


def write_deck(path: Path, mesh: Path, clamped: Sequence[int], load_node: int,
               load_n: float, n_modes: int) -> Path:
    path.write_text(DECK.format(mesh=mesh.name, clamped=",".join(map(str, clamped)),
                                load_node=load_node, load_n=load_n, n_modes=n_modes),
                    encoding="utf-8")
    return path


def run_ccx(deck: Path, timeout_s: int = 900) -> Path:
    subprocess.run(["ccx", "-i", deck.stem], cwd=deck.parent, check=True,
                   timeout=timeout_s, capture_output=True, text=True)
    frd = deck.with_suffix(".frd")
    if not frd.exists():
        raise RuntimeError(f"ccx produced no .frd for {deck.name}")
    return frd


_STRESS_RE = re.compile(r"^\s*-\s*4\s+(\d+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)")


@dataclass
class Results:
    von_mises: dict[int, float] = field(default_factory=dict)
    buckles: list[float] = field(default_factory=list)


def parse_frd(frd: Path) -> Results:
    """-4 blocks carry the stress tensor per node; eigenvalues live in .dat."""
    res = Results()
    for line in frd.read_text(encoding="utf-8", errors="replace").splitlines():
        m = _STRESS_RE.match(line)
        if not m:
            continue
        sxx, syy, szz = (float(m.group(i)) for i in (2, 3, 4))
        res.von_mises[int(m.group(1))] = math.sqrt(
            0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2))
    dat = frd.with_suffix(".dat")
    if dat.exists():
        for line in dat.read_text(encoding="utf-8", errors="replace").splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[0].isdigit():
                try:
                    res.buckles.append(float(parts[1]))
                except ValueError:
                    pass
    return res


@dataclass(frozen=True)
class Report:
    sigma_vm_peak_mpa: float
    kt: float
    buckle_multipliers: tuple[float, ...]
    n_nodes: int
    passed: bool

    def as_dict(self) -> dict:
        return {"sigma_vm_peak_mpa": round(self.sigma_vm_peak_mpa, 3),
                "kt": round(self.kt, 3),
                "buckle_multipliers": [round(b, 4) for b in self.buckle_multipliers],
                "n_nodes": self.n_nodes, "passed": self.passed}


def grade(step: Path, nominal_mpa: float, yield_mpa: float = 276.0,
          min_buckle: float = 1.5, kt_limit: float = 3.0) -> Report:
    """Mesh, solve, extract, gate."""
    with tempfile.TemporaryDirectory(prefix="fea_") as tmp:
        work = Path(tmp)
        mesh, n_nodes = mesh_step(step, work)
        deck = write_deck(work / "job.inp", mesh, clamped=(1, 2, 3),
                          load_node=4, load_n=1000.0, n_modes=4)
        res = parse_frd(run_ccx(deck))
    peak = max(res.von_mises.values(), default=0.0)
    kt = peak / nominal_mpa if nominal_mpa > 0 else float("inf")
    ok = (peak < yield_mpa and kt <= kt_limit and bool(res.buckles)
          and min(res.buckles) >= min_buckle)
    return Report(peak, kt, tuple(res.buckles), n_nodes, ok)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Headless STEP structural FEA grader")
    ap.add_argument("step", type=Path)
    ap.add_argument("--nominal-mpa", type=float, required=True)
    ap.add_argument("--yield-mpa", type=float, default=276.0)
    a = ap.parse_args()
    print(grade(a.step, a.nominal_mpa, a.yield_mpa).as_dict())
