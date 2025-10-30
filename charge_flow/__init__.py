"""
charge_flow — atomic charge-flow field analysis & visualization.

`charge_flow` is a Python-based toolkit for analyzing and visualizing
atomic charge-flow fields derived from per-atom charge data (e.g. Bader
analysis). It combines an .xyz structure with per-atom charges to show how
electron density flows between atoms, both globally and locally.

Typical usage (CLI):
    python -m charge_flow --xyz STRUCTURE.xyz --bader ACF.dat

This will launch the default 3D visualization:
- arrows anchored at atoms,
- arrow length encodes charge-flow magnitude,
- direction shows where charge is flowing.

Main capabilities:
- Full-system 3D charge-flow field ("field" mode), with arrows or spheres.
- Per-atom inspection ("atom" mode) for selected indices.
- 2D projected views onto Cartesian planes (xy, yz, xz, pca) or Miller-index
  planes via --plane-normal / --slab.
- Global / cumulative resultant vectors:
  * --show-total-resultant      : one arrow = sum of all flows.
  * --show-snowball             : interactive radius slider that accumulates
                                  contributions outward from the center of mass.

Command-line driven:
All behavior is controlled via CLI flags (no need to edit the code). See
    python -m charge_flow -h
for full help, including:
    --mode field / --mode atom
    --view 3d / --view 2d / --view both
    --display arrows / --display spheres
    --scale-3d, --sphere-scale
    --plane / --plane-normal / --slab
    --fig-size / --fig-lims
    --idw-radius / --idw-power
    --show-total-resultant / --show-snowball

Environment / installation (recommended workflow):
1. Clone the repository and enter it.
2. Create and activate a conda env (using conda-forge):
       conda create -n charge_flow python=3.10 numpy scipy matplotlib \\
                    pyvista vtk networkx rdkit -c conda-forge
       conda activate charge_flow
3. Install the package into that environment:
       pip install .
   (or `pip install -e .` for development)

Internal modules:
- cli.py       : argument parsing and user-facing command-line interface.
- core.py      : IO, graph construction, linear system for charge flow, and
                 computation of per-atom vectors / resultants.
- viz_3d.py    : 3D visualization (arrows, spheres, interactive snowball slider).
- viz_2d.py    : 2D projection plots, plane slicing, IDW background maps.
- exports.py   : helpers to export modified .xyz files, relabeled atoms, etc.
- utils.py     : shared helpers (geometry, math, color, etc.).
- surfaces.py  : surface- and hull-related helpers (internal; not all features
                 are exposed in the public CLI yet).

Version: 1.0.0
"""

__all__ = [
    "cli",
    "core",
    "viz_2d",
    "viz_3d",
    "surfaces",
    "exports",
    "utils",
]

__version__ = "1.0.0"
