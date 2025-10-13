
"""
charge_transfer — refactored toolkit for charge-flow analysis & visualization.

Main entry point (CLI):  python -m charge_transfer  --help
or                       python -m charge_transfer --xyz file.xyz --bader ACF.dat

This package provides:
- Robust CLI (cli.py) with sensible defaults and friendly help.
- Core algorithms (core.py): IO, graph, linear system, global vector field.
- 2D viz (viz_2d.py): projected field with streamlines/contours, auto cation/anion pivots.
- 3D viz (viz_3d.py): arrows/spheres, plane overlay, picking, resultants, snowball slider.
- Surfaces (surfaces.py): convex hull / alpha surfaces, flux accumulation, unwrap UV maps.
- Exports (exports.py): modified .xyz with labels, KIND helper (for CP2K etc.).
- Utils (utils.py): small helpers shared by modules.

Version: 0.4.0
"""
__all__ = [
    'cli', 'core', 'viz_2d', 'viz_3d', 'surfaces', 'exports', 'utils'
]

__version__ = "0.4.0"