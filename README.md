# charge_flow Package Documentation

## Overview
The `charge_flow` package is a visualization and analysis tool for exploring charge flow between atoms in an atomistic structure. It combines per-atom charge data (e.g. from Bader analysis) with structural coordinates to produce 3D or 2D visualizations of charge flow magnitudes and directions.

## Installation

You can install and run `charge_flow` in an isolated Conda environment to ensure all dependencies work correctly. The package relies on scientific Python libraries and visualization backends such as `pyvista` and `vtk`.

### Step 1. Download / clone the repository
```bash
git clone https://github.com/jllusarc/charge_flow.git
cd charge_flow
```

### Step 2. Create and activate a conda environment (recommended)

Create an environment manually using conda-forge
```bash
conda create -n charge_flow python=3.10 numpy scipy matplotlib pyvista vtk networkx rdkit -c conda-forge
conda activate charge_flow
```

### Step 3. Install the package into the environment
Standard install:
```bash
pip install .
```

Editable/development install:
```bash
pip install -e .
```

### Step 4. Run the package
```bash
python -m charge_flow --xyz your_structure.xyz --bader ACF.dat
```

## Basic Usage
```bash
python -m charge_flow --xyz STRUCTURE.xyz --bader ACF.dat
```
This command launches the 3D visualization showing charge-flow vectors (arrows). Larger arrows represent larger charge-flow magnitudes; smaller arrows represent smaller ones.

---

## Command Reference

### 1. General Help
```bash
python -m charge_flow -h
```
Displays all available options.

### 2. Field Mode (default)
```bash
python -m charge_flow --xyz file.xyz --bader ACF.dat --mode field
```
Displays the full 3D vector field.

#### Arrow Display
```bash
--display arrows
--scale-3d N
```
- Adjusts arrow size (`N=2` default).
- `N=0` will cause the program to fail.

#### Total Resultant Vector
```bash
--show-total-resultant
--resultant-scale N
--resultant-color COLOR
```
Shows a resultant arrow at the model’s center of mass, summing all charge-flow magnitudes. Scale and color are customizable.

#### Snowball Resultant
```bash
--show-snowball
--snowball-scale N
--snowball-color COLOR
```
Displays an interactive cumulative resultant vector, growing with the sphere radius from the center of mass.

#### Sphere Display
```bash
--display spheres
--sphere-scale N
```
Shows spheres representing charge-flow magnitudes. `N=0` hides spheres (wireframe only), `N=1` enlarges them.

---

### 3. Atom Mode
```bash
python -m charge_flow --xyz file.xyz --bader ACF.dat --mode atom
```
Without selecting atoms, shows only the wireframe.

#### Select Atoms
```bash
--atoms a1 a2 a3 ... aN
```
Shows charge-flow vectors for the selected atoms (1-based indexing).
If any index exceeds total atoms, an error appears.

---

### 4. View Selection
```bash
--view 3d
--view 2d
--view both
```
- `3d`: Default interactive mode.
- `2d`: 2D projected map with arrows and dots for atoms.
- `both`: Combines both visualizations.

---

### 5. 2D Plane Options

#### Cartesian Planes
```bash
--plane {xy, yz, xz, pca, ...}
```
Sets projection plane (origin at center of mass). Must be combined with `--view 2d`.

#### Miller Indices
```bash
--plane-normal h k l
--plane-origin x y z
--slab n
```
Defines projection using Miller indices and optional custom origin or slab thickness (in Å).
Recommended slab thickness: `n=1` Å.

---

### 6. Figure Customization (2D)
```bash
--fig-size X Y
--fig-lims x1min x1max x2min x2max
```
Controls figure dimensions and axis limits in projected coordinates.

---

### 7. Background Interpolation (2D)
```bash
--idw-radius N
--idw-power N
```
Controls the inverse distance weighting parameters for smooth background mapping.

---

## Notes
- Default view is equivalent to `--mode field --display arrows --view 3d`.
- `--plane` requires `--view 2d` to take effect.
- Atom indices are 1-based and must be valid within the structure size.

---

## Quick Examples

### Default 3D Field
```bash
python -m charge_flow --xyz model.xyz --bader ACF.dat
```

### With Total Resultant
```bash
python -m charge_flow --xyz model.xyz --bader ACF.dat --show-total-resultant
```

### Atom Mode
```bash
python -m charge_flow --xyz model.xyz --bader ACF.dat --mode atom --atoms 5 12 33
```

### 2D Projection
```bash
python -m charge_flow --xyz model.xyz --bader ACF.dat --view 2d --plane xy
```

### Miller Plane Projection
```bash
python -m charge_flow --xyz model.xyz --bader ACF.dat --view 2d --plane-normal 1 1 0 --slab 1
```
