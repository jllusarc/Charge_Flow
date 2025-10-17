# charge_flow --- User Guide & Tutorial

------------------------------------------------------------------------

## 1) What this tool does (high level)

This CLI takes an **XYZ structure** and **Bader charge file** and computes **bond‑resolved charge flows** that are then aggregated into per‑atom **resultant vectors**. You can:
- Visualize a **2D projected field** (arrows, optional streamlines, divergence/magnitude maps).
- (Stub) Request **3D** spheres/arrows in PyVista (the current implementation is placeholder-only).
- Compute and unwrap **surface flux maps** onto a 2D image.
- Select/filter atoms by environment, degree, element, slab, or explicit indices.
- Export a **modified XYZ** (with renamed symbols for selected atoms) and a **CP2K KIND**s file.

Core ingredients under the hood:
- **Bond graph** from RDKit covalent radii (with optional pair cutoffs).
- Global linear system **A x = b** assigning bond flows so that node conservation matches Bader charge differences.
- Per‑atom local flow vectors → per‑atom resultant vector; projection to chosen plane for plotting.

------------------------------------------------------------------------

## 2) Installation

Conda is recommended:

```bash
conda install -c conda-forge numpy networkx scipy matplotlib pyvista rdkit
```

> Python dependencies visible in the package code: `numpy`, `networkx`, `scipy`, `matplotlib`, `pyvista`, `rdkit`.

------------------------------------------------------------------------

## 3) Quick start

```bash
# From the folder that contains the package folder `charge_transfer/`
python -m charge_transfer --xyz structure.xyz --bader ACF.dat --view 2d --contour mag --stream --save-2d field.png
```

Typical defaults:
- View is **2D** by default (`--view 2d`).
- Display mode for 3D would be `arrows` unless you set `--display spheres` (but 3D is currently a stub).

------------------------------------------------------------------------

## 4) Required inputs

### 4.1 XYZ structure

- Standard `.xyz` with optional two-line header. The loader tolerates comment headers; if it can't parse the first line as an integer, it just reads all atom lines as `Symbol x y z`.

### 4.2 Bader charges (`--bader`)

- The reader expects lines like:

  ```
  <index> ... ... ... <charge>
  ```

  where `<index>` is **1‑based** in the file and gets converted to 0‑based internally. The 5th token on each line is taken as the **charge**.

> Pro tip: Ensure the **atom count** and **ordering** match between `structure.xyz` and `ACF.dat`.

------------------------------------------------------------------------

## 5) How bonds and the linear system are built

### 5.1 Bond guessing

- Bonds are guessed with RDKit covalent radii; a pair **i--j** is bonded if the distance is below `scale * (r_i + r_j)`.
- You can **override** specific pairs with `--pair-cutoff A-B:dist` (Å). You may pass this flag multiple times.

### 5.2 Global system A x = b

- Unknowns **x**: one variable per **undirected bond**.
- Rows are atoms; conservation equations enforce that the sum of signed bond flows equals each atom's **Δq** from Bader.
- Sign convention in the construction is consistent across all atoms so that solving the system yields a self‑consistent set of **bond flows**.

### 5.3 Per‑atom field

- For each atom, local bond flows are converted into **vectors** (from neighbor direction and flow magnitude) and summed into a **resultant**.
- The collection of resultants forms the **global vector field** used by the 2D plotter.

------------------------------------------------------------------------

## 6) Atom selection & filtering

You can reduce the set of atoms used when building/visualizing the field:

- `--exclude-elements EL1 EL2 ...` → skip these elements entirely.
- `--exclude-coordination d1 d2 ...` → skip atoms **whose degree** (in the computed graph) is in this list.
- `--target-environments` strings like `In:As:3,Se:1` → select atoms whose **central** element and **neighbor element counts** match.
- `--atoms i j k ...` → explicitly select 1‑based indices (useful with `--mode atom`).

> When a *target filter* is provided, only those atoms are considered when forming the per‑atom field.

------------------------------------------------------------------------

## 7) Modes

- `--mode field` (default): build the full per‑atom vector field and plot or export.
- `--mode atom` + `--atoms ...`: compute local flows and resultants **only** for the specified atom indices (1‑based).
- `--mode surface-flux`: build a surface (e.g., convex hull) and accumulate the **normal component of flow** onto the surface, then unwrap to a 2D map.
- `--mode interface-flux`: reserved for future multi‑region interfaces.

------------------------------------------------------------------------

## 8) 2D visualization (matplotlib)

**General projection controls**  
- `--plane {xy,yz,xz,pca}` (default `xy`). `pca` uses the top two PCA axes of the kept points.  
- Or specify an arbitrary plane with `--plane-normal nx ny nz` and optional `--plane-origin ox oy oz`.  
- **Slab selection**: limit atoms by distance to the plane normal with `--slab HALF_THICK` (Å), or `--slab-q q` (quantile 0..1).

**Figure & map**  
- `--figsize W H` (inches), `--fig-lims xmin xmax ymin ymax`.  
- Background map: magnitude (`--contour mag`), divergence (`--contour div`), or disabled (`--contour none`).  
- Map styling: `--contour-levels N`, `--contour-alpha a`, `--map-alpha a`, `--map-cmap <cmap>`.  
- Grid interpolation: `--grid-res NX NY`, with IDW settings `--idw-radius R`, `--idw-power p`.

**Arrows & roles**  
- Arrows are **uniform color** by default (`--arrow-color`).  
- Automatically detects **cation/anion/neutral** roles from charge.  
- Override roles: `--ion-role In:cation --ion-role As:anion`.  
- Geometry: `--arrow-length-gain g`, `--arrow-offset o`, `--dot-radius r`, `--arrow-headlength Lh`, `--arrow-headaxislength Lah`.  
- Hide atom dots: `--no-scatter`.

**Streams**  
- `--stream` overlays streamlines with `--stream-density d` and `--stream-arrowsize s`.

**Saving**  
- `--save-2d path.png` saves the figure (interactive otherwise).

------------------------------------------------------------------------

## 9) 3D visualization (PyVista)

The entry point `viz_3d.plot_global_vector_field_pyvista(...)` is presently a **stub**. The CLI accepts:

```bash
--view 3d or --view both
--display {arrows,spheres}
--scale-3d, --sphere-scale, --show-total-resultant, --show-snowball
```

…but these have **no effect** until 3D plotting is implemented.

------------------------------------------------------------------------

## 10) Surface flux & unwrapping

- Build a surface from points (convex hull or alpha shape). Workflow:
  1. Build mesh from positions;
  2. Accumulate the **normal component** of the global vector field;
  3. **Unwrap** to a 2D map and save as an image.

Example:

```bash
python -m charge_transfer --xyz structure.xyz --bader ACF.dat --mode surface-flux --unwrap-save flux.png
```

------------------------------------------------------------------------

## 11) Exports & tagging

- **Modified XYZ**:

  ```bash
  --export-modified-xyz sel.xyz --rename-elements In:In1 As:As1
  ```

  - Only **selected atoms** are renamed.
  - Non‑selected atoms keep original symbols.
  - Header atom count remains unchanged.

- **KIND file** (CP2K convenience):

  ```bash
  --build-kinds KINDS.inp
  ```

------------------------------------------------------------------------

## 12) Full CLI reference (flags)

All command-line options are listed in detail in the repository documentation.

------------------------------------------------------------------------

## 13) Example runs

**A. PCA plane projection**
```bash
python -m charge_transfer --xyz structure.xyz --bader ACF.dat --view 2d --plane pca --stream --contour div --save-2d field_div.png
```

**B. Surface atoms only**
```bash
python -m charge_transfer --xyz structure.xyz --bader ACF.dat --exclude-coordination 4 5 --save-2d field_surface.png
```

**C. Environment selection + export**
```bash
python -m charge_transfer --xyz structure.xyz --bader ACF.dat --target-environments In:As:3,Se:1 --export-modified-xyz motif.xyz --rename-elements In:In1 As:As1
```

**D. Atom-specific mode**
```bash
python -m charge_transfer --xyz structure.xyz --bader ACF.dat --mode atom --atoms 12 57 81 --save-2d atoms.png
```

------------------------------------------------------------------------

## 14) Troubleshooting

- **Empty figure** → Too strict filters; loosen them.  
- **Weird bonding** → Adjust `--pair-cutoff`.  
- **Bader mismatch** → Ensure 1-based indexing and identical atom ordering.  
- **3D not appearing** → Stub; use 2D.  
- **Noisy divergence map** → Increase grid resolution or radius.

------------------------------------------------------------------------

## 15) Developer notes

- `viz_3d.py`: placeholder for 3D plotting (PyVista).  
- `surfaces.py`: convex hull, alpha shape, unwrap utilities.  
- `exports.py`: preserves header atom count; renames selected atoms only.  
- `core.py`: core I/O, RDKit bonding, slab filter, PCA plane, linear solver.

------------------------------------------------------------------------

## 16) Minimal reproducible example

```bash
python -m charge_transfer --xyz structure.xyz --bader ACF.dat --plane xy --contour mag --save-2d out.png
```

