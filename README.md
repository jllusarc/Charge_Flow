# ⚡ charge_flow — Easy Start Guide

*A fast-start companion to the full [User Manual](charge_flow_tutorial_manual.md)*

---

## 🧬 What it does

`charge_transfer` reads an **XYZ atomic structure** and a **Bader charge file (`ACF.dat`)** to compute **bond-resolved charge-flow vectors** between atoms.
It can visualize **2D projected fields**, (stub) **3D arrows/spheres**, and export **labeled XYZ** or **CP2K KIND files**.

---

## ⚙️ Installation

Install dependencies via conda:

```bash
conda install -c conda-forge numpy networkx scipy matplotlib pyvista rdkit
```

Then clone this repository and run from its root folder.

---

## 🚀 Minimal run

```bash
python -m charge_transfer --xyz structure.xyz --bader ACF.dat
```

This opens the **3D viewer** (stub implementation).
By default, arrows indicate **charge-flow magnitude and direction** — longer arrows mean stronger flow.

---

## 🖥️ Quick examples

| Goal                                                        | Command                                                         |
| ----------------------------------------------------------- | --------------------------------------------------------------- |
| Show 3D vector field (default)                              | `python -m charge_transfer --xyz structure.xyz --bader ACF.dat` |
| Adjust arrow size                                           | `--scale-3d 1.5`                                                |
| Show total resultant vector                                 | `--show-total-resultant`                                        |
| Color and scale resultant                                   | `--resultant-color red --resultant-scale 2`                     |
| Enable “snowball” cumulative resultant (interactive slider) | `--show-snowball`                                               |
| Adjust snowball scale or color                              | `--snowball-scale 1.5 --snowball-color blue`                    |
| Display spheres instead of arrows                           | `--display spheres`                                             |
| Adjust sphere size                                          | `--sphere-scale 0.5`                                            |
| Show only specific atoms                                    | `--mode atom --atoms 5 10 22`                                   |
| 2D projected field (recommended view)                       | `--view 2d`                                                     |
| Choose projection plane                                     | `--view 2d --plane yz`                                          |
| Set plane by normal vector                                  | `--view 2d --plane-normal 1 0 0`                                |
| Define plane thickness                                      | `--view 2d --slab 1`                                            |
| Adjust figure size                                          | `--view 2d --figsize 8 7`                                       |
| Tune interpolation radius/intensity                         | `--idw-radius 3 --idw-power 2`                                  |

> ⚠️ `--mode surface-flux` and `--mode interface-flux` are placeholders and currently **not functional**.

---

## 📊 2D Output Example

```bash
python -m charge_transfer --xyz structure.xyz --bader ACF.dat \
  --view 2d --contour mag --stream --save-2d field.png
```

→ Produces a vector-field map (`field.png`) projected on the **xy-plane**.

---

## 🖾 Export options

```bash
# Export selected atoms as labeled XYZ
--export-modified-xyz sel.xyz --rename-elements In:In1 As:As1

# Create a CP2K KINDS input
--build-kinds KINDS.inp
```

---

## 🧠 Tips & Troubleshooting

* Empty figure? → Your filters excluded all atoms. Loosen `--slab` or `--target` options.
* Mismatched atoms? → Ensure `ACF.dat` and `.xyz` have **identical atom order**.
* 3D missing? → Still a stub — use `--view 2d` instead.
* Noisy divergence map? → Increase `--grid-res` or `--idw-radius`.

---

## 📚 Further Reading

For full flag documentation, theory, and advanced workflows, see:
👉 **[`charge_flow_tutorial_manual.md`](charge_flow_tutorial_manual.md)**
