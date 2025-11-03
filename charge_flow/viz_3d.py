import numpy as np
import pyvista as pv
from .core import (
    surface_atoms_by_element_degree,
    shell_depths_from_surface,
    compute_snowball_cumulative,
)

def add_structure_wireframe(
    plotter: pv.Plotter,
    positions: np.ndarray,
    bonds,
    line_radius: float = 0.02,
    opacity: float = 0.25,
):
    """Simple cylinder wireframe for bonds."""
    for i, j in bonds:
        start, end = positions[i], positions[j]
        vec = end - start
        cyl = pv.Cylinder(
            center=(start + end) / 2,
            direction=vec,
            radius=line_radius,
            height=np.linalg.norm(vec),
            resolution=12,
        )
        plotter.add_mesh(cyl, color="gray", opacity=opacity)


def _total_resultant_arrow(p: pv.Plotter, positions, global_field, scale=2.0, color="cyan"):
    """Draw a single global resultant arrow at the structure’s center of mass."""
    if not global_field:
        return
    R = np.sum([e["vector"] for e in global_field], axis=0)
    C = positions.mean(axis=0)
    arr = pv.Arrow(
        start=C,
        direction=R * scale,
        tip_length=0.2 * scale,
        tip_radius=0.08 * scale,
        shaft_radius=0.04 * scale,
    )
    p.add_mesh(arr, color=color)


def _snowball_slider(
    p: pv.Plotter,
    G,
    symbols,
    positions,
    global_field,
    scale=2.0,
    color="yellow",
):
    """Interactive cumulative (snowball) resultant via a shell-depth slider."""
    surface = surface_atoms_by_element_degree(G, symbols)
    depth = shell_depths_from_surface(G, surface)
    entries = compute_snowball_cumulative(global_field, depth)
    if not entries:
        return

    shell_vals = [e["shell"] for e in entries]
    vecs = [e["vec"] for e in entries]

    def cb(val):
        k = int(round(val))
        idx = shell_vals.index(k) if k in shell_vals else 0
        v = vecs[idx]
        try:  # remove previous snowball arrow (if any)
            p.remove_actor("snowball")
        except Exception:
            pass
        arr = pv.Arrow(
            start=positions.mean(axis=0),
            direction=v * scale,
            tip_length=0.2 * scale,
            tip_radius=0.08 * scale,
            shaft_radius=0.04 * scale,
        )
        p.add_mesh(arr, color=color, name="snowball")
        p.render()

    p.add_slider_widget(
        cb,
        rng=[min(shell_vals), max(shell_vals)],
        value=min(shell_vals),
        title="Snowball shell",
    )


def plot_global_vector_field_pyvista(
    positions,
    global_field,
    bonds,
    *,
    scale_factor: float = 2.0,
    colormap: str = "inferno",
    cmin=None,
    cmax=None,
    display_mode: str = "arrows",
    show: bool = True,
    symbols=None,
    sphere_scale: float = 0.2,
    sphere_base_radius: float = 5.0,
    show_total: bool = False,
    total_scale: float = 2.0,
    total_color: str = "cyan",
    G=None,
    show_snowball: bool = False,
    snowball_scale: float = 2.0,
    snowball_color: str = "yellow",
    plane=None,
    plane_normal=None,
    plane_origin=None,
    slab=None,
    slab_q=None,
    plane_alpha: float = 0.25,
    plane_color: str = "cyan",
):
    """Main 3D plot: spheres or arrows, optional resultant and snowball, and robust picking."""
    p = pv.Plotter()
    add_structure_wireframe(p, positions, bonds)

    kept = [e["atom_index"] for e in global_field]
    mags = np.array([e["magnitude"] for e in global_field])

    if mags.size == 0:
        if show:
            p.show()
        return p

    vmin = cmin if cmin is not None else float(mags.min())
    vmax = cmax if cmax is not None else float(mags.max())

    if display_mode == "spheres":
        cloud = pv.PolyData(positions[kept])
        cloud["magnitudes"] = mags
        sphere = pv.Sphere(radius=sphere_base_radius)
        glyphs = cloud.glyph(
            scale="magnitudes",
            geom=sphere,
            orient=False,
            factor=float(sphere_scale),
        )
        p.add_mesh(glyphs, scalars="magnitudes", cmap=colormap, clim=[vmin, vmax])
    else:
        mb = pv.MultiBlock()
        for e in global_field:
            idx = e["atom_index"]
            v = e["vector"] * scale_factor
            m = e["magnitude"]
            if m == 0:
                continue
            thick = 0.02 * scale_factor + 0.13 * scale_factor * (m / vmax if vmax > 0 else 0)
            arr = pv.Arrow(
                start=positions[idx],
                direction=v,
                tip_length=0.2 * scale_factor,
                tip_radius=1.5 * thick,
                shaft_radius=thick,
            )
            arr.point_data["Flow Magnitude"] = np.full(arr.n_points, m)
            mb.append(arr)
        arrows = mb.combine()
        p.add_mesh(arrows, scalars="Flow Magnitude", cmap=colormap, clim=[vmin, vmax])

    # ----------------------------- Robust point picking -----------------------------
    kept_pos = positions[kept]

    def _on_pick(*args):
        """
        PyVista may pass (picker) or (picker, event). This callback:
        1) uses plotter.picked_point if available
        2) otherwise tries picker.GetPickPosition()
        3) falls back to plotter.pick_mouse_position()
        """
        pt = None

        # Preferred: the last picked 3D point stored by PyVista
        picked = getattr(p, "picked_point", None)
        if picked is not None:
            try:
                pt = np.array(picked, dtype=float)
            except Exception:
                pt = None

        # If a picker object is passed, use it
        if pt is None and len(args) >= 1:
            picker = args[0]
            getpos = getattr(picker, "GetPickPosition", None)
            if callable(getpos):
                try:
                    pt = np.array(getpos(), dtype=float)
                except Exception:
                    pt = None

        # Fallback: ray-cast from mouse position
        if pt is None:
            try:
                pt = np.array(p.pick_mouse_position(), dtype=float)
            except Exception:
                return  # nothing we can do

        # Snap to nearest kept atom and report
        d = np.linalg.norm(kept_pos - pt[None, :], axis=1)
        k = int(np.argmin(d))
        atom_idx = kept[k]
        mag = mags[k]
        sym = symbols[atom_idx] if symbols is not None else "Atom"
        print(f"{sym} ({atom_idx+1}) : {mag:.4f} e⁻")

        # Visual highlight for the picked atom (replace previous, if any)
        try:
            p.remove_actor("picked_atom")
        except Exception:
            pass
        highlight = pv.Sphere(radius=0.4, center=positions[atom_idx])
        p.add_mesh(
            highlight,
            color="lime",
            opacity=0.6,
            name="picked_atom",
            reset_camera=False,
        )
        p.render()

    # New API: prefer picker-based pipeline; avoid deprecated `use_mesh`
    p.enable_point_picking(
        callback=_on_pick,
        left_clicking=False,   # set False if you prefer right click
        show_message=False,
        use_picker=True,      # <- replaces deprecated `use_mesh`
    )
    # -------------------------------------------------------------------------------

    if show_total:
        _total_resultant_arrow(p, positions, global_field, scale=total_scale, color=total_color)

    if show_snowball and G is not None:
        _snowball_slider(
            p,
            G,
            symbols,
            positions,
            global_field,
            scale=snowball_scale,
            color=snowball_color,
        )

    p.show_axes()
    if show:
        p.show()
    return p


def plot_atom_mode(p, positions, bonds, symbols, flows_dict, scale=2.0):
    """Per-atom visualization: local flows around selected centers."""
    p = p or pv.Plotter()
    add_structure_wireframe(p, positions, bonds)

    for center, flows in flows_dict.items():
        for fl in flows:
            j = fl["to"] - 1
            q = fl["charge_flow"]
            start = positions[center]
            direction = positions[center] - positions[j]
            direction = direction / (np.linalg.norm(direction) + 1e-12)
            v = (abs(q) * direction) * (scale if q < 0 else -scale)
            arr = pv.Arrow(
                start=start,
                direction=v,
                tip_length=0.2 * scale,
                tip_radius=0.08 * scale,
                shaft_radius=0.04 * scale,
            )
            color = "cyan" if q < 0 else "magenta"
            p.add_mesh(arr, color=color)

        s = pv.Sphere(radius=0.2 * scale, center=positions[center])
        p.add_mesh(s, color="yellow")

    p.show_axes()
    p.show()
    return p
