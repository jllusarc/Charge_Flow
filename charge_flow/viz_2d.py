import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap, ScalarMappable
from typing import Optional, Dict, List

from .utils import pca_basis, fixed_basis, infer_role_from_entry


def _project_points_vectors(positions, global_field, plane='xy', plane_normal=None, plane_origin=None):
    pos = np.asarray(positions, float)
    kept = [e['atom_index'] for e in global_field]
    if not kept:
        return np.zeros((0, 2)), np.zeros((0, 2)), np.zeros((0,)), [], None

    pts = pos[kept]

    if plane_normal is not None:
        n = np.asarray(plane_normal, float)
        C = np.asarray(plane_origin if plane_origin is not None else pts.mean(axis=0), float)
        ref = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        e1 = np.cross(n, ref)
        e1 = e1 / (np.linalg.norm(e1) + 1e-12)
        e2 = np.cross(n, e1)
        e2 = e2 / (np.linalg.norm(e2) + 1e-12)
    else:
        if plane == 'pca':
            C, e1, e2, n = pca_basis(pts)
        else:
            _, e1, e2, n = fixed_basis(plane)
            C = pts.mean(axis=0)

    R = pts - C
    x2 = R @ e1
    y2 = R @ e2
    P2 = np.c_[x2, y2]

    V = np.array([e['vector'] for e in global_field], float)
    V_in_plane = V - np.outer(V @ n, n)
    Vx = V_in_plane @ e1
    Vy = V_in_plane @ e2
    V2 = np.c_[Vx, Vy]

    mags = np.asarray([e['magnitude'] for e in global_field], float)
    return P2, V2, mags, kept, (C, e1, e2, n)


def _idw(P2, V2, U, V, radius=None, power=2.0):
    x = P2[:, 0]; y = P2[:, 1]
    vx = V2[:, 0]; vy = V2[:, 1]
    Xf = U.ravel(); Yf = V.ravel()
    VXf = np.zeros_like(Xf); VYf = np.zeros_like(Yf)
    r2_max = None if radius is None else float(radius) ** 2

    for k, (xx, yy) in enumerate(zip(Xf, Yf)):
        dx = xx - x; dy = yy - y; d2 = dx * dx + dy * dy
        zero = np.where(d2 == 0.0)[0]
        if zero.size:
            i0 = zero[0]; VXf[k] = vx[i0]; VYf[k] = vy[i0]; continue

        if r2_max is not None:
            m = d2 <= r2_max
            if not np.any(m):
                continue
            d2 = d2[m]; vx_sel = vx[m]; vy_sel = vy[m]
        else:
            vx_sel = vx; vy_sel = vy

        w = 1.0 / (np.power(d2, power / 2.0) + 1e-12)
        wsum = np.sum(w)
        VXf[k] = np.sum(w * vx_sel) / wsum
        VYf[k] = np.sum(w * vy_sel) / wsum

    return VXf.reshape(U.shape), VYf.reshape(V.shape)


def _finite_div(U, V, VX, VY):
    du = np.mean(np.diff(U[0, :])) if U.shape[1] > 1 else 1.0
    dv = np.mean(np.diff(V[:, 0])) if V.shape[0] > 1 else 1.0
    dVx_du = np.zeros_like(VX); dVy_dv = np.zeros_like(VY)
    dVx_du[:, 1:-1] = (VX[:, 2:] - VX[:, :-2]) / (2.0 * du)
    dVy_dv[1:-1, :] = (VY[2:, :] - VY[:-2, :]) / (2.0 * dv)
    dVx_du[:, 0] = (VX[:, 1] - VX[:, 0]) / du
    dVx_du[:, -1] = (VX[:, -1] - VX[:, -2]) / du
    dVy_dv[0, :] = (VY[1, :] - VY[0, :]) / dv
    dVy_dv[-1, :] = (VY[-1, :] - VY[-2, :]) / dv
    return dVx_du + dVy_dv


def plot_global_vector_field_2d(
    positions,
    global_field,
    *,
    plane='xy',
    plane_normal=None,
    plane_origin=None,
    cmin=None,
    cmax=None,
    figsize=(8, 7),
    save_path=None,
    show=True,
    symbols=None,
    max_arrows=600,
    base_len=0.35,
    scatter=True,
    grid_res=(300, 300),
    idw_radius=None,
    idw_power=2.0,
    contour='none',
    contour_levels=15,
    contour_alpha=0.55,
    map_alpha=0.35,
    map_cmap='inferno',
    stream=False,
    stream_density=1.4,
    stream_arrowsize=1.0,
    fig_lims=None,
    arrow_offset=0.0,
    dot_radius=0.0,
    headlength=9,
    headaxislength=7,
):
    """Projected vector field; auto cation/anion pivots from entry['delta_q']/element, plus visual offset."""
    P2, V2, mags, kept, basis = _project_points_vectors(
        positions, global_field, plane=plane, plane_normal=plane_normal, plane_origin=plane_origin
    )
    fig, ax = plt.subplots(figsize=figsize)

    if len(kept) == 0:
        if show:
            plt.show()
        return fig

    vmin = cmin if cmin is not None else float(mags.min())
    vmax = cmax if cmax is not None else float(mags.max())

    x, y = P2[:, 0], P2[:, 1]
    if fig_lims is None:
        m = 0.05
        xmin, xmax = x.min(), x.max(); xr = xmax - xmin; xmin -= xr * m; xmax += xr * m
        ymin, ymax = y.min(), y.max(); yr = ymax - ymin; ymin -= yr * m; ymax += yr * m
    else:
        xmin, xmax, ymin, ymax = fig_lims

    uu = np.linspace(xmin, xmax, grid_res[0])
    vv = np.linspace(ymin, ymax, grid_res[1])
    U, V = np.meshgrid(uu, vv)

    if idw_radius is None:
        diag = np.hypot(xmax - xmin, ymax - ymin)
        idw_radius = 0.17 * diag

    VX, VY = _idw(P2, V2, U, V, radius=idw_radius, power=idw_power)
    M = np.hypot(VX, VY)

    cmap = get_cmap(map_cmap)
    norm = Normalize(vmin=vmin, vmax=vmax)

    # Background magnitude map
    im = ax.imshow(
        M,
        extent=[xmin, xmax, ymin, ymax],
        origin='lower',
        alpha=map_alpha,
        cmap=map_cmap,
        vmin=vmin,
        vmax=vmax,
        aspect='auto',
    )

    # Optional contours
    if contour == 'mag':
        CS = ax.contour(U, V, M, levels=contour_levels, alpha=contour_alpha)
        ax.clabel(CS, inline=True, fontsize=8)
    elif contour == 'div':
        DIV = _finite_div(U, V, VX, VY)
        CS = ax.contour(U, V, DIV, levels=contour_levels, alpha=contour_alpha)
        ax.clabel(CS, inline=True, fontsize=8)

    # Optional streamlines
    if stream:
        ax.streamplot(U, V, VX, VY, density=stream_density, arrowsize=stream_arrowsize)

    # Downsample for quiver
    n = len(P2)
    step = max(1, int(np.ceil(n / max_arrows)))
    idxs_all = np.arange(n)[::step]

    # role per atom
    roles = []
    for k in range(n):
        entry = global_field[k]
        sym = entry.get('element', symbols[kept[k]] if symbols is not None else None)
        roles.append(infer_role_from_entry(entry, sym))

    def quiver_subset(idxs, pivot):
        if idxs.size == 0:
            return None
        P = P2[idxs].copy()
        Vv = V2[idxs]
        mag = mags[idxs]
        Uv = Vv / (np.linalg.norm(Vv, axis=1, keepdims=True) + 1e-12)
        off = dot_radius if dot_radius > 0 else arrow_offset
        if off > 0:
            if pivot == 'tail':
                P = P + Uv * off
            elif pivot == 'tip':
                P = P - Uv * off
        return ax.quiver(
            P[:, 0], P[:, 1], Vv[:, 0], Vv[:, 1], mag,
            cmap=map_cmap, norm=norm,
            scale=1 / base_len, width=0.003,
            headlength=headlength, headaxislength=headaxislength,
            pivot=pivot
        )

    idx_cation  = idxs_all[[i for i in range(idxs_all.size) if roles[idxs_all[i]] == 'cation']]
    idx_anion   = idxs_all[[i for i in range(idxs_all.size) if roles[idxs_all[i]] == 'anion']]
    idx_neutral = idxs_all[[i for i in range(idxs_all.size) if roles[idxs_all[i]] not in ('anion', 'cation')]]

    Qc = quiver_subset(idx_cation, 'tail')
    Qa = quiver_subset(idx_anion,  'tip')
    Qn = quiver_subset(idx_neutral, 'mid')

    # Colorbar bound to this axes (use ScalarMappable when quiver colors by array)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])  # required by Matplotlib for ScalarMappable-only colorbars
    fig.colorbar(sm, ax=ax, label='|Flow| (rel.)')

    # atoms overlay (if requested)
    if scatter:
        ax.scatter(P2[:, 0], P2[:, 1], s=6, facecolors='white', edgecolors='k', linewidths=0.2, alpha=0.8)

    ax.set_xlim([xmin, xmax]); ax.set_ylim([ymin, ymax])
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x1 (angstroms)'); ax.set_ylabel('x2 (angstroms)')
    ax.set_title('Projected charge-flow field')

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300)
    if show:
        plt.show()
    return fig
