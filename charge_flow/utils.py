
from __future__ import annotations
import numpy as np
from typing import Tuple, Dict, List, Iterable, Optional

# Default periodic-chemistry heuristics for roles
DEFAULT_ANIONS  = {'F','Cl','Br','I','O','S','Se','Te','N','P','As','Sb','Bi'}
DEFAULT_CATIONS = {
    'Li','Na','K','Rb','Cs','Fr','Be','Mg','Ca','Sr','Ba','Ra',
    'Al','Ga','In','Tl','B','Sc','Y','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
    'Ti','Zr','Hf','V','Nb','Ta','Cr','Mo','W','Mn','Tc','Re','Fe','Ru','Os','Co','Rh','Ir','Ni','Pd','Pt','Cu','Ag','Au','Zn','Cd','Hg'
}

def infer_role_from_entry(entry: dict, symbol: str | None) -> str:
    """Infer ion role from entry/meta data; fallback to periodic defaults."""
    # Prefer explicit per-atom data if present
    for key in ('electron_gain','gains_electrons','is_acceptor'):
        if key in entry:
            return 'anion' if bool(entry[key]) else 'cation'
    for key in ('delta_e','dq','delta_q','net_charge','charge_diff'):
        if key in entry:
            try:
                v = float(entry[key])
                if v < 0: return 'anion'
                if v > 0: return 'cation'
            except Exception:
                pass
    # Fallback
    if symbol in DEFAULT_ANIONS:  return 'anion'
    if symbol in DEFAULT_CATIONS: return 'cation'
    return 'neutral'

def pca_basis(points: np.ndarray):
    X = np.asarray(points, float)
    C = X.mean(axis=0)
    _, _, Vt = np.linalg.svd(X - C, full_matrices=False)
    e1, e2, n = Vt[0], Vt[1], Vt[2]
    return C, e1/np.linalg.norm(e1), e2/np.linalg.norm(e2), n/np.linalg.norm(n)

def fixed_basis(plane: str):
    if plane == 'xy': return np.zeros(3), np.array([1.,0.,0.]), np.array([0.,1.,0.]), np.array([0.,0.,1.])
    if plane == 'yz': return np.zeros(3), np.array([0.,1.,0.]), np.array([0.,0.,1.]), np.array([1.,0.,0.])
    if plane == 'xz': return np.zeros(3), np.array([1.,0.,0.]), np.array([0.,0.,1.]), np.array([0.,1.,0.])
    raise ValueError("plane must be 'xy'|'yz'|'xz'")

def plane_from_args(plane: str, plane_normal, plane_origin, points: np.ndarray):
    pts = np.asarray(points, float)
    if plane_normal is not None:
        n = np.asarray(plane_normal, float)
        n = n / (np.linalg.norm(n) + 1e-12)
        helper = np.array([1.0,0.0,0.0]) if abs(n[0]) < 0.9 else np.array([0.0,1.0,0.0])
        e1 = helper - n * np.dot(helper, n); e1 /= (np.linalg.norm(e1)+1e-12)
        e2 = np.cross(n, e1); e2 /= (np.linalg.norm(e2)+1e-12)
        C  = np.asarray(plane_origin, float) if plane_origin is not None else pts.mean(axis=0)
        return C, e1, e2, n
    if plane == 'pca': return pca_basis(pts)
    _, e1, e2, n = fixed_basis(plane)
    C = pts.mean(axis=0)
    return C, e1, e2, n

def filter_by_slab(positions: np.ndarray, C: np.ndarray, n: np.ndarray,
                   half_thickness: float | None, quantile: float | None):
    """Return boolean mask of atoms inside the slab |(r-C)·n| <= h or <= q-quantile."""
    r = np.asarray(positions, float)
    d = np.abs((r - C) @ n)
    if half_thickness is not None:
        return d <= float(half_thickness)
    if quantile is not None:
        thr = np.quantile(d, float(quantile))
        return d <= thr
    return np.ones(len(r), dtype=bool)

def parse_pair_cutoffs(items: list[str] | None):
    pair_cut = {}
    for s in (items or []):
        try:
            ab, d = s.split(':'); a, b = ab.split('-'); pair_cut[tuple(sorted((a, b)))] = float(d)
        except Exception:
            print(f"[warn] Ignoring malformed --pair-cutoff '{s}' (use A-B:dist)")
    return pair_cut

def parse_kv_pairs(items: list[str] | None):
    out = {}
    for it in (items or []):
        if ':' in it:
            k,v = it.split(':',1); out[k]=v
    return out

def parse_target_environment(env_str: str):
    """Format: Central:El1:n1,El2:n2  -> returns ('Central', {'El1':n1, ...})"""
    central_part, neighbors_part = env_str.split(':', 1)
    neighbor_dict = {}
    for item in neighbors_part.split(','):
        el, n = item.split(':')
        neighbor_dict[el] = int(n)
    return central_part, neighbor_dict