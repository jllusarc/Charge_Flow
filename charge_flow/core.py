
from __future__ import annotations
import numpy as np
import networkx as nx
from rdkit import Chem
from collections import Counter, defaultdict, deque
from typing import List, Dict, Tuple, Iterable, Optional

from .utils import plane_from_args, filter_by_slab, parse_target_environment

def load_xyz(filepath: str):
    symbols, positions = [], []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    # handle XYZ with header (first 2 lines) or plain
    start = 2 if len(lines) >= 2 and lines[0].strip().isdigit() else 0
    for line in lines[start:]:
        parts = line.strip().split()
        if len(parts) >= 4:
            symbols.append(parts[0])
            positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.array(positions, float)

def load_bader(filepath: str):
    """Parse ACF.dat: index, x, y, z, charge ... (we use column 5)."""
    bader = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 5 and parts[0].isdigit():
                idx = int(parts[0]) - 1
                charge = float(parts[4])
                bader[idx] = charge
    return bader

def guess_bonds_rdkit(symbols, positions, scale=1.2, pair_cutoffs=None):
    """Guess bonds by covalent radii sum (Cordero) with optional per-pair overrides A-B:dist."""
    custom_radii = {
        'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
        'P': 1.07, 'S': 1.05, 'Cl': 1.02,
        'As': 1.19, 'Se': 1.20, 'Br': 1.20,
        'Ga': 1.26, 'In': 1.42, 'Zn': 1.22, 'Cd': 1.44,
        'Sb': 1.39, 'Te': 1.38, 'Cs': 2.35, 'Pb': 1.46, 'N': 0.71, 'Hf': 1.58
    }
    pt = Chem.GetPeriodicTable()
    bonds = set()
    N = len(symbols)
    for i in range(N):
        for j in range(i+1, N):
            pair = tuple(sorted((symbols[i], symbols[j])))
            d = np.linalg.norm(positions[i] - positions[j])
            if pair_cutoffs and pair in pair_cutoffs:
                if d < pair_cutoffs[pair]: bonds.add((i, j))
                continue
            r1 = custom_radii.get(symbols[i], pt.GetRcovalent(symbols[i]))
            r2 = custom_radii.get(symbols[j], pt.GetRcovalent(symbols[j]))
            if d < scale * (r1 + r2):
                bonds.add((i, j))
    return list(bonds)

def build_graph_and_data(symbols, positions, bonds, bader):
    """Build NetworkX graph and per-atom data including delta_q if valence is known."""
    valence_Z = {'In': 13, 'Ga': 13, 'Zn': 12, 'Cd': 12, 'Hf': 12, 'O': 6, 'S': 6, 'Se': 6, 'Te': 6,
                 'P': 5, 'As': 5, 'Sb': 5, 'Cl': 7, 'Br': 7, 'Cs': 9, 'Pb': 4, 'N': 5, 'H':1, 'C': 4}
    G = nx.Graph()
    atom_data = {}
    for i, (el, pos) in enumerate(zip(symbols, positions)):
        G.add_node(i, element=el, position=pos)
        if i in bader:
            if el in valence_Z:
                Z = valence_Z[el]
                dq = Z - bader[i]
            else:
                dq = 0.0
            atom_data[i] = {'element': el, 'bader': bader[i], 'delta_q': float(dq)}
    for i, j in bonds:
        G.add_edge(i, j)
    return G, atom_data

def build_global_system(atom_data, G, bond_pairs):
    """Construct A x = b where x are signed flows on bonds, rows are atom conservation eqs."""
    valid_atoms = list(atom_data.keys())
    valid_bonds = [tuple(sorted((i, j))) for (i, j) in bond_pairs if i in atom_data and j in atom_data]
    bond_index = {bond: k for k, bond in enumerate(valid_bonds)}
    atom_to_row = {a: r for r, a in enumerate(valid_atoms)}
    A = np.zeros((len(valid_atoms), len(valid_bonds)), float)
    b = np.zeros(len(valid_atoms), float)
    for atom, row in atom_to_row.items():
        b[row] = atom_data[atom]['delta_q']
        for nb in G.neighbors(atom):
            if nb in atom_to_row:
                col = bond_index[tuple(sorted((atom, nb)))]
                sign = 1.0 if atom < nb else -1.0
                A[row, col] += sign
    return valid_atoms, valid_bonds, A, b

def extract_local_flows(atom_index, valid_bonds, bond_charges, G):
    flows = []
    for (i, j), idx in zip(valid_bonds, range(len(valid_bonds))):
        if atom_index in (i, j):
            other = j if i == atom_index else i
            sign = 1 if i == atom_index else -1
            qij = sign * bond_charges[idx]
            flows.append({
                "from": atom_index + 1, "to": other + 1,
                "element": G.nodes[other]['element'],
                "source_element": G.nodes[atom_index]['element'],
                "charge_flow": float(qij)
            })
    return flows

def compute_charge_flow_vectors(atom_index, flows, positions):
    """Arrow sign: points toward acceptor (electron gain)."""
    vectors = []
    for fl in flows:
        src = fl['from'] - 1; tgt = fl['to'] - 1
        q = fl['charge_flow']
        mag = abs(q)
        nb = tgt if src == atom_index else src
        direction = positions[atom_index] - positions[nb]
        nrm = np.linalg.norm(direction)
        unit = direction / nrm if nrm != 0 else np.zeros(3)
        arrow_vec = (mag * unit) if q < 0 else (-mag * unit)
        vectors.append({
            'start': positions[atom_index], 'vector': arrow_vec, 'magnitude': mag,
            'neighbor_index': nb, 'charge_flow': q
        })
    return vectors

def compute_resultant_flow_vector(flow_vectors):
    if not flow_vectors: return np.zeros(3), 0.0
    tot = np.sum([v['vector'] for v in flow_vectors], axis=0)
    return tot, float(np.linalg.norm(tot))

def compute_global_vector_field(valid_atoms, valid_bonds, x, G, positions, symbols,
                                exclude_elements=None, exclude_coord=None, target_filter=None, atom_data=None):
    """Compute resultant vector & magnitude per atom; include element & delta_q for viz logic."""
    exclude_elements = set(exclude_elements or [])
    exclude_coord = set(exclude_coord or [])
    target_set = set(target_filter or [])

    field = []
    for a in valid_atoms:
        if symbols[a] in exclude_elements: continue
        if exclude_coord and (G.degree[a] in exclude_coord): continue
        if target_set and (a not in target_set): continue
        flows = extract_local_flows(a, valid_bonds, x, G)
        vecs = compute_charge_flow_vectors(a, flows, positions)
        R, M = compute_resultant_flow_vector(vecs)
        dq = 0.0
        if atom_data and a in atom_data: dq = float(atom_data[a].get('delta_q', 0.0))
        field.append({'atom_index': a, 'vector': R, 'magnitude': M,
                      'element': symbols[a], 'delta_q': dq})
    return field

def match_environment(G, symbols, target_central, target_neighbors):
    from collections import Counter
    matched_atoms = []
    for node in G.nodes:
        if symbols[node] != target_central: continue
        neighbors = [symbols[n] for n in G.neighbors(node)]
        if Counter(neighbors) == Counter(target_neighbors):
            matched_atoms.append(node); matched_atoms.extend(G.neighbors(node))
    return list(set(matched_atoms))

def select_indices_by_environments(G, symbols, env_strings):
    selected = []
    if not env_strings: return selected
    for env_str in env_strings:
        try:
            central, neighbor_dict = parse_target_environment(env_str)
            neighbor_list = sum([[k] * v for k, v in neighbor_dict.items()], [])
            selected.extend(match_environment(G, symbols, central, neighbor_list))
        except Exception as e:
            print(f"[warn] Bad --target-environments item '{env_str}': {e}")
    return sorted(set(selected))

def plane_from_args_public(plane, plane_normal, plane_origin, positions):
    return plane_from_args(plane, plane_normal, plane_origin, positions)

def filter_by_slab_public(positions, C, n, half_thickness=None, quantile=None):
    return filter_by_slab(positions, C, n, half_thickness, quantile)

def compute_total_resultant_from_atoms(global_field):
    if not global_field: return np.zeros(3), 0.0
    R = np.sum([e['vector'] for e in global_field], axis=0)
    return R, float(np.linalg.norm(R))

# Surface depth helpers for snowball
def surface_atoms_by_element_degree(G, symbols):
    from collections import defaultdict
    max_deg = defaultdict(int)
    for i in G.nodes:
        el = symbols[i]; max_deg[el] = max(max_deg[el], G.degree[i])
    return [i for i in G.nodes if G.degree[i] < max_deg[symbols[i]]]

def shell_depths_from_surface(G, surface_nodes):
    from collections import deque
    depth = {i: np.inf for i in G.nodes}
    q = deque()
    for s in surface_nodes:
        depth[s] = 0; q.append(s)
    while q:
        u = q.popleft()
        for v in G.neighbors(u):
            if depth[v] == np.inf:
                depth[v] = depth[u] + 1; q.append(v)
    return depth

def compute_snowball_cumulative(global_field, depth_map):
    if not global_field: return []
    all_depths = [int(depth_map.get(e['atom_index'], 0)) for e in global_field]
    max_shell = int(max(all_depths)) if all_depths else 0
    entries = []
    for k in range(max_shell, -1, -1):
        vec = np.sum([e['vector'] for e in global_field if int(depth_map.get(e['atom_index'], 0)) >= k], axis=0)
        entries.append({'shell': k, 'vec': vec, 'mag': float(np.linalg.norm(vec))})
    return entries
