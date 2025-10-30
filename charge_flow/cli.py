from typing import Optional, Sequence
"""
CLI for charge-flow analysis.

- Compute per-atom resultant charge-flow vectors
- Visualize in 2D/3D
- Export tagged XYZ (coordination / environment / magnitude)
- Emit CP2K &KIND blocks for labels
"""


from __future__ import annotations
import argparse, os
import numpy as np

from .core import (
    load_xyz, load_bader, guess_bonds_rdkit, build_graph_and_data,
    build_global_system, compute_global_vector_field, extract_local_flows,
    select_indices_by_environments, plane_from_args_public, filter_by_slab_public
)
from .exports import (
    export_xyz_with_labels, tag_by_coordination, tag_by_environment, tag_by_magnitude,
    write_kind_sections_for_xyz, kind_sections_from_labels
)
from .utils import parse_pair_cutoffs
from . import viz_2d, viz_3d

def parse_args(argv: Optional[Sequence[str]] = None):
    p = argparse.ArgumentParser(description='Charge flow analysis')
    # Required IO
    p.add_argument('--xyz', required=True, help='XYZ structure file')
    p.add_argument('--bader', required=True, help='Bader ACF.dat file')
    # Selection and graph
    p.add_argument('--pair-cutoff', action='append', default=[], help="Override bond cutoff A-B:dist (Å). Repeatable.")
    p.add_argument('--exclude-elements', nargs='*', default=[], help="Elements to exclude from the field.")
    p.add_argument('--exclude-coordination', nargs='*', type=int, default=[], help="Exclude atoms with these coordinations.")
    p.add_argument('--target-environments', nargs='*', default=[], help="e.g. In:As:3,Se:1 (repeatable).")
    p.add_argument('--atoms', nargs='*', type=int, default=[], help="1-based indices for --mode atom.")

    # Modes (surface-flux/interface-flux removed)
    p.add_argument('--mode', choices=['field','atom'], default='field')
    p.add_argument('--view', choices=['2d','3d','both'], default='3d')
    p.add_argument('--display', choices=['arrows','spheres'], default='arrows')

    # 3D options
    p.add_argument('--scale-3d', type=float, default=2.0)
    p.add_argument('--sphere-scale', type=float, default=0.2)
    p.add_argument('--show-total-resultant', action='store_true')
    p.add_argument('--resultant-scale', type=float, default=2.0)
    p.add_argument('--resultant-color', default='cyan')
    p.add_argument('--show-snowball', action='store_true')
    p.add_argument('--snowball-scale', type=float, default=2.0)
    p.add_argument('--snowball-color', default='yellow')

    # Plane/slab
    p.add_argument('--plane', choices=['xy','yz','xz','pca'], default='xy')
    p.add_argument('--plane-normal', nargs=3, type=float)
    p.add_argument('--plane-origin', nargs=3, type=float)
    p.add_argument('--slab', type=float, help="Half-thickness Å")
    p.add_argument('--slab-q', type=float, help="Quantile (0..1) by distance to plane")

    # 2D options
    p.add_argument('--figsize', nargs=2, type=float, default=[8,7], help="2D figure size (W H in inches)")
    p.add_argument('--fig-lims', nargs=4, type=float, help="2D window (xmin xmax ymin ymax)")
    p.add_argument('--save-2d')
    p.add_argument('--grid-res', nargs=2, type=int, default=[300,300])
    p.add_argument('--idw-radius', type=float)
    p.add_argument('--idw-power', type=float, default=2.0)
    p.add_argument('--contour', choices=['none','mag','div'], default='none')
    p.add_argument('--contour-levels', type=int, default=15)
    p.add_argument('--contour-alpha', type=float, default=0.55)
    p.add_argument('--map-alpha', type=float, default=0.35)
    p.add_argument('--map-cmap', default='inferno')
    p.add_argument('--stream', action='store_true')
    p.add_argument('--stream-density', type=float, default=1.4)
    p.add_argument('--stream-arrowsize', type=float, default=1.0)
    p.add_argument('--no-scatter', action='store_true')
    p.add_argument('--color-min', type=float); p.add_argument('--color-max', type=float)
    p.add_argument('--arrow-offset', type=float, default=0.0)
    p.add_argument('--dot-radius', type=float, default=0.0)
    p.add_argument('--arrow-headlength', type=float, default=9.0)
    p.add_argument('--arrow-headaxislength', type=float, default=7.0)

    # Tagging + labeled XYZ + KINDs
    p.add_argument('--export-tags-xyz', help='Export an XYZ using per-atom labels (from --tag-*).')
    p.add_argument('--tag-coordination', action='store_true', help='Use labels like In4, As3, ...')
    p.add_argument('--tag-environment', action='store_true', help='Use labels like In2As2Se, ...')
    p.add_argument('--tag-magnitude', nargs='*', help='Bins name:min:max, e.g. low:0:0.2 mid:0.2:0.5 high:0.5:2.0')
    p.add_argument('--kinds-stdout', action='store_true', help='Also print the KIND sections to stdout.')
    p.add_argument('--no-kinds-file', action='store_true', help='Do not write *_kinds.txt automatically.')

    return p.parse_args(argv)

def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    symbols, positions = load_xyz(args.xyz)
    bader = load_bader(args.bader)

    pair_cut = parse_pair_cutoffs(args.pair_cutoff)
    bonds = guess_bonds_rdkit(symbols, positions, pair_cutoffs=pair_cut)

    G, atom_data = build_graph_and_data(symbols, positions, bonds, bader)
    valid_atoms, valid_bonds, A, b = build_global_system(atom_data, G, bonds)
    x, *_ = np.linalg.lstsq(A, b, rcond=None)

    target_filter = select_indices_by_environments(G, symbols, args.target_environments) if args.target_environments else None

    field_full = compute_global_vector_field(valid_atoms, valid_bonds, x, G, positions, symbols,
                                             exclude_elements=args.exclude_elements,
                                             exclude_coord=args.exclude_coordination,
                                             target_filter=target_filter,
                                             atom_data=atom_data)

    # optional visualization slab filter (does not affect math)
    if args.slab is not None or args.slab_q is not None:
        C,e1,e2,n = plane_from_args_public(args.plane, args.plane_normal, args.plane_origin, positions)
        mask = filter_by_slab_public(positions, C, n, half_thickness=args.slab, quantile=args.slab_q)
        field = [e for e in field_full if mask[e['atom_index']]]
    else:
        field = field_full

    # Tagging export workflow
    if args.export_tags_xyz:
        base_out = os.path.splitext(os.path.basename(args.export_tags_xyz))[0]
        labels = list(symbols)
        stats = {}

        if args.tag_magnitude:
            labels, stats = tag_by_magnitude(symbols, field, args.tag_magnitude, base_out=base_out)
        elif args.tag_environment:
            labels, stats = tag_by_environment(symbols, G, field, base_out=base_out)
        elif args.tag_coordination:
            labels, stats = tag_by_coordination(symbols, G, field, base_out=base_out)
        else:
            print("[warn] --export-tags-xyz was set but no --tag-* option provided; writing original symbols.")
            labels = list(symbols)

        export_xyz_with_labels(args.export_tags_xyz, labels, positions)

        # KIND sections for labels (ready to copy/paste)
        if labels != list(symbols) and not args.no_kinds_file:
            write_kind_sections_for_xyz(labels, args.export_tags_xyz)
        if args.kinds_stdout:
            print("\n# ---- CP2K KIND sections (copy/paste) ----")
            print(kind_sections_from_labels(labels))

    # Modes
    if args.mode == 'atom':
        centers = [i-1 for i in args.atoms] if args.atoms else []
        flows_dict = {c: extract_local_flows(c, valid_bonds, x, G) for c in centers}
        viz_3d.plot_atom_mode(None, positions, bonds, symbols, flows_dict, scale=args.scale_3d)
        return

    # field (default) view
    if args.view in ('3d','both'):
        viz_3d.plot_global_vector_field_pyvista(
            positions, field, bonds, show=True, symbols=symbols,
            display_mode=args.display, scale_factor=args.scale_3d, sphere_scale=args.sphere_scale,
            cmin=args.color_min, cmax=args.color_max,
            show_total=args.show_total_resultant, total_scale=args.resultant_scale, total_color=args.resultant_color,
            G=G if args.show_snowball else None, show_snowball=args.show_snowball,
            plane=args.plane, plane_normal=args.plane_normal, plane_origin=args.plane_origin,
            slab=args.slab, slab_q=args.slab_q
        )
    if args.view in ('2d','both'):
        viz_2d.plot_global_vector_field_2d(
            positions, field, plane=args.plane, plane_normal=args.plane_normal, plane_origin=args.plane_origin,
            cmin=args.color_min, cmax=args.color_max, save_path=args.save_2d, show=True, symbols=symbols,
            grid_res=tuple(args.grid_res), idw_radius=args.idw_radius, idw_power=args.idw_power,
            contour=args.contour, contour_levels=args.contour_levels, contour_alpha=args.contour_alpha,
            map_alpha=args.map_alpha, map_cmap=args.map_cmap, stream=args.stream, stream_density=args.stream_density,
            stream_arrowsize=args.stream_arrowsize, figsize=tuple(args.figsize),
            fig_lims=tuple(args.fig_lims) if args.fig_lims else None,
            arrow_offset=args.arrow_offset, dot_radius=args.dot_radius,
            headlength=int(args.arrow_headlength), headaxislength=int(args.arrow_headaxislength)
        )

if __name__ == '__main__':
    main()

    return 0
