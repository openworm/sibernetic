"""Sibernetic configuration loader for the Metal substrate.

Parses configuration/<scenario> files (the same ones OpenCL consumes)
and writes binary buffers in the format our Metal CLI ops expect.

Sections handled:
  [position]         — N rows of (x, y, z, type)
                       type 1.x = liquid, 2.x = elastic, 3.x = boundary
  [velocity]         — N rows of (vx, vy, vz, _)
  [connection]       — n_elastic × 32 rows of (jd, rij0, val1, val2)
                       jd = neighbor index (-1 = empty slot); rij0 = rest length
  [membranes]        — N_membrane rows of (i, j, k) elastic-particle
                       triangle vertices (global indices)
  [particleMemIndex] — n_elastic × 7 lines of one int each:
                       per-elastic-particle list of incident membrane
                       indices (-1 sentinel for unused slots)

Output mapping for our Metal substrate:
  active particles  ← (liquid + elastic), in original index order
  static particles  ← (boundary), in original index order
  bonds             ← decoded elastic connections (each bond once, i<j),
                      where i, j are indices INTO THE ACTIVE BUFFER
  membrane_tris     ← int32 [N_membrane, 3], particle indices remapped
                      to the ACTIVE buffer's local indexing
  pmem_index        ← int32 [n_elastic, 7], same indexing as the
                      OpenCL particleMembranesList; -1 means "no
                      membrane here"

Run as a script:
    .venv/bin/python src/metal_diff/load_config.py demo1 [--out /tmp/demo1]
"""
import argparse
import os
import re

import numpy as np


def parse_config(path: str) -> dict:
    """Read a Sibernetic configuration file into a dict of arrays.

    Returns:
        {'box': [x0,x1,y0,y1,z0,z1],
         'pos': np.ndarray [N,4]  (x,y,z,type),
         'vel': np.ndarray [N,4],
         'connections': np.ndarray [n_elastic*32, 4]  (jd, rij0, val1, val2),
         'membranes': np.ndarray [N_membrane, 3] int32,
         'pmem_index': np.ndarray [n_elastic*7] int32}
    """
    sections = {'box': [], 'pos': [], 'vel': [], 'connections': [],
                'membranes': [], 'pmem_index': []}
    section_keys = {
        '[simulation box]': 'box',
        '[position]': 'pos',
        '[velocity]': 'vel',
        '[connection]': 'connections',
        '[membranes]': 'membranes',
        '[particleMemIndex]': 'pmem_index',
        '[end]': None,
    }
    cur = None
    pre_header_floats = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s in section_keys:
                cur = section_keys[s]
                if cur is None:
                    break
                continue
            if cur is None:
                # Some configs (e.g. configuration/test/one_sprig_test) emit
                # box bounds as 6 bare numbers BEFORE [position] without a
                # [simulation box] header. Collect them as box defaults.
                try:
                    pre_header_floats.append(float(s))
                except ValueError:
                    pass
                continue
            sections[cur].append(s)
    # If the file had no [simulation box] but had bare leading numbers, treat
    # the first 6 as the box bounds.
    if not sections['box'] and len(pre_header_floats) >= 6:
        sections['box'] = [str(v) for v in pre_header_floats[:6]]

    box = [float(x) for x in sections['box']]
    pos = np.array(
        [[float(x) for x in line.split()] for line in sections['pos']],
        dtype=np.float32,
    )
    vel = np.array(
        [[float(x) for x in line.split()] for line in sections['vel']],
        dtype=np.float32,
    )
    conn_lines = sections['connections']
    if conn_lines:
        connections = np.array(
            [[float(x) for x in line.split()] for line in conn_lines],
            dtype=np.float32,
        )
    else:
        connections = np.zeros((0, 4), dtype=np.float32)
    if sections['membranes']:
        membranes = np.array(
            [[int(x) for x in line.split()] for line in sections['membranes']],
            dtype=np.int32,
        )
    else:
        membranes = np.zeros((0, 3), dtype=np.int32)
    if sections['pmem_index']:
        # Each line is a single int (-1 sentinel for unused slot).
        pmem_index = np.array(
            [int(line.split()[0]) for line in sections['pmem_index']],
            dtype=np.int32,
        )
    else:
        pmem_index = np.zeros(0, dtype=np.int32)
    return {'box': box, 'pos': pos, 'vel': vel, 'connections': connections,
            'membranes': membranes, 'pmem_index': pmem_index}


def split_by_type(pos: np.ndarray, vel: np.ndarray):
    """Partition particles by type → (active_idx, static_idx).

    Active = liquid (type 1.x) + elastic (type 2.x); these particles
    move under physics. Static = boundary (type 3.x); these particles
    are frozen and only act as neighbors for forces / density.

    Returns:
        active_idx (np.ndarray[N_active]), static_idx (np.ndarray[N_static])
        each indexing into the original [N] arrays.
    """
    types = pos[:, 3]
    is_active = types < 3.0    # liquid + elastic, NOT boundary
    active_idx = np.where(is_active)[0]
    static_idx = np.where(~is_active)[0]
    return active_idx, static_idx


def decode_bonds(connections: np.ndarray, n_elastic: int,
                 active_idx: np.ndarray, all_pos: np.ndarray = None):
    """Decode the 32-slot connection table into elastic-elastic bonds and
    elastic→boundary anchors.

    Each elastic particle i_elastic (0..n_elastic-1) has 32 connection
    slots. Slot k of particle i is at connection row i*32 + k. A slot
    is "empty" if jd <= 0 (typically -1 in the file).

    The neighbor index `jd` is in GLOBAL indexing. Three cases:
    - jd indexes another active particle (elastic) → bond
    - jd indexes a boundary particle (type ≥ 3) → anchor
    - jd <= 0 or out of range → skip

    Returns:
        bonds:   structured ndarray (i, j, rest, pad) — elastic-elastic
                 bonds, deduped (i < j).
        anchors: structured ndarray (i, ax, ay, az, rest, pad, pad, pad)
                 — elastic-to-boundary anchors. `i` is the elastic
                 active-local index; (ax, ay, az) is the FIXED boundary
                 position (the anchor target); rest is the spring rest
                 length. NOT deduped (one row per anchor connection).
                 Pass `all_pos[:, :3]` (shape [N_total, 3]) so we can
                 look up boundary positions by global index.
    """
    bonds_set = set()
    bonds_rest = {}
    anchors_list = []   # list of (local_i, ax, ay, az, rest)

    global_to_local = -np.ones(int(active_idx.max()) + 1, dtype=np.int64)
    for local, glob in enumerate(active_idx):
        if glob < len(global_to_local):
            global_to_local[glob] = local

    for i_elastic in range(n_elastic):
        global_i = i_elastic
        local_i = int(global_to_local[global_i]) if global_i < len(global_to_local) else -1
        if local_i < 0:
            continue
        for slot in range(32):
            row = i_elastic * 32 + slot
            if row >= len(connections):
                break
            jd = connections[row, 0]
            if jd <= 0:
                continue
            global_j = int(jd)
            if all_pos is not None and 0 <= global_j < len(all_pos):
                # If the target is a boundary particle, treat as anchor.
                # all_pos[:, 3] holds the type field.
                target_type = all_pos[global_j, 3] if all_pos.shape[1] > 3 else None
            else:
                target_type = None
            if target_type is not None and target_type >= 3.0:
                # Boundary anchor — capture its position (it's frozen).
                ax, ay, az = (float(all_pos[global_j, 0]),
                              float(all_pos[global_j, 1]),
                              float(all_pos[global_j, 2]))
                rest_len = float(connections[row, 1])
                anchors_list.append((local_i, ax, ay, az, rest_len))
                continue
            # Otherwise it should be an active (elastic) particle.
            if global_j < 0 or global_j >= len(global_to_local):
                continue
            local_j = int(global_to_local[global_j])
            if local_j < 0:
                continue
            a, b = (local_i, local_j) if local_i < local_j else (local_j, local_i)
            if (a, b) in bonds_set:
                continue
            bonds_set.add((a, b))
            bonds_rest[(a, b)] = float(connections[row, 1])

    bonds = np.zeros(len(bonds_set),
                     dtype=[('i', np.int32), ('j', np.int32),
                            ('rest', np.float32), ('pad', np.float32)])
    for k, (i, j) in enumerate(sorted(bonds_set)):
        bonds[k] = (i, j, bonds_rest[(i, j)], 0.0)

    # Anchors: 8-float-aligned struct (i + 1 pad int + 3 anchor + 1 rest + 2 pad = 32B).
    anchors = np.zeros(len(anchors_list),
                       dtype=[('i', np.int32), ('pad0', np.int32),
                              ('ax', np.float32), ('ay', np.float32),
                              ('az', np.float32), ('rest', np.float32),
                              ('pad1', np.float32), ('pad2', np.float32)])
    for k, (i, ax, ay, az, rest) in enumerate(anchors_list):
        anchors[k] = (i, 0, ax, ay, az, rest, 0.0, 0.0)
    return bonds, anchors


def build_static_grid(positions: np.ndarray, h: float):
    """Build a uniform spatial hash grid over the (frozen) static particles.

    Boundary positions don't change so we build the grid once per scenario.

    Returns:
        sorted_pos:  [n, 3] fp32 — positions reordered so particles in the
                     same cell are contiguous in memory
        cell_start:  [n_cells + 1] int32 — particles in cell c are at
                     sorted_pos[cell_start[c] : cell_start[c+1]]
        grid_dim:    [3] int32 — cells per axis
        grid_origin: [3] fp32 — world position of cell (0,0,0)
    """
    box_min = positions.min(axis=0) - 0.1
    box_max = positions.max(axis=0) + 0.1
    grid_dim = np.ceil((box_max - box_min) / h).astype(np.int32)
    n_cells = int(grid_dim.prod())

    cells = np.floor((positions - box_min) / h).astype(np.int32)
    cell_ids = (cells[:, 0]
                + cells[:, 1] * grid_dim[0]
                + cells[:, 2] * grid_dim[0] * grid_dim[1])

    perm = np.argsort(cell_ids, kind='stable')
    sorted_pos = positions[perm].astype(np.float32)
    sorted_ids = cell_ids[perm]

    # cell_start[c] = number of particles with cell_id < c
    counts = np.bincount(sorted_ids, minlength=n_cells)
    cell_start = np.zeros(n_cells + 1, dtype=np.int32)
    cell_start[1:] = np.cumsum(counts)
    return sorted_pos, cell_start, grid_dim, box_min.astype(np.float32)


def load_to_metal_buffers(scenario: str, out_dir: str = "/tmp",
                            h: float = 1.67):
    """End-to-end: parse config + write binary buffers Metal CLI consumes.

    `h` is the SPH smoothing radius — needed for the spatial grid build.

    Returns a dict with the parameters needed for the CLI (n_active,
    n_static, n_bonds, paths, grid info) so callers can spawn sib_metal
    directly.
    """
    config_path = os.path.join("configuration", scenario)
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"no configuration/{scenario}")

    cfg = parse_config(config_path)
    n_total = len(cfg['pos'])

    active_idx, static_idx = split_by_type(cfg['pos'], cfg['vel'])
    n_active = len(active_idx)
    n_static = len(static_idx)
    n_elastic = int(((cfg['pos'][:, 3] >= 2.0) & (cfg['pos'][:, 3] < 3.0)).sum())

    pos_active = cfg['pos'][active_idx, :3].astype(np.float32)
    vel_active = cfg['vel'][active_idx, :3].astype(np.float32)
    pos_static = cfg['pos'][static_idx, :3].astype(np.float32)

    bonds, anchors = decode_bonds(cfg['connections'], n_elastic, active_idx,
                                    all_pos=cfg['pos'])

    # Membrane topology: triangle particle IDs are GLOBAL in the config;
    # remap to local active-buffer indexing (same convention as bonds).
    # OpenCL keeps them global because its position buffer is also
    # global; our active-only buffer needs the remap. demo1's elastic
    # particles come first in the global array, so global_to_local is
    # the identity for elastic indices, but we don't rely on that.
    membranes_global = cfg['membranes']
    if len(membranes_global):
        global_to_local = -np.ones(int(active_idx.max()) + 1, dtype=np.int64)
        for local, glob in enumerate(active_idx):
            if glob < len(global_to_local):
                global_to_local[glob] = local
        membrane_tris = np.zeros_like(membranes_global)
        for ti in range(len(membranes_global)):
            for vi in range(3):
                g = int(membranes_global[ti, vi])
                if 0 <= g < len(global_to_local):
                    membrane_tris[ti, vi] = global_to_local[g]
                else:
                    membrane_tris[ti, vi] = -1
    else:
        membrane_tris = np.zeros((0, 3), dtype=np.int32)
    n_membranes = len(membrane_tris)

    # Per-particle membrane index list — (n_elastic * 7) int32, stored
    # row-major so pmem_index[elastic_i * 7 + slot] is the membrane
    # index (or -1) for that slot. Indexing matches OpenCL exactly.
    pmem_index = cfg['pmem_index'].astype(np.int32)
    expected_pmem_len = n_elastic * 7
    if len(pmem_index) and len(pmem_index) != expected_pmem_len:
        # Some configs may have trailing -1 lines truncated; pad.
        if len(pmem_index) < expected_pmem_len:
            pad = np.full(expected_pmem_len - len(pmem_index), -1, dtype=np.int32)
            pmem_index = np.concatenate([pmem_index, pad])
        else:
            pmem_index = pmem_index[:expected_pmem_len]

    # Build static-particle spatial grid (one-time per scenario).
    sorted_static, cell_start, grid_dim, grid_origin = build_static_grid(
        pos_static, h)
    n_cells = int(grid_dim.prod())

    os.makedirs(out_dir, exist_ok=True)
    paths = {
        'pos_active':    os.path.join(out_dir, f"{scenario}_pos_active.bin"),
        'vel_active':    os.path.join(out_dir, f"{scenario}_vel_active.bin"),
        'pos_static':    os.path.join(out_dir, f"{scenario}_pos_static.bin"),
        'bonds':         os.path.join(out_dir, f"{scenario}_bonds.bin"),
        'anchors':       os.path.join(out_dir, f"{scenario}_anchors.bin"),
        'sorted_static': os.path.join(out_dir, f"{scenario}_sorted_static.bin"),
        'cell_start':    os.path.join(out_dir, f"{scenario}_cell_start.bin"),
        'membranes':     os.path.join(out_dir, f"{scenario}_membranes.bin"),
        'pmem_index':    os.path.join(out_dir, f"{scenario}_pmem_index.bin"),
    }
    pos_active.tofile(paths['pos_active'])
    vel_active.tofile(paths['vel_active'])
    pos_static.tofile(paths['pos_static'])
    bonds.tofile(paths['bonds'])
    anchors.tofile(paths['anchors'])
    sorted_static.tofile(paths['sorted_static'])
    cell_start.tofile(paths['cell_start'])
    membrane_tris.astype(np.int32).tofile(paths['membranes'])
    pmem_index.astype(np.int32).tofile(paths['pmem_index'])

    return {
        'scenario': scenario,
        'n_total': n_total,
        'n_active': n_active,
        'n_static': n_static,
        'n_elastic': n_elastic,
        'n_liquid': n_active - n_elastic,
        'n_bonds': len(bonds),
        'n_anchors': len(anchors),
        'n_membranes': n_membranes,
        'n_cells': n_cells,
        'grid_dim': grid_dim.tolist(),
        'grid_origin': grid_origin.tolist(),
        'h': h,
        'box': cfg['box'],
        'paths': paths,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('scenario', help="e.g. demo1, worm")
    ap.add_argument('--out', default="/tmp", help="output directory for binary buffers")
    args = ap.parse_args()

    info = load_to_metal_buffers(args.scenario, args.out)
    print(f"Loaded configuration/{info['scenario']}")
    print(f"  Total particles: {info['n_total']}")
    print(f"    elastic:  {info['n_elastic']}")
    print(f"    liquid:   {info['n_liquid']}")
    print(f"    boundary: {info['n_static']}")
    print(f"  Active (elastic+liquid): {info['n_active']}")
    print(f"  Bonds (deduped):         {info['n_bonds']}")
    print(f"  Boundary anchors:        {info['n_anchors']}")
    print(f"  Membrane triangles:      {info['n_membranes']}")
    print(f"  Sim box: {info['box']}")
    print(f"  Wrote binary buffers to:")
    for k, v in info['paths'].items():
        sz = os.path.getsize(v)
        print(f"    {k}: {v}  ({sz} bytes)")


if __name__ == "__main__":
    main()
