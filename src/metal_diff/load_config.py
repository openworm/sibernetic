"""Sibernetic configuration loader for the Metal substrate.

Parses configuration/<scenario> files (the same ones OpenCL consumes)
and writes binary buffers in the format our Metal CLI ops expect.

Sections handled:
  [position]    — N rows of (x, y, z, type)
                  type 1.x = liquid, 2.x = elastic, 3.x = boundary
  [velocity]    — N rows of (vx, vy, vz, _)
  [connection]  — n_elastic × 32 rows of (jd, rij0, val1, val2)
                  jd = neighbor index (-1 = empty slot); rij0 = rest length

Output mapping for our Metal substrate:
  active particles  ← (liquid + elastic), in original index order
  static particles  ← (boundary), in original index order
  bonds             ← decoded elastic connections (each bond once, i<j),
                      where i, j are indices INTO THE ACTIVE BUFFER

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
         'connections': np.ndarray [n_elastic*32, 4]  (jd, rij0, val1, val2)}
    """
    sections = {'box': [], 'pos': [], 'vel': [], 'connections': []}
    section_keys = {
        '[simulation box]': 'box',
        '[position]': 'pos',
        '[velocity]': 'vel',
        '[connection]': 'connections',
        '[membranes]': '_membranes',     # ignored for now
        '[particleMemIndex]': '_pmemidx', # ignored for now
        '[end]': None,
    }
    cur = None
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
            if cur is None or cur.startswith('_'):
                continue
            sections[cur].append(s)

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
    return {'box': box, 'pos': pos, 'vel': vel, 'connections': connections}


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
                 active_idx: np.ndarray) -> np.ndarray:
    """Decode the 32-slot connection table into a deduplicated bonds list.

    Each elastic particle i_elastic (0..n_elastic-1) has 32 connection
    slots. Slot k of particle i is at connection row i*32 + k. A slot
    is "empty" if jd <= 0 (typically -1 in the file).

    The neighbor index `jd` is in GLOBAL indexing (over all particles).
    We need to remap to the ACTIVE buffer's local indexing.

    Returns:
        bonds: structured ndarray with fields (i, j, rest, pad) — the
        same format `xpbd_step` and `step_bond_fwd` already consume.
        Each bond appears once (i < j).
    """
    bonds_set = set()
    bonds_rest = {}
    # Map global → local (active_idx position)
    global_to_local = -np.ones(int(active_idx.max()) + 1, dtype=np.int64)
    for local, glob in enumerate(active_idx):
        if glob < len(global_to_local):
            global_to_local[glob] = local

    for i_elastic in range(n_elastic):
        # In the config, elastic particles are in some position within
        # the global array. The CONNECTION block is indexed by elastic
        # particle, so i_elastic = 0..n_elastic-1 corresponds to the
        # FIRST n_elastic active particles in the original order? Or
        # the first n_elastic of TYPE 2.x?
        # Sibernetic's loader (owConfigProperty.cpp) lays out particles
        # in the order they appear in the file, and elastic particles
        # come FIRST in our demo1 (218 elastic, then 125 liquid, then
        # boundary). So elastic particle i_elastic = global index i_elastic.
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
    return bonds


def load_to_metal_buffers(scenario: str, out_dir: str = "/tmp"):
    """End-to-end: parse config + write binary buffers Metal CLI consumes.

    Returns a dict with the parameters needed for the CLI (n_active,
    n_static, n_bonds, paths) so callers can spawn sib_metal directly.
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

    bonds = decode_bonds(cfg['connections'], n_elastic, active_idx)

    os.makedirs(out_dir, exist_ok=True)
    paths = {
        'pos_active': os.path.join(out_dir, f"{scenario}_pos_active.bin"),
        'vel_active': os.path.join(out_dir, f"{scenario}_vel_active.bin"),
        'pos_static': os.path.join(out_dir, f"{scenario}_pos_static.bin"),
        'bonds':      os.path.join(out_dir, f"{scenario}_bonds.bin"),
    }
    pos_active.tofile(paths['pos_active'])
    vel_active.tofile(paths['vel_active'])
    pos_static.tofile(paths['pos_static'])
    bonds.tofile(paths['bonds'])

    return {
        'scenario': scenario,
        'n_total': n_total,
        'n_active': n_active,
        'n_static': n_static,
        'n_elastic': n_elastic,
        'n_liquid': n_active - n_elastic,
        'n_bonds': len(bonds),
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
    print(f"  Sim box: {info['box']}")
    print(f"  Wrote binary buffers to:")
    for k, v in info['paths'].items():
        sz = os.path.getsize(v)
        print(f"    {k}: {v}  ({sz} bytes)")


if __name__ == "__main__":
    main()
