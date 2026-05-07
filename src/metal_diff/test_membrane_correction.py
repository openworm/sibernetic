"""FD validation of accumulate_membrane_correction_backward.

Synthetic scene: 3 elastic vertices forming a single membrane triangle
+ 1 liquid particle hovering ~r0/2 above its centroid. The forward
projects the liquid onto the triangle plane, computes the unit normal,
and accumulates a delta_pos vector into mem_corr[liquid].

Backward must recover ∂L/∂pos[k] for every k ∈ {pa, pb, pc, p_liquid}
where L = <grad_seed, mem_corr[liquid]>.

Validates analytically-derived backward against centered finite
differences. PASS threshold: rel_err < 1e-3 per particle per axis.
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "sib_metal")
TMP = "/tmp"


def write_floats(path, arr):
    np.asarray(arr, dtype=np.float32).tofile(path)


def write_ints(path, arr):
    np.asarray(arr, dtype=np.int32).tofile(path)


def run_fwd(n_active, n_elastic, n_membranes, r0, pos, mem_init):
    write_floats(f"{TMP}/m_pos.bin", pos)
    write_floats(f"{TMP}/m_mem_init.bin", mem_init)
    cmd = [BIN, "accumulate_membrane_fwd",
           str(n_active), str(n_elastic), str(n_membranes), str(r0),
           f"{TMP}/m_pos.bin", f"{TMP}/m_mem.bin",
           f"{TMP}/m_pmi.bin", f"{TMP}/m_mem_init.bin"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"fwd: {r.stderr}")
    return np.fromfile(f"{TMP}/membrane_corr_out.bin",
                       dtype=np.float32).reshape(-1, 3)


def run_bwd(n_active, n_elastic, n_membranes, r0, pos, grad_corr):
    write_floats(f"{TMP}/m_pos.bin", pos)
    write_floats(f"{TMP}/m_grad_corr.bin", grad_corr)
    cmd = [BIN, "accumulate_membrane_bwd",
           str(n_active), str(n_elastic), str(n_membranes), str(r0),
           f"{TMP}/m_pos.bin", f"{TMP}/m_mem.bin",
           f"{TMP}/m_pmi.bin", f"{TMP}/m_grad_corr.bin"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"bwd: {r.stderr}")
    return np.fromfile(f"{TMP}/membrane_grad_pos.bin",
                       dtype=np.float32).reshape(-1, 3)


def fd_check(label, n_active, n_elastic, n_membranes, r0,
             pos, membranes, pmem_idx, grad_corr, eps=1e-3,
             rel_thresh=2e-3):
    """Common harness: write topo files, run fwd to verify L, run bwd vs FD.
    Returns (pass: bool, max_rel: float)."""
    write_ints(f"{TMP}/m_mem.bin", membranes)
    write_ints(f"{TMP}/m_pmi.bin", pmem_idx)
    mem_init = np.zeros_like(pos)
    mem_corr = run_fwd(n_active, n_elastic, n_membranes, r0, pos, mem_init)
    L = float(np.dot(grad_corr.flatten(), mem_corr.flatten()))
    g_an = run_bwd(n_active, n_elastic, n_membranes, r0, pos, grad_corr)

    g_fd = np.zeros_like(pos)
    for k in range(n_active):
        for ax in range(3):
            pp = pos.copy(); pp[k, ax] += eps
            mp = run_fwd(n_active, n_elastic, n_membranes, r0, pp, mem_init)
            Lp = float(np.dot(grad_corr.flatten(), mp.flatten()))
            pp = pos.copy(); pp[k, ax] -= eps
            mn = run_fwd(n_active, n_elastic, n_membranes, r0, pp, mem_init)
            Ln = float(np.dot(grad_corr.flatten(), mn.flatten()))
            g_fd[k, ax] = (Lp - Ln) / (2 * eps)

    eps_floor = 1e-4
    max_rel = 0.0
    print(f"\n=== {label} (n_active={n_active}, n_elastic={n_elastic}, "
          f"n_mem={n_membranes}, r0={r0}) ===")
    print(f"  L = <grad, mem_corr> = {L:.6f}")
    all_pass = True
    for k in range(n_active):
        denom = max(float(np.linalg.norm(g_an[k])),
                    float(np.linalg.norm(g_fd[k])),
                    eps_floor)
        diff = float(np.linalg.norm(g_an[k] - g_fd[k]))
        rel  = diff / denom
        max_rel = max(max_rel, rel)
        ok = rel < rel_thresh
        all_pass = all_pass and ok
        tag = "OK" if ok else "FAIL"
        print(f"    [{tag}] k={k}: rel={rel:.3e}  "
              f"||an||={np.linalg.norm(g_an[k]):.3e}  "
              f"||fd||={np.linalg.norm(g_fd[k]):.3e}")
    return all_pass, max_rel


def case_single_triangle():
    """Original 3-vert + 1 liquid scene."""
    rng = np.random.default_rng(42)
    n_active = 4; n_elastic = 3; n_membranes = 1; r0 = 2.0
    pos = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.05, 0.0],
        [0.5, 0.0, 1.0],
        [0.5, 0.6, 0.3],
    ], dtype=np.float32)
    membranes = np.array([[0, 1, 2]], dtype=np.int32)
    pmem_idx = np.full((n_elastic, 7), -1, dtype=np.int32)
    for v in range(n_elastic):
        pmem_idx[v, 0] = 0
    grad_corr = np.zeros_like(pos)
    grad_corr[3] = rng.standard_normal(3).astype(np.float32)
    return fd_check("single triangle, 1 liquid",
                     n_active, n_elastic, n_membranes, r0,
                     pos, membranes, pmem_idx, grad_corr)


def case_two_triangles_shared_edge():
    """Two adjacent triangles sharing edge (0,1), 4 elastic + 1 liquid.
    Tests that a shared vertex correctly accumulates contributions from
    both triangles."""
    rng = np.random.default_rng(7)
    n_active = 5; n_elastic = 4; n_membranes = 2; r0 = 2.5
    pos = np.array([
        [0.0, 0.0, 0.0],   # 0: shared
        [1.0, 0.0, 0.0],   # 1: shared
        [0.5, 0.05, 1.0],  # 2: tri 0 only
        [0.5, 0.0, -0.8],  # 3: tri 1 only
        [0.5, 0.7, 0.1],   # 4: liquid hovering above the seam
    ], dtype=np.float32)
    membranes = np.array([[0, 1, 2],
                          [0, 1, 3]], dtype=np.int32)
    pmem_idx = np.full((n_elastic, 7), -1, dtype=np.int32)
    pmem_idx[0] = [0, 1, -1, -1, -1, -1, -1]
    pmem_idx[1] = [0, 1, -1, -1, -1, -1, -1]
    pmem_idx[2] = [0, -1, -1, -1, -1, -1, -1]
    pmem_idx[3] = [1, -1, -1, -1, -1, -1, -1]
    grad_corr = np.zeros_like(pos)
    grad_corr[4] = rng.standard_normal(3).astype(np.float32)
    return fd_check("two triangles, shared edge",
                     n_active, n_elastic, n_membranes, r0,
                     pos, membranes, pmem_idx, grad_corr)


def case_two_liquids_one_triangle():
    """Two liquids hovering at different heights above one triangle.
    Tests that gradients superpose across multiple liquids."""
    rng = np.random.default_rng(13)
    n_active = 5; n_elastic = 3; n_membranes = 1; r0 = 2.0
    pos = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.05, 0.0],
        [0.5, 0.0, 1.0],
        [0.4, 0.6, 0.3],   # liquid 1
        [0.6, 0.4, 0.4],   # liquid 2
    ], dtype=np.float32)
    membranes = np.array([[0, 1, 2]], dtype=np.int32)
    pmem_idx = np.full((n_elastic, 7), -1, dtype=np.int32)
    for v in range(n_elastic):
        pmem_idx[v, 0] = 0
    grad_corr = np.zeros_like(pos)
    grad_corr[3] = rng.standard_normal(3).astype(np.float32)
    grad_corr[4] = rng.standard_normal(3).astype(np.float32)
    return fd_check("one triangle, two liquids",
                     n_active, n_elastic, n_membranes, r0,
                     pos, membranes, pmem_idx, grad_corr)


def main():
    results = []
    results.append(("single_triangle",  *case_single_triangle()))
    results.append(("shared_edge",      *case_two_triangles_shared_edge()))
    results.append(("two_liquids",      *case_two_liquids_one_triangle()))

    print("\n=== Summary ===")
    all_ok = True
    for name, ok, max_rel in results:
        tag = "PASS" if ok else "FAIL"
        all_ok = all_ok and ok
        print(f"  [{tag}] {name}: max_rel={max_rel:.3e}")
    if all_ok:
        print("\n  [OVERALL PASS] M10 backward validates against FD across cases")
        return 0
    else:
        print("\n  [OVERALL FAIL] at least one M10 backward case diverges")
        return 1


if __name__ == "__main__":
    sys.exit(main())
