"""Validate the CUDA wpoly6_inplace kernel against a numpy reference.

Run after building the binary:
    src/cuda_diff/build.bat       (Windows)
    src/cuda_diff/build.sh        (Linux)
    python src/cuda_diff/test_wpoly6_cuda.py

Mirrors src/metal_diff/test_dist.py::test_wpoly6 — same numpy reference,
same tolerance — but Windows-portable (uses tempfile, finds .exe or no-ext
binary). Local-only: not part of any PR yet.
"""
import math
import os
import platform
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY_NAME = "sib_cuda.exe" if platform.system() == "Windows" else "sib_cuda"
BINARY = os.path.join(HERE, BINARY_NAME)
TMP = tempfile.gettempdir()


def run_wpoly6(r2: np.ndarray, h: float, iters: int = 1):
    inout_path = os.path.join(TMP, "sib_cuda_wpoly6.bin")
    r2.astype(np.float32).tofile(inout_path)
    t0 = time.perf_counter()
    args = [BINARY, "wpoly6_inplace", str(r2.size), str(h), inout_path]
    if iters > 1:
        args.append(str(iters))
    subprocess.run(args, check=True)
    wall = time.perf_counter() - t0
    out = np.fromfile(inout_path, dtype=np.float32).reshape(r2.shape)
    return out, wall


def numpy_wpoly6(r2: np.ndarray, h: float):
    h2 = h * h
    poly6_const = 315.0 / (64.0 * math.pi * h ** 9)
    diff = h2 - r2
    out = poly6_const * diff ** 3
    out = np.where(r2 >= h2, 0.0, out)
    return out.astype(np.float32)


def test_wpoly6(name, n_active, n_static, h, tolerance=1e-3):
    np.random.seed(0)
    active = np.random.randn(n_active, 3).astype(np.float32) * (h * 0.7)
    static_p = np.random.randn(n_static, 3).astype(np.float32) * (h * 0.7)
    diff = active[:, None, :] - static_p[None, :, :]
    r2 = np.sum(diff * diff, axis=-1).astype(np.float32)

    cuda_W, cuda_wall = run_wpoly6(r2, h)
    expected_W = numpy_wpoly6(r2, h)

    err_max = float(np.abs(cuda_W - expected_W).max())
    err_mean = float(np.abs(cuda_W - expected_W).mean())
    nonzero_frac = float((cuda_W > 0).sum()) / cuda_W.size

    print(f"=== {name} (n={r2.size}, h={h}) ===")
    print(f"  Nonzero W fraction: {nonzero_frac:.1%}")
    print(f"  Max error:          {err_max:.6e}")
    print(f"  Mean error:         {err_mean:.6e}")
    print(f"  CUDA wall (incl IO): {cuda_wall*1000:.1f} ms")
    passed = err_max < tolerance
    print(f"  {'[PASS]' if passed else '[FAIL]'}")
    print()
    return passed


def main():
    if not os.path.exists(BINARY):
        print(f"sib_cuda binary missing at {BINARY}")
        print(f"Build first: {os.path.join(HERE, 'build.bat' if platform.system() == 'Windows' else 'build.sh')}")
        return 1

    h = 5.0
    all_pass = True
    all_pass &= test_wpoly6("M6.1 small (50x80)",                50,  80, h)
    all_pass &= test_wpoly6("M6.1 mid (343x500)",               343, 500, h)
    all_pass &= test_wpoly6("M6.1 demo1-realistic (343x17498)", 343, 17498, h)

    # Steady-state perf check (iters=200) on the demo1-scale buffer.
    np.random.seed(0)
    active = np.random.randn(343, 3).astype(np.float32) * (h * 0.7)
    static_p = np.random.randn(17498, 3).astype(np.float32) * (h * 0.7)
    diff = active[:, None, :] - static_p[None, :, :]
    r2 = np.sum(diff * diff, axis=-1).astype(np.float32)
    print("=== M6.1 steady-state perf (343x17498, 200 iters) ===")
    _, _ = run_wpoly6(r2, h, iters=200)
    print()

    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
