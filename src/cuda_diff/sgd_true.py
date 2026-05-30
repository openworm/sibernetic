"""TRUE analytic-gradient SGD wrapper for the native CUDA substrate.

Each iteration runs ``xpbd_full_fwd`` at the current (rho, K, v, alpha)
to produce a saved state file plus the final trajectory frame, computes
a scalar loss + ∂L/∂x_final, runs ``xpbd_full_bwd`` to get the analytic
parameter gradients, and applies an Adam update in log-space. Gradients
are EXACT (modulo the ~17% stop-gradient through density on visc_K —
see ``test_param_grads.py``).

Cost: K=50000 sim steps → ~90s fwd + ~8min bwd ≈ 10min per iter on
RTX 2070.

Implementation note: this script is a thin wrapper around
``src/metal_diff/_sgd_true_impl.py`` — the substantive SGD/Adam/loss
code is shared between the Metal and CUDA substrates. Only ``BIN`` and
``TMP`` differ.
"""
import os
import platform
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE,
                   "sib_cuda.exe" if platform.system() == "Windows" else "sib_cuda")
TMP = tempfile.gettempdir()

# Import the substrate's load_to_metal_buffers via the cuda_diff shim
# (which re-exports the canonical metal_diff version).
sys.path.insert(0, HERE)
from load_config import load_to_metal_buffers  # noqa: E402

# Pull the shared SGD impl from metal_diff/.
_METAL_DIFF = os.path.normpath(os.path.join(HERE, os.pardir, "metal_diff"))
if _METAL_DIFF not in sys.path:
    sys.path.insert(0, _METAL_DIFF)
from _sgd_true_impl import main as _impl_main  # noqa: E402


if __name__ == "__main__":
    sys.exit(_impl_main(BIN, TMP, load_to_metal_buffers))
