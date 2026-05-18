"""TRUE analytic-gradient SGD wrapper for the native Metal substrate.

Each iteration runs ``xpbd_full_fwd`` at the current (rho, K, v, alpha)
to produce a saved state file plus the final trajectory frame, computes
a scalar loss + ∂L/∂x_final, runs ``xpbd_full_bwd`` to get the analytic
parameter gradients, and applies an Adam update in log-space. Gradients
are EXACT (modulo the ~17% stop-gradient through density on visc_K —
see ``test_param_grads.py``).

Cost: K=50000 sim steps → ~90s fwd + ~8min bwd ≈ 10min per iter on
Apple M-series.

Implementation note: this script is a thin wrapper around
``_sgd_true_impl.py`` — the substantive SGD/Adam/loss code is shared
with ``src/cuda_diff/sgd_true.py``. Only ``BIN`` and ``TMP`` differ.
``TMP`` here is hardcoded to ``/tmp`` (preserving the historical
Metal-substrate scratch location); set TMPDIR if you want it
elsewhere via ``tempfile.gettempdir()``.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "sib_metal")
TMP = "/tmp"

sys.path.insert(0, HERE)
from load_config import load_to_metal_buffers  # noqa: E402
from _sgd_true_impl import main as _impl_main  # noqa: E402


if __name__ == "__main__":
    sys.exit(_impl_main(BIN, TMP, load_to_metal_buffers))
