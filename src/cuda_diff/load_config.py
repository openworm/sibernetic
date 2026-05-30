"""Sibernetic config loader for the native CUDA substrate.

Thin re-export of the canonical parser at
``src/metal_diff/load_config.py``. The config format is
substrate-agnostic, so there is no reason to maintain two copies; this
shim keeps the existing import idiom

    sys.path.insert(0, HERE)
    from load_config import load_to_metal_buffers, parse_config

working from anywhere under ``src/cuda_diff/`` while the actual parser
lives in one place.

Why ``importlib`` instead of a regular ``from load_config import *``:
when a caller has already added ``src/cuda_diff/`` to ``sys.path`` and
imported ``load_config``, *this* module ends up in ``sys.modules`` under
the name ``load_config``. A plain ``from load_config import *`` from
inside this file would then re-bind the partially-initialised shim to
itself, not the metal_diff version. Loading the metal_diff file under a
unique module name via ``importlib.util.spec_from_file_location``
sidesteps the cache collision.
"""
import importlib.util as _importlib_util
import os as _os

_THIS_DIR = _os.path.dirname(_os.path.abspath(__file__))
_METAL_LOAD_CONFIG = _os.path.normpath(
    _os.path.join(_THIS_DIR, _os.pardir, "metal_diff", "load_config.py"))

_spec = _importlib_util.spec_from_file_location(
    "_sib_load_config_metal_canonical", _METAL_LOAD_CONFIG)
if _spec is None or _spec.loader is None:
    raise ImportError(
        f"cuda_diff/load_config.py shim could not locate the canonical "
        f"parser at {_METAL_LOAD_CONFIG}. If you moved or renamed "
        f"metal_diff/load_config.py, update this path.")
_module = _importlib_util.module_from_spec(_spec)
_spec.loader.exec_module(_module)

# Re-export the public surface so `from load_config import X` works the
# same way as before. Pulling explicit names (not `*`) keeps static
# analysers happy and documents what cuda_diff actually consumes.
parse_config = _module.parse_config
load_to_metal_buffers = _module.load_to_metal_buffers

# Anything else the canonical module exposes also becomes available, so
# callers that grew new imports don't silently break when the shim is
# in place.
for _name in dir(_module):
    if not _name.startswith("_") and _name not in globals():
        globals()[_name] = getattr(_module, _name)

del _importlib_util, _os, _THIS_DIR, _METAL_LOAD_CONFIG, _spec, _module, _name
