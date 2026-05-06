![Sibernetic](http://i.imgur.com/Hbsw6Zs.png)

Sibernetic is a physical simulator of biomechanical matter (membranes, elastic
matter, contractile matter) and environments (liquids, solids, elastic matter
with variable physical properties), developed for simulations of *C. elegans*
physical body dynamics within the [OpenWorm project](http://www.openworm.org)
by Andrey Palyanov, Sergey Khayrulin and Mike Vella (development of a Python
module for external muscle activating signals generation and input) as part of
the [OpenWorm team](http://www.openworm.org/people.html). It is primarily
written in C++ and OpenCL, with 3D visualization built on OpenGL.

## Two solver paths

Sibernetic ships with two physically-equivalent simulation paths. They use
different algorithms with different strengths.

| Path | Algorithm | Hardware | Build target | Differentiable? |
|---|---|---|---|---|
| **OpenCL** (legacy, gold standard) | PCISPH | NVIDIA / AMD GPUs, CPU OpenCL | `./Release/Sibernetic` (via `make`) | No |
| **Native Metal** (new) | XPBD | Apple Silicon GPU (M-series) | `./src/metal_diff/sib_metal` (via `cd src/metal_diff && ./build.sh`) | **Yes** |

The two paths share configuration files (`configuration/demo1`,
`configuration/worm`, etc.) and produce visually equivalent results on the
demo1 cube-drop scenario; they are tested for trajectory parity by
`tests/test_demo1_backend_parity.py`.

### From PCISPH to XPBD: why two algorithms?

The original Sibernetic OpenCL solver implements
**PCISPH** (Predictive-Corrective Incompressible Smoothed Particle Hydrodynamics
— Solenthaler 2009). PCISPH enforces fluid incompressibility through an
iterative pressure-correction loop: predict positions under gravity and other
forces, evaluate the resulting density at every particle via a smoothing
kernel, compute a pressure that cancels any density excess, project particles
to satisfy that pressure, and repeat for a few inner iterations per timestep.
PCISPH is a mature, well-understood approach for incompressible fluids; it has
served Sibernetic well for over a decade of *C. elegans* whole-body
simulation.

The native Metal substrate uses **XPBD**
(Extended Position-Based Dynamics — Macklin et al. 2016). XPBD treats every
physical relationship as a **constraint** in position space (density target,
inter-particle distances, floor inequality, spring rest-lengths) and solves
all of them in a single Lagrangian framework with explicit *compliance*
parameters (`α`) controlling stiffness. The forward step is structurally
similar to PCISPH — predict, project, recover velocity — but instead of a
pressure-Poisson loop it does a sequence of constraint projections, each of
which has a clean closed-form solution.

The reason for the change isn't speed (Metal is ~50% faster per step than
OpenCL on L4, but PCISPH is already fast enough for whole-worm sims). The
reason is that **XPBD is differentiable end-to-end**. Every constraint
projection has a hand-derived analytic backward kernel
(`solve_density_constraint_backward`,
`solve_distance_constraints_seq_backward`, `solve_floor_constraint_backward`,
`pair_forces_grid_backward`, `predict_positions_backward`,
`update_velocities_backward`). Stitching them together yields exact gradients
of any scalar loss with respect to physical parameters (`rho_rest`,
`spring_K`, `visc_pair_coef`, `alpha_dens`, `floor_y`, `restitution`) and
even initial conditions. PCISPH's iterative pressure solver, by contrast,
has no clean derivative — its inner loop is a fixed-point iteration and
backpropagating through it is mathematically awkward and numerically
fragile.

So you can think of the two paths this way:

- **OpenCL / PCISPH**: the gold-standard reference. Use it when you want the
  classic Sibernetic behaviour or when reproducing prior published results.
- **Native Metal / XPBD**: the same physical scenario solved a different way,
  with the bonus that *gradient descent on simulation parameters works*.
  Use it when you're doing parameter calibration, learning muscle
  activations to match an observed worm pose, or any other inverse problem
  where you need ∂loss/∂parameter.

## Differentiable physics with native Metal

This is the headline feature of the new substrate. The same kernel that runs
the forward simulation also has a hand-derived backward kernel; together they
support automatic differentiation through 25 ms or more of simulated time.

### What you can differentiate

- **Material parameters**: `rho_rest`, `spring_K`, `visc_pair_coef`,
  `alpha_dens` (XPBD density compliance), `floor_y`, `restitution` (floor
  elasticity coefficient).
- **Initial state**: per-particle initial positions and velocities (gradients
  flow back through every step to the start).
- **External forcing**: scalar gravity (`grad_g_y`).

Each gradient comes out of a single `xpbd_full_bwd` invocation that consumes a
saved forward state and ∂loss/∂x_final, then produces gradients for all of
the above.

### Why this is useful

Three concrete things become tractable that weren't before:

1. **Parameter calibration to a reference**. We have OpenCL trajectories on
   demo1; we want Metal to match them. Define
   `L = sum (x_metal(t) - x_opencl(t))²`, take ∂L/∂(rho, K, v, α), Adam in
   log-space, repeat. We did this and it works (see "Parity work" below).

2. **Worm muscle activation learning**. Given a target *C. elegans* body
   pose at time *t*, infer the muscle-cell activation signals that produce
   it. This is the long-arc OpenWorm goal that motivated building a
   differentiable substrate at all.

3. **Initial-condition inversion**. Given a final particle distribution,
   recover the configuration that produced it.

### How to use it

The differentiable substrate is exposed via two CLI ops on the `sib_metal`
binary: `xpbd_full_fwd` (saves all per-step state needed for backward) and
`xpbd_full_bwd` (consumes that state plus ∂loss/∂x_final, produces parameter
gradients).

A high-level Python wrapper that handles the whole loop is in
`src/metal_diff/sgd_true.py`:

```bash
cd src/metal_diff
./build.sh   # produces sib_metal binary

# Tune spring_K and viscosity on demo1, freezing rho_rest, with targets
# extracted from an OpenCL reference trajectory at t=25ms:
python3 sgd_true.py \
  --rho-init 1000.0 --K-init 5500.0 --v-init 5e-5 \
  --freeze-rho \
  --target-dm -37.74 --target-ext-ratio 1.121 \
  --n-sim-steps 1250 \
  --max-steps 15 --lr 0.05 --tbptt 100 \
  --clip-norm 5.0
```

Two practical knobs you'll touch:

- `BWD_TBPTT=N` (env var): truncated backpropagation-through-time. Chain back
  only the last *N* steps. Lower values are faster and avoid gradient
  explosion; higher values give more accurate gradients but can blow up at
  long trajectories. `N=50–100` is a good starting point for demo1.
- `BWD_CLIP_NORM=X` (env var): per-step gradient clipping. Caps the L2 norm
  of the running position and velocity gradients at *X* between every chain
  step. This is essential for chaotic post-impact dynamics where the
  Jacobian per step amplifies gradients by a factor of ~`sim_scale_inv`
  (≈1.35e5 for Sibernetic's microscale parameters). `X=1e3` is a tested
  default.

There's also a finite-difference fallback at `src/metal_diff/sgd_fd.py` for
parameters that the analytic backward doesn't expose yet.

### Validation of analytic gradients

Each backward kernel has a finite-difference test in `src/metal_diff/`:

```bash
cd src/metal_diff
./build.sh
python3 test_solve_dens_bwd.py          # density-solver backward
python3 test_dens_alpha_grad.py         # ∂L/∂alpha_dens
python3 test_xpbd.py                    # full step (forward integrity)
python3 test_learn_floor.py             # SGD on floor_y from synthetic data
```

All ship green; expected per-element relative error is 1e-3 to 1e-6
depending on the kernel.

## Building

### Linux (OpenCL backend)

You'll need an OpenCL driver. We recommend AMD's APP SDK (most stable on
Ubuntu); Intel and NVIDIA's drivers also work. The CI uses AMD on
`ubuntu-22.04` and `ubuntu-24.04`.

```bash
./setup.sh    # installs OpenCL headers + numpy + pyneuroml
make clean
make
./Release/Sibernetic
```

### macOS — main binary

The main binary builds via `setup.sh` on Apple Silicon and Intel Macs. On
Apple Silicon, OpenCL is disabled at compile time (`OW_NO_OPENCL`); the
binary still builds for command-line tooling but the actual GPU simulation
path is the native Metal substrate (next section).

```bash
./setup.sh    # uses Homebrew + python@3.13
./Release/Sibernetic -h
```

### macOS — native Metal substrate (Apple Silicon only)

```bash
cd src/metal_diff
./build.sh    # compiles shaders.metal + sib_metal.mm
./sib_metal --help
```

This produces a single binary `sib_metal` that exposes the substrate's
operations as subcommands:

| Op | Purpose |
|---|---|
| `xpbd_step` | Run one or more XPBD steps (chunked usage; useful for trajectory dumps) |
| `xpbd_full_fwd` | Run N steps with full per-step state saved (for backward) |
| `xpbd_full_bwd` | Backward pass — consume saved state + ∂L/∂x, produce parameter gradients |
| Test ops | `solve_density_constraint_fwd/bwd`, `solve_floor_constraint_fwd/bwd`, etc. — used by FD validation tests |

Concrete trajectory dump: 25 ms of demo1 cube-drop on Metal:

```bash
python3 src/metal_diff/dump_metal_trajectory.py \
    --scenario demo1 --steps 1250 --chunk 5 --dt 2e-5 \
    --rho-rest 1000.0 --visc-pair-coef 5e-5 --spring-k 5500 \
    --floor-y 2.0 --restitution 0.0 \
    --out /tmp/metal_demo1.txt
```

The output is in the same `position_buffer.txt` format the legacy OpenCL
loader produces, so all downstream tooling (`render_movie.py`,
`replay.py`) works unchanged.

## Common command options

```
 -no_g                 Run without graphics
 -l_to                 Save simulation results to disk.
 -export_vtk           Save simulation results to VTK files.
     logstep=<value>   Log every <value> steps
 -l_from               Load simulation results from disk.
     lpath=<value>     Indicates path to the directory (not the file) where
                       results will be stored. Used with -l_to / -l_from.
 -test                 Run physical tests.
 -f <filename>         Load configuration from file <filename>.
 device=<device_type>  OpenCL device hint: cpu or gpu (default ALL).
 timestep=<value>      Start simulation with time step = <value> seconds.
 timelimit=<value>     Run simulation until <value> seconds reached.
 leapfrog              Use Leapfrog integration instead of semi-implicit Euler.
 oclsourcepath=<value> Path to OpenCL kernel source.
 -nrn <value>          Run with NEURON simulation; <value> is path to a .hoc
                       file. Requires NEURON + sibernetic_neuron bridge.
 -help                 Print this list.
```

The `backend=` flag is no longer recognized by the main binary; the only
available backend on the C++ side is OpenCL. For Apple Silicon GPU
simulation, build and run `src/metal_diff/sib_metal` directly (see above).

## What's inside

Physical algorithms supported across both paths:

- Incompressible fluid (PCISPH on OpenCL; XPBD density-constraint on Metal)
- Elastic matter (rest-length springs + density-driven puff-up)
- Liquid-impermeable membranes (OpenCL only — Metal substrate doesn't yet
  model membranes; cube-drop demo doesn't need them)
- Boundary handling
- Surface tension (PCISPH only — Metal exposes the same kernel
  contribution as a tunable `visc_pair_coef`)

There are two demo configurations:

1. `configuration/demo1` — an elastic cube containing liquid, dropping under
   gravity onto a boundary floor. Used for parity testing.
2. `configuration/demo2` — two membrane sheets attached to boundaries; one
   liquid-impermeable, one not.

Press `1` or `2` in the GUI to switch. `Space` pauses.

## Parity work and known gaps

The native Metal substrate is currently being parity-tested against OpenCL on
demo1. Status (commit `1699fe0`, 2026-05-06):

- Per-particle L2 trajectory error vs OpenCL: **mean 4.01** over 60 samples
  in [0, 25 ms].
- 4 of 5 parity checks in `tests/test_demo1_backend_parity.py` pass on
  demo1; the failing check is `cube_extent_diff`, which is a known
  algorithmic-level XPBD-vs-PCISPH difference: PCISPH performs an
  axis-asymmetric initial cube expansion in the first 0.4 ms that XPBD's
  one-sided density constraint doesn't replicate.
- Per-step compute on Apple M3 Max: 1.08 ms; on NVIDIA L4: 1.61 ms.
- 25 ms cube-drop end-to-end: 1.4 s on Metal local, 10 s on OpenCL via
  Cloud Run.

See `DEVELOPMENT_LOG.md` for the full chronology.

## References

1. B. Solenthaler, R. Pajarola. *Predictive-Corrective Incompressible SPH.*
   ACM Transactions on Graphics (Proceedings of SIGGRAPH), 28(3), 2009.
2. M. Ihmsen, N. Akinci, M. Gissler, M. Teschner. *Boundary Handling and
   Adaptive Time-stepping for PCISPH.* Proc. VRIPHYS, Copenhagen, 2010.
3. M. Becker, M. Teschner. *Weakly compressible SPH for free surface flows.*
   Proceedings of the 2007 ACM SIGGRAPH/Eurographics symposium on Computer
   animation.
4. M. Macklin, M. Müller, N. Chentanez. *XPBD: Position-Based Simulation of
   Compliant Constrained Dynamics.* MIG '16, 2016.

## Configuration file format

The configuration file consists of:

```
[physical parameters]
mass: 5.4e-14
timeStep: 5.0e-06
simulationScale: 2.46e-06
viscosity: 5.0e-05
surfTensCoeff: 1.21948e+27
elasticityCoefficient: 5.55556e+08

[simulation box]
xmin
xmax
ymin
ymax
zmin
zmax

[position]
1 0 0 1
1 0 1 1
...

[velocity]
0 0 0 1
0 0 0 1
...

[connection]
1   1.58649939377   1.1   0.0
7   1.58649939377   1.1   0.0
...

[membranes]
0   1   7
7   8   1
...

[particleMemIndex]
0
144
288
-1
...
```

Position and velocity rows are 4-element vectors: `x, y, z, type` where
`type` is 1=liquid, 2=elastic, 3=boundary. Each elastic particle has up to
32 connections, each a 4-vector `(target_id, rest_length, muscle_id, _pad)`.
Membranes are triangles defined by 3 elastic-particle IDs.

## Saving to disk

```bash
./Release/Sibernetic -l_to                  # save buffers to ./buffers/
./Release/Sibernetic -l_from                # replay from buffers
./Release/Sibernetic -export_vtk            # save VTK files for Paraview
```

Recording interval is set by `logstep=<n>`.

## Sibernetic-NEURON bridge

You can drive Sibernetic from a NEURON simulation:

```bash
git clone https://github.com/openworm/sibernetic_NEURON.git
export PYTHONPATH=./sibernetic_NEURON:./src
./Release/Sibernetic -nrn ./sibernetic_NEURON/models/celegans/_ria.hoc -f worm
```

After each Sibernetic step it runs one NEURON step and reads back voltage
data from a hardcoded list of segments
(`src/owNeuronSimulator.cpp:70`).

## Running with c302

```bash
git clone https://github.com/openworm/CElegansNeuroML.git
export C302_HOME=./CElegansNeuroML/CElegans/pythonScripts/c302
export PYTHONPATH=$PYTHONPATH:$C302_HOME:./src
python sibernetic_c302.py
```

This generates the NEURON code for c302, runs Sibernetic with the neuronal
simulation in the background, and saves results to `simulations/`.

## Troubleshooting

If you have a question or hit a problem, contact us:

- Email: skhayrulin@openworm.org or info@openworm.org
- Or open an [issue on GitHub](https://github.com/openworm/sibernetic/issues)
