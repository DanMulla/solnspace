# Solution spaces of shoulder muscle activity

Code, models, inputs and outputs for the simulation study:

> **Musculoskeletal simulations reveal how kinematic degrees of freedom and effort shape shoulder neuromuscular control**


## Contents

- [Overview](#overview)
- [Replicating the study](#replicating-the-study)
- [Requirements](#requirements)
- [Repository layout](#repository-layout)
- [The linear problem](#the-linear-problem)
- [Running on macOS (Apple silicon)](#running-on-macos-apple-silicon)

## Overview

The study characterises the space of muscle activation patterns that can satisfy
the mechanical demands of a static shoulder exertion, and asks how that space is shaped by
(i) the number of kinematic degrees of freedom (DOF) the model is required to balance and
(ii) an effort-based motor control criterion.

The pipeline has four stages:

1. **Inverse dynamics** (Python + OpenSim) — pose the shoulder model, apply a hand load,
   and extract the net joint moments, muscle moment arms, maximum isometric forces, and the
   glenohumeral (GH) joint reaction force geometry.
2. **Solution space sampling** (MATLAB) — assemble the mechanical equilibrium requirements
   as an underdetermined linear problem, find the minimum- and maximum-effort solutions by
   optimisation, then sample the feasible set with a Riemannian Hamiltonian Monte Carlo
   (RHMC) algorithm. Repeated for four models of increasing DOF (3, 5, 8, 11).
3. **Summary statistics** (R) — activation ranges and quantiles by muscle group; distance
   between "low effort" solutions and the minimum-effort solution.
4. **Figures** (R) — activation distributions (box-and-whisker) and individual activation
   patterns (flipped parallel-coordinates plots).

Each stage reads the previous stage's files from disk, so stages can be re-run independently.

## Requirements

| Stage | Software |
| --- | --- |
| Inverse dynamics | Python 3 with the **OpenSim** Python API, plus `numpy`, `pandas`, `matplotlib`, `openpyxl` |
| Sampling | **MATLAB** with the Optimization Toolbox (`fmincon`) and **PolytopeSamplerMatlab** |
| Analysis / figures | **R** with `conflicted`, `gdata`, `gridExtra`, `tidyverse`, `viridis`, `rstudioapi` |

- OpenSim Python API setup: <https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53085346/Scripting+in+Python>
- OpenSim API reference: <https://simtk.org/api_docs/opensim/api_docs/>
- PolytopeSamplerMatlab: <https://github.com/ConstrainedSampler/PolytopeSamplerMatlab>

  The sampler is GPL-3.0 licensed and is **not** redistributed here. Clone it into the
  repository root before running Stage 2:

  ```bash
  git clone https://github.com/ConstrainedSampler/PolytopeSamplerMatlab.git
  ```

  `feasiblesolns.m` adds it to the MATLAB path via `addpath(genpath(...))` on the
  repository folder, so no further configuration is needed.

On Apple silicon see [Running on macOS](#running-on-macos-apple-silicon) before the first run.

A working environment can be created with:

```bash
conda create -n solnspace -c opensim-org -c conda-forge python=3.11 "opensim=4.4" numpy pandas matplotlib openpyxl
```

## Repository layout

```
solnspace/
├── inversedynamics.py          Stage 1 driver: loops over all tasks
├── utils.py                    Stage 1 helper functions (OpenSim wrappers, .mot writers,
│                               scapular coordinate system, GH JRF geometry, plotting)
├── feasiblesolns.m             Stage 2: builds the linear problem, optimises, samples
├── build_polytopesampler_macos.sh   Builds the sampler's MEX files on Apple silicon
├── dataanalysis.R              Stage 3: summary statistics across all tasks
├── dataviz.R                   Stage 4: per-task publication figures
├── models/shoulder/            OpenSim models + the model-reduction script
├── inputs/                     Task definitions and OpenSim setup/template files
├── outputs/task01 … task13/    Created by the pipeline, one folder per task
│                               (not published; regenerate by running Stages 1 and 2)
├── statistics/                 Stage 3 summary .csv files (the values reported in the paper)
└── PolytopeSamplerMatlab/      RHMC sampler (not redistributed; clone separately)
```

## The linear problem

For a rigid linked-segment system in static equilibrium with `m` balanced kinematic DOF and
`n` actuators, the mechanical equilibrium requirement is

```
A_eq · a = M          A_eq = [ R·F | +T | −T ]
```

where

- `R` is the `m × 42` matrix of muscle moment arms at the balanced DOF,
- `F` is the `42 × 42` diagonal matrix of maximum isometric (tendon) forces at the simulated
  pose, so `R·F` is the maximum moment each muscle can generate at each DOF,
- `T` is the `m × m` diagonal matrix of reserve actuator strengths,
- `M` is the `m × 1` vector of net joint moments from inverse dynamics,
- `a` is the `n × 1` vector of actuator activations, subject to `0 ≤ a ≤ 1`.

Reserve strengths are **10 Nm at each SC and AC DOF** and **1 Nm at each GH, HU and RU DOF**.

**Glenohumeral stability** is imposed as eight linear inequalities approximating the glenoid
dislocation boundary (Halder et al. 2001; Dickerson et al. 2007; Blache et al. 2016):

The joint reaction force is the resultant of every non-contact force on the humerus — the
muscle pull, which is linear in the activations, plus the contribution of gravity and the
external hand load, which is constant. Writing the cone as `U · F_total ≤ 0` and separating
the two:

```
U · (R'·a + f_load)  ≤  0        →        (U·R')·a  ≤  −U·f_load
└────┬────┘                                └──┬──┘      └───┬───┘
  F_total                                   Aineq          bineq
```

so `Aineq = U·R'` and `bineq = −U·f_load`. `R = [Fx Fz Fy]` holds the per-actuator GH joint
reaction force components in the modified scapular frame (`Fx` is the glenoid-normal /
compression axis; reserve actuators contribute zero), and `U` holds the eight
shear-to-compression threshold triples, starting superiorly and proceeding clockwise when
looking into the glenoid.

**Effort** is reported two ways, both normalised by the actuator count `n`:

```
J₁ = Σ aᵢ  / n        (linear)
J₂ = Σ aᵢ² / n        (non-linear)
```

**Two sampling passes** are run per model:

1. *Activation landscape* — approximately uniform sampling over the whole feasible polytope.
   `nsamples = 10^4` is requested.
2. *Effort landscape* — the linear effort range between the minimum- and maximum-effort
   solutions is divided into 100 equal intervals, and `nsamples/100 = 100` solutions are
   requested from each by adding the pair of inequalities `J₁ ∈ [eₖ, eₖ₊₁]`. Intervals are
   defined on the **linear** effort term because the sampler requires linear constraints,
   while the interval end points are the linear effort of the solutions that minimise and
   maximise the **quadratic** term.

RHMC does not return exactly the requested number of samples. In this repository the
activation-landscape files contain roughly 8,400–10,600 solutions and the effort-landscape
files roughly 17,100–17,900.

## Replicating the study

Start to finish, on a clean machine. Expect roughly two hours, most of it Stage 2.

1. **Install the environment.**

   ```bash
   conda create -n solnspace -c opensim-org -c conda-forge python=3.11 "opensim=4.4" numpy pandas matplotlib openpyxl
   conda activate solnspace
   ```

   You also need MATLAB with the Optimization Toolbox, and R with `conflicted`, `gdata`,
   `gridExtra`, `tidyverse` and `viridis`. On Apple silicon, build the sampler's MEX files
   first — see [Running on macOS](#running-on-macos-apple-silicon).

2. **Check the sampler loads.** In MATLAB, from the repository root:

   ```matlab
   addpath(genpath(pwd)); compile_solver(0); compile_solver(4)
   ```

   This should return silently. If it raises *"No C++ mex compiler is available"*, see the
   macOS section.

3. **Stage 1 — inverse dynamics** (about 3 minutes for all 13 tasks):

   ```bash
   python inversedynamics.py
   ```

   Writes five `.csv` files per task into `outputs/<task>/`. Sanity check: `<task>-Fiso.csv`
   must be identical across all 13 tasks, because every task uses the same pose.

4. **Stage 2 — sampling** (about 2 minutes per task per model). In MATLAB, once per task:

   ```matlab
   feasiblesolns      % prompts for a task name, e.g. task03
   ```

   Repeat for `task01` through `task13`. Each run loops over the four models internally and
   writes to `outputs/<task>/solutions/`.

5. **Stage 3 — statistics.** Open `dataanalysis.R` in RStudio and source it. Writes the
   summary `.csv` files to `statistics/`.

6. **Stage 4 — figures.** Open `dataviz.R` in RStudio and source it once. It loops over all
   13 tasks and the three effort percentiles (0.01, 0.05, 0.10) internally, writing the figures
   into each task's `outputs/<task>/solutions/` folder. The distribution and sample plots are
   written only at the 0.05 percentile; the percentile comparison plot is written for each.

**Reproducibility.** `feasiblesolns.m` sets `rngseed = 2024` and derives an explicit seed for
every sampler call, so re-running on the same machine reproduces the sample files exactly.
Results may *not* bit-reproducible across platforms: the x86 and arm64 builds of the sampler's
`PackedChol` take different code paths (AVX2 versus a portable scalar fallback), so trajectories
diverge from rounding even under an identical seed. Distributions agree; individual solutions do
not. Change `rngseed` to draw an independent sample set.

## Running on macOS (Apple silicon)

The pipeline was developed on Windows. Two changes were needed to make
`feasiblesolns.m` run on an Apple silicon Mac. Both apply to the PolytopeSamplerMatlab
clone rather than to this repository, so they must be made after cloning it; neither
affects the Windows or Linux build.

**1. `PackedChol` MEX binaries.** PolytopeSamplerMatlab compiles its `PackedChol` solver
on first use, via `compile_solver.m` calling MATLAB's `mex`. On macOS this fails with

```
Error using compile_solver (line 29)
No C++ mex compiler is available.
```

even when Xcode Command Line Tools and a working `clang++` are installed. The reason is
that MATLAB's C++ mex configuration (`$MATLABROOT/bin/maca64/mexopts/clang++_maca64.xml`)
resolves `XCODE_AGREED_VERSION` by reading Xcode's own preference keys
(`IDEXcodeVersionForAgreedToGMLicense`), which only exist with a full Xcode install whose
licence has been accepted. Command Line Tools alone can never satisfy that check, so
MATLAB reports no compiler and `compile_solver` aborts. A second problem follows: the
`PolytopeSamplerMatlab/bin/` folder ships prebuilt binaries for Windows (`.mexw64`),
Linux (`.mexa64`) and Intel macOS (`.mexmaci64`), but MATLAB R2024b on Apple silicon is a
native arm64 application and looks for `.mexmaca64`, so it finds nothing to load and
falls through to compiling.

Either fix works:

- **Build the MEX files directly** with the Command Line Tools compiler, using the flags
  from MATLAB's own mexopts XML:

  ```bash
  bash build_polytopesampler_macos.sh
  ```

  This writes `PackedChol{0,1,4}arm.mexmaca64` into `PolytopeSamplerMatlab/bin/`. Once
  they exist, `compile_solver` finds them (`exist(name) == 3`), probes them, and never
  invokes `mex`. Pass a different MATLAB path as the first argument if not using
  `/Applications/MATLAB_R2024b.app`. The binaries are committed, so this only needs
  re-running after a MATLAB upgrade or if they are deleted.

- **Or install full Xcode** from the App Store, run `sudo xcodebuild -license accept`,
  then `mex -setup C++` in MATLAB, and let `compile_solver` build them itself. This is
  the route MathWorks supports; the script above just avoids the download.

**2. `immintrin.h` include guard.** `PolytopeSamplerMatlab/code/solver/PackedCSparse/FloatArray.h`
included `<immintrin.h>` unconditionally. That header is x86-only and fails to compile on
arm64. The x86 intrinsics in that file are already confined to an `#ifdef __AVX2__` block
with a complete portable fallback, so the include is now guarded:

```cpp
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#include <immintrin.h>
#endif
```

`_M_X64` is defined by MSVC, so Windows x64 builds still get the header and the AVX2 path
exactly as before. Current upstream PolytopeSamplerMatlab carries an equivalent guard, so a recent clone
may already contain it; check the file before editing.

Verify the result in MATLAB:

```matlab
addpath(genpath('/path/to/solnspace'))
compile_solver(0); compile_solver(4);   % should return silently
```

Note that `MexSolver.solverName()` appends `arm` to the solver name when
`sysctl -n machdep.cpu.brand_string` contains `Apple`, which is why the files are named
`PackedChol0arm` rather than `PackedChol0`. There is no MEX-free fallback:
`Solver.m` always constructs a `MexSolver`, and `MatlabSolver` uses one internally for its
exact solves, so these binaries are required, not merely an optimisation.

The remaining macOS obstacles were Windows path separators and hard-coded absolute paths,
now corrected; those affected Stage 1 (`inversedynamics.py`), not the sampling stage.