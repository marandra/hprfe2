# hprfe2

Python package for managing the full lifecycle of **HPR-FE²**
(High-Performance Reduced Finite Element²) multiscale simulations on
HPC clusters: from parameterized case generation through offline
training to result collection, field reconstruction, and packaging.

Companion code to the
[MultiscaleROMApplication](https://github.com/marandra/MultiscaleROMApplication)
module of [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos).

## What is HPR-FE²?

HPR-FE² is a hyper-reduced order modelling technique for two-scale
finite-element analysis of materials with heterogeneous microstructure.
The full FE² method couples a macroscopic finite-element problem to an
independent microscopic FE problem (a Representative Volume Element,
RVE) at every macroscopic integration point — accurate but
computationally intractable for industrial cases.

HPR-FE² approximates the RVE response in two complementary ways:

1. **Reduced order modelling (ROM)** — the RVE displacement field is
   projected onto a low-dimensional basis obtained from a Singular
   Value Decomposition (SVD) of strain snapshots collected during an
   offline training stage.
2. **Hyper-reduced integration** — the internal force and tangent
   stiffness integrals over the RVE are evaluated only at a small set
   of optimally selected sampling points, avoiding a full loop over
   Gauss points.

The combined effect is up to 4 orders of magnitude of speed-up over a
standard FE² solver while keeping the accuracy loss below 1%, enabling
industrial-scale multiscale simulations.

## What `hprfe2` does

`hprfe2` is the orchestration / data-engineering layer around the
solver. It is responsible for the parts of an HPR-FE² study that live
*outside* the FEM solver itself:

- **Case generation** — read a parameter file and generate a set of
  ready-to-run simulation case directories (mesh, materials, BCs,
  solver settings, job scripts).
- **Offline training** — drive snapshot acquisition runs, then build
  the reduced basis via Randomized SVD or Block SVD
  (NumPy / scikit-learn), and write the trained ROM to disk in a
  format the solver can consume.
- **HPC job orchestration** — submit Slurm jobs for case batches,
  monitor status, retry failed jobs, and collect logs.
- **Result collection and packaging** — gather per-case outputs into
  a single dataset (HDF5 + metadata) with full traceability of input
  parameters.
- **Field reconstruction** — reconstruct full-field RVE responses from
  the reduced solution at user-selected macroscopic points for
  visualization and post-processing.

## Why a separate package

Keeping the orchestration logic in Python — outside the C++ solver —
makes the workflow:

- **Reproducible**: every study is a parameter file + a git revision of
  the package.
- **Testable**: pure-Python parts of the pipeline are covered by
  `pytest`.
- **Extensible**: alternative training algorithms or HPC backends can
  be plugged in without touching the solver.
- **Documented**: API reference and tutorials are built with Sphinx.

## Installation

```bash
git clone https://github.com/marandra/hprfe2
cd hprfe2
pip install -e .
```

Requires Python 3.8+, NumPy, SciPy, scikit-learn, h5py.
A working install of Kratos Multiphysics with
[MultiscaleROMApplication](https://github.com/marandra/MultiscaleROMApplication)
is required to actually run simulations; the package itself only
depends on Python.

## Quick start

```bash
# Generate a batch of cases from a parameter file
hprfe2 cases generate study.yaml --out cases/

# Run the offline training stage on the generated snapshots
hprfe2 train --cases cases/ --rom-out rom/

# Submit the parametric study to Slurm and collect results
hprfe2 run --cases cases/ --rom rom/ --backend slurm

# Reconstruct full RVE fields at selected macro points
hprfe2 reconstruct --results results.h5 --points points.csv
```

See `docs/` for the full CLI reference, configuration schema, and a
worked example.

## Architecture

```
hprfe2/
├── cli.py            # Click-based CLI entry points
├── cases/            # Case generation from parameter files
├── train/            # Offline SVD / Block-SVD training
├── orchestrate/      # HPC backends (local, Slurm), job state machine
├── io/               # HDF5 result aggregation, metadata
├── reconstruct/      # Field reconstruction from reduced bases
└── tests/            # pytest test suite
```

## Status

This repository preserves the version used in the published HPR-FE²
studies (see *Citing* below). The public package is maintained for
reproducibility of the published work.

## Citing

If you use `hprfe2` or HPR-FE² in academic work, please cite:

> M. Raschi, O. Lloberas, A. Huespe, J. Oliver.
> *High performance technique for multiscale finite element (HPR-FE²):
> towards industrial multiscale FE software.*
> Computer Methods in Applied Mechanics and Engineering, 2021.

## License

See `LICENSE`.
