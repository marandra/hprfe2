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
- **Reduced-basis construction** — build the reduced bases from
  solution snapshots via SVD (NumPy / scikit-learn) and select the
  hyper-reduced integration-point (cubature) sets.
- **HPC job scripts** — generate Slurm array launcher scripts for the
  sampling and validation case batches, ready to submit on the cluster.
- **Packaging** — pack the bases and datasets into HDF5 for the Kratos
  solver to consume.
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

Requires Python 3.6+ and `docopt`, `numpy`, `h5py`, `scikit-learn`,
`meshio`.
A working install of Kratos Multiphysics with
[MultiscaleROMApplication](https://github.com/marandra/MultiscaleROMApplication)
is required to actually run simulations; the package itself only
depends on Python.

## Quick start

```bash
# Create a default configuration file in the project root
hprfe2 config

# Deploy the sampling case structure + Slurm launcher scripts from a template case
hprfe2 deploy

# (run the sampling jobs on the cluster with the generated Slurm array scripts)

# Build reduced bases, integration-point sets and reconstruction datasets
hprfe2 generate

# Generate single-integration-point multiscale validation cases
hprfe2 validate
```

See `docs/` for the full CLI reference, configuration schema, and a
worked example.

## Architecture

```
hprfe2/
├── hprfe2            # docopt CLI entry point (config / deploy / generate / validate)
├── common.py         # configuration, project paths, HDF5 resource I/O
├── sampling.py       # deploy sampling cases + Slurm array launcher scripts
├── bases.py          # reduced-basis construction via SVD (NumPy / scikit-learn)
├── roc.py            # integration-point (cubature) selection for hyper-reduction
├── pack.py           # pack bases / datasets into HDF5 for the solver
├── reconstruction.py # full-field RVE reconstruction (displacement, stress, damage)
└── multiscale.py     # single-IP multiscale validation cases
tests/                # pytest suite
docs/                 # Sphinx documentation
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
