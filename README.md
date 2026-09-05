# ConformationalFluctuation

Molecular dynamics simulation and analysis code for investigating how conformational fluctuations of intrinsically disordered proteins (IDPs) are coupled to intermolecular interactions and biomolecular condensate properties.

This repository contains the simulation setup, data-driven coarse-graining workflow, and analysis scripts used to connect **single-chain conformational ensembles** with **multichain condensate structure and dynamics**, with the FUS low-complexity domain (FUS-LC) as the primary model system.

## Overview

The workflow combines atomistic molecular dynamics, residue-level coarse-grained simulations, and a maximum-entropy-based correction of coarse-grained interactions.

In broad terms:

1. **All-atom simulations** provide reference conformational and contact statistics.
2. A conventional **HPS residue-level coarse-grained model** provides the baseline coarse-grained description.
3. A **maximum-entropy / data-driven refinement procedure** iteratively adjusts residue-type interaction terms so that the coarse-grained model reproduces selected contact-probability statistics from the reference ensemble.
4. **Single-chain and multichain simulations** are then analyzed to determine how different levels of conformational fluctuation affect molecular dimensions, intermolecular contacts, dynamics, network organization, and condensate structure.

The repository is organized primarily as research code supporting this workflow rather than as a general-purpose software package.

## Repository Structure

```text
ConformationalFluctuation/
├── all_atom/      # Analysis of atomistic MD trajectories
├── hps/           # Baseline HPS coarse-grained simulations and analysis
├── max_entr/      # Maximum-entropy/data-driven CG refinement and single-chain analysis
├── mul_chn/       # Multichain/condensate simulations and analysis
└── coeff_tsf/     # Analysis and comparison of interaction coefficients
```

### `all_atom/`

Analysis scripts for atomistic simulations and reference conformational ensembles.

Representative analyses include:

- radius of gyration
- RMSD
- molecular shape
- residue-residue distances
- contact probabilities
- reduced/contact-based representations
- conformational-angle analysis

These results provide reference structural statistics for evaluating and refining coarse-grained models.

### `hps/`

Scripts associated with the conventional **hydrophobicity scale (HPS)** residue-level coarse-grained model.

The `hps/init/` directory contains scripts for generating LAMMPS data, parameter, and input files for FUS-LC, including conventional simulations and replica-exchange molecular dynamics (REMD).

Additional directories contain analyses of:

- radius of gyration
- shape
- residue-residue distances
- contact probabilities
- trajectory processing

### `max_entr/`

Implementation and analysis associated with the **maximum-entropy/data-driven refinement** of the coarse-grained model.

The refinement procedure uses contact-probability observables obtained from the reference ensemble and iteratively updates residue-type interaction coefficients. The covariance matrix of the selected observables is used to determine the coefficient correction at each iteration.

The main refinement workflow is located under:

```text
max_entr/init/type_mom_prob_cov_hps/
```

Important scripts include:

- `max_entr_init.py` — initialization of the refinement procedure
- `max_entr_file.py` — preparation of simulation/input files
- `max_entr_pc.py` — calculation of contact-probability statistics
- `max_entr_coeff.py` — update of data-driven interaction coefficients
- `max_entr.sh` — iterative workflow driver

The remaining directories contain analyses of the refined single-chain models, including:

- radius of gyration and conformational fluctuations
- molecular shape
- distance distributions
- contact probabilities
- contact lifetimes
- autocorrelation functions
- diffusion-related quantities
- Fourier-transform-based analyses
- potential-energy functions

### `mul_chn/`

Simulation analysis for **multichain systems and biomolecular condensates**.

The analyses include:

- intermolecular contacts
- contact lifetimes
- molecular diffusion
- chain clustering
- intermolecular network connectivity
- chain orientation and orientational order
- residue distributions
- density and spatial distributions
- chain dimensions
- end-to-end distances
- local structure across the condensed phase and interface
- stress-related quantities
- trajectory processing

This part of the repository is used to connect differences in single-chain conformational fluctuation with collective properties of the condensed phase.

### `coeff_tsf/`

Scripts for analyzing the data-driven residue-type interaction coefficients and their relationships with underlying HPS parameters such as residue size, hydrophobicity, and charge.

## Requirements

The code was developed for molecular-dynamics analysis in a research/HPC environment. The exact requirements vary between scripts, but the main dependencies include:

- Python 3
- NumPy
- SciPy
- Matplotlib
- scikit-learn
- LAMMPS

Some workflows also use shell scripts and VMD-based trajectory processing.

Several Python scripts import locally developed helper functions, for example:

```python
from pyw import pyw
from xtc import xtc_rd
from tot_len import tot_len
```

and currently contain machine-specific paths such as:

```python
sys.path.append(...)
```

Therefore, before running the scripts on another machine, these paths should be changed and the corresponding helper functions should be made available in the Python environment.

## Basic Usage

Clone the repository:

```bash
git clone https://github.com/CiboFeng/ConformationalFluctuation.git
cd ConformationalFluctuation
```

Most scripts are intended to be run from their corresponding analysis directories because file names, relative paths, and system-specific parameters are defined directly in the scripts.

For scripts using command-line arguments, help information can generally be inspected with:

```bash
python script_name.py -h
```

For example, the maximum-entropy coefficient update script accepts input coefficient files, trajectory-statistics directories, sequence information, reference contact probabilities, and output file names through command-line arguments.

Simulation parameters, trajectory locations, system names, and analysis settings should be checked inside each script before execution.

## Data-Driven Coarse-Graining

The central refinement strategy in this repository introduces additional residue-type-dependent interaction terms to a baseline coarse-grained model.

For a set of observables \( \mathbf{P} \), the interaction coefficients are iteratively corrected using the discrepancy between the simulated and reference contact-probability statistics together with their covariance matrix. Schematically,

\[
\Delta \boldsymbol{\alpha}
\propto
\mathbf{C}^{-1}
\left(
\mathbf{P}_{\mathrm{sim}}
-
\mathbf{P}_{\mathrm{ref}}
\right),
\]

where \( \mathbf{C} \) is the covariance matrix of the selected observables and \( \boldsymbol{\alpha} \) denotes the data-driven interaction coefficients.

The purpose of this procedure is to construct coarse-grained models that retain the efficiency required for large-scale simulations while reproducing selected conformational statistics of the higher-resolution reference ensemble.

## Notes on Reproducibility

This repository contains the scripts used in the original research workflow. Some scripts therefore include:

- system-specific file names
- absolute paths from the original workstation/HPC environment
- predefined protein sequences and simulation parameters
- intermediate plotting and testing scripts
- generated figures alongside source scripts

Users wishing to reproduce or extend the calculations should inspect these settings before execution.

Large raw molecular-dynamics trajectories are not stored in this repository.

## Citation

If you use this code or the associated methodology in published work, please cite the corresponding publication describing the conformational-fluctuation models and simulations.

The full citation will be added here when available.

## Contact

For questions about the repository, please open a GitHub issue or contact the repository owner through GitHub.
