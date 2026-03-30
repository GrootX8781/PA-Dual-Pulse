# Polarization-Agile Dual-Pulse Radar Simulation

This repository contains the MATLAB implementation for the simulation framework used in the paper on polarization-agile dual-pulse radar against deceptive jamming. The code focuses on the polarization-domain signal model, pulse-wise estimation, task-oriented transmit polarization design, and the corresponding robustness evaluations.

## Overview

The implemented framework follows a two-pulse logic.

In pulse 1, the radar uses an estimation-oriented polarization trajectory to probe the scene and obtain the information required for subsequent design, including target-scattering-related quantities and jammer-polarization-related quantities.

In pulse 2, the transmit polarization is redesigned in a task-oriented manner using the information extracted from pulse 1. The receive side then applies constrained processing to improve the output SINR and the overall anti-jamming performance.

The simulations in this repository cover the following aspects:

- the main SINR and gain evaluation versus JSR,
- the pulse-1 jammer-polarization estimation error behavior,
- several robustness studies under model mismatch or implementation perturbations,
- and one visualization example for the polarization-trajectory and response behavior.

## Folder Structure

    FUNCTIONS/
    SCRIPTS/

- `FUNCTIONS/` contains the supporting functions used by the simulation scripts.
- `SCRIPTS/` contains the main experiment scripts and plotting utilities.

## Main Scripts

### `run_main.m`

Main entry for the SINR and gain simulation.

This script generates the core performance data of the proposed method and the comparison methods under different JSR settings. It is the main script for reproducing the SINR-gain results reported in the paper.

### `pjerr_test.m`

Pulse-1 jammer-polarization estimation error experiment.

This script evaluates the estimation behavior in the first pulse and produces the corresponding error results.

### `pjest_plot.m`

Visualization example.

This script provides a representative visualization example of the estimated or designed polarization behavior, corresponding to the illustrative figure in the paper in the style of Fig. 5.

### `robust_phase.m`

Robustness experiment under phase perturbation or phase mismatch.

### `robust_pj.m`

Robustness experiment with jammer-polarization-related perturbation or estimation deviation.

### `robust_S.m`

Robustness experiment with target-scattering-related mismatch.

### `config.m`

Parameter configuration script used by the simulation scripts.

### `configplot.m`

Common plotting configuration script.

## Data File

### `FUNCTIONS/psm.mat`

This file stores the scattering-related data required by the simulations. Please keep it in the expected location relative to the scripts.

## Quick Start

1. Open MATLAB.
2. Add the project folders to the MATLAB path.
3. Enter the `SCRIPTS/` folder.
4. Run the script corresponding to the experiment of interest.

A typical order is:

    run_main
    pjerr_test
    pjest_plot
    robust_phase
    robust_pj
    robust_S

If your local MATLAB setup does not automatically include subfolders in the search path, run the following command from the project root before executing the scripts:

    addpath(genpath(pwd))

## Implementation Logic

The code is organized according to the main signal-processing flow of the paper.

First, the waveform and polarization trajectory are generated. Then, the target echo and jamming-contaminated response are constructed under the adopted fully polarimetric signal model. Based on the pulse-1 observation, the relevant quantities for the second pulse are estimated. These estimates are then used to form the pulse-2 task-oriented polarization design and the corresponding constrained receive processing. Finally, the scripts evaluate the output SINR, gain, estimation error, and robustness metrics under different simulation settings.

In this implementation, the experiment scripts in `SCRIPTS/` are kept relatively explicit so that each experiment can be run independently, while the common computational components are placed in `FUNCTIONS/`.

## Notes

- The scripts are intended to be run independently according to the target experiment.
- Figures and numerical results may depend on Monte Carlo settings and random initialization if random seeds are not fixed in the current script version.
- For consistent reproduction, it is recommended to keep the folder structure unchanged.
- The plotting style is controlled mainly through `configplot.m`.

## MATLAB Environment

The code is written in MATLAB. A standard recent MATLAB version is recommended.

If your local environment requires additional toolboxes for numerical linear algebra or plotting, please make sure they are installed before running the scripts.

## Citation

If you use this code in your research, please cite the corresponding paper.