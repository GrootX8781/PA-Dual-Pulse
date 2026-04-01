# PA-Dual-Pulse

MATLAB simulation code for the paper:

**Polarization-Agile Trajectory Optimization for Deception Jamming Suppression in Dual-Pulse Cooperative Polarimetric Radar**

## Overview

This repository contains the MATLAB code used to reproduce the main simulation results in the paper.

The implemented framework follows a dual-pulse cooperative logic:
- Pulse 1 is used for estimation-oriented probing.
- Pulse 2 is used for task-oriented redesign and suppression.

## Repository Structure

- `SCRIPTS/`  
  Main scripts for reproducing the paper results.
- `FUNCTIONS/`  
  Supporting functions used by the main scripts.
- `README.md`  
  Repository description.

## Requirements

- MATLAB
- CVX

## Usage

1. Open MATLAB and set the project root as the current folder.
2. Run a main script in `SCRIPTS/`, such as `main_fig3.m`, `main_fig4.m`, etc.
3. For runtime evaluation, run `runtime.m`.

## Configuration

- `config.m` is the main configuration file.
- `configplot.m` is used for plotting configuration.
- Before running the scripts, please manually define the target PSM in `config.m`.
- The PSM should be given as a `2x2` matrix and stored in the variable `S`.

## Notes

- The `main_fig*.m` files are the main simulation scripts.
- The files in `FUNCTIONS/` are called internally by the scripts and usually do not need to be modified.
- Please make sure the MATLAB path includes both `SCRIPTS/` and `FUNCTIONS/`.

## Citation

If you use this repository in your research, please cite the corresponding paper.