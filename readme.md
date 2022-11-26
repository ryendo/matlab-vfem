# Rigorous shape derivative of eigenvalues estimation (alpha version)

Ryoki ENDO

First version: 2022/08/24

## Running environment

- The codes run on latest MATLAB environment along with INTLAB toolbox (http://www.ti3.tu-harburg.de/rump/intlab/).
- Since Octave does not supply **ldl** function for indefinite symmetric matrix, the rigorous evaluation is not available for Octave.

## Configuration

Edit **my_intlab_mode_config.m** to configure the computing environment.

```MATLAB
%The path of INTLAB toolbox and initialization.
addpath('Intlab_V12')

%The path of the library of verified eigenvalue estimation for matrix.
addpath('verified_eig_estimation')

%The path of the codes for switch between verified computing and approximate computing.
addpath('mode_swith_interface')

%The path of the library of FEMs on triangular domains
addpath('functions')
startintlab;

global INTERVAL_MODE; %The value can be 0 or 1. See readme.md

%INTERVAL_MODE=1; for rigorous computing based on interval arithmetic.
%INTERVAL_MODE=0; for approximate computing with rounding error inside.
INTERVAL_MODE=1;

```

Remeber to run "my_intlab_mode_config" before other codes.

## How to use it?

For the purpose of rigorous computation, please run the code in MATLAB along with the INTLAB toolbox.

The code has two running mode: approximate evaluation and verified evaluation.

To swith between each other, please set the value of **INTERVAL_MODE**.

- INTERVAL_MODE=0: approximate computation mode.
- INTERVAL_MODE=1: verified computation mode. It takes longer time that approximate computation mode. So be careful with the mesh size, which should not be too small. 
