%The path of INTLAB toolbox and initialization.
addpath('/Users/raphael/Documents/calc/MATLAB/lib/Intlab_V12')

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
INTERVAL_MODE=0;
