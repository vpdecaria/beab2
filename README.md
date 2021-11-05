# MOOSE-IMEX-12
This is an implementation of an adaptive stepsize and or Implicit/Explicit (IMEX) timestepping scheme for the incompressible Navier-Stokes equations. The method is called MOOSE-IMEX-12. MOOSE stands for Multiple Order One Solve Embedded, because there is only one linear solve each time step, but two solutions of different orders of accuracy are produced. The '12' is because approximations are of orders one and two.

This code was used to produce the IMEX results in 
"An embedded variable step IMEX scheme for the incompressible Navier–Stokes equations", Victor DeCaria and Michael Schneier, CMAMEm 2021. It is based on a code I used earlier to produce numerical results for the linearly-implicit version in 
"A Time-Accurate, Adaptive Discretization for Fluid Flow Problems", Victor Decaria, William Layton and Haiyun Zhao, IJNAM, 2020. If you find this code useful in your research, it would be much appreciated if you cited the appropriate paper.

We use FEniCS for the finite element discretization.

This is an embedded timestepping method based on backward Euler/Adams Bashforth 2 (BE-AB2), which is where the name of the repository comes from.

BEAB2.py is the main script. With FEniCS installed, it is run with

python3 BEAB2.py

This will run the code with default parameters.

To use MOOSE-IMEX-12, run with the vo (for variable order flag)

python3 BEAB2.py --vo 12

Specify the starting stepsize with -k (default is 1e-6) and the tolerance with -t (default is 1e-3).
Give it your own problem input deck with -p, or --problem. The default is taylor_green_problem. For example,

python3 BEAB2.py --vo 12 -k 0.1 -t 1e-5 --problem taylor_green_problem

This runs MOOSE-IMEX-12 with a starting stepsize of 0.1, a tolerance of 1e-5, and specifies the Taylor Green vortex problem. taylor_green_problem is a python file, but DO NOT put .py at the end.

# Directories

## tests

Contains bash scripts that automate tests BEAB2.py. Should be run from top directory, using

bash tests/nameofscript.sh

## problems

Contains modules that define the domain, boundary conditions, and body forces for different tests.


