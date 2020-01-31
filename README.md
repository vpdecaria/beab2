# MOOSE-IMEX-12
This is an implementation of an adaptive stepsize and or Implicit/Explicit (IMEX) timestepping scheme for the incompressible Navier-Stokes equations. The method is called MOOSE-IMEX-12. MOOSE stands for Multiple Order One Solve Embedded, because there is only one linear solve each time step, but two solutions of different orders of accuracy are produced. The '12' is because approximations of order one and two are created each timestep.

We use FEniCS for the finite element discretization.

This is an embedded timestepping method based on backward Euler/Adams Bashforth 2 (BE-AB2), which is where the name of the repository comes from.

BEAB2.py is the main script. With FEniCS installed, it is run with

python3 BEAB2.py

This will run the code with default parameters.

# Directories

## tests

Contains bash scripts that automate tests BEAB2.py. Should be run from top directory, using

bash tests/nameofscript.sh

## problems

Contains modules that define the domain, boundary conditions, and body forces for different tests.


