# beab2
Backward Euler/Adams Bashforth 2. An adaptive, implicit/explicit, timestepping method for the Navier-Stokes equation using FEniCS.

BEAB2.py is the main script. With FEniCS installed, it is run with

python3 BEAB2.py

This will run the code with default parameters.

# Directories

## tests

Contains bash scripts that automate tests BEAB2.py. Should be run from top directory, using

bash tests/nameofscript.sh

## problems

Contains modules that define the domain for different test.


