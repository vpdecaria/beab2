This directory contains scripts (not data) to make the plots found in
"An embedded variable step IMEX scheme for the incompressible Navier-Stokes equations."
- Victor DeCaria, Michael Schneier.

This readme explains how to make specific plots founds in the paper.
From the directory that contains this readme.txt, copy and 
paste the commands below a title into the terminal to make these plots,
assuming that you ran the tests to create the data.

"Nonadaptive velocity/pressure error"
cd constant_step/
bash make_plots.sh

"Adaptive velocity/pressure error"
cd convergence_tol_cylinder_harder/
bash make_plots.sh

"Adaptive vs nonadaptive velocity/pressure error"
cd convergence_adapt_vs_constant/
bash make_plots.sh

"Adaptive vs nonadaptive velocity norms with same number of stokes solves"
norm_adapt_vs_nonadapt.py