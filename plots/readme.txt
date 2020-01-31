This directory contains scripts (not data) to make the plots found in
"An embedded variable step IMEX scheme for the incompressible Navier-Stokes equations."
- Victor DeCaria, Michael Schneier.

This readme explains how to make specific plots founds in the paper.
From the directory that contains this readme.txt, copy and 
paste the commands below a title into the terminal to make these plots.


This assumes that you ran the tests yourself.
If you do not run the tests yourself to create the data, our copy is in 
errors/data_for_paper. From the directory with this readme file, run the following 
two commands in a terminal. Then you should be able to generate the plots.

cp ../errors/data_for_paper/convergenceTestWRTtoleranceHarder/*.txt ../errors/convergenceTestWRTtoleranceHarder/
cp ../errors/data_for_paper/convergenceTestConstantStep/*.txt ../errors/convergenceTestConstantStep/

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