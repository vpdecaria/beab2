#Meant to be run from one directory level up.
let mpi_threads=4

for vo in 12 # 2 1
do
	for stepsize in '0.023100616016427104' # '1.5517241379310345' '0.36' '0.20361990950226244' '0.1323529411764706' '0.0774526678141136' '0.043478260869565216'  #'1e-7'
	do
		mpirun -np $mpi_threads python3 BEAB2.py\
			   -p cutoff_problem \
			   --error \
			   -o errors/convergenceTestWRTtoleranceHarder/order-$vo-tol-$stepsize-constantstep.txt \
			   --vo 2 \
			   -k $stepsize\
			   --constant
	done

done