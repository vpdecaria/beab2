#Meant to be run from one directory level up.
let mpi_threads=1

for vo in 1 12 2
do
	for tolerance in '1' '1e-1' '1e-2' '1e-3' '1e-4' '1e-5' '1e-6' '1e-7'
	do
		mpirun -np $mpi_threads python3 BEAB2.py\
			   -p cylinder_with_exact_solution \
			   --error \
			   -o errors/convergenceTestWRTtolerance/order-$vo-tol-$tolerance.txt \
			   --vo $vo \
			   -t $tolerance
	done

done