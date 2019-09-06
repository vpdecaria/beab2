#Mean to be run from one directory level up.
let mpi_threads=4

for vo in 1
do
	for extrap in 'fe' 'ab2'
	do
		for dt in '0.2' '0.1' '0.05' '0.025' '0.0125' '0.00625'
		do
			mpirun -np $mpi_threads python3 BEAB2.py --constant -p taylor_green_problemP2P1disc -k $dt --error -o errors/convergenceTestConstantStepP2P1disc/order-$vo-extrap-$extrap-dt-$dt.txt --vo $vo --solver lu --extrap $extrap
		done
	done

done


for vo in 2
do
	for extrap in 'ab2'
	do
		for dt in '0.2' '0.1' '0.05' '0.025' '0.0125' '0.00625'
		do
			mpirun -np $mpi_threads python3 BEAB2.py --constant -p taylor_green_problemP2P1disc -k $dt --error -o errors/convergenceTestConstantStepP2P1disc/order-$vo-extrap-$extrap-dt-$dt.txt --vo $vo --solver lu --extrap $extrap
		done
	done

done

