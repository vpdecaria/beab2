"""
This is a not fully implicit Euler code for solving the Navier-Stokes equations.
The user can optionally use an unconditionally stable linearly implicit method,
or a conditionally stable implicit-explicit method.

The code is also adaptive by default, but adaptivity can be turned off with the
--constant flag.

The user can supply their own problem via the -p or --problem tag. Omit the '.py'
at the end of the file. These problem files are like an input deck and must
contain certain information, but they are python files.
"""

from __future__ import print_function
from fenics import *
import numpy as np

from types import ModuleType
import lib.timefilters as tf
import time as pytimer
import argparse
from subprocess import check_output


comm = MPI.comm_world
mpiRank = MPI.rank(comm)

#Define a print function to be used for parallel runs
def mpi_print(output):
	if mpiRank==0:
		print(output)
	else:
		pass

#------------------------------------ CONSOLE INPUT ------------------------------------------------

parser = argparse.ArgumentParser(description='An adaptive, linearly implicit/implicit-explicit \
	euler code.')

parser.add_argument(\
'-f','--filters', help ="Run the test using up to the filter number specified \
	by this argument. An input of 0 is unmodified Backward Euler. The default selection is\
	4, which allows for the calculation of the estimator of the fourth order method.",\
	 type = int,default = 1)
parser.add_argument(\
'-o','--output', help ="The file name for the text file containing the errors\
	 with respect to delta t. Default file name is date-time.txt",type =str, default = \
	 "errors/temp/" + str(check_output(['date','+%Y-%m-%d_%R:%S'])).strip()+".txt")
parser.add_argument(\
'-t','--tolerance', help ="The tolerance used to adapt the step size. \
	Right now, it is just based on error committed per time step, but can be made more \
	sophisticated later.",type =np.float64, default = 1e-3)
parser.add_argument(\
'-r','--ratio', help ="The maximum step size ratio allowed, where step size \
	ratio is (step size at time n divided by step size at time n minus 1) The default value is 2.",\
	type =np.float64, default = 2)
parser.add_argument(\
'--constant', help ="Enabling this will disable adaptivity altogether.",\
	action="store_true")
parser.add_argument(\
'-k','--startingStepSize', help ="The initial step size taken.",\
	type =np.float64, default = 0.000001)
parser.add_argument(\
'--forcefilter', help ="Forces the exclusive use of the filter \
	passed to argument f.",action="store_true")
parser.add_argument(\
'-p','--problem', help ="The name of a user created python module that contains\
	all the information necessary to run a specific test problem, including mesh, \
	boundary conditions, body forces, exact solutions, etc. There is a certain syntax that I need\
					to specify eventually. ",type =str, default = 'taylor_green_problem')

parser.add_argument(\
'--paraview', help ="Output name for pvd and vtu files",type =str, default = "pvd/tmp")
parser.add_argument(\
'--parfreq', help ="Frequency with respect to delta t to take paraview snapshots.",\
	 type = int,default = 1000000)

parser.add_argument(\
'--error', help ="Evaluate Error norms.",action="store_true")

parser.add_argument(\
'--vo', help ="Which orders to use, such as 2, 23, 3, 34, 234", type = int,default = 1)


parser.add_argument(\
'--nopics', help ="Don't write paraview output.",action="store_true")

parser.add_argument(\
'--semi', help ="Use the more stable, but costlier, linearly implicit version.",action="store_true")

parser.add_argument(\
'--extrap', help ="Which IMEX scheme to use. Choices are ab2 and fe.",default = "ab2")
parser.add_argument(\
'--coldstart', help ="In the absense of a nice solution to initialize with, this will take the\
	 first three steps with BEFE.",action="store_true")
parser.add_argument(\
'-s','--solver', help ="specify which solver to use. Defauts to the solve command",\
	type=str,default="solve")
parser.add_argument(\
'--stepseq', help ="This feature is intended to force some possibly pathological stepsize sequence.\
					 Possible options include...",type =str, default = "false")

#----------------------------------- PARSING THE CONSOLE INPUT -----------------------------------												
args = parser.parse_args()						
solver_type=args.solver 	
stepsequence=args.stepseq
coldstart=args.coldstart
use_semi_implicit= args.semi
imex_extrapolation=args.extrap
calculate_errors = args.error
writePVD = args.nopics
orders_to_use = args.vo
#This if statement will trigger if any value is entered besides false or 0
#Step sequence has format
#firstratio_secondratio....._finalratio For example, to double every step, 
#stepsequence = '2'. 
#To double then halve
#stepsequence = '2-0.5'.      More complicated sequences can be encoded this way.
if(stepsequence=='false'):
	stepsequence=False
if(stepsequence):
	if(orders_to_use == 12):
		error("CUSTOM STEPSIZE SEQUENCE NOT INTENDED TO BE USED WITH VARIABLE ORDER.")
	else:	
		stepsequence=stepsequence.strip().split('-')
		#Convert to floats
		stepsequence=[float(i) for i in stepsequence]

#Paraview arguments
pviewOut = args.paraview
writePVD = True
if(args.nopics):
	writePVD = False


paraview_frequency = args.parfreq

exec('from problems.'+str(args.problem) +' import *')
#from problem import *

print("Using filter number "+ str(orders_to_use))

constantStepSize = args.constant
if(constantStepSize):
	print("Constant step size of " +str(args.startingStepSize)+ ". Adaptivity is turned off.")
	
print("Printing to file " + args.output)

tolerance = args.tolerance

maxRatio = args.ratio

forcefilter = args.forcefilter
#-------------------------------- END PARSING CONSOLE INPUT --------------------------------------

#
#Initialize quantities related to adaptivity
EstVector = np.ones(2)*1e54
tempSolutions = ["Intermediate","Filtered","Solution","Values"]
safetyFactor = 0.9
numOfFailures = 0
minStepSize = 1e-12

errorfName = args.output
errorfile = open(errorfName, 'w')
output = "T final =" + str(T) +", Filters Used "+ str(orders_to_use)  +'\n'
errorfile.write(output)

dt = np.float64(args.startingStepSize)
dt_n=dt
dt_nM1=dt_n
k   = dt
K = np.ones(4)*dt

#vector of times. Earlier times at start of array

#Determine how many steps the method is (one step or multistep)
if (orders_to_use == 12 or orders_to_use == 2):
	if not constantStepSize:
		total_num_steps = 3
	elif constantStepSize:
		total_num_steps = 2
if (orders_to_use == 1):
	if constantStepSize:
		if imex_extrapolation == 'fe':
			total_num_steps = 1
		elif imex_extrapolation == 'ab2':
			total_num_steps = 2
	if not constantStepSize:
		if imex_extrapolation == 'ab2':
			total_num_steps = 2

Ts = np.array([-j*dt for j in range(total_num_steps,-1,-1)])

#----------------------- DECLARE FINITE ELEMENT FUNCTIONS -----------------------------------------

v,q  = TestFunctions(W)

# Define functions
w   = TrialFunction(W)  # current solution

w_ = Function(W)
(u_,p_) = split(w_)

w_n  = Function(W)  # solution from previous converged step
w_extrap  = Function(W)  # solution from previous converged step

# Split mixed functions
u,  p  = split(w)

w_n = Function(W)
w_nM1=Function(W)
w_nM2=Function(W)

p_temp = Function(W.sub(1).collapse())

null_vector=Function(W.sub(1).collapse())
assign(null_vector,interpolate(Constant(1.0),W.sub(1).collapse()))

#------------------------ CREATE INITIAL CONDITIONS -----------------------------------------------
assign(w_n.sub(0),interpolate(get_u_exact(0),W.sub(0).collapse()))
assign(w_n.sub(1),interpolate(get_p_exact(0),W.sub(1).collapse()))

assign(w_nM1.sub(0),interpolate(get_u_exact(-dt),W.sub(0).collapse()))
assign(w_nM1.sub(1),interpolate(get_p_exact(-dt),W.sub(1).collapse()))

assign(w_nM2.sub(0),interpolate(get_u_exact(-2*dt),W.sub(0).collapse()))
assign(w_nM2.sub(1),interpolate(get_p_exact(-2*dt),W.sub(1).collapse()))

u_n,p_n=split(w_n)
u_nM1=split(w_nM1)[0]
u_nM2=split(w_nM2)[0]
#-----------------------  END CREATE INITIAL CONDITONS --------------------------------------------

#initialize temp solutions
error_1st_order=Function(W)
error_2nd_order=Function(W)

y_2 = Function(W)

#TIME STEPPING
output = "\n Time, StepSize, FilterUsed, NormU, NormP ,NormUExact, NormPExact, Uerror, Perror"
errorfile.write(output)

#Compute measure of domain
measure = assemble(Constant(1.0)*dx(mesh))

# Define convection term

def convect(u,v,w):
	return dot(dot(u, nabla_grad(v)), w)*dx

# Define skew-symmeterized forn
def b(u,v,w):
	return convect(u,v,w)+0.5*div(u)*dot(v,w)*dx

#####  INITIALIZE ERROR QUANTITIES ####
L2_error = 0                          #
                                      #
l2L2_error = 0                        #
l2L2 = 0                              #
l2L2_error_pressure = 0               #
l2L2_pressure = 0                     #
l2H1 = 0                              #
l2H1_error = 0                        #
                                      #
exact_vel_norm_sq = 0                 #
exact_pres_norm_sq = 0                #
p_L2_error = 0                        #
total_error_time = 0                  #
#######################################

#LAST INITIALIZATIONS BEFORE TIMESTEPPING LOOP
t = 0
tOld = 0
step_counter = 0

if writePVD:
	file = File(pviewOut + ".pvd")


###########################    CONFIGURE KRYLOV SOLVER    ##########################################

solver = KrylovSolver('bicgstab', "amg")
#solver = KrylovSolver('cg', "amg")

#solver = KrylovSolver('gmres', "amg")

#I'm not sure what PETScKrylovSolver was. Maybe related to null space stuff?
#solver = PETScKrylovSolver("gmres","amg")

#solver.parameters["relative_tolerance"] = 1e-8
solver.parameters["error_on_nonconvergence"] = False
solver.parameters["maximum_iterations"] = 100

null_space = VectorSpaceBasis([null_vector.vector()])

loop_timer = pytimer.time()
###########################################################
#                                                         #
#        LL         OOOOO      OOOOO     PPPPP            #
#        LL        O     O    O     O    P    P           #
#        LL        O     O    O     O    P    P           #
#        LL        O     O    O     O    PPPPP            #
#        LL        O     O    O     O    P                #
#        LL        O     O    O     O    P                #
#        LLLLLLLLL  OOOOO      OOOOO     P                #
#                                                         #
###########################################################

while (tOld < T-1e-15):
	if(mpiRank == 0):
		print("Current dt - ", dt)

	t = tOld + dt
	k=Constant(dt)

	Ts[total_num_steps] = t
	
	# Update body force
	f = get_f(t)
	
	################   COMPUTE EXTRAPOLATION ######################
	omega_n = dt/dt_n
	omega_nM1 = dt_n/dt_nM1
	extrap_order = 4

	if(imex_extrapolation=="fe" or (coldstart and step_counter<=2)):
		w_extrap.vector()[:] = w_n.vector().get_local()
		print("USING FE EXTRAPOLATION")
	elif(imex_extrapolation=="ab2"):
		w_extrap.vector()[:] = (1+omega_n)*w_n.vector().get_local()\
			 - omega_n*w_nM1.vector().get_local()
		print("USING AB2 EXTRAPOLATION")
	#w_extrap.vector()[:]=((1+omega_n)*(1+omega_nM1*(1+omega_n)))/(1+omega_n)\
	#			*w_n.vector().get_local()\
	#			-omega_n*(1+omega_nM1*(1+omega_n))\
	#			*w_nM1.vector().get_local()\
	#			+(omega_nM1**2*omega_n*(1+omega_n))/(1+omega_nM1)\
	#			*w_nM2.vector().get_local()
	#w_.vector()[:] = w_n[total_num_steps].vector()


	###############   END COMPUTE EXTRAPOLATION   ###################
	

	##############   DEFINE FUNCTIONAL   #########3########       
	if(not use_semi_implicit):
		F = dot((u-u_n)/k,v)*dx                           \
			+b(w_extrap.sub(0),w_extrap.sub(0),v)         \
			+nu*inner(nabla_grad(u), nabla_grad(v))*dx    \
			-p*div(v)*dx                                  \
			+ div(u)*q*dx                                 \
			-dot(f,v)*dx  
	else:#---------USE the semi implicit method-----------
		F = dot((u-u_n)/k,v)*dx                           \
			+  b(w_extrap.sub(0),u,v)                     \
			+nu*inner(nabla_grad(u), nabla_grad(v))*dx    \
			-p*div(v)*dx                                  \
			+ div(u)*q*dx                                 \
			-dot(f,v)*dx 
	#######################################################

	## ASSEMBLE MATRIX IF NECESSARY ##
	if(dt != dt_n or step_counter==0 or use_semi_implicit):
		#The 'dt!=dt_n' condition checks if the formulation changed, and requires a new matrix.
		#Using semi_implicit requires a new matrix every time.
		A = assemble(lhs(F))

	bcs = get_bcs(t)
	b_rhs = assemble(rhs(F))            ## APPLY BOUNDARY CONDITIONS ##
	for bc in bcs:
		bc.apply(A)
		bc.apply(b_rhs)
	
	
	###############    SOLVE THE LINEAR SYSTEM    #################
	
	#KRYLOV SOLVER
	if(solver_type=="krylov"):
		if(dt != dt_n or step_counter==0):
			solver.set_operator(A)
			#the null space thing makes it not work?
			#as_backend_type(A).set_nullspace(null_space)
			print('HI')

		num_krylov_iterations = solver.solve(w_.vector(),b_rhs)

		print("Krylov solver converged in ", num_krylov_iterations)
	

	elif(solver_type=="solve"):
		#DEFAULT SOLVE (Whatever 'solve' happens to be in this version of fenics)
		#In single thread mode, seems to be umfpack.
		solve(A,w_.vector(),b_rhs)

	elif(solver_type== "lu"):
		#LU SOLVER - Don't use this if using semi-implicit or if dt often
		if(dt != dt_n or step_counter==0):
			#Configure LU SOLVER
			print("Performing LU factorization.")
			lu_solver = LUSolver(A)
			print("Finished setting up LU solver.")
		lu_solver.solve(w_.vector() ,b_rhs)
	else:
		error("User argument is not a valid linear solver.")

	#Make pressure mean zero
	p_temp =  w_.split(True)[1]
	p_temp.vector()[:] = p_temp.vector().get_local() - \
		assemble(p_temp/measure*dx)*np.ones_like(p_temp.vector().get_local())
	assign(w_.sub(1), p_temp)
	
	print(EstVector)

	#------------------------------- APPLY TIME FILTERS  ------------------------------------------

	if (not constantStepSize or orders_to_use==2 or orders_to_use == 12):
		#Form error estimator
		error_1st_order.vector()[:] = omega_n/(1.+2*omega_n)*(w_.vector().get_local()\
			-w_extrap.vector().get_local())
		if(orders_to_use == 1 or orders_to_use == 12):
			EstVector[0] = norm(error_1st_order.sub(0),'L2',mesh)

		if(orders_to_use == 2 or orders_to_use == 12):
			#Use the second order filter
			y_2.vector()[:] = w_.vector().get_local() - error_1st_order.vector().get_local()

			if(not constantStepSize):
				error_2nd_order.vector()[:] = \
					(omega_n*omega_nM1*(1+omega_n))/(1+2*omega_n + \
					omega_nM1*(1+4*omega_n+3*omega_n**2))\
					*(y_2.vector().get_local() \
					-((1+omega_n)*(1+omega_nM1*(1+omega_n)))/(1+omega_nM1)\
					*w_n.vector().get_local()\
					+omega_n*(1+omega_nM1*(1+omega_n))\
					*w_nM1.vector().get_local()\
					-(omega_nM1**2*omega_n*(1+omega_n))/(1+omega_nM1)\
					*w_nM2.vector().get_local())
				EstVector[1]=norm(error_2nd_order.sub(0),'L2',mesh)


	#----------------------------  CHOOSE NEXT STEPSIZE  ------------------------------------------

	#If using a predetermined stepsize ratio sequence...
	if(stepsequence):

		knp1=dt*stepsequence[step_counter % len(stepsequence)]

	elif not constantStepSize and not (coldstart and step_counter<=2):
		#The LARGE_NUM below is a hacky way to exclude solutions that don't satisfy
		#the tolerance from consideration for picking the next step size.
	
		LARGE_NUM=1e40
		TempEstVector = EstVector + ~(EstVector< tolerance)*LARGE_NUM
		[knp1,J] = tf.pickSolutionMaxK(TempEstVector,tolerance,dt,safetyFactor,[1,2])
		knp1 = np.max([np.max([np.min([knp1,maxRatio*dt]),dt/2]),minStepSize])
		#Force them all to end at the same time
		knp1 = np.min([knp1,T-t])
	elif (coldstart and step_counter<=2):
		J=0
		knp1=dt
	#-----------------  EITHER ADVANCE SOLUTION OR REJECT AND RECALCULATE ------------------

	if(stepsequence or constantStepSize or EstVector[J] < tolerance\
	or np.abs(knp1 - minStepSize) < 1e-10 or (coldstart and step_counter<=2)):

		#        PICK WHICH SOLUTION TO ADVANCE IN TIME WITH
		#--------------------------------------------------------------
		if(orders_to_use == 1):      #IF ORDER IS CONSTANT            #
			J=0                      #THEN SET J TO CHOOSE            #
		elif(orders_to_use == 2):    #THAT ORDER HERE.                #
			J=1                      #ORDER = J+1                     #
        #--------------------------------------------------------------
		if(J==0):                    #J = 0. NO NEED TO DO            #
			pass                     #ANYTHING. OTHERWISE, SET EQUAL  #
		elif(J==1):                  #TO FILTERED SOLUTION y_2.       #
			assign(w_.sub(0), y_2.sub(0))                             #
		#--------------------------------------------------------------
			
		#-------------- PRINT OUT SOLUTION INFORMATION AT THIS TIME LEVEL ------------------------
		#------------------ CALCULATE ERRORS IF USER SPECIFIED TO  -------------------------------
		if(mpiRank == 0):
			print("Using order ",J+1)
			print('At Time t = %.6f' % (t))
		if(calculate_errors):
			error_timer_begin=pytimer.time()
			u_exact = get_u_exact(t)
			u_exact_interpolated = interpolate(u_exact,W.sub(0).collapse())
			L2_error = errornorm(u_exact_interpolated, w_.sub(0),degree_rise = 0)
			
			p_exact = get_p_exact(t)
			p_exact_interpolated = interpolate(p_exact,W.sub(1).collapse())
			p_L2_error = errornorm(p_exact_interpolated, w_.sub(1),degree_rise = 0)

			#_______________Update Temporal Error________________________
			exact_vel_norm_sq =norm(u_exact,'L2',mesh)**2               #
			l2L2 += exact_vel_norm_sq*dt                    #VELOCITY   #
			l2L2_error += L2_error**2*dt                                #
			#-----------------------------------------------------------#
			exact_pres_norm_sq = norm(p_exact,'L2',mesh)**2             #
			l2L2_pressure += exact_pres_norm_sq*dt          #PRESSURE   #
			l2L2_error_pressure += p_L2_error**2*dt                     #
			#___________________________________________________________|

			if(mpiRank == 0):
				print('t = %.2f: L2_error = %.3g' % (t, L2_error))
				print('t = %.2f: L2_error_pressure = %.3g' % (t, p_L2_error))
			total_error_time += pytimer.time()-error_timer_begin


		normU = norm(w_.sub(0),'L2',mesh)
		normP = norm(w_.sub(1),'L2',mesh)
		if(mpiRank == 0):
			print('t = %.2f: Norm U = %.3g' % (t, normU))
			print('t = %.2f: Norm P = %.3g' % (t, normP))
		

		output = "\n" + str(t) + "," + str(dt) + "," + str(J) + "," +str(normU)+ "," + \
		str(normP)+ ","+str(exact_vel_norm_sq**0.5)+ ","+ str(exact_pres_norm_sq**0.5) \
		+ "," +str(L2_error)+ ","+str(p_L2_error)
		errorfile.write(output)
		
		#____________________________________________________
		#            REINITIALIZE VARIABLES AND              |
		#       STEPSIZE TO GO BACK TO BEGINNING OF LOOP     |
		#____________________________________________________|
		w_nM2.vector()[:]=w_nM1.vector().get_local()
		w_nM1.vector()[:]=w_n.vector().get_local()
		w_n.vector()[:] = w_.vector().get_local()

		tOld += dt
		
		for j in range(total_num_steps):
			Ts[j] = Ts[j+1]
		if not constantStepSize:
			
			dt_nM1=dt_n
			dt_n=dt
			dt = knp1

			k = knp1

			
		if(step_counter%paraview_frequency ==0 and writePVD):
			file << (w_.sub(0),t)

		step_counter+=1

	#---------  ELSE STATEMENT BELOW IS USED IF SOLUTION IS REJECTED ---------------
	#---------  PICK NEW STEPSIZE, GO BACK TO BEGINNING OF LOOP --------------------
	#---------  DO NOT ADVANCE TIME ------------------------------------------------
	else:
		print("Failed: Halving Time Step")
		#Update Estimates
		saftey_factor_failed = 0.7
		[knp1,junk] = tf.pickSolutionMaxK(EstVector,tolerance,dt,saftey_factor_failed,[1,2])
		#dt = np.max([dt/2.,knp1])
		dt = knp1
		numOfFailures += 1
		#Force all simulations to end at the final time T
		knp1 = np.min([knp1,T-t])

##################################################################################
#                                                                                #
#    EEEEEEEE  N     N   DDDDD          L         OOOOO      OOOOO     PPPPP     #                  
#    E         NN    N   D    D         L        O     O    O     O    P    P    #  
#    E         N N   N   D     D        L        O     O    O     O    P    P    # 
#    EEEEEEEE  N  N  N   D     D        L        O     O    O     O    PPPPP     #           
#    E         N   N N   D     D        L        O     O    O     O    P         #
#    E         N    NN   D    D         L        O     O    O     O    P         # 
#    EEEEEEEE  N     N   DDDDD          LLLLLLLL  OOOOO      OOOOO     P         #  
#                                                                                #                     
################################################################################## 

elapsed_time = pytimer.time()-loop_timer
print("Main loop took ",elapsed_time ," seconds.")
print("Spent ",total_error_time ," seconds calculating errors.")

#------------------------------   CALCULATE ERRORS  -----------------------------------------------
#-------------------------------   WRITE TO FILES     ---------------------------------------------

errorfile.write("\nEND\ntolerance, l2L2 error, l2H1 error, l2L2 Pressure error, \
	Number Of Rejections, Starting Step Size, Elapsed Time")
if calculate_errors:
	relative_l2L2_error = np.sqrt(l2L2_error)/np.sqrt(l2L2)
	relative_l2L2_error_pressure = np.sqrt(l2L2_error_pressure)/np.sqrt(l2L2_pressure)
	relative_l2H1_error = 'NOT CALCULATED'
else:
	relative_l2L2_error = 'Not applicable'
	relative_l2L2_error_pressure = 'Not applicable'
	relative_l2H1_error = 'Not applicable'
print("l2L2 error:",relative_l2L2_error)
print("l2L2P error:",relative_l2L2_error_pressure)

#Write errors to file
output = "\n"  +str(tolerance)            + "," + str(relative_l2L2_error)          + "," \
			   +str(relative_l2H1_error)  + "," + str(relative_l2L2_error_pressure) + "," \
			   +str(numOfFailures)        + "," + str(args.startingStepSize)        + "," \
			   +str(elapsed_time)

errorfile.write(output)
errorfile.close()

