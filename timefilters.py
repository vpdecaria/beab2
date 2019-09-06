"""
Time Filters Module - contain helper functions for the vsvo backward euler time stepping scheme.
"""


import numpy as np
#import timefiltercoefficient as tfc


def nu(filterNum, tau):
	"""
	Computes the various filter parameters
	tau_n = tau[4]
	.
	.
	.
	tau_nm4 = tau[0]
	"""
	switcher = {
		1:tfc.nu1(tau),
		2:tfc.nu2(tau),
		3:tfc.nu3(tau),
		4:tfc.nu4(tau),
	}
	return switcher.get(filterNum, "No such parameter")
	
def kappa(filterNum,tau,U):
	switcher = {
		1:tfc.kappa1(tau,U),
		2:tfc.kappa2(tau,U),
		3:tfc.kappa3(tau,U),
		4:tfc.kappa4(tau,U),
	}
	return switcher.get(filterNum, "No such parameter")

def estimator(filterNum, K, U):
	"""
	U is the vector of the previous 5 solutions and most recent one
	U[5] = unplus1
	...
	U[0] = unminus4
	This function returns the filterNum(th) estimator that is to be subtracted from the  previously calculated solutionsolution
	"""
	tau = K[1:5]/K[0:4]
	switcher = {
		1: nu(1,tau)/2.*kappa(1,tau,U),
		2: nu(2,tau)/2.*kappa(2,tau,U),
		3: nu(3,tau)/2.*kappa(3,tau,U),
		4: nu(4,tau)/2.*kappa(4,tau,U),
	}
	return switcher.get(filterNum, "no such filter")

def filter(filterNum,K,U):
	"""
	Returns the filtered solution, and the error estimator for the unfiltered solution
	filterNum is between 1-4.
	K is a vector of the five most recent step sizes
	U is a vector of the six most recently calculated solutions
	"""
	errorEst = estimator(filterNum, K, U)
	return [U[5] - errorEst, errorEst]


def estDeltaT(k_n,errEst,tol,order):
	"""
	k_n - The most recent timestep taken
	errEst - Some norm of the error estimator
	order - the order of the method that ErrEst is supposed to be estimating
	the local truncation error in.
	tol - Some prescribed tolerance.
	"""
	return  k_n*(tol/(errEst))**(1./(order+1));
	

def pickSolutionMaxK(Est,tolerance,k_n,safetyFactor,orders):
	"""
	Y - Vector of solutions [Y(1),Y(2),...], where the 
	index of Y represents the order of the method used to compute the
	solution.
	Est - Array of estimators for Y, in the same length as the matrix Y.
	tolerance - Some desired tolerance prescribed by user
	k_n - the most recent step size
	safetyFactor - Allows a margin of error for the estimated step size
	Pick between zero and one
	OUTPUT
	knPlus1 - The next step size
	j - The index of Y(j) with the largest allowable step size based on 
	the estimated error and the tolerance.
	
	This has been modified from the backward euler code so that it can 
	account for other methods
	"""
	n = len(Est)
	K_nPlus1 = np.zeros(n)
	for j, est in enumerate(Est):
		K_nPlus1[j] = safetyFactor*estDeltaT(k_n,est,tolerance,orders[j]);
	return [np.amax(K_nPlus1),np.argmax(K_nPlus1)]
