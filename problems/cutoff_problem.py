"""
This is a fully coupled, fully implicity Backward Euler Code to solve NSE in a periodic box with an exact solution.

An exact solution to the Taylor Green Vortex in 2D is given by
u = (F(t)cos(x)sin(y), -F(t)cos(y)sin(x))
p = -(1/4)F(t)^2 (cos(2x) + cos(2y))

No boundary conditions because periodic domain.
"""

from __future__ import print_function
from fenics import *
import numpy as np
import matplotlib.pyplot as plt

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):

	# Left boundary is "target domain" G
	def inside(self, x, on_boundary):
		# return True if on left or bottom boundary AND NOT on one of the two corners (0, 2pi) and (2pi, 0)
		return bool((near(x[0], 0) or near(x[1], 0)) and 
                (not ((near(x[0], 0) and near(x[1], 2*np.pi)) or 
                        (near(x[0], 2*np.pi) and near(x[1], 0)))) and on_boundary)

	def map(self, x, y):
		if near(x[0], 2*np.pi) and near(x[1], 2*np.pi):
			y[0] = x[0] - 2*np.pi
			y[1] = x[1] - 2*np.pi
		elif near(x[0], 2*np.pi):
			y[0] = x[0] - 2*np.pi
			y[1] = x[1]
		elif near(x[1], 2*np.pi):
			y[0] = x[0]
			y[1] = x[1] - 2*np.pi
		else:
			y[0] = -1000
			y[1] = -1000
            
periodic_bc = PeriodicBoundary()

# Define boundary conditions
def boundary(x, on_boundary):
	return on_boundary

#Basic Problem info
nu = 1.0/1.0             # kinematic viscosity
Re = 1.0 / nu
T = 45.           # final time
penalty = 1e-6

# Create mesh and define function spaces
N = 50  #number of nodes
mesh = RectangleMesh(Point(0,0),Point(2*np.pi,2*np.pi), N, N,"left/right")  #Even number ensures no degenerate nodes.

V = VectorElement("Lagrange", mesh.ufl_cell(), 3)
Q = FiniteElement("Lagrange", mesh.ufl_cell(), 2)
TH = V * Q

W = FunctionSpace(mesh, TH,constrained_domain=periodic_bc)
dw = TrialFunction(W)

bcs = []

# Define expressions used in variational forms
n   = FacetNormal(mesh)
scale = 2
translation = 1-1./scale
translation = 0

def cutoff(t):
	if t <=0:
		return 0
	else:
		return np.exp(-1./((10*t)**10))
	
def cutoffPrime(t):
	if t <=0:
		return 0
	else:
		return np.exp(-1./((10*t)**10))*100/(10*t)**11
	

	
def flatQuickStepUp(t):
	return cutoff(t-9)
def flatQuickStepDown(t):
	return 1- cutoff(t-9)
def flatQuickStepUpPrime(t):
	return cutoffPrime(t-9)
def flatQuickStepDownPrime(t):
	return - cutoffPrime(t-9)

#def flatQuickStepUp(t):
#	return mollifierStep(t-9)
#def flatQuickStepDown(t):
##	return 1- mollifierStep(t-9)
def smoothSteps(t):
	if t <= 0:
		return 0
	elif t % 20 <= 10 - 1e-14:
		return flatQuickStepUp(t%10)
	else:
		return flatQuickStepDown(t%10)
def smoothStepsPrime(t):
	if t <= 0:
		return 0
	elif t % 20 <= 10 - 1e-14:
		return flatQuickStepUpPrime(t%10)
	else:
		return flatQuickStepDownPrime(t%10)
	


#t = np.linspace(-1e-1,40,1000000)


def discSteps(t):
	if t <= 0:
		return 0
	elif t % 20 <= 10:
		return 0
	else:
		return 1

def FF(t):
	return smoothSteps(t)
#def FF_prime(t):
#	return (FF(t+1e-5) -FF(t-1e-5))/2e-5
def FF_prime(t):
	return smoothStepsPrime(t)
#y = [2*nu*FF(T) + FF_prime(T) for T in t]


#t = np.linspace(-1e-1,40,1000000)
#y = [FF(T) for T in t]
#y = [smoothSteps(T) for T in t]
#plt.plot(t,y)
#plt.show()

#This is the time dependent portion of the rhs
def RHS_time(t):
	return 2*nu*FF(t) + FF_prime(t)

#Create Functions that returns the rhs at every time
def get_f(t):
	return Expression(('F*cos(x[0])*sin(x[1])','-F*cos(x[1])*sin(x[0])'),t=t,F = RHS_time(t),degree = 2)
#Create functions that return updated exact solutions at every time
def get_u_exact(t):
	return Expression(('F*cos(x[0])*sin(x[1])','-F*cos(x[1])*sin(x[0])'),t=t,F=FF(t),degree = 2)
def get_p_exact(t):
	return Expression('-0.25*F*F*(cos(2*x[0])+cos(2*x[1]))',t = t, F=FF(t),degree =1 )

def get_bcs(t):
	return []
"""
def get_f(t):
	if(0<= t%1 <= 1-2./scale):
		return Expression(('0','0'),degree = 2)
	elif(1-2./scale < t%1 <= translation):
		return whatever
	else:
		return Expression(('(exp)*cos(x[0])*sin(x[1])','-exp(-2 *Nu *t)*sin(x[0])*cos(x[1])'),t=0.0,Nu=nu,degree = 2)
"""
