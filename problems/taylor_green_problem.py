from __future__ import print_function
from fenics import *
import numpy as np

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
T = 4          # final time
penalty = 1e-6

# Create mesh and define function spaces
N = 50  #number of nodes
mesh = RectangleMesh(Point(0,0),Point(2*np.pi,2*np.pi), N, N,"left/right")  #Even number ensures no degenerate nodes.

V = VectorElement("Lagrange", mesh.ufl_cell(), 3)
Q = FiniteElement("Lagrange", mesh.ufl_cell(), 2)
#V = VectorElement("CR", mesh.ufl_cell(), 1)
#Q = FiniteElement("DG", mesh.ufl_cell(), 0)
TH = V * Q

W = FunctionSpace(mesh, TH,constrained_domain=periodic_bc)
dw = TrialFunction(W)

QuadSpace = FunctionSpace(mesh, FiniteElement("Lagrange", mesh.ufl_cell(), 2),constrained_domain=periodic_bc)

bcs = []

# Define expressions used in variational forms
n   = FacetNormal(mesh)

class pressure_point_boundary(SubDomain):
	def inside(self,x,on_boundary):
		return (on_boundary and near(x[0],0) and near(x[1],0))
bc1 = DirichletBC(W.sub(1), (0), pressure_point_boundary(),"pointwise")

#TAYLOR GREEN

def FF(t):
	return np.exp(-2*nu*t)
def FF_prime(t):
	return -2*nu*np.exp(-2*nu*t)

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
def get_p_exact_ema(t):
	return Expression('-0.25*F*F*(cos(2*x[0])+cos(2*x[1])) - 0.5*(pow(F*cos(x[0])*sin(x[1]),2) + pow(-F*cos(x[1])*sin(x[0]),2))',t = t, F=FF(t),degree =1 )
def get_bcs(t):
	return [bc1]

print(W.dim())

