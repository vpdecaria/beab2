from __future__ import print_function
from fenics import *
import numpy as np
from mshr import *
import matplotlib.pyplot as plt

#Basic Problem info
nu = 1.0             # kinematic viscosity
Re = 1.0 / nu
T = 1         # final time
#Mesh
#mesh = Mesh('cyliner.finer.mesh.xml.gz')

N=400
bigCircle = Circle(Point(0,0),1,N)
penalty = 1e-6

mesh = generate_mesh(bigCircle,10)
#mesh = generate_mesh(domain,20)
plot(mesh)
plt.show()


# Define function spaces

V = VectorElement("Lagrange", mesh.ufl_cell(), 2)
Q = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = V * Q

W = FunctionSpace(mesh, TH)
dw = TrialFunction(W)

#Define BC
class noslip_boundary(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary
		
class pressure_point_boundary(SubDomain):
	def inside(self,x,on_boundary):
		return (on_boundary and near(x[0],1) and near(x[1],0))

bc0 = DirichletBC(W.sub(0), (0,0), noslip_boundary())

bc1 = DirichletBC(W.sub(1), (0), pressure_point_boundary(),"pointwise")

bc = [bc0,bc1]

#Create Functions that returns the rhs at every time
def get_f(t):
	return Expression(('pi*x[1]*(-1 + pow(x[0],2) + pow(x[1],2))*cos(pi*t) - sin(pi*t)*(-1 + 8*nu*x[1] + x[0]*pow(-1 + pow(x[0],2) + pow(x[1],2),2)*sin(pi*t))'\
	,				  '-pi*x[0]*(-1 + pow(x[0],2) + pow(x[1],2))*cos(pi*t) - sin(pi*t)*(-1 - 8*nu*x[0] + x[1]*pow(-1 + pow(x[0],2) + pow(x[1],2),2)*sin(pi*t))'),t=t,nu=nu, degree = 2)
	#return Expression(('pow(0.41,-2)*1*6*x[1]*(0.41-x[1])','0'),degree = 2)
#Create functions that return dummy exact solutions at every time
def get_u_exact(t):
	return Expression(('-sin(pi*t)*x[1]*(1-x[0]*x[0] - x[1]*x[1])','sin(pi*t)*x[0]*(1-x[0]*x[0]-x[1]*x[1])'),t=t,degree = 2)
def get_p_exact(t):
	return Expression('sin(pi*t)*(x[0] + x[1])',t=t,degree =1 )
def get_bcs(t):
	bc = [bc0,bc1]
	return bc


print("DOF " + str(W.dim()))
