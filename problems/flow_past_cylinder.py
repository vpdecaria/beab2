from __future__ import print_function
from fenics import *
import numpy as np
from mshr import *

#Basic Problem info
nu = 1.0/1000.0             # kinematic viscosity
Re = 1.0 / nu
T = 8.1           # final time

N=2

#Mesh
#mesh = Mesh('cylinder.coarser.mesh.xml.gz')

segments_small = 100
resolution     = 30   #Worked well
resolution     = 90

rectangleHeight = .41
rectangleLength = 2.2

BigRectangle = Rectangle(Point(0,0), Point(rectangleLength, rectangleHeight))
smallCircle = Circle(Point(0.2,0.2), 0.05, segments_small)

domain = BigRectangle - smallCircle

mesh = generate_mesh(domain,resolution)

# Define function spaces

V = VectorElement("Lagrange", mesh.ufl_cell(), 2)
Q = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = V * Q

W = FunctionSpace(mesh, TH)
dw = TrialFunction(W)

#Define BC
class noslip_boundary(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and \
                    (x[1] > rectangleHeight - DOLFIN_EPS or x[1] < 0 + DOLFIN_EPS  or \
                    (x[0]-.2)**2 + (x[1]-.2)**2 < 0.0025 + 3*DOLFIN_EPS)
    
class outflow_boundary(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and x[0] > rectangleLength - 40*DOLFIN_EPS
        
class inflow_boundary(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and x[0] < DOLFIN_EPS

class pressure_point_boundary(SubDomain):
	def inside(self,x,on_boundary):
		return (on_boundary and near(x[0],0) and near(x[1],0))

mesh_points = mesh.coordinates()
class OriginPoint(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[0], mesh_points[0,0])  and  near(x[1], mesh_points[0,1]))
        
Inflow = Expression(('pow(0.41,-2)*6*x[1]*(0.41-x[1])','0'),degree = 2)

bc0        = DirichletBC(W.sub(0), (0,0), noslip_boundary())
bc1        = DirichletBC(W.sub(0), Inflow, inflow_boundary())
bc2        = DirichletBC(W.sub(0), Inflow, outflow_boundary())

bcpressure = DirichletBC(W.sub(1), (0), OriginPoint(),"pointwise")

bc = [bc0,bc1,bc2,bcpressure]

#Create Functions that returns the rhs at every time
def get_f(t):
    return Expression(('0','0'),degree = 2)
#Create functions that return dummy exact solutions at every time
def get_u_exact(t):
    return Expression(('0','0'),degree = 2)
def get_p_exact(t):
    return Expression( '0',     degree = 1 )
def get_bcs(t):
    bc0    = DirichletBC(W.sub(0), (0,0), noslip_boundary())
    Inflow = Expression(('pow(0.41,-2)*sin(pi*t/8.0)*6*x[1]*(0.41-x[1])','0'),t=t,degree = 2)
    bc1    = DirichletBC(W.sub(0), Inflow, inflow_boundary())
    bc2    = DirichletBC(W.sub(0), Inflow, outflow_boundary())
    bcpressure = DirichletBC(W.sub(1), (0), OriginPoint(),"pointwise")
    bc     = [bc0,bc1,bc2,bcpressure]
    return bc

print("DOF " + str(W.dim()))
