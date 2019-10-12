from __future__ import print_function
from fenics import *
import numpy as np
from mshr import *
import matplotlib.pyplot as plt

#Basic Problem info
nu = 1.0             # kinematic viscosity
T = 20               # final time

#_________________________________DEFINE THE OSCILLATIONS IN TIME___________________________________
def smooth_bridge(t):
	return np.exp(-np.exp(-1./(1-t)**2)/t**2)

period         = 10.0
peak           = 1.0
impulse_length = 5.0
def amplitude(t):
	if t%period<=impulse_length-1.0:
		return 1.0
	elif t%period<=impulse_length:
		theta = impulse_length-t%period
		return smooth_bridge(theta)
	elif t%period<=period-1.0:
		return 0.0
	else:
		theta = period-t%period
		return smooth_bridge(1-theta)
#Larger period should be multiple of period
largerPeriod=2*period
def posAndNegOscillations(t):
	if t%largerPeriod<=largerPeriod/2:
		return -amplitude(t)+1
	else:
		return amplitude(t)-1
def Amplitude(t):
	return posAndNegOscillations(t+period-1)
def leftCircleAmplitude(t):
	return posAndNegOscillations(t+1./4.*period-1)
def rightCircleAmplitude(t):
	return posAndNegOscillations(t+3./4.*period-1)


#_____________________________ DEFINE OSCILLATIONS IN TIME # 2 _____________________________________
TOLERANCE = 1e-14

def w_osc(t):
	horizontalScale = np.pi
	if np.sin(t*horizontalScale) > TOLERANCE:
		return np.exp(1)*np.exp(-1/np.sin(t*horizontalScale))
	else:
		return 0
	
def g_osc(t):
	horizontalScale = 1.
	if np.cos(t*horizontalScale) > TOLERANCE:
		return np.exp(1)*np.exp(-1/np.cos(t*horizontalScale))
	else:
		return 0
def wPrime(t):
	horizontalScale = np.pi
	if np.sin(t*horizontalScale) > TOLERANCE:
		return horizontalScale*np.cos(horizontalScale*t)/np.sin(horizontalScale*t)**2*w_osc(t)
	else:
		return 0
	
def gPrime(t):
	horizontalScale = 1.
	if np.cos(t*horizontalScale) > TOLERANCE:
		return -horizontalScale*np.sin(horizontalScale*t)/np.cos(horizontalScale*t)**2*g_osc(t)
	else:
		return 0

def amp(t):
	return np.array([0.5*w_osc(t) + 0.5*(1- g_osc(t))])

def amp_prime(t):
	return 0.5*wPrime(t) - 0.5*gPrime(t)

x = np.linspace(0,20,10000)
y = [Amplitude(t) for t in x]
y = [amp(t) for t in x]


#plt.plot(x,y)
#plt.show()
#exit()

#_______________________________________ DEFINE MESH _______________________________________________

N=1000
bigCircle = Circle(Point(0,0),1,N)

#mesh = generate_mesh(bigCircle,5) #about 7k dof P2p1, 16k p3P2
#mesh = generate_mesh(bigCircle,10) #about 7k dof P2p1, 16k p3P2
mesh = generate_mesh(bigCircle,5) #about 7k dof P2p1, 16k p3P2
#mesh = generate_mesh(bigCircle,20) # about 11k
#plot(mesh)
#plt.show()

# Define function spaces

#The exact solution is a cubic velocity and linear pressure. This element choice confines         
#the error mainly to the temporal discretization, and the error at the boundary of the circular
#domain. We took many edges along the boundary to diminish the boundary error.
V = VectorElement("Lagrange", mesh.ufl_cell(), 3)
Q = FiniteElement("Lagrange", mesh.ufl_cell(), 2)
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
#def get_f(t):
#	return Expression(('pi*x[1]*(-1 + pow(x[0],2) + pow(x[1],2))*cos(pi*t) \
#		                    - sin(pi*t)*(-1 + 8*nu*x[1] \
#		                    + x[0]*pow(-1 + pow(x[0],2) + pow(x[1],2),2)*sin(pi*t))'\
#	,				  '-pi*x[0]*(-1 + pow(x[0],2) + pow(x[1],2))*cos(pi*t) \
#	                        - sin(pi*t)*(-1 - 8*nu*x[0] \
#	                        + x[1]*pow(-1 + pow(x[0],2) + pow(x[1],2),2)*sin(pi*t))')\
#	                  ,t=t,nu=nu, degree = 2)
def get_f(t):
	amp_temp = amp(t)
	amp_prime_temp = amp_prime(t)
	return Expression(('x[1]*(-1 + pow(x[0],2) + pow(x[1],2))*amp_prime \
		                    - amp*(-1 + 8*nu*x[1] \
		                    + x[0]*pow(-1 + pow(x[0],2) + pow(x[1],2),2)*amp)'\
		              ,\
     				  '-x[0]*(-1 + pow(x[0],2) + pow(x[1],2))*amp_prime \
	                        - amp*(-1 - 8*nu*x[0] \
	                        + x[1]*pow(-1 + pow(x[0],2) + pow(x[1],2),2)*amp)')\
	                  ,t=t,nu=nu, amp = amp_temp, amp_prime = amp_prime_temp,degree = 3)
	#return Expression(('pow(0.41,-2)*1*6*x[1]*(0.41-x[1])','0'),degree = 2)
#Create functions that return dummy exact solutions at every time
#def get_u_exact(t):
#	return Expression(('-sin(pi*t)*x[1]*(1-x[0]*x[0] - x[1]*x[1])','sin(pi*t)*x[0]*(1-x[0]*x[0]-x[1]*x[1])'),t=t,degree = 2)
#def get_p_exact(t):
#	return Expression('sin(pi*t)*(x[0] + x[1])',t=t,degree =1 )
def get_u_exact(t):
	amp_temp = amp(t)
	return Expression(('-amp*x[1]*(1-x[0]*x[0] - x[1]*x[1])',\
		                'amp*x[0]*(1-x[0]*x[0] - x[1]*x[1])'),t=t,
						 amp = amp_temp, amp_prime = amp_prime(t),degree = 3)
def get_p_exact(t):
	amp_temp = amp(t)
	return Expression('amp*(x[0] + x[1])',t=t,amp = amp_temp, amp_prime = amp_prime(t),degree = 1)

def get_bcs(t):
	bc = [bc0,bc1]
	return bc


print("DOF " + str(W.dim()))
