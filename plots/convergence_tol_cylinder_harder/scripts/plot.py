import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

#A driver to access simulation results

markers = ['v','o','s','^','+']

def getData(filename):

	f = open(filename,'r')

	#First 3 lines are header information
	f.readline()
	f.readline()
	f.readline()
	
	dt = np.zeros(7)
	
	l2L2 = np.zeros_like(dt)
	wall_time = np.zeros_like(dt)
	l2L2Pressure = np.zeros_like(dt)
	stokes_solves = np.zeros_like(dt)
	for j,line in enumerate(f):
		data = line.strip().split(',')
		data = [float(i) for i in data]
		dt[j] = data[0]
		l2L2[j] = data[1]
		wall_time[j] = data[4]
		l2L2Pressure[j] = data[3]
		stokes_solves[j] = data[6]


	f.close()
	return [dt,l2L2,wall_time,l2L2Pressure,stokes_solves]

#Plotting Method

#t is the independent variable.
#methodData is a list of data from however many methods. Each data should be the same size as dt.
#title - string, title of plot
#xLabel - string, x-axis label
#lineType - list of strings to indicate type of lines we want drawn.
#labels - list of strings for titles in legend
def plotCompareMethods(dt, methodData,xLabel,yLabel,title,lineType,labels,markers,variable):
	numberOfMethods = len(methodData)  #How many methods are we comparing?
	if(lineType == []):
		for j, data in enumerate(methodData):
			plt.loglog(dt[j],data,label = labels[j])
	else:
		for j, data in enumerate(methodData):
			plt.loglog(dt[j],data,lineType[j],label = labels[j],marker = markers[j])
	plt.xlabel(xLabel,fontsize = 18)
	plt.ylabel(yLabel,fontsize=18)
	plt.title(title,fontsize = 18)
	#plt.xlim([50e-4,dt[0]+0.02])
	#plt.ylim([1e-6,1.7e-1])

	#velocity
	if(variable == "v"):  
		x_l1=337
		x_u1=13000
		x_l2=200
		x_u2=2000
		multiplier2=3e1
		multiplier1=7e-1
	elif (variable == "p"):#pressure
		x_l1=337
		x_u1=13000
		x_l2=200
		x_u2=2000
		multiplier2=8.5e1
		multiplier1=2.4e0    
	plt.loglog([x_l1,x_u1],[multiplier1*x_l1**(-1),multiplier1*x_u1**(-1)],'k:',label = 'slope -1',linewidth=2)
	plt.loglog([x_l2,x_u2],[multiplier2*x_l2**(-2),multiplier2*x_u2**(-2)],'k--',label = 'slope -2')
	plt.legend(fontsize = 14)
	
	plt.tight_layout()
	
	
	plt.show()

[dt,l2L21fe,wall_time1fe,l2L21fePressure,stokes_solves_1] = getData('order-1.txt')	
[dt,l2L21ab2,wall_time1ab2,l2L21ab2Pressure,stokes_solves_12] = getData('order-12.txt')	
#[dt,l2L22fe,wall_time2fe,l2L22fePressure] = getData('order-2.txt')	
[dt,l2L22ab2,wall_time2ab2,l2L22ab2Pressure,stokes_solves_2] = getData('order-2.txt')	

#l2L2

xLabel = 'Stokes solves'
yLabel = r'Relative $\ell^2(0,T;L^2(\Omega))$ error'
#title = r'Velocity Error'
title='Adaptive velocity error'
#lineType = ['k','k--','k-.','k.-']
lineType = []
lineType = ['k','b','g','r.-']
lineType = ['k','b','r','k']
markers = ['v','o','s','^']
legend= ['VSS BE-AB2',\
		'MOOSE-IMEX-12']


plotCompareMethods([stokes_solves_1[2:],stokes_solves_12[2:]],[l2L21fe[2:],l2L21ab2[2:]],xLabel,yLabel,title,lineType,legend,markers,'v')

#l2L2 pressure

xLabel = 'Stokes solves'
yLabel = r'Relative $\ell^2(0,T;L^2(\Omega))$ error'
#title = r'Pressure Error'
title='Adaptive pressure error'
lineType = ['k','k--','k-.','k.-']
lineType = ['k','k--','k-.']
lineType = ['k','b','r','k','m-.','c:','g']
markers = ['v','o','s','^','P','*','|']
legend= ['VSS BE-AB2',\
		'MOOSE-IMEX-12']
plotCompareMethods([stokes_solves_1[2:],stokes_solves_12[2:]],[l2L21fePressure[2:],l2L21ab2Pressure[2:]],xLabel,yLabel,title,lineType,legend,markers,'p')
