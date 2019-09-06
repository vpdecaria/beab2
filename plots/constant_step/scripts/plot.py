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
	
	dt = np.zeros(6)
	
	l2L2 = np.zeros_like(dt)
	l2H1 = np.zeros_like(dt)
	l2L2Pressure = np.zeros_like(dt)

	for j,line in enumerate(f):
		data = line.strip().split(',')
		data = [float(i) for i in data]
		dt[j] = data[0]
		l2L2[j] = data[1]
		l2H1[j] = data[2]
		l2L2Pressure[j] = data[3]


	f.close()
	return [dt,l2L2,l2H1,l2L2Pressure]

#Plotting Method

#t is the independent variable.
#methodData is a list of data from however many methods. Each data should be the same size as dt.
#title - string, title of plot
#xLabel - string, x-axis label
#lineType - list of strings to indicate type of lines we want drawn.
#labels - list of strings for titles in legend
def plotCompareMethods(dt, methodData,xLabel,yLabel,title,lineType,labels,markers):
	numberOfMethods = len(methodData)  #How many methods are we comparing?
	if(lineType == []):
		for j, data in enumerate(methodData):
			plt.loglog(dt,data,label = labels[j])
	else:
		for j, data in enumerate(methodData):
			plt.loglog(dt,data,lineType[j],label = labels[j],marker = markers[j])
	plt.xlabel(xLabel,fontsize = 18)
	plt.ylabel(yLabel,fontsize=18)
	plt.title(title,fontsize = 18)
	#plt.xlim([50e-4,dt[0]+0.02])
	#plt.ylim([1e-6,1.7e-1])
	x_l=0.025
	x_u=0.2
	#plt.loglog([x_l,x_u],[.5*x_l**2,.5*x_u**2],'k--',label = 'slope 2')
	plt.loglog([x_l,x_u],[4*x_l**2,4*x_u**2],'k--',label = 'slope 2')
	#plt.loglog([x_l,x_u],[x_l**3,x_u**3],'k-.',label = 'slope 3')
	#plt.loglog([x_l,x_u],[3*x_l**4,3*x_u**4],'k:',label = 'slope 4')
	plt.legend(loc = 4)
	
	
	
	
	plt.show()

[dt,l2L21fe,l2H11fe,l2L21fePressure] = getData('order-1-extrap-fe.txt')	
[dt,l2L21ab2,l2H11ab2,l2L21ab2Pressure] = getData('order-1-extrap-ab2.txt')	
[dt,l2L22fe,l2H12fe,l2L22fePressure] = getData('order-2-extrap-fe.txt')	
[dt,l2L22ab2,l2H12ab2,l2L22ab2Pressure] = getData('order-2-extrap-ab2.txt')	

#l2L2

xLabel = r'$\Delta t$'
yLabel = r'$\|u-u_h\|_{l2L2}/\|u\|_{l2L2}$'
#title = r'Velocity Error'
title='Nonadaptive Velocity Error'
#lineType = ['k','k--','k-.','k.-']
lineType = []
lineType = ['k','b','g','r.-']
lineType = ['k','b','r','k']
markers = ['v','o','s','^']
legend= ['BEFE',\
		'BEAB2',\
		'BEAB2+F']
plotCompareMethods(dt,[l2L21fe,l2L21ab2,l2L22ab2],xLabel,yLabel,title,lineType,legend,markers)

#l2L2 pressure

xLabel = r'k'
yLabel = r'$\|p-p_h\|_{l2L2}/\|p\|_{l2L2}$'
#title = r'Pressure Error'
title='Nonadaptive Pressure Error'
lineType = ['k','k--','k-.','k.-']
lineType = ['k','k--','k-.']
lineType = ['k','b','r','k','m-.','c:','g']
markers = ['v','o','s','^','P','*','|']
"""
legend= ['BDF3-Stab, m = '+'{0:.3f}'.format(slopel2L20),\
		'BDF3, m = '+'{0:.3f}'.format(slopel2L203),\
		'BDF3-4, m = '+'{0:.3f}'.format(slopel2L204)]
"""
legend= ['BEFE',\
		'BEAB2',\
		'BEAB2+F']
plotCompareMethods(dt,[l2L21fePressure,l2L21ab2Pressure,l2L22ab2Pressure],xLabel,yLabel,title,lineType,legend,markers)