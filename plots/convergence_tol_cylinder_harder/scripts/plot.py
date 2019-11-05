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

	for j,line in enumerate(f):
		data = line.strip().split(',')
		data = [float(i) for i in data]
		dt[j] = data[0]
		l2L2[j] = data[1]
		wall_time[j] = data[4]
		l2L2Pressure[j] = data[3]


	f.close()
	return [dt,l2L2,wall_time,l2L2Pressure]

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
			plt.loglog(dt[j],data,label = labels[j])
	else:
		for j, data in enumerate(methodData):
			plt.loglog(dt[j],data,lineType[j],label = labels[j],marker = markers[j])
	plt.xlabel(xLabel,fontsize = 18)
	plt.ylabel(yLabel,fontsize=18)
	plt.title(title,fontsize = 18)
	#plt.xlim([50e-4,dt[0]+0.02])
	#plt.ylim([1e-6,1.7e-1])
	x_l1=250
	x_u1=6500
	x_l2=28
	x_u2=1200
	multiplier2=2e2
	multiplier1=3e0
	plt.loglog([x_l1,x_u1],[multiplier1*x_l1**(-1),multiplier1*x_u1**(-1)],'k:',label = 'slope -1',linewidth=2)
	plt.loglog([x_l2,x_u2],[multiplier2*x_l2**(-2),multiplier2*x_u2**(-2)],'k--',label = 'slope -2')
	#plt.loglog([x_l,x_u],[x_l**3,x_u**3],'k-.',label = 'slope 3')
	#plt.loglog([x_l,x_u],[3*x_l**4,3*x_u**4],'k:',label = 'slope 4')
	plt.legend()
	
	
	
	
	plt.show()

[dt,l2L21fe,wall_time1fe,l2L21fePressure] = getData('order-1.txt')	
[dt,l2L21ab2,wall_time1ab2,l2L21ab2Pressure] = getData('order-12.txt')	
#[dt,l2L22fe,wall_time2fe,l2L22fePressure] = getData('order-2.txt')	
[dt,l2L22ab2,wall_time2ab2,l2L22ab2Pressure] = getData('order-2.txt')	

#l2L2

xLabel = 'runtime (s)'
yLabel = r'$\|u-u_h\|_{l2L2}/\|u\|_{l2L2}$'
#title = r'Velocity Error'
title='Adaptive Velocity Error'
#lineType = ['k','k--','k-.','k.-']
lineType = []
lineType = ['k','b','g','r.-']
lineType = ['k','b','r','k']
markers = ['v','o','s','^']
legend= ['BEFE',\
		'VSVO 12',\
		'BEAB2+F']


plotCompareMethods([wall_time1fe,wall_time1ab2,wall_time2ab2],[l2L21fe,l2L21ab2,l2L22ab2],xLabel,yLabel,title,lineType,legend,markers)

#l2L2 pressure

xLabel = 'runtime (s)'
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
		'VSVO 12',\
		'BEAB2+F']
plotCompareMethods([wall_time1fe,wall_time1ab2,wall_time2ab2],[l2L21fePressure,l2L21ab2Pressure,l2L22ab2Pressure],xLabel,yLabel,title,lineType,legend,markers)