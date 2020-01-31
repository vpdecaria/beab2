import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

#A driver to access simulation results
def getData(filename):

	f = open(filename,'r')

	#First 3 lines are header information
	f.readline()
	f.readline()
	f.readline()

	t = np.array([])
	k = np.array([])
	filters = np.array([])
	normU = np.array([])
	normP = np.array([])
	normUExact = np.array([])
	normPExact = np.array([])
	Uerror = np.array([])
	Perror=np.array([])

	for j,line in enumerate(f):
		if line == 'END\n':
			break
		data = line.strip().split(',')
		data = [float(i) for i in data]
		t = np.append(t,data[0])
		k = np.append(k, data[1])
		filters = np.append(filters,data[2])
		normU = np.append(normU,data[3])
		normP = np.append(normP,data[4])
		normUExact = np.append(normUExact,data[5])
		normPExact = np.append(normPExact,data[6])
		Uerror = np.append(Uerror,data[7])
		Perror=np.append(Perror,data[8])


	f.close()
	return [t,k,filters,normU,normP,normUExact,normPExact,Uerror,Perror]

#Plotting Method

#t is the independent variable.
#methodData is a list of data from however many methods. Each data should be the same size as dt.
#title - string, title of plot
#xLabel - string, x-axis label
#lineType - list of strings to indicate type of lines we want drawn.
#labels - list of strings for titles in legend
def plotCompareMethods(dt, methodData,xLabel,title,lineType,labels):
	numberOfMethods = len(methodData)  #How many methods are we comparing?
	if(lineType == []):
		for j, data in enumerate(methodData):
			plt.loglog(dt,data,label = labels[j])
	else:
		for j, data in enumerate(methodData):
			plt.loglog(dt,data,lineType[j],label = labels[j])
	plt.xlabel(xLabel,fontsize = 18)
	plt.title(title,fontsize = 18)
	plt.xlim([25e-4,dt[0]+0.02])
	plt.legend(loc = 4)
	
	plt.show()

#gamma  = 1e0
[t,k,filters,normU,normP,normUExact,normPExact,Uerror,Perror] = getData('../errors/convergenceTestWRTtoleranceHarder/order-12-tol-1e-2.txt')	
[t1,k1,filters1,normU1,normP1,normUExact1,normPExact1,Uerror1,Perror1] = getData('../errors/convergenceTestWRTtoleranceHarder/order-12-constantstep-tol-1e-2.txt')	

plt.figure(1)
tFine = np.linspace(0,45,10000)
#plt.show()
plt.plot(t,normU,marker = 'x',label=r'MOOSE-IMEX-12, $TOL=10^{-2}$')
plt.plot(t1,normU1,'--',label='BE-AB2+F - nonadaptive')
#plt.plot(t,normUExact)


def cutoff(t):
	if t <=0:
		return 0
	else:
		return np.exp(-1./((10*t)**10))

	
def flatQuickStepUp(t):
	return cutoff(t-9)
def flatQuickStepDown(t):
	return 1- cutoff(t-9)

def smoothSteps(t):
	if t % 20 <= 10 - 1e-10:
		return flatQuickStepUp(t%10)
	else:
		return flatQuickStepDown(t%10)

	
exact = [np.sqrt(2)*np.pi*smoothSteps(T) for T in tFine]
plt.plot(tFine,exact,'k',label=r'Exact $||u||$')
plt.legend(fontsize = 14)
plt.ylabel(r'$||u||$',fontsize = 18)
plt.xlabel('t',fontsize = 18)
plt.xlim([0,45])
plt.title('Adaptive vs nonadaptive velocity norms\n with same number of Stokes solves',fontsize=16)
plt.show()

plt.figure(1)
tFine = np.linspace(0,45,10000)
#plt.show()
plt.plot(t,normP,marker = 'x',label=r'MOOSE-IMEX-12, $TOL=10^{-2}$')
plt.plot(t1,normP1,'--',label='BE-AB2+F - nonadaptive')
#plt.plot(t,normUExact)


def cutoff(t):
	if t <=0:
		return 0
	else:
		return np.exp(-1./((10*t)**10))

	
def flatQuickStepUp(t):
	return cutoff(t-9)
def flatQuickStepDown(t):
	return 1- cutoff(t-9)

def smoothSteps(t):
	if t % 20 <= 10 - 1e-10:
		return flatQuickStepUp(t%10)
	else:
		return flatQuickStepDown(t%10)

	
exact = [np.sqrt(2)*np.pi*smoothSteps(T) for T in tFine]
plt.plot(t,normPExact,'k',label=r'Exact $||u||$')
plt.legend(fontsize = 14)
plt.ylabel(r'$||p||$',fontsize = 18)
plt.xlabel('t',fontsize = 18)
plt.xlim([0,45])
plt.title('Adaptive vs nonadaptive pressure norms\n with same number of Stokes solves',fontsize=16)
plt.show()
