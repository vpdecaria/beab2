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
[t,k,filters,normU,normP,normUExact,normPExact,Uerror,Perror] = getData('1e-3second.txt')	
[t1,k1,filters1,normU1,normP1,normUExact1,normPExact1,Uerror1,Perror1] = getData('1e-4second-constant-step.txt')	

plt.figure(1)

plt.subplot(211)
plt.title('VSVO versus Constant Stepsize Constant Order')
plt.plot(t,filters+1)
plt.ylabel('Order')
plt.yticks([1,2])
plt.xlim([0,45])
#plt.plot(t1,filters1)
#
#plt.plot(t,TotalModelDissipation1e6)
#plt.plot(t,KineticEnergy1e6+ TotalModelDissipation1e6)
#plt.plot(t,TotalModelDissipation1e10)
#plt.plot(t,TotalModelDissipation1e2)
plt.subplot(212)
tFine = np.linspace(0,45,10000)
#plt.show()
plt.plot(t,normU,marker = 'x',label=r'VSVO-12, $TOL=10^{-3}$, 342 steps')
plt.plot(t1,normU1,'--',label='second order - nonadaptive, 535 steps')
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
plt.legend()
plt.ylabel(r'$||u||$')
plt.xlabel('t')
plt.xlim([0,45])
plt.show()

"""
time = 78
gammas = [1.0,1e2,1e4,1e6,1e8,1e10]
ModelDissipation1sec = [ModelDissipation1e0[time], ModelDissipation1e2[time],\
	ModelDissipation1e4[time],ModelDissipation1e6[time],\
	ModelDissipation1e8[time],ModelDissipation1e10[time]]
	
TotalModelDissipation1sec = [TotalModelDissipation1e0[time], TotalModelDissipation1e2[time],\
	TotalModelDissipation1e4[time],TotalModelDissipation1e6[time],\
	TotalModelDissipation1e8[time],TotalModelDissipation1e10[time]]

plt.semilogx(gammas,TotalModelDissipation1sec)
plt.title("Cumulative Dissipation at time " + str(t[time]), fontsize = 18)
plt.xlabel(r"\gamma")
plt.show()
"""
