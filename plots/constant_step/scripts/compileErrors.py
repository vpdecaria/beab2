import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

#This function handles the case that casting to
#float is impossible.
def to_float(input_string):
	try:
		new_float=float(input_string)
		return new_float
	except:
		return float('nan')

#A driver to access simulation results
def getWorkData(filename):

	f = open(filename,'r')

	#First 3 lines are header information
	f.readline()
	f.readline()
	f.readline()

	counter = 0

	while True:
		line = f.readline()
		if line == 'END\n':
			break
		counter += 1
	f.readline()
	line = f.readline()
	data = line.strip().split(',')
	data = [to_float(i) for i in data]
	l2L2 = data[1]
	l2H1 = data[2]
	l2L2Pressure = data[3]
	dt = data[5]
	
	
	f.close()
	return [counter,l2L2,l2H1,l2L2Pressure,dt]

#Plotting Method

def compileErrors(fName,dtList):
	counter = np.array([])
	tol = np.array([])
	l2L2 = np.array([])
	l2H1 = np.array([])
	l2L2Pressure = np.array([])
	dts = np.array([])
	rejections = np.array([])
	for dt in dtList:
		[countertemp,l2L2temp,l2H1temp,l2L2Pressuretemp,currentDt] = getWorkData(fName + str(dt)+'.txt')
		counter = np.append(counter,countertemp)
		l2L2 = np.append(l2L2,l2L2temp)
		l2H1 = np.append(l2H1,l2H1temp)
		l2L2Pressure = np.append(l2L2Pressure,l2L2Pressuretemp)
		dts = np.append(dts,currentDt)
	return [counter,l2L2,l2H1,l2L2Pressure,dts]

def writeErrorsToFile(fName,dtList,outputfname):
	[counter,l2L2,l2H1,l2L2Pressure,dts] = compileErrors(fName,dtList)
	f = open(outputfname,'w')
	f.write('Junk\n')
	f.write('dt, l2L2 error, l2H1 error, l2L2 Pressure error \n')
	for i, dt in enumerate(dts):
		output = "\n" + str(dt) + "," + str(l2L2[i]) + "," +str(l2H1[i])  + "," + str(l2L2Pressure[i])
		f.write(output)
		
	f.close()
		
#t is the independent variable.
#methodData is a list of data from however many methods. Each data should be the same size as dt.
#title - string, title of plot
#xLabel - string, x-axis label
#lineType - list of strings to indicate type of lines we want drawn.
#labels - list of strings for titles in legend


dtList = [0.2,0.1,0.05, 0.025,0.0125,0.00625]

writeErrorsToFile('../../../errors/convergenceTestConstantStep/order-1-extrap-fe-dt-',dtList,'order-1-extrap-fe.txt')
writeErrorsToFile('../../../errors/convergenceTestConstantStep/order-1-extrap-ab2-dt-',dtList,'order-1-extrap-ab2.txt')
writeErrorsToFile('../../../errors/convergenceTestConstantStep/order-2-extrap-fe-dt-',dtList,'order-2-extrap-fe.txt')
writeErrorsToFile('../../../errors/convergenceTestConstantStep/order-2-extrap-ab2-dt-',dtList,'order-2-extrap-ab2.txt')