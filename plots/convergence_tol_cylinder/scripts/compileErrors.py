import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

#A driver to access simulation results

#This function handles the case that casting to
#float is impossible.
def to_float(input_string):
	try:
		new_float=float(input_string)
		return new_float
	except:
		return float('nan')

def getWorkData(filename):

	f = open(filename,'r')

	#First 3 lines are header information
	f.readline()
	f.readline()
	f.readline()

	counter = 0

	dt_average = 0

	while True:
		line = f.readline()

		if line == 'END\n':
			break
		data = line.strip().split(',')
		data = [to_float(i) for i in data]
		dt_average += data[1]

		counter += 1
	f.readline()
	line = f.readline()
	data = line.strip().split(',')
	data = [to_float(i) for i in data]
	l2L2 = data[1]
	l2H1 = data[2]
	l2L2Pressure = data[3]
	tol = data[0]
	time = data[6]
	dt_average/=counter
	
	
	f.close()
	return [counter,l2L2,l2H1,l2L2Pressure,tol,time,dt_average]

#Plotting Method

def compileErrors(fName,tolList):
	counter = np.array([])
	tol = np.array([])
	l2L2 = np.array([])
	l2H1 = np.array([])
	l2L2Pressure = np.array([])
	tols = np.array([])
	rejections = np.array([])
	times = np.array([])
	dt_aves = np.array([])
	for tol in tolList:
		[countertemp,l2L2temp,l2H1temp,l2L2Pressuretemp,currentTol,time,dt_ave] = \
			getWorkData(fName + str(tol)+'.txt')
		counter = np.append(counter,countertemp)
		l2L2 = np.append(l2L2,l2L2temp)
		l2H1 = np.append(l2H1,l2H1temp)
		l2L2Pressure = np.append(l2L2Pressure,l2L2Pressuretemp)
		tols = np.append(tols,currentTol)
		times = np.append(times,time)
		dt_aves = np.append(dt_aves,dt_ave)
	return [counter,l2L2,l2H1,l2L2Pressure,tols,times,dt_aves]

def writeErrorsToFile(fName,tolList,outputfname):
	[counter,l2L2,l2H1,l2L2Pressure,tols, times,dt_aves] = compileErrors(fName,tolList)
	f = open(outputfname,'w')
	f.write('Junk\n')
	f.write('tolerance, l2L2 error, l2H1 error, l2L2 Pressure error, run time, dt average \n')
	for i, tol in enumerate(tols):
		output = "\n" + str(tol) + "," + str(l2L2[i]) + "," +str(l2H1[i])  + "," \
			+ str(l2L2Pressure[i]) + ',' + str(times[i]) + ',' + str(dt_aves[i])
		f.write(output)
		
	f.close()
		
#t is the independent variable.
#methodData is a list of data from however many methods. Each data should be the same size as dt.
#title - string, title of plot
#xLabel - string, x-axis label
#lineType - list of strings to indicate type of lines we want drawn.
#labels - list of strings for titles in legend


TolList = ['1', '1e-1', '1e-2', '1e-3' ,'1e-4' ,'1e-5' ,'1e-6', '1e-7']

TolList = ['1', '1e-1', '1e-2', '1e-3' ,'1e-4' ,'1e-5' ,'1e-6']


writeErrorsToFile('../../../errors/convergenceTestWRTtolerance/order-1-tol-',TolList,'order-1.txt')
writeErrorsToFile('../../../errors/convergenceTestWRTtolerance/order-12-tol-',TolList,'order-12.txt')
writeErrorsToFile('../../../errors/convergenceTestWRTtolerance/order-2-tol-',TolList,'order-2.txt')