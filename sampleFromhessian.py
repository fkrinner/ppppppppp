#!/usr/bin/python
# sampleFromhessian.py
# Created: 2018-12-06 13:33:51.387882
# Author: Fabian Krinner
import os, sys
import numpy as np
import numpy.linalg as la


def loadFile(inFileName):
	retVal = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			retVal.append([float(v) for v in line.split()])
	return np.array(retVal)

def checkHessian(matrix, numLim = 1.e-5):
	if not matrix.shape[0] == matrix.shape[1]:
		raise ValueError("Matrix is not quadratic")
	dim = len(matrix)
	for i in range(dim):
		for j in range(dim):
			if abs(matrix[i,j] - matrix[j,i]) > numLim:
				raise ValueError("Matrix is not symmetric")

def main():
	inFileName  = sys.argv[1]
	outFileName = inFileName.replace('parameterHessian','sample')

	nSample = 10000
	if len(sys.argv) > 2:
		nSample = int(sys.argv[2])

	hessian    = loadFile(inFileName)
	checkHessian(hessian)
	coma = la.inv(hessian)
	dim  = len(coma)

	mean = np.zeros(len(hessian))

	with open(outFileName,'w') as outFile:
		for i in range(nSample):
			smpl = np.random.multivariate_normal(mean, coma)
			for v in smpl:
				outFile.write(str(v) + ' ')
			outFile.write('\n')

	print "Sample file '" + outFileName + "' written"
	print "Now run './Dost.exe "+ inFileName.replace("BELLE_fit_results/hessians","BELLE_fit_results").replace("parameterHessian","fit") + " -coma'"


if __name__ == "__main__":
	main()
