#!/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Python_ultra/Python-2.7.10/bin/python
# regularize_integral_matrix.py
# Created: 2018-07-20 15:35:27.287781
# Author: Fabian Krinner
import os, sys
import numpy as np
import numpy.linalg as la

def toComplex(s):
	chunks = s[1:~0].split(',')
	return float(chunks[0]) + 1.j*float(chunks[1])

def toString(val):
	return '('+str(val.real)+','+str(val.imag)+')'

def parseMatrixFile(inFileName):
	matrix = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			lineVals = []
			chunks = line.split()
			for c in chunks:
				lineVals.append(toComplex(c))
			matrix.append(lineVals)
	return np.asarray(matrix, dtype = complex)

def writeMatrixFile(outFileName, matrix):
	with open(outFileName, 'w') as outFile:
		for i in range(len(matrix)):
			for j in range(len(matrix)):
				outFile.write(toString(matrix[i,j])+' ')
			outFile.write('\n')

def regularize(inFileName, outFileName):
	matrix = parseMatrixFile(inFileName)
	val, vec = la.eig(matrix)
	indices = []
	for i,v in enumerate(val):
		if v.real < 0.:
			print v," !!!"
			indices.append(i)
	for i in indices:
		for j in range(len(val)):
			for k in range(len(val)):
				matrix[j,k] -= 2*val[i]*vec[j,i]*vec[k,i]
	writeMatrixFile(outFileName, matrix)

def main():
	inFileName  = "./build/ac_integral_bg_model.dat"
	outFileName = "./build/ac_integral_bg_model_regular.dat"
	regularize(inFileName, outFileName)
	inFileName  = "./build/ps_integral_bg_model.dat"
	outFileName = "./build/ps_integral_bg_model_regular.dat"
	regularize(inFileName, outFileName)

if __name__ == "__main__":
	main()
