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

def regularize(inFileNames, outFileName):
	matrix = parseMatrixFile(inFileNames[0])
	for i in range(1,len(inFileNames)):
		matrix += parseMatrixFile(inFileNames[i]) 
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
	import ROOT

	inFileNames  = ["./build/ac_integral_finer_binning.dat","./build/ac_model_integral.dat"]
	inFileNames  = ["./build/ps_integral_finer_binning.dat","./build/ps_model_integral.dat"]
	mm = parseMatrixFile(inFileNames[0])
	mm += parseMatrixFile(inFileNames[1])
	delta = 103
	hist = ROOT.TH1D('h','h',100, .9, 1.1)
	count = 0
	for i in range(delta):
		for j in range(delta):
			if mm[i,j] != 0.:
				rat = mm[i+delta,j+delta]/mm[i,j]
				hist.Fill(rat)
				count += 1
			if mm[i+delta,j] != 0.:
				rat = mm[i,j+delta]/mm[i+delta,j]
				hist.Fill(rat)
				count += 1
			if mm[i,j+delta] != 0.:
				rat = mm[i+delta,j]/mm[i,j+delta]
				hist.Fill(rat)
				count += 1
	print count
	hist.Draw()
	raw_input()


	outFileName = "./build/ac_integral_bg_model_regular_merged.dat"
#	regularize(inFileNames, outFileName)
	inFileNames  = ["./build/ps_integral_finer_binning.dat","./build/ps_model_integral.dat"]
	outFileName = "./build/ps_integral_bg_model_regular_merged.dat"
#	regularize(inFileNames, outFileName)

if __name__ == "__main__":
	main()
