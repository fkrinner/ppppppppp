import numpy as np
import numpy.linalg as la
import os

def toComplex(string):
	vals = string[1:-1].split(',')
	return float(vals[0]) + 1.j*float(vals[1])

def pyComplex(string):
	return complex(string)

class integral:
	def __init__(self, inFileName, conversionMode = toComplex):
		if not self.loadFile(inFileName, conversionMode = conversionMode):
			raise IOError("could not load integral file")
		self.normalized = False
		self.hasEigen   = False
		self.numLim     = 1.e-10

	def loadFile(self, fileName, conversionMode = toComplex):
		if not os.path.isfile(fileName):
			print "File '" + fileName + "' does not exist"
			return False
		matrix = []
		with open(fileName, 'r') as inin:
			for line in inin.readlines():
				vals = [conversionMode(chunk) for chunk in line.split()]
				matrix.append(vals)
		for line in matrix:
			if not len(line) == len(matrix):
				print "Matrix is not quadratic"
				return False
		self.nAmpl = len(matrix)
		self.integralMatrix = np.asarray(matrix, dtype = complex)
		return True

	def removeIndices(self, indices):
		if isinstance(indices, int):
			indices = [indices]
		if not isinstance(indices, list):
			raise ValueError("Indices is not a list")
		self.hasEigen = False
		newIndexMap = []
		for i in range(self.nAmpl):
			if not i in indices:
				newIndexMap.append(i)
		self.nAmpl = len(newIndexMap)
		newMatrix = np.zeros((self.nAmpl, self.nAmpl), dtype = complex)
		for i,I in enumerate(newIndexMap):
			for j,J in enumerate(newIndexMap):
				newMatrix[i,j] = self.integralMatrix[I,J]
		self.integralMatrix = newMatrix

	def eigen(self):
		self.wasNormalized = self.normalized
		val, self.vec = la.eig(self.integralMatrix)
		vals = []
		for v in val:
			if abs(v.imag) > self.numLim:
				raise ValueError("Eigenvalue not real " + str(v))
			vals.append(v.real)
		self.val = np.asarray(vals)
		self.hasEigen = True

	def getSmallVectors(self, maxVal):
		if not self.hasEigen:
			raise RuntimeError("Eigensystem not calculated")
		vectors = []
		vals    = []
		for i in range(self.nAmpl):
			if self.val[i] < maxVal and not self.val[i] == 0.:
				vec = np.zeros((self.nAmpl), dtype = complex)
				print self.val[i]
				vals.append(self.val[i])
				for j in range(self.nAmpl):
					vec[j] = self.vec[j,i]
				vectors.append(vec)
		return vals, vectors

	def getBiggestVector(self):
		if not self.hasEigen:
			raise RuntimeError("Eigensystem not calculated")
		vectors = [[]]
		maxVal = 0.
		maxI   = 0
		for i in range(self.nAmpl):
			if self.val[i] > maxVal:
				maxVal = self.val[i]
		print maxVal
		vals = [maxVal]
		for i in range(self.nAmpl):
			vectors[0].append(self.vec[i,maxI])
		return vals, vectors		


	def norm(self):
		if self.normalized:
			return
		self.norms = np.zeros((self.nAmpl))
		for i in range(self.nAmpl):
			self.norms[i] = self.integralMatrix[i,i].real		
		for i in range(self.nAmpl):
			for j in range(self.nAmpl):
				if i == j:
					continue
				norm = (self.integralMatrix[i,i]*self.integralMatrix[j,j])**.5
				if norm == 0.+0.j:
					continue
				self.integralMatrix[i,j] /= norm
		for i in range(self.nAmpl):
			if not self.integralMatrix[i,i] == 0.:
				self.integralMatrix[i,i] = 1.+0.j
		self.normalized = True

	def unnorm(self):
		if not self.normalized:
			return
		for i in range(self.nAmpl):
			for j in range(self.nAmpl):
				self.integralMatrix[i,j] *= (self.norms[i] * self.norms[j])**.5
		self.normalized = False

def main():
	fileName = "./build/integral.dat"
#	fileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/integralBinStudy/resultat/intMatrix_15696_104772"
	inte     = integral(fileName, toComplex)
#	inte.removeIndices(0)
#	inte.norm()
	inte.eigen()
	vals,vecs = inte.getSmallVectors(0.005)
	for i in range(len(vecs[0])):
		print vecs[0][i]
	print vals


if __name__ == "__main__":
	main()
