#!/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Python_ultra/Python-2.7.10/bin/python
# backgroundModelFromDalitz.py
# Created: 2018-08-08 14:49:29.042090
# Author: Fabian Krinner
import os, sys
from rootfabi import root_open

from math import sin, cos, pi

import numpy as np
import numpy.linalg as la
import scipy.special as sps

from analyzeBELLE import loadConstants

mD0 = 1.86483;
#mD0 = 1.905;

mKs = 0.49761;
mPi = 0.13957;
mel = 0.0005109989461;

def get_ABC(vals, bins, basisFunctions):
	nF = len(basisFunctions)
	A  = np.zeros((nF,nF))
	B  = np.zeros((nF))
	C  = 0.
	for i, centers in enumerate(bins):
		v = vals[i]
		e = vals[i]**.5
		C += (v/e)**2
		fVals = [f(*centers) for f in basisFunctions]
		for i in range(nF):
			B[i] -= 2*fVals[i]*v/e**2
			for j in range(nF):
				A[i,j] += fVals[i]*fVals[j]/e**2
	return A,B,C

class legendre_with_transformation:
	def __init__(self, degX, xMin, xMax, degY, yMin, yMax):
		self.degX   = degX
		self.xMin   = xMin
		self.xMax   = xMax
		self.xRange = xMax - xMin

		self.degY   = degY
		self.yMin   = yMin
		self.yMax   = yMax
		self.yRange = yMax - yMin

		self.Xp = sps.legendre(self.degX)
		self.Yp =  sps.legendre(self.degY)
	def __call__(self, mpp2, mKp2):

		x,y,j = getCosTmPP(mpp2, mKp2)
		X = 2*(x-self.xMin)/self.xRange - 1.
		Y = 2*(y-self.yMin)/self.yRange - 1.
		return self.Xp(X)*self.Yp(Y)

	def __str__(self):
		return "legendre_with_transformation(" + str(self.degX) + "," + str(self.xMin) + "," + str(self.xMax) + ")(" + str(self.degY) + "," + str(self.yMin) + "," + str(self.yMax) + ")"

	def getCoefficients(self):
		retVal = [[0.]*(self.degY+1) for _ in range(self.degX+1)]
		for i in range(self.degX+1):
			for j in range(self.degY+1):
				retVal[i][j] = self.Xp.c[~i]*self.Yp.c[~j]
		return retVal

class legendre:
	def __init__(self, degX, xMin, xMax, degY, yMin, yMax):
		self.degX   = degX
		self.xMin   = xMin
		self.xMax   = xMax
		self.xRange = xMax - xMin

		self.degY   = degY
		self.yMin   = yMin
		self.yMax   = yMax
		self.yRange = yMax - yMin

		self.Xp = sps.legendre(self.degX)
		self.Yp = sps.legendre(self.degY)

	def __call__(self, x,y):
		X = 2*(x-self.xMin)/self.xRange - 1.
		Y = 2*(y-self.yMin)/self.yRange - 1.
		return self.Xp(X)*self.Yp(Y)

	def __str__(self):
		return "legendre(" + str(self.degX) + "," + str(self.xMin) + "," + str(self.xMax) + ")(" + str(self.degY) + "," + str(self.yMin) + "," + str(self.yMax) + ")"

	def getCoefficients(self):
		retVal = [[0.]*(self.degY+1) for _ in range(self.degX+1)]
		for i in range(self.degX+1):
			for j in range(self.degY+1):
				retVal[i][j] = self.Xp.c[~i]*self.Yp.c[~j]
		return retVal

class chebyt:
	def __init__(self, degX, xMin, xMax, degY, yMin, yMax):
		self.degX   = degX
		self.xMin   = xMin
		self.xMax   = xMax
		self.xRange = xMax - xMin

		self.degY   = degY
		self.yMin   = yMin
		self.yMax   = yMax
		self.yRange = yMax - yMin

		self.Xp = sps.chebyt(self.degX)
		self.Yp = sps.chebyt(self.degY)

	def __call__(self, x,y):
		X = 2*(x-self.xMin)/self.xRange - 1.
		Y = 2*(y-self.yMin)/self.yRange - 1.
		return self.Xp(X)*self.Yp(Y)

	def __str__(self):
		return "legendre(" + str(self.degX) + "," + str(self.xMin) + "," + str(self.xMax) + ")(" + str(self.degY) + "," + str(self.yMin) + "," + str(self.yMax) + ")"

	def getCoefficients(self):
		retVal = [[0.]*(self.degY+1) for _ in range(self.degX+1)]
		for i in range(self.degX+1):
			for j in range(self.degY+1):
				retVal[i][j] = self.Xp.c[~i]*self.Yp.c[~j]
		return retVal


def cut(st, A,B,C):
	a = np.zeros((len(st),len(st)))
	b = np.zeros((len(st)))
	for i,I in enumerate(st):
		b[i] = B[I]
		for j,J in enumerate(st):
			a[i,j] = A[I,J]
	return a,b,C

def getCosTmPP(mpp2, mKp2):
	mpp = mpp2**.5

	p1pK = (mKp2 - mPi*mPi - mKs*mKs)/2
	Epi  = mpp2**.5/2.
	EKs  = (mD0*mD0 - mpp2 - mKs*mKs)/4/Epi
	pPi  = (Epi*Epi - mPi*mPi)**.5
	pKs  = (EKs*EKs - mKs*mKs)**.5
	cosT = (Epi*EKs - p1pK)/pPi/pKs

	DmppDmpp2  = 1./(2*mpp2**.5)
	DcosTDmKp2 = -2*cosT/(mD0**2 - 2*mKp2 + mKs**2 + 2 * mPi**2 - mpp2)
#	print "derivatives:",DmppDmpp2, DcosTDmKp2
	jac = DmppDmpp2 * DcosTDmKp2

	return cosT, mpp, jac

def c2ndfScan(masses, values, fSet):
	bestChi2NDF = float("inf")
	A,B,C       = get_ABC(values, masses, fSet)
	bestSet     = []
	while True:
		if len(bestSet) == len(fSet):
			print "best Chi2/DNF with ALL functions"
			break
		bestChi2inThisSize = float("inf")
		bestI = None
		for i in range(len(fSet)):
			if i in bestSet:
				continue
			newSet  = bestSet + [i]
			a,b,c   = cut(newSet,A,B,C)
			cff     =-np.dot(la.inv(a+a.T), b)
			chi2ndf = np.dot(cff, np.dot(a,cff)) + np.dot(cff,b) + C
			if chi2ndf < bestChi2inThisSize:
				bestChi2inThisSize = chi2ndf
				bestI = i
		newChi2NDF = bestChi2inThisSize/(len(masses) - len(newSet))
		if newChi2NDF < bestChi2NDF:
			bestSet.append(bestI)
			bestChi2NDF = newChi2NDF
			print "Appending",bestI,"to get",bestChi2NDF
		else:
			print "The best set is",bestSet,'||',len(bestSet)
			break
	return bestSet

def main():
	inFileName = "./build/binned_dalitzs.root"
	histName2  = "dalitz_1.815000-1.915000"
#	histName   = "cosTp_vs_m2pp"
	masses     = [ ]
	values     = [ ]
	with root_open(inFileName, "READ") as inFile:
		dalitz = inFile.Get(histName2)

		for i in range(dalitz.GetNbinsX()):
			x = dalitz.GetXaxis().GetBinCenter(i+1)
			for j in range(dalitz.GetNbinsY()):
				val = dalitz.GetBinContent(i+1, j+1)
				if val > 0.:
					y = dalitz.GetYaxis().GetBinCenter(j+1)
					try:
						c,m,j = getCosTmPP(x,y)
					except ValueError:
						continue
					if abs(c) > 1.:
						continue
					masses.append((x,y))
					values.append(val)

#		xMin = -1.0 
#		xMax =  1.0 
#		yMin =  0.26
#		yMax =  1.44913767462

		xMin = min([m[0] for m in masses])
		xMax = max([m[0] for m in masses])
		yMin = min([m[1] for m in masses])
		yMax = max([m[1] for m in masses])

#		for m in masses:
#			yMin = min(yMin, m[1])
#			yMax = max(yMax, m[1])

		print xMin, "< x <",xMax
		print yMin, "< y <",yMax

		fSet   = []
		maxDim = 20
		for i in range(maxDim):
			if i == 1:
				continue
			for j in range(maxDim):
				if j == 1:
					continue
#				fSet.append(legendre_with_transformation(i,xMin,xMax,j,yMin,yMax))
#				fSet.append(legendre(i, xMin, xMax, j, yMin, yMax))

		bestSet = c2ndfScan(masses, values, fSet)
		return

		bestModel = [0, 1, 30, 2, 20, 5, 44, 23, 71, 21, 40, 3, 4, 14, 6,
		             24, 17, 80, 54, 9, 79, 50, 47, 29, 97, 77, 78, 90,
		             13, 31, 96, 98, 32, 27, 45, 62, 18, 95, 48, 26, 91,
		             38, 82, 41, 65, 74, 28, 58, 70, 51, 64, 16, 81, 61,
		             84, 8, 92, 36, 7, 25, 34, 60, 76, 66, 52, 39, 15, 
		             75, 33, 19, 67, 83, 72, 46, 53, 99, 94, 63, 22, 73,
		             59, 85, 86, 37, 88, 49, 87, 68, 56, 93, 55, 43, 42] ## len([...]) = 93

		print len(bestModel)

		bestFset = []
		for i in bestModel:
			bestFset.append(fSet[i])

		maxDegX = 0
		maxDegY = 0
		for f in bestFset:
			maxDegX = max(maxDegX, f.degX) 
			maxDegY = max(maxDegY, f.degY)

		overAllCoeffs = [[0.]*(maxDegY+1) for _ in range(maxDegX+1)]

		A,B,C  = get_ABC(values, masses, bestFset)
		coeffs = -np.dot(la.inv(A+A.T), B)
		for n,f in enumerate(bestFset):
			legendre_coeffs = f.getCoefficients()
			for i in range(len(legendre_coeffs)):
				for j in range(len(legendre_coeffs[i])):
					overAllCoeffs[i][j] += coeffs[n] * legendre_coeffs[i][j]
		with open("./build/dalitzBackground_version2.pol", 'w') as outFile:
			for i in range(maxDegX+1):
				for j in range(maxDegY+1):
					coeff = overAllCoeffs[i][j]
					if not coeff == 0.:
						outFile.write(str(i) + ' ' + str(j) + ' ' + str(coeff) + '\n')

#		chi2   = np.dot(coeffs, np.dot(A,coeffs)) + np.dot(coeffs,B) + C
#		print "chi2/NDF =",chi2/(len(masses) - len(bestFset))
#
		return 

		histTheo = hist.Clone()
		histTheo.Reset()
		histTheo.SetTitle("Hist theo cTm")
		histTheo.SetName("Hist_theo_cTm")

		for c,centers in enumerate(masses):
			i   = hist.GetXaxis().FindBin(centers[0])
			j   = hist.GetYaxis().FindBin(centers[1])
			dat = hist.GetBinContent(i,j)
			if dat == 0.:
				continue
			val = 0.
			for f,func in enumerate(bestFset):
				val += coeffs[f]*func(*centers)
			histTheo.SetBinContent(i,j,val)

		c2 = 0.
		ndf = 0
		jacPlot  = dalitz.Clone()
		dalModel = dalitz.Clone()

		jacPlot.Reset()
		jacPlot.SetTitle("Jacobian")
		jacPlot.SetName("Jacobian")

		dalModel.Reset()
		dalModel.SetTitle("Dalitz model")
		dalModel.SetName("Dalitz_model")

		for i in range(dalitz.GetNbinsX()):
			mpp2 = dalitz.GetXaxis().GetBinCenter(i+1)
			for j in range(dalitz.GetNbinsY()):
				dat = dalitz.GetBinContent(i+1,j+1)
				if dat == 0.:
					continue
				ndf += 1
				mKp2 = dalitz.GetYaxis().GetBinCenter(j+1)
				try:
					cosT, mpp, jac = getCosTmPP(mpp2, mKp2)
					print mpp2, mKp2, cosT, mpp
				except ValueError:
					print "ERRR at",mpp2, mKp2
					continue
				if abs(cosT) > 1:
					continue

				val = 0.
				try:
					for f,func in enumerate(bestFset):
		#				val += coeffs[f]*func(cosT,mpp)
						val += coeffs[f]*func(mpp2, mKp2)
					dalModel.SetBinContent(i+1,j+1, val)
				except ValueError:
					continue

				with root_open("polotossos.root","RECREATE"):
					hist.Write()
					dalitz.Write()
					histTheo.Write()
#					jacPlot.Write()
					dalModel.Write()

		dalitz.Draw("colz")

		print "Im dalitz", c2/(ndf - len(bestModel))

#		histTheo.Draw("colz")
#		raw_input()

if __name__ == "__main__":
	main()

