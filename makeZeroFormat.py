import sys
sys.path.append("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/zeroModeFitter_2D")
import numpy as np
import numpy.linalg as la
from integral import toComplex, integral
import pyRootPwa
from rootfabi import root_open
import utils

binning3Pi = [.5]
while binning3Pi[-1] < 2.5:
	binning3Pi.append(binning3Pi[-1] + .04)
binning3Pi = np.asarray(binning3Pi, dtype = np.float64)

def loadBinning(inFileName):
	binning = []
	with open(inFileName, 'r') as inin:
		for line in inin.readlines():
			binning += [float(chunk)**.5 for chunk in line.split()]
	return np.asarray(binning)


def loadRealMatrixFile(fileName):
	matrix = []
	with open(fileName, 'r') as inin:
		for line in inin.readlines():
			vals = [float(chunk) for chunk in line.split()]
			matrix.append(vals)
	return np.asarray(matrix)

def loadAmplitudeFile(fileName):
	ampls = []
	with open(fileName, 'r') as inin:
		for line in inin.readlines():
			ampls.append(toComplex(line.strip()))
	return np.asarray(ampls, dtype = complex)

def projectOutPhaseDirection(matrix, amplitudes):
	phaseDirection = []
	norm           = 0.
	dim            = len(matrix)
	for a in amplitudes:
		phaseDirection.append(-a.imag)
		phaseDirection.append( a.real)
		norm += abs(a)**2
	norm **= .5
	for i in range(dim):
		phaseDirection[i] /= norm
	projector = np.zeros((dim,dim))
	for i in range(dim):
		projector[i,i] += 1.
		for j in range(dim):
			projector[i,j] -= phaseDirection[i] * phaseDirection[j]
	return np.dot(projector, np.dot(matrix, projector))			


def main():
	bin          = 34
	amplFileName = "./build/amplitudes.dat"
	hessFileName = "./build/hessian.dat"
	inteFileName = "./build/integral.dat"
	binFfileName = "./build/binningF0.dat"
	binRfileName = "./build/binningRho.dat"

	binningF0  = loadBinning(binFfileName)
	binningRho = loadBinning(binRfileName)

	nBins = len(binningF0) + len(binningRho) - 2

	sectors = {"D0[pi,pi]0++PiS" : (0,len(binningF0)-1), "D0[pi,pi]1--PiP":(len(binningF0)-1,len(binningF0)+len(binningRho)-2)}

	hessian = loadRealMatrixFile(hessFileName)
	ampl    = loadAmplitudeFile(amplFileName)

	inte = integral(inteFileName)
	inte.removeIndices(0)
	inte.norm()
	inte.eigen()
	norm     = inte.norms
	sev = inte.getSmallVectors(0.005)
	zeroMode = sev[0]
	

	comaName = "COMA_0_"+str(bin)
#	COMA     = la.inv(hessian)
	COMA     = utils.pinv(hessian)
	COMA     = projectOutPhaseDirection(COMA, ampl)

	with root_open("D0PiPiPi.root", "RECREATE") as outFile:
		COMAhist = pyRootPwa.ROOT.TH2D(comaName, comaName, 2*nBins, 0.,1., 2*nBins, 0.,1.)
		for i in range(2*nBins):
#			COMAhist.SetBinContent(i+1, i+1, 1.) # first two are fixed isobar
			for j in range(2*nBins):
				COMAhist.SetBinContent(i+1, j+1, COMA[i+2, j+2]) # first two are fixed isobar
				pass
		COMAhist.Write()

		eigenHist = pyRootPwa.ROOT.TH1D("eigen0_0", "D0[pi,pi]0++PiS<|>D0[pi,pi]1--PiP", 50, .5 , 2.5)
		eigenHist.SetBinContent(bin+1, 0.000908109182548)
		eigenHist.Write()

		zeroHist = pyRootPwa.ROOT.TH2D("zero0_0","D0[pi,pi]0++PiS<|>D0[pi,pi]1--PiP", 50, .5 , 2.5, nBins, 0.,1.)
		for i in range(nBins):
			zeroHist.SetBinContent(bin+1, i+1, zeroMode[i].real)
		zeroHist.Write()

		for sector in sectors:
			if "1--" in sector:
				isobarBinning = binningRho
			elif "0++" in sector:
				isobarBinning = binningF0
			else:
				"Sector not valid: '"+sector+"'"
			histIntens = pyRootPwa.ROOT.TH2D(sector+"_0_intens", sector+"_0_intens", len(binning3Pi)-1, binning3Pi, len(isobarBinning)-1, isobarBinning)
			histReal   = pyRootPwa.ROOT.TH2D(sector+"_0_real"  , sector+"_0_real"  , len(binning3Pi)-1, binning3Pi, len(isobarBinning)-1, isobarBinning)
			histImag   = pyRootPwa.ROOT.TH2D(sector+"_0_imag"  , sector+"_0_imag"  , len(binning3Pi)-1, binning3Pi, len(isobarBinning)-1, isobarBinning)
			histNorm   = pyRootPwa.ROOT.TH2D(sector+"_0_norm"  , sector+"_0_norm"  , len(binning3Pi)-1, binning3Pi, len(isobarBinning)-1, isobarBinning)
			histIndex  = pyRootPwa.ROOT.TH2D(sector+"_0_index" , sector+"_0_index" , len(binning3Pi)-1, binning3Pi, len(isobarBinning)-1, isobarBinning)
			for b,B in enumerate(range(sectors[sector][0],sectors[sector][1])):
				histIntens.SetBinContent(bin+1, b+1, abs(ampl[B+1])**2)
				histReal.SetBinContent(  bin+1, b+1, ampl[B+1].real)
				histImag.SetBinContent(  bin+1, b+1, ampl[B+1].imag)
				histNorm.SetBinContent(  bin+1, b+1, norm[B])
				histIndex.SetBinContent( bin+1, b+1, B)
			histIntens.Write()
			histReal.Write()
			histImag.Write()
			histNorm.Write()
			histIndex.Write()


if __name__ == "__main__":
	main()
