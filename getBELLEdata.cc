#include "getBELLEdata.h"
#include "constants.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>

void BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge(std::string infilename, std::string outfilename) {
	TFile *file0 = new TFile(infilename.c_str(),"READ");
	TTree *tree0 = (TTree*)file0->Get("h1");

	TFile *filenew = new TFile(outfilename.c_str(),"RECREATE");

    	std::string cut_string_selection = "p_bestcand_random==1 && 0.1504<Dstar_massdiff && Dstar_massdiff<0.16 && abs(D0_mass-1.865)<0.015 && ( ( (e_expno==43 || e_expno==53 || e_expno==67 || e_expno==67 || e_expno==69 || e_expno==71) && Dstar_cm_ptot>3.1 ) || (Dstar_cm_ptot>2.5) )";

        TCut cut( cut_string_selection.c_str() );

	TTree *treenew = tree0->CopyTree( cut );
	
	treenew->Write();
	filenew->Close();
	file0->Close();
}

std::vector<std::vector<double> > getBELLEevents(std::string inFileName, int SP_sign, bool SIGNSWITCH) {
	if (SP_sign*SP_sign*SP_sign != SP_sign) {
		std::cout << "SP_sign (soft pion) can only be 1, -1, or 0 (for positive, negative, or both)" << std::endl;
		throw;
	}
	TFile *file0   = new TFile(inFileName.c_str(),"READ");
	TTree *treenew = (TTree*)file0->Get("h1");

	size_t nPoints = treenew->GetEntries();
	std::vector<std::vector<double> > retVal(nPoints, std::vector<double>(3,mD0*mD0));

	TBranch* branch = treenew->GetBranch("h1");

	TLeaf* L_m01 = branch-> GetLeaf("D0_massconstrained_mass_squared01");
	TLeaf* L_m02 = branch-> GetLeaf("D0_massconstrained_mass_squared02");
	TLeaf* L_m12 = branch-> GetLeaf("D0_massconstrained_mass_squared12");
	TLeaf* L_spc = branch-> GetLeaf("softpion_charge");

	size_t count = 0;
	for (int i = 0, N = treenew->GetEntries(); i < N; ++i) {
		treenew->GetEntry(i);
		if (L_m01->GetValue(0) == 0. || L_m02->GetValue(0) == 0.|| L_m12->GetValue(0) == 0.) {
			continue;
		}

		if ((L_spc->GetValue(0) - SP_sign)*SP_sign != 0.) {
			continue;
		}
		if (L_spc->GetValue(0) == 1. or not SIGNSWITCH) {
			retVal[count][1] = L_m01->GetValue(0);
			retVal[count][2] = L_m12->GetValue(0);
		} else {
			retVal[count][1] = L_m02->GetValue(0);
			retVal[count][2] = L_m12->GetValue(0);
		}
		++count;	
	}
	std::cout << count << " events have been found for SP_charge of " << SP_sign << std::endl;
	file0->Close();
	retVal.resize(count);
	return retVal;
}

TH2D makeDalitzPlot(std::string name, const std::vector<std::vector<double> > data, double sMin, double sMax, double nBin) {
	TH2D hist(name.c_str(), name.c_str(), nBin, sMin, sMax, nBin, sMin, sMax);
	for (size_t e = 0; e < data.size(); ++e) {
		hist.Fill(data[e][1], data[e][2]);
	}
	return hist;
}


/*
int main() {
	std::string inFileName = "/nfs/freenas/tuph/e18/data/belle/d0_kpipi/forDanielGreenwald/data/data_official_skim_KSPIPI.root";
	std::string outFileName = "./BELLE_data.root";
//	BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge(inFileName, outFileName);
	TFile* outFile = new TFile("./histograms.root", "RECREATE");
	TH2D hist = makeDalitzPlot("both_true" ,getEvents(outFileName,  0, true ));
	outFile->cd();
	hist.Write();
	hist = makeDalitzPlot("pos_true"  ,getEvents(outFileName,  1, true ));
	outFile->cd();
	hist.Write();
	hist = makeDalitzPlot("neg_true"  ,getEvents(outFileName, -1, true ));
	outFile->cd();
	hist.Write();
	hist = makeDalitzPlot("both_false",getEvents(outFileName,  0, false));
	outFile->cd();
	hist.Write();
	hist = makeDalitzPlot("pos_false" ,getEvents(outFileName,  1, false));
	outFile->cd();
	hist.Write();
	hist = makeDalitzPlot("neg_false" ,getEvents(outFileName, -1, false));
	outFile->cd();
	hist.Write();
	outFile->Close();
	return 0;
}
*/
