#include <iostream>
#include <fstream>
#include "BELLE_efficiency.h"
#include <vector>
#include "TFile.h"
#include "TH2D.h"
#include "getBELLEdata.h"
int main() {
/*
	TFile* outFile = new TFile("effciency.root", "RECREATE");
	TH2D*  hist = new TH2D("eff", "eff", 200, 0., 4., 200., 0., 4.);
	for (int i = 0; i < 200; ++i) {
		double sx = hist->GetXaxis()->GetBinCenter(i+1);
		for (int j = 0; j < 200; ++j) {
			double sy = hist->GetYaxis()->GetBinCenter(j+1);
			std::vector<double> kin = {sx,sy};
			double eff = Efficiency(kin);
			hist->SetBinContent(i+1,j+1,eff);
		}
	}
	hist->Write();
	outFile->Close();
*/
	std::string inFileName  = "/nfs/freenas/tuph/e18/data/belle/d0_kpipi/forDanielGreenwald/data/data_official_skim_KSPIPI.root";
//	std::string outFileName = "./BELLE_data_eCut.root";
//	std::string outFileName = "./BELLE_data_antiEcut.root";
//	std::string belleData   = "./BELLE_data.root";
	std::string outFileNameSignal         = "./BELLE_data.root";
	std::string outFileNameDeltaMSideband = "./BELLE_deltaMsideband.root";
	std::string outFileNameMDsideband     = "./BELLE_mDsideband.root";
	std::string outFileNameBothSidebands  = "./BELLE_bothSidebands.root";

//	BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge_andElectronCut(inFileName, outFileName);

//	BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi(inFileName, outFileNameSignal        , false, false);
//	BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi(inFileName, outFileNameDeltaMSideband, true , false);
//	BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi(inFileName, outFileNameMDsideband    , false, true );
//	BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi(inFileName, outFileNameBothSidebands , true , true );

	makeIDplot(outFileNameSignal        , "idPlots_signal.root");
	makeIDplot(outFileNameDeltaMSideband, "idPlots_deltaMsideband.root");
	makeIDplot(outFileNameMDsideband    , "idPlots_mDsideband.root");
	makeIDplot(outFileNameBothSidebands , "idPlots_bothSidebands.root");
	return 0;
}
