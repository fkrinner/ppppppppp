#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include <string>

void BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge(std::string infilename, std::string outfilename) {
	TFile *file0 = new TFile(infilename.c_str(),"READ");
	TTree *tree0 = (TTree*)file0->Get("h1");

	TFile *filenew = new TFile(outfilename.c_str(),"RECREATE");

    	std::string cut_string_selection = "p_bestcand_random==1 && 0.1504<Dstar_massdiff && Dstar_massdiff<0.16 && abs(D0_mass-1.865)<0.015 && ( ( (e_expno==43 || e_expno==53 || e_expno==67 || e_expno==67 || e_expno==69 || e_expno==71) && Dstar_cm_ptot>3.1 ) || (Dstar_cm_ptot>2.5) )";

        TCut cut( cut_string_selection.c_str() );

	TTree *treenew = tree0->CopyTree( cut );
	
	treenew->Write();
	filenew->Close();
}


int main() {
	std::string inFileName = "/nfs/freenas/tuph/e18/data/belle/d0_kpipi/forDanielGreenwald/data/data_official_skim_KSPIPI.root";
	std::string outFileName = "./BELLE_data.root";
	BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge(inFileName, outFileName);
}
