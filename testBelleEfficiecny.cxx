#include <iostream>
#include <fstream>
#include "BELLE_efficiency.h"
#include <vector>
#include "TFile.h"
#include "TH2D.h"
int main() {
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
	return 0;
}
