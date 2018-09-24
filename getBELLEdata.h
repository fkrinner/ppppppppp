#ifndef GETBELLEDATA
#define GETBELLEDATA
#include<vector>
#include<string>
#include"TH2D.h"
void BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge(std::string infilename, std::string outfilename);
void BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge_andElectronCut(std::string infilename, std::string outfilename);
void BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi(const std::string inFileName, const std::string outFileName, bool massDiffSideBand = false, bool mDsideBand = false, int mDsbRegion = 0);

std::vector<std::vector<double> > getBELLEevents(std::string inFileName, int SP_sign, bool SIGNSWITCH = true);

void makeDmassBinnedDalitzs(const std::string inFileName, const std::string outFileName, std::vector<std::pair<double, double> > binning, int SP_sign, bool SIGNSWITCH = true, int nBinsX = 100, int nBinsY = 100);

TH2D makeDalitzPlot(std::string name, const std::vector<std::vector<double> > data, double sMin = 0., double sMax = 4., double nBin = 200);

void makeIDplot(const std::string inFileName, const std::string outFileName);
#endif//GETBELLEDATA
