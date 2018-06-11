#ifndef GETBELLEDATA
#define GETBELLEDATA
#include<vector>
#include<string>
#include"TH2D.h"
void BELLE_apply_selection_final_Dstar0ToD0piplus_D0ToKspipi_massdiffsidebandlarge(std::string infilename, std::string outfilename);

std::vector<std::vector<double> > getBELLEevents(std::string inFileName, int SP_sign, bool SIGNSWITCH = true);

TH2D makeDalitzPlot(std::string name, const std::vector<std::vector<double> > data, double sMin = 0., double sMax = 4., double nBin = 200);
#endif//GETBELLEDATA
