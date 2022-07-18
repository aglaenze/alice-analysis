#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>

#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TDatime.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"


using namespace std;

int DrawXSections() {
    
    const int n = 3;
    
    Double_t massCenter[n] = {1.25, 1.75, 2.25};
    Double_t massInterval[n] = {0.25, 0.25, 0.25};
    
    Double_t rapCenter[n] = {3.25, 2.875, 3.625};
    Double_t rapInterval[n] = {0.25, 0.25, 0.25};
    
    // Cross sections
    Double_t xSec1[n] = {5100, 2980, 2210};
    Double_t xSecError1[n] = {530, 460, 220};
    
    Double_t xSec2[n] = {1960, 1110, 816};
    Double_t xSecError2[n] = {190, 140, 86};
    
    Double_t xSec3[n] = {638, 357, 265};
    Double_t xSecError3[n] = {115, 105, 48};
    
    TGraphErrors* gr1 = new TGraphErrors(n, rapCenter, xSec1, rapInterval, xSecError1);
    TGraphErrors* gr2 = new TGraphErrors(n, rapCenter, xSec2, rapInterval, xSecError2);
    TGraphErrors* gr3 = new TGraphErrors(n, rapCenter, xSec3, rapInterval, xSecError3);
    
    TCanvas* cv = new TCanvas("cv", "cv", 300, 300);
    
    gr1->Draw("AP");
    gr2->Draw("P same");
    gr3->Draw("P same");
    
    cv->SaveAs("Plots/xSections/pdf");
    
    return 0;
}
