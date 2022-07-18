#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>

#include "TString.h"
#include "TDatime.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"
#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>
#include <TCut.h>
#include <TTree.h>
#include <TLegend.h>

using namespace std;



int DrawEffFunctionOfMass() {
    
    const int n = 5;

    Double_t massMean[n] = {1.25, 1.75, 2.25, 2.75, 3.25};
    Double_t massInterval[n] = {0.25, 0.25, 0.25, 0.25, 0.25};
    
    Double_t accEffValues[n] = {0.0165553, 0.0304076, 0.0409487, 0.0459062, 0.048367};
    Double_t accEffErrors[n] = {0.000135916, 0.000318241, 0.000559454, 0.000835111, 0.00115493};
    
    TGraphErrors* gr = new TGraphErrors(n, massMean, accEffValues, massInterval, accEffErrors);
    gr->SetTitle("A x Eff = f(M_{#mu#mu})");
    gr->GetHistogram()->GetXaxis()->SetTitle("M_{#mu#mu} (GeV)");
    gr->GetHistogram()->GetYaxis()->SetTitle("A x Eff");
    
    gr->GetHistogram()->GetYaxis()->SetRangeUser(0,0.06);
    
    TCanvas* cv = new TCanvas();
    gr->Draw("AP");
    cv->SaveAs("Plots/EffVsMass.pdf");
    return 0;
}


