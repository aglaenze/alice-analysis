#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>

#include "TSystem.h"
#include <TROOT.h>
#include <TMath.h>
#include "TString.h"
#include "TDatime.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TColor.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "RooWorkspace.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooChi2Var.h"

#include "../Include/FitUtils.C"

using namespace RooFit;
using namespace std;



void MakePlots(RooWorkspace* ws, string period, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax, bool drawPulls, bool logScale, bool exp, bool exclusiveOnly) {

    gStyle->SetTitleSize(0.04);
    gStyle->SetTitleFontSize(0.05);
    gStyle->SetTitleSize(0.04, "XY");
    gStyle->SetLabelSize(.04, "XY");
    
    // get what we need of the workspace and fit
    RooRealVar* m = ws->var("fTrkTrkM");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    
    // Define mass frame
    RooPlot* mframe = m->frame(Title("Plot of invariant mass"));
    int binNum = 50;
    if (exclusiveOnly) binNum = 30;
    data->plotOn(mframe, Binning(binNum));
    //data->plotOn(mframe);
    
    
    // for naming output files
    string cutType = "";
    if (!useCuts) cutType = "-nocuts";
    string suf;
    if (exp) suf = "exp";
    else suf = "powerlaw";
    if (exclusiveOnly) suf += "-exclusive-only";

    TCanvas* c2 = new TCanvas("plot","fit",600,500) ;

    // Draw mass plot
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;

    // Add numbers
    TLegend* lgd = new TLegend(0.6, 0.7, 0.9, 0.9);
    lgd->SetMargin(0.05);
    lgd->SetTextSize(0.04);
    lgd->AddEntry((TObject*)0, Form("N_{data} = %.0f", data->sumEntries()), "");

    mframe->Draw();
    lgd->Draw();
    
    c2->SaveAs(Form("Plots/%s/Mass-plot-2D%s-%.1f-%.1f-%s-%.1f-%.1f.pdf", period.c_str(), cutType.c_str(), mMin, mMax, suf.c_str(), abs(rapMax), abs(rapMin)));
}


void DataMass(string rootfilePath, string rootfilePathMC, vector<string> periods = {"LHC16r", "LHC16s"}) {
    
    Double_t mMin = 2., mMax = 3.5, ptMin = 0., ptMax = 3.5;
    Double_t rapMin = -4, rapMax = -2.5;
    bool useCuts = true, exp = false, exclusiveOnly = false;
    bool logScale = false, drawPulls = false;
    
    bool noInclusive = false;
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.07);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetTitleSize(.07);
    gStyle->SetTextSize(.07);
    gStyle->SetLabelSize(.05, "XY");
    //gStyle->SetMarkerSize(0.5);
    //gStyle->SetMarkerStyle(20);
    
    gROOT->ProcessLine(".L ../Include/ExtendedCrystalBall.cxx+") ;
    gSystem->Load("./../Include/ExtendedCrystalBall_cxx.so") ;
    
    const int nPeriod = periods.size();
    
    for (int k = 0; k<nPeriod; k++) {
        string period = periods[k];
        if ( ! Initiate(period, mMin, mMax, ptMin, ptMax, rapMin, rapMax, useCuts, logScale, drawPulls, exp, exclusiveOnly)) {cout << "Something wrong at initialisation"; return;}
        //return;
        if (exclusiveOnly) {
            ptMax = 1.2;
            /*
             ptMin = 0.2;
             mMin = 2.8;
             mMax = 3.3;
             */
        }
        else {
            if (period == "LHC16s") ptMax = 3;
        }
        
        // Define cuts
        std::list<TCut> mCutList = DefineCuts(period, exclusiveOnly);
        TCut mCut = "";
        
        for (std::list <TCut>::iterator iter = mCutList.begin(); iter != mCutList.end(); ++iter) {mCut += *iter;}
        if (!useCuts) mCut = "";
        
        // Open the file
        TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()),"READ");
        
        // Connect to the tree
        TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
        //Create a new workspace to manage the project
        RooWorkspace* wspace = new RooWorkspace("myJpsi");
        //return;
        ImportDataSet(wspace, fAnaTree, mCut, mMin, mMax, ptMin, ptMax, rapMin, rapMax);
        //wspace->Print();
        MakePlots(wspace, period, useCuts, mMin, mMax, ptMin, ptMax, rapMin, rapMax, drawPulls, logScale, exp, exclusiveOnly);
        
        return;
    }
}


