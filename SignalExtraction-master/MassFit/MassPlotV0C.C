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

Double_t mLimitPsi2s = 3.65;



void AddModel(RooWorkspace* ws, string rootfilePath, string rootfilePathMC, string period, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax, bool exp, bool exclusiveOnly) {
    // Define 2D model
    // First define fits in mass
    // Second define fits in pt
    // Then make the product to get 2D PDFs
    
    LoadMassFitFunctions(ws, period);
    
    // Now collect m pdf
    RooAbsPdf* pdfJpsi = ws->pdf("jpsi");
    pdfJpsi->SetName("pdfJpsi");
    RooAbsPdf* pdfTwoGamma = ws->pdf("bkg");
    pdfTwoGamma->SetName("pdfTwoGamma");
    
    // Load yields
    RooRealVar* yieldJpsi = new RooRealVar("yieldJpsi", "yieldJpsi", 1000, 10, 5000);
    RooRealVar* yieldTwoGamma = new RooRealVar("yieldTwoGamma", "yieldTwoGamma", 1000, 10, 5000);
    
    // Assemble all components in sets
    RooArgList* pdfList = new RooArgList(*pdfJpsi, *pdfTwoGamma);
    RooArgList* yieldList = new RooArgList(*yieldJpsi, *yieldTwoGamma);
    
    // Create fit model
    RooAbsPdf* fitModel = new RooAddPdf("model", "model", *pdfList, *yieldList, kFALSE);
    
    ws->import(*fitModel);
}


void MakePlots(RooWorkspace* ws, string period, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax, bool drawPulls, bool logScale, bool exp, bool exclusiveOnly, int i) {
    
    //int ptBinNumber = int(15*(ptMax-ptMin));
    //int ptBinNumber = int(10*(ptMax-ptMin));
    //int mBinNumber = int(20*(mMax-mMin));
    //int ptBinNumber = 20;
    
    
    // get what we need of the workspace and fit
    RooRealVar* m = ws->var("fTrkTrkM");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    RooAbsPdf* fitModel = ws->pdf("model");
    
    
    // Fit data
    RooFitResult* r = fitModel->fitTo(*data, Extended(), Minos(true), Strategy(1), Save());
    //RooFitResult* r = fitModel->fitTo(*h_data, Extended(), Minos(true), Strategy(1), Save());
    //RooFitResult* r = fitModel->fitTo(*data, Minos(false), Strategy(0), Save());    // be quick (for testing)
    //RooFitResult* r = nullptr;
    
    
    RooAbsPdf* pdfJpsi = ws->pdf("pdfJpsi");
    RooAbsPdf* pdfTwoGamma = ws->pdf("pdfTwoGamma");
    
    RooRealVar* yieldJpsi = ws->var("yieldJpsi");
    RooRealVar* yieldTwoGamma = ws->var("yieldTwoGamma");
    
    // Define mass frame
    RooPlot* mframe = m->frame(Title("Fit of invariant mass"));
    int binNum = 50;
    if (exclusiveOnly) binNum = 30;
    data->plotOn(mframe, Binning(binNum));
    //data->plotOn(mframe);
    fitModel->plotOn(mframe, Name("model"), LineColor(kRed), LineWidth(2));
    fitModel->plotOn(mframe,Name("pdfJpsi"),Components(*pdfJpsi),LineStyle(kDashed), LineColor(kOrange+8), LineWidth(2));
    fitModel->plotOn(mframe,Name("pdfTwoGamma"),Components(*pdfTwoGamma),LineStyle(kDashed), LineColor(kMagenta), LineWidth(2));
    
    
    // for naming output files
    string cutType = "";
    if (!useCuts) cutType = "-nocuts";
    string suf;
    if (exp) suf = "exp";
    else suf = "powerlaw";
    if (exclusiveOnly) suf += "-exclusive-only";
    
    
    
    
    TLegend* legend2;
    if (exclusiveOnly) legend2 = new TLegend(0.48, 0.7, 1, 1);
    else legend2 = new TLegend(0.55, 0.65, 1, 1);
    legend2->SetFillColor(kWhite);
    legend2->SetLineColor(kWhite);
    legend2->AddEntry(mframe->findObject("pdfJpsi"), "J/#Psi","L");
    legend2->AddEntry(mframe->findObject("pdfTwoGamma"), "#gamma#gamma","L");
    legend2->AddEntry(mframe->findObject("model"),"sum","L");
    
    TCanvas* c2 = new TCanvas("2Dplot","2D fit",600,500) ;
    // Draw mass plot
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    
    
    // Add number of J/Psis
    TLatex* txtJpsi = new TLatex(2.6, 0.75*mframe->GetMaximum(), Form("#bf{N_{J/#Psi} = %.0f #pm %.0f}", yieldJpsi->getVal(), yieldJpsi->getError()));
    txtJpsi->SetTextSize(0.045);
    mframe->addObject(txtJpsi) ;
    TLatex* txtJpsi3 = new TLatex(2.6, 0.65*mframe->GetMaximum(), Form("#bf{N_{#gamma#gamma} = %.0f #pm %.0f}", yieldTwoGamma->getVal(), yieldTwoGamma->getError()));
    txtJpsi3->SetTextSize(0.045);
    mframe->addObject(txtJpsi3) ;
    mframe->SetTitle("");
    mframe->SetMaximum(mframe->GetMaximum()*1.1);
    mframe->Draw();
    
    
    c2->SaveAs(Form("Plots/%s/V0C-%d/Fit-simple-2D%s-%.1f-%.1f-%s-%.1f-%.1f.pdf", period.c_str(), i, cutType.c_str(), mMin, mMax, suf.c_str(), abs(rapMax), abs(rapMin)));
}

void GetYields(RooWorkspace* ws, string period, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax, bool drawPulls, bool logScale, bool exp, bool exclusiveOnly, double &jpsiNum, double &ggNum) {
    
    //int ptBinNumber = int(15*(ptMax-ptMin));
    //int ptBinNumber = int(10*(ptMax-ptMin));
    //int mBinNumber = int(20*(mMax-mMin));
    //int ptBinNumber = 20;
    
    
    // get what we need of the workspace and fit
    RooRealVar* m = ws->var("fTrkTrkM");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    RooAbsPdf* fitModel = ws->pdf("model");
    
    
    // Fit data
    RooFitResult* r = fitModel->fitTo(*data, Extended(), Minos(true), Strategy(1), Save());
    //RooFitResult* r = fitModel->fitTo(*h_data, Extended(), Minos(true), Strategy(1), Save());
    //RooFitResult* r = fitModel->fitTo(*data, Minos(false), Strategy(0), Save());    // be quick (for testing)
    //RooFitResult* r = nullptr;
    
    
    RooAbsPdf* pdfJpsi = ws->pdf("pdfJpsi");
    RooAbsPdf* pdfTwoGamma = ws->pdf("pdfTwoGamma");
    
    RooRealVar* yieldJpsi = ws->var("yieldJpsi");
    RooRealVar* yieldTwoGamma = ws->var("yieldTwoGamma");
    jpsiNum = yieldJpsi->getVal();
    ggNum = yieldTwoGamma->getVal();
    return;
}

void MassPlotV0C(string rootfilePath, string rootfilePathMC, vector<string> periods = {"LHC16r", "LHC16s"}) {
    
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
        //ptMin = 0.5;
        //ptMax = 2;
        int nBins = 20;
        TH1F* hist = new TH1F("hist", "Number of J/#psi with nV0C < nMatched + N", nBins+1, -0.5, nBins+0.5);
        TH1F* hist2 = new TH1F("hist2", "Number of #gamma#gamma with nV0C < nMatched + N", nBins+1, -0.5, nBins+0.5);
        for (int i = 0; i<nBins+1; i++) {
            // Define cuts
            std::list<TCut> mCutList = DefineCutsBis(period, exclusiveOnly);
            TCut mCut = "";
            
            for (std::list <TCut>::iterator iter = mCutList.begin(); iter != mCutList.end(); ++iter) {mCut += *iter;}
            mCut += Form("fNV0C < fV0CNMatched + %d", i+1);
            if (!useCuts) mCut = "";
            
            // Open the file
            TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()),"READ");
            
            // Connect to the tree
            TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
            //Create a new workspace to manage the project
            RooWorkspace* wspace = new RooWorkspace("myJpsi");
            //return;
            ImportDataSet(wspace, fAnaTree, mCut, mMin, mMax, ptMin, ptMax, rapMin, rapMax);
            AddModel(wspace, rootfilePath, rootfilePathMC, period, mMin, mMax, ptMin, ptMax, rapMin, rapMax, exp, exclusiveOnly);
            //return;
            //wspace->Print();
            MakePlots(wspace, period, useCuts, mMin, mMax, ptMin, ptMax, rapMin, rapMax, drawPulls, logScale, exp, exclusiveOnly, i);
            double jpsiNum = 0, ggNum = 0;
            GetYields(wspace, period, useCuts, mMin, mMax, ptMin, ptMax, rapMin, rapMax, drawPulls, logScale, exp, exclusiveOnly, jpsiNum, ggNum);
            hist->Fill(i, jpsiNum);
            hist2->Fill(i, ggNum);
        
        }
        hist->SetMinimum(0);
        //hist->SetMaximum(4000);
        hist->SetXTitle("nV0C <= nMatched + N");
        hist2->SetMinimum(0);
        //hist2->SetMaximum(4000);
        hist2->SetXTitle("nV0C <= nMatched + N");
        
        // Draw
        TCanvas* cv = new TCanvas("cv", "cv", 800, 400);
        cv->Divide(2);
        cv->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        hist->Draw("hist");
        cv->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        hist2->Draw("hist");
        
        string suffix = "";
        if (exclusiveOnly) suffix = "-exc";
        cv->SaveAs(Form("Plots/%s/Numbers-V0C%s.pdf", period.c_str(), suffix.c_str()));
        
        return;
    }
}


