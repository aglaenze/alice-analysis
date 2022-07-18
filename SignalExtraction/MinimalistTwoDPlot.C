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

#include "Include/FitUtils.C"
#include "Include/Initiate.C"

using namespace RooFit;
using namespace std;

Double_t mLimitPsi2s = 3.65;


void AddModel(RooWorkspace* ws, string rootfilePath, string rootfilePathMC, string period, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax, bool exp, bool exclusiveOnly, bool noInclusive = true) {
    // Define 2D model
    // First define fits in mass
    // Second define fits in pt
    // Then make the product to get 2D PDFs
    
    LoadMassFitFunctions(ws, period);
    LoadJpsiPtFitFunctions(ws, rootfilePathMC, period, ptMin, ptMax, exp, exclusiveOnly);
    LoadBkgPtFitFunctions(ws, rootfilePath, rootfilePathMC, period, mMin, mMax, ptMin, ptMax, exclusiveOnly);
    
    // Now collect m pdf
    RooAbsPdf* jpsi = ws->pdf("jpsi");
    RooAbsPdf* psi2s = ws->pdf("psi2s");
    RooAbsPdf* bkg = ws->pdf("bkg");
    
    // Then pt pdf
    RooAbsPdf* ptJpsiExclusive = ws->pdf("ptJpsiExclusive");
    RooAbsPdf* ptJpsiDissociative = ws->pdf("ptJpsiDissociative");
    RooAbsPdf* ptJpsiGammaPb = ws->pdf("ptJpsiGammaPb");
    RooAbsPdf* ptJpsiInclusive = ws->pdf("ptJpsiInclusive");
    RooAbsPdf* ptPsi2s = ws->pdf("ptPsi2s");
    RooAbsPdf* ptTwoGamma = ws->pdf("ptTwoGamma");
    RooAbsPdf* ptBackground = ws->pdf("ptBackground");
    
    // Last step:
    // Product fit(m) x fit(pt)
    RooProdPdf* pdfJpsiExclusive = new RooProdPdf("pdfJpsiExclusive","jpsi*ptJpsiExclusive",RooArgList(*jpsi,*ptJpsiExclusive));
    RooProdPdf* pdfJpsiDissociative = new RooProdPdf("pdfJpsiDissociative","jpsi*ptJpsiDissociative",RooArgList(*jpsi,*ptJpsiDissociative));
    RooProdPdf* pdfJpsiGammaPb = new RooProdPdf("pdfJpsiGammaPb","jpsi*ptJpsiGammaPb",RooArgList(*jpsi,*ptJpsiGammaPb));
    RooProdPdf* pdfJpsiInclusive = new RooProdPdf("pdfJpsiInclusive","jpsi*ptJpsiInclusive",RooArgList(*jpsi,*ptJpsiInclusive));
    RooProdPdf* pdfPsi2s = new RooProdPdf("pdfPsi2s","psi2s*ptPsi2s",RooArgList(*psi2s,*ptPsi2s));
    RooProdPdf* pdfTwoGamma = new RooProdPdf("pdfTwoGamma","bkg*ptTwoGamma",RooArgList(*bkg,*ptTwoGamma));
    RooProdPdf* pdfBackground = new RooProdPdf("pdfBackground","bkg*ptBackground",RooArgList(*bkg,*ptBackground));
    
    // Load yields
    LoadYields(ws, period, rapMin, rapMax, mMin, mMax, ptMin, ptMax, exclusiveOnly);
    RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
    RooRealVar* yieldJpsiDissociative = ws->var("yieldJpsiDissociative");
    RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
    RooRealVar* yieldJpsiInclusive = ws->var("yieldJpsiInclusive");
    RooRealVar* yieldPsi2s = ws->var("yieldPsi2s");
    RooRealVar* yieldTwoGamma = ws->var("yieldTwoGamma");
    RooRealVar* yieldBkg = ws->var("yieldBkg");
    
    // Assemble all components in sets
    RooArgList* pdfList = new RooArgList(*pdfJpsiExclusive, *pdfJpsiGammaPb, *pdfTwoGamma, *pdfBackground);
    RooArgList* yieldList = new RooArgList(*yieldJpsiExclusive, *yieldJpsiGammaPb, *yieldTwoGamma, *yieldBkg);
    
    if (!exclusiveOnly) {
        pdfList->add(*pdfJpsiDissociative);
        yieldList->add(*yieldJpsiDissociative);
    }
    
    // Create fit model
    RooAbsPdf* fitModel = new RooAddPdf("model", "model", *pdfList, *yieldList, kFALSE);
    
    ws->import(*fitModel);
}

void MakePlots(RooWorkspace* ws, string period, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax, bool drawPulls, bool logScale, bool exp, bool exclusiveOnly, double& chi2, bool &noInclusive) {
    
    int ptBinNumber = int(15*(ptMax-ptMin));
    
    //get what we need of the workspace
    RooRealVar* m = ws->var("fTrkTrkM");
    RooRealVar* pt = ws->var("fTrkTrkPt");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    RooAbsPdf* fitModel = ws->pdf("model");
    
    // Fit data
    RooFitResult* r = fitModel->fitTo(*data, Extended(), Minos(true), Strategy(1), Save());
    //RooFitResult* r = fitModel->fitTo(*data, Minos(false), Strategy(0), Save());    // be quick (for testing)
    //RooFitResult* r = nullptr;
    
    RooAbsPdf* pdfJpsiExclusive = ws->pdf("pdfJpsiExclusive");
    RooAbsPdf* pdfJpsiDissociative = ws->pdf("pdfJpsiDissociative");
    RooAbsPdf* pdfJpsiGammaPb = ws->pdf("pdfJpsiGammaPb");
    RooAbsPdf* pdfTwoGamma = ws->pdf("pdfTwoGamma");
    RooAbsPdf* pdfBackground = ws->pdf("pdfBackground");
    RooAbsPdf* pdfJpsiInclusive = ws->pdf("pdfJpsiInclusive");
    RooAbsPdf* pdfPsi2s = ws->pdf("pdfPsi2s");
    
    RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
    RooRealVar* yieldJpsiDissociative = ws->var("yieldJpsiDissociative");
    RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
    RooRealVar* yieldTwoGamma = ws->var("yieldTwoGamma");
    RooRealVar* yieldBkg = ws->var("yieldBkg");
    RooRealVar* yieldJpsiInclusive = ws->var("yieldJpsiInclusive");
    RooRealVar* yieldPsi2s = ws->var("yieldPsi2s");
    

    RooRealVar* bExc = ws->var("bExc");
    RooRealVar* bDiss = ws->var("bDiss");
    RooRealVar* nDiss = nullptr;
    if (!exp) nDiss = ws->var("nDiss");
    
    RooRealVar* meanGg = ws->var("mean_gg");
    RooRealVar* sigmaGg = ws->var("sigma_gg");
    
    // Define mass frame
    RooPlot* mframe = m->frame(Title("Fit of invariant mass"));
    int mbins = int((mMax-mMin)*50);
    data->plotOn(mframe, Binning(mbins));
    fitModel->plotOn(mframe, Name("sum"), LineColor(kRed), LineWidth(2));
    fitModel->plotOn(mframe,Name("pdfJpsiExclusive"),Components(*pdfJpsiExclusive),LineStyle(kDashed), LineColor(kOrange+8), LineWidth(2));
    if (!exclusiveOnly) { fitModel->plotOn(mframe,Name("pdfJpsiDissociative"),Components(*pdfJpsiDissociative),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(2));
    }
    fitModel->plotOn(mframe,Name("pdfBackground"),Components(*pdfBackground),LineStyle(kDashed), LineColor(kBlack), LineWidth(2));
    fitModel->plotOn(mframe,Name("pdfJpsiGammaPb"),Components(*pdfJpsiGammaPb),LineStyle(kDashed), LineColor(4), LineWidth(2));
    if (mMax > mLimitPsi2s) fitModel->plotOn(mframe,Name("pdfPsi2s"),Components(*pdfPsi2s),LineStyle(kDashed), LineColor(6), LineWidth(2));
    fitModel->plotOn(mframe,Name("pdfTwoGamma"),Components(*pdfTwoGamma),LineStyle(kDashed), LineColor(kMagenta), LineWidth(2));
    
    if (!noInclusive) fitModel->plotOn(mframe,Name("pdfJpsiInclusive"),Components(*pdfJpsiInclusive),LineStyle(kDashed), LineColor(kMagenta+2), LineWidth(2));
    
    Double_t ptRangeMax = ptMax;
    if (!exclusiveOnly) ptRangeMax = 3;
    // Define pt frame
    vector<RooPlot*> ptframes= {pt->frame(Title("Fit of p_{T} (log scale and full p_{T} range)")), pt->frame(Title("Fit of p_{T} (zoom)"))};
    for (int k = 0; k<(int)ptframes.size(); k++) {
    //for (int k = 0; k< 1; k++) {
        if (exclusiveOnly) ptframes[k]->SetTitle("Fit of p_{T}");
        data->plotOn(ptframes[k], Binning(ptBinNumber));
        fitModel->plotOn(ptframes[k], Name("sum"), LineColor(kRed), LineWidth(2));
        fitModel->plotOn(ptframes[k],Name("pdfJpsiExclusive"),Components(*pdfJpsiExclusive),LineStyle(kDashed), LineColor(kOrange+8), LineWidth(2));
        if (!exclusiveOnly) fitModel->plotOn(ptframes[k],Name("pdfJpsiDissociative"),Components(*pdfJpsiDissociative),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(2));
        fitModel->plotOn(ptframes[k],Name("pdfJpsiGammaPb"),Components(*pdfJpsiGammaPb),LineStyle(kDashed), LineColor(4), LineWidth(2));
        if (mMax > mLimitPsi2s) fitModel->plotOn(ptframes[k],Name("pdfPsi2s"),Components(*pdfPsi2s),LineStyle(kDashed), LineColor(6), LineWidth(2));
        fitModel->plotOn(ptframes[k],Name("pdfTwoGamma"),Components(*pdfTwoGamma),LineStyle(kDashed), LineColor(kMagenta), LineWidth(2));
        fitModel->plotOn(ptframes[k],Name("pdfBackground"),Components(*pdfBackground),LineStyle(kDashed), LineColor(kBlack), LineWidth(2));
        if (!noInclusive) {
            fitModel->plotOn(ptframes[k],Name("pdfJpsiInclusive"),Components(*pdfJpsiInclusive),LineStyle(kDashed), LineColor(kMagenta+2), LineWidth(2));
        }
        ptframes[k]->GetXaxis()->SetRangeUser(0, ptRangeMax);
    }
    //mframe->GetXaxis()->SetRangeUser(0, ptRangeMax);
    
    TCanvas* cv = new TCanvas("cv", "cv", 600, 300);
    cv->Divide(2);
    cv->cd(1);
    mframe->Draw();
    cv->cd(2);
    ptframes[0]->Draw();
    
    cv->SaveAs("Plots/LHC16r/test.pdf");
    return;
    
    
    
}

void TwoDPlot(string rootfilePath, string rootfilePathMC, vector<string> periods = {"LHC16r", "LHC16s"}) {
    
    Double_t mMin = 2., mMax = 3.5, ptMin = 0., ptMax = 3.5;
    Double_t rapMin = -4, rapMax = -2.5;
    bool useCuts = true, exp = false, exclusiveOnly = false;
    bool logScale = false, drawPulls = false;
    
    bool noInclusive = true;
    
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.07);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetTitleSize(.07);
    gStyle->SetTextSize(.07);
    gStyle->SetLabelSize(.05, "XY");
    //gStyle->SetMarkerSize(0.5);
    //gStyle->SetMarkerStyle(20);
    
    gROOT->ProcessLine(".L Include/ExtendedCrystalBall.cxx+") ;
    gSystem->Load("./Include/ExtendedCrystalBall_cxx.so") ;
    
    const int nPeriod = periods.size();
    
    for (int k = 0; k<nPeriod; k++) {
        string period = periods[k];
        if ( ! Initiate(period, mMin, mMax, ptMin, ptMax, rapMin, rapMax, useCuts, logScale, drawPulls, exp, exclusiveOnly)) {cout << "Something wrong at initialisation"; return;}
        
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
        ImportDataSet(wspace, fAnaTree, mCut, mMin, mMax, ptMin, ptMax, rapMin, rapMax);
        AddModel(wspace, rootfilePath, rootfilePathMC, period, mMin, mMax, ptMin, ptMax, rapMin, rapMax, exp, exclusiveOnly, noInclusive);
        wspace->Print();
        double chi2;
        MakePlots(wspace, period, useCuts, mMin, mMax, ptMin, ptMax, rapMin, rapMax, drawPulls, logScale, exp, exclusiveOnly, chi2, noInclusive);

    }
}


