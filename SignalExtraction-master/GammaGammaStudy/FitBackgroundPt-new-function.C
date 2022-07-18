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

#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPlot.h"

#include "../Include/Initiate.C"
#include "../Include/Utils.C"

using namespace RooFit;
using namespace std;

void GetPtHistMC(RooWorkspace* ws, string rootfilePathMC, string period, TString process, Double_t ptMin = 0, Double_t ptMax = 4., TCut cutMc = ""){
    TFile *fSimu = new TFile(Form("%s/AnalysisResults_%s_MC_", rootfilePathMC.c_str(), period.c_str()) + process + ".root","READ");
    TTree* fAnaTree = (TTree*)fSimu->Get("MyTask/fAnaTree");
    
    RooRealVar pt = *ws->var("fTrkTrkPt");
    Int_t nBins = int((ptMax-ptMin)*100);
    TH1F* histPt = new TH1F("hPt"+process, "hPt"+process, nBins, ptMin, ptMax);
    fAnaTree->Draw("fTrkTrkPt>>hPt"+process, cutMc);
    RooDataHist* dTemplatePt = new RooDataHist("ptHist"+process,"ptHist"+process, RooArgList(pt),histPt);
    RooHistPdf* ptPdf = new RooHistPdf("pt"+process, "pt"+process, pt, *dTemplatePt);
    ws->import(*ptPdf);
}


//____________________________________
void MakePlots(RooWorkspace *ws, string rootfilePath, string rootfilePathMc, string period, bool exclusiveOnly, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax){
    
    bool plot = true;
    int ptBinNumber = int(10*(ptMax-ptMin));
    
    string suffix = "";
    if (exclusiveOnly) suffix = "-exclusive";
    
    bool twoComp = false;
    
    if (!exclusiveOnly) twoComp = true;
    
    
    //get what we need of the workspace
    RooRealVar* pt = ws->var("fTrkTrkPt");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    
    RooRealVar* yieldTwoGamma = new RooRealVar("yieldTwoGamma","yieldTwoGamma", 300, 0., 2.e3);
    RooRealVar* yieldBkg = new RooRealVar("yieldBkg","yieldBkg", 100, 0., 2.e3);
    
    /*
     // With MC
     TCut cutMc = Form("fTrkTrkM > %f && fTrkTrkM < %f", mMin, mMax);
     GetPtHistMC(ws, rootfilePathMc, period, "kTwoGammaToMuLow", ptMin, ptMax, cutMc);
     RooAbsPdf* ptTwoGamma = ws->pdf("ptkTwoGammaToMuLow");
     ptTwoGamma->SetName("ptTwoGamma");
     */
    // With Landau
    /*
     RooRealVar mean("mean","mean",0.06, 0.01,0.1);
     RooRealVar sigma("sigma","sigma", 0.03, 0.01, 0.1); //3.43,3.73);
     */
    // Landau
    double muValue = 0, sigmaValue = 0;
    double pt0Value = 0, nIncValue = 0, bIncValue = 0;
    //if (! GetParametersLandau(period, muValue, sigmaValue, pt0Value, nIncValue) ) {cout << "Parameters not loaded" << endl; return;}
    if (! GetParametersLandau(period, muValue, sigmaValue, bIncValue, nIncValue) ) {cout << "Parameters not loaded" << endl; return;}
    
    RooRealVar* mean = new RooRealVar("mean_gg","mean_gg", muValue, 0.04, 0.2);
    RooRealVar* sigma = new RooRealVar("sigma_gg","sigma_gg", sigmaValue, 0.005, 0.1); //3.43,3.73);
    
    if (!exclusiveOnly) {
        mean->setVal(muValue);
        mean->setConstant();
        sigma->setVal(sigmaValue);
        sigma->setConstant();
    }
    // m.setConstant(kTRUE);
    // sigma->setConstant(kTRUE);
    RooRealVar *bGg = new RooRealVar("bGg","bGg", 100, 0, 200);
    RooRealVar *nGg = new RooRealVar("nGg","nGg", 2, 0, 10);
    if (!exclusiveOnly) {
        bGg->setVal(200);
        bGg->setConstant();
        nGg->setVal(1.83431);
        nGg->setConstant();
    }
    RooAbsPdf* ptTwoGamma = new RooGenericPdf("ptTwoGamma","ptTwoGamma","(2*fTrkTrkPt*(1.+(fTrkTrkPt**2)*(bGg/nGg))**(-nGg))", RooArgSet(*pt, *bGg, *nGg));
    

     RooRealVar *bInc = new RooRealVar("bInc","bInc", 4, 1, 6);
     RooRealVar *nInc = new RooRealVar("nInc","nInc", 2, 1, 6);

    RooGenericPdf* ptBkg = new RooGenericPdf("ptBackground","ptBackground","(2*fTrkTrkPt*(1.+(fTrkTrkPt**2)*(bInc/nInc))**(-nInc))", RooArgSet(*pt, *bInc, *nInc));
    
    
    RooAbsPdf* model;
    Double_t fitRangeMax = ptMax;
    if (twoComp) {
        model = new RooAddPdf("fit", "fit", RooArgList(*ptTwoGamma, *ptBkg), RooArgList(*yieldTwoGamma, *yieldBkg), kFALSE);
        model->fitTo(*data, Extended(), Minos(true), Strategy(1));
    }
    else {
        model = new RooAddPdf("fit", "fit", RooArgList(*ptTwoGamma), RooArgList(*yieldTwoGamma), kFALSE);
        fitRangeMax = 0.38;
        model->fitTo(*data, Range(0, fitRangeMax), Extended(), Minos(true), Strategy(1));
    }
    
    
    if (plot) {
        RooPlot* frame = pt->frame();
        int nBinsPt = 100;
        if (exclusiveOnly) nBinsPt = 7;
        Double_t nBinsPt2 = int(ptMax/fitRangeMax*nBinsPt);
        data->plotOn(frame, Binning(nBinsPt2));
        //data->plotOn(frame);
        
        model->plotOn(frame);
        model->plotOn(frame, Components(*model), Name("sum"), LineStyle(kDashed), LineColor(kRed));
        model->plotOn(frame, Components(*ptTwoGamma), Name("ptTwoGamma"), LineStyle(kDashed), LineColor(kBlue));
        if (twoComp) model->plotOn(frame, Components(*ptBkg), Name("ptBkg"), LineStyle(kDashed), LineColor(kGreen));
        
        double textSize = 0.04;
        double xPos = ptMin+(1-ptMin)*0.4;
        double yMax0 = frame->GetMaximum();
        
        // write parameters info
        TLegend* lgdParam = new TLegend(0.65, 0.55, 0.9, 0.75);
        lgdParam->SetMargin(0.05);
        lgdParam->SetTextSize(0.03);
        lgdParam->AddEntry((TObject*)0,Form("%.1f #pm %.1f #gamma#gamma", yieldTwoGamma->getVal(), yieldTwoGamma->getError()), "");
        
        lgdParam->AddEntry((TObject*)0, Form("b_{gg} = %.2f #pm %.2f", bGg->getVal(), bGg->getError()), "");
        lgdParam->AddEntry((TObject*)0, Form("n_{gg} = %.2f #pm %.2f", nGg->getVal(), nGg->getError()), "");
        if (twoComp) {
            lgdParam->AddEntry((TObject*)0, Form("%.1f #pm %.1f bkg", yieldBkg->getVal(), yieldBkg->getError()), "");
            lgdParam->AddEntry((TObject*)0, Form("b_{inc} = %.1f #pm %.1f", bInc->getVal(), bInc->getError() ), "");
            lgdParam->AddEntry((TObject*)0, Form("n_{inc} = %.1f #pm %.1f", nInc->getVal(), nInc->getError()), "");
        }
        
        //make some plots
        TCanvas* cv = new TCanvas() ;
        gPad->SetLeftMargin(0.15) ;
        gPad->SetBottomMargin(0.15) ;
        //gPad->SetLogy();
        
        // Chi2
        // Chi2 will be calculated until pT = 1 so that there would be no empty bin
        Double_t ptMaxchi2 = 1;
        Int_t nBinsChi2 = int(ptMaxchi2/fitRangeMax*nBinsPt);
        TH1* histBkg = data->createHistogram("hPtBackground",*pt, Binning(nBinsChi2, 0, ptMaxchi2));
        RooDataHist* roo_hist_bkg = new RooDataHist("roo_hist", "roo_hist", RooArgSet(*pt), histBkg);
        TH1* histBkgModel = model->createHistogram("hPtBackgroundModel",*pt, Binning(nBinsChi2, 0, ptMaxchi2));
        RooDataHist* roo_hist_bkg_model = new RooDataHist("roo_hist_model", "roo_hist_model", RooArgSet(*pt), histBkgModel);
        RooHistPdf* roo_hist_pdf = new RooHistPdf("roo_hist_pdf", "roo_hist_pdf", RooArgSet(*pt), *roo_hist_bkg_model);
        int nFitParam = model->getParameters(data)->selectByAttrib("Constant",kFALSE)->getSize();
        //int nDof = nBinsPt - nFitParam;
        int nDof = nBinsChi2 - nFitParam;
        //RooChi2Var* chi2Var2 = new RooChi2Var("chi2","chi2", *model, *roo_hist_bkg, DataError(RooAbsData::Poisson));
        RooChi2Var* chi2Var2 = new RooChi2Var("chi2","chi2", *roo_hist_pdf, *roo_hist_bkg, DataError(RooAbsData::Poisson));
        //Double_t chi2M = frame->chiSquare( "sum", "h_data", nFitParam);
        Double_t chi2M = chi2Var2->getVal();
        //TLatex* txtChi = new TLatex(xPos*1.2,0.5*yMax0,Form("#bf{#chi^{2}/ndf = %.2f / %d = %.2f}", chi2M*nFitParam, nFitParam, chi2M) );
        TLatex* txtChi = new TLatex(xPos,0.5*yMax0,Form("#bf{#chi^{2}/ndf = %.1f/%d = %.2f}", chi2M, nDof, chi2M/nDof) );
        //frame->addObject(txtChi) ;
        frame->chiSquare() ;
        
        // Other draw options
        frame->SetTitle("");
        frame->Draw();
        
        TLegend* lgd;
        if (exclusiveOnly) lgd = new TLegend(0.73, 0.8, 0.9, 0.9);
        else lgd = new TLegend(0.7, 0.75, 0.9, 0.9);
        lgd->SetTextSize(0.045);
        lgd->AddEntry(frame->findObject("ptTwoGamma"), "#gamma #gamma", "L");
        if (twoComp) {
            lgd->AddEntry(frame->findObject("ptBkg"), "other bkg", "L");
            lgd->AddEntry(frame->findObject("sum"), "sum", "L");
        }
        lgd->Draw();
        
        lgdParam->Draw();
        
        bool eps = true;
        
        // write y and pt info
        TLegend* lgdInfo = new TLegend(0.3, 0.73, 0.6, 0.9);
        lgdInfo->SetMargin(0.05);
        lgdInfo->SetTextSize(0.04);
        lgdInfo->AddEntry((TObject*)0, Form("%.2f < y < %.2f", -rapMax, -rapMin), "");
        lgdInfo->AddEntry((TObject*)0, Form("%.1f < M < %.1f GeV/c^{2}", mMin, mMax), "");
        lgdInfo->AddEntry((TObject*)0, Form("%.1f < p_{T} < %.1f GeV/c", ptMin, ptMax), "");
        lgdInfo->Draw();
        
        TString figureName = Form("Plots/%s/%.1f-%.1f/Landau-%.1f-%.1f%s", period.c_str(), mMin, mMax, abs(rapMax), abs(rapMin), suffix.c_str());
        if (test) figureName+= "-test";
        cv->SaveAs(figureName+".pdf");
        if (eps) cv->SaveAs(figureName+".eps");
        
        TCanvas* cv2 = new TCanvas();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        //gPad->SetLogy();
        frame->GetXaxis()->SetRangeUser(0, 1);
        frame->Draw();
        lgd->Draw();
        lgdInfo->Draw();
        
        cv->SaveAs(figureName+"-zoom.pdf");
        if (eps) cv->SaveAs(figureName+"-zoom.eps");
    }
    ws->import(*bGg);
    ws->import(*nGg);
    ws->import(*yieldTwoGamma);
    ws->import(*yieldBkg);
    return;
}



void FitBackgroundPt(string rootfilePath = "", string rootfilePathMc = "", std::vector<string> periods = {"LHC16r", "LHC16s"}, Double_t mMin = 2.2, Double_t mMax = 2.7, Double_t ptMin = 0, Double_t ptMax = 3, Double_t rapMin = -4, Double_t rapMax = -2.5, bool exclusiveOnly = false) {
    
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetTitleSize(.05);
    gStyle->SetLabelSize(.05, "XY");
    //gStyle->SetMarkerSize(0.5);
    //gStyle->SetMarkerStyle(20);
    
    bool useCuts = true;
    
    const int nPeriod = periods.size();
    for (int k = 0; k<nPeriod; k++) {
        string period = periods[k];
        
        std::list<TCut> mCutList = DefineCuts(period, exclusiveOnly);
        // Define cuts
        TCut mCut = "";
        
        for (std::list <TCut>::iterator iter = mCutList.begin(); iter != mCutList.end(); ++iter) {mCut += *iter;}
        
        for (std::list <TCut>::iterator iter = mCutList.begin(); iter != mCutList.end(); ++iter) {
            cout << *iter << endl;
        }
        // Open the file
        TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()), "READ");
        
        // Connect to the tree
        TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
        //Create a new workspace to manage the project
        RooWorkspace* wspace = new RooWorkspace("GammaGamma");
        ImportDataSet(wspace, fAnaTree, mCut, mMin, mMax, ptMin, ptMax, rapMin, rapMax);
        wspace->Print();
        double chi2 = 0;
        MakePlots(wspace, rootfilePath, rootfilePathMc, period, exclusiveOnly, useCuts, mMin, mMax, ptMin, ptMax, rapMin, rapMax);
        delete wspace;
        
    }
    
    
}


