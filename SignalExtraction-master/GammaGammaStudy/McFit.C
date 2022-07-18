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
    //Int_t nBins = int((ptMax-ptMin)*100);
    Int_t nBins = 60;
    TH1F* histPt = new TH1F("hPt"+process, "hPt"+process, nBins, ptMin, ptMax);
    fAnaTree->Draw("fTrkTrkPt>>hPt"+process, cutMc);
    RooDataHist* dTemplatePt = new RooDataHist("ptHist"+process,"ptHist"+process, RooArgList(pt),histPt);
    RooHistPdf* ptPdf = new RooHistPdf("pt"+process, "pt"+process, pt, *dTemplatePt);
    ws->import(*ptPdf);
}


//____________________________________
void MakePlots(RooWorkspace *ws, string rootfilePath, string rootfilePathMc, string period, bool exclusiveOnly, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax){
    
    int ptBinNumber = int(10*(ptMax-ptMin));
    
    string suffix = "";
    if (exclusiveOnly) suffix = "-exclusive";
    
    //get what we need of the workspace
    RooRealVar* pt = ws->var("fTrkTrkPt");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    
    RooRealVar* yieldTwoGamma = new RooRealVar("yieldTwoGamma","yieldTwoGamma", 1000, 10., 20000);
    
    // With MC
    TCut cutMc = Form("fTrkTrkM > %f && fTrkTrkM < %f", mMin, mMax);
    GetPtHistMC(ws, rootfilePathMc, period, "kTwoGammaToMuLow", ptMin, ptMax, cutMc);
    RooAbsPdf* ptTwoGamma = ws->pdf("ptkTwoGammaToMuLow");
    ptTwoGamma->SetName("ptTwoGamma");
    
    RooAbsPdf* model = new RooAddPdf("fit", "fit", RooArgList(*ptTwoGamma), RooArgList(*yieldTwoGamma));
    //model->fitTo(*data, Extended(), Minos(true), Strategy(1));
    model->fitTo(*data);
    
    RooPlot* frame = pt->frame();
    int nBinsPt = 60;
    data->plotOn(frame, Binning(nBinsPt));
    //frame->GetXaxis()->SetRangeUser(0., 0.4);
    //data->plotOn(frame);
    
    model->plotOn(frame);
    model->plotOn(frame, Components(*model), Name("sum"), LineStyle(kDashed), LineColor(kRed));
    model->plotOn(frame, Components(*ptTwoGamma), Name("ptTwoGamma"), LineStyle(kDashed), LineColor(kBlue));
    
    double textSize = 0.04;
    double xPos = ptMin+(1-ptMin)*0.4;
    double yMax0 = frame->GetMaximum();
    
    
    // write parameters info
    double xLgd = 0.65;
    double yLgd = 0.33;
    double dx = 0.3;
    double dy = 0.18;
    TLegend* lgdParam = new TLegend(xLgd, yLgd, xLgd+dx, yLgd+dy);
    lgdParam->SetMargin(0.05);
    lgdParam->SetTextSize(0.037);
    lgdParam->AddEntry((TObject*)0,Form("%.1f #pm %.1f #gamma#gamma", yieldTwoGamma->getVal(), yieldTwoGamma->getError()), "");
    
    
    //make some plots
    TCanvas* cv = new TCanvas() ;
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    //gPad->SetLogy();
    
    // chi2
    int nFitParam = model->getParameters(data)->selectByAttrib("Constant",kFALSE)->getSize();
    double chi2 = frame->chiSquare("ptTwoGamma", "h_data", nFitParam);       // prints chi2/nBins
    int nDof = nBinsPt*0.4/1.2-nFitParam;
    double chi2Real = chi2*nBinsPt;
    
    cout << endl << endl << "chi2/ndf = " << chi2Real << "/" << nDof << " = " << chi2Real/nDof << endl;
    frame->Print();
    
    

    // Save histograms
    TH1* histBkg = data->createHistogram("hPtBackground",*pt, Binning(nBinsPt, 0, ptMax));
    TH1* histModel = model->createHistogram("hPtModel",*pt, Binning(nBinsPt, 0, ptMax));
    for (int i = 0; i < nBinsPt; i++) {
        histModel->SetBinError(i+1, 0);
    }
    histModel->SetOption("hist");
    histModel->Scale(107./histModel->GetMaximum());
    //histModel->Scale((yieldTwoGamma->getVal())/histModel->Integral());
    
    /*
    // normalize model histograms
    histModel->Scale(1./histModel->Integral());
    double factor = (double)nBinsModel/nBinsPt2;
    histModel->Scale((yieldBkg->getVal()+yieldTwoGamma->getVal())/histModel->Integral()*factor);
    histGgModel->Scale((double)yieldTwoGamma->getVal()/histGgModel->Integral()*factor);
    histExtraBkgModel->Scale((double)yieldBkg->getVal()/histExtraBkgModel->Integral()*factor);
     */
    
    /*
    // chi2 with other method
    RooChi2Var* chi2Var2 = new RooChi2Var("chi2","chi2", *ptTwoGamma, *data, DataError(RooAbsData::Poisson));
    Double_t chi2Var22 = chi2Var2->getVal();
    RooDataHist* roo_hist_bkg = new RooDataHist("roo_hist", "roo_hist", RooArgSet(*pt), histBkg);
    RooDataHist* roo_hist_bkg_model = new RooDataHist("roo_hist_model", "roo_hist_model", RooArgSet(*pt), histModel);
     
     */
    
     

     
     //cout << endl << endl << "yieldBkg->getVal()+yieldTwoGamma->getVal()" << yieldBkg->getVal()+yieldTwoGamma->getVal() << endl;
     // write them
     TFile* f = new TFile("mc-fit.root", "recreate");
     histBkg->Write();
     histModel->Write();
     f->Close();
    
    
    // Other draw options
    frame->SetTitle("");
    frame->Draw();
    
    TLegend* lgd;
    if (exclusiveOnly) lgd = new TLegend(0.68, 0.8, 0.9, 0.9);
    else lgd = new TLegend(0.67, 0.75, 0.9, 0.9);
    lgd->SetTextSize(0.04);
    lgd->AddEntry(frame->findObject("sum"), "Sum", "L");
    lgd->AddEntry(frame->findObject("ptTwoGamma"), "#gamma#gamma", "L");
    lgd->Draw();
    
    lgdParam->Draw();
    
    // write y and pt info
    TLegend* lgdInfo = new TLegend(0.3, 0.73, 0.6, 0.9);
    lgdInfo->SetMargin(0.05);
    lgdInfo->SetTextSize(0.04);
    lgdInfo->AddEntry((TObject*)0, Form("%.2f < y < %.2f", -rapMax, -rapMin), "");
    lgdInfo->AddEntry((TObject*)0, Form("%.1f < M < %.1f GeV/c^{2}", mMin, mMax), "");
    lgdInfo->AddEntry((TObject*)0, Form("%.1f < p_{T} < %.1f GeV/c", ptMin, ptMax), "");
    lgdInfo->Draw();
    
    
    cv->SaveAs(Form("Plots/%s/%.1f-%.1f/mc-fit-%.1f-%.1f%s.pdf", period.c_str(), mMin, mMax, abs(rapMax), abs(rapMin), suffix.c_str()));
}



void McFit() {
    
    Double_t mMin = 1.0, mMax = 2.5;
    Double_t ptMin = 0, ptMax = 1.2;
    Double_t rapMin = -4, rapMax = -2.5;
    bool exclusiveOnly = true;
    
    ptMax = 0.4;
    
    std::vector<string> periods = {"LHC16r"};
    string rootfilePath = "/Volumes/Transcend/rootFiles-pPb/std";
    string rootfilePathMc = "/Volumes/Transcend/rootFiles-pPb/MC-std";
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetTitleSize(.05);
    gStyle->SetLabelSize(.05, "XY");
    //gStyle->SetMarkerSize(0.5);
    //gStyle->SetMarkerStyle(20);
    
    
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
        MakePlots(wspace, rootfilePath, rootfilePathMc, period, exclusiveOnly, mMin, mMax, ptMin, ptMax, rapMin, rapMax);
        //cleanup
        delete wspace;
        
    }
    
    
}


