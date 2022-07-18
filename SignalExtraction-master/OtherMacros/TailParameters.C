#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>
#include <TROOT.h>

#include "TSystem.h"
#include "TString.h"
#include "TDatime.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TColor.h"
#include <TROOT.h>
#include <TMath.h>
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
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"

#include "../Include/ExtendedCrystalBall.h"

using namespace RooFit;
using namespace std;

Double_t mMin = 2.5;
Double_t mMax = 3.5;

RooPlot* DrawCB(string rootfilePath, string period, TString process, bool log, bool drawText) {
    
    
    TFile *fSimu = new TFile(Form("%s/AnalysisResults_%s_MC_", rootfilePath.c_str(), period.c_str()) + process + ".root","READ");
    TTree* fAnaTree = (TTree*)fSimu->Get("MyTask/fAnaTree");
    
    Double_t mMin = 2.5, mMax = 3.5;
    Double_t ptMin = 0, ptMax = 1.2;
    Double_t rapMin = -4.0, rapMax = -2.5;

    
    TCut cutMc = Form("fTrkTrkM > %f && fTrkTrkM < %f && fTrkTrkY > %f && fTrkTrkY < %f && fTrkTrkPt > %f && fTrkTrkPt < %f", mMin, mMax, rapMin, rapMax, ptMin, ptMax);
    //TCut cutMc = Form("fTrkTrkM > %.1f && fTrkTrkM < %.1f ", mMin, mMax);

    
    RooRealVar* pt = new RooRealVar("fTrkTrkPt","p_{T} (GeV/c)", ptMin, ptMax);
    RooRealVar* m = new RooRealVar("fTrkTrkM","M(GeV/c^{2})", mMin, mMax);
    RooRealVar* rap = new RooRealVar("fTrkTrkY","Y", rapMin, rapMax);
    RooDataSet* data = new RooDataSet("data","data",RooArgSet(*pt, *m, *rap),Import(*fAnaTree), Cut(cutMc));

    
    // Import binned data (otherwise it never ends)
    int nBins = 500;
    TH1* histData = data->createHistogram("hMass",*m, Binning(nBins, mMin, mMax));
    RooDataHist* dTemplateM = new RooDataHist("mHist"+process,"mHist"+process, RooArgList(*m),histData);
    
    
    double x = 3.1, xmin = 2.9, xmax = 3.3;
    if (process == "kIncohPsi2sToMu") {x = 3.6, xmin = 3.5, xmax = 3.7;}
    
    RooRealVar mean_jpsi("mean_jpsi","mean_jpsi", x, xmin, xmax);
    RooRealVar sigma_jpsi("sigma_jpsi","sigma_jpsi", 0.05,0,0.5);
    RooRealVar alpha_jpsi_L("alpha_jpsi_L","alpha_jpsi_L", 1, 0, 10);
    RooRealVar n_jpsi_L("n_jpsi_L","n_jpsi_L", 2, 0, 10);
    RooRealVar alpha_jpsi_R("alpha_jpsi_R","alpha_jpsi_R", 3, 0, 10);
    RooRealVar n_jpsi_R("n_jpsi_R","n_jpsi_R", 3, 0, 10);
    
    ExtendedCrystalBall *jpsi = new ExtendedCrystalBall("jpsi","crystal ball PDF", *m,
                                                        mean_jpsi, sigma_jpsi, alpha_jpsi_L,
                                                        n_jpsi_L, alpha_jpsi_R, n_jpsi_R);
    RooRealVar yield1("yield"+process,"yield"+process, 1.e6, 10, 1.e7);
    RooArgList* pdfList = new RooArgList(*jpsi);
    RooArgList yieldList = RooArgList(yield1);
    // Create fit model
    RooAbsPdf* model = new RooAddPdf("model", "model", *pdfList, yieldList, kFALSE);
    //model->fitTo(*data,Minos(true),Strategy(2));
    RooFitResult* r = model->fitTo(*dTemplateM, Minos(true), Strategy(2), Save());
    
    TString txtt1 = Form("#bf{#alpha_{L} = %.3f #pm %.3f }", alpha_jpsi_L.getVal(), alpha_jpsi_L.getError());
    TString txtt2 = Form("#bf{n_{L} = %.3f #pm %.3f }", n_jpsi_L.getVal(), n_jpsi_L.getError());
    TString txtt3 = Form("#bf{#alpha_{R} = %.3f #pm %.3f}", alpha_jpsi_R.getVal(), alpha_jpsi_R.getError());
    TString txtt4 = Form("#bf{n_{R} = %.3f #pm %.3f}", n_jpsi_R.getVal(), n_jpsi_R.getError());
    
    TString txtt5 = Form("#bf{#mu = %.3f #pm %.3f}", mean_jpsi.getVal(), mean_jpsi.getError());
    TString txtt6 = Form("#bf{#sigma = %.3f #pm %.3f}", sigma_jpsi.getVal(), sigma_jpsi.getError());

    cout << txtt1 << endl;
    cout << txtt2 << endl;
    cout << txtt3 << endl;
    cout << txtt4 << endl;
    cout << txtt5 << endl;
    cout << txtt6 << endl;
    
    TString titleName = "";
    if (process == "kIncohJpsiToMu") titleName = "J/#Psi";
    else if (process == "kIncohPsi2sToMu") titleName = "#Psi(2s)";
    RooPlot* mframe = m->frame(Title("Fit of invariant mass of the "+titleName));
    dTemplateM->plotOn(mframe);
    model->plotOn(mframe, Name("sum"), LineColor(kRed));
    
    double yMax = mframe->GetMaximum();
    double y1 = 0.9*yMax, y2 = 0.8*yMax, y3 = 0.7*yMax, y4 = 0.6*yMax, y5 = 0.5*yMax, y6 = 0.4*yMax;
    if (log) {y1 = yMax/pow(2,1); y2 = yMax/pow(2,2); y3 = yMax/pow(2,3); y4 = yMax/pow(2,4); y5 = yMax/pow(2,5); y6 = yMax/pow(2,6);}
    double xPos = mMin + (mMax-mMin)*0.05;
    cout << endl << endl << "xPos = " << xPos << endl << endl << endl;
    if (process == "kIncohPsi2sToMu") xPos = mMin + (mMax-mMin)*0.1;
    TLatex* txt1 = new TLatex(xPos,y1, txtt1);
    TLatex* txt2 = new TLatex(xPos,y2, txtt2);
    TLatex* txt3 = new TLatex(xPos,y3, txtt3);
    TLatex* txt4 = new TLatex(xPos,y4, txtt4);
    //TLatex* txt5 = new TLatex(xPos,y5,Form("#bf{# candidates = %.1f #pm %.1f}", yield1.getVal(), yield1.getError()));
    TLatex* txt5 = new TLatex(xPos,y5, txtt5);
    TLatex* txt6 = new TLatex(xPos,y6, txtt6);
    if (drawText) {
        mframe->addObject(txt1);
        mframe->addObject(txt2);
        mframe->addObject(txt3);
        mframe->addObject(txt4);
        //mframe->addObject(txt5) ;
        mframe->addObject(txt5);
        mframe->addObject(txt6);
    }
    
    // Save histograms
    TH1* histModel = model->createHistogram("hModel",*m, Binning(nBins, mMin, mMax));
    for (int i = 0; i < nBins; i++) {
        histModel->SetBinError(i+1, 0);
    }
    histModel->SetOption("hist");
    histModel->Scale((yield1.getVal())/histModel->Integral());
    
    // write them
    TFile* f = new TFile("cb2.root", "recreate");
    histData->Write();
    histModel->Write();
    f->Close();
    
    
    return mframe;
}

void GetNCandidatesAtLowMass(string rootfilePath, string period, TString process) {
    
    
    TFile *fSimu = new TFile(Form("%s/AnalysisResults_%s_MC_", rootfilePath.c_str(), period.c_str()) + process + ".root","READ");
    TTree* fAnaTree = (TTree*)fSimu->Get("MyTask/fAnaTree");
    TTree* fGenTree = (TTree*)fSimu->Get("MyTask/fGenTree");
    
    double mLimit = 2.5;
    
    TH1F* h1 = new TH1F("h1", "h1", 1000, 0, mLimit);
    fGenTree->Draw("fMCTrkTrkM>>h1", Form("fMCTrkTrkM<%.1f", mLimit));
    cout << "Generated: number of JPsi below " << mLimit << "GeV = " << h1->GetEntries() << endl;
    TH1F* h2 = new TH1F("h2", "h2", 1000, 0, 10);
    fGenTree->Draw("fMCTrkTrkM>>h2");
    cout << "Generated: total number of JPsi = " << h2->GetEntries() << endl;
    
    TH1F* h3 = new TH1F("h3", "h3", 1000, 0, mLimit);
    fAnaTree->Draw("fTrkTrkM>>h3", Form("fTrkTrkM<%.1f", mLimit));
    cout << "Reconstructed: number of JPsi below " << mLimit << "GeV = " << h3->GetEntries() << endl;
    TH1F* h4 = new TH1F("h4", "h4", 1000, 0, 10);
    fAnaTree->Draw("fTrkTrkM>>h4");
    cout << "Reconstructed: total number of JPsi = " << h4->GetEntries() << endl;
    
}


void ExtractB(string rootfilePath, string period) {
    
    gStyle->SetOptStat(0);
    
    TFile *fSimu = new TFile(Form("%s/AnalysisResults_%s_MC_kIncohJpsiToMu.root", rootfilePath.c_str(), period.c_str()),"READ");
    TTree* fSimuTree = (TTree*)fSimu->Get("MyTask/fAnaTree");
    
    RooRealVar pt("fTrkTrkPt","p_{T} (GeV/c2)", 0, 4);
    RooArgSet variables(pt);
    RooDataSet* data = new RooDataSet("data","data",variables,Import(*fSimuTree));
    
    RooRealVar bExc("bExc","bExc", 6, 2, 10);
    RooAbsPdf* jpsiExclusive = new RooGenericPdf("jpsiExc","exclusive jPsi PDF","(fTrkTrkPt*exp(-bExc*(fTrkTrkPt**2)))",RooArgSet(pt,bExc)) ;
    jpsiExclusive->SetName("jpsiExc");
    RooRealVar fsigJpsiExc("fsigJpsiExc","fsigJpsiExc",1000,0.,1.e8);
    
    RooAbsPdf* model = new RooAddPdf("ptfit", "ptfit", RooArgList(*jpsiExclusive), RooArgList(fsigJpsiExc), kFALSE);
    
    model->fitTo(*data, Extended(), Minos(true), Strategy(2));
    
    RooPlot* ptframe = pt.frame(Title("Fit of pt of exclusive J/Psi"));
    model->plotOn(ptframe, Name("MC"), LineColor(kRed));
    
    double yMax = ptframe->GetMaximum();
    TLatex* txtExc = new TLatex(1.5,0.45*yMax,Form("b_{exc} = %.2f #pm %.2f", bExc.getVal(), bExc.getError()));
    ptframe->addObject(txtExc) ;
    
    ptframe->Draw();
}

void TailParameters() {
    
    string rootfilePath="/Volumes/Transcend/rootFiles-pPb/MC-std/";
    std::vector<string> periods = {"LHC16r"};
    bool logScale = false;
    
    gStyle->SetOptStat(0);
    gStyle->SetTextSize(0.045);
    
    const int nPeriod = periods.size();
    
    if (false) {
        TCanvas* cv = new TCanvas("cv","cv",800,800) ;
        cv->Divide(1,2);
        for (int k = 0; k<nPeriod;k++) {
            string period = periods[k];
            cv->cd(k+1);
            ExtractB(rootfilePath, period);
        }
        cv->SaveAs("Plots/ExcBparam.pdf");
    }
    
    
    gROOT->ProcessLine(".L ../Include/ExtendedCrystalBall.cxx+") ;
    gSystem->Load("./../Include/ExtendedCrystalBall_cxx.so") ;
    
    for (int k = 0; k<nPeriod;k++) {
        string period = periods[k];
        
        RooPlot* mframeJpsi = DrawCB(rootfilePath, period, "kIncohJpsiToMu", logScale, true);
        //RooPlot* mframePsi2s = DrawCB(rootfilePath, period, "kIncohPsi2sToMu", logScale);
        RooPlot* mframeJpsi2 = DrawCB(rootfilePath, period, "kIncohJpsiToMu", logScale, false);
        
        TCanvas* cv2 = new TCanvas("cv","cv",800,400) ;
        cv2->Divide(2,1);
        
        cv2->cd(1);
        if (logScale) gPad->SetLogy();
        gPad->SetLeftMargin(0.15) ;
        mframeJpsi->Draw();
        cv2->cd(2);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15) ;
        mframeJpsi2->Draw();
        //mframePsi2s->Draw();
        
        cv2->SaveAs(Form("Plots/2sidedCB-%s.pdf", period.c_str()));
        GetNCandidatesAtLowMass(rootfilePath, period, "kIncohJpsiToMu");
    }
}

