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

//____________________________________
void MakePlots(RooWorkspace *ws, string rootfilePath, string rootfilePathMc, string period, bool exclusiveOnly, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax){
    
    // get MC data
    TCut cutMc = Form("fTrkTrkM > %f && fTrkTrkM < %f && fTrkTrkY > %f && fTrkTrkY < %f ", mMin, mMax, rapMin, rapMax);
    //TCut cutMc = Form("fTrkTrkM > %.1f && fTrkTrkM < %.1f ", mMin, mMax);
    TFile *fSimu = new TFile(Form("%s/AnalysisResults_%s_MC_kTwoGammaToMuLow.root", rootfilePathMc.c_str(), period.c_str()),"READ");
    TTree* fAnaTree = (TTree*)fSimu->Get("MyTask/fAnaTree");
    
    RooRealVar* pt = new RooRealVar("fTrkTrkPt","p_{T} (GeV/c)", ptMin, ptMax);
    RooRealVar* m = new RooRealVar("fTrkTrkM","M(GeV/c^{2})", mMin, mMax);
    RooRealVar* rap = new RooRealVar("fTrkTrkY","Y", rapMin, rapMax);
    RooDataSet* data = new RooDataSet("data","data",RooArgSet(*pt, *m, *rap),Import(*fAnaTree), Cut(cutMc));
    TH1* h = new TH1F("h", "h", 1, 0, 0.38);
    TH1* h2 = new TH1F("h2", "h2", 1, 0, 1.2);
    fAnaTree->Draw("fTrkTrkPt>>h", cutMc);
    fAnaTree->Draw("fTrkTrkPt>>h2", cutMc);
    
    std::cout << std::endl << std::endl << "Entries before 0.38 GeV/c = " << h->Integral() << std::endl;
    std::cout << "Entries before 1.2 GeV/c = " << h2->Integral() << std::endl << std::endl << std::endl;
    
    int nBinsPt = 1200;
    double fitRangeMax = ptMax;
    //double fitRangeMax = 0.38;
    
    
    RooRealVar* yieldTwoGamma = new RooRealVar("yieldTwoGamma","yieldTwoGamma", 300, 0., 2.e5);
    RooRealVar* yieldBkg = new RooRealVar("yieldBkg","yieldBkg", 100, 0., 2.e3);
    
    
    // With Landau
    RooRealVar* mean = new RooRealVar("mean","mean",0.06, 0.01,0.1);
    RooRealVar* sigma = new RooRealVar("sigma","sigma", 0.03, 0.01, 0.1); //3.43,3.73);
    RooAbsPdf* ptTwoGamma = new RooLandau("ptTwoGamma","ptTwoGamma", *pt, *mean, *sigma);
    
    mean->setVal(0.056);
    mean->setConstant();
    sigma->setVal(0.023);
    sigma->setConstant();
    
    RooRealVar *bInc = new RooRealVar("bInc","bInc", 4, 1, 6);
    RooRealVar *nInc = new RooRealVar("nInc","nInc", 2, 1, 6);
    RooGenericPdf* ptBkg = new RooGenericPdf("ptBackground","ptBackground","(2*fTrkTrkPt*(1.+(fTrkTrkPt**2)*(bInc/nInc))**(-nInc))", RooArgSet(*pt, *bInc, *nInc));
    
    
    //RooAbsPdf* model = new RooAddPdf("fit", "fit", RooArgList(*ptTwoGamma), RooArgList(*yieldTwoGamma), kFALSE);
    RooAbsPdf* model = new RooAddPdf("fit", "fit", RooArgList(*ptTwoGamma, *ptBkg), RooArgList(*yieldTwoGamma, *yieldBkg), kFALSE);
    model->fitTo(*data, Range(0, fitRangeMax), Extended(), Minos(true), Strategy(1));
    
    
    RooPlot* frame = pt->frame();
    data->plotOn(frame, Binning(nBinsPt));
    //data->plotOn(frame);
    
    model->plotOn(frame);
    model->plotOn(frame, Components(*model), Name("sum"), LineStyle(kDashed), LineColor(kRed));
    model->plotOn(frame, Components(*ptTwoGamma), Name("ptTwoGamma"), LineStyle(kDashed), LineColor(kBlue));
    model->plotOn(frame, Components(*ptBkg), Name("ptBkg"), LineStyle(kDashed), LineColor(kGreen));
    
    double textSize = 0.04;
    double xPos = ptMin+(1-ptMin)*0.4;
    double yMax0 = frame->GetMaximum();
    /*
     //double xPos = ptMin+(ptMax-ptMin)*0.3;
     //if (mMax < mLimitPsi2s) xPos = mMin+(mMax-mMin)/3;
     TLatex* txt0 = new TLatex(xPos,0.7*yMax0,Form("#bf{%.1f #pm %.1f #gamma#gamma}", yieldTwoGamma->getVal(), yieldTwoGamma->getError()));
     txt0->SetTextSize(textSize);
     
     // Fit parameters
     TLatex* txt2 = new TLatex(xPos,0.62*yMax0,Form("#bf{#mu(Landau) = %.3f #pm %.3f}", mean->getVal(), mean->getError() ));
     txt2->SetTextSize(textSize);
     TLatex* txt3 = new TLatex(xPos,0.54*yMax0,Form("#bf{#sigma(Landau) = %.3f #pm %.3f}", sigma->getVal(), sigma->getError()));
     txt3->SetTextSize(textSize);
     frame->addObject(txt2);
     frame->addObject(txt3);
     
     TLatex* txt1 = new TLatex(xPos,0.46*yMax0,Form("#bf{%.1f #pm %.1f bkg}", yieldBkg->getVal(), yieldBkg->getError()));
     txt1->SetTextSize(textSize);
     frame->addObject(txt0);
     if (twoComp) {frame->addObject(txt1);}
     //TLatex* txt4 = new TLatex(xPos,0.4*yMax0,Form("#bf{p_{T0} = %.2f #pm %.2f}", pt0->getVal(), pt0->getError() ));
     TLatex* txt4 = new TLatex(xPos,0.38*yMax0,Form("#bf{b_{inc} = %.1f #pm %.1f}", bInc->getVal(), bInc->getError() ));
     txt4->SetTextSize(textSize);
     TLatex* txt5 = new TLatex(xPos,0.3*yMax0,Form("#bf{n_{inc} = %.1f #pm %.1f}", nInc->getVal(), nInc->getError()));
     txt5->SetTextSize(textSize);
     if (twoComp) {frame->addObject(txt4); frame->addObject(txt5);}
     */
    // write parameters info
    double xLgd = 0.65;
    double yLgd = 0.33;
    double dx = 0.3;
    double dy = 0.18;
    TLegend* lgdParam = new TLegend(xLgd, yLgd, xLgd+dx, yLgd+dy);
    lgdParam->SetMargin(0.05);
    lgdParam->SetTextSize(0.037);
    lgdParam->AddEntry((TObject*)0,Form("%.1f #pm %.1f #gamma#gamma", yieldTwoGamma->getVal(), yieldTwoGamma->getError()), "");
    lgdParam->AddEntry((TObject*)0, Form("#mu(Landau) = %.3f #pm %.3f", mean->getVal(), mean->getError()), "");
    lgdParam->AddEntry((TObject*)0, Form("#sigma(Landau) = %.3f #pm %.3f", sigma->getVal(), sigma->getError()), "");
    
    
    //make some plots
    TCanvas* cv = new TCanvas() ;
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    //gPad->SetLogy();
    
    /*
     // Chi2
     // Chi2 will be calculated until pT = 1 so that there would be no empty bin
     Double_t ptMaxchi2 = 1;
     Int_t nBinsChi2 = int(ptMaxchi2/fitRangeMax*nBinsPt);
     //TH1* histBkg = data->createHistogram("hPtBackground",*pt, Binning(nBinsChi2, 0, ptMaxchi2));
     TH1* histBkg = data->createHistogram("hPtBackground",*pt, Binning(nBinsPt2, 0, ptMax));
     RooDataHist* roo_hist_bkg = new RooDataHist("roo_hist", "roo_hist", RooArgSet(*pt), histBkg);
     //TH1* histBkgModel = model->createHistogram("hPtBackgroundModel",*pt, Binning(nBinsChi2, 0, ptMaxchi2));
     int nBinsModel = 1000;
     TH1* histModel = model->createHistogram("hPtModel",*pt, Binning(nBinsModel, 0, ptMax));
     TH1* histGgModel = ptTwoGamma->createHistogram("hPtGgModel",*pt, Binning(nBinsModel, 0, ptMax));
     TH1* histExtraBkgModel = ptBkg->createHistogram("hPtExtraModel",*pt, Binning(nBinsModel, 0, ptMax));
     
     RooDataHist* roo_hist_bkg_model = new RooDataHist("roo_hist_model", "roo_hist_model", RooArgSet(*pt), histModel);
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
     frame->addObject(txtChi) ;
     frame->chiSquare() ;
     */
    
    /*
    // normalize model histograms
    histModel->Scale(1./histModel->Integral());
    double factor = (double)nBinsModel/nBinsPt2;
    histModel->Scale((yieldBkg->getVal()+yieldTwoGamma->getVal())/histModel->Integral()*factor);
    histGgModel->Scale((double)yieldTwoGamma->getVal()/histGgModel->Integral()*factor);
    histExtraBkgModel->Scale((double)yieldBkg->getVal()/histExtraBkgModel->Integral()*factor);
    
    histModel->SetOption("hist");
    histGgModel->SetOption("hist");
    histExtraBkgModel->SetOption("hist");
    
    for (int i = 0; i < nBinsModel; i++) {
        histModel->SetBinError(i+1, 0);
        histGgModel->SetBinError(i+1, 0);
        histExtraBkgModel->SetBinError(i+1, 0);
    }
    
    //cout << endl << endl << "yieldBkg->getVal()+yieldTwoGamma->getVal()" << yieldBkg->getVal()+yieldTwoGamma->getVal() << endl;
    // write them
    TFile* f = new TFile("mc-landau.root", "recreate");
    histBkg->Write();
    histModel->Write();
    histExtraBkgModel->Write();
    f->Close();
     */
    
    
    
    // Save histograms
    TH1* histBkg = data->createHistogram("hPtBackground",*pt, Binning(nBinsPt, 0, ptMax));
    TH1* histModel = model->createHistogram("hPtModel",*pt, Binning(nBinsPt, 0, ptMax));
    for (int i = 0; i < nBinsPt; i++) {
        histModel->SetBinError(i+1, 0);
    }
    histModel->SetOption("hist");
    histModel->Scale((yieldTwoGamma->getVal())/histModel->Integral());

    // write them
    TFile* f = new TFile("mc-landau-fit.root", "recreate");
    histBkg->Write();
    histModel->Write();
    f->Close();
    
    
    // Other draw options
    frame->SetTitle("");
    frame->Draw();
    
    TLegend* lgd;
    lgd = new TLegend(0.68, 0.8, 0.9, 0.9);
    lgd->SetTextSize(0.04);
    lgd->AddEntry(frame->findObject("sum"), "sum", "L");
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
    
    gPad->SetLogy();
    
    cv->SaveAs(Form("Plots/%s/%.1f-%.1f/Landau-mc-%.1f-%.1f-test.pdf", period.c_str(), mMin, mMax, abs(rapMax), abs(rapMin)));
    
    return;
}



void LandauOnMc() {
    
    Double_t mMin = 1.0, mMax = 2.5;
    //Double_t ptMin = 0, ptMax = 1.2;
    Double_t ptMin = 0, ptMax = 3.0;
    Double_t rapMin = -4, rapMax = -2.5;
    
    string rootfilePath = "/Volumes/Transcend/rootFiles-pPb/std";
    string rootfilePathMc = "/Volumes/Transcend/rootFiles-pPb/MC-std";
    std::vector<string> periods = {"LHC16r"};
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetTitleSize(.05);
    gStyle->SetLabelSize(.05, "XY");
    //gStyle->SetMarkerSize(0.5);
    //gStyle->SetMarkerStyle(20);
    
    bool exclusiveOnly = true;
    
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


