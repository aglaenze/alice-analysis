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

void WriteResults(RooWorkspace* ws, string period, Double_t rapMin, Double_t rapMax, Double_t mMin, Double_t mMax) {
    
    RooRealVar* mean_gg = ws->var("mean_gg");
    RooRealVar* sigma_gg = ws->var("sigma_gg");
    
    RooRealVar* yieldTwoGamma = ws->var("yieldTwoGamma");
    RooRealVar* yieldBkg = ws->var("yieldBkg");
    
    // Write the values in a text file
    string file = Form("output/%s/%.1f-%.1f/output-gg-%.1f-%.1f.txt", period.c_str(), mMin, mMax, -rapMax, -rapMin);
    ofstream monFlux(file.c_str(), ios::app); // append to the file
    
    if (monFlux) {
        // First parameters
        monFlux << mean_gg->getVal() <<  " " << mean_gg->getError() << " " << sigma_gg->getVal() <<  " " << sigma_gg->getError() << endl;
        // Then yields
        monFlux << yieldTwoGamma->getVal() <<  " " << yieldTwoGamma->getError() <<  " " << yieldBkg->getVal() <<  " " << yieldBkg->getError() << endl;
    }
    else {
        cout << "\n\nERREUR: Impossible d'ouvrir le fichier.\n\n" << endl;
    }
}

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
    
    RooRealVar* yieldTwoGamma = new RooRealVar("yieldTwoGamma","yieldTwoGamma", 200, 0., 2.e3);
    RooRealVar* yieldBkg = new RooRealVar("yieldBkg","yieldBkg", 500, 0., 2.e3);
    
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
    RooLandau* ptTwoGamma = new RooLandau("ptTwoGamma","ptTwoGamma", *pt, *mean, *sigma);
    
    
    //RooRealVar *pt0 = new RooRealVar("pt0","pt0", 0.4, 0.2, 3);
    RooRealVar *bInc = new RooRealVar("bInc","bInc", 5, 0, 8);
    RooRealVar *nInc = new RooRealVar("nInc","nInc", 3, 1, 8);
    //RooGenericPdf* ptBkg = new RooGenericPdf("ptBkg","Inclusive bkg PDF","fTrkTrkPt/((1.+(fTrkTrkPt/pt0)**2)**nInc)", RooArgSet(*pt, *pt0, *nInc)) ;
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
    

    //cout << "va ici" << endl;
    //return;
    //if (plot) {
        RooPlot* frame = pt->frame();
        int nBinsPt = 100;
        if (exclusiveOnly) nBinsPt = 20;
        Double_t nBinsPt2 = int(ptMax/fitRangeMax*nBinsPt);
        data->plotOn(frame, Binning(nBinsPt2));
        //data->plotOn(frame);
        
        model->plotOn(frame);
        model->plotOn(frame, Components(*model), Name("sum"), LineStyle(kDashed), LineColor(kRed));
        model->plotOn(frame, Components(*ptTwoGamma), Name("ptTwoGamma"), LineStyle(kDashed), LineColor(kBlue));
        if (twoComp) model->plotOn(frame, Components(*ptBkg), Name("ptBkg"), LineStyle(kDashed), LineColor(kGreen));
        
        double yMax0 = frame->GetMaximum();
        //double xPos = ptMin+(ptMax-ptMin)*0.3;
        double xPos = ptMin+(1-ptMin)*0.4;
        //if (mMax < mLimitPsi2s) xPos = mMin+(mMax-mMin)/3;
        TLatex* txt0 = new TLatex(xPos,0.9*yMax0,Form("#bf{%.1f #pm %.1f #gamma#gamma}", yieldTwoGamma->getVal(), yieldTwoGamma->getError()));
        TLatex* txt1 = new TLatex(xPos,0.8*yMax0,Form("#bf{%.1f #pm %.1f bkg}", yieldBkg->getVal(), yieldBkg->getError()));
        frame->addObject(txt0);
        if (twoComp) {frame->addObject(txt1);}
        
        // Fit parameters
        TLatex* txt2 = new TLatex(xPos,0.7*yMax0,Form("#bf{#mu(Landau) = %.3f #pm %.3f}", mean->getVal(), mean->getError() ));
        TLatex* txt3 = new TLatex(xPos,0.6*yMax0,Form("#bf{#sigma(Landau) = %.3f #pm %.3f}", sigma->getVal(), sigma->getError()));
        frame->addObject(txt2);
        frame->addObject(txt3);
        //TLatex* txt4 = new TLatex(xPos,0.4*yMax0,Form("#bf{p_{T0} = %.2f #pm %.2f}", pt0->getVal(), pt0->getError() ));
        TLatex* txt4 = new TLatex(xPos,0.4*yMax0,Form("#bf{b_{inc} = %.2f #pm %.2f}", bInc->getVal(), bInc->getError() ));
        TLatex* txt5 = new TLatex(xPos,0.3*yMax0,Form("#bf{n_{inc} = %.2f #pm %.2f}", nInc->getVal(), nInc->getError()));
        if (twoComp) {frame->addObject(txt4); frame->addObject(txt5);}
        
    //cout << "on est là " << endl;
    //return;
        //make some plots
        TCanvas* cv = new TCanvas() ;
        gPad->SetLeftMargin(0.15) ;
        gPad->SetBottomMargin(0.15) ;
        
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
        
    //cout << "check encore " << endl;
    //return;
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
        cv->SaveAs(Form("Plots/%s/%.1f-%.1f/Landau-%.1f-%.1f%s.pdf", period.c_str(), mMin, mMax, abs(rapMax), abs(rapMin), suffix.c_str()));
        
        TCanvas* cv2 = new TCanvas();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        frame->GetXaxis()->SetRangeUser(0, 1);
        frame->Draw();
        lgd->Draw();
        cv->SaveAs(Form("Plots/%s/%.1f-%.1f/Landau-%.1f-%.1f%s-zoom.pdf", period.c_str(), mMin, mMax, abs(rapMax), abs(rapMin), suffix.c_str()));
    //}
    ws->import(*mean);
    ws->import(*sigma);
    ws->import(*yieldTwoGamma);
    ws->import(*yieldBkg);
    return;
}



void FitBackgroundPt2(string rootfilePath = "", string rootfilePathMc = "", std::vector<string> periods = {"LHC16r", "LHC16s"}, Double_t mMin = 2.2, Double_t mMax = 2.7, Double_t ptMin = 0, Double_t ptMax = 3, Double_t rapMin = -4, Double_t rapMax = -2.5, bool exclusiveOnly = false) {
    
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetTitleSize(.05);
    gStyle->SetLabelSize(.05, "XY");
    //gStyle->SetMarkerSize(0.5);
    //gStyle->SetMarkerStyle(20);
    
    bool useCuts = true;
    if (exclusiveOnly) ptMax = 1.2;
    
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
        if (!exclusiveOnly) WriteResults(wspace, period, rapMin, rapMax, mMin, mMax);
        //cleanup
        delete wspace;
        
    }
    
    
}


