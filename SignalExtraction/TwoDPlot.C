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
#include "RooMinimizer.h"
#include "RooChi2Var.h"

#include "Include/FitUtils.C"

using namespace RooFit;
using namespace std;

Double_t mLimitPsi2s = 3.65;

void WriteResults(RooWorkspace* ws, string period, Double_t rapMin, Double_t rapMax, bool exp, double chi2, bool noInclusive) {
    
    RooRealVar* a1 = ws->var("a1");
    RooRealVar* bExc = ws->var("bExc");
    RooRealVar* bDiss = ws->var("bDiss");
    RooRealVar* nDiss = nullptr;
    if (!exp) nDiss = ws->var("nDiss");
    else {nDiss = new RooRealVar("nDiss", "nDiss", 0); nDiss->setConstant();}
    RooRealVar* nInc = nullptr;
    RooRealVar* pt0 = nullptr;
    if (noInclusive) {
        nInc = new RooRealVar("nInc", "nInc", 0); nInc->setConstant();
        pt0 = new RooRealVar("pt0", "pt0", 0); pt0->setConstant();
    }
    else {
        nInc = ws->var("nInc");
        pt0 = ws->var("pt0");
    }
    RooRealVar* mean_gg = ws->var("mean_gg");
    RooRealVar* sigma_gg = ws->var("sigma_gg");
    
    RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
    RooRealVar* yieldJpsiDissociative = ws->var("yieldJpsiDissociative");
    RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
    RooRealVar* yieldJpsiInclusive = ws->var("yieldJpsiInclusive");
    RooRealVar* yieldPsi2s = ws->var("yieldPsi2s");
    RooRealVar* yieldTwoGamma = ws->var("yieldTwoGamma");
    RooRealVar* yieldBkg = ws->var("yieldBkg");
    RooRealVar* r_diss_exc = ws->var("r_diss_exc");
    
    // Write the values in a text file
    string file = "output-"+period;
    if (exp) file += "-exp";
    else file += "-powerlaw";
    file += Form("-%.1f-%.1f", -rapMax, -rapMin);
    file += ".txt";
    ofstream monFlux(file.c_str(), ios::app); // append to the file
    
    if (monFlux) {
        // First parameters
        monFlux << a1->getVal() <<  " " << a1->getError() <<  " " << bExc->getVal() <<  " " << bExc->getError() <<  " " << bDiss->getVal() <<  " " << bDiss->getError() <<  " " << nDiss->getVal() <<  " " << nDiss->getError() <<  " " << nInc->getVal() <<  " " << nInc->getError() <<  " " << pt0->getVal() <<  " " << pt0->getError() << " " << mean_gg->getVal() <<  " " << mean_gg->getError() << " " << sigma_gg->getVal() <<  " " << sigma_gg->getError() << " " << chi2 << endl;
        // Then yields
        monFlux << yieldJpsiExclusive->getVal() <<  " " << yieldJpsiExclusive->getError() <<  " " << yieldJpsiDissociative->getVal() <<  " " << yieldJpsiDissociative->getError() <<  " " << yieldJpsiGammaPb->getVal() <<  " " << yieldJpsiGammaPb->getError() <<  " " << yieldJpsiInclusive->getVal() <<  " " << yieldJpsiInclusive->getError() <<  " " << yieldTwoGamma->getVal() <<  " " << yieldTwoGamma->getError() <<  " " << yieldBkg->getVal() <<  " " << yieldBkg->getError() <<  " " << r_diss_exc->getVal() <<  " " << r_diss_exc->getError() << endl;
    }
    else {
        cout << "\n\nERREUR: Impossible d'ouvrir le fichier.\n\n" << endl;
    }
    
    cout << endl << endl << "ratio = " << r_diss_exc->getVal() << " pm " << r_diss_exc->getError() << endl;
}


void AddModel(RooWorkspace* ws, string rootfilePath, string rootfilePathMC, string period, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax, bool exp, bool exclusiveOnly, bool noInclusive = true) {
    // Define 2D model
    // First define fits in mass
    // Second define fits in pt
    // Then make the product to get 2D PDFs
    
    LoadMassFitFunctions(ws, period);
    LoadJpsiPtFitFunctions(ws, rootfilePathMC, period, ptMin, ptMax, exp, exclusiveOnly);
    LoadBkgPtFitFunctions(ws, rootfilePath, rootfilePathMC, period, exclusiveOnly);
    
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
    //RooArgList* pdfList = new RooArgList(*jpsi, *bkg);
    RooArgList* yieldList = new RooArgList(*yieldJpsiExclusive, *yieldJpsiGammaPb, *yieldTwoGamma, *yieldBkg);
    
    if (mMax > mLimitPsi2s) {
        pdfList->add(*pdfPsi2s);
        yieldList->add(*yieldPsi2s);
    }
    if (!exclusiveOnly) {
        pdfList->add(*pdfJpsiDissociative);
        yieldList->add(*yieldJpsiDissociative);
    }
    if (!noInclusive) {
        pdfList->add(*pdfJpsiInclusive);
        yieldList->add(*yieldJpsiInclusive);
    }
    
    // Create fit model
    RooAbsPdf* fitModel = new RooAddPdf("model", "model", *pdfList, *yieldList, kFALSE);
    
    ws->import(*fitModel);
}

std::vector <int> FindBinNumbers(RooWorkspace* ws) {
    
    // get what we need
    RooRealVar* m = ws->var("fTrkTrkM");
    RooRealVar* pt = ws->var("fTrkTrkPt");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    RooAbsPdf* fitModel = ws->pdf("model");
    
    /*
     RooRealVar mCopy = *m;
     RooRealVar ptCopy = *pt;
     mCopy.setRange(2.7, 3.3),
     ptCopy.setRange(0,2);
     */
    //m->setRange(2.7, 3.3),
    //pt->setRange(0,2);
    
    int nFitParam = fitModel->getParameters(data)->selectByAttrib("Constant",kFALSE)->getSize();
    
    int mBinNumber = 1000;
    int ptBinNumber = 0;
    int ptBinNumberMax = 1000;
    
    int minBinM = 2;
    int minBinPt = 4;
    // Find maximum value of ptBinNumber
    TH1* hist = data->createHistogram("fTrkTrkM,fTrkTrkPt", minBinM, ptBinNumberMax);
    while (hist->GetMinimum() < 1) {
        ptBinNumberMax--;
        delete hist;
        hist = data->createHistogram("fTrkTrkM,fTrkTrkPt", minBinM, ptBinNumberMax);
    }
    
    cout << "max ptBinNumber = " << ptBinNumberMax << " to have at least " << minBinM << " bins in M" << endl;
    
    
    std::vector <int> ptBinVector = {};
    std::vector <int> mBinVector = {};
    std::vector <int> productVector = {};
    std::vector <double> chi2Vec = {};
    int ndf = 0;
    for (int k = minBinPt; k < ptBinNumberMax+1; k++) {
        ptBinNumber = k;
        hist = data->createHistogram("fTrkTrkM,fTrkTrkPt", mBinNumber, ptBinNumber);
        while (hist->GetMinimum() < 1 && mBinNumber >= 2) {
            mBinNumber--;
            delete hist;
            hist = data->createHistogram("fTrkTrkM,fTrkTrkPt", mBinNumber, ptBinNumber);
        }
        
        // Compute chi2
        ndf = mBinNumber*ptBinNumber-nFitParam;
        delete hist;
        if (ndf < 1) continue;
        if (mBinNumber < minBinM) continue;
        //hist = data->createHistogram("fTrkTrkM,fTrkTrkPt", mBinNumber+1, ptBinNumber+1);
        hist = data->createHistogram("fTrkTrkM,fTrkTrkPt", mBinNumber, ptBinNumber);
        RooDataHist* roo_hist = new RooDataHist("roo_hist", "roo_hist", RooArgSet(*m,*pt), hist);
        
        // Construct histo PDF
        TH1* hist2 = fitModel->createHistogram("fTrkTrkM,fTrkTrkPt", mBinNumber, ptBinNumber);
        RooDataHist* hist_pdf = new RooDataHist("roo_hist2", "roo_hist2", RooArgSet(*m,*pt), hist2);
        RooHistPdf* roo_hist_pdf = new RooHistPdf("roo_hist_pdf", "roo_hist_pdf", RooArgSet(*m,*pt), *hist_pdf);
        
        /*
         TH1* h1 = roo_hist->createHistogram("fTrkTrkM,fTrkTrkPt", mBinNumber, ptBinNumber);
         TH1* h2 = roo_hist_pdf->createHistogram("fTrkTrkM,fTrkTrkPt", mBinNumber, ptBinNumber);
         
         h1->Scale(1);
         h2->Scale(h1->Integral()/h2->Integral());
         
         if (k==5) {
         TCanvas* cv = new TCanvas("", "", 600, 300);
         cv->Divide(2);
         cv->cd(1);
         gPad->SetLeftMargin(0.15);
         gPad->SetBottomMargin(0.15);
         gPad->SetRightMargin(0.15);
         h1->Draw("colz");
         //hist->Draw("colz");
         cv->cd(2);
         gPad->SetLeftMargin(0.15);
         gPad->SetBottomMargin(0.15);
         gPad->SetRightMargin(0.15);
         h2->Draw("colz");
         //hist2->Draw("colz");
         cv->SaveAs("Plots/testSub2.pdf");
         }
         
         RooDataHist* dataHist1 = new RooDataHist("roo_hist1", "roo_hist1", RooArgSet(*m,*pt), h1);
         RooDataHist* dataHist2 = new RooDataHist("roo_hist2", "roo_hist2", RooArgSet(*m,*pt), h2);
         RooHistPdf* modelHist2 = new RooHistPdf("roo_hist_pdf2", "roo_hist_pdf2", RooArgSet(*m,*pt), *dataHist2);
         */
        /*
         TCanvas* cv = new TCanvas();
         gPad->SetLeftMargin(0.15);
         gPad->SetBottomMargin(0.15);
         gPad->SetRightMargin(0.15);
         hist->Draw("colz");
         cv->SaveAs(Form("Plots/hist-chi-%d.pdf", k-minBinPt+1));
         */
        
        
        //RooChi2Var* chi2Var = new RooChi2Var("chi2","chi2", *fitModel, *roo_hist, DataError(RooAbsData::Poisson));
        RooChi2Var* chi2Var = new RooChi2Var("chi2","chi2", *roo_hist_pdf, *roo_hist, DataError(RooAbsData::Poisson));
        mBinVector.push_back(mBinNumber);
        ptBinVector.push_back(ptBinNumber);
        productVector.push_back(mBinNumber*ptBinNumber);
        chi2Vec.push_back((double)chi2Var->getVal()/ndf);
        
        std::cout << endl << "ptBinNumber = " << ptBinNumber << endl;
        std::cout << "mBinNumber = " << mBinNumber << endl;
        std::cout << "chi2Var = " << (double)chi2Var->getVal()/ndf << std::endl << endl;
        
        ptBinNumber--;
        delete hist;
    }
    
    int maxProductIndex = std::max_element(productVector.begin(),productVector.end()) - productVector.begin();
    int minChi2Index = std::min_element(chi2Vec.begin(),chi2Vec.end()) - chi2Vec.begin();
    /*
     std::cout << "Maximum  number of bins = " << productVector[maxProductIndex] << endl;
     std::cout << "Number of bins in M = " << mBinVector[maxProductIndex] << endl;
     std::cout << "Number of bins in Pt = " << ptBinVector[maxProductIndex] << endl;
     
     std::cout << "Min Chi2 = " << chi2Vec[minChi2Index] << endl;
     std::cout << "Total number of bins = " << productVector[minChi2Index] << endl;
     std::cout << "Number of bins in M = " << mBinVector[minChi2Index] << endl;
     std::cout << "Number of bins in Pt = " << ptBinVector[minChi2Index] << endl;
     */
    
    std::vector<int> binNumbersVec = {};
    binNumbersVec.push_back(mBinVector[minChi2Index]);
    binNumbersVec.push_back(ptBinVector[minChi2Index]);
    
    return binNumbersVec;
    
}

void MakePlots(RooWorkspace* ws, string period, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax, bool drawPulls, bool logScale, bool exp, bool exclusiveOnly, double& chi2, bool &noInclusive) {
    
    //int ptBinNumber = int(15*(ptMax-ptMin));
    //int ptBinNumber = int(10*(ptMax-ptMin));
    //int mBinNumber = int(20*(mMax-mMin));
    //int ptBinNumber = 20;
    
    bool eps = true;    // save in .eps
    
    // get what we need of the workspace and fit
    RooRealVar* m = ws->var("fTrkTrkM");
    RooRealVar* pt = ws->var("fTrkTrkPt");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    RooAbsPdf* fitModel = ws->pdf("model");
    
    
    // Fit data
    RooFitResult* r = fitModel->fitTo(*data, Extended(), Minos(true), Strategy(1), Save());
    //RooFitResult* r = fitModel->fitTo(*h_data, Extended(), Minos(true), Strategy(1), Save());
    //RooFitResult* r = fitModel->fitTo(*data, Minos(false), Strategy(0), Save());    // be quick (for testing)
    //RooFitResult* r = nullptr;
    
    // Find chi2 and bin numbers
    int nFitParam = fitModel->getParameters(data)->selectByAttrib("Constant",kFALSE)->getSize();
    
    bool computeChi2 = false;
    if (computeChi2) {
        // To search for bin numbers
        std::vector <int> binNumbersVec = FindBinNumbers(ws);
        int mBinNumber = binNumbersVec[0];
        int ptBinNumber = binNumbersVec[1];
        std::cout << " mBinNumber = " << mBinNumber << std::endl;
        std::cout << " ptBinNumber = " << ptBinNumber << std::endl;
        //return;
        
        
        /*
         // To force bin numbers
         int mBinNumber = 5;
         int ptBinNumber = 6;
         TH1* hist = data->createHistogram("fTrkTrkM,fTrkTrkPt", mBinNumber, ptBinNumber);
         while (hist->GetMinimum() < 1) {
         mBinNumber--;
         delete hist;
         hist = data->createHistogram("fTrkTrkM,fTrkTrkPt", mBinNumber, ptBinNumber);
         if (hist->GetMinimum() < 1) {
         ptBinNumber--;
         delete hist;
         hist = data->createHistogram("fTrkTrkM,fTrkTrkPt", mBinNumber, ptBinNumber);
         }
         }
         */
        
        TH1* hist_data2 = data->createHistogram("fTrkTrkM,fTrkTrkPt", mBinNumber, ptBinNumber);
        RooDataHist* h_data2 = new RooDataHist("h_data2", "h_data2", RooArgSet(*m,*pt), hist_data2);
        
        int ndf = ndf = mBinNumber*ptBinNumber-nFitParam;
        // Chi2
        TH1* hist_model2 = fitModel->createHistogram("fTrkTrkM,fTrkTrkPt", mBinNumber, ptBinNumber);
        RooDataHist* h_model2 = new RooDataHist("h_model2", "h_model", RooArgSet(*m,*pt), hist_model2);
        RooHistPdf* roo_hist_pdf = new RooHistPdf("roo_hist_pdf", "roo_hist_pdf", RooArgSet(*m,*pt), *h_model2);
        
        //RooChi2Var* chi2Var = new RooChi2Var("chi2","chi2", *fitModel, *h_data2, DataError(RooAbsData::Poisson));  // IntegrateBins(0.01)
        RooChi2Var* chi2Var = new RooChi2Var("chi2","chi2", *roo_hist_pdf, *h_data2, DataError(RooAbsData::Poisson));  // IntegrateBins(0.01)
        std::cout << "chi2/ndf (with histogram) = " << chi2Var->getVal() << "/" << ndf << " = " << (double)chi2Var->getVal()/ndf << endl;
        
        std::cout << "Maximum  number of bins = " << mBinNumber*ptBinNumber << endl;
        std::cout << "Number of bins in M = " << mBinNumber << endl;
        std::cout << "Number of bins in Pt = " << ptBinNumber << endl;
        // Plot things
        
        TCanvas* cv0 = new TCanvas();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        TH2F* hTest = (TH2F*)hist_data2->Clone();
        hTest->Add(hist_model2, -1);
        hTest->Draw("colz");
        //hist_data2->Draw("colz");
        cv0->SaveAs("Plots/test.pdf");
        //return;
    }
    
    /*
     if (mBinNumber < 2 || ptBinNumber < 2) {
     std::cout << "Problem! Number of bins in 1D < 2!" << endl;
     return;
     }
     */
    
    
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
    
    /*
     if (!noInclusive) {
     if (yieldJpsiInclusive->getVal()<5) noInclusive = true;
     }
     */
    
    RooRealVar* bExc = ws->var("bExc");
    RooRealVar* bDiss = ws->var("bDiss");
    RooRealVar* nDiss = nullptr;
    if (!exp) nDiss = ws->var("nDiss");
    
    RooRealVar* meanGg = ws->var("mean_gg");
    RooRealVar* sigmaGg = ws->var("sigma_gg");
    
    // Define mass frame
    RooPlot* mframe = m->frame(Title("Fit of invariant mass"));
    int binNum = 50;
    if (exclusiveOnly) binNum = 30;
    data->plotOn(mframe, Binning(binNum));
    //data->plotOn(mframe);
    fitModel->plotOn(mframe, Name("sum"), LineColor(kRed), LineWidth(2));
    fitModel->plotOn(mframe,Name("pdfJpsiExclusive"),Components(*pdfJpsiExclusive),LineStyle(kDashed), LineColor(kOrange+8), LineWidth(2));
    if (!exclusiveOnly) { fitModel->plotOn(mframe,Name("pdfJpsiDissociative"),Components(*pdfJpsiDissociative),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(2));
    }
    fitModel->plotOn(mframe,Name("pdfBackground"),Components(*pdfBackground),LineStyle(kDashed), LineColor(kBlack), LineWidth(2));
    fitModel->plotOn(mframe,Name("pdfJpsiGammaPb"),Components(*pdfJpsiGammaPb),LineStyle(kDashed), LineColor(4), LineWidth(2));
    if (mMax > mLimitPsi2s) fitModel->plotOn(mframe,Name("pdfPsi2s"),Components(*pdfPsi2s),LineStyle(kDashed), LineColor(6), LineWidth(2));
    fitModel->plotOn(mframe,Name("pdfTwoGamma"),Components(*pdfTwoGamma),LineStyle(kDashed), LineColor(kMagenta), LineWidth(2));
    
    if (!noInclusive) fitModel->plotOn(mframe,Name("pdfJpsiInclusive"),Components(*pdfJpsiInclusive),LineStyle(kDashed), LineColor(kMagenta+2), LineWidth(2));
    
    // Define pt frame
    vector<RooPlot*> ptframes= {pt->frame(Title("Fit of p_{T} (log scale and full p_{T} range)")), pt->frame(Title("Fit of p_{T} (zoom)"))};
    for (int k = 0; k<(int)ptframes.size(); k++) {
        //for (int k = 0; k< 1; k++) {
        if (exclusiveOnly) ptframes[k]->SetTitle("Fit of p_{T}");
        data->plotOn(ptframes[k], Binning(50));
        //data->plotOn(ptframes[k]);
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
    }
    
    /*
     TCanvas* cv = new TCanvas("cv", "cv", 600, 300);
     cv->Divide(2);
     cv->cd(1);
     mframe->Draw();
     cv->cd(2);
     ptframes[0]->Draw();
     
     cv->SaveAs("Plots/LHC16r/test.pdf");
     return;
     */
    
    
    /*
     Double_t chi2M = mframe->chiSquare( "sum", "h_data", nFitParam);
     std::cout << "chi2/ndf (method mframe->chiSquare) = " << chi2M << endl;
     //return;
     */
    
    TCanvas* c1 = new TCanvas("2Dplot","2D fit",1800,1000) ;
    //TCanvas* c1 = new TCanvas("2Dplot","2D fit",800,300) ;
    if (drawPulls) c1->Divide(4,3) ;
    else {c1->Divide(4,2) ; c1->SetCanvasSize(1800, 600);}
    
    std::cout << std::endl << std::endl << "Number of degrees of freedom = " << nFitParam << std::endl << std::endl << std::endl;
    // Mass Plot
    c1->cd(2) ;
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    if (logScale) gPad->SetLogy() ;
    //mframe->addObject(txtChi) ;
    //mframe->chiSquare() ;
    /*
     TLatex* txtbinNumM = new TLatex(2.6, 0.5*mframe->GetMaximum(), Form("#bf{binNum = %d}", mBinNumber));
     mframe->addObject(txtbinNumM);
     */
    mframe->Draw();
    //std::cout << std::endl << "chi2 = " << chi2 << std::endl;
    
    // pt plot
    c1->cd(3) ;
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    //if (logScale) gPad->SetLogy() ;
    double yMax2 = ptframes[0]->GetMaximum();
    double y1 = 0.75*yMax2, y2 = 0.65*yMax2, y3 = 0.57*yMax2, y4 = 0.5*yMax2;
    if (logScale) {y1 = yMax2/pow(2.,1), y2 = yMax2/pow(2.,2), y3 = yMax2/pow(2.,2.8), y4 = yMax2/pow(2.,4);}
    double ptRangeMax=ptMax;
    ptframes[0]->GetXaxis()->SetRangeUser(0, ptRangeMax);
    ptframes[0]->Draw();
    ptframes[0]->Print();
    
    // Zoom pt
    c1->cd(4);
    //gPad->SetLogy();
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    if (!exclusiveOnly) ptRangeMax = 3;
    double xText = ptMin+(ptRangeMax-ptMin)*1/2;
    //chi2 = ptframes[0]->chiSquare("sum", "h_data", nFitParam);
    //TLatex* txtChi2 = new TLatex(xText, y4,Form("#chi^{2}/ndf = %.1f / %d = %.3f", chi2*nFitParam, nFitParam, chi2));
    //TLatex* txtChi2 = new TLatex(xText, y4, chi2String);
    TLatex* txtExc = new TLatex(xText, y1, Form("b_{exc} = %.2f #pm %.2f", bExc->getVal(), bExc->getError()));
    ptframes[1]->addObject(txtExc) ;
    //ptframes[1]->addObject(txtChi2) ;
    ptframes[1]->GetXaxis()->SetRangeUser(0, ptRangeMax);
    //TLatex* txtbinNumPt = new TLatex(0.5, 0.8*ptframes[1]->GetMaximum(), Form("#bf{binNum = %d}", ptBinNumber));
    //ptframes[1]->addObject(txtbinNumPt);
    ptframes[1]->Draw();
    
    /*
     // Quality plots: (data-fit)/sigma
     c1->cd(6);
     gPad->SetLeftMargin(0.15) ;
     //gPad->SetTopMargin(0.15) ;
     data->plotOn(mframe, Binning(50));
     fitModel->plotOn(mframe, Binning(50));
     RooHist* hpullM = mframe->pullHist();
     hpullM->GetXaxis()->SetRangeUser(mMin, mMax);
     hpullM->SetTitle("(data - fit)/#sigma for m distribution");
     hpullM->Draw("");
     
     c1->cd(7);
     gPad->SetLeftMargin(0.15) ;
     data->plotOn(ptframes[k], Binning(ptBinNumber));
     fitModel->plotOn(ptframes[0], Binning(ptBinNumber));
     RooHist* hpullPt = ptframes[0]->pullHist();
     //hpullPt->GetXaxis()->SetRangeUser(ptMin, ptMax);
     hpullPt->GetXaxis()->SetRangeUser(0,ptRangeMax);
     hpullPt->SetTitle("(data - fit)/#sigma for p_{T} distribution");
     hpullPt->Draw("");
     */
    
    // for naming output files
    string cutType = "";
    if (!useCuts) cutType = "-nocuts";
    string suf;
    if (exp) suf = "exp";
    else suf = "powerlaw";
    if (exclusiveOnly) suf += "-exclusive-only";
    
    //Correlation coefficients
    ofstream myfile;
    TString fileNameCor = Form("correlation-coeff-%s-%.1f-%.1f-%s-%.1f-%.1f.txt", period.c_str(), mMin, mMax, suf.c_str(), abs(rapMax), abs(rapMin));
    const string fileNameCorString = fileNameCor.Data();
    myfile.open(fileNameCorString);
    
    c1->cd(6);
    gStyle->SetTextSize(0.05);
    
    RooArgList argList = r->floatParsInit();
    //return;
    bool cor = false;
    int l = 0;
    vector <TText*> txtVec = {};
    for (int i = 0; i<nFitParam; i++) {
        RooAbsArg* arg1 = argList.at(i);
        for (int j = i+1; j < nFitParam; j++) {
            RooAbsArg* arg2 = argList.at(j);
            double correl = r->correlation(*arg1, *arg2);
            if (abs(correl) > 0.3) {
                cor=true;
                if (abs(correl) > 0.7) gStyle->SetTextColor(kRed);
                else if (abs(correl) > 0.5) gStyle->SetTextColor(kOrange+7);
                else gStyle->SetTextColor(kBlack);
                //cout << "\n\n\n" << Form("Correlation between %s and %s = %f",  arg1->GetName(), arg2->GetName(), correl) << "\n\n\n" << endl;
                l++;
                int l2 = l;
                if (l>=10) {l2=l-10; c1->cd(7);}
                else c1->cd(6);
                TText* txxx = new TText(0.1, 0.9-l2*0.07, Form("%s and %s = %.2f",  arg1->GetName(), arg2->GetName(), correl) );
                //TText txxx(0.2, 0.3, Form("Correlation between %d and %d", i, j) );
                txtVec.push_back(txxx);
                txxx->Draw();
            }
            string lineCor = Form("%s and %s = %.4f",  arg1->GetName(), arg2->GetName(), correl);
            myfile << lineCor << "\n";
        }
    }
    myfile.close();
    //return;
    
    if (cor) {
        c1->cd(6);
        TText* txxx = new TText(0.1, 0.9, "Correlation between:" );
        txxx->Draw();
    }
    
    gStyle->SetTextSize(0.07);
    gStyle->SetTextColor(kBlack);
    
    // Quality plots: data-fit
    if (drawPulls) {
        c1->cd(10);
        gPad->SetLeftMargin(0.15) ;
        RooHist* hresidM = mframe->residHist();
        hresidM->GetXaxis()->SetRangeUser(mMin, mMax);
        hresidM->SetTitle("Residuals (data - fit) for m distribution");
        hresidM->Draw("");
        
        c1->cd(11);
        gPad->SetLeftMargin(0.15) ;
        RooHist* hresidPt = ptframes[0]->residHist();
        //hresidPt->GetXaxis()->SetRangeUser(ptMin, ptMax);
        hresidPt->GetXaxis()->SetRangeUser(0, 4);
        hresidPt->SetTitle("Residuals (data - fit) for p_{T} distribution");
        hresidPt->Draw("");
    }
    //return;
    // Legend in a subcanvas
    c1->cd(1);
    TLegend* legend = new TLegend(0.1, 0.3, 0.9, 0.9);
    legend->SetFillColor(kWhite);
    legend->SetLineColor(kWhite);
    legend->AddEntry(ptframes[0]->findObject("pdfJpsiExclusive"), "Exclusive J/#Psi","L");
    if (!exclusiveOnly) {
        legend->AddEntry(ptframes[0]->findObject("pdfJpsiDissociative"), "J/#Psi with dissociative p","L");
    }
    if (!noInclusive) legend->AddEntry(ptframes[0]->findObject("pdfJpsiInclusive"), "Inclusive events","L");
    legend->AddEntry(ptframes[0]->findObject("pdfJpsiGammaPb"), "#gamma + Pb","L");
    if (mMax > mLimitPsi2s) legend->AddEntry(ptframes[0]->findObject("pdfPsi2s"), "Psi(2s)","L");
    legend->AddEntry(ptframes[0]->findObject("pdfTwoGamma"), "#gamma#gamma #rightarrow #mu^{+} #mu^{-}","L");
    legend->AddEntry(ptframes[0]->findObject("pdfBackground"), "extra background","L");
    legend->AddEntry(ptframes[0]->findObject("sum"),"sum","L");
    legend->Draw();
    //return;
    // Number of different contributions in a subcanvas
    c1->cd(5);
    
    double yBkg = 0.5;
    // Write number of candidates
    //TLatex* txt1 = new TLatex(3.3,0.9*yMax,Form("Coherent Jpsi : %.1f", yieldCohJpsi.getVal()));
    TLatex* txt2 = new TLatex(0.2,0.9,Form("Exclusive J/#Psi : %.1f #pm %.1f", yieldJpsiExclusive->getVal(), yieldJpsiExclusive->getError()));
    if (!exclusiveOnly) {
        TLatex* txt3 = new TLatex(0.2,0.8,Form("Dissociative J/#Psi : %.1f #pm %.1f", yieldJpsiDissociative->getVal(), yieldJpsiDissociative->getError()));
        txt3->Draw();
    }
    TLatex* txt9 = new TLatex(0.2,yBkg-0.1,Form("Other background : %.1f #pm %.1f", yieldBkg->getVal(), yieldBkg->getError()));
    txt9->Draw();
    TLatex* txt4 = new TLatex(0.2,0.7,Form("#gamma-Pb J/#Psi : %.1f #pm %.1f", yieldJpsiGammaPb->getVal(), yieldJpsiGammaPb->getError()));
    if (mMax > mLimitPsi2s) {
        TLatex* txt5 = new TLatex(0.2,yBkg,Form("#Psi(2s) : %.1f #pm %.1f", yieldPsi2s->getVal(), yieldPsi2s->getError()));
        txt5->Draw();
        yBkg = 0.4;
    }
    TLatex* txt6 = new TLatex(0.2,yBkg,Form("#gamma#gamma #rightarrow #mu^{+} #mu^{-} : %.1f #pm %.1f", yieldTwoGamma->getVal(), yieldTwoGamma->getError()));
    txt2->Draw(); txt4->Draw(); txt6->Draw();
    if (!noInclusive) {
        TLatex* txt4bis = new TLatex(0.2,0.6,Form("Inclusive J/#Psi : %.1f #pm %.1f", yieldJpsiInclusive->getVal(), yieldJpsiInclusive->getError()));
        txt4bis->Draw();
    }
    /*
     // Chi2
     //double xChi = mMin + (mMax-mMin)*0.7;
     //Double_t chi2M = mframe->chiSquare( "sum", "h_data", nFitParam);
     //TLatex* txtChi = new TLatex(xChi, 0.5*mframe->GetMaximum(),Form("#bf{#chi^{2}/ndf = %.2f / %d = %.2f}", chi2M*nFitParam, nFitParam, chi2M) );
     TString chi2String = Form("#bf{#chi^{2}/ndf = %.1f/%d = %.2f}", chi2Var->getVal(), ndf, (double)chi2Var->getVal()/ndf);
     TLatex* txtChi = new TLatex(0.2, 0.6, chi2String);
     if (noInclusive) txtChi->Draw();
     */
    
    // Compute ratios N_diss/N_exc and N_(gamma-Pb)/N_exc
    string percent = "%";
    RooFormulaVar r_diss_exc0("r_diss_exc0","yieldJpsiDissociative/yieldJpsiExclusive",RooArgSet(*yieldJpsiDissociative, *yieldJpsiExclusive));
    RooRealVar *r_diss_exc = new RooRealVar("r_diss_exc", "r_diss_exc", r_diss_exc0.getVal());
    r_diss_exc->setError(r_diss_exc0.getPropagatedError(*r));
    ws->import(*r_diss_exc);
    if (!exclusiveOnly) {
        TLatex* txt7 = new TLatex(0.2,0.3,Form("N_{diss}/N_{exc} = %.2f #pm %.2f %s", r_diss_exc0.getVal()*100, r_diss_exc0.getPropagatedError(*r)*100, percent.c_str()));
        txt7->Draw();
    }
    RooFormulaVar r_gammaPb_exc("r_gammaPb_exc","yieldJpsiGammaPb/yieldJpsiExclusive",RooArgSet(*yieldJpsiGammaPb, *yieldJpsiExclusive));
    TLatex* txt8 = new TLatex(0.2,0.2,Form("N_{#gamma-Pb}/N_{exc} = %.2f #pm %.2f %s", r_gammaPb_exc.getVal()*100, r_gammaPb_exc.getPropagatedError(*r)*100, percent.c_str()));
    txt8->Draw();
    
    RooRealVar *r_diss_exc2 = new RooRealVar("r_diss_exc2", "r_diss_exc2", r_gammaPb_exc.getVal());
    r_diss_exc2->setError(r_gammaPb_exc.getPropagatedError(*r));
    if (mMin < 3.2 && mMax > 3) ws->import(*r_diss_exc2);
    
    c1->cd(8);
    if (false) {
        //if (!noInclusive) {
        RooRealVar* nInc = ws->var("nInc");
        RooRealVar* pt0 = ws->var("pt0");
        TLatex* tInc0 = new TLatex(0.1, 0.9, "Inclusive: #frac{p_{T}}{(1 + (p_{T}/p_{0})^{2})^{n}} with");
        TLatex* tInc1 = new TLatex(0.1, 0.75, Form("n_{inc} = %.2f #pm %.2f", nInc->getVal(), nInc->getError()));
        TLatex* tInc2 = new TLatex(0.1, 0.65, Form("p_{0} = %.2f #pm %.2f", pt0->getVal(), pt0->getError()));
        tInc0->Draw(); tInc1->Draw(); tInc2->Draw();
    }
    else {
        TLatex* tGamma = new TLatex(0.1, 0.9, "Landau parameters for #gamma#gamma:");
        TLatex* tGamma1 = new TLatex(0.1, 0.75, Form("#mu = %.4f #pm %.4f", meanGg->getVal(), meanGg->getError()));
        TLatex* tGamma2 = new TLatex(0.1, 0.65, Form("#sigma = %.4f #pm %.4f", sigmaGg->getVal(), sigmaGg->getError()));
        tGamma->Draw(); tGamma1->Draw(); tGamma2->Draw();
    }
    
    
    TLatex* tDiss0 = nullptr;
    TLatex* tDiss1 = new TLatex(0.1, 0.35,Form("b_{diss} = %.2f #pm %.2f", bDiss->getVal(), bDiss->getError()));
    if (exp) {
        tDiss0 = new TLatex(0.1, 0.45, "Dissociative: p_{T}*exp(-b_{diss}*p_{T}^{2}) with");
    }
    else {
        tDiss0 = new TLatex(0.1, 0.45, "Dissociative: p_{T}*(1+p_{T}^{2}*b_{diss}/n_{diss})^{-n_{diss}} with");
        TLatex* tDiss2 = new TLatex(0.1, 0.25,Form("n_{diss} = %.2f #pm %.2f", nDiss->getVal(), nDiss->getError()));
        if (!exclusiveOnly) tDiss2->Draw();
    }
    if (!exclusiveOnly) {
        tDiss0->Draw(); tDiss1->Draw();
    }
    
    
    // save plot
    c1->SaveAs(Form("Plots/%s/Fit-2D%s-%.1f-%.1f-%s-%.1f-%.1f.pdf", period.c_str(), cutType.c_str(), mMin, mMax, suf.c_str(), abs(rapMax), abs(rapMin)));
    if (eps) c1->SaveAs(Form("Plots/%s/Fit-2D%s-%.1f-%.1f-%s-%.1f-%.1f.eps", period.c_str(), cutType.c_str(), mMin, mMax, suf.c_str(), abs(rapMax), abs(rapMin)));
    
    RooRealVar* pt0 = ws->var("pt0");
    RooRealVar* nInc = ws->var("nInc");
    
    /*
     cout << "\n\n\n\npt0 = " << pt0->getVal() << endl;
     cout << "nInc = " << nInc->getVal() << "\n\n\n" << endl;
     */
    //return;
    TLegend* legend2;
    if (exclusiveOnly) legend2 = new TLegend(0.48, 0.7, 1, 1);
    else legend2 = new TLegend(0.55, 0.65, 1, 1);
    legend2->SetFillColor(kWhite);
    legend2->SetLineColor(kWhite);
    legend2->AddEntry(ptframes[0]->findObject("pdfJpsiExclusive"), "Exclusive J/#Psi","L");
    if (!exclusiveOnly) {legend2->AddEntry(ptframes[0]->findObject("pdfJpsiDissociative"), "J/#Psi with dissociative p","L");}
    if (!noInclusive) legend2->AddEntry(ptframes[0]->findObject("pdfJpsiInclusive"), "Inclusive events","L");
    legend2->AddEntry(ptframes[0]->findObject("pdfJpsiGammaPb"), "#gamma + Pb","L");
    if (mMax > mLimitPsi2s) legend2->AddEntry(ptframes[0]->findObject("pdfPsi2s"), "Psi(2s)","L");
    legend2->AddEntry(ptframes[0]->findObject("pdfTwoGamma"), "#gamma#gamma #rightarrow #mu^{+} #mu^{-}","L");
    legend2->AddEntry(ptframes[0]->findObject("pdfBackground"), "extra background","L");
    legend2->AddEntry(ptframes[0]->findObject("sum"),"sum","L");
    
    //return;
    
    TCanvas* c2 = new TCanvas("2Dplot","2D fit",1200,500) ;
    c2->Divide(2);
    // Draw mass plot
    c2->cd(1);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    //TLatex* txtChiSimple = new TLatex(2.55, 0.87*mframe->GetMaximum(),Form("#bf{#chi^{2}/ndf = %.2f / %d = %.2f}", chi2M*nFitParam, nFitParam, chi2M) );
    //TLatex* txtChiSimple = new TLatex(2.6, 0.87*mframe->GetMaximum(),Form("#bf{#chi^{2}/ndf = %.2f}", chi2M) );
    /*
     TLatex* txtChiSimple = new TLatex(2.55, mframe->GetMaximum(),chi2String);
     txtChiSimple->SetTextSize(0.045);
     */
    mframe->remove();
    //mframe->addObject(txtChiSimple) ;
    // Add number of J/Psis
    TLatex* txtJpsi = new TLatex(2.6, 0.72*mframe->GetMaximum(), Form("#bf{N_{exc J/#Psi} = %.0f #pm %.0f}", yieldJpsiExclusive->getVal(), yieldJpsiExclusive->getError()));
    txtJpsi->SetTextSize(0.045);
    mframe->addObject(txtJpsi) ;
    TLatex* txtJpsi2 = new TLatex(2.6, 0.62*mframe->GetMaximum(), Form("#bf{N_{diss J/#Psi} = %.0f #pm %.0f}", yieldJpsiDissociative->getVal(), yieldJpsiDissociative->getError()));
    txtJpsi2->SetTextSize(0.045);
    if (!exclusiveOnly) mframe->addObject(txtJpsi2) ;
    TLatex* txtJpsi3 = new TLatex(2.6, 0.52*mframe->GetMaximum(), Form("#bf{N_{#gamma#gamma} = %.0f #pm %.0f}", yieldTwoGamma->getVal(), yieldTwoGamma->getError()));
    txtJpsi3->SetTextSize(0.045);
    mframe->addObject(txtJpsi3) ;
    TLatex* txtJpsi4 = new TLatex(2.6, 0.42*mframe->GetMaximum(), Form("#bf{N_{bkg} = %.0f #pm %.0f}", yieldBkg->getVal(), yieldBkg->getError()));
    txtJpsi4->SetTextSize(0.045);
    mframe->addObject(txtJpsi4) ;
    mframe->SetTitle("");
    mframe->SetMaximum(mframe->GetMaximum()*1.1);
    /*
     TLatex* txtRap = new TLatex(2.6, 0.93*mframe->GetMaximum(), Form("#bf{%.2f < y < %.2f, %.1f < p_{T} < %.1f GeV/c}", -rapMax, -rapMin, ptMin, ptMax));
     TLatex* txtMass = new TLatex(2.6, 0.85*mframe->GetMaximum(), Form("#bf{%.1f < M < %.1f GeV/c^{2}}", mMin, mMax));
     txtRap->SetTextSize(0.042);
     txtMass->SetTextSize(0.042);
     mframe->addObject(txtRap) ;
     mframe->addObject(txtMass) ;
     */
    
    mframe->Draw();
    //return;
    // write y and pt info
    TLegend* lgdInfo = new TLegend(0.15, 0.73, 0.45, 0.9);
    lgdInfo->SetMargin(0.05);
    lgdInfo->SetTextSize(0.04);
    lgdInfo->AddEntry((TObject*)0, Form("%.2f < y < %.2f", -rapMax, -rapMin), "");
    lgdInfo->AddEntry((TObject*)0, Form("%.1f < M < %.1f GeV/c^{2}", mMin, mMax), "");
    lgdInfo->AddEntry((TObject*)0, Form("%.1f < p_{T} < %.1f GeV/c", ptMin, ptMax), "");
    lgdInfo->Draw();
    
    
    // Draw pT plot
    c2->cd(2);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    
    ptframes[1]->remove();
    ptframes[1]->remove();
    ptframes[1]->SetTitle("");
    TLatex* txtBexc = nullptr;
    //TLatex* txtChiPt = nullptr;
    if (exclusiveOnly) {
        txtBexc = new TLatex(0.6, 0.7*ptframes[1]->GetMaximum(), Form("#bf{b_{exc} = %.2f #pm %.2f}", bExc->getVal(), bExc->getError()));
        //txtChiPt = new TLatex(0.64, 0.55*ptframes[1]->GetMaximum(),Form("#bf{#chi^{2}/ndf = %.2f / %d = %.2f}", chi2*nFitParam, nFitParam, chi2) );
        //txtChiPt = new TLatex(0.7, 0.55*ptframes[1]->GetMaximum(),Form("#bf{#chi^{2}/ndf = %.2f}", chi2) );
        //txtChiPt = new TLatex(0.6, 0.6*ptframes[1]->GetMaximum(), chi2String );
    }
    else {
        txtBexc = new TLatex(1.2, 0.6*ptframes[1]->GetMaximum(), Form("#bf{b_{exc} = %.2f #pm %.2f}", bExc->getVal(), bExc->getError()));
        //txtChiPt = new TLatex(1.4, 0.4*ptframes[1]->GetMaximum(),Form("#bf{#chi^{2}/ndf = %.2f / %d = %.2f}", chi2*nFitParam, nFitParam, chi2) );
        //txtChiPt = new TLatex(1.5, 0.4*ptframes[1]->GetMaximum(),Form("#bf{#chi^{2}/ndf = %.2f}", chi2) );
        //txtChiPt = new TLatex(1.2, 0.5*ptframes[1]->GetMaximum(), chi2String);
    }
    txtBexc->SetTextSize(0.045);
    ptframes[1]->SetMaximum(ptframes[1]->GetMaximum()*1.1);
    ptframes[1]->addObject(txtBexc);
    //txtChiPt->SetTextSize(0.045);
    //ptframes[1]->addObject(txtChiPt);
    
    ptframes[1]->Draw();
    legend2->Draw();
    
    // write y and pt info
    
    
    c2->SaveAs(Form("Plots/%s/Fit-simple-2D%s-%.1f-%.1f-%s-%.1f-%.1f.pdf", period.c_str(), cutType.c_str(), mMin, mMax, suf.c_str(), abs(rapMax), abs(rapMin)));
    if (eps) c2->SaveAs(Form("Plots/%s/Fit-simple-2D%s-%.1f-%.1f-%s-%.1f-%.1f.eps", period.c_str(), cutType.c_str(), mMin, mMax, suf.c_str(), abs(rapMax), abs(rapMin)));
    
    /*
     std::cout << "chi2/ndf (with histogram) = " << chi2Var->getVal() << "/" << ndf << " = " << (double)chi2Var->getVal()/ndf << endl;
     std::cout << "mBinNumber = " << mBinNumber << endl;
     std::cout << "ptBinNumber = " << ptBinNumber << endl;
     */
    
}

TH1* createHistoFromPdf(RooAbsPdf* pdf, double yield, const char * name, RooRealVar* var, int nBins, double min, double max, int nBinsIni) {
    TH1* h = pdf->createHistogram(name,*var, Binning(nBins, min, max));
    double factor = (double)nBins/nBinsIni;
    h->Scale(1./h->Integral());
    //cout << endl << "name = " << name << endl;
    //cout << "h->Integral() = " << h->Integral() << endl;
    //h->Scale((double)factor/h->Integral());
    h->Scale(yield*(double)factor/h->Integral());
    for (int i = 0; i < nBins; i++) {h->SetBinError(i+1, 0);}
    h->SetOption("hist");
    return h;
}

void ExportHist(RooWorkspace* ws, bool exclusiveOnly) {
    
    string suffix = "";
    if (exclusiveOnly) suffix = "-exc";
    // First get quantities needed
    
    RooDataSet* data = (RooDataSet*) ws->data("data");
    RooRealVar* m = ws->var("fTrkTrkM");
    RooRealVar* pt = ws->var("fTrkTrkPt");
    
    RooAbsPdf* fitModel = ws->pdf("model");
    
    RooAbsPdf* pdfJpsiExclusive = ws->pdf("pdfJpsiExclusive");
    RooAbsPdf* pdfJpsiDissociative = nullptr;
    if (!exclusiveOnly) pdfJpsiDissociative = ws->pdf("pdfJpsiDissociative");
    RooAbsPdf* pdfJpsiGammaPb = ws->pdf("pdfJpsiGammaPb");
    RooAbsPdf* pdfTwoGamma = ws->pdf("pdfTwoGamma");
    RooAbsPdf* pdfBackground = ws->pdf("pdfBackground");
    /*
     RooAbsPdf* pdfJpsiInclusive = ws->pdf("pdfJpsiInclusive");
     RooAbsPdf* pdfPsi2s = ws->pdf("pdfPsi2s");
     */
    
    RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
    RooRealVar* yieldJpsiDissociative = new RooRealVar("yieldJpsiDissociative", "yieldJpsiDissociative", 0);
    if (!exclusiveOnly) yieldJpsiDissociative = ws->var("yieldJpsiDissociative");
    RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
    RooRealVar* yieldTwoGamma = ws->var("yieldTwoGamma");
    RooRealVar* yieldBkg = ws->var("yieldBkg");
    /*
     RooRealVar* yieldJpsiInclusive = ws->var("yieldJpsiInclusive");
     RooRealVar* yieldPsi2s = ws->var("yieldPsi2s");
     */
    
    // number of bins
    int nBinsM = 50;
    int nBinsPt = 50;
    if (exclusiveOnly) {
        //nBinsM = 100;
        nBinsPt = 100;
    }
    double ptMin = 0.00000000000001, ptMax = 3.;
    double mMin = 2.5, mMax = 3.5;
    
    // Create histograms for data
    TH1* hMData = data->createHistogram("hMData",*m, Binning(nBinsM, mMin, mMax));
    TH1* hPtData = data->createHistogram("hPtData",*pt, Binning(nBinsPt, ptMin, ptMax));
    
    // Create histograms for model
    int nBinsHistPdf = int((ptMax-ptMin)*500)+1;
    //int nBinsModel=1000;
    int nBinsModel=nBinsHistPdf;
    // First sum
    double totalInt = yieldJpsiExclusive->getVal() + yieldJpsiDissociative->getVal() + yieldJpsiGammaPb->getVal() + yieldTwoGamma->getVal() + yieldBkg->getVal();
    cout << endl << endl << "totalInt = " << totalInt << endl << endl;
    TH1* hMModel = createHistoFromPdf(fitModel, totalInt, "hMModel", m, nBinsModel, mMin, mMax, nBinsM);
    TH1* hPtModel = createHistoFromPdf(fitModel, totalInt, "hPtModel", pt, nBinsHistPdf, ptMin, ptMax, nBinsPt);
    
    // Now components
    
    TH1* hMJpsiExclusive = createHistoFromPdf(pdfJpsiExclusive, yieldJpsiExclusive->getVal(), "hMJpsiExclusive", m, nBinsModel, mMin, mMax, nBinsM);
    TH1* hMJpsiDissociative = nullptr;
    if (!exclusiveOnly) hMJpsiDissociative = createHistoFromPdf(pdfJpsiDissociative, yieldJpsiDissociative->getVal(), "hMJpsiDissociative", m, nBinsModel, mMin, mMax, nBinsM);
    TH1* hMJpsiGammaPb = createHistoFromPdf(pdfJpsiGammaPb, yieldJpsiGammaPb->getVal(), "hMJpsiGammaPb", m, nBinsModel, mMin, mMax, nBinsM);
    TH1* hMTwoGamma = createHistoFromPdf(pdfTwoGamma, yieldTwoGamma->getVal(), "hMTwoGamma", m, nBinsModel, mMin, mMax, nBinsM);
    TH1* hMBackground = createHistoFromPdf(pdfBackground, yieldBkg->getVal(), "hMBackground", m, nBinsModel, mMin, mMax, nBinsM);
    /*
     TH1* hMJpsiInclusive = createHistoFromPdf(pdfJpsiInclusive, yieldJpsiInclusive->getVal(), "hMJpsiInclusive", m, nBinsModel, mMin, mMax, nBinsM);
     TH1* hMPsi2s = createHistoFromPdf(pdfPsi2s, yieldPsi2s->getVal(), "hMPsi2s", m, nBinsModel, mMin, mMax, nBinsM);
     */
    
    TH1* hPtJpsiExclusive = createHistoFromPdf(pdfJpsiExclusive, yieldJpsiExclusive->getVal(), "hPtJpsiExclusive", pt, nBinsModel, ptMin, ptMax, nBinsPt);
    TH1* hPtJpsiDissociative = nullptr;
    if (!exclusiveOnly) hPtJpsiDissociative = createHistoFromPdf(pdfJpsiDissociative, yieldJpsiDissociative->getVal(), "hPtJpsiDissociative", pt, nBinsModel, ptMin, ptMax, nBinsPt);
    TH1* hPtJpsiGammaPb = createHistoFromPdf(pdfJpsiGammaPb, yieldJpsiGammaPb->getVal(), "hPtJpsiGammaPb", pt, nBinsHistPdf, ptMin, ptMax, nBinsPt);
    TH1* hPtTwoGamma = createHistoFromPdf(pdfTwoGamma, yieldTwoGamma->getVal(), "hPtTwoGamma", pt, nBinsModel, ptMin, ptMax, nBinsPt);
    TH1* hPtBackground = createHistoFromPdf(pdfBackground, yieldBkg->getVal(), "hPtBackground", pt, nBinsModel, ptMin, ptMax, nBinsPt);
    /*
     TH1* hPtJpsiInclusive = createHistoFromPdf(pdfJpsiInclusive, yieldJpsiInclusive->getVal(), "hPtJpsiInclusive", pt, nBinsModel, ptMin, ptMax, nBinsPt);
     TH1* hPtPsi2s = createHistoFromPdf(pdfPsi2s, yieldPsi2s->getVal(), "hPtPsi2s", pt, nBinsModel, ptMin, ptMax, nBinsPt);
     */
    
    /*
    TH1* hPtModel = createHistoFromPdf(pdfJpsiExclusive, yieldJpsiExclusive->getVal(), "hPtModel", pt, nBinsModel, ptMin, ptMax, nBinsPt);
    hPtModel->Add(hPtJpsiDissociative);
    hPtModel->Add(hPtJpsiGammaPb);
    hPtModel->Add(hPtTwoGamma);
    hPtModel->Add(hPtBackground);
     */
    
    // Write them
    TFile* fMass = new TFile(Form("mass%s.root", suffix.c_str()), "recreate");
    hMData->Write();
    hMModel->Write();
    hMJpsiExclusive->Write();
    if (!exclusiveOnly) hMJpsiDissociative->Write();
    hMJpsiGammaPb->Write();
    hMTwoGamma->Write();
    hMBackground->Write();
    /*
     hMJpsiInclusive->Write();
     hMPsi2s->Write();
     */
    fMass->Close();
    
    
    TFile* fPt = new TFile(Form("pt%s.root", suffix.c_str()), "recreate");
    hPtData->Write();
    hPtModel->Write();
    hPtJpsiExclusive->Write();
    if (!exclusiveOnly) hPtJpsiDissociative->Write();
    hPtJpsiGammaPb->Write();
    hPtTwoGamma->Write();
    hPtBackground->Write();
    /*
     hPtJpsiInclusive->Write();
     hPtPsi2s->Write();
     */
    
    fPt->Close();
    
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
    
    /*
    gROOT->ProcessLine(".L Include/ExtendedCrystalBall.cxx+") ;
    gSystem->Load("./Include/ExtendedCrystalBall_cxx.so") ;
     */
    
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
        
        if (noInclusive) std::cout << "No inclusive contribution" << std::endl;
        else std::cout << "There is an inclusive contribution" << std::endl;
        if (!exclusiveOnly) {
            WriteResults(wspace, period, rapMin, rapMax, exp, chi2, noInclusive);
        }
        ExportHist(wspace, exclusiveOnly);
    }
    
}
