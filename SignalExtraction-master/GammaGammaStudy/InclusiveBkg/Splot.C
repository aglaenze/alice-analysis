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
#include "RooLandau.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooChi2Var.h"

#include "../../Include/Initiate.C"
#include "../../Include/Utils.C"

using namespace RooFit;
using namespace std;

void AddModel(RooWorkspace* ws, string period) {
    // Define model
    
    //get what we need of the workspace
    RooRealVar* pt = ws->var("fTrkTrkPt");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    
    RooRealVar* yieldTwoGamma = new RooRealVar("yieldTwoGamma","yieldTwoGamma", 200, 0., 2.e3);
    RooRealVar* yieldBkg = new RooRealVar("yieldBkg","yieldBkg", 500, 0., 2.e3);
    
    // Landau
    double muValue = 0, sigmaValue = 0;
    double pt0Value = 0, nIncValue = 0, bIncValue = 0;
    //if (! GetParametersLandau(period, muValue, sigmaValue, pt0Value, nIncValue) ) {cout << "Parameters not loaded" << endl; return;}
    if (! GetParametersLandau(period, muValue, sigmaValue, bIncValue, nIncValue) ) {cout << "Parameters not loaded" << endl; return;}
    
    RooRealVar* mean = new RooRealVar("mean_gg","mean_gg", muValue);
    RooRealVar* sigma = new RooRealVar("sigma_gg","sigma_gg", sigmaValue); //3.43,3.73);
    
    RooLandau* ptTwoGamma = new RooLandau("ptTwoGamma","ptTwoGamma", *pt, *mean, *sigma);
    
    //RooRealVar *pt0 = new RooRealVar("pt0","pt0", 0.4, 0.2, 3);
    RooRealVar *bInc = new RooRealVar("bInc","bInc", 5, 0, 8);
    RooRealVar *nInc = new RooRealVar("nInc","nInc", 3, 1, 8);
    RooGenericPdf* ptBkg = new RooGenericPdf("ptBackground","ptBackground","(2*fTrkTrkPt*(1.+(fTrkTrkPt**2)*(bInc/nInc))**(-nInc))", RooArgSet(*pt, *bInc, *nInc));
    
    
    RooAbsPdf* model = new RooAddPdf("fit", "fit", RooArgList(*ptTwoGamma, *ptBkg), RooArgList(*yieldTwoGamma, *yieldBkg), kFALSE);
    
    ws->import(*model);
}

void DoSPlot(RooWorkspace* ws) {
    
    RooDataSet* data = (RooDataSet*) ws->data("data");
    RooRealVar* pt = ws->var("fTrkTrkPt");
    RooRealVar* m = ws->var("fTrkTrkM");
    
    RooAbsPdf* model = ws->pdf("fit");
    
    RooRealVar* yieldTwoGamma = ws->var("yieldTwoGamma");
    RooRealVar* yieldBkg = ws->var("yieldBkg");
    
    model->fitTo(*data, Extended(), Minos(true), Strategy(1));
    //model->fitTo(*data, Extended(), Minos(true), Strategy(2), SumW2Error(true));
    //The sPlot technique requires that we fix the parameters
    // of the model that are not yields after doing the fit
    RooRealVar* mean_gg = ws->var("mean_gg");
    RooRealVar* sigma_gg = ws->var("sigma_gg");
    RooRealVar* bInc = ws->var("bInc");
    RooRealVar* nInc = ws->var("nInc");
    
    mean_gg->setConstant();
    sigma_gg->setConstant();
    bInc->setConstant();
    nInc->setConstant();
    
    RooMsgService::instance().setSilentMode(true);
    
    //Now we use the SPlot class to add SWeight to our data set
    // based on our model and our yield variables
    
    RooStats::SPlot * sData = new RooStats::SPlot("sData","splot", *data, model, RooArgList(*yieldTwoGamma, *yieldBkg) );
    
    //Check Sweight properties
    std::cout << "Check SWeights: " << std::endl;
    
    std::cout << std::endl << "Yield of Gamma-Gamma is "
    << yieldTwoGamma->getVal() << ". From sWeights it is "
    << sData->GetYieldFromSWeight("yieldTwoGamma") << std::endl;
    
    
    std::cout << std::endl << "Yield of bkg is "
    << yieldBkg->getVal() << ". From sWeights it is "
    << sData->GetYieldFromSWeight("yieldBkg") << std::endl;
    
    for (Int_t i=0; i < 10; i ++ )
    {
        std::cout << "Gamma-Gamma Weight "<< sData->GetSWeight(i, "yieldTwoGamma")
        << " bkg Yield "<< sData->GetSWeight(i, "yieldBkg")
        << " Total Weight "<< sData->GetSumOfEventSWeight(i)
        << std::endl;
    }
    
    //import the new data set with Sweight
    
    std::cout << "import new dataset with sWeight" << std::endl;
    ws->import(*data, Rename("dataWithSWeights"));
}


//____________________________________
void MakePlots(RooWorkspace *ws, string rootfilePath, string rootfilePathMc, string period, bool exp, bool exclusiveOnly, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax){
    
    
    int ptBinNumber = int(10*(ptMax-ptMin));
    
    
    //get what we need of the workspace
    RooAbsPdf* model = ws->pdf("fit");
    RooAbsPdf* ptTwoGamma = ws->pdf("ptTwoGamma");
    RooAbsPdf* ptBackground = ws->pdf("ptBackground");
    
    RooRealVar* m = ws->var("fTrkTrkM");
    RooRealVar* pt = ws->var("fTrkTrkPt");
    
    RooDataSet* sdata = (RooDataSet*) ws->data("dataWithSWeights");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    
    RooRealVar* yieldTwoGamma = ws->var("yieldTwoGamma");
    RooRealVar* yieldBkg = ws->var("yieldBkg");
    
    
    vector <RooPlot*> frames = {};
    for (int i =0; i<2; i++) {
        frames.push_back(pt->frame());
        sdata->plotOn(frames[i], Binning(ptBinNumber));
        //h_sdata->plotOn(frame, Name("h_sdata"), LineColor(kWhite));
        model->plotOn(frames[i], LineWidth(2));
        model->plotOn(frames[i], Components(*model), LineStyle(kDashed), LineColor(kBlue), Name("sum"), LineWidth(2));
        model->plotOn(frames[i], Components(*ptTwoGamma), LineStyle(kDashed), LineColor(kRed), Name("ptTwoGamma"), LineWidth(2));
        model->plotOn(frames[i], Components(*ptBackground), LineStyle(kDashed), LineColor(kGreen), Name("ptBackground"), LineWidth(2));
        frames[i]->SetTitle("Fit to model to discriminating variable");
    }
    frames[1]->SetTitle("Fit to model to discriminating variable (zoom)");
    TLegend* lgd1 = new TLegend(0.71, 0.7, 0.9, 0.9);
    lgd1->AddEntry(frames[0]->findObject("ptTwoGamma"), "#gamma#gamma", "L");
    lgd1->AddEntry(frames[0]->findObject("ptBackground"), "Bkg", "L");
    lgd1->AddEntry(frames[0]->findObject("sum"), "Sum", "L");
    
    double yMax0 = frames[0]->GetMaximum();
    double xPos = 0.5;
    TLatex* txt0 = new TLatex(xPos,0.8*yMax0,Form("#bf{N_{#gamma#gamma} = %.0f #pm %.0f}", yieldTwoGamma->getVal(), yieldTwoGamma->getError()));
    TLatex* txt1 = new TLatex(xPos,0.7*yMax0,Form("#bf{N_{bkg} = %.0f #pm %.0f}", yieldBkg->getVal(), yieldBkg->getError()));
    frames[0]->addObject(txt0) ;
    frames[0]->addObject(txt1) ;
    
    frames[1]->GetXaxis()->SetRangeUser(0,1);
    
    // draw plots
    TCanvas* cv = new TCanvas("splot","splot", 600, 900) ;
    cv->Divide(2,3);
    // Draw pt fit
    cv->cd(1);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    frames[0]->Draw();
    lgd1->Draw();
    cv->cd(2);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    frames[1]->Draw();
    lgd1->Draw();
    
    
    // Draw quality plots
    cv->cd(3);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    sdata->plotOn(frames[0]);
    model->plotOn(frames[0]);
    RooHist* hpull = frames[0]->pullHist();
    hpull->SetTitle("(data - fit)/#sigma of the p_{T} fit");
    hpull->Draw("");
    
    cv->cd(4);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    sdata->plotOn(frames[1]);
    model->plotOn(frames[1]);
    RooHist* hpull2 = frames[1]->pullHist();
    hpull2->GetXaxis()->SetRangeUser(0,1);
    hpull2->SetTitle("(data - fit)/#sigma of the p_{T} fit (zoom)");
    hpull2->Draw("");
    
    cv->cd(5);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    RooHist* hresid = frames[0]->residHist();
    hresid->SetTitle("Residuals (data - fit) of the p_{T} fit");
    hresid->Draw("");
    
    cv->cd(6);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    RooHist* hresid2 = frames[1]->residHist();
    hresid2->GetXaxis()->SetRangeUser(0,1);
    hresid2->SetTitle("Residuals (data - fit) of the p_{T} fit (zoom)");
    hresid2->Draw("");
    
    // Save plots
    cv->SaveAs(Form("Plots/%s/splot-discrimation-%.1f-%.1f-%.1f-%.1f.pdf", period.c_str(), mMin, mMax, abs(rapMax), abs(rapMin)));
    
    // The SPlot class adds a new variable that has the name of the corresponding
    // yield + "_sw".
    // create weighted data set for JPsi
    RooDataSet * dataw_gg = new RooDataSet(sdata->GetName(),sdata->GetTitle(),sdata,*sdata->get(), 0, "yieldTwoGamma_sw") ;
    // create weighted data set for Background
    RooDataSet * dataw_bkg = new RooDataSet(sdata->GetName(),sdata->GetTitle(),sdata,*sdata->get(), 0, "yieldBkg_sw") ;
    
    // scattered plots of gamma-gamma and bkg
    TCanvas* cv2 = new TCanvas("cv2", "cv2", 600, 300);
    cv2->Divide(2);
    // First gamma-gamma
    cv2->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    
    TH1* histGg = dataw_gg->createHistogram("fTrkTrkM,fTrkTrkPt", 10, ptBinNumber);
    histGg->SetTitle("#gamma#gamma");
    histGg->GetXaxis()->SetMaxDigits(2);
    histGg->Draw("colz");
    
    // Second bkg
    cv2->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    TH1* histBkg = dataw_bkg->createHistogram("fTrkTrkM,fTrkTrkPt", 10, ptBinNumber);
    histBkg->SetTitle("Inclusive continuum");
    histBkg->GetXaxis()->SetMaxDigits(2);
    histBkg->Draw("colz");
    
    // Save plots
    cv2->SaveAs(Form("Plots/%s/splot-scattered-%.1f-%.1f-%.1f-%.1f.pdf", period.c_str(), mMin, mMax, abs(rapMax), abs(rapMin)));
    
}

//____________________________________
void DrawWeights(RooWorkspace *ws, string period, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax){
    
    //get what we need of the workspace
    RooRealVar* m = ws->var("fTrkTrkM");
    RooRealVar* pt = ws->var("fTrkTrkPt");
    
    RooDataSet* sData = (RooDataSet*) ws->data("dataWithSWeights");
    
    //std::cout << std::endl << std::endl << std::endl << std::endl;
    const int nEntries = (int) sData->sumEntries();
    Double_t massForWeights[nEntries], ptForWeights[nEntries], ggWeight[nEntries], bkgWeight[nEntries];
    
    for (Int_t i=0; i < (Int_t)nEntries; i ++ ) {
        massForWeights[i] = sData->get(i)->getRealValue("fTrkTrkM");
        ptForWeights[i] = sData->get(i)->getRealValue("fTrkTrkPt");
        
        ggWeight[i] = sData->get(i)->getRealValue("yieldTwoGamma_sw");
        bkgWeight[i] = sData->get(i)->getRealValue("yieldBkg_sw");
        
    }
    //std::cout << std::endl << std::endl << std::endl << std::endl;
    
    // TGraphs of sWeights = f(m) or sWeights = f(pt)
    TGraph* gWeightsMassGg = new TGraph(nEntries, massForWeights, ggWeight);
    TGraph* gWeightsPtGg = new TGraph(nEntries, ptForWeights, ggWeight);
    TGraph* gWeightsMassBkg = new TGraph(nEntries, massForWeights, bkgWeight);
    TGraph* gWeightsPtBkg = new TGraph(nEntries, ptForWeights, bkgWeight);
    
    gWeightsMassGg->SetTitle("#gamma#gamma  sWeights = f(m)");
    gWeightsPtGg->SetTitle("#gamma#gamma  sWeights = f(pt)");
    gWeightsMassBkg->SetTitle("Background sWeights = f(m)");
    gWeightsPtBkg->SetTitle("Background sWeights = f(pt)");
    
    gWeightsMassGg->GetXaxis()->SetTitle("m");
    gWeightsPtGg->GetXaxis()->SetTitle("pt");
    gWeightsMassBkg->GetXaxis()->SetTitle("m");
    gWeightsPtBkg->GetXaxis()->SetTitle("pt");
    
    gWeightsMassGg->GetYaxis()->SetTitle("#gamma#gamma sWeights");
    gWeightsPtGg->GetYaxis()->SetTitle("#gamma#gamma sWeights");
    gWeightsMassBkg->GetYaxis()->SetTitle("Bkg sWeights");
    gWeightsPtBkg->GetYaxis()->SetTitle("Bkg sWeights");
    
    TCanvas* cv = new TCanvas("sweights","sweights", 600, 900) ;
    cv->Divide(2,3);
    
    cv->cd(1);
    gPad->SetLeftMargin(0.2);
    gWeightsMassGg->Draw("A*");
    
    cv->cd(2);
    gPad->SetLeftMargin(0.2);
    gWeightsPtGg->Draw("A*");
    
    cv->cd(3);
    gPad->SetLeftMargin(0.2);
    gWeightsMassBkg->Draw("A*");
    
    cv->cd(4);
    gPad->SetLeftMargin(0.2);
    gWeightsPtBkg->Draw("A*");
    
    // 2d histograms
    TH2F* hWeightsBkg2d = new TH2F("hWeightsBkg2d", "Background sWeights)", 50, mMin, mMax, 100, ptMin, ptMax);
    TH2F* hWeightsGg2d = new TH2F("hWeightsGg2d", "#gamma#gamma sWeights", 50, mMin, mMax, 100, ptMin, ptMax);
    for (int i = 0; i<nEntries; i++) {
        hWeightsBkg2d->Fill(massForWeights[i], ptForWeights[i], bkgWeight[i]);
        hWeightsGg2d->Fill(massForWeights[i], ptForWeights[i], ggWeight[i]);
    }
    cv->cd(5);
    gPad->SetLeftMargin(0.2);
    hWeightsBkg2d->GetXaxis()->SetTitle("Mass");
    hWeightsBkg2d->GetYaxis()->SetTitle("Pt");
    //hWeightsBkg2d->Draw("CONTZ");
    hWeightsBkg2d->Draw("COLZ");
    
    cv->cd(6);
    gPad->SetLeftMargin(0.2);
    hWeightsGg2d->GetXaxis()->SetTitle("Mass");
    hWeightsGg2d->GetYaxis()->SetTitle("Pt");
    hWeightsGg2d->Draw("COLZ");
    
    // Save plots
    cv->SaveAs(Form("Plots/%s/splot-weights-%.1f-%.1f.pdf", period.c_str(), mMin, mMax));
}

void Splot(string rootfilePath = "", string rootfilePathMc = "", std::vector<string> periods = {"LHC16r", "LHC16s"}) {
    
    Double_t mMin = 2., mMax = 3.5, ptMin = 0., ptMax = 3.5;
    Double_t rapMin = -4., rapMax = -2.5;
    bool useCuts = true, exp = false, exclusiveOnly = false;
    bool logScale = false, drawPulls = false;
    
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
        if ( ! Initiate(period, mMin, mMax, ptMin, ptMax, rapMin, rapMax, useCuts, logScale, drawPulls, exp, exclusiveOnly)) {cout << "Something wrong at initialisation"; return;}
        
        if (exclusiveOnly) {
            ptMax = 1.2;
        }
        /*
         else {
         if (period == "LHC16s") ptMax = 3;
         }
         */
        
        std::list<TCut> mCutList = DefineCuts(period, exclusiveOnly);
        // Define cuts
        TCut mCut = "";
        
        for (std::list <TCut>::iterator iter = mCutList.begin(); iter != mCutList.end(); ++iter) {mCut += *iter;}
        if (!useCuts) mCut = "";
        
        // Open the file
        TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()), "READ");
        
        // Connect to the tree
        TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
        //Create a new workspace to manage the project
        RooWorkspace* wspace = new RooWorkspace("myJpsi");
        ImportDataSet(wspace, fAnaTree, mCut, mMin, mMax, ptMin, ptMax, rapMin, rapMax);
        AddModel(wspace, period);
        wspace->Print();
        
        DoSPlot(wspace);
        MakePlots(wspace, rootfilePath, rootfilePathMc, period, exp, exclusiveOnly, useCuts, mMin, mMax, ptMin, ptMax, rapMin, rapMax);
        DrawWeights(wspace, period, mMin, mMax, ptMin, ptMax);
        
        //cleanup
        delete wspace;
    }
    
    
}


