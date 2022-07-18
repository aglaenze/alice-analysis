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

#include "Include/FitUtils.C"

using namespace RooFit;
using namespace std;

//Double_t mLimitPsi2s = 3.65;
Double_t mLimitPsi2s = 10;


void SubtractBkg(TH1* hPtMc , TH1* hBackground, TH1* hSub) {
    
    int nBins = hBackground->GetNbinsX();
    
    double scale = (double)hBackground->GetMaximum()/hPtMc->GetMaximum();
    //hPtMc->Scale(scale);
    
    // Create the subtraction histo between packground data and MC (two gamma)
    //TH1F* hSub = new TH1F("histSub", "hist sub", nBins, xMin, xMax);
    for (int l = 0; l<nBins; l++) {
        int binBkg = hBackground->GetBinContent(l);
        int binMc = hPtMc->GetBinContent(l);
        if (binBkg > binMc) hSub->SetBinContent(l, binBkg - binMc);
        else hSub->SetBinContent(l, 0);
    }
}

void SmoothHistLog(TH1* hPtBackground, TH1* hPtBkgSmooth, int avgSize = 3) {
    double negValue = -10;    // in case the value is negative, pb log(neg) so set log to this negvalue
    int ptBinNumber = hPtBkgSmooth->GetNbinsX();
    for (int k = 0; k<ptBinNumber; k++) {
        int binNum = k+1;
        double smoothVal = 0;
        int nSmoothBins = 1;
        if (hPtBackground->GetBinContent(binNum) > 0) {smoothVal = TMath::Log(hPtBackground->GetBinContent(binNum));}
        else smoothVal = -10;    // by default say there's 10^-2 background events (originally) in this bin
        double size = binNum*0.8;    // number of bins to smooth
        if (size > avgSize) size = avgSize;
        for (int j = 1; j < (int)size; j++ ) {
            if (binNum-j > 0) {
                nSmoothBins++;
                double neighbourVal = hPtBackground->GetBinContent(binNum-j);
                if (neighbourVal>0) {smoothVal+= TMath::Log(neighbourVal);}
                else smoothVal+= negValue;
            }
            if (binNum+j < ptBinNumber) {
                nSmoothBins++;
                double neighbourVal = hPtBackground->GetBinContent(binNum+j);
                if (neighbourVal>0) {smoothVal+= TMath::Log(hPtBackground->GetBinContent(binNum+j));}
                else smoothVal+= negValue;
            }
        }
        smoothVal /= (double)nSmoothBins;
        smoothVal = TMath::Exp(smoothVal);
        Double_t x = hPtBackground->GetXaxis()->GetBinCenter(binNum);
        hPtBkgSmooth->Fill(x, smoothVal);
        //std::cout << x << " " << smoothVal << std::endl;
    }
    hPtBkgSmooth->SetOption("hist");
}

void SmoothHist(TH1* hPtBackground, TH1* hPtBkgSmooth, int avgSize = 3) {
    int ptBinNumber = hPtBkgSmooth->GetNbinsX();
    for (int k = 0; k<ptBinNumber; k++) {
        int binNum = k+1;
        double smoothVal = hPtBackground->GetBinContent(binNum);
        int nSmoothBins = 1;
        double size = binNum*0.8;    // number of bins to smooth
        if (size > avgSize) size = avgSize;
        for (int j = 1; j < (int)size; j++ ) {
            if (binNum-j > 0) {
                nSmoothBins++;
                smoothVal+= hPtBackground->GetBinContent(binNum-j);
            }
            if (binNum+j < ptBinNumber) {
                nSmoothBins++;
                smoothVal+= hPtBackground->GetBinContent(binNum+j);
            }
        }
        smoothVal /= (double)nSmoothBins;
        Double_t x = hPtBackground->GetXaxis()->GetBinCenter(binNum);
        hPtBkgSmooth->Fill(x, smoothVal);
        //std::cout << x << " " << smoothVal << std::endl;
    }
    hPtBkgSmooth->Rebin(2);
    hPtBkgSmooth->Scale(1./2);
    hPtBkgSmooth->SetOption("hist");
}

void AddModel(RooWorkspace* ws, string period, Double_t mMin, Double_t mMax) {
    // Define model
    
    //RooRealVar m("fTrkTrkM","M_{#mu#mu} (GeV/c2)",2.5,3.5);
    RooRealVar m = *ws->var("fTrkTrkM");
    
    LoadMassFitFunctions(ws, period);
    // Now collect m pdf
    RooAbsPdf* jpsi = ws->pdf("jpsi");
    RooAbsPdf* psi2s = ws->pdf("psi2s");
    RooAbsPdf* bkg = ws->pdf("bkg");
    
    //RooRealVar fsig("fsig","signalPhi",0.1,0.,1.);
    RooRealVar fsigJpsi("fsigJpsi","signalJPsi",1000,0.,1.e4);
    RooRealVar fsigPsi2s("fsigPsi2s","signalPsi2s",100,0.,1.e3);
    RooRealVar fbkg("fbkg","fbkg",500,0.,1.e7);
    
    if (mMin > 3.2 || mMax < 3) {fsigJpsi.setVal(0); fsigJpsi.setConstant();}
    
    RooAbsPdf* model;
    if (mMax > mLimitPsi2s) model = new RooAddPdf("mfit", "mfit", RooArgList(*jpsi, *psi2s, *bkg), RooArgList(fsigJpsi, fsigPsi2s, fbkg), kFALSE);
    else model = new RooAddPdf("mfit", "mfit", RooArgList(*jpsi, *bkg), RooArgList(fsigJpsi, fbkg), kFALSE);
    
    ws->import(*model);
}

void AddPtJpsiModel(RooWorkspace* ws, string rootfilePathMC, string period, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax, bool exp, bool exclusiveOnly) {
    // Define fitting model for pt for J/Psi signal only
    // Dissociative contribution and exclusive contribution
    
    RooRealVar pt = *ws->var("fTrkTrkPt");
    LoadJpsiPtFitFunctions(ws, rootfilePathMC, period, ptMin, ptMax, exp, exclusiveOnly);
    
    // Then pt pdf
    RooAbsPdf* ptJpsiExclusive = ws->pdf("ptJpsiExclusive");
    RooAbsPdf* ptJpsiDissociative = ws->pdf("ptJpsiDissociative");
    RooAbsPdf* ptJpsiGammaPb = ws->pdf("ptJpsiGammaPb");
    RooAbsPdf* ptJpsiInclusive = ws->pdf("ptJpsiInclusive");
    
    // Load yields
    LoadYields(ws, period, rapMin, rapMax, mMin, mMax, ptMin, ptMax, exclusiveOnly);
    RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
    RooRealVar* yieldJpsiDissociative = ws->var("yieldJpsiDissociative");
    RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
    RooRealVar* yieldJpsiInclusive = ws->var("yieldJpsiInclusive");
    
    //RooAbsPdf* model = new RooAddPdf("ptfit", "ptfit", RooArgList(*ptJpsiExclusive, *ptJpsiDissociative, *ptJpsiGammaPb, *ptJpsiInclusive), RooArgList(*yieldJpsiExclusive, *yieldJpsiDissociative, *yieldJpsiGammaPb , *yieldJpsiInclusive), kFALSE);
    RooAbsPdf* model = new RooAddPdf("ptfit", "ptfit", RooArgList(*ptJpsiExclusive, *ptJpsiGammaPb), RooArgList(*yieldJpsiExclusive, *yieldJpsiGammaPb), kFALSE);
    
    ws->import(*model);
}

void DoSPlot(RooWorkspace* ws, Double_t mMax) {
    
    RooAbsPdf* model = ws->pdf("mfit");
    RooRealVar* jPsiYield = ws->var("fsigJpsi");
    RooRealVar* psi2sYield = ws->var("fsigPsi2s");
    RooRealVar* bkgYield = ws->var("fbkg");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    
    //model->fitTo(*data, Extended(), Minos(true), Strategy(2), SumW2Error(true));
    model->fitTo(*data, Extended(), Minos(true), Strategy(2));

    //The sPlot technique requires that we fix the parameters
    // of the model that are not yields after doing the fit
    // J/Psi
    RooRealVar* mean_jpsi = ws->var("mean_jpsi");
    RooRealVar* sigma_jpsi = ws->var("sigma_jpsi");
    //RooRealVar* alpha_jpsi = ws->var("alpha_jpsi");
    //RooRealVar* n_jpsi = ws->var("n_jpsi");
    RooRealVar* alpha_jpsi_L = ws->var("alpha_jpsi_L");
    RooRealVar* n_jpsi_L = ws->var("n_jpsi_L");
    RooRealVar* alpha_jpsi_R = ws->var("alpha_jpsi_R");
    RooRealVar* n_jpsi_R = ws->var("n_jpsi_R");
    
    // Psi(2s)
    RooRealVar *mean_psi2s, *sigma_psi2s, *alpha_psi2s_L, *n_psi2s_L, *alpha_psi2s_R, *n_psi2s_R;
    if (mMax > mLimitPsi2s) {
        mean_psi2s = ws->var("mean_psi2s");
        sigma_psi2s = ws->var("sigma_psi2s");
        alpha_psi2s_L = ws->var("alpha_psi2s_L");
        n_psi2s_L = ws->var("n_psi2s_L");
        alpha_psi2s_R = ws->var("alpha_psi2s_R");
        n_psi2s_R = ws->var("n_psi2s_R");
        alpha_psi2s_L->setConstant();
        n_psi2s_L->setConstant();
        alpha_psi2s_R->setConstant();
        n_psi2s_R->setConstant();
    }
    
    // Background
    RooRealVar* a1 = ws->var("a1");
    
    mean_jpsi->setConstant();
    sigma_jpsi->setConstant();
    alpha_jpsi_L->setConstant();
    n_jpsi_L->setConstant();
    alpha_jpsi_R->setConstant();
    n_jpsi_R->setConstant();
    
    a1->setConstant();
    
    RooMsgService::instance().setSilentMode(true);
    
    //Now we use the SPlot class to add SWeight to our data set
    // based on our model and our yield variables
    
    RooStats::SPlot * sData;
    if (mMax > mLimitPsi2s) sData = new RooStats::SPlot("sData","splot", *data, model, RooArgList(*jPsiYield, *psi2sYield, *bkgYield) );
    else sData = new RooStats::SPlot("sData","splot", *data, model, RooArgList(*jPsiYield, *bkgYield) );
    
    //Check Sweight properties
    std::cout << "Check SWeights: " << std::endl;
    
    std::cout << std::endl << "Yield of JPsi is "
    << jPsiYield->getVal() << ". From sWeights it is "
    << sData->GetYieldFromSWeight("fsigJpsi") << std::endl;
    
    
    std::cout << std::endl << "Yield of bkg is "
    << bkgYield->getVal() << ". From sWeights it is "
    << sData->GetYieldFromSWeight("fbkg") << std::endl;
    
    if (mMax > mLimitPsi2s) {
        std::cout << std::endl << "Yield of Psi(2s) is "
        << psi2sYield->getVal() << ". From sWeights it is "
        << sData->GetYieldFromSWeight("fsigPsi2s") << std::endl;
        for (Int_t i=0; i < 10; i ++ )
        {
            std::cout << "JPsi Weight "<< sData->GetSWeight(i, "fsigJpsi")
            << " psi(2s) Yield "<< sData->GetSWeight(i, "fsigPsi2s")
            << " bkg Yield "<< sData->GetSWeight(i, "fbkg")
            << " Total Weight "<< sData->GetSumOfEventSWeight(i)
            << std::endl;
        }
    }
    else {
        for (Int_t i=0; i < 10; i ++ )
        {
            std::cout << "JPsi Weight "<< sData->GetSWeight(i, "fsigJpsi")
            << " bkg Yield "<< sData->GetSWeight(i, "fbkg")
            << " Total Weight "<< sData->GetSumOfEventSWeight(i)
            << std::endl;
        }
    }
    
    //import the new data set with Sweight
    
    std::cout << "import new dataset with sWeight" << std::endl;
    ws->import(*data, Rename("dataWithSWeights"));
    
    std::cout << "Sigma = " << sigma_jpsi->getVal() << "\n\n\n\n\n" << std::endl;
}


//____________________________________
void MakePlots(RooWorkspace *ws, string rootfilePath, string rootfilePathMc, string period, bool exp, bool exclusiveOnly, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax){
    
    
    int ptBinNumber = int(10*(ptMax-ptMin));
    
    //make some plots
    TCanvas* cv = new TCanvas("splot","splot", 800, 800) ;
    cv->Divide(3,3);
    
    //get what we need of the workspace
    RooAbsPdf* model = ws->pdf("mfit");
    RooAbsPdf* jPsiModel = ws->pdf("jpsi");
    RooAbsPdf* bkgModel = ws->pdf("bkg");
    
    RooRealVar* m = ws->var("fTrkTrkM");
    RooRealVar* pt = ws->var("fTrkTrkPt");
    
    RooDataSet* sdata = (RooDataSet*) ws->data("dataWithSWeights");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    
    RooRealVar* jPsiYield = ws->var("fsigJpsi");
    RooRealVar* psi2sYield = ws->var("fsigPsi2s");
    RooRealVar* bkgYield = ws->var("fbkg");
    
    int nBins = 20;
    TH1* hist_sdata = sdata->createHistogram("fTrkTrkM", nBins);
    RooDataHist* h_sdata = new RooDataHist("h_sdata", "h_sdata", RooArgSet(*m), hist_sdata);
    
    // Draw mass fit
    cv->cd(1);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    RooPlot* frame = m->frame();
    sdata->plotOn(frame, Binning(50));
    //h_sdata->plotOn(frame, Name("h_sdata"), LineColor(kWhite));
    model->plotOn(frame, LineWidth(2));
    model->plotOn(frame, Components(*model), LineStyle(kDashed), LineColor(kBlue), Name("sum"), LineWidth(2));
    model->plotOn(frame, Components(*jPsiModel), LineStyle(kDashed), LineColor(kRed), Name("jpsi"), LineWidth(2));
    model->plotOn(frame, Components(*bkgModel), LineStyle(kDashed), LineColor(kGreen), Name("bkg"), LineWidth(2));
    frame->SetTitle("Fit to model to discriminating variable");
    double yMax0 = frame->GetMaximum();
    double xPos = mMin+(mMax-mMin)*0.06;
    //if (mMax < mLimitPsi2s) xPos = mMin+(mMax-mMin)/3;
    TLatex* txt0 = new TLatex(xPos,0.7*yMax0,Form("#bf{N_{J/#Psi} = %.0f #pm %.0f}", jPsiYield->getVal(), jPsiYield->getError()));
    TLatex* txt1 = new TLatex(xPos,0.6*yMax0,Form("#bf{N_{bkg} = %.0f #pm %.0f}", bkgYield->getVal(), bkgYield->getError()));
    frame->addObject(txt0) ;
    frame->addObject(txt1) ;
    if (mMax > mLimitPsi2s) {
        RooAbsPdf* psi2sModel = ws->pdf("psi2s");
        model->plotOn(frame, Components(*psi2sModel), LineStyle(kDashed), LineColor(kOrange));
        TText* txtPsi2s = new TText((mMin+mMax+0.1)/2,0.45*yMax0,Form("#bf{%.1f #pm %.1f Psi(2s)}", psi2sYield->getVal(), psi2sYield->getError()));
        frame->addObject(txtPsi2s) ;
    }
    int nDofM = model->getParameters(data)->selectByAttrib("Constant",kFALSE)->getSize();
    //cout << "nDOF = " << nDofM << endl;
    //return;
    //RooDataHist* h_dataWithSWeights = (RooDataHist*)frame->findObject("h_dataWithSWeights");
    //TH1* h_dataWithSWeights = (TH1*)frame->findObject("h_dataWithSWeights");
    /*
     TCanvas* cv0 = new TCanvas();
     //h_dataWithSWeights->Draw();
     h_data->Draw();
     cv0->SaveAs("Plots/test2.pdf");
     */
    //return;
    //Int_t ndf = h_dataWithSWeights->numEntries()-nDofM;
    Int_t ndf = nBins-nDofM;
    cout << " nDofM = " << nDofM << endl;
    cout << " ndf = " << ndf << endl;
    //cout << " entries in hist = " << h_dataWithSWeights->GetEntries() << endl;
    
    //double chi2M = frame->chiSquare("sum", "h_sdata", nDofM);
    RooChi2Var* chi2Var = new RooChi2Var("chi2","chi2", *model, *h_sdata, DataError(RooAbsData::Poisson));
    double chi2M = chi2Var->getVal();
    TLatex* txtChi2 = new TLatex(xPos, yMax0,Form("#bf{#chi^{2}/ndf = %.1f/%d = %.2f}", chi2M, ndf, chi2M/ndf));
    //frame->addObject(txtChi2);
    frame->SetMaximum(1.1*frame->GetMaximum());
    frame->Draw();
    cout << "chi2M = " << chi2M << endl;
    //return;
    
    TLegend* lgd1 = new TLegend(0.71, 0.7, 0.9, 0.9);
    lgd1->AddEntry(frame->findObject("jpsi"), "J/#Psi", "L");
    lgd1->AddEntry(frame->findObject("bkg"), "Bkg", "L");
    lgd1->AddEntry(frame->findObject("sum"), "Sum", "L");
    lgd1->Draw();
    
    // write y and pt info
    TLegend* lgdInfo = new TLegend(0.15, 0.73, 0.5, 0.9);
    lgdInfo->SetMargin(0.05);
    lgdInfo->SetTextSize(0.04);
    lgdInfo->AddEntry((TObject*)0, Form("%.2f < y < %.2f", -rapMax, -rapMin), "");
    lgdInfo->AddEntry((TObject*)0, Form("%.1f < M < %.1f GeV/c^{2}", mMin, mMax), "");
    lgdInfo->AddEntry((TObject*)0, Form("%.1f < p_{T} < %.1f GeV/c", ptMin, ptMax), "");
    lgdInfo->Draw();
    
    
    // Draw quality plots
    cv->cd(4);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    sdata->plotOn(frame, Binning(ptBinNumber));
    model->plotOn(frame, LineWidth(2));
    RooHist* hpull = frame->pullHist();
    hpull->GetXaxis()->SetRangeUser(mMin, mMax);
    hpull->SetTitle("(data - fit)/#sigma of the mass fit");
    hpull->Draw("");
    cv->cd(5);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    RooHist* hresid = frame->residHist();
    hresid->GetXaxis()->SetRangeUser(mMin, mMax);
    hresid->SetTitle("Residuals (data - fit) of the mass fit");
    hresid->Draw("");
    
    // The SPlot class adds a new variable that has the name of the corresponding
    // yield + "_sw".
    // create weighted data set for JPsi
    RooDataSet * dataw_jpsi = new RooDataSet(sdata->GetName(),sdata->GetTitle(),sdata,*sdata->get(),0,"fsigJpsi_sw") ;
    //RooDataSet * dataw_jpsi = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),"fTrkTrkM < 3.2 && fTrkTrkM > 2.8","fsigJpsi_sw") ;
    // create weighted data set for Background
    RooDataSet * dataw_bkg = nullptr;
    if (mMax > 3.1 && mMin < 3.1 ) dataw_bkg = new RooDataSet(sdata->GetName(),sdata->GetTitle(),sdata,*sdata->get(),0,"fbkg_sw") ;
    else dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0) ;
    
    
    // Draw J/Psi pt with weights
    cv->cd(2);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    RooPlot* frameJpsi = pt->frame() ;
    int binningNum = 15*int(ptMax-ptMin);
    
    double fitRangeMax = 0.38;
    //double fitRangeMax = ptMax;
    int nBinsPt = 30;
    Double_t newPtMax = 1.2;
    int nBinsPt2 = int( newPtMax/fitRangeMax*nBinsPt )-1;
    
    //dataw_jpsi->plotOn(frameJpsi, Name("jPsiData"), DataError(RooAbsData::SumW2), Binning(binningNum) ) ;
    dataw_jpsi->plotOn(frameJpsi, Name("jPsiData"), DataError(RooAbsData::SumW2), Binning(nBinsPt2) ) ;
    frameJpsi->SetTitle("p_{T} distribution for J/#Psi with weights");
    

    
    /*
     // Write parameters on plot
     RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
     RooRealVar* yieldJpsiDissociative = ws->var("yieldJpsiDissociative");
     RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
     RooRealVar* yieldJpsiInclusive = ws->var("yieldJpsiInclusive");
     TLatex* txt2 = new TLatex(xPosJpsi,0.85*yMax,Form("%.1f #pm %.1f exclusive J/Psi", yieldJpsiExclusive->getVal(), yieldJpsiExclusive->getError()));
     frameJpsi->addObject(txt2) ;
     TLatex* txt = new TLatex(xPosJpsi,0.77*yMax,Form("%.1f #pm %.1f dissociative J/Psi", yieldJpsiDissociative->getVal(), yieldJpsiDissociative->getError()));
     frameJpsi->addObject(txt) ;
     TLatex* txt4 = new TLatex(xPosJpsi,0.69*yMax,Form("%.1f #pm %.1f #gamma-Pb J/Psi", yieldJpsiGammaPb->getVal(), yieldJpsiGammaPb->getError()));
     frameJpsi->addObject(txt4) ;
     TLatex* txt3 = new TLatex(xPosJpsi,0.61*yMax,Form("%.1f #pm %.1f inclusive J/Psi", yieldJpsiInclusive->getVal(), yieldJpsiInclusive->getError()));
     if (!exclusiveOnly && period == "LHC16r") frameJpsi->addObject(txt3) ;
     */
    
    /*
     // Write fit function parameters on the plot
     double xPosJpsi2 = xPosJpsi*1.2;
     RooRealVar* bExc = ws->var("bExc");
     TLatex* txtExc = new TLatex(xPosJpsi2,0.45*yMax,Form("b_{exc} = %.2f #pm %.2f", bExc->getVal(), bExc->getError()));
     frameJpsi->addObject(txtExc) ;
     RooRealVar* bDiss = ws->var("bDiss");
     
     TLatex* txtDiss = new TLatex(xPosJpsi2,0.35*yMax,Form("b_{diss} = %.2f #pm %.2f", bDiss->getVal(), bDiss->getError()));
     frameJpsi->addObject(txtDiss) ;
     
     if (!exp) {
     RooRealVar* nDiss = ws->var("nDiss");
     TLatex* txtDiss2 = new TLatex(xPosJpsi2,0.25*yMax,Form("n_{diss} = %.2f #pm %.2f", nDiss->getVal(), nDiss->getError()));
     frameJpsi->addObject(txtDiss2) ;
     }
     RooRealVar* pt0 = ws->var("pt0");
     RooRealVar* nInc = ws->var("nInc");
     TLatex* txtInc1 = new TLatex(xPosJpsi2,0.15*yMax,Form("p_{0} = %.2f #pm %.2f", pt0->getVal(), pt0->getError()));
     TLatex* txtInc2 = new TLatex(xPosJpsi2,0.05*yMax,Form("n_{inc} = %.2f #pm %.2f", nInc->getVal(), nInc->getError()));
     if (!exclusiveOnly && period == "LHC16r") {frameJpsi->addObject(txtInc1) ; frameJpsi->addObject(txtInc2) ;}
     */
    
    frameJpsi->Draw() ;
    
    // Draw Psi(2s) pt with sweights
    if (mMax > mLimitPsi2s) {
        cv->cd(6);
        gPad->SetLeftMargin(0.15) ;
        gPad->SetBottomMargin(0.15) ;
        // create weighted data set for Psi(2s)
        RooDataSet * dataw_psi2s = new RooDataSet(sdata->GetName(),sdata->GetTitle(),sdata,*sdata->get(),0,"fsigPsi2s_sw") ;
        RooPlot* framePsi2s = pt->frame() ;
        dataw_psi2s->plotOn(framePsi2s, DataError(RooAbsData::SumW2), Binning(50) ) ;
        framePsi2s->SetTitle("p_{T} distribution for #Psi(2s) with weights");
        framePsi2s->Draw() ;
    }
    
    // Draw background pt with sweights
    cv->cd(3);
    gPad->SetLeftMargin(0.15) ;
    gPad->SetBottomMargin(0.15) ;
    RooPlot* frameBkg = pt->frame() ;
    frameBkg->SetTitle("p_{T} distribution for background with weights");
    //frameBkg->GetXaxis()->SetRangeUser(0,3);
    
    TLegend* legend = new TLegend(0.4, 0.6, 0.9, 0.9);
    legend->SetTextSize(0.05);
    legend->SetFillColor(kWhite);
    legend->SetLineColor(kWhite);
    
    
    //dataw_bkg->plotOn(frameBkg, DataError(RooAbsData::SumW2), Binning(binningNum)) ;
    dataw_bkg->plotOn(frameBkg, DataError(RooAbsData::SumW2), Binning(nBinsPt2)) ;
    
    // Draw fit of pt
    // Add MC templates for (gammma gamma to mu mu) for comparison
    /*
     TCut cutMc = Form("fTrkTrkM > %f && fTrkTrkM < %f", mMin, mMax);
     GetPtHistMC(ws, rootfilePathMc, period, "kTwoGammaToMuLow", ptMin, ptMax, cutMc);
     GetV0Template(ws, rootfilePath, period);
     RooAbsPdf* pdfV0 = ws->pdf("ptV0C7");
     */
    // Landau
    RooRealVar muGg("mean_gg","mean_gg", 0.09, 0.05, 0.13); //3.43,3.73);
    RooRealVar sigmaGg("sigma_gg","sigma_gg",0.015, 0.008, 0.10); //3.43,3.73);
    
    if (!exclusiveOnly) {
        muGg.setVal(0.0799);
        sigmaGg.setVal(0.0351);
        muGg.setConstant();
        sigmaGg.setConstant();
    }
    
    bool twoComp = true;
    if (exclusiveOnly) twoComp = false;
    
    RooLandau* pdfTwoGamma = new RooLandau("ptTwoGamma","ptTwoGamma", *pt, muGg, sigmaGg);
    //RooAbsPdf* pdfTwoGamma = ws->pdf("ptkTwoGammaToMuLow");
    
    //RooRealVar *pt0 = new RooRealVar("pt0","pt0", 0.7, 0.25, 5);
    RooRealVar *bInc = new RooRealVar("bInc","bInc", 2, 0.1, 7);
    RooRealVar *nInc = new RooRealVar("nInc","nInc", 3.5, 1, 7);
    //RooGenericPdf* pdfExtraBkg = new RooGenericPdf("ptBkg","Inclusive #gamma#gamma PDF","fTrkTrkPt/((1.+(fTrkTrkPt/pt0)**2)**nInc)", RooArgSet(*pt, *pt0, *nInc)) ;
    RooGenericPdf* pdfExtraBkg = new RooGenericPdf("ptBkg","Inclusive #gamma#gamma PDF","(2*fTrkTrkPt*(1.+(fTrkTrkPt**2)*(bInc/nInc))**(-nInc))", RooArgSet(*pt, *bInc, *nInc)) ;
    
    RooRealVar* yieldTwoGamma = new RooRealVar("yieldTwoGamma","yieldTwoGamma", 100, 0, 5000);
    RooRealVar* yieldExtraBkg = new RooRealVar("yieldExtraBkg","yieldExtraBkg", 0, 0, 2000);
    
    RooAbsPdf* fitModel = nullptr;
    
    if (twoComp) {
        fitModel = new RooAddPdf("fit", "fit", RooArgList(*pdfTwoGamma, *pdfExtraBkg), RooArgList(*yieldTwoGamma, *yieldExtraBkg), kFALSE);
        fitModel->fitTo(*dataw_bkg, Extended(), Minos(true), Strategy(1));
    }
    else {
        fitModel = new RooAddPdf("fit", "fit", RooArgList(*pdfTwoGamma), RooArgList(*yieldTwoGamma), kFALSE);
        fitModel->fitTo(*dataw_bkg, Range(0, fitRangeMax), Extended(), Minos(true), Strategy(1));
        //fitModel->fitTo(*roo_hist_bkg, Extended(), Minos(true), Strategy(1));
    }
    frameBkg->GetXaxis()->SetRangeUser(0, 1.2);
    frameBkg->Draw() ;
    frameBkg->Print("v");
    
    //legend->Draw("same");
    
    // Plot isolation for QCD component.
    // Eg. plot all events weighted by the sWeight for the QCD component.
    // The SPlot class adds a new variable that has the name of the corresponding
    // yield + "_sw".
    
    
    // Save plots
    bool eps = true;
    string cutType = "";
    if (!useCuts) cutType = "-nocuts";
    string suffix = "";
    if (exp) suffix = "exp";
    else suffix = "powerlaw";
    string suffix2 = "";
    if (exclusiveOnly) suffix2 = "-exclusive-only";
    cv->SaveAs(Form("Plots/%s/Splot%s-%.1f-%.1f-%.1f-%.1f-%s%s.pdf", period.c_str(), cutType.c_str(), mMin, mMax, abs(rapMax), abs(rapMin), suffix.c_str(), suffix2.c_str()));
    if (eps) cv->SaveAs(Form("Plots/%s/Splot%s-%.1f-%.1f-%.1f-%.1f-%s%s.eps", period.c_str(), cutType.c_str(), mMin, mMax, abs(rapMax), abs(rapMin), suffix.c_str(), suffix2.c_str()));
    
    // New Canvas to fit background pT
    TCanvas* cv2 = new TCanvas();
    RooPlot* frameBkg2 = pt->frame() ;
    frameBkg2->SetTitle("p_{T} distribution for background with weights");
    dataw_bkg->plotOn(frameBkg2, Name("bkgData"), DataError(RooAbsData::SumW2), Binning(nBinsPt2) ) ;
    //dataw_bkg->plotOn(frameBkg2, Name("bkgData"), DataError(RooAbsData::SumW2), MarkerColor(kRed)) ;
    //roo_hist_bkg->plotOn(frameBkg2, Name("bkgData2"), DataError(RooAbsData::SumW2)) ;
    fitModel->plotOn(frameBkg2, Name("sumBkg"), LineWidth(2), LineColor(kBlue));
    if (twoComp) {
        fitModel->plotOn(frameBkg2, Components(*pdfTwoGamma), LineStyle(kDashed), LineColor(kBlue), LineWidth(2));
        fitModel->plotOn(frameBkg2, Components(*pdfExtraBkg), LineStyle(kDashed), LineColor(kGreen), LineWidth(2));
    }
    Double_t xPosBkg = newPtMax*0.4;
    //Double_t xPosBkg = 0.2;
    TLatex* txtGammaGamma = new TLatex(xPosBkg, 0.65*frameBkg2->GetMaximum(), Form("#bf{N_{#gamma#gamma} = %.1f #pm %.1f}", yieldTwoGamma->getVal(), yieldTwoGamma->getError()));
    TLatex* txtGammaGamma1 = new TLatex(xPosBkg, 0.55*frameBkg2->GetMaximum(), Form("#bf{#mu_{Landau} = %.3f #pm %.3f}", muGg.getVal(), muGg.getError()));
    TLatex* txtGammaGamma2 = new TLatex(xPosBkg, 0.45*frameBkg2->GetMaximum(), Form("#bf{#sigma_{Landau} = %.3f #pm %.3f}", sigmaGg.getVal(), sigmaGg.getError()));
    
    //TLatex* txtBkg1 = new TLatex(xPosBkg, 0.43*frameBkg2->GetMaximum(), Form("#bf{p_{T0} = %.2f #pm %.2f}", pt0->getVal(), pt0->getError()));
    if (twoComp) {
        TLatex* txtBkg = new TLatex(xPosBkg, 0.35*frameBkg2->GetMaximum(), Form("#bf{N_{Bkg} = %.1f #pm %.1f}", yieldExtraBkg->getVal(), yieldExtraBkg->getError()));
        TLatex* txtBkg1 = new TLatex(xPosBkg, 0.25*frameBkg2->GetMaximum(), Form("#bf{b_{inc} = %.2f #pm %.2f}", bInc->getVal(), bInc->getError()));
        TLatex* txtBkg2 = new TLatex(xPosBkg, 0.15*frameBkg2->GetMaximum(), Form("#bf{n_{inc} = %.2f #pm %.2f}", nInc->getVal(), nInc->getError()));
        frameBkg2->addObject(txtBkg); frameBkg2->addObject(txtBkg1); frameBkg2->addObject(txtBkg2);
    }
    
    frameBkg2->addObject(txtGammaGamma); frameBkg2->addObject(txtGammaGamma1); frameBkg2->addObject(txtGammaGamma2);
    frameBkg2->GetXaxis()->SetRangeUser(0, newPtMax);
    
    // Chi2
    /*
     TH1* histModel = data->createHistogram("fTrkTrkM", nBinsPt);
     histBkg->SetBins(int(nBinsPt/2), 0, fitRangeMax);
     RooDataHist* roo_hist_bkg = new RooDataHist("roo_hist", "roo_hist", RooArgSet(*m), histBkg);
     */
    //int nDofBkg = fitModel->getParameters(dataw_bkg)->selectByAttrib("Constant",kFALSE)->getSize();
    TH1* histBkg = dataw_bkg->createHistogram("hPtBackground2",*pt, Binning(nBinsPt, 1.e-6, fitRangeMax));
    RooDataHist* roo_hist_bkg = new RooDataHist("roo_hist", "roo_hist", RooArgSet(*pt), histBkg);
    int nDofBkg = histBkg->GetNbinsX()-fitModel->getParameters(dataw_bkg)->selectByAttrib("Constant",kFALSE)->getSize();
    //double yMax = frameJpsi->GetMaximum();
    RooChi2Var* chi2Var2 = new RooChi2Var("chi2","chi2", *fitModel, *roo_hist_bkg, DataError(RooAbsData::Poisson));
    TLatex* txtChiBkg = new TLatex(xPosBkg, 0.55*frameBkg2->GetMaximum(),Form("#bf{#chi^{2}/ndf = %.1f/ %d = %.2f}", chi2Var2->getVal(), nDofBkg, chi2Var2->getVal()/nDofBkg));
    //frameBkg2->addObject(txtChiBkg);
    
    frameBkg2->Draw() ;
    
    // write y and pt info
    TLegend* lgdInfo2 = new TLegend(0.6, 0.73, 0.9, 0.9);
    lgdInfo2->SetMargin(0.05);
    lgdInfo2->SetTextSize(0.04);
    lgdInfo2->AddEntry((TObject*)0, Form("%.2f < y < %.2f", -rapMax, -rapMin), "");
    lgdInfo2->AddEntry((TObject*)0, Form("%.1f < M < %.1f GeV/c^{2}", mMin, mMax), "");
    lgdInfo2->AddEntry((TObject*)0, Form("%.1f < p_{T} < %.1f GeV/c", ptMin, ptMax), "");
    lgdInfo2->Draw();
    
    cv2->SaveAs(Form("Plots/%s/sPlot-bkg%s-%.1f-%.1f-%.1f-%.1f-%s%s.pdf", period.c_str(), cutType.c_str(), mMin, mMax, abs(rapMax), abs(rapMin), suffix.c_str(), suffix2.c_str()));
    if (eps) cv2->SaveAs(Form("Plots/%s/sPlot-bkg%s-%.1f-%.1f-%.1f-%.1f-%s%s.eps", period.c_str(), cutType.c_str(), mMin, mMax, abs(rapMax), abs(rapMin), suffix.c_str(), suffix2.c_str()));
    
    /*
     // Look at V0CCounts
     cv->cd(2);
     gPad->SetLeftMargin(0.15) ;
     RooRealVar* v0Ccounts = ws->var("fV0CCounts");
     RooPlot* frame6 = v0Ccounts->frame() ;
     //dataw_jpsi->plotOn(frame6, DataError(RooAbsData::SumW2) ) ;
     data->plotOn(frame6, DataError(RooAbsData::SumW2) ) ;
     
     frame6->SetTitle("V0C Counts");
     frame6->Draw() ;
     
     // Look at 2D hist of V0CCounts et pt
     TCanvas* c4 = new TCanvas("ptV0", "ptV0", 600, 600);
     TH1* hh_data = dataw_jpsi->createHistogram("fTrkTrkPt,fV0CCounts",40,10) ;
     gPad->SetLeftMargin(0.15) ; gPad->SetRightMargin(0.15) ; hh_data->Draw("colz") ;
     
     c4->SaveAs(Form("Plots/Pt-V0-2D-%s.pdf", period.c_str()));
     */
    
    /*
     // create new pdf to have the full range
     RooRealVar* pt2 = ws->var("fTrkTrkPt");
     RooRealVar muGg2("muGg2", "muGg2", muGg.getVal());
     RooRealVar sigmaGg2("sigmaGg2", "sigmaGg2",sigmaGg.getVal());
     RooDataSet* dataw_bkg2 = new RooDataSet(sdata->GetName(),sdata->GetTitle(),sdata,*sdata->get(),0,"fbkg_sw") ;
     RooLandau* pdfTwoGamma2 = new RooLandau("ptTwoGamma2","ptTwoGamma2", *pt, muGg2, sigmaGg2);
     pdfTwoGamma2->fitTo(*dataw_bkg2, Range(0, ptMax), Strategy(0));
     */
    
    // new canvas for jpsi
    
    // Fit pt distrib of J/Psi with exclusive / dissociative components
    RooAbsPdf* ptmodel = ws->pdf("ptfit");
    RooAbsPdf* ptJpsiExclusive = ws->pdf("ptJpsiExclusive");
    RooAbsPdf* ptJpsiDissociative = ws->pdf("ptJpsiDissociative");
    RooAbsPdf* ptJpsiGammaPb = ws->pdf("ptJpsiGammaPb");
    RooAbsPdf* ptJpsiInclusive = ws->pdf("ptJpsiInclusive");
    ptMax = 3.0;
    ptmodel->fitTo(*dataw_jpsi, Range(0, ptMax), Extended(), Minos(true), Strategy(1));
    
    TCanvas* cv3 = new TCanvas();
    RooPlot* frameJpsi2 = pt->frame() ;
    frameJpsi2->SetTitle("p_{T} distribution for jpsi with weights");
    //dataw_jpsi->plotOn(frameJpsi2, Name("jpsiData"), DataError(RooAbsData::SumW2), Binning(nBinsPt2) ) ;
    dataw_jpsi->plotOn(frameJpsi2, Name("jpsiData"), DataError(RooAbsData::SumW2)) ;
    //dataw_bkg->plotOn(frameBkg2, Name("bkgData"), DataError(RooAbsData::SumW2), MarkerColor(kRed)) ;
    //roo_hist_bkg->plotOn(frameBkg2, Name("bkgData2"), DataError(RooAbsData::SumW2)) ;
    ptmodel->plotOn(frameJpsi2, Name("jpsimodel"), LineWidth(2), LineColor(kBlue));
    ptmodel->plotOn(frameJpsi2, Components(*ptJpsiExclusive), LineStyle(kDashed), LineColor(kRed));
    ptmodel->plotOn(frameJpsi2, Components(*ptJpsiGammaPb), LineStyle(kDashed), LineColor(kGray));
    frameJpsi2->GetXaxis()->SetRangeUser(0, ptMax);
    
    // Write chi2
    double xPosJpsi = ptMax/3;
    int nDof = ptmodel->getParameters(dataw_jpsi)->selectByAttrib("Constant",kFALSE)->getSize();
    double yMax = frameJpsi2->GetMaximum();
    //TLatex* txtChi = new TLatex(xPosJpsi/2, 0.93*yMax,Form("#chi^{2}/ndf = %.3f", frameJpsi->chiSquare( "ptfit", "jPsiData", nDof)));
    //frameJpsi2->addObject(txtChi);
    cout << endl << endl << endl << endl;
    
    frameJpsi2->Draw();
    
    // write legend with parameters
    RooRealVar* bExc = ws->var("bExc");
    RooRealVar* yieldJpsiExclusive = ws->var("yieldJpsiExclusive");
    RooRealVar* yieldJpsiGammaPb = ws->var("yieldJpsiGammaPb");
    TLegend* lgd3 = new TLegend(0.7, 0.7, 0.9, 0.9);
    lgd3->AddEntry((TObject*)0, Form("N_{J/#psi} = %.1f #pm %.1f", yieldJpsiExclusive->getVal(), yieldJpsiExclusive->getError()), "");
    lgd3->AddEntry((TObject*)0, Form("b_{exc} = %.2f #pm %.2f", bExc->getVal(), bExc->getError()), "");
    lgd3->AddEntry((TObject*)0, Form("N_{#gamma-Pb} = %.1f #pm %.1f", yieldJpsiGammaPb->getVal(), yieldJpsiGammaPb->getError()), "");
    lgd3->Draw();
    
    cv3->SaveAs(Form("Plots/%s/sPlot-jpsi%s-%.1f-%.1f-%.1f-%.1f-%s%s.pdf", period.c_str(), cutType.c_str(), mMin, mMax, abs(rapMax), abs(rapMin), suffix.c_str(), suffix2.c_str()));
    
    // And save pt shapes
    TFile* fSplot = new TFile("splot.root","RECREATE");
    int nBinsModel = 1000;
    int nBinsM = 50;
    
    // mass histo
    TH1* hMass = data->createHistogram("hMass", *m, Binning(nBinsM, mMin, mMax));
    TH1* hModel = model->createHistogram("hModel", *m, Binning(nBinsModel, mMin, mMax));
    TH1* hJpsi = jPsiModel->createHistogram("hJpsi", *m, Binning(nBinsModel, mMin, mMax));
    TH1* hBkg = bkgModel->createHistogram("hBkg", *m, Binning(nBinsModel, mMin, mMax));
    
    // pt histos
    TH1* hPtJpsi = dataw_jpsi->createHistogram("hPtJpsi", *pt, Binning(nBinsPt2, 0, ptMax));
    TH1* hPtBackground = dataw_bkg->createHistogram("hPtBackground", *pt, Binning(nBinsPt2, 0, ptMax));
    TH1* hPtTwoGamma = pdfTwoGamma->createHistogram("hPtTwoGamma", *pt, Binning(nBinsModel, 0, fitRangeMax));
    
    hModel->Scale((double)nBinsModel/nBinsM*(jPsiYield->getVal()+bkgYield->getVal())/hModel->Integral());
    hJpsi->Scale((double)nBinsModel/nBinsM*jPsiYield->getVal()/hJpsi->Integral());
    hBkg->Scale((double)nBinsModel/nBinsM*bkgYield->getVal()/hBkg->Integral());
    hPtTwoGamma->Scale((double)nBinsModel/(nBinsPt2*fitRangeMax/ptMax)*yieldTwoGamma->getVal()/hPtTwoGamma->Integral());
    
    for (int i = 0; i<nBinsModel; i++) {
        hModel->SetBinError(i+1, 0);
        hJpsi->SetBinError(i+1, 0);
        hBkg->SetBinError(i+1, 0);
        hPtTwoGamma->SetBinError(i+1, 0);
    }
    
    hModel->SetOption("hist");
    hJpsi->SetOption("hist");
    hBkg->SetOption("hist");
    hPtTwoGamma->SetOption("hist");
    
    hMass->Write();
    hModel->Write();
    hJpsi->Write();
    hBkg->Write();
    hPtJpsi->Write();
    hPtBackground->Write();
    hPtTwoGamma->Write();
    
    fSplot->Close();
    
    // Draw final plots
    TCanvas* cvFinal = new TCanvas("cvFinal", "cvFinal", 1500, 400);
    cvFinal->Divide(3);
    cvFinal->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    frame->SetTitle("");
    frame->Draw();
    lgdInfo->Draw();
    
    cvFinal->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    frameJpsi->SetTitle("");
    frameJpsi->Draw();
    
    cvFinal->cd(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    frameBkg2->SetTitle("");
    frameBkg2->Draw();
    
    cvFinal->SaveAs(Form("Plots/%s/Final-sPlot-bkg%s-%.1f-%.1f-%.1f-%.1f-%s%s.pdf", period.c_str(), cutType.c_str(), mMin, mMax, abs(rapMax), abs(rapMin), suffix.c_str(), suffix2.c_str()));
    if (eps) cvFinal->SaveAs(Form("Plots/%s/Final-sPlot-bkg%s-%.1f-%.1f-%.1f-%.1f-%s%s.eps", period.c_str(), cutType.c_str(), mMin, mMax, abs(rapMax), abs(rapMin), suffix.c_str(), suffix2.c_str()));
    
    
}


void DrawBkgPt(string rootfilePath, string period, bool exclusiveOnly) {
    
    gStyle->SetTextSize(.04);
    
    string suffix = "";
    if (exclusiveOnly) suffix = "-exclusive-only";
    TFile* fTemplates = new TFile(Form("%s/sPlotTemplates-%s%s.root", rootfilePath.c_str(), period.c_str(), suffix.c_str()),"READ");
    
    TH1* hPtBackground = (TH1*)fTemplates->Get("hPtBackground__fTrkTrkPt");
    hPtBackground->SetTitle("Background p_{T} distribution with sPlot");
    TH1* hSubtraction = (TH1*)fTemplates->Get("hSubtraction");
    TH1* hSubSmooth = (TH1*)fTemplates->Get("hSubSmooth");
    TH1* hPtMc = (TH1*)fTemplates->Get("hPtMc__fTrkTrkPt");
    //TH1* hPtMc = (TH1*)fTemplates->Get("hPtMc");
    TCanvas* cv3 = new TCanvas();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hPtBackground->SetLineColor(kBlack);
    hPtBackground->SetLineWidth(3);
    hPtBackground->SetMinimum(0);
    hPtBackground->SetMaximum(1.1*hPtMc->GetMaximum());
    hPtBackground->Draw("hist");
    hPtMc->SetLineColor(kRed);
    hPtMc->SetLineWidth(2);
    hSubtraction->SetLineColor(kRed+2);
    hSubtraction->SetLineWidth(2);
    hSubSmooth->SetLineColor(kBlue);
    hSubSmooth->SetLineWidth(3);
    hPtMc->Draw("hist same");
    hSubtraction->Draw("hist same");
    hSubSmooth->Draw("hist same");
    /*
     hPtBkgSmooth->SetLineColor(kRed);
     hPtBkgSmooth->Draw("hist same");
     */
    
    // Write number of candidates
    Double_t yMax = hPtBackground->GetMaximum();
    TLatex* t1 = new TLatex(2, 0.2*yMax, Form("Number of bkg = %.1f", hPtBackground->Integral()));
    TLatex* t2 = new TLatex(2, 0.3*yMax, Form("Number of #gamma#gamma = %.1f", hPtMc->Integral()));
    t1->Draw(); t2->Draw();
    
    TLegend* lgd = new TLegend(0.5, 0.5, 0.9, 0.9);
    lgd->AddEntry(hPtBackground, "Data with sPlot", "l");
    lgd->AddEntry(hPtMc, "MC #gamma#gamma", "l");
    lgd->AddEntry(hSubtraction, "Subtraction: Data - MC", "l");
    lgd->AddEntry(hSubSmooth, "Smoothed subtraction", "l");
    lgd->Draw();
    cv3->SaveAs(Form("Plots/%s/pt-splotsmooth%s.pdf", period.c_str(), suffix.c_str()));
}

//____________________________________
void DrawWeights(RooWorkspace *ws){
    
    
    //get what we need of the workspace
    RooRealVar* m = ws->var("fTrkTrkM");
    RooRealVar* pt = ws->var("fTrkTrkPt");
    
    RooDataSet* sData = (RooDataSet*) ws->data("dataWithSWeights");
    
    //std::cout << std::endl << std::endl << std::endl << std::endl;
    const int nEntries = (int) sData->sumEntries();
    Double_t massForWeights[nEntries], ptForWeights[nEntries], jPsiWeight[nEntries], bkgWeight[nEntries];
    
    for (Int_t i=0; i < (Int_t)nEntries; i ++ ) {
        massForWeights[i] = sData->get(i)->getRealValue("fTrkTrkM");
        ptForWeights[i] = sData->get(i)->getRealValue("fTrkTrkPt");
        
        jPsiWeight[i] = sData->get(i)->getRealValue("fsigJpsi_sw");
        bkgWeight[i] = sData->get(i)->getRealValue("fbkg_sw");
        
    }
    //std::cout << std::endl << std::endl << std::endl << std::endl;
    
    // TGraphs of sWeights = f(m) or sWeights = f(pt)
    TGraph* gWeightsMJpsi = new TGraph(nEntries,massForWeights, jPsiWeight);
    TGraph* gWeightsPtJpsi = new TGraph(nEntries,ptForWeights, jPsiWeight);
    TGraph* gWeightsMBkg = new TGraph(nEntries,massForWeights, bkgWeight);
    TGraph* gWeightsPtBkg = new TGraph(nEntries,ptForWeights, bkgWeight);
    
    gWeightsMJpsi->SetTitle("JPsi sWeights = f(m)");
    gWeightsPtJpsi->SetTitle("JPsi sWeights = f(pt)");
    gWeightsMBkg->SetTitle("Background sWeights = f(m)");
    gWeightsPtBkg->SetTitle("Background sWeights = f(pt)");
    
    gWeightsMJpsi->GetXaxis()->SetTitle("m");
    gWeightsPtJpsi->GetXaxis()->SetTitle("pt");
    gWeightsMBkg->GetXaxis()->SetTitle("m");
    gWeightsPtBkg->GetXaxis()->SetTitle("pt");
    
    gWeightsMJpsi->GetYaxis()->SetTitle("J/Psi sWeights");
    gWeightsPtJpsi->GetYaxis()->SetTitle("J/Psi sWeights");
    gWeightsMBkg->GetYaxis()->SetTitle("Bkg sWeights");
    gWeightsPtBkg->GetYaxis()->SetTitle("Bkg sWeights");
    
    /*
     TCanvas* cv = new TCanvas("sweights","sweights", 800, 1200) ;
     cv->Divide(2,3);
     
     cv->cd(1);
     gPad->SetLeftMargin(0.2);
     gWeightsMJpsi->Draw("A*");
     
     cv->cd(2);
     gPad->SetLeftMargin(0.2);
     gWeightsPtJpsi->Draw("A*");
     
     cv->cd(3);
     gPad->SetLeftMargin(0.2);
     gWeightsMBkg->Draw("A*");
     
     cv->cd(4);
     gPad->SetLeftMargin(0.2);
     gWeightsPtBkg->Draw("A*");
     
     // 2d histograms
     TH2F* hWeightsBkg2d = new TH2F("hWeightsBkg2d", "Background sWeights)", 50, mMin, mMax, 100, ptMin, 4.5);
     TH2F* hWeightsJpsi2d = new TH2F("hWeightsJpsi2d", "J/Psi sWeights", 50, mMin, mMax, 100, ptMin, 4.5);
     for (int i = 0; i<nEntries; i++) {
     hWeightsBkg2d->Fill(massForWeights[i], ptForWeights[i], bkgWeight[i]);
     hWeightsJpsi2d->Fill(massForWeights[i], ptForWeights[i], jPsiWeight[i]);
     }
     cv->cd(5);
     gPad->SetLeftMargin(0.2);
     hWeightsBkg2d->GetXaxis()->SetTitle("Mass");
     hWeightsBkg2d->GetYaxis()->SetTitle("Pt");
     //hWeightsBkg2d->Draw("CONTZ");
     hWeightsBkg2d->Draw("COLZ");
     
     cv->cd(6);
     gPad->SetLeftMargin(0.2);
     hWeightsJpsi2d->GetXaxis()->SetTitle("Mass");
     hWeightsJpsi2d->GetYaxis()->SetTitle("Pt");
     hWeightsJpsi2d->Draw("COLZ");
     
     // Save plots
     string cutType = "";
     if (!useCuts) cutType = "-nocuts";
     cv->SaveAs(Form("Plots/Splot-weights-%s%s-%.1f-%.1f.pdf", period.c_str(), cutType.c_str(), mMin, mMax));
     */
}

void Splot(string rootfilePath = "", string rootfilePathMc = "", std::vector<string> periods = {"LHC16r", "LHC16s"}) {
    
    Double_t mMin = 2., mMax = 3.5, ptMin = 0., ptMax = 3.5;
    Double_t rapMin = -4., rapMax = -2.5;
    bool useCuts = true, exp = false, exclusiveOnly = false;
    bool logScale = false, drawPulls = false;
    
    /*
    gROOT->ProcessLine(".L Include/ExtendedCrystalBall.cxx+") ;
    gSystem->Load("./Include/ExtendedCrystalBall_cxx.so") ;
     */
    
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
        AddModel(wspace, period, mMin, mMax);
        AddPtJpsiModel(wspace, rootfilePathMc, period, mMin, mMax, ptMin, ptMax, rapMin, rapMax, exp, exclusiveOnly);
        wspace->Print();
        
        DoSPlot(wspace, mMax);

        MakePlots(wspace, rootfilePath, rootfilePathMc, period, exp, exclusiveOnly, useCuts, mMin, mMax, ptMin, ptMax, rapMin, rapMax);
        //DrawWeights(wspace);
        
        //cleanup
        delete wspace;
        //DrawBkgPt(rootfilePath, period, exclusiveOnly);
    }
    
    
}


