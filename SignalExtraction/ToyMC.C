#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <stdio.h>
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
#include "TCanvas.h"
#include "TAxis.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBasket.h"
#include "Riostream.h"

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
#include "RooHistPdf.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPlot.h"

using namespace RooFit;
using namespace std;

Double_t mMin = 2.5, mMax = 3.5;
Double_t ptMin = 0., ptMax = 1.2;
//int ptBinNumber = int(10*ptMax);
int ptBinNumber = 120;


std::map<std::string, Int_t> nEventsSignal, nEventsBkg;

void InitiateNumbers(std::map<std::string, Int_t> &nEventsSignal, std::map<std::string, Int_t> &nEventsBkg) {
    nEventsSignal["LHC16r"] = 1300;
    nEventsSignal["LHC16s"] = 600;
    nEventsBkg["LHC16r"] = 500;
    nEventsBkg["LHC16s"] = 100;
}


template <typename T>
bool contains(std::vector<T> & listOfElements, const T& element)
{
    // Find the iterator if element in list
    auto it = std::find(listOfElements.begin(), listOfElements.end(), element);
    //return if iterator points to end or not. It points to end then it means element
    // does not exists in list
    return it != listOfElements.end();
}


void AddModel(RooWorkspace* ws) {
    // Define model
    
    //Crystal ball for the J/psi
    
    RooRealVar m = *ws->var("m");
    
    RooRealVar mean_jpsi("mean_jpsi","mean_jpsi",3.1,2.9,3.2);
    RooRealVar sigma_jpsi("sigma_jpsi","sigma_jpsi",0.07,0,0.15);
    RooRealVar alpha_jpsi("alpha_jpsi","alpha_jpsi",1, 0, 5);
    RooRealVar n_jpsi("n_jpsi","n_jpsi",5, 0, 10);
    RooCBShape *jpsi = new RooCBShape("jpsi","crystal ball PDF", m, mean_jpsi, sigma_jpsi,alpha_jpsi,n_jpsi);
    
    //Exponential background
    //RooRealVar a1("a1","a1",-1.,-5,-0.);
    // Constant function
    //RooRealVar a1("a1","a1",0);
    RooRealVar a1("a1","a1",-1, -5, 1);
    RooExponential *bkg = new RooExponential("exp","exp",m,a1);
    
    //RooRealVar fsig("fsig","signalPhi",0.1,0.,1.);
    RooRealVar fsigJPsi("fsigJPsi","signalJPsi",5.e3,1.e1,1.e5);
    RooRealVar fbkg("fbkg","fbkg",5.e3,1.e1,1.e5);
    
    RooAbsPdf* model = new RooAddPdf("mfit", "mfit", RooArgList(*jpsi, *bkg), RooArgList(fsigJPsi, fbkg), kFALSE);
    
    ws->import(*model);
}

void GetPtHistMC(RooWorkspace* ws, std::string rootfilePathMC, std::string period){
    TFile *fSimu = new TFile(Form("%s/toy_MC_%s.root", rootfilePathMC.c_str(), period.c_str()),"READ");
    TTree* tSignal = (TTree*)fSimu->Get("tSignal");
    TTree* tBkg = (TTree*)fSimu->Get("tBkg");
    
    RooRealVar ptSignal("pt","Dimuon p_{T} (GeV/c)",0,ptMax);
    RooRealVar ptBkg("pt","Dimuon p_{T} (GeV/c)",0,ptMax);
    RooDataSet* dataSignal = new RooDataSet("data","data",RooArgSet(ptSignal),Import(*tSignal));
    RooDataSet* dataBkg = new RooDataSet("data","data",RooArgSet(ptBkg),Import(*tBkg));
    
    TH1F* histSignalPt = new TH1F("hPtSignal", "hPtSignal", ptBinNumber, 0, ptMax);
    TH1F* histBkgPt = new TH1F("hPtBkg", "hPtBkg", ptBinNumber, 0, ptMax);
    tSignal->Draw("pt>>hPtSignal");
    tBkg->Draw("pt>>hPtBkg");
    RooDataHist* ptSignalHist = new RooDataHist("ptSignalHist","ptSignalHist", RooArgList(ptSignal),histSignalPt);
    RooDataHist* ptBkgHist = new RooDataHist("ptBkgHist","ptBkgHist", RooArgList(ptBkg),histBkgPt);
    RooHistPdf* ptSignalPdf = new RooHistPdf("ptSignalPdf", "ptSignalPdf", ptSignal, *ptSignalHist);
    RooHistPdf* ptBkgPdf = new RooHistPdf("ptBkgPdf", "ptBkgPdf", ptBkg, *ptBkgHist);
    
    ws->import(*ptBkgHist);
    ws->import(*ptSignalPdf);
    ws->import(*ptBkgPdf);
    
}

void DoSPlot(RooWorkspace* ws) {
    
    RooAbsPdf* model = ws->pdf("mfit");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    RooRealVar* jPsiYield = ws->var("fsigJPsi");
    RooRealVar* bkgYield = ws->var("fbkg");
    
    model->fitTo(*data, Extended(), Minos(true), Strategy(2));
    
    //The sPlot technique requires that we fix the parameters
    // of the model that are not yields after doing the fit
    RooRealVar* mean_jpsi = ws->var("mean_jpsi");
    RooRealVar* sigma_jpsi = ws->var("sigma_jpsi");
    RooRealVar* alpha_jpsi = ws->var("alpha_jpsi");
    RooRealVar* n_jpsi = ws->var("n_jpsi");
    
    RooRealVar* a1 = ws->var("a1");
    
    mean_jpsi->setConstant();
    sigma_jpsi->setConstant();
    alpha_jpsi->setConstant();
    n_jpsi->setConstant();
    a1->setConstant();
    
    RooMsgService::instance().setSilentMode(true);
    
    //Now we use the SPlot class to add SWeight to our data set
    // based on our model and our yield variables
    
    RooStats::SPlot * sData = new RooStats::SPlot("sData","splot", *data, model, RooArgList(*jPsiYield,*bkgYield) );
    
    //Check Sweight properties
    std::cout << "Check SWeights: " << std::endl;
    
    std::cout << std::endl << "Yield of JPsi is "
    << jPsiYield->getVal() << ". From sWeights it is "
    << sData->GetYieldFromSWeight("fsigJPsi") << std::endl;
    
    std::cout << std::endl << "Yield of bkg is "
    << bkgYield->getVal() << ". From sWeights it is "
    << sData->GetYieldFromSWeight("fbkg") << std::endl;
    
    for (Int_t i=0; i < 10; i ++ )
    {
        std::cout << "JPsi Weight "<< sData->GetSWeight(i, "fsigJPsi")
        << " bkg Yield "<< sData->GetSWeight(i, "fbkg")
        << " Total Weight "<< sData->GetSumOfEventSWeight(i)
        << std::endl;
    }
    
    //import the new data set with Sweight
    
    std::cout << "import new dataset with sWeight" << std::endl;
    ws->import(*data, Rename("dataWithSWeights"));
    
}


void MakePlots(RooWorkspace *ws, std::string rootfilePathMC, std::string period, bool drawPulls, bool write, std::vector <TH1*> histVec) {
    
    
    //get what we need of the workspace
    RooAbsPdf* model = ws->pdf("mfit");
    RooAbsPdf* jPsiModel = ws->pdf("jpsi");
    RooAbsPdf* bkgModel = ws->pdf("exp");
    
    RooAbsPdf* ptSignalPdf = ws->pdf("ptSignalPdf");
    RooAbsPdf* ptBkgPdf = ws->pdf("ptBkgPdf");
    
    RooRealVar* m = ws->var("m");
    RooRealVar* pt = ws->var("pt");
    
    RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");
    
    RooDataHist* ptBkgHist = (RooDataHist*) ws->data("ptBkgHist");
    
    // Get weighted data
    // The SPlot class adds a new variable that has the name of the corresponding
    // yield + "_sw".
    RooDataSet *dataw_jpsi = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"fsigJPsi_sw") ;
    RooDataSet *dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"fbkg_sw") ;
    
    // Extract pt templates to compare them with reconstructed data
    RooRealVar yieldSignal("yieldSignal","yieldSignal",nEventsSignal[period]);
    RooRealVar yieldBkg("yieldBkg","yieldBkg",nEventsBkg[period]);
    RooAbsPdf* ptSignalModel = new RooAddPdf("ptSignalMC", "ptSignalMC", RooArgList(*ptSignalPdf), RooArgList(yieldSignal), kFALSE);
    RooAbsPdf* ptBkgModel = new RooAddPdf("ptBkgMC", "ptBkgMC", RooArgList(*ptBkgPdf), RooArgList(yieldBkg), kFALSE);
    //ptSignalModel->fitTo(*dataw_jpsi, Extended(), Minos(true), Strategy(2));
    //ptBkgModel->fitTo(*dataw_bkg, Extended(), Minos(true), Strategy(2));
    
    RooRealVar* jPsiYield = ws->var("fsigJPsi");
    RooRealVar* bkgYield = ws->var("fbkg");
    
    // Ecrire les valeurs dans un fichier
    if (write) {
        string const nomFichier(Form("Candidate-Numbers/sPlot-Numbers-%s.txt", period.c_str()));
        ofstream monFlux(nomFichier.c_str(), ios::app);
        
        if (monFlux) {
            monFlux << jPsiYield->getVal() << " ";
            monFlux << bkgYield->getVal() << endl;
        }
        else {cout << "File not opened" << std::endl;}
    }
    
    
    
    
    if (!write) {
        //if (true) {
        //make some plots
        TCanvas* cv = new TCanvas("splot","splot",900,900) ;
        if (drawPulls) cv->Divide(3, 3);
        else {cv->Divide(3, 2); cv->SetCanvasSize(900, 600);}
        
        // First draw mass fit
        cv->cd(1);
        gPad->SetLeftMargin(0.15);
        RooPlot* frame = m->frame();
        data->plotOn(frame);
        model->plotOn(frame);
        model->plotOn(frame, Components(*model), LineStyle(kDashed), LineColor(kRed));
        model->plotOn(frame, Components(*jPsiModel), LineStyle(kDashed), LineColor(kRed));
        model->plotOn(frame, Components(*bkgModel), LineStyle(kDashed), LineColor(kGreen));
        frame->SetTitle("Fit to model to discriminating variable");
        double yMax0 = frame->GetMaximum();
        double x0 = mMin+(mMax-mMin)*0.2;
        TText* txt0 = new TText(x0,0.6*yMax0,Form("%.1f J/Psi", jPsiYield->getVal()));
        TText* txt1 = new TText(x0,0.5*yMax0,Form("%.1f Bkg", bkgYield->getVal()));
        //txt3->SetTextSize(0.05) ;
        frame->addObject(txt0) ;
        frame->addObject(txt1) ;
        frame->Draw();
        
        // Draw pt with weights
        // For signal
        cv->cd(2);
        gPad->SetLeftMargin(0.15);
        RooPlot* frame2 = pt->frame() ;
        dataw_jpsi->plotOn(frame2, DataError(RooAbsData::SumW2) ) ;
        ptSignalModel->plotOn(frame2, LineStyle(kDashed), LineColor(kRed));
        frame2->SetTitle("p_{T} distribution for J/#Psi with weights");
        frame2->GetXaxis()->SetRangeUser(0,2);
        frame2->Draw() ;
        TLegend* lgdJpsi = new TLegend(0.3, 0.7, 0.9, 0.9);
        lgdJpsi->AddEntry(frame2->findObject("h_dataWithSWeights"), "reconstructed signal", "LP");
        lgdJpsi->AddEntry(frame2->findObject("ptSignalMC_Norm[pt]"), "original distribution", "L");
        lgdJpsi->Draw();
        
        
        // Add 2 landau fits for bkg pt
        RooRealVar* mean = new RooRealVar("mean","mean",0.06, 0.01,0.1);
        RooRealVar* sigma = new RooRealVar("sigma","sigma", 0.03, 0.01, 0.1); //3.43,3.73);
        RooRealVar* yieldTwoGamma = new RooRealVar("yieldTwoGamma","yieldTwoGamma", 300, 10., 1000);
        RooAbsPdf* ptTwoGamma = new RooLandau("ptTwoGamma","ptTwoGamma", *pt, *mean, *sigma);
        RooAbsPdf* model1 = new RooAddPdf("fit1", "fit1", RooArgList(*ptTwoGamma), RooArgList(*yieldTwoGamma), kFALSE);
        
        RooRealVar* mean2 = new RooRealVar("mean2","mean2",0.06, 0.01,0.1);
        RooRealVar* sigma2 = new RooRealVar("sigma2","sigma2", 0.03, 0.01, 0.1); //3.43,3.73);
        RooRealVar* yieldTwoGamma2 = new RooRealVar("yieldTwoGamma2","yieldTwoGamma2", 300, 10., 1000);
        RooAbsPdf* ptTwoGamma2 = new RooLandau("ptTwoGamma2","ptTwoGamma2", *pt, *mean2, *sigma2);
        RooAbsPdf* model2 = new RooAddPdf("fit2", "fit2", RooArgList(*ptTwoGamma2), RooArgList(*yieldTwoGamma2), kFALSE);
        
        
        model1->fitTo(*dataw_bkg, Extended(), Range(0, 0.38), Minos(true), Strategy(1));     // distribution reconstruite avec sPlot
        model2->fitTo(*ptBkgHist, Extended(), Range(0, 0.38), Minos(true), Strategy(1));    // distribution initiale
        
        // import Landau distributions
        ws->import(*model1);
        ws->import(*model2);
        ws->import(*yieldTwoGamma);
        ws->import(*yieldTwoGamma2);
        
        std::cout << std::endl << "reconstructed" << std::endl;
        std::cout << "mu1 = " << Form("%.4f",mean->getVal()) << " pm " << Form("%.4f",mean->getError()) << std::endl;
        std::cout << "sigma1 = " << Form("%.4f",sigma->getVal()) << " pm " << Form("%.4f",sigma->getError()) << std::endl;
        std::cout << std::endl << "initial distribution" << std::endl;
        std::cout << "mu2 = " << Form("%.4f",mean2->getVal()) << " pm " << Form("%.4f",mean2->getError()) << std::endl;
        std::cout << "sigma1 = " << Form("%.4f",sigma2->getVal()) << " pm " << Form("%.4f",sigma2->getError()) << std::endl << std::endl;
        
        //std::vector <TH1*> histVec = {hMuRec, hMuIni, hSigRec, hSigIni};
        histVec[0]->Fill(mean->getVal()-mean2->getVal());
        histVec[1]->Fill(sigma->getVal()-sigma2->getVal());
        
        
        // And for background
        cv->cd(3);
        gPad->SetLeftMargin(0.15);
        RooPlot* frame3 = pt->frame() ;
        dataw_bkg->plotOn(frame3,DataError(RooAbsData::SumW2) ) ;
        ptBkgModel->plotOn(frame3, LineStyle(kDashed), LineColor(kRed));
        // plot Landau distributions
        //plot landau distributions
        model1->plotOn(frame3, LineColor(kBlue));
        model2->plotOn(frame3, LineColor(kGreen));
        frame3->SetTitle("p_{T} distribution for #gamma#gamma #rightarrow #mu#mu with weights");
        frame3->GetXaxis()->SetRangeUser(0,2);
        frame3->Draw() ;
        frame3->Print("v");
        TLegend* lgdBkg = new TLegend(0.3, 0.7, 0.9, 0.9);
        lgdBkg->AddEntry(frame3->findObject("h_dataWithSWeights"), "reconstructed signal", "LP");
        lgdBkg->AddEntry(frame3->findObject("ptBkgMC_Norm[pt]"), "original distribution", "L");
        lgdBkg->Draw();
        //return;
        
        
        // Then finally, draw quality plots
        // For the mass fit
        cv->cd(4);
        gPad->SetLeftMargin(0.15);
        data->plotOn(frame, Binning(120));
        model->plotOn(frame, Binning(120));
        RooHist* hpull = frame->pullHist();
        hpull->GetXaxis()->SetRangeUser(mMin, mMax);
        hpull->SetTitle("(data - fit)/#sigma of the mass fit");
        hpull->Draw("");
        cv->cd(7);
        gPad->SetLeftMargin(0.15);
        RooHist* hresid = frame->residHist();
        hresid->GetXaxis()->SetRangeUser(mMin, mMax);
        hresid->SetTitle("Residuals (data - fit) of the mass fit");
        hresid->Draw("");
        
        // Draw quality histograms
        // First (data-fit)/sigma
        // For signal pt
        cv->cd(5);
        gPad->SetLeftMargin(0.15);
        dataw_jpsi->plotOn(frame2, Binning(120));
        ptSignalModel->plotOn(frame2, Binning(120));
        RooHist* hpull2 = frame2->pullHist();
        hpull2->GetXaxis()->SetRangeUser(ptMin, ptMax);
        hpull2->SetTitle("(reconstructed - original)/#sigma of J/#Psi p_{T}");
        hpull2->Draw("");
        
        // For background pt
        cv->cd(6);
        gPad->SetLeftMargin(0.15);
        dataw_bkg->plotOn(frame3, Binning(120));
        ptBkgModel->plotOn(frame3, Binning(120));
        RooHist* hpull3 = frame3->pullHist();
        hpull3->GetXaxis()->SetRangeUser(ptMin, ptMax);
        hpull3->SetTitle("(reconstructed - original)/#sigma of bkg p_{T}");
        hpull3->Draw("");
        
        // Draw pull histograms = data-fit
        if (drawPulls) {
            cv->cd(8);
            gPad->SetLeftMargin(0.15);
            RooHist* hresid2 = frame2->residHist();
            hresid2->GetXaxis()->SetRangeUser(ptMin, ptMax);
            hresid2->SetTitle("Residuals (reconstructed - original) of J/#Psi p_{T}");
            hresid2->Draw("");
            
            cv->cd(9);
            gPad->SetLeftMargin(0.15);
            RooHist* hresid3 = frame3->residHist();
            hresid3->GetXaxis()->SetRangeUser(ptMin, ptMax);
            hresid3->SetTitle("Residuals (reconstructed - original) of bkg p_{T}");
            hresid3->Draw("");
        }
        
        cv->SaveAs(Form("Plots/toyMC-splot-%s-%d-%d.pdf", period.c_str(), nEventsBkg[period], nEventsSignal[period]));
    }
    
    
    // Now save pt shapes
    TFile* fNewTemplates = new TFile(Form("%s/sPlotTemplates-%s.root", rootfilePathMC.c_str(), period.c_str()),"RECREATE");
    
    TH1* hPtExclusive = dataw_jpsi->createHistogram("hPtExclusive",*pt, Binning(ptBinNumber, 0, ptMax)) ;
    hPtExclusive->Write();
    
    // Import original pt shape
    TFile *fSimu = new TFile(Form("%s/toy_MC_%s.root", rootfilePathMC.c_str(), period.c_str()),"READ");
    TTree* tBkg = (TTree*)fSimu->Get("tBkg");
    TH1F* hPtBkgOriginal = new TH1F("hPtBkgOriginal", "hPtBkgOriginal", ptBinNumber, 0, ptMax);
    tBkg->Draw("pt>>hPtBkgOriginal");
    
    fNewTemplates->cd();
    hPtBkgOriginal->Write();
    TH1* hPtBackground = dataw_bkg->createHistogram("hPtBackground",*pt, Binning(ptBinNumber, 0, ptMax)) ;
    hPtBackground->Write();
    
    // Create a smooth histogram (the size of the smooting depends on pt, should be logarithmic so that the smoothing is a lot more agressive for high pt)
    TH1* hPtBkgSmooth = new TH1F("hPtBkgSmooth", "smoothed background pt template", ptBinNumber, 0, ptMax);
    double negValue = -10;    // in case the value is negative, pb log(neg) so set log to this negvalue
    for (int k = 0; k<ptBinNumber; k++) {
        int binNum = k+1;
        double smoothVal = 0;
        int nSmoothBins = 1;
        if (hPtBackground->GetBinContent(binNum) > 0) {smoothVal = TMath::Log(hPtBackground->GetBinContent(binNum));}
        else smoothVal = -10;    // by default say there's 10^-2 background events (originally) in this bin
        double size = binNum*0.8;    // number of bins to smooth
        if (size > 10) size = 10.;
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
    
    fNewTemplates->cd();
    hPtBkgSmooth->Write();
    fNewTemplates->Close();
    fSimu->Close();
    
    ws->import(*dataw_jpsi, Rename("dataJpsi"));
    ws->import(*dataw_bkg, Rename("dataBkg"));
    
    ws->import(*jPsiModel);
    ws->import(*bkgModel);
    
    
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

void ExportHist(RooWorkspace* ws) {
    
    // First get quantities needed
    RooDataSet* data = (RooDataSet*) ws->data("data");
    //RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");
    RooDataSet* dataw_jpsi = (RooDataSet*) ws->data("dataJpsi");
    RooDataSet* dataw_bkg = (RooDataSet*) ws->data("dataBkg");
    
    RooRealVar* m = ws->var("m");
    RooRealVar* pt = ws->var("pt");
    
    
    // get yields
    RooRealVar* jPsiYield = ws->var("fsigJPsi");
    RooRealVar* bkgYield = ws->var("fbkg");
    
    // get model for mass
    RooAbsPdf* model = ws->pdf("mfit");
    RooAbsPdf* jPsiModel = ws->pdf("jpsi");
    RooAbsPdf* bkgModel = ws->pdf("exp");
    
    // gets histos for mass
    int nBinsMass = 50;
    int nBinsModel = 1000;
    double mMin = 2.5, mMax = 3.5;
    TH1* hMass = data->createHistogram("hMass",*m, Binning(nBinsMass, mMin, mMax));
    TH1* hMassModel = createHistoFromPdf(model, jPsiYield->getVal()+bkgYield->getVal(), "hMassModel", m, nBinsModel, mMin, mMax, nBinsMass);
    TH1* hMassJpsi = createHistoFromPdf(jPsiModel, jPsiYield->getVal(), "hMassJpsi", m, nBinsModel, mMin, mMax, nBinsMass);
    TH1* hMassBkg = createHistoFromPdf(bkgModel, bkgYield->getVal(), "hMassBkg", m, nBinsMass, mMin, mMax, nBinsMass);
    
    // Get reconstructed pt distributions
    int nBinsPt = 120;
    double ptMax = 1.2;
    int nModel = 1200;
    TH1* hPtBkg = dataw_bkg->createHistogram("hPtBkg",*pt, Binning(nBinsPt, 0, ptMax));
    TH1* hPtJpsi = dataw_jpsi->createHistogram("hPtJpsi",*pt, Binning(nBinsPt, 0, ptMax));
    
    // Get original pt ditribution
    RooAbsPdf* ptSignalPdf = ws->pdf("ptSignalPdf");
    RooAbsPdf* ptBkgPdf = ws->pdf("ptBkgPdf");
    
    TH1* hPtJpsiOriginal = createHistoFromPdf(ptSignalPdf, jPsiYield->getVal(), "hPtJpsiOriginal", pt, nBinsPt, ptMin, ptMax, nBinsPt);
    TH1* hPtBkgOriginal = createHistoFromPdf(ptBkgPdf, bkgYield->getVal(), "hPtBkgOriginal", pt, nBinsPt, ptMin, ptMax, nBinsPt);
    
    // Get Landau distributions
    RooRealVar* yieldTwoGamma1 = ws->var("yieldTwoGamma");
    RooRealVar* yieldTwoGamma2 = ws->var("yieldTwoGamma2");
    RooAbsPdf* ptTwoGamma1 = ws->pdf("ptTwoGamma");
    RooAbsPdf* ptTwoGamma2 = ws->pdf("ptTwoGamma2");
    
    TH1* hPtLandauOriginal = createHistoFromPdf(ptTwoGamma2, yieldTwoGamma2->getVal(), "hPtLandauOriginal", pt, nModel, ptMin, ptMax, nBinsPt);
    TH1* hPtLandauReconstructed = createHistoFromPdf(ptTwoGamma1, yieldTwoGamma1->getVal(), "hPtLandauReconstructed", pt, nModel, ptMin, ptMax, nBinsPt);
    
    // Write them
    TFile* f = new TFile("mc-toy.root", "recreate");
    // write mass distributions
    hMass->Write();
    hMassModel->Write();
    hMassJpsi->Write();
    hMassBkg->Write();
    // write pt distributions
    hPtBkgOriginal->Write();
    hPtBkg->Write();
    hPtJpsiOriginal->Write();
    hPtJpsi->Write();
    // write Landau distributions
    hPtLandauOriginal->Write();
    hPtLandauReconstructed->Write();
    f->Close();
    
}



bool CreateTrees(std::string rootfilePathMC, std::string period, bool draw = true) {
    
    
    // 1ere étape : générer des "fausses" données en utilisant les MC, les enregistrer dans des histogrammes
    
    // MC files
    TFile *fSignal = new TFile(Form("%s/AnalysisResults_%s_MC_kIncohJpsiToMu.root", rootfilePathMC.c_str(), period.c_str()),"READ");
    TTree* fSignalTree = (TTree*)fSignal->Get("MyTask/fAnaTree");
    
    TFile *fBkg = new TFile(Form("%s/AnalysisResults_%s_MC_kTwoGammaToMuLow.root", rootfilePathMC.c_str(), period.c_str()),"READ");
    TTree* fBkgTree = (TTree*)fBkg->Get("MyTask/fAnaTree");
    
    Double_t mSignal = 0, mBkg = 0, ptSignal = 0, ptBkg = 0;
    fSignalTree->SetBranchAddress("fTrkTrkM", &mSignal);
    fSignalTree->SetBranchAddress("fTrkTrkPt", &ptSignal);
    fBkgTree->SetBranchAddress("fTrkTrkM", &mBkg);
    fBkgTree->SetBranchAddress("fTrkTrkPt", &ptBkg);
    fBkgTree->SetBranchStatus("*",0);    // Select the branches to look at
    fBkgTree->SetBranchStatus("fTrkTrkM",1);    // Select the branches to look at
    fBkgTree->SetBranchStatus("fTrkTrkPt",1);    // Select the branches to look at
    
    const int nSignal = fSignalTree->GetEntries();
    const int nBkg = fBkgTree->GetEntries();
    
    std::vector<int> signalIndices = {}, bkgIndices = {};
    
    std::cout << "There are " << nSignal << " signal entries" << std::endl;
    std::cout << "There are " << nBkg << " background entries" << std::endl;
    
    
    // new file
    TFile *fAna = new TFile(Form("%s/toy_MC_%s.root", rootfilePathMC.c_str(), period.c_str()),"RECREATE");
    
    // new TTrees
    TTree* tSignal = new TTree("tSignal", "tSignal");
    tSignal->Branch("m", &mSignal, "m/D");
    tSignal->Branch("pt", &ptSignal, "pt/D");
    TTree* tBkg = new TTree("tBkg", "tBkg");
    tBkg->Branch("m", &mBkg, "m/D");
    tBkg->Branch("pt", &ptBkg, "pt/D");
    Double_t mMixed = 0, ptMixed = 0;
    TTree* tMixed = new TTree("tMixed", "tMixed");
    tMixed->Branch("m", &mMixed, "m/D");
    tMixed->Branch("pt", &ptMixed, "pt/D");
    
    // get jpsi events
    //for (int i = 0; i<nSignal; i++) {
    int nPbSig = 0;
    //std::cout << "nSignal = " << nSignal << endl;
    //return 0;
    srand (time(NULL));
    for (int i = 0; i<nEventsSignal[period]; i++) {
        int r = rand() % nSignal + 1;
        //std::cout << "r = " << r << std::endl;
        if ((int)signalIndices.size() == nSignal) {std::cout << "not enough stats! bye" << std::endl; return false;}
        if (contains(signalIndices,r)) {i--; continue;}
        signalIndices.push_back(r);
        fSignalTree->GetEntry(r);
        if (mSignal < mMin || mSignal > mMax || ptSignal < ptMin || ptSignal > ptMax) {
            nPbSig++;
            //std::cout << mSignal << std::endl;
            //std::cout << "event (signal) ignored because not in the right range of pt or m" << std::endl;
            i--;
            continue;}
        mMixed = mSignal;
        ptMixed = ptSignal;
        
        tSignal->Fill();
        tMixed->Fill();
        /*
         histM->Fill(mSignal);
         histPt->Fill(ptSignal);
         histSignalM->Fill(mSignal);
         histSignalPt->Fill(ptSignal);
         */
    }
    
    // Create a histogram that is constant or exp
    int nBins=int(1000*(mMax-mMin));
    TH1F* histFunction = new TH1F("histFunction", "histFunction", nBins, mMin, mMax);
    for (int i = 0; i<nBins; i++) {
        double value = 1;
        //double value = TMath::Exp(0.001*i);
        histFunction->SetBinContent(i, value);
    }
    
    
    int nPbBkg = 0;
    //for (int i = 0; i<nBkg; i++) {
    for (int i = 0; i<nEventsBkg[period]; i++) {
        
        int r = rand() % nBkg;
        if ((int)bkgIndices.size() == nBkg) {std::cout << "not enough stats! bye" << std::endl; return false;}
        if (contains(bkgIndices,r)) {i--; continue;}
        bkgIndices.push_back(r);
        fBkgTree->GetEntry(r);
        if (mBkg < mMin || mBkg > mMax || ptBkg < ptMin || ptBkg > ptMax) {
            nPbBkg++;
            //std::cout << mBkg << std::endl;
            //std::cout << "event (background) ignored because not in the right range of pt or m" << std::endl;
            i--;
            continue;}
        //mBkg = histFunction->GetRandom();    // when not from MC root file
        mMixed = mBkg;
        ptMixed = ptBkg;
        
        /*
         histM->Fill(mBkg);
         histPt->Fill(ptBkg);
         histBkgM->Fill(mBkg);
         histBkgPt->Fill(ptBkg);
         */
        tBkg->Fill();
        tMixed->Fill();
    }
    
    std::cout << nPbSig << " signal events ignored because not in the right range of pt or m" << std::endl;
    std::cout << nPbBkg << " background events ignored because not in the right range of pt or m" << std::endl;
    
    
    if (draw) {
        TCanvas* cv = new TCanvas("ToyMC","Toy MC",600,300) ;
        //TCanvas* cv = new TCanvas("2Dplot","2D fit",800,300) ;
        TH1F* histM = new TH1F("histM", "M distribution", 100, mMin, mMax);
        TH1F* histPt = new TH1F("histPt", "Pt distribution", 100, ptMin, ptMax);
        
        TH1F* histSignalM = new TH1F("histSignalM", "M Signal distribution", 100, mMin, mMax);
        TH1F* histSignalPt = new TH1F("histSignalPt", "Pt Signal distribution", 100, ptMin, ptMax);
        TH1F* histBkgM = new TH1F("histBkgM", "M Bkg distribution", 100, mMin, mMax);
        TH1F* histBkgPt = new TH1F("histBkgPt", "Pt Bkg distribution", 100, ptMin, ptMax);
        // m distributions
        cv->Divide(3,2) ;
        cv->cd(1);
        tSignal->Draw("m>>histSignalM");
        histSignalM->Draw();
        cv->cd(2);
        tBkg->Draw("m>>histBkgM");
        histBkgM->Draw();
        cv->cd(3);
        tMixed->Draw("m>>histM");
        histM->Draw();
        
        // pt distributions
        cv->cd(4);
        tSignal->Draw("pt>>histSignalPt");
        histSignalPt->GetXaxis()->SetRangeUser(0,2);
        histSignalPt->Draw();
        cv->cd(5);
        tBkg->Draw("pt>>histBkgPt");
        histBkgPt->GetXaxis()->SetRangeUser(0,2);
        histBkgPt->Draw();
        cv->cd(6);
        tMixed->Draw("pt>>histPt");
        histPt->GetXaxis()->SetRangeUser(0,2);
        histPt->Draw();
        
        cv->SaveAs(Form("Plots/ToyMC-%s-newtest.pdf", period.c_str()));
    }
    
    // Finally, write TTrees
    fAna->cd();
    tSignal->Write();
    tBkg->Write();
    tMixed->Write();
    fAna->Close();
    
    
    return true;
}

// 2e étape : faire un sPlot pour extraire ma distribution en pt
void ToyMC(std::string rootfilePathMC = "/Volumes/Transcend/rootFiles-pPb/MC-std", bool write = false) {
    
    bool drawPulls = false;
    //std::vector <std::string> periods = {"LHC16r", "LHC16s"};
    std::vector <std::string> periods = {"LHC16r"};
    const int nPeriods = periods.size();
    
    InitiateNumbers(nEventsSignal, nEventsBkg);
    
    for (int k = 0; k<nPeriods; k++) {
        std::string period = periods[k];
        
        TH1* hMu = new TH1F("hMu", "hMu", 20, -0.005, 0.005);
        TH1* hSigma = new TH1F("hSigma", "hSigma", 20, -0.005, 0.005);
        
        std::vector <TH1*> histVec = {};
        histVec.push_back(hMu);
        histVec.push_back(hSigma);
        
        for (int l = 0; l<1; l++)
        {
            
            
            bool newTrees = CreateTrees(rootfilePathMC, period, !write);    // false = don't draw
            if (!newTrees) {std::cout << "Trees not created! bye" << std::endl; return;}
            
            // Open the file
            TFile *fAna = new TFile(Form("%s/toy_MC_%s.root", rootfilePathMC.c_str(), period.c_str()), "READ");
            // Connect to the tree
            TTree* t = (TTree*)fAna->Get("tMixed");
            
            // New working space
            RooWorkspace* wspace = new RooWorkspace("myJpsi");
            // Import data
            RooRealVar m("m","M_{#mu#mu} (GeV/c2)", mMin, mMax);
            RooRealVar pt("pt","Dimuon p_{T} (GeV/c)", ptMin, ptMax);
            RooArgSet variables(m, pt);
            RooDataSet* data = new RooDataSet("data","data",variables,Import(*t));
            wspace->import(*data, Rename("data"));
            data->Print();
            
            // Add model for m distribution
            AddModel(wspace);
            
            wspace->Print();
            
            DoSPlot(wspace);
            
            GetPtHistMC(wspace, rootfilePathMC, period);
            MakePlots(wspace, rootfilePathMC, period, drawPulls, write, histVec);
            ExportHist(wspace);
            //cleanup
            delete wspace;
            
        }
        
        double size = 0.02;
        gStyle->SetTitleSize(size);
        gStyle->SetTitleSize(size, "XY");
        gStyle->SetLabelSize(size,"xyz");
        gStyle->SetTextFont(42);
        gStyle->SetTextSize(size);
        // Draw hists
        TCanvas* cv = new TCanvas("cv", "cv", 600, 300);
        cv->Divide(2);
        cv->cd(1);
        histVec[0]->SetTitle("#mu_{rec}-#mu_{ini}");
        histVec[0]->GetXaxis()->SetMaxDigits(4);
        histVec[0]->Draw("hist");
        cv->cd(2);
        histVec[1]->SetTitle("#sigma_{rec}-#sigma_{ini}");
        histVec[1]->GetXaxis()->SetMaxDigits(4);
        histVec[1]->Draw("hist");
        //cv->SaveAs("hists.pdf");
    }
}
