#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <sstream>
#include <cmath>

#include <TROOT.h>
#include <TStyle.h>
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TText.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TLatex.h"


using namespace std;

Double_t Square(Double_t x) {
    return x*x;
}

//____________________________________________
TGraphErrors* CreateTGraph( Int_t size, const Double_t* x, const Double_t* y, const Double_t* xErr, const Double_t* yErr )
{
    TGraphErrors* tg = new TGraphErrors();
    for( Int_t i = 0; i < size; ++i )
    {
        tg->SetPoint( i, x[i], y[i] );
        tg->SetPointError( i, xErr[i], yErr[i] );
    }
    
    return tg;
}

//____________________________________________
TGraphErrors* CreateTGraph( const std::vector<Double_t>& x, const std::vector<Double_t>& y, const std::vector<Double_t>& xErr, const std::vector<Double_t>& yErr )
{ return CreateTGraph( x.size(), &x[0], &y[0], &xErr[0], &yErr[0] ); }

//____________________________________________
TGraphErrors* CreateTGraph( const std::vector<Double_t>& x, const std::vector<Double_t>& y, Double_t xErr, Double_t yErr )
{ return CreateTGraph( x.size(), &x[0], &y[0], &xErr, &yErr ); }

bool MakeTree(bool exp, string period, Double_t rapMin, Double_t rapMax, Double_t bExcMean, Double_t bExcError) {
    
    string suffix = "";
    if (exp) suffix = "exp";
    else suffix = "powerlaw";
    
    Double_t a1, a1Err, bExc, bExcErr, bDiss, bDissErr, nDiss, nDissErr, nInc, nIncErr, pt0, pt0Err, meanGg, meanGgErr, sigmaGg, sigmaGgErr, chi2;
    Double_t jpsiExc, jpsiExcErr, jpsiDiss, jpsiDissErr, jpsiGammaPb, jpsiGammaPbErr, jpsiInc, jpsiIncErr, twoGamma, twoGammaErr, bkg, bkgErr, ratio, ratioErr;
    
    // Here build a tree
    TTree* t = new TTree("tFitResults", "Parameters resulting of the fit");
    t->Branch("a1", &a1);
    t->Branch("a1Err", &a1Err);
    t->Branch("bExc", &bExc);
    t->Branch("bExcErr", &bExcErr);
    t->Branch("bDiss", &bDiss);
    t->Branch("bDissErr", &bDissErr);
    t->Branch("nDiss", &nDiss);
    t->Branch("nDissErr", &nDissErr);
    t->Branch("nInc", &nInc);
    t->Branch("nIncErr", &nIncErr);
    t->Branch("pt0", &pt0);
    t->Branch("pt0Err", &pt0Err);
    t->Branch("meanGg", &meanGg);
    t->Branch("meanGgErr", &meanGgErr);
    t->Branch("sigmaGg", &sigmaGg);
    t->Branch("sigmaGgErr", &sigmaGgErr);
    t->Branch("chi2", &chi2);
    
    t->Branch("jpsiExc", &jpsiExc);
    t->Branch("jpsiExcErr", &jpsiExcErr);
    t->Branch("jpsiDiss", &jpsiDiss);
    t->Branch("jpsiDissErr", &jpsiDissErr);
    t->Branch("jpsiGammaPb", &jpsiGammaPb);
    t->Branch("jpsiGammaPbErr", &jpsiGammaPbErr);
    t->Branch("jpsiInc", &jpsiInc);
    t->Branch("jpsiIncErr", &jpsiIncErr);
    t->Branch("twoGamma", &twoGamma);
    t->Branch("twoGammaErr", &twoGammaErr);
    t->Branch("bkg", &bkg);
    t->Branch("bkgErr", &bkgErr);
    t->Branch("ratio", &ratio);
    t->Branch("ratioErr", &ratioErr);
    
    string filename = "output/output-" + period + "-" + suffix + Form("-%.1f-%.1f", -rapMax, -rapMin) + ".txt";
    ifstream file(filename, ios::in);
    if (file) {
        string line;
        while(getline(file, line)) {
            stringstream stream(line);
            stream >> a1 >> a1Err >> bExc >> bExcErr >> bDiss >> bDissErr >> nDiss >> nDissErr >> nInc >> nIncErr >> pt0 >> pt0Err >> meanGg >> meanGgErr >> sigmaGg >> sigmaGgErr >> chi2;
            getline(file, line);
            stringstream stream2(line);
            stream2 >> jpsiExc >> jpsiExcErr >> jpsiDiss >> jpsiDissErr >> jpsiGammaPb >> jpsiGammaPbErr >> jpsiInc >> jpsiIncErr >> twoGamma >> twoGammaErr >> bkg >> bkgErr >> ratio >> ratioErr;
            if (bExc > bExcMean-3*bExcError && bExc < bExcMean+3*bExcError ) {
                t->Fill();
            }
        }
        t->SaveAs(Form("output/%s/FitResults-%s-%.1f-%.1f.root", period.c_str(), suffix.c_str(), -rapMax, -rapMin));
        file.close();
    }
    
    else {
        cout << "Error: not possible to open " << filename << " file in reading mode" << endl;
        return false;
    }
    
    return true;
    
}


void PlotResults(bool exp, string period, Double_t rapMin, Double_t rapMax) {
    
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(0.3);
    
    
    string suffix = "";
    if (exp) suffix = "exp";
    else suffix = "powerlaw";
    
    Double_t a1, a1Err, bExc, bExcErr, bDiss, bDissErr, nDiss, nDissErr, nInc, nIncErr, pt0, pt0Err, meanGg, meanGgErr, sigmaGg, sigmaGgErr, chi2;
    Double_t jpsiExc, jpsiExcErr, jpsiDiss, jpsiDissErr, jpsiGammaPb, jpsiGammaPbErr, jpsiInc, jpsiIncErr, twoGamma, twoGammaErr, bkg, bkgErr, ratio, ratioErr;
    
    vector<Double_t> a1Vec, a1ErrVec, bExcVec, bExcErrVec, bDissVec, bDissErrVec, nDissVec, nDissErrVec, nIncVec, nIncErrVec, pt0Vec, pt0ErrVec, meanGgVec, meanGgErrVec, sigmaGgVec, sigmaGgErrVec, chi2Vec;
    vector<Double_t> jpsiExcVec, jpsiExcErrVec, jpsiDissVec, jpsiDissErrVec, jpsiGammaPbVec, jpsiGammaPbErrVec, jpsiIncVec, jpsiIncErrVec, twoGammaVec, twoGammaErrVec, bkgVec, bkgErrVec, ratioVec, ratioErrVec;
    
    TFile* f = new TFile(Form("output/%s/FitResults-%s-%.1f-%.1f.root", period.c_str(), suffix.c_str(), -rapMax, -rapMin), "READ");
    
    // Here build a tree
    TTree* t = (TTree*)f->Get("tFitResults");
    t->SetBranchAddress("a1", &a1);
    t->SetBranchAddress("a1Err", &a1Err);
    t->SetBranchAddress("bExc", &bExc);
    t->SetBranchAddress("bExcErr", &bExcErr);
    t->SetBranchAddress("bDiss", &bDiss);
    t->SetBranchAddress("bDissErr", &bDissErr);
    t->SetBranchAddress("nDiss", &nDiss);
    t->SetBranchAddress("nDissErr", &nDissErr);
    t->SetBranchAddress("nInc", &nInc);
    t->SetBranchAddress("nIncErr", &nIncErr);
    t->SetBranchAddress("pt0", &pt0);
    t->SetBranchAddress("pt0Err", &pt0Err);
    t->SetBranchAddress("meanGg", &meanGg);
    t->SetBranchAddress("meanGgErr", &meanGgErr);
    t->SetBranchAddress("sigmaGg", &sigmaGg);
    t->SetBranchAddress("sigmaGgErr", &sigmaGgErr);
    t->SetBranchAddress("chi2", &chi2);
    t->SetBranchAddress("jpsiExc", &jpsiExc);
    t->SetBranchAddress("jpsiExcErr", &jpsiExcErr);
    t->SetBranchAddress("jpsiDiss", &jpsiDiss);
    t->SetBranchAddress("jpsiDissErr", &jpsiDissErr);
    t->SetBranchAddress("jpsiGammaPb", &jpsiGammaPb);
    t->SetBranchAddress("jpsiGammaPbErr", &jpsiGammaPbErr);
    t->SetBranchAddress("jpsiInc", &jpsiInc);
    t->SetBranchAddress("jpsiIncErr", &jpsiIncErr);
    t->SetBranchAddress("twoGamma", &twoGamma);
    t->SetBranchAddress("twoGammaErr", &twoGammaErr);
    t->SetBranchAddress("bkg", &bkg);
    t->SetBranchAddress("bkgErr", &bkgErr);
    t->SetBranchAddress("ratio", &ratio);
    t->SetBranchAddress("ratioErr", &ratioErr);
    
    const int n = t->GetEntries();
    for (int i = 0; i<n; i++) {
        t->GetEntry(i);
        a1Vec.push_back(a1);
        a1ErrVec.push_back(a1Err);
        bExcVec.push_back(bExc);
        bExcErrVec.push_back(bExcErr);
        bDissVec.push_back(bDiss);
        bDissErrVec.push_back(bDissErr);
        nDissVec.push_back(nDiss);
        nDissErrVec.push_back(nDissErr);
        nIncVec.push_back(nInc);
        nIncErrVec.push_back(nIncErr);
        pt0Vec.push_back(pt0);
        pt0ErrVec.push_back(pt0Err);
        meanGgVec.push_back(meanGg);
        meanGgErrVec.push_back(meanGgErr);
        sigmaGgVec.push_back(sigmaGg);
        sigmaGgErrVec.push_back(sigmaGgErr);
        chi2Vec.push_back(chi2);
        jpsiExcVec.push_back(jpsiExc);
        jpsiExcErrVec.push_back(jpsiExcErr);
        jpsiDissVec.push_back(jpsiDiss);
        jpsiDissErrVec.push_back(jpsiDissErr);
        jpsiGammaPbVec.push_back(jpsiGammaPb);
        jpsiGammaPbErrVec.push_back(jpsiGammaPbErr);
        jpsiIncVec.push_back(jpsiInc);
        jpsiIncErrVec.push_back(jpsiIncErr);
        twoGammaVec.push_back(twoGamma);
        twoGammaErrVec.push_back(twoGammaErr);
        bkgVec.push_back(bkg);
        bkgErrVec.push_back(bkgErr);
        ratioVec.push_back(ratio);
        ratioErrVec.push_back(ratioErr);
        
        if (nDiss == 0) exp = true;
    }
    
    // Draw in a canvas
    // First parameters
    TCanvas* cv = new TCanvas("cv", "cv", 500, 800);
    cv->Divide(2, 4);
    TGraphErrors* gr1 = CreateTGraph(bExcVec, a1Vec, bExcErrVec, a1ErrVec);
    gr1->SetTitle("a1 = f(b_{exc})");
    TGraphErrors* gr2 = CreateTGraph(bExcVec, bDissVec, bExcErrVec, bDissErrVec);
    gr2->SetTitle("b_{diss} = f(b_{exc})");
    TGraphErrors* gr3 = CreateTGraph(bExcVec, nDissVec, bExcErrVec, nDissErrVec);
    gr3->SetTitle("n_{diss} = f(b_{exc})");
    TGraphErrors* gr4 = CreateTGraph(bExcVec, pt0Vec, bExcErrVec, pt0ErrVec);
    gr4->SetTitle("p_{T0} = f(b_{exc})");
    TGraphErrors* gr5 = CreateTGraph(bExcVec, nIncVec, bExcErrVec, nIncErrVec);
    gr5->SetTitle("n_{inc} = f(b_{exc})");
    TGraphErrors* gr6 = CreateTGraph(bExcVec, meanGgVec, bExcErrVec, meanGgErrVec);
    gr6->SetTitle("#mu_{Landau} = f(b_{exc})");
    TGraphErrors* gr7 = CreateTGraph(bExcVec, sigmaGgVec, bExcErrVec, sigmaGgErrVec);
    gr7->SetTitle("#sigma_{Landau} = f(b_{exc})");
    TGraphErrors* gr8 = CreateTGraph(bExcVec, chi2Vec, bExcErrVec, bExcErrVec);
    gr8->SetTitle("#chi^{2}= f(b_{exc})");
    
    
    cv->cd(1);
    gr1->Draw("AP");
    cv->cd(2);
    gr2->Draw("AP");
    cv->cd(3);
    gr3->Draw("AP");
    cv->cd(4);
    gr4->Draw("AP");
    cv->cd(5);
    gr5->Draw("AP");
    cv->cd(6);
    gr6->SetMinimum(0.02);
    gr6->SetMaximum(0.1);
    gr6->Draw("AP");
    cv->cd(7);
    gr7->SetMinimum(0);
    gr7->SetMaximum(0.2);
    gr7->Draw("AP");
    cv->cd(8);
    gr8->Draw("AP");
    
    cv->SaveAs(Form("Plots/%s/Param-%s-%.1f-%.1f.pdf", period.c_str(), suffix.c_str(), -rapMax, -rapMin));
    
    //Then yields
    TCanvas* cv2 = new TCanvas("cv2", "cv2", 400, 800);
    cv2->Divide(2, 4);
    TGraphErrors* grr1 = CreateTGraph(bExcVec, jpsiExcVec, bExcErrVec, jpsiExcErrVec);
    grr1->SetTitle("# exc J/#Psi = f(b_{exc})");
    grr1->GetYaxis()->SetRange(0,3000);
    TGraphErrors* grr2 = CreateTGraph(bExcVec, jpsiDissVec, bExcErrVec, jpsiDissErrVec);
    grr2->SetTitle("# diss J/#Psi = f(b_{exc})");
    TGraphErrors* grr3 = CreateTGraph(bExcVec, jpsiGammaPbVec, bExcErrVec, jpsiGammaPbErrVec);
    grr3->SetTitle("# #gamma-Pb J/#Psi = f(b_{exc})");
    TGraphErrors* grr4 = CreateTGraph(bExcVec, jpsiIncVec, bExcErrVec, jpsiIncErrVec);
    grr4->SetTitle("# inc J/#Psi = f(b_{exc})");
    TGraphErrors* grr5 = CreateTGraph(bExcVec, twoGammaVec, bExcErrVec, twoGammaVec);
    grr5->SetTitle("# #gamma#gamma = f(b_{exc})");
    TGraphErrors* grr6 = CreateTGraph(bExcVec, bkgVec, bExcErrVec, bkgErrVec);
    grr6->SetTitle("# bkg = f(b_{exc})");
    TGraphErrors* grr7 = CreateTGraph(bExcVec, ratioVec, bExcErrVec, ratioErrVec);
    grr7->SetTitle("N_{diss}/N_{exc} = f(b_{exc})");
    
    
    cv2->cd(1);
    grr1->Draw("AP");
    cv2->cd(2);
    grr2->Draw("AP");
    cv2->cd(3);
    grr3->Draw("AP");
    cv2->cd(4);
    grr4->Draw("AP");
    cv2->cd(5);
    grr5->Draw("AP");
    cv2->cd(6);
    grr6->Draw("AP");
    cv2->cd(7);
    grr7->Draw("AP");
    
    cv2->SaveAs(Form("Plots/%s/Yields-%s-%.1f-%.1f.pdf", period.c_str(), suffix.c_str(), -rapMax, -rapMin));
    f->Close();
}

void SystematicsComputation(bool exp, string period, Double_t rapMin, Double_t rapMax, double bExcMean, double bExcSigma) {
    
    
    string suffix = "";
    if (exp) suffix = "exp";
    else suffix = "powerlaw";
    
    Double_t a1, a1Err, bExc, bExcErr, bDiss, bDissErr, nDiss, nDissErr, nInc, nIncErr, pt0, pt0Err, chi2;
    Double_t jpsiExc, jpsiExcErr, jpsiDiss, jpsiDissErr, jpsiGammaPb, jpsiGammaPbErr, jpsiInc, jpsiIncErr, twoGamma, twoGammaErr, bkg, bkgErr, ratio, ratioErr;
    
    vector<Double_t> a1Vec, a1ErrVec, bExcVec, bExcErrVec, bDissVec, bDissErrVec, nDissVec, nDissErrVec, nIncVec, nIncErrVec, pt0Vec, pt0ErrVec, chi2Vec;
    vector<Double_t> jpsiExcVec, jpsiExcErrVec, jpsiDissVec, jpsiDissErrVec, jpsiGammaPbVec, jpsiGammaPbErrVec, jpsiIncVec, jpsiIncErrVec, twoGammaVec, twoGammaErrVec, bkgVec, bkgErrVec, ratioVec, ratioErrVec;
    
    TFile* f = new TFile(Form("output/%s/FitResults-%s-%.1f-%.1f.root", period.c_str(), suffix.c_str(), -rapMax, -rapMin), "READ");
    
    // Here build a tree
    TTree* t = (TTree*)f->Get("tFitResults");
    t->SetBranchAddress("a1", &a1);
    t->SetBranchAddress("a1Err", &a1Err);
    t->SetBranchAddress("bExc", &bExc);
    t->SetBranchAddress("bExcErr", &bExcErr);
    t->SetBranchAddress("bDiss", &bDiss);
    t->SetBranchAddress("bDissErr", &bDissErr);
    t->SetBranchAddress("nDiss", &nDiss);
    t->SetBranchAddress("nDissErr", &nDissErr);
    t->SetBranchAddress("nInc", &nInc);
    t->SetBranchAddress("nIncErr", &nIncErr);
    t->SetBranchAddress("pt0", &pt0);
    t->SetBranchAddress("chi2", &chi2);
    t->SetBranchAddress("pt0Err", &pt0Err);
    t->SetBranchAddress("jpsiExc", &jpsiExc);
    t->SetBranchAddress("jpsiExcErr", &jpsiExcErr);
    t->SetBranchAddress("jpsiDiss", &jpsiDiss);
    t->SetBranchAddress("jpsiDissErr", &jpsiDissErr);
    t->SetBranchAddress("jpsiGammaPb", &jpsiGammaPb);
    t->SetBranchAddress("jpsiGammaPbErr", &jpsiGammaPbErr);
    t->SetBranchAddress("jpsiInc", &jpsiInc);
    t->SetBranchAddress("jpsiIncErr", &jpsiIncErr);
    t->SetBranchAddress("twoGamma", &twoGamma);
    t->SetBranchAddress("twoGammaErr", &twoGammaErr);
    t->SetBranchAddress("bkg", &bkg);
    t->SetBranchAddress("bkgErr", &bkgErr);
    t->SetBranchAddress("ratio", &ratio);
    t->SetBranchAddress("ratioErr", &ratioErr);
    
    const int n = t->GetEntries();
    for (int i = 0; i<n; i++) {
        t->GetEntry(i);
        a1Vec.push_back(a1);
        a1ErrVec.push_back(a1Err);
        bExcVec.push_back(bExc);
        bExcErrVec.push_back(bExcErr);
        bDissVec.push_back(bDiss);
        bDissErrVec.push_back(bDissErr);
        nDissVec.push_back(nDiss);
        nDissErrVec.push_back(nDissErr);
        nIncVec.push_back(nInc);
        nIncErrVec.push_back(nIncErr);
        pt0Vec.push_back(pt0);
        chi2Vec.push_back(chi2);
        pt0ErrVec.push_back(pt0Err);
        jpsiExcVec.push_back(jpsiExc);
        jpsiExcErrVec.push_back(jpsiExcErr);
        jpsiDissVec.push_back(jpsiDiss);
        jpsiDissErrVec.push_back(jpsiDissErr);
        jpsiGammaPbVec.push_back(jpsiGammaPb);
        jpsiGammaPbErrVec.push_back(jpsiGammaPbErr);
        jpsiIncVec.push_back(jpsiInc);
        jpsiIncErrVec.push_back(jpsiIncErr);
        twoGammaVec.push_back(twoGamma);
        twoGammaErrVec.push_back(twoGammaErr);
        bkgVec.push_back(bkg);
        bkgErrVec.push_back(bkgErr);
        ratioVec.push_back(ratio);
        ratioErrVec.push_back(ratioErr);
        
        if (nDiss == 0) exp = true;
    }
    
    // First create the weights on each measurement
    vector<double> weights = {};
    double sum = 0;
    for (int i = 0; i<n; i++) {
        double value = TMath::Exp(-Square(bExcMean-bExcVec[i])/(2*Square(bExcSigma)));
        weights.push_back(value);
        sum += value;
    }
    for (int i = 0; i<n; i++) {
        weights[i] /= sum;
    }
    
    // Second compute the mean value
    Double_t jpsiExcMean = 0;
    Double_t jpsiDissMean = 0;
    Double_t ggMean = 0;
    Double_t ratioMean = 0;
    for (int i = 0; i<n; i++) {
        jpsiExcMean += weights[i]* jpsiExcVec[i];
        jpsiDissMean += weights[i]* jpsiDissVec[i];
        ggMean += weights[i]* twoGammaVec[i];
        ratioMean += weights[i]* ratioVec[i];
    }
    
    // Third compute the uncertainty
    Double_t jpsiExcError = 0;
    Double_t jpsiDissError = 0;
    Double_t ggError = 0;
    Double_t ratioError = 0;
    for (int i = 0; i<n; i++) {
        jpsiExcError += weights[i]* Square(jpsiExcVec[i]-jpsiExcMean);
        jpsiDissError += weights[i]* Square(jpsiDissVec[i]-jpsiDissMean);
        ggError += weights[i]* Square(twoGammaVec[i]-ggMean);
        ratioError += weights[i]* Square(ratioVec[i]-ratioMean);
    }
    
    jpsiExcError = TMath::Sqrt((double)n/(n-1)*jpsiExcError);
    jpsiDissError = TMath::Sqrt((double)n/(n-1)*jpsiDissError);
    ggError = TMath::Sqrt((double)n/(n-1)*ggError);
    ratioError = TMath::Sqrt((double)n/(n-1)*ratioError);
    
    cout << endl << endl;
    cout << "Number of exclusive J/Psi = " << jpsiExcMean << " #pm " << jpsiExcError << endl;
    cout << "Number of dissociative J/Psi = " << jpsiDissMean << " #pm " << jpsiDissError << endl;
    cout << "Number of gamma-gamma = " << ggMean << " #pm " << ggError << endl;
    cout << "Ratio diss / exc = " << ratioMean << " #pm " << ratioError << endl;
    cout << endl << endl;
    
    vector<double> nullVec = {};
    for (int j = 0; j<(int)weights.size(); j++) {nullVec.push_back(0);}
    
    TCanvas* cv = new TCanvas("cv", "cv", 400, 600);
    cv->Divide(2,3);
    TGraphErrors* gr1 = CreateTGraph(bExcVec, weights, bExcErrVec, nullVec);
    gr1->SetTitle("weights = f(b_{exc})");
    gr1->GetXaxis()->SetTitle("b_{exc}");
    gr1->GetYaxis()->SetTitle("weights");
    TGraphErrors* gr2 = CreateTGraph(bExcVec, jpsiExcVec, bExcErrVec, jpsiExcErrVec);
    gr2->SetTitle("N_{exc} = f(b_{exc})");
    gr2->GetYaxis()->SetTitle("N_{exc}");
    gr2->GetXaxis()->SetTitle("b_{exc}");
    TGraphErrors* gr3 = CreateTGraph(bExcVec, jpsiDissVec, bExcErrVec, jpsiDissErrVec);
    gr3->SetTitle("N_{diss} = f(b_{exc})");
    gr3->GetYaxis()->SetTitle("N_{diss}");
    gr3->GetXaxis()->SetTitle("b_{exc}");
    TGraphErrors* gr4 = CreateTGraph(bExcVec, twoGammaVec, bExcErrVec, twoGammaErrVec);
    gr4->SetTitle("N_{#gamma#gamma} = f(b_{exc})");
    gr4->GetYaxis()->SetTitle("N_{#gamma#gamma}");
    gr4->GetXaxis()->SetTitle("b_{exc}");
    TGraphErrors* gr6 = CreateTGraph(bExcVec, ratioVec, bExcErrVec, ratioErrVec);
    gr4->SetTitle("N_{diss}/N_{exc} = f(b_{exc})");
    gr4->GetYaxis()->SetTitle("N_{diss}/N_{exc}");
    gr4->GetXaxis()->SetTitle("b_{exc}");
    
    cv->cd(1);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.15);
    gr1->Draw();
    cv->cd(2);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.15);
    gr2->Draw();
    cv->cd(3);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.15);
    gr3->Draw();
    cv->cd(4);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.15);
    gr4->Draw();
    
    cv->cd(6);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.15);
    gr6->Draw();
    
    
    cv->cd(5);
    TLatex* tExc = new TLatex(0.1, 0.8, Form("#bf{<N_{exc}> = %.1f #pm %.1f}", jpsiExcMean, jpsiExcError));
    TLatex* tDiss = new TLatex(0.1, 0.7, Form("#bf{<N_{diss}> = %.1f #pm %.1f}", jpsiDissMean, jpsiDissError));
    TLatex* tTwoGamma = new TLatex(0.1, 0.6, Form("#bf{<N_{#gamma#gamma}> = %.1f #pm %.1f}", ggMean, ggError));
    TLatex* tRatio = new TLatex(0.1, 0.5, Form("#bf{<N_{diss}/N_{exc}> = %.2f #pm %.2f}", ratioMean, ratioError));
    tExc->Draw();
    tDiss->Draw();
    tTwoGamma->Draw();
    tRatio->Draw();
    
    cv->SaveAs(Form("Plots/%s/SystError-%s-%.1f-%.1f.pdf", period.c_str(), suffix.c_str(), -rapMax, -rapMin));
    f->Close();
    
}


void FitResultAnalysis(bool exp = false, Int_t rapMode = 1) {
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetTitleSize(.05);
    gStyle->SetLabelSize(.05, "XY");
    
    vector<string> periods = {"LHC16r", "LHC16s"};
    //vector<string> periods = {"LHC16s"};
    vector<double> bExcMeanVec, bExcErrorVec = {};
    
    
    
    Double_t rapMin = 0, rapMax = 0;
    if (rapMode == 1) {
        rapMin = -4.;
        rapMax = -2.5;
        bExcMeanVec.push_back(3.62);
        bExcErrorVec.push_back(0.14);
        bExcMeanVec.push_back(5.08);
        bExcErrorVec.push_back(0.24);
    }
    else if (rapMode == 2) {
        rapMin = -4.;
        rapMax = -3.25;
        bExcMeanVec.push_back(3.38);
        bExcErrorVec.push_back(0.17);
        bExcMeanVec.push_back(4.51);
        bExcErrorVec.push_back(0.26);
    }
    else if (rapMode == 3) {
        rapMin = -3.25;
        rapMax = -2.5;
        bExcMeanVec.push_back(3.86);
        bExcErrorVec.push_back(0.20);
        bExcMeanVec.push_back(5.36);
        bExcErrorVec.push_back(0.30);
    }
    else {
        cout << "Unknown rapidity mode, terminating" << endl;
        return;
    }
    
    /*
     vector<double> bExcMeanVec = {4.16, 5.93};
     vector<double> bExcErrorVec = {0.17, 0.36};
     
     vector<double> bExcMeanVec = {4.16, 5.52};
     vector<double> bExcErrorVec = {0.26, 0.52};
     
     vector<double> bExcMeanVec = {3.86, 6.31};
     vector<double> bExcErrorVec = {0.19, 0.40};
     */
    
    const int nPeriods = periods.size();
    for (int k = 0; k<nPeriods; k++) {
        string period = periods[k];
        
        bool init = MakeTree(exp, period, rapMin, rapMax, bExcMeanVec[k], bExcErrorVec[k]);
        if (!init) return;
        
        PlotResults(exp, period, rapMin, rapMax);
        
        SystematicsComputation(exp, period, rapMin, rapMax, bExcMeanVec[k], bExcErrorVec[k]);
    }
    
}
