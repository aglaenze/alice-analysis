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

Bool_t contains(vector<Double_t> vec, Double_t num) {
    for (int i = 0; i < (int) vec.size(); i++) {
        if (num > 0.99*vec[i] && num < 1.01*vec[i]) return true;
    }
    return false;
}

//____________________________________________
TGraphErrors* CreateTGraph( const std::vector<Double_t>& x, const std::vector<Double_t>& y, const std::vector<Double_t>& xErr, const std::vector<Double_t>& yErr )
{ return CreateTGraph( x.size(), &x[0], &y[0], &xErr[0], &yErr[0] ); }

//____________________________________________
TGraphErrors* CreateTGraph( const std::vector<Double_t>& x, const std::vector<Double_t>& y, Double_t xErr, Double_t yErr )
{ return CreateTGraph( x.size(), &x[0], &y[0], &xErr, &yErr ); }

bool MakeTree(string period, Double_t rapMin, Double_t rapMax, Double_t massMin, Double_t massMax, Double_t mu, Double_t muErr, Double_t sigma, Double_t sigmaErr, Int_t& nMu, Int_t& nSigma) {
    
    Double_t meanGg, meanGgErr, sigmaGg, sigmaGgErr;
    Double_t twoGamma, twoGammaErr, bkg, bkgErr;
    
    // Here build a tree
    TTree* t = new TTree("tFitResults", "Parameters resulting of the fit");
    t->Branch("meanGg", &meanGg);
    t->Branch("meanGgErr", &meanGgErr);
    t->Branch("sigmaGg", &sigmaGg);
    t->Branch("sigmaGgErr", &sigmaGgErr);
    
    t->Branch("twoGamma", &twoGamma);
    t->Branch("twoGammaErr", &twoGammaErr);
    t->Branch("bkg", &bkg);
    t->Branch("bkgErr", &bkgErr);
    
    vector<Double_t> muList = {}, sigmaList = {};
    string filename = Form("output/%s/%.1f-%.1f/output-gg-%.1f-%.1f.txt", period.c_str(), massMin, massMax, -rapMax, -rapMin);
    ifstream file(filename, ios::in);
    if (file) {
        string line;
        while(getline(file, line)) {
            stringstream stream(line);
            stream >> meanGg >> meanGgErr >> sigmaGg >> sigmaGgErr;
            getline(file, line);
            stringstream stream2(line);
            stream2 >> twoGamma >> twoGammaErr >> bkg >> bkgErr;
            if ((meanGg > mu-3*muErr && meanGg < mu+3*muErr) && (sigmaGg > sigma-3*sigmaErr && sigmaGg < sigma+3*sigmaErr) ) {
                t->Fill();
                if (!contains(muList, meanGg)) {muList.push_back(meanGg);}
                if (!contains(sigmaList, sigmaGg)) {sigmaList.push_back(sigmaGg);}
            }
        }
        t->SaveAs(Form("output/%s/%.1f-%.1f/FitResults-%.1f-%.1f.root", period.c_str(), massMin, massMax, -rapMax, -rapMin));
        file.close();
        nMu = (int)muList.size();
        nSigma = (int)sigmaList.size();
    }
    
    else {
        cout << "Error: not possible to open " << filename << " file in reading mode" << endl;
        return false;
    }
    
    return true;
    
}

/*
void PlotResults(string period, Double_t rapMin, Double_t rapMax, Double_t massMin, Double_t massMax) {
    
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(0.3);
    
    Double_t meanGg, meanGgErr, sigmaGg, sigmaGgErr;
    Double_t twoGamma, twoGammaErr, bkg, bkgErr;
    
    vector<Double_t> meanGgVec, meanGgErrVec, sigmaGgVec, sigmaGgErrVec;
    vector<Double_t> twoGammaVec, twoGammaErrVec, bkgVec, bkgErrVec;
    
    TFile* f = new TFile(Form("output/%s/%.1f-%.1f/FitResults-%.1f-%.1f.root", period.c_str(), massMin, massMax, -rapMax, -rapMin), "READ");
    
    // Here build a tree
    TTree* t = (TTree*)f->Get("tFitResults");
    t->SetBranchAddress("meanGg", &meanGg);
    t->SetBranchAddress("meanGgErr", &meanGgErr);
    t->SetBranchAddress("sigmaGg", &sigmaGg);
    t->SetBranchAddress("sigmaGgErr", &sigmaGgErr);
    
    t->SetBranchAddress("twoGamma", &twoGamma);
    t->SetBranchAddress("twoGammaErr", &twoGammaErr);
    t->SetBranchAddress("bkg", &bkg);
    t->SetBranchAddress("bkgErr", &bkgErr);
    
    const int n = t->GetEntries();
    for (int i = 0; i<n; i++) {
        t->GetEntry(i);
        meanGgVec.push_back(meanGg);
        meanGgErrVec.push_back(meanGgErr);
        sigmaGgVec.push_back(sigmaGg);
        sigmaGgErrVec.push_back(sigmaGgErr);
        
        twoGammaVec.push_back(twoGamma);
        twoGammaErrVec.push_back(twoGammaErr);
        bkgVec.push_back(bkg);
        bkgErrVec.push_back(bkgErr);
        
    }
    
    
    //Then yields
    TCanvas* cv2 = new TCanvas("cv2", "cv2", 600, 300);
    cv2->Divide(2);
    TGraphErrors* grr1 = CreateTGraph(meanGgVec, twoGammaVec, meanGgErrVec, twoGammaErrVec);
    grr1->SetTitle("# #gamma#gamma = f(#mu)");
    TGraphErrors* grr2 = CreateTGraph(meanGgVec, bkgVec, meanGgErrVec, bkgErrVec);
    grr2->SetTitle("# bkg = f(#mu)");
    
    
    cv2->cd(1);
    grr1->Draw("AP");
    cv2->cd(2);
    grr2->Draw("AP");
    
    
    cv2->SaveAs(Form("Plots/%s/%.1f-%.1f/Yields-%.1f-%.1f.pdf", period.c_str(), massMin, massMax, -rapMax, -rapMin));
    f->Close();
}
 */

void SystematicsComputation(string period, Double_t rapMin, Double_t rapMax, Double_t massMin, Double_t massMax, double ggMean, double ggMeanError, double ggSigma, double ggSigmaError, int nMu, int nSigma, double r) {
    
    
    Double_t meanGg, meanGgErr, sigmaGg, sigmaGgErr;
    Double_t twoGamma, twoGammaErr, bkg, bkgErr;
    
    vector<vector<Double_t>> meanGgVec, meanGgErrVec, sigmaGgVec, sigmaGgErrVec;
    vector<vector<Double_t>> twoGammaVec, twoGammaErrVec, bkgVec, bkgErrVec;
    
    TFile* f = new TFile(Form("output/%s/%.1f-%.1f/FitResults-%.1f-%.1f.root", period.c_str(), massMin, massMax, -rapMax, -rapMin), "READ");
    
    // Here build a tree
    TTree* t = (TTree*)f->Get("tFitResults");
    
    t->SetBranchAddress("meanGg", &meanGg);
    t->SetBranchAddress("meanGgErr", &meanGgErr);
    t->SetBranchAddress("sigmaGg", &sigmaGg);
    t->SetBranchAddress("sigmaGgErr", &sigmaGgErr);
    
    t->SetBranchAddress("twoGamma", &twoGamma);
    t->SetBranchAddress("twoGammaErr", &twoGammaErr);
    t->SetBranchAddress("bkg", &bkg);
    t->SetBranchAddress("bkgErr", &bkgErr);

    for (int i = 0; i<nMu; i++) {
        meanGgVec.push_back({});
        meanGgErrVec.push_back({});
        sigmaGgVec.push_back({});
        sigmaGgErrVec.push_back({});
        
        twoGammaVec.push_back({});
        twoGammaErrVec.push_back({});
        bkgVec.push_back({});
        bkgErrVec.push_back({});
        for (int j = 0; j<nSigma; j++) {
            t->GetEntry(nMu*i+j);
            meanGgVec[i].push_back(meanGg);
            meanGgErrVec[i].push_back(meanGgErr);
            sigmaGgVec[i].push_back(sigmaGg);
            sigmaGgErrVec[i].push_back(sigmaGgErr);
            
            twoGammaVec[i].push_back(twoGamma);
            twoGammaErrVec[i].push_back(twoGammaErr);
            bkgVec[i].push_back(bkg);
            bkgErrVec[i].push_back(bkgErr);
            
        }
    }
    
    /*
    cout << endl << endl << "nMu = " << nMu << endl << endl;
    cout << endl << endl << "nSigma = " << nSigma << endl << endl;
    for (int i = 0; i<nMu; i++) {
        for (int j = 0; j < nSigma; j ++ ) {
            cout << "mu = " << meanGgVec[i][j] << endl;
            cout << "sigma = " << sigmaGgVec[i][j] << endl;
        }
        cout << endl;
    }
    */
    cout << endl << endl;
    double muStep = (meanGgVec[nMu-1][nSigma-1]-meanGgVec[0][0])/nMu;
    //double muStep = ((*max_element(meanGgVec.begin(), meanGgVec.end())-meanGgVec[0][0]))/nMu;
    double muMin = meanGgVec[0][0]-muStep/2;
    double muMax = meanGgVec[nMu-1][nSigma-1]+muStep/2;
    double sigmaStep = (sigmaGgVec[nMu-1][nSigma-1]-sigmaGgVec[0][0])/nSigma;
    double sigmaMin = sigmaGgVec[0][0]-sigmaStep/2;
    double sigmaMax = sigmaGgVec[nMu-1][nSigma-1]+sigmaStep/2;
    
    cout << "mu step = " << muStep << endl;
    cout << "sigma step = " << sigmaStep << endl;
    
    // First create the weights on each measurement
    vector<vector<double>> weights = {};
    double sum = 0;
    for (int i = 0; i<nMu; i++) {
        weights.push_back({});
        for (int j = 0; j < nSigma; j ++ ) {
            double term1 = Square((ggMean-meanGgVec[i][j])/ggMeanError);
            double term2 = Square((ggSigma-sigmaGgVec[i][j])/ggSigmaError);
            double term3 = (ggMean-meanGgVec[i][j])*(ggSigma-sigmaGgVec[i][j])/(ggMeanError*ggSigmaError);
            double value = TMath::Exp(-1./2*(term1 + term2 - 2*r *term3)/(1-r*r));
            double norm = 2*TMath::Pi()*ggMeanError*ggSigmaError*TMath::Sqrt(1-r*r);
//            cout << endl << "term1 = " << term1 << endl;
//            cout << "term2 = " << term2 << endl;
//            cout << "term3 = " << term3 << endl;
//            cout << "norm = " << norm << endl;
//            cout << "value = " << value << endl;
//            cout << "value/norm = " << value/norm << endl << endl;
            weights[i].push_back(value/norm);
            sum += value/norm;
        }
    }
    for (int i = 0; i< nMu; i++) {
        for (int j = 0; j < nSigma; j ++ ) {
            weights[i][j] /= sum;
        }
    }

    // Second compute the mean value
    Double_t weightedGgNum = 0;
    Double_t weightMin = 1, weightMax = 0;
    for (int i = 0; i< nMu; i++) {
        for (int j = 0; j < nSigma; j ++ ) {
            weightedGgNum += weights[i][j]*twoGammaVec[i][j];
            if (weightMin > weights[i][j]) weightMin = weights[i][j];
            if (weightMax < weights[i][j]) weightMax = weights[i][j];
        }
    }
    
    // Third compute the uncertainty
    Double_t weightedGgNumError = 0;
    for (int i = 0; i < nMu; i++) {
        for (int j = 0; j < nSigma; j ++ ) {
            weightedGgNumError += weights[i][j]* Square(twoGammaVec[i][j]-weightedGgNum);
        }
    }
    
    weightedGgNumError = TMath::Sqrt((double)nMu*nSigma/(nMu*nSigma-1)*weightedGgNumError);
    
    cout << endl << endl;
    cout << "Number of gamma-gamma = " << weightedGgNum << " #pm " << weightedGgNumError << endl;
    cout << endl << endl;
    
    TCanvas* cv = new TCanvas("cv", "cv", 600, 600);
    cv->Divide(2,2);
    /*
     vector<double> nullVec = {};
     for (int j = 0; j<(int)weights.size(); j++) {nullVec.push_back(0);}
     TGraphErrors* gr1 = CreateTGraph(meanGgVec, weights, meanGgErrVec, nullVec);
     gr1->SetTitle("weights = f(#mu)");
     gr1->GetYaxis()->SetTitle("weights");
     gr1->GetHistogram()->SetMinimum(0);
     gr1->GetXaxis()->SetTitle("#mu");
     TGraphErrors* gr2 = CreateTGraph(sigmaGgVec, weights, sigmaGgErrVec, nullVec);
     gr2->SetTitle("weights = f(#sigma)");
     gr2->GetYaxis()->SetTitle("weights");
     gr2->GetHistogram()->SetMinimum(0);
     gr2->GetXaxis()->SetTitle("#sigma");
     
     TGraphErrors* gr3 = CreateTGraph(meanGgVec, twoGammaVec, meanGgErrVec, twoGammaErrVec);
     gr3->SetTitle("N_{#gamma#gamma} = f(#mu)");
     gr3->GetYaxis()->SetTitle("N_{#gamma#gamma}");
     gr3->GetXaxis()->SetTitle("#mu");
     */

    TH2F* hWeights = new TH2F("hWeights", "weights", nMu, muMin, muMax, nSigma, sigmaMin, sigmaMax);
    for (int i = 0; i < nMu; i++) {
        for (int j = 0; j < nSigma; j ++ ) {
            hWeights->Fill(meanGgVec[i][j], sigmaGgVec[i][j], weights[i][j]);
            //hWeights->SetBinContent(i+1, j+1, weights[i][j]);
            //hWeights->Fill(meanGgVec[i][j], sigmaGgVec[i][j], 1);
        }
    }
    hWeights->GetXaxis()->SetTitle("#mu");
    hWeights->GetYaxis()->SetTitle("#sigma");
    
    TH2F* hNumGg = new TH2F("hNumGg", "N_{#gamma#gamma}", nMu, muMin, muMax, nSigma, sigmaMin, sigmaMax);
    Double_t numMin = 1000, numMax = 0;
    for (int i = 0; i < nMu; i++) {
        for (int j = 0; j < nSigma; j ++ ) {
            hNumGg->Fill(meanGgVec[i][j], sigmaGgVec[i][j], twoGammaVec[i][j]);
            //hNumGg->SetBinContent(i+1, j+1, twoGammaVec[i][j]);
            if (numMin > twoGammaVec[i][j]) numMin = twoGammaVec[i][j];
            if (numMax < twoGammaVec[i][j]) numMax = twoGammaVec[i][j];
        }
    }
    hNumGg->GetXaxis()->SetTitle("#mu");
    hNumGg->GetYaxis()->SetTitle("#sigma");
    
    cv->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    //gr1->Draw("APC");
    //hWeights->GetZaxis()->SetColzRange(hWeights->GetMinimum(),hWeights->GetMaximum());
    hWeights->GetZaxis()->SetRangeUser(weightMin, weightMax);
    hWeights->GetXaxis()->SetMaxDigits(3);
    hWeights->Draw("colz");
    //hWeights->Draw("");
    cout << "weightMin = " << weightMin << endl;
    cout << "weightMax = " << weightMax << endl;
    
    cv->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    //gr2->Draw("APC");
    hNumGg->GetZaxis()->SetRangeUser(numMin, numMax);
    hNumGg->GetXaxis()->SetMaxDigits(3);
    hNumGg->Draw("colz");
    
    cv->cd(4);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    //gr3->Draw();
    
    cv->cd(3);
    string percent = "%";
    TLatex* tTwoGamma = new TLatex(0.1, 0.6, Form("#bf{<N_{#gamma#gamma}> = %.1f #pm %.1f}", weightedGgNum, weightedGgNumError));
    tTwoGamma->Draw();
    TLatex* tRelUncertainty = new TLatex(0.1, 0.45, Form("#bf{#frac{#Delta <N_{#gamma#gamma}>}{<N_{#gamma#gamma}>} = %.1f %s}", weightedGgNumError/weightedGgNum*100, percent.c_str()));
    tTwoGamma->Draw();
    tRelUncertainty->Draw();
    
    cv->SaveAs(Form("Plots/%s/%.1f-%.1f/SystError-%.1f-%.1f.pdf", period.c_str(), massMin, massMax, -rapMax, -rapMin));
    f->Close();
    
}


void FitResultAnalysis(Int_t rapMode = 1, Int_t massMode = 1) {
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleSize(.05);
    gStyle->SetTitleXSize(.04);
    gStyle->SetTitleYSize(.04);
    gStyle->SetLabelSize(.04, "XY");
    
    //vector<string> periods = {"LHC16r", "LHC16s"};
    vector<string> periods = {"LHC16r"};
    vector<double> meanVec, meanErrorVec = {};
    vector<double> sigmaVec, sigmaErrorVec = {};
    vector<double> corrVec = {};
    
    Double_t massMin = 0, massMax = 0;
    Double_t rapMin = 0, rapMax = 0;
    
    if (massMode == 1) {
        massMin = 1.0;
        massMax = 1.5;
        
        if (rapMode == 1) {
            rapMin = -4.;
            rapMax = -2.5;
            meanVec.push_back(0.067);
            meanErrorVec.push_back(0.003);
            sigmaVec.push_back(0.030);
            sigmaErrorVec.push_back(0.002);
            corrVec.push_back(0.52);
        }
        else if (rapMode == 2) {
            rapMin = -4.;
            rapMax = -3.25;
            meanVec.push_back(0.070);
            meanErrorVec.push_back(0.003);
            sigmaVec.push_back(0.032);
            sigmaErrorVec.push_back(0.002);
            corrVec.push_back(0.51);
        }
        else if (rapMode == 3) {
            rapMin = -3.25;
            rapMax = -2.5;
            meanVec.push_back(0.054);
            meanErrorVec.push_back(0.005);
            sigmaVec.push_back(0.023);
            sigmaErrorVec.push_back(0.003);
            corrVec.push_back(0.58);
        }
        else {
            cout << "Unknown rapidity mode, terminating" << endl;
            return;
        }
    }
    else if (massMode == 2) {
        massMin = 1.5;
        massMax = 2.0;
        if (rapMode == 1) {
            rapMin = -4.;
            rapMax = -2.5;
            meanVec.push_back(0.077);
            meanErrorVec.push_back(0.004);
            sigmaVec.push_back(0.035);
            sigmaErrorVec.push_back(0.003);
            corrVec.push_back(0.52);
        }
        else if (rapMode == 2) {
            rapMin = -4.;
            rapMax = -3.25;
            meanVec.push_back(0.076);
            meanErrorVec.push_back(0.005);
            sigmaVec.push_back(0.034);
            sigmaErrorVec.push_back(0.003);
            corrVec.push_back(0.56);
        }
        else if (rapMode == 3) {
            rapMin = -3.25;
            rapMax = -2.5;
            meanVec.push_back(0.078);
            meanErrorVec.push_back(0.008);
            sigmaVec.push_back(0.039);
            sigmaErrorVec.push_back(0.006);
            corrVec.push_back(0.44);
        }
        else {
            cout << "Unknown rapidity mode, terminating" << endl;
            return;
        }
    }
    else if (massMode == 3) {
        massMin = 2.0;
        massMax = 2.5;
        
        if (rapMode == 1) {
            rapMin = -4.;
            rapMax = -2.5;
            meanVec.push_back(0.089);
            meanErrorVec.push_back(0.006);
            sigmaVec.push_back(0.036);
            sigmaErrorVec.push_back(0.004);
            corrVec.push_back(0.64);
        }
        else if (rapMode == 2) {
            rapMin = -4.;
            rapMax = -3.25;
            meanVec.push_back(0.093);
            meanErrorVec.push_back(0.008);
            sigmaVec.push_back(0.032);
            sigmaErrorVec.push_back(0.004);
            corrVec.push_back(0.72);
        }
        else if (rapMode == 3) {
            rapMin = -3.25;
            rapMax = -2.5;
            meanVec.push_back(0.083);
            meanErrorVec.push_back(0.009);
            sigmaVec.push_back(0.038);
            sigmaErrorVec.push_back(0.006);
            corrVec.push_back(0.54);
        }
        else {
            cout << "Unknown rapidity mode, terminating" << endl;
            return;
        }
    }
    else {
        cout << "Unknown mass mode, terminating" << endl;
        return;
    }
    
    /*
     vector<double> ggMeanVec = {4.16, 5.93};
     vector<double> ggMeanErrorVec = {0.17, 0.36};
     
     vector<double> ggMeanVec = {4.16, 5.52};
     vector<double> ggMeanErrorVec = {0.26, 0.52};
     
     vector<double> ggMeanVec = {3.86, 6.31};
     vector<double> ggMeanErrorVec = {0.19, 0.40};
     */
    
    const int nPeriods = periods.size();
    for (int k = 0; k<nPeriods; k++) {
        string period = periods[k];
        double r = corrVec[k];
        
        int nMu = 0, nSigma = 0;
        bool init = MakeTree(period, rapMin, rapMax, massMin, massMax, meanVec[k], meanErrorVec[k], sigmaVec[k], sigmaErrorVec[k], nMu, nSigma);
        if (!init) return;
        //PlotResults(period, rapMin, rapMax, massMin, massMax);
        SystematicsComputation(period, rapMin, rapMax, massMin, massMax, meanVec[k], meanErrorVec[k], sigmaVec[k], sigmaErrorVec[k], nMu, nSigma, r);
    }
    
}
