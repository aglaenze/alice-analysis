#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>

#include "TString.h"
#include "TDatime.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"
#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>
#include <TFitResultPtr.h>

using namespace std;

template <typename T>
bool contains(std::list<T> & listOfElements, const T & element)
{
    // Find the iterator if element in list
    auto it = std::find(listOfElements.begin(), listOfElements.end(), element);
    //return if iterator points to end or not. It points to end then it means element
    // does not exists in list
    return it != listOfElements.end();
}

Double_t FitLin(Double_t* x, Double_t* par) {
    return  par[0]*x[0];
}

struct RunInfo {
    Int_t runNumber=0;
    Double_t mu=0;
    Double_t mu2=0;
    Int_t ctrueCounter=0;
    Int_t cutsCounter[7]={0,0,0,0,0,0,0};
    Bool_t trdActive = kFALSE;
};

void WriteTree(TString path, string period, vector<TString> cutNames) {
    
    // import data
    //TFile* fAna = TFile::Open("AnalysisResults.root", "READ");
    
    TString rootFileName = path + Form("AnalysisResults_%s_mu.root", period.c_str());
    
    TFile* fAna = TFile::Open(rootFileName, "READ");
    TTree* fAnaTree = (TTree*) fAna->Get("MyTask/AnaTree");
    
    // the two following lists come from Guillermo's presentation
    std::list<Int_t> selectedRuns;
    std::list<Int_t> trdRuns;
    if (period == "LHC16r") {
        //selectedRuns = {265589, 265594, 265596, 265607, 265669, 265691, 265694, 265742, 265744, 265746, 265754, 265756, 265785, 265787, 266074, 266076, 266081, 266084, 266085, 266086, 266117, 266299, 266300, 266304, 266305, 266312, 266316, 266318, 265696, 265697, 265698, 265700, 265701, 265709, 265788, 265789, 265792, 265795, 265797, 265840, 266022, 266023, 266025, 266034, 266208, 266234, 266235, 266296, 266190, 266193, 266196, 266197};
        selectedRuns = {265589, 265594, 265596, 265607, 265669, 265691, 265694, 265697, 265698, 265700, 265701, 265709, 265713, 265714, 265740, 265741, 265742, 265744, 265746, 265754, 265756, 265785, 265787, 265788, 265789, 265792, 265795, 265797, 265840, 266022, 266023, 266025, 266034, 266074, 266076, 266081, 266084, 266085, 266086, 266117, 266187, 266189, 266190, 266193, 266196, 266197, 266208, 266234, 266235, 266296, 266299, 266300, 266304, 266305, 266312, 266316, 266318};
        trdRuns = {265589, 265594, 265596, 265607, 265669, 265696, 265697, 265698, 265700, 265701, 265709, 265713, 265714, 265741, 265742, 265744, 265746, 265754, 265756, 265788, 265789, 265792, 265795, 265797, 266034, 266074, 266076, 266081, 266084, 266085, 266086, 266117, 266187, 266189, 266190, 266193, 266196, 266197, 266208, 266296, 266299, 266300, 266304, 266305, 266316, 266318 };
    }
    else if (period == "LHC16s") {
        //selectedRuns = {266405, 266437, 266438, 266439, 266441, 266470, 266472, 266533, 266534, 266539, 266543, 266549, 266584, 266587, 266657, 266658, 266659, 266665, 266668, 266669, 266674, 266857, 266878, 266880, 266882, 266883, 266885, 266886, 266998, 267020, 267022, 267062, 267063, 267067, 267070, 266479, 266480, 266487, 266514, 266516, 266518, 266588, 266591, 266593, 266595, 266613, 266614, 266676, 266702, 266703, 266706, 266708, 266775, 266912, 266915, 266940, 266942, 266943, 266944, 267072, 267077, 267109, 267110, 267130, 267131, 266520, 266522, 266523, 266525, 266615, 266618, 266621, 266630, 266776, 266800, 266805, 266807, 266988, 266993, 266994, 266997};
        selectedRuns = {266439, 266441, 266472, 266479, 266480, 266487, 266514, 266516, 266518, 266520, 266522, 266523, 266525, 266533, 266534, 266539, 266543, 266549, 266584, 266587, 266588, 266591, 266593, 266595, 266613, 266614, 266618, 266621, 266630, 266657, 266658, 266659, 266665, 266668, 266669, 266674, 266676, 266702, 266703, 266706, 266708, 266775, 266776, 266800, 266805, 266807, 266857, 266878, 266880, 266882, 266883, 266885, 266886, 266912, 266915, 266940, 266942, 266943, 266944, 266988, 266993, 266994, 266997, 266998, 267020, 267022, 267062, 267063, 267067, 267070, 267072, 267077, 267109, 267110, 267130, 267131};
        trdRuns = {266439, 266441, 266487, 266543, 266549, 266630, 266882, 266883, 266885, 266886, 266944, 266993, 266994, 266997, 266998, 267070, 267072, 267077, 267110};
    }
    
    //return;
    // Data from the anaTree
    Int_t fRunNum;
    Double_t fMu;
    Bool_t fCtrue;
    Int_t nCuts = cutNames.size();
    Bool_t cutList[nCuts];
    fAnaTree->SetBranchAddress("fRunNum",&fRunNum);
    fAnaTree->SetBranchAddress("fCtrue",&fCtrue);
    fAnaTree->SetBranchAddress("fMu",&fMu);
    for (int i = 0; i< nCuts; i++) fAnaTree->SetBranchAddress(cutNames[i], &cutList[i]);
    
    std::map<string, RunInfo > allRunData;
    
    // Browse over all data to collect run numbers
    Int_t nEvents = fAnaTree->GetEntries();
    //nEvents = 100000;
    for (int k = 0; k < nEvents; k++) {
        //for (int k = 0; k < 20; k++) {
        fAnaTree->GetEntry(k);
        if (!contains(selectedRuns, fRunNum)) continue;
        std::string numString = std::to_string(fRunNum);
        if (allRunData.find(numString) == allRunData.end()) { // run was not found
            RunInfo newRun;
            newRun.runNumber = fRunNum;
            allRunData.insert(std::make_pair(numString, newRun));
        }
    }
    
    // open the tree that contains the real mu
    TFile* fTrending = TFile::Open(path + "trending_merged.root", "READ");
    TTree* fTrendingTree = (TTree*) fTrending->Get("trending");
    
    // data from trending tree
    Int_t nEntries2 = fTrendingTree->GetEntries();
    Int_t fRunNum2;
    Double_t mu2;
    fTrendingTree->SetBranchAddress("mu",&mu2);
    fTrendingTree->SetBranchAddress("run",&fRunNum2);
    for (int k = 0; k < nEntries2; k++) {
        fTrendingTree->GetEntry(k);
        std::string numString = std::to_string(fRunNum2);
        // only select runs among the ones compiled in the analysis
        if(allRunData.find(numString) == allRunData.end()) continue;
        allRunData[numString].mu2=mu2;
    }
    
    
    std::vector <TH1*> cutsHistTrd;
    std::vector <TH1*> cutsHistNoTrd;
    for (int i = 0; i< nCuts; i++) {
        //TH1D* h = new TH1D(cutNames[i], cutNames[i], 100, 0, 0.08);
        cutsHistTrd.push_back(new TH1D(cutNames[i]+"Trd", cutNames[i]+"Trd", 2000, 0, 0.08));
        cutsHistNoTrd.push_back(new TH1D(cutNames[i]+"NoTrd", cutNames[i]+"NoTrd", 2000, 0, 0.08));
    }
    TH1* ctrueHist = new TH1D("CTrue", "CTrue", 2000, 0, 0.08);
    // collect selection info
    for (int k = 0; k < nEvents; k++) {
        //for (int k = 0; k < 20; k++) {
        fAnaTree->GetEntry(k);
        std::string numString = std::to_string(fRunNum);
        if (allRunData.find(numString) == allRunData.end()) continue;// run was not found
        allRunData[numString].mu=fMu;
        allRunData[numString].trdActive=contains(trdRuns, fRunNum);
        if (fCtrue) {
            allRunData[numString].ctrueCounter++;
            for (int i = 0; i< nCuts; i++) {
                if (cutList[i]) {
                    allRunData[numString].cutsCounter[i]++;
                    if (allRunData[numString].trdActive) cutsHistTrd[i]->Fill(allRunData[numString].mu2);
                    else cutsHistNoTrd[i]->Fill(allRunData[numString].mu2);
                }
            }
            ctrueHist->Fill(allRunData[numString].mu2);
        }
        //else if (!fCtrue && cutList[i]) std::cout << "\n\n\nproblem!\n\n" << std::endl; // should never happen because the cuts are obtained only in ctrue events
        
    }
    
    //now write data per run in a new tree
    
    TString muAnalysisFile = path + Form("MuAnalysis_%s.root", period.c_str());
    TFile* f = new TFile(muAnalysisFile,"RECREATE");
    TTree* newTree = new TTree("muAnalysisTree","Tree with cuts and mu parameter");
    
    RunInfo runA;
    Double_t cutsCounterNorm[nCuts];
    newTree->Branch("runNum",&runA.runNumber,"runNum/I");
    newTree->Branch("mu",&runA.mu,"mu/D");
    newTree->Branch("mu2",&runA.mu2,"mu2/D");
    newTree->Branch("fCtrue",&runA.ctrueCounter,"fCtrue/I");
    newTree->Branch("fTrdActive",&runA.trdActive,"fTrdActive/O");
    for (int i = 0; i< nCuts; i++){
        newTree->Branch(cutNames[i]+"andfCtrue", &runA.cutsCounter[i], cutNames[i]+"andfCtrue/I");
        newTree->Branch(cutNames[i]+"andfCtrueNorm", &cutsCounterNorm[i], cutNames[i]+"andfCtrueNorm/D");
    }
    
    // normalize the selection
    for(std::map<string,RunInfo>::iterator iter = allRunData.begin(); iter != allRunData.end(); ++iter)
    {
        string key =  iter->first;
        runA = iter->second;
        //runA = allRunData[key];
        for (int i = 0; i< nCuts; i++){
            if (runA.ctrueCounter!=0) cutsCounterNorm[i]= (double)runA.cutsCounter[i]/runA.ctrueCounter;
            else cutsCounterNorm[i]= 0;
        }
        newTree->Fill();
    }
    newTree->Write();
    ctrueHist->Write();
    for (int i = 0; i< nCuts; i++) {cutsHistTrd[i]->Write(); cutsHistNoTrd[i]->Write();}
    f->Close();
    return;
}




//void PlotMu(Int_t period = 0) {
void PlotMu() {
    
    gStyle->SetOptStat("ne");
    gStyle->SetOptFit(1);
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetTitleSize(.05);
    gStyle->SetLabelSize(.05, "XY");
    gStyle->SetMarkerSize(0.5);
    gStyle->SetMarkerStyle(20);
    
    TString path = "/Users/aglaenzer/alice/analysis-alice/CTRUE/rootFiles/";
    
    vector<string> periods = {"LHC16r", "LHC16s"};
    for (int l = 0; l < (int)periods.size(); l++) {
        string period = periods[l];
        //WriteTree(path, period);
        //return;
        
        TString muAnalysisFile =  Form("MuAnalysis_%s", period.c_str());
        TFile* fNewTree = TFile::Open(path + muAnalysisFile+".root", "READ");
        //TTree* newTree = (TTree*) fNewTree->Get("muAnalysisTree");
        
        TString tag = "fCut14";
        vector<TString> cutNames = {tag+"on", tag+"offDAA", tag+"offDVA", tag+"offDVC", tag+"offZN", tag+"off", tag};
        const int nCuts = cutNames.size();
        
        TCanvas* cv = new TCanvas("cv", "cv", 600, 1200);
        cv->Divide(2, 4);
        
        TH1D* hCtrue = (TH1D*)fNewTree->Get("CTrue");
        for (int i = 0; i< nCuts; i++) {
            TH1D* hCutTrd = (TH1D*)fNewTree->Get(cutNames[i]+"Trd");
            TH1D* hCutNoTrd = (TH1D*)fNewTree->Get(cutNames[i]+"NoTrd");
            // Define the ratio plot
            TH1D *hCutTrdNorm = (TH1D*)hCutTrd->Clone(cutNames[i]+"TrdNorm");
            TH1D *hCutNoTrdNorm = (TH1D*)hCutNoTrd->Clone(cutNames[i]+"NoTrdNorm");
            hCutTrdNorm->Sumw2();
            hCutNoTrdNorm->Sumw2();
            hCutTrdNorm->Divide(hCtrue);
            hCutNoTrdNorm->Divide(hCtrue);
            hCutNoTrdNorm->SetMinimum(0.);  // Define Y ..
            if (period == "LHC16r") hCutNoTrdNorm->SetMaximum(0.4); // .. range
            else hCutNoTrdNorm->SetMaximum(1); // .. range
            hCutTrdNorm->SetMarkerColor(kBlue);
            hCutTrdNorm->SetLineColor(kBlue);
            hCutNoTrdNorm->SetMarkerColor(kRed);
            hCutNoTrdNorm->SetLineColor(kRed);
            cv->cd(i+1);
            //hCutTrdNorm->Draw("ep");
            //hCutTrdNorm->Draw("hist");
            //TF1* hfit = new TF1("pol1","pol1", 0.04, 0.06);
            //hCutNoTrdNorm->Fit(hfit,"","", 0.04, 0.078);
            //hCutNoTrdNorm->Fit("pol1","","", 0.04, 0.078);
            
            if (period == "LHC16s") {
                TString tag2 = "fCut23";
                vector<TString> cutNames2 = {tag2+"on", tag2+"offDAA", tag2+"offDVA", tag2+"offDVC", tag2+"offZN", tag2+"off", tag2};
                hCutNoTrdNorm->SetTitle(cutNames2[i]);
            }
            else {
                hCutNoTrdNorm->SetTitle(cutNames[i]);
            }
            
            Double_t fitRangeMin = 0.04, fitRangeMax = 0.078;
            TF1 *f1 = new TF1("f1", FitLin, fitRangeMin, fitRangeMax, 1);
            hCutNoTrdNorm->Fit(f1, "0", "0", fitRangeMin, fitRangeMax);
            hCutNoTrdNorm->Draw("ep");
            f1->Draw("same");
            
            gPad->Update();
            TPaveStats *st = (TPaveStats*)hCutNoTrdNorm->FindObject("stats");
            st->SetOptStat(0);
            st->SetOptFit(1);
            TList *listOfLines = st->GetListOfLines();
            // Remove lines // ça marche pas...
            TText *p0 = st->GetLineWith("p0");
            //listOfLines->Remove(p0);
            
            st->Draw();
            hCutTrdNorm->Draw("same");
            
        }
        
        
        cv->SaveAs("Plots/"+muAnalysisFile+".pdf");
    }
    
    return;
}

