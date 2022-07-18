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
#include <TCut.h>
#include <TTree.h>
#include <TLegend.h>

using namespace std;


bool contains(std::map<Int_t, Int_t> &myMap, const Int_t &element)
{
    for(std::map<Int_t, Int_t>::iterator iter = myMap.begin(); iter != myMap.end(); ++iter)
    {
        Int_t k =  iter->first;
        if (element == iter->first) return true;
    }
    return false;
}

int GetTotalEff(string process, string period, TCut cut, TCut cutMc, Int_t& nGen, Int_t& nRec) {
    
    TString path = "/Volumes/Transcend/rootFiles-pPb/MC-std/";
    TString fileName = Form("%s_MC_%s.root", period.c_str(), process.c_str());
    TFile *f = new TFile(path+ "AnalysisResults_"+fileName,"read");
    
    bool draw = true;
    
    TTree* fAnaTree = (TTree*)f->Get("MyTask/fAnaTree");
    TTree* fGenTree = (TTree*)f->Get("MyTask/fGenTree");
    
    
    TCanvas* c0 = new TCanvas();    // To avoid the warning
    
    TH1F* h1 = new TH1F("h1","", 100, -10, 1000);
    TH1F* h2 = new TH1F("h2","", 100, -10, 1000);
    fAnaTree->Draw("fTrkTrkM>>h1", cut);
    fGenTree->Draw("fMCTrkTrkM>>h2", cutMc);
    
    nRec = h1->GetEntries();
    nGen = h2->GetEntries();
    
    delete h1;
    delete h2;
    
    return 0;
}



int DrawEffPt(string process = "kIncohJpsiToMu", vector<string> periods = {"LHC16r"}, Double_t mMin = 2.5, Double_t mMax = 3.5, Double_t rapMin = 2.5, Double_t rapMax = 4.0) {
    
    // Cuts for generated MC
    TCut rapCutMc = Form("fMCTrkTrkY < -%f && fMCTrkTrkY > -%f", rapMin, rapMax);
    //TCut rapCutMc = Form("fMCTrkTrkY < -%f && fMCTrkTrkY > -%f", -5.0, 5.0);
    TCut unlikeSignCutMc = "(fMCTrkQ1 < 0 && fMCTrkQ2 > 0) || (fMCTrkQ1 > 0 && fMCTrkQ2 < 0)";
    // Don't use the mass cut in generated events
    TCut massCutMc = Form("fMCTrkTrkM < %f && fMCTrkTrkM > %f", mMax, mMin);
    
    // Cuts for reconstructed
    TCut rapCut = Form("fTrkTrkY < -%f && fTrkTrkY > -%f", rapMin, rapMax);
    TCut triggerCut;    // defined in the for loop
    TCut unlikeSignCut = "(fTrkQ1 < 0 && fTrkQ2 > 0) || (fTrkQ1 > 0 && fTrkQ2 < 0)";
    TCut massCut = Form("fTrkTrkM < %f && fTrkTrkM > %f", mMax, mMin);
    
    for (int l = 0; l < (int)periods.size(); l++) {
        string period = periods[l];
        if (period == "LHC16r") triggerCut = "0<(fL0inputs&(1<<5))";
        else if (period == "LHC16s") triggerCut = "0<(fL0inputs&(1<<13))";
        else {cout << "problem" << endl; return -1;}
        cout << endl;
        
        cout << "Going through period " << period << " and process " << process << endl;
        
        Int_t nGen;
        Int_t nRec;
        
        // First: events with no cut
        TCut cutMc = "";
        TCut cut = triggerCut;
        
        // Second: events within the rapidity cut
        cut += rapCut;
        cutMc += rapCutMc;
        
        // Third: unlike sign events
        cut += unlikeSignCut;
        cutMc += unlikeSignCutMc;
        
        // Fourth: events within the given mass range
        cut += massCut;
        cutMc += massCutMc;
        GetTotalEff(process, periods[l], cut, cutMc, nGen, nRec);

        
        float eff = (float)nRec/nGen;
        
        float effExc = 0, effDiss = 0;

        // Fifth: histogram as a function of pt
        Int_t nBins = 10;
        double weightsExc[nBins];
        double weightsDiss[nBins];
        double weightsExcSum = 0, weightsDissSum = 0;
        Double_t ptMin = 0, ptMax = 3;
        Double_t step = (ptMax - ptMin) / nBins;
        TH1F* histPt = new TH1F("histPt", "AxE = f(p_{T})", nBins, ptMin, ptMax);
        for (int i = 0; i < nBins; i++) {
            //TCut newCut = cut;
            TCut newCut = cut + Form("fTrkTrkPt>%.1f && fTrkTrkPt<%.1f", ptMin+i*step, ptMin+(i+1)*step);
            TCut newCutMc = cutMc + Form("fMCTrkTrkPt>%.1f && fMCTrkTrkPt<%.1f", ptMin+i*step, ptMin+(i+1)*step);
            //TCut newCutMc = cutMc;
            GetTotalEff(process, periods[l], newCut, newCutMc, nGen, nRec);
            if (nGen == 0) continue;

     
            cout << "Efficiency = " << (float)nRec/nGen << " between " << ptMin+i*step << " and " << ptMin+(i+1)*step << endl;
            //histPt->Fill(ptMin+(i+0.5)*step, (float)nRec/nGen);
            histPt->SetBinContent(i+1, (float)nRec/nGen);
            histPt->SetBinError(i+1, TMath::Sqrt((float)1/nGen));
            
            // compute weights
            double pT = ptMin+(i+0.5)*step;
            double bExc = 3.62;
            weightsExc[i] = pT*TMath::Exp(-bExc*pT*pT);
            weightsExcSum += weightsExc[i];
            double bDiss = 0.77;
            double nDiss = 6.19;
            weightsDiss[i] = pT/pow(1+pT*pT+(bDiss/nDiss), nDiss);
            weightsDissSum += weightsDiss[i];
            
            effExc += weightsExc[i]*(float)nRec/nGen;
            effDiss += weightsDiss[i]*(float)nRec/nGen;
        }
        effExc /= weightsExcSum;
        effDiss /= weightsDissSum;
        std::cout << "effExc = " << effExc << std::endl;
        std::cout << "effDiss = " << effDiss << std::endl;
        
        histPt->SaveAs(Form("test-%s.root", period.c_str()));
        TCanvas* cv = new TCanvas();
        histPt->Draw();
        cv->SaveAs(Form("Plots/Efficiency-vs-pt-%s.pdf", period.c_str()));
        
    }
    return 0;
}


