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
#include <TTree.h>
#include <TGraph.h>
#include <TF1.h>
#include <TPaveStats.h>
#include <TText.h>

/*
 Dans cette macro je veux obtenir la fraction des cuts en fonction de mu, sachant que mu est défini par run
 Le plus simple serait de faire un histogramme par cut, rempli à chaque fois, puis le diviser par l'autre
 USE TIMING INFO INSTEAD OF ZN ENERGY
 */


using namespace std;

vector<Int_t> GetTrdList(string period) {
    vector <Int_t> trdRuns16r = {265589, 265594, 265596, 265607, 265669, 265696, 265697, 265698, 265700, 265701, 265709, 265713, 265714, 265741, 265742, 265744, 265746, 265754, 265756, 265788, 265789, 265792, 265795, 265797, 266034, 266074, 266076, 266081, 266084, 266085, 266086, 266117, 266187, 266189, 266190, 266193, 266196, 266197, 266208, 266296, 266299, 266300, 266304, 266305, 266316, 266318};
    vector <Int_t> trdRuns16s = {266439, 266441, 266487, 266543, 266549, 266630, 266882, 266883, 266885, 266886, 266944, 266993, 266994, 266997, 266998, 267070, 267072, 267077, 267110};
    if (period == "LHC16r")
        return trdRuns16r;
    else if (period == "LHC16s")
        return trdRuns16s;
    else {
        cout << "Wrong period" << endl;
        return {};
    }
}

Double_t FitLin(Double_t* x, Double_t* par) {
    return par[0]*x[0];
}

Double_t GetMu(map<Int_t, Double_t> mymap, Int_t runNum) {
    
    map<Int_t, Double_t>::iterator it = mymap.begin();
    while (it != mymap.end()) {
        if (it->first == runNum) {
            return it->second;
        }
        it++;
    }
    cout << "run not found" << endl;
    return -1;
}

Bool_t Contains(vector <Int_t> runList, Int_t run) {
    for (int l = 0; l < (int)runList.size(); l++) {
        if (runList[l] == run) return true;
    }
    return false;
}

map<Int_t, Double_t> GetMuMap(TString path, vector <Int_t> runList) {
    /* tree to get mu */
    map <Int_t, Double_t> mymap = {};
    TFile* fTrending = TFile::Open(path + "trending_merged.root", "READ");
    TTree* fTrendingTree = (TTree*) fTrending->Get("trending");
    Double_t mu;
    Int_t run = -1;
    fTrendingTree->SetBranchAddress("mu",&mu);
    fTrendingTree->SetBranchAddress("run",&run);
    for (int i = 0; i < fTrendingTree->GetEntries(); i++) {
        fTrendingTree->GetEntry(i);
        if (Contains(runList, run)) mymap.insert ( pair<Int_t, Double_t>(run,mu) );
    }
    return mymap;
}

/*
 1e étage : récupérer tous les histos
 */
void GetCutsHists(TString path, string period, bool exc = true) {
    
    TCut emptyEvents = "", fCuton = "", fCutoffAD = "", fCutoffV0 = "", fCutoffADp = "", fCutoffZN = "", fCutoff = "", fCut = "";
    TCut fZN = "";
    int nBins = 1000;
    
    TFile* f = new TFile(path+Form("AnalysisResults_%s_mu.root", period.c_str()));
    cout << "Reading " + path + Form("AnalysisResults_%s_mu.root", period.c_str()) << endl;
    
    TTree* fAnaTree = (TTree*)f->Get("MyTask/fAnaTree");
    /* Get info from TTree */
    Int_t fRunNum;
    fAnaTree->SetBranchAddress("fRunNum", &fRunNum);
    
    Bool_t f0VBA(0), f0VGA(0), f0UBA(0), f0VBC(0), f0UGC(0), f0UBC(0), f0VC5(0), f0SH2(0);
    fAnaTree->SetBranchAddress("f0VBA", &f0VBA);
    fAnaTree->SetBranchAddress("f0VGA", &f0VGA);
    fAnaTree->SetBranchAddress("f0UBA", &f0UBA);
    fAnaTree->SetBranchAddress("f0VBC", &f0VBC);
    fAnaTree->SetBranchAddress("f0UGC", &f0UGC);
    fAnaTree->SetBranchAddress("f0UBC", &f0UBC);
    fAnaTree->SetBranchAddress("f0VC5", &f0VC5);
    fAnaTree->SetBranchAddress("f0SH2", &f0SH2);
    
    Double_t fZNCEnergy, fZNAEnergy;
    Bool_t fZNAgoodTiming, fZNCgoodTiming, fZNAgoodTiming2, fZNCgoodTiming2;
    Int_t fV0ADecision, fV0CDecision, fADADecision, fADCDecision;
    fAnaTree->SetBranchAddress("fZNCEnergy", &fZNCEnergy);
    fAnaTree->SetBranchAddress("fZNAEnergy", &fZNAEnergy);
    fAnaTree->SetBranchAddress("fZNAgoodTiming", &fZNAgoodTiming);
    fAnaTree->SetBranchAddress("fZNCgoodTiming", &fZNCgoodTiming);
    fAnaTree->SetBranchAddress("fZNAgoodTiming2", &fZNAgoodTiming2);
    fAnaTree->SetBranchAddress("fZNCgoodTiming2", &fZNCgoodTiming2);
    fAnaTree->SetBranchAddress("fV0ADecision", &fV0ADecision);
    fAnaTree->SetBranchAddress("fV0CDecision", &fV0CDecision);
    fAnaTree->SetBranchAddress("fADADecision", &fADADecision);
    fAnaTree->SetBranchAddress("fADCDecision", &fADCDecision);
    
    vector<Int_t> runList = {};
    for (int i = 0; i < fAnaTree->GetEntries(); i++) {
        fAnaTree->GetEntry(i);
        if (!Contains(runList, fRunNum)) {
            runList.push_back(fRunNum);
        }
    }
    
    // Get list of mu
    map<Int_t, Double_t> mumap = GetMuMap(path, runList);
    cout << "Size of the map = " << (int)mumap.size() << endl;
    
    // Definition of cuts
    if (period == "LHC16r") {
        // study the behaviour of the online tigger
        emptyEvents = "fADCDecision == 0 && fV0CDecision == 0 && !fZNCgoodTiming";
        fZN = "fZNAgoodTiming" && emptyEvents;
        /*
         fCuton = emptyEvents && "f0UBA || f0VBA";
         
         // study the ADA decision
         fCutoffAD = fCuton && "(fADADecision == 1)";
         
         // study the V0 decision
         fCutoffV0 = !(fCuton || fCutoffAD ) && "(fV0ADecision == 1)";
         
         if (exc) {
         // study the AD decision in the proton side
         fCutoffADp = !(fCuton || fCutoffAD || fCutoffV0) && "(fADCDecision == 1)";
         // study the ZN decision on both sides
         fCutoffZN = !(fCuton || fCutoffAD || fCutoffV0 || fCutoffADp) && "(fZNAgoodTiming || fZNCgoodTiming)"; //(fZNAEnergy>fZNAEnergyLimit || fZNCEnergy>fZNCEnergyLimit);
         }
         
         else {
         fCutoffADp = ""; // no cut in the proton side for the loose selection
         // study the ZN decision on Pb side only
         fCutoffZN = !(fCuton || fCutoffAD || fCutoffV0 || fCutoffADp) && "(fZNAgoodTiming)";
         }
         */
    }
    
    /*
     else if (period == "LHC16s") {
     fCuton = "f0UBC || f0UGC || f0VBA || f0VGA || f0SH2 || f0VC5";
     
     fCutoffAD = !fCuton && "(fADCDecision == 1)";
     fCutoffV0 = !(fCuton || fCutoffAD ) && "(fV0ADecision == 1)";
     if (exc ) {
     fCutoffADp = !(fCuton || fCutoffAD || fCutoffV0) && "(fADADecision == 1)";
     fCutoffZN = !(fCuton || fCutoffAD || fCutoffV0 || fCutoffADp) && "(fZNAgoodTiming || fZNCgoodTiming)";
     }
     else {
     fCutoffADp = ""; // no cut in the proton side for the loose selection
     fCutoffZN = !(fCuton || fCutoffAD || fCutoffV0 || fCutoffADp) && "(fZNCgoodTiming)";
     }
     }
     */
    else return;
    /*
     // study all offline vetoes (exclusive selection)
     fCutoff = !(fCuton) && (fCutoffAD || fCutoffV0 || fCutoffADp || fCutoffZN);
     // all vetoes at the same time
     fCut = fCuton || fCutoffAD || fCutoffV0 || fCutoffADp || fCutoffZN;
     */
    
    //vector<TCut> cutList = {emptyEvents, fCuton, fCutoffAD, fCutoffV0, fCutoffADp, fCutoffZN, fCutoff, fCut};
    vector<TCut> cutList = {emptyEvents, fZN};
    const int nCuts = (int)cutList.size();
    const int nn = (int)runList.size();
    
    vector<Int_t> runTrdList = GetTrdList(period);
    vector<Double_t> muVec = {}, muVecTrd = {};
    vector<Double_t> proportions[nCuts], proportionsTrd[nCuts];
    vector<Double_t> propErrors[nCuts], propErrorsTrd[nCuts];
    Int_t nEntries = 0, nTotal = 0;
    for (int i = 0; i < nn; i++) {
        TCut runCut = Form("fRunNum == %d", runList[i]);
        if (Contains(runTrdList, runList[i])) {muVecTrd.push_back(GetMu(mumap, runList[i]));}
        else {muVec.push_back(GetMu(mumap, runList[i]));}
        for (int j = 0; j<nCuts; j++) {
            TH1F* hist = new TH1F("hist", "hist", nBins, 0, 0.08);
            fAnaTree->Draw("fRunNum>>hist", cutList[j] && runCut);
            nEntries = hist->GetEntries();
            if (j==0) nTotal = nEntries;
            // Compute errors using binomial formula
            // variance = npq
            double p = (double)nEntries / nTotal;
            double q = (1-p);
            if (Contains(runTrdList, runList[i])) {
                proportionsTrd[j].push_back(p);
                propErrorsTrd[j].push_back(p*p*q);
            }
            else {
                proportions[j].push_back(p);
                propErrors[j].push_back(p*p*q);
            }
            
            delete hist;
        }
    }
    
    vector <TGraphErrors*> grList = {};
    vector <TGraphErrors*> grListTrd = {};
    const int nProp = (int)proportions[0].size();
    const int nPropTrd = (int)proportionsTrd[0].size();
    Double_t li1[nProp], li2[nProp], li2Err[nProp];
    Double_t li1Trd[nPropTrd], li2Trd[nPropTrd], li2ErrTrd[nPropTrd];
    // Create TGaphErrors
    for (int j = 0; j<nCuts; j++) {
        for (int i = 0; i < nProp; i++) {
            li1[i] = muVec[i];
            li2[i] = proportions[j][i];
            li2Err[i] = propErrors[j][i];
        }
        TGraphErrors* gr = new TGraphErrors(nProp, li1, li2, 0, li2Err);
        grList.push_back(gr);
        for (int i = 0; i < nPropTrd; i++) {
            li1Trd[i] = muVecTrd[i];
            li2Trd[i] = proportionsTrd[j][i];
            li2ErrTrd[i] = propErrorsTrd[j][i];
        }
        TGraphErrors* grTrd = new TGraphErrors(nPropTrd, li1Trd, li2Trd, 0, li2ErrTrd);
        grListTrd.push_back(grTrd);
    }
    
    // Draw them
    int n2 = (int)nCuts/2;
    if ((double)nCuts/2 > n2) n2++;
    TCanvas* cv = new TCanvas("cv", "cv", 600, 300*n2);
    cv->Divide(2, n2);
    
    for (int j = 0; j<nCuts; j++) {
        cv->cd(j+1);
        grList[j]->SetMinimum(0);
        grList[j]->SetMaximum(1.2);
        /*
         grList[j]->SetMaximum(0.4);
         if (j==0) grList[j]->SetMaximum(1.2);
         */
        grList[j]->GetXaxis()->SetLimits(0.0,0.08);
        grList[j]->SetTitle(cutList[j]);
        grList[j]->SetMarkerColor(kRed);
        grList[j]->Draw("AP");
        grListTrd[j]->SetMarkerColor(kBlue);
        grListTrd[j]->Draw("P same");
        
        // Fit
        if (j>0) {
            Double_t fitRangeMin = 0.04, fitRangeMax = 0.078;
            TF1 *f1 = new TF1("f1", FitLin, fitRangeMin, fitRangeMax, 1);
            grList[j]->Fit(f1, "0", "0", fitRangeMin, fitRangeMax);
            f1->Draw("same");
            
            gPad->Update();
            TPaveStats *st = (TPaveStats*)grList[j]->FindObject("stats");
            st->SetOptStat(0);
            st->SetOptFit(1);
            TList *listOfLines = st->GetListOfLines();
            // Remove lines
            TText *p0 = st->GetLineWith("p0");
            //listOfLines->Remove(p0);
            st->Draw();
        }
    }
    
    if (exc) cv->SaveAs(Form("Plots/cuts-hists-%s-exc.pdf", period.c_str()));
    else cv->SaveAs(Form("Plots/cuts-hists-%s-loose.pdf", period.c_str()));
}


void Ana() {
    
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
    
    //vector<string> periods = {"LHC16r", "LHC16s"};
    vector<string> periods = {"LHC16r"};
    for (int l = 0; l < (int)periods.size(); l++) {
        string period = periods[l];
        GetCutsHists(path, period);
    }
    return;
}

