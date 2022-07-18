#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>

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
#include "RooPlot.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"

#include "../Include/Initiate.C"

using namespace std;


void SaveMassAndPt(int testCase = 3) {  // 0 no selection, 1 standard selection, 2 exclusive selection, 3 reconstructed MC
    
    bool exclusiveOnly = false;
    bool noCut = true;
    
    Double_t mMin = 1., mMax = 2.5, ptMin = 0., ptMax = 3;
    Double_t rapMin = -4, rapMax = -2.5;
    
    string rootfilePathData = "/Volumes/Transcend/rootFiles-pPb/std";
    string rootfilePathMc = "/Volumes/Transcend/rootFiles-pPb/MC-std";
    string rootfilePath = rootfilePathData;
    
    if (testCase == 1) noCut = false;
    else if (testCase == 2) {
        noCut = false;
        exclusiveOnly = true;
        ptMax = 1.2;
    }
    else if (testCase == 3) {   // MC
        noCut = true;
        rootfilePath = rootfilePathMc;
    }
    
    //vector<string> periods = {"LHC16r", "LHC16s"};
    vector<string> periods = {"LHC16r"};
    

    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.07);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetTitleSize(.07);
    gStyle->SetTextSize(.07);
    gStyle->SetLabelSize(.05, "XY");
    //gStyle->SetMarkerSize(0.5);
    //gStyle->SetMarkerStyle(20);

    const int nPeriod = periods.size();
    
    for (int k = 0; k<nPeriod; k++) {
        string period = periods[k];
        if (exclusiveOnly) {ptMax = 1.2;}
        string fileName = Form("AnalysisResults_%s.root", period.c_str());
        if (testCase == 3) {
            //fileName = Form("AnalysisResults_%s_MC_kIncohJpsiToMu.root", period.c_str());
            fileName = Form("AnalysisResults_%s_MC_kTwoGammaToMuLow.root", period.c_str());
        }
        
        // Define cuts
        std::list<TCut> mCutList = DefineCuts(period, exclusiveOnly);
        TCut mCut = "";
        
        for (std::list <TCut>::iterator iter = mCutList.begin(); iter != mCutList.end(); ++iter) {mCut += *iter;}
        
        // Open the file
        TFile *fAna = new TFile(Form("%s/%s", rootfilePath.c_str(), fileName.c_str()),"READ");
        
        TString histnameM = "mhist-sel";
        TString histnamePt = "pthist-sel";
        if (exclusiveOnly) {
            histnameM += "-exc";
            histnamePt += "-exc";
        }
        
        if (noCut) {
            mCut = "";
            histnameM = "mhist";
            histnamePt = "pthist";
        }
        if (testCase == 3) {
            histnameM += "-mc";
            histnamePt += "-mc";
        }
        mCut += Form("fTrkTrkM > %.1f && fTrkTrkM < %1f", mMin, mMax);
        mCut += Form("fTrkTrkPt > %.1f && fTrkTrkPt < %1f", ptMin, ptMax);
        
        int nBinsMass = 250;
        int nBinsPt = int((ptMax-ptMin)*100);
        if (testCase == 3) {
            nBinsMass = 2500;
            nBinsPt = int((ptMax-ptMin)*1000);
        }
        
        //int nBins = 125;
        TH1* hM = new TH1F(histnameM, histnameM, nBinsMass, mMin, mMax);
        TH1* hPt = new TH1F(histnamePt, histnamePt, nBinsPt, ptMin, ptMax);
        
        // Connect to the tree
        TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
        fAnaTree->Draw("fTrkTrkM>>"+histnameM, mCut, "hist");
        fAnaTree->Draw("fTrkTrkPt>>"+histnamePt, mCut, "hist");
        
        hM->SaveAs(histnameM+".root");
        hPt->SaveAs(histnamePt+".root");
        
        /*
        //Create a new workspace to manage the project
        RooWorkspace* wspace = new RooWorkspace("myJpsi");
        ImportDataSet(wspace, fAnaTree, mCut, mMin, mMax, ptMin, ptMax, rapMin, rapMax);
         */
        
    }
}
