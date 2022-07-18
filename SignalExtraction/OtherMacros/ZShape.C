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


void ZShape() {
	gStyle->SetOptStat(0);
	
	Double_t mMin = 2.8;
	Double_t mMax = 3.3;
	vector<string> periods = {"LHC16r", "LHC16s"};
    
    bool excOnly = false;
	
	for (int i = 0; i<(int)periods.size(); i++) {
		string period = periods[i];
		
		std::list<TCut> mCutList = DefineCuts(period, excOnly);
		TCut mCut = "";
		for (std::list <TCut>::iterator iter = mCutList.begin(); iter != mCutList.end(); ++iter) {mCut += *iter;}
		//mCut += "fTrkTrkM>1.7 && fTrkTrkM<2.2";
		mCut += Form("fTrkTrkM>%f && fTrkTrkM<%f", mMin, mMax);
		mCut += "fTrkTrkPt>0 && fTrkTrkPt<3";
		
		// Get data file
		TFile *fAna = new TFile(Form("../../p-Pb-2016/rootFiles/std/AnalysisResults_%s.root", period.c_str()),"READ");
		TTree* t = (TTree*)fAna->Get("MyTask/fAnaTree");
        
        TH1F* hist1 = new TH1F("hist1", "hist1", 300, 0, 1.2);
        TH1F* hist2 = new TH1F("hist2", "hist2", 300, 0, 1.2);
		
		TCanvas* cv = new TCanvas("cv", "cv", 400, 600);
		cv->Divide(1,2);
		cv->cd(1);
        hist1->SetTitle("Z distribution for the selected data with 2.8 < M_{#mu#mu} < 3.3");
		t->Draw("fTrkTrkZ>>hist1", mCut);
		cv->cd(2);
        hist2->SetTitle("Z distribution for data with M_{#mu#mu} < 2.5");
		t->Draw("fTrkTrkZ>>hist2", "fTrkTrkM < 2.5");

        
		
		cv->SaveAs(Form("Plots/z-shape-%s.pdf", period.c_str()));
	}
	
}
