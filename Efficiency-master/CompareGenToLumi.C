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
#include <TChain.h>

// AliROOT
#include "AliCDBManager.h"
#include "AliTriggerScalers.h"
#include "AliTriggerRunScalers.h"
#include "AliTimeStamp.h"
#include "AliTriggerScalersRecord.h"
#include "AliTriggerConfiguration.h"
#include "AliLHCData.h"
#include "AliTriggerClass.h"
#include "AliTriggerBCMask.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"

#include "../Luminosity/LuminosityCalculationStandAlone.C"

using namespace std;


int DrawLumiRunByRun(string process, string period) {
	
	vector<Int_t> runList;
	TString className;
	string trg;
	if (period == "LHC16r") {
		runList = LHC16rRuns;
		className = "CMUP14-B-NOPF-MUFAST";
		trg = "MUP14";
	}
	else if (period == "LHC16s") {
		runList = LHC16sRuns;
		className = "CMUP23-B-NOPF-MUFAST";
		trg = "MUP23";
	}
	else {
		cout << "Wrong period" << endl;
		return -1;
	}
	
	TString path = "../p-Pb-2016/rootFiles/MC-std/";
	TString fileName = Form("%s_MC_%s.root", period.c_str(), process.c_str());
	TFile *f = new TFile(path+ "AnalysisResults_"+fileName,"read");
	
	bool draw = true;
	
	TTree* fAnaTree = (TTree*)f->Get("MyTask/fAnaTree");
	TTree* fGenTree = (TTree*)f->Get("MyTask/fGenTree");

	int nRuns = (int)runList.size();
	Int_t firstRun = runList[0];
	Int_t lastRun = runList[nRuns-1];
	
	// Create the histogram of efficiencies
	TH1F* hGen = new TH1F("hGen", "hGen", 3, 0, 3);
	hGen->SetStats(0);
	//hGen->SetFillColor(38);
	hGen->LabelsDeflate();
	hGen->SetLineColor(kRed);
	
	TH1F* hLumi = new TH1F("hLumi", "hLumi", 3, 0, 3);
	hLumi->SetStats(0);
	//hLumi->SetFillColor(38);
	hLumi->LabelsDeflate();
	hLumi->SetLineColor(kBlue);
	
	vector<ULong64_t > seenEvents = {};
	for (int k = 0; k<nRuns; k++) seenEvents.push_back(0);
	if (period == "LHC16r") SeenEvents(0, seenEvents);
	else if (period == "LHC16s") SeenEvents(1, seenEvents);
	// Here you put histogram with fired triggers per run in your dataset
	TH1D* hTrg = new TH1D(Form("h%strg",trg.c_str()),"",nRuns,0,nRuns);
	
	for (Int_t irun=0;irun<nRuns;irun++){
		// Randomly filling this histo for testing purposes
		//    hMUP6trg->Fill(irun,gRandom->Integer(1000000));
		Printf("seen %llu",seenEvents[irun]);
		hTrg->Fill(irun,seenEvents[irun]);
	}
	const Int_t nrunsmax = 1000;


	TChain* t = new TChain("trending");
	//t->AddFile("trending.root");
	t->AddFile("trending_merged.root");
	t->LoadTree(0);
	TObjArray* classes = new TObjArray();
	Double_t  lumi_seen[nrunsmax] = {0};
	Double_t  class_lumi[nrunsmax] = {0};
	Double_t  class_ds[nrunsmax] = {0};
	ULong64_t  class_l2a[nrunsmax] = {0};
	Int_t run;
	Double_t mu = 0;
	t->SetBranchAddress("mu",&mu);
	t->SetBranchAddress("run",&run);
	t->SetBranchAddress("classes",&classes);
	t->SetBranchAddress("lumi_seen",&lumi_seen);
	t->SetBranchAddress("class_lumi",&class_lumi);
	t->SetBranchAddress("class_ds",&class_ds);
	t->SetBranchAddress("class_l2a",&class_l2a);
	t->BuildIndex("run");
	
	vector<Int_t> missingRuns = {};
	TCanvas* c0 = new TCanvas();
	for (int i = 0; i<nRuns; i++) {
		Int_t r = runList[i];
		char* srun = Form("%i",r);
		cout << "Run " << srun << endl;
		TH1F* h2 = new TH1F("h2","", lastRun-firstRun+1, firstRun, lastRun);
		fGenTree->Draw("fRunNum>>h2", Form("fRunNum == %d", r));
		if (h2->GetEntries() == 0) {
			delete h2;
			missingRuns.push_back(r);
			continue;
		}

		// Lumi now
		t->GetEntryWithIndex(r);
		AliTriggerClass* cl = (AliTriggerClass*) classes->FindObject(className.Data());
		if (!cl) {cout << "class not found in trending.root" << endl; continue;}
		Int_t iclass = classes->IndexOf(cl);
		Printf("%i %i %s",run, iclass, cl->GetName());
		Double_t l2a = (Double_t) class_l2a[iclass];
		if (l2a == 0) continue;
		Printf("l2a a triggers %.llu %f ",class_l2a[iclass], hTrg->GetBinContent(i+1));
		hLumi->Fill(srun,class_lumi[iclass]);
		
		hGen->Fill( srun , (double)h2->GetEntries()/class_lumi[iclass]);
		delete h2;
	}
	
	cout << endl;
	cout << "Missing runs: " << endl;
	for (int k = 0; k<(int)missingRuns.size(); k++) {cout << missingRuns[k] << endl;}
	cout << endl;
	
	TLegend* lgd = new TLegend(0.5, 0.7, 0.9, 0.9);
	lgd->AddEntry(hGen, "Generated events (MC) / Luminosity", "l");
	//lgd->AddEntry(hLumi, "Luminosity", "l");
	
	// Draw it in a canvas
	TCanvas* cv = new TCanvas("cv", "cv", 800, 300);
	hGen->Draw("hist");
	//hLumi->Draw("hist same");
	lgd->Draw();
	cv->SaveAs(Form("Plots/GeneratedVsLumi-%s-%s.pdf", period.c_str(), process.c_str()));
	
	return 0;
}

int CompareGenToLumi(string process, vector<string> periods = {"LHC16r", "LHC16s"}) {

	for (int l = 0; l < (int)periods.size(); l++) {
		string period = periods[l];
		cout << endl;
		cout << "Going through period " << period << " and process " << process << endl;
		DrawLumiRunByRun(process, period);
	}
	return 0;
}


