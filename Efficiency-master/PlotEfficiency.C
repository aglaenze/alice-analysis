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


// LHC16r
const vector<Int_t> LHC16rRuns = {265594, 265596, 265607, 265691, 265694, 265697, 265698, 265700, 265701, 265709, 265713, 265714, 265740, 265741, 265742, 265744, 265746, 265754, 265756, 265785, 265787, 265788, 265789, 265792, 265795, 265797, 265840, 266022, 266023, 266025, 266034, 266074, 266076, 266081, 266084, 266085, 266086, 266117, 266187, 266189, 266190, 266193, 266196, 266197, 266208, 266234, 266235, 266296, 266299, 266300, 266304, 266305, 266312, 266316, 266318};

// LHC16s
const vector<Int_t> LHC16sRuns  = {266439, 266441, 266472, 266479, 266480, 266487, 266514, 266516, 266518, 266520, 266522, 266523, 266525, 266533, 266534, 266539, 266543, 266549, 266584, 266587, 266588, 266591, 266593, 266595, 266613, 266614, 266618, 266621, 266630, 266657, 266658, 266659, 266665, 266668, 266669, 266674, 266676, 266702, 266703, 266706, 266708, 266775, 266776, 266800, 266805, 266807, 266857, 266878, 266880, 266882, 266883, 266885, 266886, 266912, 266915, 266940, 266942, 266943, 266944, 266988, 266993, 266994, 266997, 266998, 267020, 267022, 267062, 267063, 267067, 267070, 267072, 267077, 267109, 267110, 267130, 267131};

bool contains(std::map<Int_t, Int_t> &myMap, const Int_t &element)
{
	for(std::map<Int_t, Int_t>::iterator iter = myMap.begin(); iter != myMap.end(); ++iter)
	{
		Int_t k =  iter->first;
		if (element == iter->first) return true;
	}
	return false;
}

int GetTotalEff(string process, string period, TCut cut, TCut cutMc) {
	
	TString path = "../p-Pb-2016/rootFiles/MC-std/";
	TString fileName = Form("%s_MC_%s.root", period.c_str(), process.c_str());
	TFile *f = new TFile(path+ "AnalysisResults_"+fileName,"read");
	
	bool draw = true;
	
	TTree* fAnaTree = (TTree*)f->Get("MyTask/fAnaTree");
	TTree* fGenTree = (TTree*)f->Get("MyTask/fGenTree");
	
	
	TCanvas* c0 = new TCanvas();	// To avoid the warning
	
	TH1F* h1 = new TH1F("h1","", 100, 0, 100);
	TH1F* h2 = new TH1F("h2","", 100, 0, 100);
	fAnaTree->Draw("fTrkTrkM>>h1", cut);
	fGenTree->Draw("fMCTrkTrkM>>h2", cutMc);
	
	Int_t num = h1->GetEntries();
	Int_t denom = h2->GetEntries();
	/*
	 cout << h1->GetEntries() << endl;
	 cout << h2->GetEntries() << endl;
	 */
	delete h1;
	delete h2;
	
	
	Double_t eff = (double)num/denom*100;
	
	cout << endl;
	cout << "Efficiency = " << eff << "%" << endl;
	cout << "Reconstructed events that passed the cuts = " << num << endl;
	cout << "Events generated within the cuts = " << denom << endl;
	
	return 0;
}


int DrawEffRunByRun(string process, string period, TCut cut, TCut cutMc, Double_t mMin, Double_t mMax) {
	
	TString path = "../p-Pb-2016/rootFiles/MC-std/";
	TString fileName = Form("%s_MC_%s.root", period.c_str(), process.c_str());
	TFile *f = new TFile(path+ "AnalysisResults_"+fileName,"read");
	
	bool draw = true;
	
	TTree* fAnaTree = (TTree*)f->Get("MyTask/fAnaTree");
	TTree* fGenTree = (TTree*)f->Get("MyTask/fGenTree");
	
	vector<Int_t> runList;
	if (period == "LHC16r") runList = LHC16rRuns;
	else if (period == "LHC16s") runList = LHC16sRuns;
	else {cout << "Wrong period" << endl; return -1;}
	
	int nRuns = (int)runList.size();
	Int_t firstRun = runList[0];
	Int_t lastRun = runList[nRuns-1];
	
	// Create the histogram of efficiencies
	TH1F* hEff = new TH1F("hEff", "hEff", 3, 0, 3);
	hEff->SetStats(0);
	hEff->SetFillColor(38);
	hEff->LabelsDeflate();
	
	vector<Int_t> missingRuns = {};
	Int_t num = 0;
	Int_t denom = 0;
	TCanvas* c0 = new TCanvas();
	for (int k = 0; k<nRuns; k++) {
		cout << "Run " << runList[k] << endl;
		TH1F* h1 = new TH1F("h1","", lastRun-firstRun+1, firstRun, lastRun);
		TH1F* h2 = new TH1F("h2","", lastRun-firstRun+1, firstRun, lastRun);
		fAnaTree->Draw("fRunNum>>h1", cut+ Form("fRunNum == %d", runList[k]));
		fGenTree->Draw("fRunNum>>h2", cutMc+ Form("fRunNum == %d", runList[k]));
		if (h2->GetEntries() == 0) {
			delete h1;
			delete h2;
			missingRuns.push_back(runList[k]);
			continue;
		}
		hEff->Fill( Form("%d", runList[k]) , (double)h1->GetEntries()/h2->GetEntries());
		num += h1->GetEntries();
		denom += h2->GetEntries();
		delete h1;
		delete h2;
	}
	
	cout << endl;
	if ((int)missingRuns.size()>0) {
		cout << "Missing runs: " << endl;
		for (int k = 0; k<(int)missingRuns.size(); k++) {cout << missingRuns[k] << endl;}
	}
	else {cout << "No missing runs" << endl;}
	cout << endl;
	
	
	TLegend* lgd = new TLegend(0.6, 0.7, 0.9, 0.9);
	string percent = "%";
	lgd->AddEntry(hEff, Form("Efficiency = %.3f %s (%.1f < M < %.1f)", (double)num/denom*100, percent.c_str(), mMin, mMax ), "l");
	
	
	// Draw it in a canvas
	TCanvas* cv = new TCanvas("cv", "cv", 800, 300);
	hEff->Draw("hist");
	lgd->Draw();
	cv->SaveAs(Form("Plots/Efficiency-%s-%s.pdf", period.c_str(), process.c_str()));
	
	return 0;
}

int PlotEfficiency(string process, vector<string> periods = {"LHC16r", "LHC16s"}, Double_t mMin = 2.9, Double_t mMax = 3.2, Double_t rapMin = 2.5, Double_t rapMax = 4.0, bool drawRunByRun = true) {
	
	// Cuts for generated MC
	TCut rapCutMc = Form("fMCTrkTrkY < -%f && fMCTrkTrkY > -%f", rapMin, rapMax);
	TCut unlikeSignCutMc = "(fMCTrkQ1 < 0 && fMCTrkQ2 > 0) || (fMCTrkQ1 > 0 && fMCTrkQ2 < 0)";
	// Don't use the mass cut in generated events
	TCut massCutMc = Form("fMCTrkTrkM < %f && fMCTrkTrkM > %f", mMax, mMin);
	TCut cutMc = rapCutMc + unlikeSignCutMc;
	
	// Cuts for reconstructed
	TCut rapCut = Form("fTrkTrkY < -%f && fTrkTrkY > -%f", rapMin, rapMax);
	TCut unlikeSignCut = "(fTrkQ1 < 0 && fTrkQ2 > 0) || (fTrkQ1 > 0 && fTrkQ2 < 0)";
	TCut triggerCut;	// defined in the for loop
	TCut massCut = Form("fTrkTrkM < %f && fTrkTrkM > %f", mMax, mMin);
	
	for (int l = 0; l < (int)periods.size(); l++) {
		string period = periods[l];
		if (period == "LHC16r") triggerCut = "0<(fL0inputs&(1<<5))";
		else if (period == "LHC16s") triggerCut = "0<(fL0inputs&(1<<13))";
		else {cout << "problem" << endl; return -1;}
		TCut cut = rapCut + unlikeSignCut + massCut + triggerCut;
		cout << endl;
		cout << "Going through period " << period << " and process " << process << endl;
		if (drawRunByRun) DrawEffRunByRun(process, period, cut, cutMc, mMin, mMax);
		GetTotalEff(process, periods[l], cut, cutMc);
	}
	return 0;
}


