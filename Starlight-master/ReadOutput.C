#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <sstream>
#include <cmath>

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
#include "TGraph.h"

using namespace std;


void TrkTrkKine(Int_t gpid[2], Double_t px[2], Double_t py[2], Double_t pz[2], Double_t& fTrkTrkPt, Double_t& fTrkTrkPhi, Double_t& fTrkTrkY, Double_t& fTrkTrkM, Double_t& fTrkPt1, Double_t& fTrkPt2, Double_t& fTrkEta1, Double_t& fTrkEta2, Double_t& fTrkPhi1, Double_t& fTrkPhi2, Double_t& fTrkQ1, Double_t& fTrkQ2) {
	
	Double_t MuonMass = 0.105658389; // GeV/c^2
	
	Int_t posIdx, negIdx;
	
	if (gpid[0] == 5 &&  gpid[1] == 6 ) {posIdx = 0; negIdx = 1;}
	else if (gpid[0] == 6 &&  gpid[1] == 5 ) {posIdx = 1; negIdx = 0;}
	else {cout << "Problem in kinematics" << endl; return;}
	
	
	// --  positive track
	TLorentzVector LV1;
	Double_t e[2];
	for (int k = 0; k<2; k++) {e[k] = sqrt( pow(px[k],2) + pow(py[k],2) + pow(pz[k],2) + pow(MuonMass,2));}
	LV1.SetPxPyPzE(px[posIdx], py[posIdx], pz[posIdx], e[posIdx]);
	//LV1.SetPtEtaPhiM(Track1->Pt(), Track1->Eta(), Track1->Phi(), TrkMass);
	// --  negative track
	TLorentzVector LV2;
	LV2.SetPxPyPzE(px[negIdx], py[negIdx], pz[negIdx], e[negIdx]);
	
	// vector of Trk+Trk
	TLorentzVector TrkTrk = LV1+LV2;
	
	// get tree variables
	fTrkTrkPt = TrkTrk.Pt();
	fTrkTrkPhi = TrkTrk.Phi();
	fTrkTrkY = TrkTrk.Rapidity();
	fTrkTrkM = TrkTrk.M();
	fTrkPt1 = LV1.Pt();
	fTrkPt2 = LV2.Pt();
	fTrkEta1 = LV1.Eta();
	fTrkEta2 = LV2.Eta();
	fTrkPhi1 = LV1.Phi();
	fTrkPhi2 = LV2.Phi();
	fTrkQ1 = +1;
	fTrkQ2 = -1;
}



int ReadOutput(string process, int year, string config) {
	
	TFile* fOut = new TFile(Form("files/%s/tree-%d-%s.root", process.c_str(), year, config.c_str()), "RECREATE");
	TTree* t = new TTree("fAnaTree", "Kinematics");
	Double_t fTrkTrkPt, fTrkTrkPhi, fTrkTrkY, fTrkTrkM, fTrkPt1, fTrkPt2, fTrkEta1, fTrkEta2, fTrkPhi1, fTrkPhi2, fTrkQ1, fTrkQ2;
	t->Branch("fTrkTrkPt", &fTrkTrkPt, "fTrkTrkPt/D");
	t->Branch("fTrkTrkPhi", &fTrkTrkPhi, "fTrkTrkPhi/D");
	t->Branch("fTrkTrkY", &fTrkTrkY, "fTrkTrkY/D");
	t->Branch("fTrkTrkM", &fTrkTrkM, "fTrkTrkM/D");
	t->Branch("fTrkPt1", &fTrkPt1, "fTrkPt1/D");
	t->Branch("fTrkPt2", &fTrkPt2, "fTrkPt2/D");
	t->Branch("fTrkEta1", &fTrkEta1, "fTrkEta1/D");
	t->Branch("fTrkEta2", &fTrkEta2, "fTrkEta2/D");
	t->Branch("fTrkPhi1", &fTrkPhi1, "fTrkPhi1/D");
	t->Branch("fTrkPhi2", &fTrkPhi2, "fTrkPhi2/D");
	t->Branch("fTrkQ1", &fTrkQ1, "fTrkQ1/D");
	t->Branch("fTrkQ2", &fTrkQ2, "fTrkQ2/D");
	
	TString filename = Form("files/%s/slight-%d-%s.out", process.c_str(), year, config.c_str());
	ifstream file(filename);
	if(file) {
		string line; //Une variable pour stocker les lignes lues
		int nVertex = 0;
		int gpid[2], nev[2];
		double px[2], py[2], pz[2];
		string firstWord;
		while(getline(file, line)) {
			/*
			if (line.find("VERTEX:") != string::npos) {
				//cout << nVertex << endl;
				nVertex++;
				
			}
			 */
			if (line.find("VERTEX:") == string::npos) {continue;}
			//cout << line << endl;
			nVertex++;
			// Get first track line
			getline(file, line);
			if (line.find("TRACK:") == string::npos) {
				cout << line << endl;
				cout << "Problem TRACK1" << endl;
				return -1;
			}
			stringstream stream(line);
			stream >> firstWord >> gpid[0] >> px[0] >> py[0] >> pz[0] >> nev[0];
			
			// Get second track line
			getline(file, line);
			if (line.find("TRACK:") == string::npos) {
				cout << line << endl;
				cout << "Problem TRACK2" << endl;
				return -1;
			}
			stringstream stream2(line);
			stream2 >> firstWord >> gpid[1] >> px[1] >> py[1] >> pz[1] >> nev[1];
			if (nev[0] != nev[1]) {
				cout << "not same event ID = " << nev[0] << " and " << nev[1] << ", return" << endl;
				return -1;}
			
			// Now compute kinematics
			TrkTrkKine(gpid, px, py, pz, fTrkTrkPt, fTrkTrkPhi, fTrkTrkY, fTrkTrkM, fTrkPt1, fTrkPt2, fTrkEta1, fTrkEta2, fTrkPhi1, fTrkPhi2, fTrkQ1, fTrkQ2);
			t->Fill();
			
		}
		//cout << line << endl;
		cout << "nVertex = " << nVertex << endl;
	}
	
	else {cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl; return -1;}
	
	fOut->cd();
	t->Write();
	fOut->Close();
	
	
	return 0;
	
	
}
