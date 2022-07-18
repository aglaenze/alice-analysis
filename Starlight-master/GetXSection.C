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

//____________________________________________
TGraph* CreateTGraph( Int_t size, const Double_t* x, const Double_t* y)
{
	TGraph* tg = new TGraph();
	for( Int_t i = 0; i < size; ++i ) {tg->SetPoint( i, x[i], y[i] );}
	return tg;
}

//____________________________________________
TGraph* CreateTGraph( const std::vector<Double_t>& x, const std::vector<Double_t>& y)
{ return CreateTGraph( x.size(), &x[0], &y[0]); }

int GetXSection(string process, int year, string config, bool pbEmitter, Double_t mMin, Double_t mMax, Double_t rapMin = 2.5, Double_t rapMax = 4) {
	
	bool twok16 = (year==2016);
	
	
	//if (!twok16 && config == "Pb-p") {rapMin = 2.6; rapMax = 3.6;}
	
	Double_t step = 0.02;
	Double_t factor = 1;
	if (process == "kCohJpsiToMu" && config != "Pb-Pb") factor = 0.001;
	
	TString filename = Form("files/%s/output-%d-%s.txt", process.c_str(), year, config.c_str());
	
	// Output variables
	double dSigAve = 0.;	// dSigma(p+Pb)/dy averaged over y
	double dSigInt = 0.;	// dSigma(p+Pb)/dy integrated over y
	double phFluxAve = 0.;	// Photon flux averaged over y
	double sigGpAve = 0.;	// Sigma (gamma-p) averaged over y
	double sigGpInt = 0.;	// Sigma (gamma-p) divided by the photon flux and integrated over y
	
	double sigTotal = 0.;
	int rapMaxInput = 0;
	
	bool start = false;
	
	ifstream file(filename);
	if(file) {
		string line; //Une variable pour stocker les lignes lues
		int i = 0;
		double y, Wgp1, phFlux1, sigma2, dsig1, Wgp2, phFlux2, sigma1, dsig2, dsig;
		vector<double> yVec = {}, phFluxVec = {}, dsigVec = {}, wgpVec = {}, sigGpVec = {};
		int yMin = 0, yMax = 0;
		while(getline(file, line)) {
			if (!start){
				if (line.find("Total cross section") != string::npos) {
					stringstream stream0(line);
					string s1, s2, s3;
					stream0 >> s1 >> s2 >> s3 >> sigTotal;
				}
				else if (line.find("maximum absolute value for rapidity") != string::npos) {
					stringstream stream1(line);
					string s1, s2, s3, s4, s5, s6;
					stream1 >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> rapMaxInput;
					//cout << line << endl;
				}
				else if (line.find("Wgp") != string::npos) {
					start= true;
					/*
					 cout << "table starts here" << endl;
					 cout << line << endl;
					 */
				}
				continue;}
			
			if (line.find("+") == string::npos) break;
			//cout << line << endl;
			stringstream stream(line);
			stream >> y >> Wgp1 >> phFlux1 >> sigma2 >> dsig1 >> Wgp2 >> phFlux2 >> sigma1 >> dsig2 >> dsig;
			yVec.push_back(y);
			//cout << y << endl;
			if ( (config == "p-Pb" && pbEmitter) || (config == "Pb-p" && !pbEmitter)) {
				//if (config == "p-Pb") {
				phFluxVec.push_back(phFlux2);
				dsigVec.push_back(dsig2);
				wgpVec.push_back(Wgp2);
				sigGpVec.push_back(sigma1*1000);
			}
			else {
				phFluxVec.push_back(phFlux1);
				dsigVec.push_back(dsig1);
				wgpVec.push_back(Wgp1);
				sigGpVec.push_back(sigma2*1000);
			}
			if (y > rapMin-step && y < rapMin+step) yMin = i;
			else if (y > rapMax-step && y < rapMax+step) yMax = i;
			i++;
			
		} // end of while getline
		if (start) {
			TCanvas* cv = new TCanvas("cv", "cv", 800, 400);
			cv->Divide(2);
			cv->cd(1);
			TGraph* grFlux = CreateTGraph(yVec, phFluxVec);
			grFlux->SetTitle("Photon flux = f(y)");
			grFlux->GetXaxis()->SetRangeUser(-5., 5);
			grFlux->GetXaxis()->SetTitle("y");
			grFlux->GetYaxis()->SetTitle("k*dn/dk");
			grFlux->Draw("AP");
			cv->cd(2);
			TGraph* grDSig = CreateTGraph(yVec, dsigVec);
			grDSig->GetXaxis()->SetRangeUser(-5.,5);
			grDSig->SetTitle("dSigma/dy = f(y)");
			grDSig->GetXaxis()->SetTitle("y");
			grDSig->GetYaxis()->SetTitle("dSigma/dy");
			grDSig->Draw("AC");
			
			cv->SaveAs(Form("files/%s/Sigma-%d-%s.pdf", process.c_str(), year, config.c_str()));
			
			
			//cout << "yMin1 = " << yMin << endl;
			//cout << "yMax1 = " << yMax << endl;
			
			int numberOfPoints = 0.;
			for (int k = yMin; k < yMax; k++) {
				//phFluxInt[l] += phFluxVec[k]*(yVec[k+1]-yVec[k]);
				dSigAve += factor*dsigVec[k];
				phFluxAve += phFluxVec[k]/factor;
				dSigInt += factor*dsigVec[k]*(yVec[k+1]-yVec[k]);
				sigGpAve += sigGpVec[k];
				numberOfPoints++;
				
				sigGpInt += factor*dsigVec[k]*(yVec[k+1]-yVec[k])/phFluxVec[k];
			}
			phFluxAve /= numberOfPoints;
			dSigAve /= numberOfPoints;
			sigGpAve /= numberOfPoints;
			
			sigGpInt /= numberOfPoints/1000.;
			
			cout << endl;
			cout << rapMin << " < y < " << rapMax << endl;
			cout << min(wgpVec[yMin], wgpVec[yMax]) << " GeV < Wgp < " << max(wgpVec[yMin], wgpVec[yMax]) << " GeV" << endl;
			
			cout << "Photon Flux = " << phFluxAve << endl;
			cout << "dSigma/dy = " << dSigAve << " \u03BCb" << endl;
			cout << "Integrated Sigma = " << dSigInt << " \u03BCb" << endl;
			if (pbEmitter) {
				cout << "Sigma (gamma-p) = " << sigGpAve << " nb" << endl;
				//cout << "Sigma (gamma-p) = " << sigGpInt << " nb" << endl;
			}
			else {
				cout << "Sigma (gamma-Pb) = " << sigGpAve << " nb" << endl;
				//cout << "Sigma (gamma-Pb) = " << sigGpInt << " nb" << endl;
			}
		}// end of if (start)
	}
	else {cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl; return -1;}
	
	//TCut rapCut = Form("fTrkTrkY > %f && fTrkTrkY < %f", rapMin, rapMax;
	TCut rapCut = Form("!(fTrkTrkY < %f || fTrkTrkY > %f)", rapMin, rapMax);
	//TCut massCut = Form("fTrkTrkM>%f && fTrkTrkM<%f", mMin, mMax);
	TCut massCut = Form("!(fTrkTrkM<%f || fTrkTrkM>%f)", mMin, mMax);
	TCut cut = rapCut + massCut;
	TString rootFileName = Form("files/%s/tree-%d-%s.root", process.c_str(), year, config.c_str());
	//cout << rootFileName << endl;
	//return 0;
	TFile* f = TFile::Open(rootFileName, "READ");
	TTree* t = (TTree*)f->Get("fAnaTree");
	Int_t nEntries = t->GetEntries();
	TH1F* hist = new TH1F("histPt", "histPt", 1000, -10, 100);
	
	t->Draw("fTrkTrkPt>>histPt", cut);
	Int_t nCut = hist->GetEntries();
	
	Double_t p = (double)nCut/nEntries;
	Double_t sigErr = sigTotal * sqrt(p*(1-p)/nEntries);
	
	// Now print the results
	cout << endl;
	cout << "Starlight says Sigma = " << sigTotal << " \u03BCb for " << -rapMaxInput << " < Y < " << rapMaxInput << endl;
	
	cout << "nCut / nEntries = " << nCut << " / " << nEntries << " = " << (double)nCut/nEntries << endl;
	if (start) cout << "Ratio of Sigmas = " << dSigInt/sigTotal << endl;
	
	cout << endl;
	cout << "So it should be Sigma (" << config << ") = " << p*sigTotal << " ± " << sigErr << " \u03BCb for " << rapMin << " < Y < " << rapMax << " and " << mMin << " GeV/c2 < M < " << mMax << " GeV/c2" << endl;
	
	return 0;
	
}
