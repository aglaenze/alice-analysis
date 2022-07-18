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
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooLandau.h"
#include "RooHistPdf.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"

#include "ExtendedCrystalBall.h"
#include "Utils.C"


using namespace RooFit;
using namespace std;



//RooDataHist* GetPtHist(string period, TString process){
void GetPtHistMC(RooWorkspace* ws, string rootfilePathMC, string period, TString process, Double_t ptMin = 0, Double_t ptMax = 4., TCut cutMc = ""){
	TFile *fSimu = new TFile(Form("%s/AnalysisResults_%s_MC_", rootfilePathMC.c_str(), period.c_str()) + process + ".root","READ");
	TTree* fAnaTree = (TTree*)fSimu->Get("MyTask/fAnaTree");
	
	RooRealVar pt = *ws->var("fTrkTrkPt");
	Int_t nBins = int((ptMax-ptMin)*100);
	TH1F* histPt = new TH1F("hPt"+process, "hPt"+process, nBins, ptMin, ptMax);
	fAnaTree->Draw("fTrkTrkPt>>hPt"+process, cutMc);
	RooDataHist* dTemplatePt = new RooDataHist("ptHist"+process,"ptHist"+process, RooArgList(pt),histPt);
	RooHistPdf* ptPdf = new RooHistPdf("pt"+process, "pt"+process, pt, *dTemplatePt);
	ws->import(*ptPdf);
}

void GetMHistMC(RooWorkspace* ws, string rootfilePathMC, string period, TString process, Double_t mMin = 2.5, Double_t mMax = 3.5, TCut cutMc = ""){
	TFile *fSimu = new TFile(Form("%s/AnalysisResults_%s_MC_", rootfilePathMC.c_str(), period.c_str()) + process + ".root","READ");
	TTree* fAnaTree = (TTree*)fSimu->Get("MyTask/fAnaTree");
	
	RooRealVar m = *ws->var("fTrkTrkM");
	
	Int_t nBins = int((mMax-mMin)*100);
	TH1F* histM = new TH1F("hM"+process, "hM"+process, nBins, mMin, mMax);
	fAnaTree->Draw("fTrkTrkM>>hM"+process, cutMc);
	RooDataHist* dTemplateM = new RooDataHist("mHist"+process,"mHist"+process, RooArgList(m),histM);
	RooHistPdf* mPdf = new RooHistPdf("m"+process, "m"+process, m, *dTemplateM);
	ws->import(*mPdf);
}


void GetV0Template(RooWorkspace* ws, string rootfilePath, string period) {
	
	// Open the file
	TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()),"READ");
	
	// Connect to the tree
	TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
	
	TCut newCut = "fV0CBBNHits>7 && fTrkTrkM<3";
	if (period == "LHC16s") newCut = "fADABBNHits>1 && fTrkTrkM<3";
	Double_t ptMin = 0, ptMax = 8;
	
	Int_t nBins = int((ptMax-ptMin)*10);
	TH1F* hist = new TH1F("hV0C7", "hV0C7", nBins, ptMin, ptMax);
	fAnaTree->Draw("fTrkTrkPt>>hV0C7", newCut);
	
	RooRealVar pt = *ws->var("fTrkTrkPt");
	RooDataHist* dTemplatePt = new RooDataHist("hV0C7","hV0C7", RooArgList(pt),hist);
	RooHistPdf* ptPdf = new RooHistPdf("ptV0C7", "ptV0C7", pt, *dTemplatePt);
	
	ws->import(*ptPdf);
	//gStyle->SetOptStat(0);
}


void LoadMassFitFunctions(RooWorkspace* ws, string period) {
	
	RooRealVar m = *ws->var("fTrkTrkM");
	RooRealVar pt = *ws->var("fTrkTrkPt");
	RooDataSet* data = (RooDataSet*) ws->data("data");
	
	// First mass PDFs
	// J/Psi peak
	// Take tails parameters from TailParameters.C
	double alphaL = 0.961839, nL = 7.521515, alphaR = 2.641260, nR = 3.325886;
	if (period == "LHC16s") {alphaL = 0.993482; nL = 6.845735; alphaR = 2.669157; nR = 3.078395;}
	RooRealVar mean_jpsi("mean_jpsi","mean_jpsi",3,3.0,3.3);
	RooRealVar sigma_jpsi("sigma_jpsi","sigma_jpsi",0.0811359, 0.05, 0.15);
	//RooRealVar sigma_jpsi("sigma_jpsi","sigma_jpsi",0.081);
	RooRealVar alpha_jpsi_L("alpha_jpsi_L","alpha_jpsi_L", alphaL);
	RooRealVar n_jpsi_L("n_jpsi_L","n_jpsi_L", nL);
	RooRealVar alpha_jpsi_R("alpha_jpsi_R","alpha_jpsi_R", alphaR);
	RooRealVar n_jpsi_R("n_jpsi_R","n_jpsi_R", nR);
	//RooCBShape *jpsi = new RooCBShape("jpsi","crystal ball PDF", m, mean_jpsi, sigma_jpsi,alpha_jpsi,n_jpsi);
	ExtendedCrystalBall *jpsi = new ExtendedCrystalBall("jpsi","crystal ball PDF", m,
														mean_jpsi, sigma_jpsi, alpha_jpsi_L,
														n_jpsi_L, alpha_jpsi_R, n_jpsi_R);
	
	ws->import(*jpsi);
	
	// Then Psi(2s)
	// Most of the time not needed
	// That is kIncohPsi2sToMu in MC data
	double alphaL2 = 1.200001, nL2 = 3.017759, alphaR2 = 2.928444, nR2 = 2.256593;
	if (period == "LHC16s") {alphaL2 = 1.192599; nL2 = 3.118819; alphaR2 = 2.927051; nR2 = 2.055075;}
	
	Double_t factorMean = 1.1902; // maybe to adjust
								  //RooRealVar sigma_psi("sigma_psi","sigma_psi",0.07,0,0.1);
								  // scaling factor Psi(2S) / JPsi
	Double_t factorSigma = 1.05; // maybe to adjust
	RooRealVar mean_scaling("mean_scaling", "", factorMean);
	RooRealVar sigma_scaling("sigma_scaling", "", factorSigma);
	RooFormulaVar mean_psi("mean_psi","mean_jpsi*mean_scaling",RooArgSet(mean_jpsi, mean_scaling));
	RooFormulaVar sigma_psi("sigma_psi", "sigma_jpsi*sigma_scaling", RooArgSet(sigma_jpsi, sigma_scaling));
	//RooRealVar sigma_psi("sigma_psi","sigma_psi",0.07, 0, 0.15);
	RooRealVar alpha_psi_L("alpha_psi_L","alpha_psi_L",alphaL2);
	RooRealVar n_psi_L("n_psi_L","n_psi_L",nL2);
	RooRealVar alpha_psi_R("alpha_psi_R","alpha_psi_R",alphaR2);
	RooRealVar n_psi_R("n_psi_R","n_psi_R",nR2);
	//RooCBShape *psi = new RooCBShape("psi","crystal ball PDF",m,mean_psi,sigma_psi,alpha_psi,n_psi);
	ExtendedCrystalBall *psi2s = new ExtendedCrystalBall("psi2s","crystal ball PDF", m,
													   mean_psi, sigma_psi, alpha_psi_L,
													   n_psi_L, alpha_psi_R, n_psi_R);
	ws->import(*psi2s);
	
	// Finally background
	
	//Exponential background for mass
	RooRealVar a1("a1","a1",-1,-6,0.1);
	RooExponential* bkg = new RooExponential("bkg","bkg",m,a1);
	/*
	 // using mass template from MC data
	 GetMHistMC(ws, rootfilePathMC, period, "kTwoGammaToMuLow", mMin, mMax);
	 RooAbsPdf* bkg = ws->pdf("mkTwoGammaToMuLow");
	 bkg->SetName("exp");
	 */
	ws->import(*bkg);
}

void LoadJpsiPtFitFunctions(RooWorkspace* ws, string rootfilePathMC, string period, Double_t ptMin, Double_t ptMax, bool exp, bool excOnly) {
	
	RooRealVar m = *ws->var("fTrkTrkM");
	RooRealVar pt = *ws->var("fTrkTrkPt");
	RooDataSet* data = (RooDataSet*) ws->data("data");

	// Contribution from exclusive J/Psis
	
	// Using H1 formula (b is free or not)
	double bExcValue = 4;
	double unused;
	if (! GetParameters(period, bExcValue, unused) ) {cout << "Parameters not loaded" << endl; return;}
	
	RooRealVar bExc("bExc","bExc", bExcValue, 2.8, 7);
	if (!excOnly) {bExc.setVal(bExcValue); bExc.setConstant();}
	
	// H1 formula
	RooGenericPdf *ptJpsiExclusive = new RooGenericPdf("ptJpsiExclusive","exclusive jPsi PDF","(2*fTrkTrkPt*exp(-bExc*(fTrkTrkPt**2)))",RooArgSet(pt,bExc)) ;
	
	/*
	 // using template pre-defined in sPlot
	 TFile* fTemplates = new TFile(Form("%s/sPlotTemplates-%s.root", rootfilePath.c_str(), period.c_str()),"READ");
	 TH1F* hPtExclusive = (TH1F*)fTemplates->Get("hPtExclusive");
	 RooDataHist* ptHistExclusive = new RooDataHist("ptHistData","ptHistData", RooArgList(pt),hPtExclusive);
	 RooHistPdf* ptJpsiExclusive = new RooHistPdf("ptJpsiExclusive", "ptExclusive", pt, *ptHistExclusive);
	 */
	/*
	 // using pt template from MC data
	 GetPtHistMC(ws, rootfilePathMC, period, "kIncohJpsiToMu");
	 RooAbsPdf* ptJpsiExclusive = ws->pdf("ptkIncohJpsiToMu");
	 ptJpsiExclusive->SetName("ptJpsiExclusive");
	 */
	ws->import(*ptJpsiExclusive);
	
	// JPsi Dissociative
	/*
	 // using template pre-defined in sPlot
	 TH1F* hptJpsiDissociative = (TH1F*)fTemplates->Get("hptJpsiDissociative");
	 RooDataHist* ptHistDissociative = new RooDataHist("ptHistData","ptHistData", RooArgList(pt),hptJpsiDissociative);
	 RooHistPdf* ptJpsiDissociative = new RooHistPdf("ptJpsiDissociative", "ptJpsiDissociative", pt, *ptHistDissociative);
	 */
	RooGenericPdf *ptJpsiDissociative = nullptr;
	RooRealVar* bDiss = nullptr;
	RooRealVar* nDiss = nullptr;
	if (exp) {
		// H1 formula (the first one, with the exp)
		//RooRealVar bDiss("bDiss","bDiss", 0.323027, 0, 2);
		bDiss = new RooRealVar("bDiss","bDiss", 0.323027, 0.23, 1.8);
		ptJpsiDissociative = new RooGenericPdf("ptJpsiDissociative","Dissociative jPsi PDF","(2*fTrkTrkPt*exp(-bDiss*(fTrkTrkPt**2)))", RooArgSet(pt,*bDiss)) ;
	}
	else {
		// H1 formula (the second one, with the power law)
		bDiss = new RooRealVar("bDiss","bDiss", 1.7, 0.5, 4.);
		nDiss = new RooRealVar("nDiss","nDiss", 3.5, 1, 8);
        /*
		if (period == "LHC16r") {
            bDiss->setVal(1.7);
            nDiss->setVal(3.58);
            bDiss->setConstant();
            nDiss->setConstant();
		}
         */
        if (period == "LHC16s") {bDiss->setRange(1, 4);}

		ptJpsiDissociative = new RooGenericPdf("ptJpsiDissociative","Dissociative jPsi PDF","(2*fTrkTrkPt*(1.+(fTrkTrkPt**2)*(bDiss/nDiss))**(-nDiss))", RooArgSet(pt, *nDiss, *bDiss)) ;
	}

	/*
	if (period == "LHC16r") {
		bDiss->setVal(0.73); bDiss->setConstant();
		nDiss->setVal(3.35); nDiss->setConstant();
	}
	 */

	ws->import(*ptJpsiDissociative);
	
	// Contribution from gamma Pb events
	// using pt template from MC data
	GetPtHistMC(ws, rootfilePathMC, period, "kCohJpsiToMu", ptMin, ptMax);
	RooAbsPdf* ptJpsiGammaPb = ws->pdf("ptkCohJpsiToMu");
	ptJpsiGammaPb->SetName("ptJpsiGammaPb");
	ws->import(*ptJpsiGammaPb);
	
	// Inclusive events
	/*
	 // From data with additional cuts
	 GetInclusiveTemplate(ws, rootfilePath, period, mMin, mMax, ptMin, ptMax);
	 RooAbsPdf* ptJpsiInclusive = ws->pdf("ptJpsiInclusive");
	 ptJpsiInclusive->SetName("ptJpsiInclusive");
	 */
	// With a new formula


	RooRealVar *pt0 = new RooRealVar("pt0","pt0", 2.5, 2, 6);
	RooRealVar *nInc = new RooRealVar("nInc","nInc", 3.5, 2, 10);
    
    //if (!excOnly && period == "LHC16r") {pt0->setVal(2.2); pt0->setConstant();}

    /*
	 RooRealVar *pt0 = new RooRealVar("pt0","pt0", 4.5);
	 RooRealVar *nInc = new RooRealVar("nInc","nInc", 3.5);
     */

	RooGenericPdf* ptJpsiInclusive = new RooGenericPdf("ptJpsiInclusive","Inclusive jPsi PDF","fTrkTrkPt/((1.+(fTrkTrkPt/pt0)**2)**nInc)", RooArgSet(pt, *pt0, *nInc)) ;
    
    /*
    // From HERA
    TFile* f = new TFile("../p-Pb-2016/rootFiles/HEPData-ins586977-v1-root.root", "READ");
    TH1D* histInclusifs2 = (TH1D*)f->Get("Table 4/Hist1D_y1");
    TH1D* histInclusifs = (TH1D*)histInclusifs2->Clone();
    histInclusifs->SetXTitle("p_{T}");
    histInclusifs->SetYTitle("");
    histInclusifs->SetTitle("");
    for (int bin=0;bin<=histInclusifs->GetNcells();++bin) {
        histInclusifs->SetBinContent(bin, sqrt(histInclusifs2->GetBinContent(bin) ));
    }
    TH1D hist = *histInclusifs;
    f->Close();
    TCanvas* cv = new TCanvas();
    hist.Draw("hist");
    cv->SaveAs("Plots/test-inclusive-template.pdf");
    RooDataHist* ptHistInclusifs = new RooDataHist("ptHistData","ptHistData", RooArgList(pt), &hist);
    RooHistPdf* ptJpsiInclusive = new RooHistPdf("ptJpsiDissociative", "ptJpsiDissociative", pt, *ptHistInclusifs);
     */
	
	ws->import(*ptJpsiInclusive);
	
	// Psi(2s)
	GetPtHistMC(ws, rootfilePathMC, period, "kIncohPsi2sToMu");
	RooAbsPdf* ptPsi2s = ws->pdf("ptkIncohPsi2sToMu");
	ptPsi2s->SetName("ptPsi2s");
	ws->import(*ptPsi2s);
}

void LoadBkgPtFitFunctions(RooWorkspace* ws, string rootfilePath, string rootfilePathMC, string period, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, bool excOnly) {
	
	RooRealVar m = *ws->var("fTrkTrkM");
	RooRealVar pt = *ws->var("fTrkTrkPt");
	RooDataSet* data = (RooDataSet*) ws->data("data");


    /*
	// Background
	// Using MC: gamma gamma to Mu Mu
	TCut cutMc = Form("fTrkTrkM > %f && fTrkTrkM < %f", mMin, mMax);
	GetPtHistMC(ws, rootfilePathMC, period, "kTwoGammaToMuLow", ptMin, ptMax, cutMc);
	RooAbsPdf* ptTwoGamma = ws->pdf("ptkTwoGammaToMuLow");
	ptTwoGamma->SetName("ptTwoGamma");
     */

    // Landau
    RooRealVar mean_gg("mean_gg","mean_gg",0.2, 0.003,0.7);
    RooRealVar sigma_gg("sigma_gg","sigma_gg", 1, 0.01,1.);
    // m.setConstant(kTRUE);
    // k.setConstant(kTRUE);
    RooLandau* ptTwoGamma = new RooLandau("ptTwoGamma","ptTwoGamma", pt, mean_gg, sigma_gg);
    
    ws->import(*ptTwoGamma);
    
    /*
     // Bkg from sidebands
     TFile* fData = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()),"READ");
     TTree* t = (TTree*)fData->Get("MyTask/fAnaTree");
     TH1F* hPtBackground = new TH1F("h", "h", int((ptMax-ptMin)*10), ptMin, ptMax);
     TCut bkgCut = "fTrkTrkM > 2 && fTrkTrkM < 2.8 && fTrkTrkPt > 0.2";
     std::list<TCut> cutList = DefineCuts(period, excOnly);
     for (std::list <TCut>::iterator iter = cutList.begin(); iter != cutList.end(); ++iter) {bkgCut += *iter;}
     t->Draw("fTrkTrkPt>>h", bkgCut);
     */
	

	// Extra background: pt distribution from background obtained with sPlot
	string suffix = "";
	if (excOnly) suffix = "-exclusive-only";
	TFile* fTemplates = new TFile(Form("%s/sPlotTemplates-%s%s.root", rootfilePath.c_str(), period.c_str(), suffix.c_str()),"READ");
	//TH1F* hPtBackground = (TH1F*)fTemplates->Get("hSubtraction");
	TH1F* hPtBackground = (TH1F*)fTemplates->Get("hSubSmooth");
    //TH1F* hPtBackground = (TH1F*)fTemplates->Get("hPtBackground__fTrkTrkPt");
     RooDataHist* ptHistBackground = new RooDataHist("ptHistData","ptHistData", RooArgList(pt), hPtBackground);
     RooHistPdf* ptBackground = new RooHistPdf("ptBackground", "ptBackground", pt, *ptHistBackground);
    
    /*
    // Extra background: function
    RooRealVar *pt0 = new RooRealVar("pt0_bkg","pt0_bkg", 0.7, 0.3, 5);
    RooRealVar *nInc = new RooRealVar("nInc_bkg","nInc_bkg", 3.5, 1, 6);

    RooGenericPdf* ptBackground = new RooGenericPdf("ptBackground","ptBackground","fTrkTrkPt/((1.+(fTrkTrkPt/pt0_bkg)**2)**nInc_bkg)", RooArgSet(pt, *pt0, *nInc)) ;
     */

	ws->import(*ptBackground);
	
	/*
	 // using sidebands
	 GetSidebandsTemplate(ws, rootfilePath, period, mMin, mMax, ptMin, ptMax);
	 RooAbsPdf* ptBackground = ws->pdf("ptSidebands");
	 ptBackground->SetName("ptBackground");
	 */
	
}


void LoadYields(RooWorkspace* ws, string period, Double_t rapMin, Double_t rapMax, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, bool excOnly) {
	
	double yieldGammaPbValue = 40;
	double unused;
	if (! GetParameters(period, unused, yieldGammaPbValue) ) {cout << "Parameters not loaded" << endl; return;}
	
	RooRealVar yieldJpsiExclusive("yieldJpsiExclusive","yieldJpsiExclusive",1200,200.,1700);
	RooRealVar yieldJpsiDissociative("yieldJpsiDissociative","yieldJpsiDissociative",1400,0.,1900);
	//RooRealVar yieldJpsiGammaPb("yieldJpsiGammaPb","yieldJpsiGammaPb", yieldGammaPbValue, yieldGammaPbValue*0.91, yieldGammaPbValue*1.09);
	RooRealVar yieldJpsiGammaPb("yieldJpsiGammaPb","yieldJpsiGammaPb", yieldGammaPbValue, yieldGammaPbValue*0, yieldGammaPbValue*1.09);
	RooRealVar yieldJpsiInclusive("yieldJpsiInclusive","yieldJpsiInclusive",800,0.,1600);
	
	RooRealVar yieldPsi2s("yieldPsi2s","yieldPsi2s",0,0.,1000);
	RooRealVar yieldTwoGamma("yieldTwoGamma","yieldTwoGamma",200,0.,300);
	RooRealVar yieldBkg("yieldBkg","yieldBkg",300,0.,1200);
	if (period == "LHC16r") {
		yieldJpsiExclusive.setRange(600*(rapMax-rapMin), 1500*(rapMax-rapMin));
		yieldJpsiExclusive.setVal(900*(rapMax-rapMin));
		//yieldJpsiExclusive.setRange(960, 1130);
		yieldJpsiDissociative.setRange(100*(rapMax-rapMin),1300*(rapMax-rapMin));
		yieldJpsiDissociative.setVal(500*(rapMax-rapMin));
		yieldJpsiInclusive.setRange(0,2000);
		yieldJpsiInclusive.setVal(600*(rapMax-rapMin));
	}
	else {
		yieldJpsiExclusive.setRange(300*(rapMax-rapMin), 1100*(rapMax-rapMin));
		yieldJpsiExclusive.setVal(700*(rapMax-rapMin));
		//yieldJpsiExclusive.setRange(580, 650);
		yieldJpsiDissociative.setRange(0,600*(rapMax-rapMin)); yieldJpsiDissociative.setVal(100*(rapMax-rapMin));
		if (!excOnly) {yieldJpsiInclusive.setVal(0); yieldJpsiInclusive.setConstant();}
	}
	 if (excOnly) {
         /*
         yieldBkg.setRange(0, 10);
         yieldBkg.setVal(0); yieldBkg.setConstant();
          */
         
		 yieldJpsiExclusive.setRange(80*(rapMax-rapMin),1200*(rapMax-rapMin));
		 yieldJpsiExclusive.setVal(800*(rapMax-rapMin));
		 //if (period == "LHC16s") {yieldJpsiDissociative.setVal(100);}
		 yieldJpsiDissociative.setRange(0,1300*(rapMax-rapMin));
		 yieldJpsiDissociative.setVal(0); yieldJpsiDissociative.setConstant();
	
		 yieldJpsiInclusive.setRange(0,1000);
		 yieldJpsiInclusive.setVal(0); yieldJpsiInclusive.setConstant();
      
		 /*
		 if (period == "LHC16s") {
			 yieldJpsiDissociative.setVal(0); yieldJpsiDissociative.setConstant();
		 }
		 else {
			 yieldJpsiDissociative.setRange(0,400*(rapMax-rapMin));
			 yieldJpsiDissociative.setVal(0);
		 }
		  */
	 }

	/*
	yieldJpsiExclusive.setRange(0, 10);
	yieldJpsiExclusive.setVal(0); yieldJpsiExclusive.setConstant();
	//yieldJpsiInclusive.setVal(0); yieldJpsiInclusive.setConstant();
	//yieldJpsiDissociative.setVal(0); yieldJpsiDissociative.setConstant();
	yieldJpsiGammaPb.setRange(0,yieldGammaPbValue);
	 */

	/*
	yieldJpsiDissociative.setRange(1300,1380);
	yieldJpsiDissociative.setVal(1300);
	 */

/*
	yieldJpsiExclusive.setRange(0, 10);
	yieldJpsiExclusive.setVal(0); yieldJpsiExclusive.setConstant();
 */
    /*
    yieldJpsiInclusive.setRange(0,1000);
    yieldJpsiInclusive.setVal(0); yieldJpsiInclusive.setConstant();
     */

    if (mMax < 3 || mMin > 3.2) {
        yieldJpsiExclusive.setRange(0,10);
        yieldJpsiExclusive.setVal(0); yieldJpsiExclusive.setConstant();
        yieldJpsiDissociative.setRange(0,10);
        yieldJpsiDissociative.setVal(0); yieldJpsiDissociative.setConstant();
        yieldJpsiInclusive.setRange(0,10);
        yieldJpsiInclusive.setVal(0); yieldJpsiInclusive.setConstant();
        yieldJpsiGammaPb.setRange(0,10);
        yieldJpsiGammaPb.setVal(0); yieldJpsiGammaPb.setConstant();
    }
    
    if (ptMin > 0.1) {
        yieldJpsiGammaPb.setRange(0,10);
        yieldJpsiGammaPb.setVal(0); yieldJpsiGammaPb.setConstant();
        yieldTwoGamma.setRange(0,10);
        yieldTwoGamma.setVal(0); yieldTwoGamma.setConstant();
    }
    
    //yieldTwoGamma.setVal(0); yieldTwoGamma.setConstant();
    yieldJpsiInclusive.setVal(0); yieldJpsiInclusive.setConstant();
	
	ws->import(yieldJpsiExclusive);
	ws->import(yieldJpsiDissociative);
	ws->import(yieldJpsiGammaPb);
	ws->import(yieldJpsiInclusive);
	ws->import(yieldPsi2s);
	ws->import(yieldTwoGamma);
	ws->import(yieldBkg);
}

