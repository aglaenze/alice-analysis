#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>

#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TDatime.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCut.h"
#include "TTree.h"

#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPlot.h"

#include "Include/FitUtils.C"
#include "Include/Initiate.C"

using namespace RooFit;
using namespace std;

//Double_t mLimitPsi2s = 3.65;
Double_t mLimitPsi2s = 10;


void AddModel(RooWorkspace* ws, string period, Double_t mMin, Double_t mMax) {
	// Define model
	
	//RooRealVar m("fTrkTrkM","M_{#mu#mu} (GeV/c2)",2.5,3.5);
	RooRealVar m = *ws->var("fTrkTrkM");
	
	LoadMassFitFunctions(ws, period);
	// Now collect m pdf
	RooAbsPdf* jpsi = ws->pdf("jpsi");
	RooAbsPdf* psi2s = ws->pdf("psi2s");
	RooAbsPdf* bkg = ws->pdf("bkg");
	
	//RooRealVar fsig("fsig","signalPhi",0.1,0.,1.);
	RooRealVar fsigJpsi("fsigJpsi","signalJPsi",1000,0.,1.e4);
	RooRealVar fsigPsi2s("fsigPsi2s","signalPsi2s",100,0.,1.e3);
	RooRealVar fbkg("fbkg","fbkg",500,0.,1.e7);
    
    if (mMin > 3.2 || mMax < 3) {fsigJpsi.setVal(0); fsigJpsi.setConstant();}
	
	RooAbsPdf* model;
	if (mMax > mLimitPsi2s) model = new RooAddPdf("mfit", "mfit", RooArgList(*jpsi, *psi2s, *bkg), RooArgList(fsigJpsi, fsigPsi2s, fbkg), kFALSE);
	else model = new RooAddPdf("mfit", "mfit", RooArgList(*jpsi, *bkg), RooArgList(fsigJpsi, fbkg), kFALSE);
	
	ws->import(*model);
    ws->import(*jpsi);
    ws->import(*jpsi);
}

void DoSPlot(RooWorkspace* ws, Double_t mMax) {
	
	RooAbsPdf* model = ws->pdf("mfit");
	RooRealVar* jPsiYield = ws->var("fsigJpsi");
	RooRealVar* psi2sYield = ws->var("fsigPsi2s");
	RooRealVar* bkgYield = ws->var("fbkg");
	RooDataSet* data = (RooDataSet*) ws->data("data");
	
	model->fitTo(*data, Extended(), Minos(true), Strategy(2));
	//The sPlot technique requires that we fix the parameters
	// of the model that are not yields after doing the fit
	// J/Psi
	RooRealVar* mean_jpsi = ws->var("mean_jpsi");
	RooRealVar* sigma_jpsi = ws->var("sigma_jpsi");
	//RooRealVar* alpha_jpsi = ws->var("alpha_jpsi");
	//RooRealVar* n_jpsi = ws->var("n_jpsi");
	RooRealVar* alpha_jpsi_L = ws->var("alpha_jpsi_L");
	RooRealVar* n_jpsi_L = ws->var("n_jpsi_L");
	RooRealVar* alpha_jpsi_R = ws->var("alpha_jpsi_R");
	RooRealVar* n_jpsi_R = ws->var("n_jpsi_R");
	
	// Psi(2s)
	RooRealVar *mean_psi2s, *sigma_psi2s, *alpha_psi2s_L, *n_psi2s_L, *alpha_psi2s_R, *n_psi2s_R;
	if (mMax > mLimitPsi2s) {
		mean_psi2s = ws->var("mean_psi2s");
		sigma_psi2s = ws->var("sigma_psi2s");
		alpha_psi2s_L = ws->var("alpha_psi2s_L");
		n_psi2s_L = ws->var("n_psi2s_L");
		alpha_psi2s_R = ws->var("alpha_psi2s_R");
		n_psi2s_R = ws->var("n_psi2s_R");
		alpha_psi2s_L->setConstant();
		n_psi2s_L->setConstant();
		alpha_psi2s_R->setConstant();
		n_psi2s_R->setConstant();
	}
	
	// Background
	RooRealVar* a1 = ws->var("a1");
	
	mean_jpsi->setConstant();
	sigma_jpsi->setConstant();
	alpha_jpsi_L->setConstant();
	n_jpsi_L->setConstant();
	alpha_jpsi_R->setConstant();
	n_jpsi_R->setConstant();
	
	a1->setConstant();
	
	RooMsgService::instance().setSilentMode(true);
	
	//Now we use the SPlot class to add SWeight to our data set
	// based on our model and our yield variables
	
	RooStats::SPlot * sData;
	if (mMax > mLimitPsi2s) sData = new RooStats::SPlot("sData","splot", *data, model, RooArgList(*jPsiYield, *psi2sYield, *bkgYield) );
	else sData = new RooStats::SPlot("sData","splot", *data, model, RooArgList(*jPsiYield, *bkgYield) );
	
	//Check Sweight properties
	std::cout << "Check SWeights: " << std::endl;
	
	std::cout << std::endl << "Yield of JPsi is "
	<< jPsiYield->getVal() << ". From sWeights it is "
	<< sData->GetYieldFromSWeight("fsigJpsi") << std::endl;
	
	
	std::cout << std::endl << "Yield of bkg is "
	<< bkgYield->getVal() << ". From sWeights it is "
	<< sData->GetYieldFromSWeight("fbkg") << std::endl;
	
	if (mMax > mLimitPsi2s) {
		std::cout << std::endl << "Yield of Psi(2s) is "
		<< psi2sYield->getVal() << ". From sWeights it is "
		<< sData->GetYieldFromSWeight("fsigPsi2s") << std::endl;
		for (Int_t i=0; i < 10; i ++ )
		{
			std::cout << "JPsi Weight "<< sData->GetSWeight(i, "fsigJpsi")
			<< " psi(2s) Yield "<< sData->GetSWeight(i, "fsigPsi2s")
			<< " bkg Yield "<< sData->GetSWeight(i, "fbkg")
			<< " Total Weight "<< sData->GetSumOfEventSWeight(i)
			<< std::endl;
		}
	}
	else {
		for (Int_t i=0; i < 10; i ++ )
		{
			std::cout << "JPsi Weight "<< sData->GetSWeight(i, "fsigJpsi")
			<< " bkg Yield "<< sData->GetSWeight(i, "fbkg")
			<< " Total Weight "<< sData->GetSumOfEventSWeight(i)
			<< std::endl;
		}
	}
	
	//import the new data set with Sweight
	
	std::cout << "import new dataset with sWeight" << std::endl;
	ws->import(*data, Rename("dataWithSWeights"));
	
	std::cout << "Sigma = " << sigma_jpsi->getVal() << "\n\n\n\n\n" << std::endl;
}


//____________________________________
void MakePlots(RooWorkspace *ws, string rootfilePath, string rootfilePathMc, string period, bool exp, bool exclusiveOnly, bool useCuts, Double_t mMin, Double_t mMax, Double_t ptMin, Double_t ptMax, Double_t rapMin, Double_t rapMax){
	
	int rapBinNumber = int(5*(rapMax-rapMin));
    int massBinNumber = 50;
	
	//make some plots
	//TCanvas* cv = new TCanvas("splot","splot", 900, 900) ;
	//cv->Divide(3,3);
    TCanvas* cv = new TCanvas("splot","splot", 900, 600) ;
    cv->Divide(3,2);
	//get what we need of the workspace
	RooAbsPdf* model = ws->pdf("mfit");
	RooAbsPdf* jPsiModel = ws->pdf("jpsi");
	RooAbsPdf* bkgModel = ws->pdf("bkg");
	
	RooRealVar* m = ws->var("fTrkTrkM");
	RooRealVar* rap = ws->var("fTrkTrkY");
	
	RooDataSet* sdata = (RooDataSet*) ws->data("dataWithSWeights");
	RooDataSet* data = (RooDataSet*) ws->data("data");
	
	RooRealVar* jPsiYield = ws->var("fsigJpsi");
	RooRealVar* psi2sYield = ws->var("fsigPsi2s");
	RooRealVar* bkgYield = ws->var("fbkg");
    //return;
	
	// Draw mass fit
	cv->cd(1);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	RooPlot* frame = m->frame();
	sdata->plotOn(frame, Binning(massBinNumber));
	model->plotOn(frame);
	model->plotOn(frame, Components(*model), LineStyle(kDashed), LineColor(kBlue));
	model->plotOn(frame, Components(*jPsiModel), LineStyle(kDashed), LineColor(kRed));
	model->plotOn(frame, Components(*bkgModel), LineStyle(kDashed), LineColor(kGreen));
	frame->SetTitle("Fit to model to discriminating variable");
	double yMax0 = frame->GetMaximum();
	double xPos = mMin+(mMax-mMin)*0.1;
	//if (mMax < mLimitPsi2s) xPos = mMin+(mMax-mMin)/3;
	TText* txt0 = new TText(xPos,0.7*yMax0,Form("%.1f J/Psi", jPsiYield->getVal()));
	TText* txt1 = new TText(xPos,0.55*yMax0,Form("%.1f bkg", bkgYield->getVal()));
	frame->addObject(txt0) ;
	frame->addObject(txt1) ;
	if (mMax > mLimitPsi2s) {
		RooAbsPdf* psi2sModel = ws->pdf("psi2s");
		model->plotOn(frame, Components(*psi2sModel), LineStyle(kDashed), LineColor(kOrange));
		TText* txtPsi2s = new TText((mMin+mMax+0.1)/2,0.45*yMax0,Form("%.1f Psi(2s)", psi2sYield->getVal()));
		frame->addObject(txtPsi2s) ;
	}
	frame->Draw();
    //return;
    
	// Draw quality plots
	cv->cd(4);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	sdata->plotOn(frame, Binning(massBinNumber));
    //sdata->plotOn(frame);
	model->plotOn(frame);
	RooHist* hpull = frame->pullHist();
	hpull->GetXaxis()->SetRangeUser(mMin, mMax);
	hpull->SetTitle("(data - fit)/#sigma of the mass fit");
	hpull->Draw("");
	cv->cd(5);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	RooHist* hresid = frame->residHist();
	hresid->GetXaxis()->SetRangeUser(mMin, mMax);
	hresid->SetTitle("Residuals (data - fit) of the mass fit");
	hresid->Draw("");
    //return;
	
	// The SPlot class adds a new variable that has the name of the corresponding
	// yield + "_sw".
	// create weighted data set for JPsi
	RooDataSet * dataw_jpsi = new RooDataSet(sdata->GetName(),sdata->GetTitle(),sdata,*sdata->get(),0,"fsigJpsi_sw") ;
	//RooDataSet * dataw_jpsi = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),"fTrkTrkM < 3.2 && fTrkTrkM > 2.8","fsigJpsi_sw") ;
	// create weighted data set for Background
	RooDataSet * dataw_bkg = nullptr;
	if (mMax > 3.1 && mMin < 3.1 ) dataw_bkg = new RooDataSet(sdata->GetName(),sdata->GetTitle(),sdata,*sdata->get(),0,"fbkg_sw") ;
	else dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0) ;
    //return;
	
	// Draw J/Psi pt with weights
	cv->cd(2);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	RooPlot* frameJpsi = rap->frame() ;
	int binningNum = 100;
	if (exclusiveOnly) {binningNum = 40;}
	dataw_jpsi->plotOn(frameJpsi, Name("jPsiData"), DataError(RooAbsData::SumW2), Binning(binningNum) ) ;
	frameJpsi->SetTitle("rapidity distribution for J/#Psi with weights");
	
	// Fit rap distrib of J/Psi with exclusive / dissociative components
	RooAbsPdf* rapModel = ws->pdf("rapfit");

    //return;
	// Write chi2
	double xPosJpsi = rapMax/3;
    /*
	int nDof = rapModel->getParameters(dataw_jpsi)->selectByAttrib("Constant",kFALSE)->getSize();
	double yMax = frameJpsi->GetMaximum();
	TLatex* txtChi = new TLatex(xPosJpsi/2, 0.93*yMax,Form("#chi^{2}/ndf = %.3f", frameJpsi->chiSquare( "rapfit", "jPsiData", nDof)));
	frameJpsi->addObject(txtChi);
    return;
     */

	frameJpsi->Draw() ;

    
	// Draw background rap with sweights
	cv->cd(3);
	gPad->SetLeftMargin(0.15) ;
	gPad->SetBottomMargin(0.15) ;
	RooPlot* frameBkg = rap->frame() ;
	frameBkg->SetTitle("rapidity distribution for background with weights");
	//frameBkg->GetXaxis()->SetRangeUser(0,3);
	
	TLegend* legend = new TLegend(0.4, 0.6, 0.9, 0.9);
	legend->SetTextSize(0.05);
	legend->SetFillColor(kWhite);
	legend->SetLineColor(kWhite);
	
	dataw_bkg->plotOn(frameBkg,DataError(RooAbsData::SumW2), Binning(binningNum)) ;
    frameBkg->Draw();
    
    // draw rapidity from MC
    cv->cd(6);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    RooRealVar rapMc("fTrkTrkY","Dimuon rap mC", rapMin, rapMax);
    RooRealVar mMc("fTrkTrkM","Dimuon m MC", mMin, mMax);
    RooRealVar ptMc("fTrkTrkPt","Dimuon pt MC", ptMin, ptMax);
    RooArgSet variablesMc(rapMc, mMc, ptMc);
    TString filenameMc = rootfilePathMc + Form("/AnalysisResults_%s_MC_kTwoGammaToMuLow.root", period.c_str());
    TFile* f = new TFile(filenameMc, "READ");
    TTree* tree = (TTree*)f->Get("MyTask/fAnaTree");
    TCut cut = "";
    RooDataSet* dataMc = new RooDataSet("data","data",variablesMc,Import(*tree),Cut(cut));
    RooPlot* frameRapMc = rapMc.frame();
    dataMc->plotOn(frameRapMc, Binning(binningNum));
    frameRapMc->SetTitle("rapidity distribution for Monte Carlo Gamma-Gamma");
    frameRapMc->Draw();
    
	
	
	// Save plots
	string cutType = "";
	if (!useCuts) cutType = "-nocuts";
	string suffix = "";
	if (exp) suffix = "exp";
	else suffix = "powerlaw";
	string suffix2 = "";
	if (exclusiveOnly) suffix2 = "-exclusive-only";
	cv->SaveAs(Form("Plots/%s/Splot%s-%.1f-%.1f-%s%s.pdf", period.c_str(), cutType.c_str(), mMin, mMax, suffix.c_str(), suffix2.c_str()));
}


//____________________________________
void DrawWeights(RooWorkspace *ws){
	
	
	//get what we need of the workspace
	RooRealVar* m = ws->var("fTrkTrkM");
	RooRealVar* rap = ws->var("fTrkTrkY");
	
	RooDataSet* sData = (RooDataSet*) ws->data("dataWithSWeights");
	
	//std::cout << std::endl << std::endl << std::endl << std::endl;
	const int nEntries = (int) sData->sumEntries();
	Double_t massForWeights[nEntries], rapForWeights[nEntries], jPsiWeight[nEntries], bkgWeight[nEntries];
	
	for (Int_t i=0; i < (Int_t)nEntries; i ++ ) {
		massForWeights[i] = sData->get(i)->getRealValue("fTrkTrkM");
		rapForWeights[i] = sData->get(i)->getRealValue("fTrkTrkY");
		
		jPsiWeight[i] = sData->get(i)->getRealValue("fsigJpsi_sw");
		bkgWeight[i] = sData->get(i)->getRealValue("fbkg_sw");
		
	}
	//std::cout << std::endl << std::endl << std::endl << std::endl;
	
	// TGraphs of sWeights = f(m) or sWeights = f(rap)
	TGraph* gWeightsMJpsi = new TGraph(nEntries,massForWeights, jPsiWeight);
	TGraph* gWeightsRapJpsi = new TGraph(nEntries,rapForWeights, jPsiWeight);
	TGraph* gWeightsMBkg = new TGraph(nEntries,massForWeights, bkgWeight);
	TGraph* gWeightsRapBkg = new TGraph(nEntries,rapForWeights, bkgWeight);
	
	gWeightsMJpsi->SetTitle("JPsi sWeights = f(m)");
	gWeightsRapJpsi->SetTitle("JPsi sWeights = f(rap)");
	gWeightsMBkg->SetTitle("Background sWeights = f(m)");
	gWeightsRapBkg->SetTitle("Background sWeights = f(rap)");
	
	gWeightsMJpsi->GetXaxis()->SetTitle("m");
	gWeightsRapJpsi->GetXaxis()->SetTitle("rap");
	gWeightsMBkg->GetXaxis()->SetTitle("m");
	gWeightsRapBkg->GetXaxis()->SetTitle("rap");
	
	gWeightsMJpsi->GetYaxis()->SetTitle("J/Psi sWeights");
	gWeightsRapJpsi->GetYaxis()->SetTitle("J/Psi sWeights");
	gWeightsMBkg->GetYaxis()->SetTitle("Bkg sWeights");
	gWeightsRapBkg->GetYaxis()->SetTitle("Bkg sWeights");
}

void SplotRap(string rootfilePath = "", string rootfilePathMc = "", std::vector<string> periods = {"LHC16r", "LHC16s"}) {
	
	Double_t mMin = 2., mMax = 3.5, ptMin = 0, ptMax = 5;
	Double_t rapMin = -4., rapMax = -2.5;
	bool useCuts = true, exp = false, exclusiveOnly = false;
	bool logScale = false, drawPulls = false;

	gROOT->ProcessLine(".L Include/ExtendedCrystalBall.cxx+") ;
	gSystem->Load("./Include/ExtendedCrystalBall_cxx.so") ;
	
	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetTitleSize(.05);
	gStyle->SetLabelSize(.05, "XY");
	//gStyle->SetMarkerSize(0.5);
	//gStyle->SetMarkerStyle(20);
	
	const int nPeriod = periods.size();
	for (int k = 0; k<nPeriod; k++) {
		string period = periods[k];
		if ( ! Initiate(period, mMin, mMax, ptMin, ptMax, rapMin, rapMax, useCuts, logScale, drawPulls, exp, exclusiveOnly)) {cout << "Something wrong at initialisation"; return;}
		

		if (exclusiveOnly) {
			ptMax = 0.8;
		}
         /*
		else {
			if (period == "LHC16s") ptMax = 3;
		}
		 */
		
		std::list<TCut> mCutList = DefineCuts(period, exclusiveOnly);
		// Define cuts
		TCut mCut = "";
		
		for (std::list <TCut>::iterator iter = mCutList.begin(); iter != mCutList.end(); ++iter) {mCut += *iter;}
		if (!useCuts) mCut = "";
		
		// Open the file
		TFile *fAna = new TFile(Form("%s/AnalysisResults_%s.root", rootfilePath.c_str(), period.c_str()), "READ");
		
		// Connect to the tree
		TTree* fAnaTree = (TTree*)fAna->Get("MyTask/fAnaTree");
		//Create a new workspace to manage the project
		RooWorkspace* wspace = new RooWorkspace("myJpsi");
		ImportDataSet(wspace, fAnaTree, mCut, mMin, mMax, ptMin, ptMax, rapMin, rapMax);
		AddModel(wspace, period, mMin, mMax);
		wspace->Print();
		
		DoSPlot(wspace, mMax);
		MakePlots(wspace, rootfilePath, rootfilePathMc, period, exp, exclusiveOnly, useCuts, mMin, mMax, ptMin, ptMax, rapMin, rapMax);
		//DrawWeights(wspace);
		
		//cleanup
		delete wspace;
	}
	
	
}


