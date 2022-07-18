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

using namespace RooFit;
using namespace std;

//vector <TString> processes = {"kCohJpsiToMu", "kIncohJpsiToMu", "kIncohPsi2sToMuPi", "kIncohPsi2sToMu", "kTwoGammaToMuLow", "kTwoGammaToMuMedium"};
//vector <TString> processes = {"kCohJpsiToMu", "kIncohJpsiToMu", "kIncohPsi2sToMu"};

//vector <TString> processes = {"kCohJpsiToMu", "kIncohJpsiToMu", "kIncohPsi2sToMuPi", "kIncohPsi2sToMu", "kTwoGammaToMuLow", "kTwoGammaToMuMedium"};
//vector <TString> processes = {"kCohJpsiToMu", "kIncohJpsiToMu", "kIncohPsi2sToMu"};

list<TCut> DefineCuts(string period, bool excOnly) {
    TCut ZNcut[2], ADcut[2], V0cut[2];
    
    int i = -1;
    if (period == "LHC16r") i = 0;
    else if (period == "LHC16s") i = 1;
    if (i == -1) {std::cout << "Wrong period in DefineCuts" << std::endl; return {};}
    
    ZNcut[0] = "!fZNAgoodTiming";    // no event in beam-beam time window [-2 ns, +2 ns] (because UPC only)
    ZNcut[1] = "!fZNCgoodTiming";
    //ZNcut = "fZNAfired==0";
    
    /*
     ZNcut[0] = "fZNAEnergy < 10 || !fZNAgoodTiming";    // no event in beam-beam time window (because UPC only)
     ZNcut[1] = "fZNCEnergy < 10 || !fZNCgoodTiming";
     */
    /*
     ZNcut[0] = "fZNAEnergy<30";    // no event in beam-beam time window (because UPC only)
     ZNcut[1] = "fZNCEnergy<30";
     */
    
    /*
     ZNcut[0] = "fZNATDC[0]>-2 && fZNATDC[0]<2";
     ZNcut[1] = "fZNCTDC[0]>-2 && fZNCTDC[0]<2";
     
     for (int j = 1; j<4; j++) {
     ZNcut[0] = ZNcut[0] || Form("fZNATDC[%d]>-2 && fZNATDC[%d]<2", j, j);
     ZNcut[1] = ZNcut[1] || Form("fZNCTDC[%d]>-2 && fZNCTDC[%d]<2", j, j);
     }
     ZNcut[0] = !ZNcut[0];
     ZNcut[1] = !ZNcut[1];
     */
    
    V0cut[0] = "fV0ADecision != 1";        // Pb does not dissociate
    /*
     // nMatch
     TCut v0nMatch0 = "fV0CCounts < 3 && !((fTrkEta1 > -3.7 && fTrkEta1 <-1.7) || (fTrkEta2 > -3.7 && fTrkEta2 <-1.7))";
     TCut v0nMatch1 = "fV0CCounts < 4 && ((fTrkEta1 > -3.7 && fTrkEta1 <-1.7) && !(fTrkEta2 > -3.7 && fTrkEta2 <-1.7))";
     TCut v0nMatch1bis = "fV0CCounts < 4 && (!(fTrkEta1 > -3.7 && fTrkEta1 <-1.7) && (fTrkEta2 > -3.7 && fTrkEta2 <-1.7))";
     TCut v0nMatch2 = "fV0CCounts < 5 && (fTrkEta1 > -3.7 && fTrkEta1 <-1.7 && fTrkEta2 > -3.7 && fTrkEta2 <-1.7)";
     V0cut[1] = v0nMatch0 || (v0nMatch1 || v0nMatch1bis) || v0nMatch2;
     */
    V0cut[1] = "fNV0C < fV0CNMatched + 3";
    //if (period == "LHC16s") V0cut[1] = "fNV0C < fV0CNMatched + 5";
    //V0cut[1] = "fNV0C < fV0CNMatched + 1";
    //V0cut[1] = "";
    
    ADcut[0] = "fADADecision != 1";    // Pb does not dissociate
    ADcut[1] = "fADCDecision != 1";
    
    //bool exclusiveOnly = false;    // if true, p is not allowed to break
    bool inclusiveOnly = false;
    bool dissociativeOnly = false;
    bool midSelection = true;
    if (excOnly) midSelection = true;
    
    TCut unlikeSignCut = "(fTrkQ1 < 0 && fTrkQ2 > 0) || (fTrkQ1 > 0 && fTrkQ2 < 0)";
    TCut sameSignCut = "(fTrkQ1 < 0 && fTrkQ2 < 0) || (fTrkQ1 > 0 && fTrkQ2 > 0)";
    list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ZNcut[i], ADcut[i], V0cut[i]};
    // test same sign muons
    //list<TCut> mCutList = {"", "fAnaType==0", sameSignCut, ZNcut[i], ADcut[i], V0cut[i]};
    //list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ADcut[i], V0cut[i], "fV0CBBNHits < 3"};
    
    if (midSelection) {
        TCut SPDcutExc = "fTracklets < 3";
        TCut V0cutExc = V0cut[0] && V0cut[1];
        //TCut V0cutExc = V0cut[i];
        mCutList.push_back(V0cutExc);
        mCutList.push_back(SPDcutExc);
    }
    
    if (excOnly) {
        TCut ZDCcutExc[2];
        ZDCcutExc[0]= "!fZNAgoodTiming && !fZNCgoodTiming2";
        ZDCcutExc[1]= "!fZNCgoodTiming && !fZNAgoodTiming2";
        //ZDCcutExc[1]= "";
        /*
         ZDCcutExc[0] = "fZNCEnergy < 5 || !fZNCgoodTiming";
         ZDCcutExc[1] = "fZNAEnergy < 5 || !fZNAgoodTiming";    // no event in beam-beam time window (because UPC only)
         */
        /*
         ZDCcutExc[0] = "fZNCEnergy < 5";
         ZDCcutExc[1] = "fZNAEnergy < 5";
         */
        //TCut ZDCcutExc = "(fZNAEnergy < 10 || !fZNAgoodTiming) && (fZNCEnergy < 10 || !fZNCgoodTiming)";
        //TCut V0cutExc = "fV0ADecision == 0 && (fV0CDecision == 0 || fV0ADecision == 1)";
        TCut ADcutExc = "fADADecision != 1 && fADCDecision != 1";    // Pb and p do not dissociate
        mCutList.push_back(ADcutExc);
        mCutList.push_back(ZDCcutExc[i]);
    }
    
    
    if (inclusiveOnly) {
        //list<TCut> cutListInclusive = {"fV0CBBNHits > 2"};
        //list<TCut> cutListInclusive = {"fTracklets > 0 || fADADecision == 1 || fADCDecision == 1 || fZNAgoodTiming || fZNCgoodTiming || fV0ADecision == 1 || fV0CDecision == 1"};
        list<TCut> cutListInclusive = {"fTracklets > 0 || fADADecision == 1 || fADCDecision == 1 || fV0ADecision == 1 || fV0CDecision == 1"};
        //list<TCut> cutListInclusive = {"fTracklets > 0"};
        for (list<TCut>::iterator it=cutListInclusive.begin(); it!=cutListInclusive.end(); it++)
        mCutList.push_back(*it);
    }
    if (dissociativeOnly) {
        ADcut[0] += "fADCDecision == 1";    // Pb does not dissociate
        ADcut[1] += "fADADecision == 1";
        TCut SPDcutExc = "fTracklets == 0";
        list<TCut> cutListDiss = {ADcut[i], SPDcutExc};
        for (list<TCut>::iterator it=cutListDiss.begin(); it!=cutListDiss.end(); it++)
        mCutList.push_back(*it);
    }
    
    //list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ZNcut[i], ADcut[i], V0cut[i]};
    //list<TCut> mCutList = {"", "fAnaType==0", ZNcut[i], ADcut[i], V0cut[i]};
    //list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ZNcut[i], ADcut[i]};
    //mCutList.push_back("fTrkTrkM>2.5 && fTrkTrkM<3.5");
    return mCutList;
}

list<TCut> DefineCutsBis(string period, bool excOnly) {
    TCut ZNcut[2], ADcut[2], V0cut[2];
    
    int i = -1;
    if (period == "LHC16r") i = 0;
    else if (period == "LHC16s") i = 1;
    if (i == -1) {std::cout << "Wrong period in DefineCuts" << std::endl; return {};}
    
    ZNcut[0] = "!fZNAgoodTiming";    // no event in beam-beam time window [-2 ns, +2 ns] (because UPC only)
    ZNcut[1] = "!fZNCgoodTiming";
    
    
    V0cut[0] = "fV0ADecision != 1";        // Pb does not dissociate
    V0cut[1] = "";
    
    ADcut[0] = "fADADecision != 1";    // Pb does not dissociate
    ADcut[1] = "fADCDecision != 1";
    
    TCut unlikeSignCut = "(fTrkQ1 < 0 && fTrkQ2 > 0) || (fTrkQ1 > 0 && fTrkQ2 < 0)";
    list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ZNcut[i], ADcut[i], V0cut[i]};
    
    TCut SPDcutExc = "fTracklets < 3";
    TCut V0cutExc = V0cut[0] && V0cut[1];
    //TCut V0cutExc = V0cut[i];
    mCutList.push_back(V0cutExc);
    mCutList.push_back(SPDcutExc);
    
    if (excOnly) {
        TCut ZDCcutExc[2];
        ZDCcutExc[0]= "!fZNAgoodTiming && !fZNCgoodTiming2";
        ZDCcutExc[1]= "!fZNCgoodTiming && !fZNAgoodTiming2";
        TCut ADcutExc = "fADADecision != 1 && fADCDecision != 1";    // Pb and p do not dissociate
        mCutList.push_back(ADcutExc);
        mCutList.push_back(ZDCcutExc[i]);
    }
    
    return mCutList;
}

list<TCut> DefineCutsTer(string period, bool excOnly) {
    TCut ZNcut[2], ADcut[2], V0cut[2];
    
    int i = -1;
    if (period == "LHC16r") i = 0;
    else if (period == "LHC16s") i = 1;
    if (i == -1) {std::cout << "Wrong period in DefineCuts" << std::endl; return {};}
    
    ZNcut[0] = "!fZNAgoodTiming";    // no event in beam-beam time window [-2 ns, +2 ns] (because UPC only)
    ZNcut[1] = "!fZNCgoodTiming";
    
    
    V0cut[0] = "fV0ADecision != 1";        // Pb does not dissociate
    V0cut[1] = "fNV0C < fV0CNMatched + 3";
    
    ADcut[0] = "fADADecision != 1";    // Pb does not dissociate
    ADcut[1] = "fADCDecision != 1";
    
    TCut unlikeSignCut = "(fTrkQ1 < 0 && fTrkQ2 > 0) || (fTrkQ1 > 0 && fTrkQ2 < 0)";
    list<TCut> mCutList = {"", "fAnaType==0", unlikeSignCut, ZNcut[i], ADcut[i], V0cut[i]};
    
    //TCut SPDcutExc = "fTracklets < 3";
    TCut SPDcutExc = "";
    TCut V0cutExc = V0cut[0] && V0cut[1];
    //TCut V0cutExc = V0cut[i];
    mCutList.push_back(V0cutExc);
    mCutList.push_back(SPDcutExc);
    
    if (excOnly) {
        TCut ZDCcutExc[2];
        ZDCcutExc[0]= "!fZNAgoodTiming && !fZNCgoodTiming2";
        ZDCcutExc[1]= "!fZNCgoodTiming && !fZNAgoodTiming2";
        TCut ADcutExc = "fADADecision != 1 && fADCDecision != 1";    // Pb and p do not dissociate
        mCutList.push_back(ADcutExc);
        mCutList.push_back(ZDCcutExc[i]);
    }
    
    return mCutList;
}


void ImportDataSet(RooWorkspace* ws, TTree* tree, TCut cut = "", Double_t mMin = 1.5, Double_t mMax = 5, Double_t ptMin = 0, Double_t ptMax = 8, Double_t rapMin = -4, Double_t rapMax = -2.5) {
    // import data
    
    RooRealVar m("fTrkTrkM","M_{#mu#mu} (GeV/c^{2})", mMin, mMax);
    RooRealVar pt("fTrkTrkPt","p_{T} (GeV/c)", ptMin, ptMax);
    RooRealVar phi("fTrkTrkPhi","Dimuon phi", -4, 4);
    RooRealVar fTrkTrkZ("fTrkTrkZ","Dimuon z", 0, 2);
    
    RooArgSet variables(m, pt, phi);
    
    //RooRealVar pt1("fTrkPt1","pt1",0,12);
    RooRealVar anaType("fAnaType","fAnaType",-2,4);
    RooRealVar fZNAgoodTiming("fZNAgoodTiming","fZNAgoodTiming",-1,2);
    RooRealVar fZNCgoodTiming("fZNCgoodTiming","fZNCgoodTiming",-1,2);
    RooRealVar fZNAgoodTiming2("fZNAgoodTiming2","fZNAgoodTiming2",-1,2);
    RooRealVar fZNCgoodTiming2("fZNCgoodTiming2","fZNCgoodTiming2",-1,2);
    
    RooRealVar fZNAEnergy("fZNAEnergy","fZNAEnergy",-1.e5,1.e5);
    RooRealVar fZNCEnergy("fZNCEnergy","fZNCEnergy",-1.e5,1.e5);
    RooRealVar fZPAEnergy("fZPAEnergy","fZPAEnergy",-1.e5,1.e5);
    RooRealVar fZPCEnergy("fZPCEnergy","fZPCEnergy",-1.e5,1.e5);
    
    // V0
    RooRealVar fV0ADecision("fV0ADecision","fV0ADecision",-1,5);
    RooRealVar fV0CDecision("fV0CDecision","fV0CDecision",-1,5);
    RooRealVar fV0ACounts("fV0ACounts","fV0ACounts",-1,1000);
    RooRealVar fV0CCounts("fV0CCounts","fV0CCounts",-1,1000);
    RooRealVar fV0CBBNHits("fV0CBBNHits","fV0CBBNHits", 0, 33);
    RooRealVar fV0CBGNHits("fV0CBGNHits","fV0CBGNHits", 0, 33);
    RooRealVar fNMatched("fV0CNMatched","fV0CNMatched", 0, 3);
    //RooRealVar fNMatched("fNMatched","fNMatched", 0, 3);
    RooRealVar fNV0C("fNV0C","fNV0C", -1, 1000);
    
    RooRealVar fADADecision("fADADecision","fADADecision",-1,5);
    RooRealVar fADCDecision("fADCDecision","fADCDecision",-1,5);
    
    RooRealVar fADABBNHits("fADABBNHits","fADABBNHits",-1,30);
    RooRealVar fADCBBNHits("fADCBBNHits","fADCBBNHits",-1,30);
    
    //RooRealVar fTrkTrkY("fTrkTrkY","fTrkTrkY",-20,20);
    RooRealVar fTrkTrkY("fTrkTrkY","fTrkTrkY", rapMin, rapMax);
    //RooRealVar fTrkTrkY("fTrkTrkY","fTrkTrkY",-3.25,-2.5);
    //RooRealVar fTrkTrkY("fTrkTrkY","fTrkTrkY",-3.2,-2.7);
    
    RooRealVar fTrkEta1("fTrkEta1","fTrkEta1", -5, 5);
    RooRealVar fTrkEta2("fTrkEta2","fTrkEta2", -5, 5);
    
    RooRealVar fTracklets("fTracklets","fTracklets",-2.0,3.e5);
    
    RooRealVar fTrkQ1("fTrkQ1","fTrkQ1",-5,5);
    RooRealVar fTrkQ2("fTrkQ2","fTrkQ2",-5,5);
    
    variables.add(anaType);
    variables.add(fZNAgoodTiming);
    variables.add(fZNCgoodTiming);
    variables.add(fZNAgoodTiming2);
    variables.add(fZNCgoodTiming2);
    
    variables.add(fZNAEnergy);
    variables.add(fZNCEnergy);
    variables.add(fZPAEnergy);
    variables.add(fZPCEnergy);
    
    variables.add(fV0ADecision);
    variables.add(fV0CDecision);
    
    variables.add(fV0ACounts);
    variables.add(fV0CCounts);
    variables.add(fV0CBBNHits);
    variables.add(fV0CBGNHits);
    variables.add(fNMatched);
    variables.add(fNV0C);
    
    variables.add(fADADecision);
    variables.add(fADCDecision);
    
    variables.add(fADABBNHits);
    variables.add(fADCBBNHits);
    
    variables.add(fTrkTrkY);
    variables.add(fTrkEta1);
    variables.add(fTrkEta2);
    
    variables.add(fTrkQ1);
    variables.add(fTrkQ2);
    variables.add(fTracklets);
    
    variables.add(fTrkTrkZ);
    
    RooDataSet* data = new RooDataSet("data","data",variables,Import(*tree),Cut(cut));
    
    data->Print();
    //return data;
    ws->import(*data, Rename("data"));
    
    Int_t nData = data->numEntries();
    Int_t nEntries = tree->GetEntries();
    std::cout << "\n\nNumber of entries = " << nData << "\nAnd entries in TTree = "<< nEntries << "\n\n" << std::endl;
}

void ImportDataSet(RooWorkspace* ws, TTree* tree, Double_t mMin = 1.5, Double_t mMax = 5, Double_t ptMin = 0, Double_t ptMax = 8) {
    ImportDataSet(ws, tree, "fTrkTrkM>0", mMin, mMax, ptMin, ptMax);
    
}
