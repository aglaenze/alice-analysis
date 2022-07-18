/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// c++ headers
#include <iostream>
#include <fstream>

// root headers
#include <TMath.h>
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <TStopwatch.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TLatex.h>
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TList.h"
#include "TChain.h"


// aliroot headers
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMuonTrackCuts.h"
#include "AliPIDResponse.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"


// my headers
#include "AliAnalysisTaskNanoJPsi2016.h"



class AliAnalysisTaskNanoJPsi2016;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskNanoJPsi2016) // classimp: necessary for root


AliAnalysisTaskNanoJPsi2016::AliAnalysisTaskNanoJPsi2016() : AliAnalysisTaskSE(),
fPeriod(0),
fIsMC(0),
fPIDResponse(0),
fAOD(0),
fOutputList(0),
fCounterH(0),
fMuonTrackCounterH(0),
fTriggerCounterFwdH(0),
fMuonTrackCuts(0x0),
fAnaTree(0),
fRunNum(0),
fOrbitNum(0),
fBunchCrossNum(0),
fL0inputs(0),
fTracklets(0),
fAnaType(-1),
fZNAfired(0),
fZNCfired(),
fZNCEnergy(0),
fZNAEnergy(0),
fZPCEnergy(0),
fZPAEnergy(0),
fV0ADecision(-10),
fV0CDecision(-10),
fV0ACounts(0),
fV0CCounts(0),
fV0AMultiplicity(0),
fV0CMultiplicity(0),
//fV0ARingMultiplicity(0),
//fV0CRingMultiplicity(0),
fV0ATime(0),
fV0CTime(0),
fV0ABBNHits(0),
fV0CBBNHits(0),
fV0ABGNHits(0),
fV0CBGNHits(0),
fV0CNMatched(0),
//fV0BB(0),
//fV0BG(0),
fADADecision(-10),
fADCDecision(-10),
fADABBNHits(0),
fADCBBNHits(0),
fADABGNHits(0),
fADCBGNHits(0),
fIR1Map(0),
fIR2Map(0),
fTrkTrkPt(0),
fTrkTrkPhi(0),
fTrkTrkY(0),
fTrkTrkM(0),
fTrkTrkEnergy(0),
fTrkTrkZ(0),
fTrkTrkZBis(0),
fTrkTrkW(0),
fTrkPt1(0),
fTrkPt2(0),
fTrkEta1(0),
fTrkEta2(0),
fTrkPhi1(0),
fTrkPhi2(0),
fTrkQ1(0),
fTrkQ2(0),
fTrkRabs1(0),
fTrkRabs2(0),
fGenTree(0),
fMCTrkTrkPt(0),
fMCTrkTrkPhi(0),
fMCTrkTrkY(0),
fMCTrkTrkM(0),
fMCTrkTrkEnergy(0),
fMCTrkTrkZ(0),
fMCTrkTrkW(0),
fMCTrkPt1(0),
fMCTrkPt2(0),
fMCTrkEta1(0),
fMCTrkEta2(0),
fMCTrkPhi1(0),
fMCTrkPhi2(0),
fMCTrkQ1(0),
fMCTrkQ2(0),
fEventLabel(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskNanoJPsi2016::AliAnalysisTaskNanoJPsi2016(const char* name) : AliAnalysisTaskSE(name),
fPeriod(0),
fIsMC(0),
fPIDResponse(0),
fAOD(0),
fOutputList(0),
fCounterH(0),
fMuonTrackCounterH(0),
fTriggerCounterFwdH(0),
fMuonTrackCuts(0x0),
fAnaTree(0),
fRunNum(0),
fOrbitNum(0),
fBunchCrossNum(0),
fL0inputs(0),
fTracklets(0),
fAnaType(-1),
fZNAfired(0),
fZNCfired(),
fZNCEnergy(0),
fZNAEnergy(0),
fZPCEnergy(0),
fZPAEnergy(0),
fV0ADecision(-10),
fV0CDecision(-10),
fV0ACounts(0),
fV0CCounts(0),
fV0AMultiplicity(0),
fV0CMultiplicity(0),
//fV0ARingMultiplicity(0),
//fV0CRingMultiplicity(0),
fV0ATime(0),
fV0CTime(0),
fV0ABBNHits(0),
fV0CBBNHits(0),
fV0ABGNHits(0),
fV0CBGNHits(0),
fV0CNMatched(0),
//fV0BB(0),
//fV0BG(0),
fADADecision(-10),
fADCDecision(-10),
fADABBNHits(0),
fADCBBNHits(0),
fADABGNHits(0),
fADCBGNHits(0),
fIR1Map(0),
fIR2Map(0),
fTrkTrkPt(0),
fTrkTrkPhi(0),
fTrkTrkY(0),
fTrkTrkM(0),
fTrkTrkEnergy(0),
fTrkTrkZ(0),
fTrkTrkZBis(0),
fTrkTrkW(0),
fTrkPt1(0),
fTrkPt2(0),
fTrkEta1(0),
fTrkEta2(0),
fTrkPhi1(0),
fTrkPhi2(0),
fTrkQ1(0),
fTrkQ2(0),
fTrkRabs1(0),
fTrkRabs2(0),
fGenTree(0),
fMCTrkTrkPt(0),
fMCTrkTrkPhi(0),
fMCTrkTrkY(0),
fMCTrkTrkM(0),
fMCTrkTrkEnergy(0),
fMCTrkTrkZ(0),
fMCTrkTrkW(0),
fMCTrkPt1(0),
fMCTrkPt2(0),
fMCTrkEta1(0),
fMCTrkEta2(0),
fMCTrkPhi1(0),
fMCTrkPhi2(0),
fMCTrkQ1(0),
fMCTrkQ2(0),
fEventLabel(0)
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
    DefineOutput(3, TTree::Class());
    
}
//_____________________________________________________________________________
AliAnalysisTaskNanoJPsi2016::~AliAnalysisTaskNanoJPsi2016()
{
    // destructor
    // liberate all allocated memory
    if(fOutputList) {delete fOutputList;}
    if(fMuonTrackCuts) {delete fMuonTrackCuts;}
    if(fAnaTree) {delete fAnaTree;}
    if(fCounterH) {delete fCounterH;}
    if(fMuonTrackCounterH) {delete fMuonTrackCounterH;}
    if(fTriggerCounterFwdH) {delete fTriggerCounterFwdH;}
    if(fPIDResponse) {delete fPIDResponse;}
    if(fGenTree) {delete fGenTree;}
    
}
//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    
    ////////////////////////////////////////
    //muon track cuts
    ////////////////////////////////////////
    fMuonTrackCuts = new AliMuonTrackCuts("StdMuonCuts", "StdMuonCuts");
    fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca | AliMuonTrackCuts::kMuMatchLpt);
    //fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca);
    fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
    fMuonTrackCuts->Print("mask");
    
    
    ////////////////////////////////////////
    //pid
    ////////////////////////////////////////
    
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    
    ////////////////////////////////////////
    //output histograms
    ////////////////////////////////////////
    
    fOutputList = new TList();          // this is a list which will contain all  histograms
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
    //  counter for events passing each cut
    fCounterH = new TH1F("fCounterH", "fCounterH", 25, -0.5, 24.5);
    fOutputList->Add(fCounterH);
    //  counter for tracks passing each cut
    fMuonTrackCounterH = new TH1F("fMuonTrackCounterH", "fMuonTrackCounterH", 10, -0.5, 9.5);
    fOutputList->Add(fMuonTrackCounterH);
    //  counter for triggers per run
    Int_t FirstRun = 1;
    Int_t LastRun = -1;
    if (fPeriod == 0 ) {FirstRun=265589;LastRun=266318;}// 16r
    else if (fPeriod == 1) {FirstRun=266405;LastRun=267131;}// 16s
    Int_t nRuns = LastRun-FirstRun+1;
    fTriggerCounterFwdH = new TH1F("fTriggerCounterFwdH", "fTriggerCounterFwdH", nRuns, FirstRun-0.5, LastRun+0.5);
    fOutputList->Add(fTriggerCounterFwdH);
    
    // post data
    PostData(1, fOutputList);
    
    ////////////////////////////////////////
    //output tree
    ////////////////////////////////////////
    fAnaTree = new TTree("fAnaTree", "fAnaTree");
    fAnaTree->Branch("fEventLabel", &fEventLabel, "fEventLabel/I");
    fAnaTree->Branch("fRunNum", &fRunNum, "fRunNum/I");
    fAnaTree->Branch("fOrbitNum", &fOrbitNum, "fOrbitNum/I");
    fAnaTree->Branch("fBunchCrossNum", &fBunchCrossNum, "fBunchCrossNum/I");
    fAnaTree->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
    fAnaTree->Branch("fAnaType", &fAnaType, "fAnaType/I");
    if (!fIsMC) {
        fAnaTree->Branch("fTracklets", &fTracklets, "fTracklets/I");
        fAnaTree->Branch("fZNAfired", &fZNAfired, "fZNAfired/O");
        fAnaTree->Branch("fZNCfired", &fZNCfired, "fZNCfired/O");
        fAnaTree->Branch("fZNATDC", &fZNATDC[0], "fZNATDC[4]/D");
        fAnaTree->Branch("fZNCTDC", &fZNCTDC[0], "fZNCTDC[4]/D");
        fAnaTree->Branch("fZNAgoodTiming", &fZNAgoodTiming, "fZNAgoodTiming/O");
        fAnaTree->Branch("fZNCgoodTiming", &fZNCgoodTiming, "fZNCgoodTiming/O");
        fAnaTree->Branch("fZNAgoodTiming2", &fZNAgoodTiming2, "fZNAgoodTiming2/O");
        fAnaTree->Branch("fZNCgoodTiming2", &fZNCgoodTiming2, "fZNCgoodTiming2/O");
        fAnaTree->Branch("fZNCEnergy", &fZNCEnergy, "fZNCEnergy/D");
        fAnaTree->Branch("fZNAEnergy", &fZNAEnergy, "fZNAEnergy/D");
        fAnaTree->Branch("fZPCEnergy", &fZPCEnergy, "fZPCEnergy/D");
        fAnaTree->Branch("fZPAEnergy", &fZPAEnergy, "fZPAEnergy/D");
        fAnaTree->Branch("fZPATDC", &fZPATDC[0], "fZPATDC[4]/D");
        fAnaTree->Branch("fZPCTDC", &fZPCTDC[0], "fZPCTDC[4]/D");
        fAnaTree->Branch("fZPAgoodTiming", &fZPAgoodTiming, "fZPAgoodTiming/O");
        fAnaTree->Branch("fZPCgoodTiming", &fZPCgoodTiming, "fZPCgoodTiming/O");
        fAnaTree->Branch("fZPAgoodTiming2", &fZPAgoodTiming2, "fZPAgoodTiming2/O");
        fAnaTree->Branch("fZPCgoodTiming2", &fZPCgoodTiming2, "fZPCgoodTiming2/O");
        fAnaTree->Branch("fV0ADecision", &fV0ADecision, "fV0ADecision/I");
        fAnaTree->Branch("fV0CDecision", &fV0CDecision, "fV0CDecision/I");
        fAnaTree->Branch("fV0ACounts", &fV0ACounts, "fV0ACounts/I");
        fAnaTree->Branch("fV0CCounts", &fV0CCounts, "fV0CCounts/I");
        fAnaTree->Branch("fNV0C", &fNV0C, "fNV0C/I");
        // Change here
        fAnaTree->Branch("fV0AMultiplicity", &fV0AMultiplicity, "fV0AMultiplicity/D");
        fAnaTree->Branch("fV0CMultiplicity", &fV0CMultiplicity, "fV0CMultiplicity/D");
        fAnaTree->Branch("fV0ARingMultiplicity", &fV0ARingMultiplicity[0], "fV0ARingMultiplicity[4]/D");
        fAnaTree->Branch("fV0CRingMultiplicity", &fV0CRingMultiplicity[0], "fV0CRingMultiplicity[4]/D");
        fAnaTree->Branch("fV0ATime", &fV0ATime, "fV0ATime/D");
        fAnaTree->Branch("fV0CTime", &fV0CTime, "fV0CTime/D");
        fAnaTree->Branch("fV0ABBNHits", &fV0ABBNHits, "fV0ABBNHits/I");
        fAnaTree->Branch("fV0CBBNHits", &fV0CBBNHits, "fV0CBBNHits/I");
        fAnaTree->Branch("fV0ABGNHits", &fV0ABGNHits, "fV0ABGNHits/I");
        fAnaTree->Branch("fV0CBGNHits", &fV0CBGNHits, "fV0CBGNHits/I");
        //fAnaTree->Branch("fV0BB", &fV0BB[0], "fV0BB[64]/O");
        fAnaTree->Branch("fV0AOnlineTrigger", &fV0AOnlineTrigger[0], "fV0AOnlineTrigger[32]/O");
        fAnaTree->Branch("fV0COnlineTrigger", &fV0COnlineTrigger[0], "fV0COnlineTrigger[32]/O");
        fAnaTree->Branch("fV0AOfflineTrigger", &fV0AOfflineTrigger[0], "fV0AOfflineTrigger[32]/O");
        fAnaTree->Branch("fV0COfflineTrigger", &fV0COfflineTrigger[0], "fV0COfflineTrigger[32]/O");
        
        fAnaTree->Branch("fV0BG", &fV0BG[0], "fV0BG[64]/O");
        fAnaTree->Branch("fV0CNMatched", &fV0CNMatched, "fV0CNMatched/I");
        fAnaTree->Branch("fV0AvgPhi", &fV0AvgPhi[0], "fV0AvgPhi[64]/D");
        fAnaTree->Branch("fV0EtaMin", &fV0EtaMin[0], "fV0EtaMin[64]/D");
        fAnaTree->Branch("fV0EtaMax", &fV0EtaMin[0], "fV0EtaMax[64]/D");
        // end of new V0 features
        fAnaTree->Branch("fADADecision", &fADADecision, "fADADecision/I");
        fAnaTree->Branch("fADCDecision", &fADCDecision, "fADCDecision/I");
        fAnaTree->Branch("fADABBNHits", &fADABBNHits, "fADABBNHits/I");
        fAnaTree->Branch("fADCBBNHits", &fADCBBNHits, "fADCBBNHits/I");
        fAnaTree->Branch("fADABGNHits", &fADABGNHits, "fADABGNHits/I");
        fAnaTree->Branch("fADCBGNHits", &fADCBGNHits, "fADCBGNHits/I");
        fAnaTree->Branch("fADBB", &fADBB[0], "fADBB[16]/O");
        fAnaTree->Branch("fADBG", &fADBG[0], "fADBG[16]/O");
        fAnaTree->Branch("fIR1Map", &fIR1Map);
        fAnaTree->Branch("fIR2Map", &fIR2Map);
    }
    fAnaTree->Branch("fTrkTrkPt", &fTrkTrkPt, "fTrkTrkPt/D");
    fAnaTree->Branch("fTrkTrkPhi", &fTrkTrkPhi, "fTrkTrkPhi/D");
    fAnaTree->Branch("fTrkTrkY", &fTrkTrkY, "fTrkTrkY/D");
    fAnaTree->Branch("fTrkTrkM", &fTrkTrkM, "fTrkTrkM/D");
    fAnaTree->Branch("fTrkTrkEnergy", &fTrkTrkEnergy, "fTrkTrkEnergy/D");
    fAnaTree->Branch("fTrkTrkZ", &fTrkTrkZ, "fTrkTrkZ/D");
    fAnaTree->Branch("fTrkTrkZBis", &fTrkTrkZBis, "fTrkTrkZBis/D");
    fAnaTree->Branch("fTrkTrkW", &fTrkTrkW, "fTrkTrkW/D");
    fAnaTree->Branch("fTrkPt1", &fTrkPt1, "fTrkPt1/D");
    fAnaTree->Branch("fTrkPt2", &fTrkPt2, "fTrkPt2/D");
    fAnaTree->Branch("fTrkEta1", &fTrkEta1, "fTrkEta1/D");
    fAnaTree->Branch("fTrkEta2", &fTrkEta2, "fTrkEta2/D");
    fAnaTree->Branch("fTrkPhi1", &fTrkPhi1, "fTrkPhi1/D");
    fAnaTree->Branch("fTrkPhi2", &fTrkPhi2, "fTrkPhi2/D");
    fAnaTree->Branch("fTrkQ1", &fTrkQ1, "fTrkQ1/D");
    fAnaTree->Branch("fTrkQ2", &fTrkQ2, "fTrkQ2/D");
    fAnaTree->Branch("fTrkRabs1", &fTrkRabs1, "fTrkRabs1/D");
    fAnaTree->Branch("fTrkRabs2", &fTrkRabs2, "fTrkRabs2/D");
    
    // post data
    PostData(2, fAnaTree);
    
    ////////////////////////////////////////
    //output tree
    ////////////////////////////////////////
    fGenTree = new TTree("fGenTree", "fGenTree");
    fGenTree->Branch("fEventLabel", &fEventLabel, "fEventLabel/I");
    fGenTree->Branch("fRunNum", &fRunNum, "fRunNum/I");
    fGenTree->Branch("fMCTrkTrkPt", &fMCTrkTrkPt, "fMCTrkTrkPt/D");
    fGenTree->Branch("fMCTrkTrkPhi", &fMCTrkTrkPhi, "fMCTrkTrkPhi/D");
    fGenTree->Branch("fMCTrkTrkY", &fMCTrkTrkY, "fMCTrkTrkY/D");
    fGenTree->Branch("fMCTrkTrkM", &fMCTrkTrkM, "fMCTrkTrkM/D");
    fGenTree->Branch("fMCTrkTrkEnergy", &fMCTrkTrkEnergy, "fMCTrkTrkEnergy/D");
    fGenTree->Branch("fMCTrkTrkZ", &fMCTrkTrkZ, "fMCTrkTrkZ/D");
    fGenTree->Branch("fMCTrkTrkW", &fMCTrkTrkW, "fMCTrkTrkW/D");
    fGenTree->Branch("fMCTrkPt1", &fMCTrkPt1, "fMCTrkPt1/D");
    fGenTree->Branch("fMCTrkPt2", &fMCTrkPt2, "fMCTrkPt2/D");
    fGenTree->Branch("fMCTrkEta1", &fMCTrkEta1, "fMCTrkEta1/D");
    fGenTree->Branch("fMCTrkEta2", &fMCTrkEta2, "fMCTrkEta2/D");
    fGenTree->Branch("fMCTrkPhi1", &fMCTrkPhi1, "fMCTrkPhi1/D");
    fGenTree->Branch("fMCTrkPhi2", &fMCTrkPhi2, "fMCTrkPhi2/D");
    fGenTree->Branch("fMCTrkQ1", &fMCTrkQ1, "fMCTrkQ1/I");
    fGenTree->Branch("fMCTrkQ2", &fMCTrkQ2, "fMCTrkQ2/I");
    // post data
    if (fIsMC) PostData(3, fGenTree);
    
}
//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016::NotifyRun()
{
    /// Set run number for cuts
    fMuonTrackCuts->SetRun(fInputHandler);
}

//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016::SetPeriod(Int_t period)
{
    // period = 0 => 2016 s, = 1 => 2016 r
    fPeriod = period;
}

//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016::SetMC(Bool_t flag)
{
    // set if MC file
    fIsMC = flag;
}
//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016::TrkTrkKine(Int_t idxTrk1, Int_t idxTrk2, Double_t TrkMass)
{
    // --  positive track
    TLorentzVector LV1;
    //AliAODTrack *PosTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(idxPosTrk[0]));
    AliAODTrack *Track1 = static_cast<AliAODTrack*>(fAOD->GetTrack(idxTrk1));
    LV1.SetPtEtaPhiM(Track1->Pt(), Track1->Eta(), Track1->Phi(), TrkMass);
    // --  negative track
    TLorentzVector LV2;
    //AliAODTrack *NegTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(idxNegTrk[0]));
    AliAODTrack *Track2 = static_cast<AliAODTrack*>(fAOD->GetTrack(idxTrk2));
    LV2.SetPtEtaPhiM(Track2->Pt(), Track2->Eta(), Track2->Phi(), TrkMass);
    
    // vector of Trk+Trk
    TrkTrk = LV2+LV1;
    
    // set tree variables
    fTrkTrkPt = TrkTrk.Pt();
    fTrkTrkPhi = TrkTrk.Phi();
    fTrkTrkY = TrkTrk.Rapidity();
    fTrkTrkM = TrkTrk.M();
    fTrkTrkEnergy = TrkTrk.Energy();
    fTrkPt1 = Track1->Pt();
    fTrkPt2 = Track2->Pt();
    fTrkEta1 = Track1->Eta();
    fTrkEta2 = Track2->Eta();
    fTrkPhi1 = Track1->Phi();
    fTrkPhi2 = Track2->Phi();
    fTrkQ1 = Track1->Charge();
    fTrkQ2 = Track2->Charge();
    fTrkRabs1 = Track1->GetRAtAbsorberEnd();
    fTrkRabs2 = Track2->GetRAtAbsorberEnd();
    
    if (fPeriod == 0) {
        TrkTrk.RotateZ(180);
    }
    
}

//_____________________________________________________________________________

void AliAnalysisTaskNanoJPsi2016::CheckTrigger(Bool_t &isTriggered)
// checks if event is triggered according to period and analysis type
{
    // Initialise: 0 = fwd, 1 = cent, 2 = semi-fwd
    
    // read trigger info
    TString trigger = fAOD->GetFiredTriggerClasses();
    
    // forward analysis
    // in 2016 r : CMUP14-B-NOPF-MUFAST = 0MSL *0VBA *0UBA
    // in 2016 s : CMUP23-B-NOPF-MUFAST = 0MUL *0UBC *0UGC *0VBA *0VGA *0SH2 *0VC5
    // central analysis
    // in 2016 r CCUP20-B-NOPF-CENTNOTRD = !VBA & !VGA & !VBC & !UBA & !UGC & !SH2 & STG & OMU
    // in 2016 s CCUP22-B-SPD2-CENTNOTRD = !VBA & !VGA & !VBC & !UBC & !UGC & !SH2 & STG & OMU
    // semi-forward analysis
    // in 2016 r CMUP15-B-NOPF-ALLNOTRD = *0VBA *0UBA *0VC5 0SMB *0SH2 0MSL
    // in 2016 s CMUP22-B-NOPF-ALLNOTRD = *0UBC *0UGC *0VBA *0VGA *0SH2 *0VC5 0MSL 0SMB
    if (fPeriod == 0 && trigger.Contains("CMUP14-B-NOPF-MUFAST")) {isTriggered = kTRUE;
    }
    else if (fPeriod == 1 && trigger.Contains("CMUP23-B-NOPF-MUFAST")) {isTriggered = kTRUE;
    }
    
    //if (trigger.Contains("CTRUE-B")) {isTriggered = kTRUE;}
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskNanoJPsi2016::TagCellVZEROC(Double_t ItsEta, Double_t ItsPhi)
{
    // coverage in eta of v0c
    // -3.7,-3.2,-2.7,-2.2,-1.7 for ring 0, 1, 2, 3 and 4, respectively
    
    Int_t i_eta = -1;
    Int_t idx   = -1;
    if      ( ItsEta>-3.7 && ItsEta<-3.2) i_eta = 0;
    else if ( ItsEta>-3.2 && ItsEta<-2.7) i_eta = 1;
    else if ( ItsEta>-2.7 && ItsEta<-2.2) i_eta = 2;
    else if ( ItsEta>-2.2 && ItsEta<-1.7) i_eta = 3;
    
    if (i_eta> -1) {
        Int_t i_phi = (Int_t) ((4.0*ItsPhi/TMath::Pi()));
        idx = i_eta*8+i_phi;
    }
    return idx;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskNanoJPsi2016::GoodCentralTrack(Int_t iTrack)
// selects good  tracks (w/o PID)
// bit 5 definition according to Michal:
//     70 < TPC clusters
//     4 > chi2 per TPC cluster
//     36 chi2 per ITS cluster
//     TPC and ITS refit
//     At least one SPD cluster
//     Tight pt dependent DCA cuts
{
    fCentralTrackCounterH->Fill(0);
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track) return kFALSE;
    fCentralTrackCounterH->Fill(1);
    if(!track->TestFilterBit(1<<5)) return kFALSE;
    fCentralTrackCounterH->Fill(2);
    if((!track->HasPointOnITSLayer(0))&&(!track->HasPointOnITSLayer(1))) return kFALSE;
    fCentralTrackCounterH->Fill(3);
    return kTRUE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskNanoJPsi2016::GoodMUONTrack(Int_t iTrack)
// selects good MUON tracks
{
    fMuonTrackCounterH->Fill(0);
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track) return kFALSE;
    fMuonTrackCounterH->Fill(1);
    if(!track->IsMuonTrack())  return kFALSE;
    fMuonTrackCounterH->Fill(2);
    if(!fMuonTrackCuts->IsSelected(track))  return kFALSE;
    fMuonTrackCounterH->Fill(3);
    if(track->GetRAtAbsorberEnd()>89.5)  return kFALSE;
    fMuonTrackCounterH->Fill(4);
    if(track->GetRAtAbsorberEnd()<17.5)  return kFALSE;
    fMuonTrackCounterH->Fill(5);
    return kTRUE;
}

Double_t Square(Double_t x) {return x*x;}

//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016::UserExec(Option_t *)
{
    ////////////////////////////////////////////
    // general info of the event
    ////////////////////////////////////////////
    Int_t iSelectionCounter = 0; // no selection applied yet
    fCounterH->Fill(iSelectionCounter); // entering UserExec (0)
    iSelectionCounter++;
    
    fEventLabel++;
    
    // get AOD event
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) {
        PostData(1, fOutputList);
        PostData(2, fAnaTree);
        return;
    }
    fCounterH->Fill(iSelectionCounter); // AOD event found (1)
    iSelectionCounter++;
    
    // get the run number
    fRunNum = fAOD->GetRunNumber();
    fOrbitNum = fAOD->GetOrbitNumber();
    fBunchCrossNum = fAOD->GetBunchCrossNumber();
    
    // is the right trigger?
    Bool_t isTriggered = kFALSE;
    // Monte Carlo
    if (fIsMC) {
        isTriggered = kTRUE;
        fMCEvent = MCEvent();
        if(fMCEvent) ProcessMCParticles();
    }
    else {
        CheckTrigger(isTriggered);
    }
    
    if (!isTriggered) {
        PostData(1, fOutputList);
        PostData(2, fAnaTree);
        return;
    }
    fCounterH->Fill(iSelectionCounter); // right trigger found (2)
    iSelectionCounter++;
    
    // trigger inputs
    fL0inputs = fAOD->GetHeader()->GetL0TriggerInputs();
    if (isTriggered) fTriggerCounterFwdH->Fill(fRunNum);
    
    ////////////////////////////////////////////
    //  find events with two good tracks
    ////////////////////////////////////////////
    
    //are there tracks at all?
    Int_t nTracks(fAOD->GetNumberOfTracks());
    if(nTracks<1) {
        PostData(1, fOutputList);
        PostData(2, fAnaTree);
        return;
    }
    fCounterH->Fill(iSelectionCounter); // At least one track (3)
    iSelectionCounter++;
    
    // loop over tracks and select good muons
    Int_t nGoodPosTrk = 0;
    Int_t nGoodNegTrk = 0;
    Int_t *idxPosTrk = new Int_t[nTracks];
    Int_t *idxNegTrk = new Int_t[nTracks];
    Int_t nGoodMUON = 0;
    
    for(Int_t iTrack(0); iTrack < nTracks; iTrack++) {
        // check if it is muon
        Bool_t isGoodMUON = GoodMUONTrack(iTrack);
        
        // if valid track
        if (isGoodMUON) {
            nGoodMUON++;    // update counters
            AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
            // set charge, increase counter and store indices
            if (track->Charge() > 0) {
                idxPosTrk[nGoodPosTrk] = iTrack;
                nGoodPosTrk++;
            } else  if (track->Charge() < 0) {
                idxNegTrk[nGoodNegTrk] = iTrack;
                nGoodNegTrk++;
            }
        }
    }
    
    ////////////////////////////////////////////
    // two track analysis
    ////////////////////////////////////////////
    
    // set type of analysis
    fAnaType = -1;
    if (nGoodMUON == 2) fAnaType = 0;
    
    // check valid analysis type
    if ( fAnaType == -1 ) {
        PostData(1, fOutputList);
        PostData(2, fAnaTree);
        return;
    }
    fCounterH->Fill(iSelectionCounter); // valid analysis (4)
    iSelectionCounter++;
    
    // We want to keep sign-like events as well
    /*
     // check unlike sign tracks
     if (!(nGoodPosTrk == 1 && nGoodNegTrk == 1)) {
     PostData(1, fAnaTree);
     PostData(2, fOutputList);
     delete [] idxNegTrk;
     delete [] idxPosTrk;
     return;
     }
     */
    Bool_t unlikeSign = kFALSE;
    if (nGoodPosTrk == 1 && nGoodNegTrk == 1) {fCounterH->Fill(iSelectionCounter); unlikeSign = true;} // exactly one positive and one negative track
    iSelectionCounter++;
    
    // set PID and compute trk+trk kinematics
    Double_t MuonMass = 0.105658; // GeV/c^2
    Double_t JpsiMass = 3.096916; // GeV/c^2
    // -- PID, fwd case
    fTrkRabs1 = fTrkRabs2 = 0;
    if (nGoodPosTrk == 1 && nGoodNegTrk == 1) {    // Opposite sign muons
        TrkTrkKine(idxPosTrk[0],idxNegTrk[0],MuonMass);
    }
    else if (nGoodPosTrk == 2 && nGoodNegTrk == 0) {
        TrkTrkKine(idxPosTrk[0],idxPosTrk[1],MuonMass);
    }
    else if (nGoodPosTrk == 0 && nGoodNegTrk == 2) {
        TrkTrkKine(idxNegTrk[0],idxNegTrk[1],MuonMass);
    }
    
    // TLorentzvector of the proton
    TLorentzVector protonP;
    protonP.SetVectM(TVector3(0, 0, 6500), 0.938);
    Double_t protonMomentum = protonP.P();
        
    Double_t rap = TrkTrk.Rapidity();
    /*
    if (fPeriod == 0 ) rap = -fTrkTrkY;
    else rap = fTrkTrkY;
     */
    fTrkTrkW = TMath::Sqrt(2*protonMomentum*TMath::Sqrt(Square(JpsiMass)+Square(fTrkTrkPt)) *TMath::Exp(-rap));
    //fTrkTrkZ = 6500*fTrkTrkEnergy * (1- TMath::Sqrt(1- (Square(fTrkTrkPt) + Square(JpsiMass))/Square(fTrkTrkEnergy)) )/ (Square(fTrkTrkW)/2);
    fTrkTrkZ = protonP* TrkTrk/(Square(fTrkTrkW)/2);     // using formula Z = (P.P(jpsi))/(P.q)
    fTrkTrkZBis = 2*protonMomentum*fTrkTrkEnergy/(Square(fTrkTrkW)/2);    // using formula Z = E(jpsi)/E(gamma)
    
    // clean up
    delete [] idxNegTrk;
    delete [] idxPosTrk;
    
    if (fIsMC) {
        fAnaTree->Fill();
        PostData(1, fOutputList);
        PostData(2, fAnaTree);
        return;
    }
    
    ////////////////////////////////////////////
    // info to determine exclusivity
    ////////////////////////////////////////////
    
    //  SPD
    fTracklets = fAOD->GetTracklets()->GetNumberOfTracklets();
    
    //  ZDC
    AliAODZDC *dataZDC = dynamic_cast<AliAODZDC*>(fAOD->GetZDCData());
    if(!dataZDC  && !fIsMC) {
        //if(!dataZDC) {
        PostData(1, fOutputList);
        PostData(2, fAnaTree);
        return;
    }
    if (unlikeSign) fCounterH->Fill(iSelectionCounter); // ZDC info is present
    iSelectionCounter++;
    fZNAfired = dataZDC->IsZNAfired();
    fZNCfired = dataZDC->IsZNCfired();
    fZNAEnergy = dataZDC->GetZNATowerEnergy()[0];
    fZNCEnergy = dataZDC->GetZNCTowerEnergy()[0];
    fZPAEnergy = dataZDC->GetZPATowerEnergy()[0];
    fZPCEnergy = dataZDC->GetZPCTowerEnergy()[0];
    fZNAgoodTiming = kFALSE; fZNCgoodTiming = kFALSE; fZNAgoodTiming2 = kFALSE; fZNCgoodTiming2 = kFALSE;
    fZPAgoodTiming = kFALSE; fZPCgoodTiming = kFALSE; fZPAgoodTiming2 = kFALSE; fZPCgoodTiming2 = kFALSE;
    for (Int_t i=0;i<4;i++) {
        fZNATDC[i] = dataZDC->GetZNATDCm(i);
        fZNCTDC[i] = dataZDC->GetZNCTDCm(i);
        fZPATDC[i] = dataZDC->GetZPATDCm(i);
        fZPCTDC[i] = dataZDC->GetZPCTDCm(i);
        
        if (fZNATDC[i] > -2. && fZNATDC[i] < 2.) fZNAgoodTiming = kTRUE;
        if (fZNCTDC[i] > -2. && fZNCTDC[i] < 2.) fZNCgoodTiming = kTRUE;
        if (fZNATDC[i] > -6. && fZNATDC[i] < 6.) fZNAgoodTiming2 = kTRUE;
        if (fZNCTDC[i] > -6. && fZNCTDC[i] < 6.) fZNCgoodTiming2 = kTRUE;
        
        if (fZPATDC[i] > -2. && fZPATDC[i] < 2.) fZPAgoodTiming = kTRUE;
        if (fZPCTDC[i] > -2. && fZPCTDC[i] < 2.) fZPCgoodTiming = kTRUE;
        if (fZPATDC[i] > -6. && fZPATDC[i] < 6.) fZPAgoodTiming2 = kTRUE;
        if (fZPCTDC[i] > -6. && fZPCTDC[i] < 6.) fZPCgoodTiming2 = kTRUE;
    }
    
    // V0
    AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
    if(!dataVZERO) {
        PostData(1, fOutputList);
        PostData(2, fAnaTree);
        return;
    }
    if (unlikeSign) fCounterH->Fill(iSelectionCounter); //  V0 info
    iSelectionCounter++;
    fV0ADecision = dataVZERO->GetV0ADecision();
    fV0CDecision = dataVZERO->GetV0CDecision();
    fV0ACounts = dataVZERO->GetNbPMV0A();
    fV0CCounts = dataVZERO->GetNbPMV0C();
    
    // Change here
    fV0AMultiplicity = dataVZERO->GetMTotV0A();
    fV0CMultiplicity = dataVZERO->GetMTotV0C();
    for (Int_t i=0;i<4;i++) {
        fV0ARingMultiplicity[i] = dataVZERO->GetMRingV0A(i);
        fV0CRingMultiplicity[i] = dataVZERO->GetMRingV0C(i);
    }
    fV0ATime = dataVZERO->GetV0ATime();
    fV0CTime = dataVZERO->GetV0CTime();
    fV0ABBNHits = 0;
    fV0CBBNHits = 0;
    fV0ABGNHits = 0;
    fV0CBGNHits = 0;
    fNV0C = 0;
    for (Int_t i = 0; i < 32; i++) {
        // online trigger
        //fV0BB[i] = dataVZERO->GetBBFlag(i);
        fV0AOnlineTrigger[i] = dataVZERO->GetBBFlag(32+i);
        fV0COnlineTrigger[i] = dataVZERO->GetBBFlag(i);
        // offline trigger
        fV0AOfflineTrigger[i] = dataVZERO->BBTriggerV0A(i);
        fV0COfflineTrigger[i] = dataVZERO->BBTriggerV0C(i);
        
        if (fV0AOnlineTrigger[i]) {
            fV0ABBNHits++;
            fV0ABBNHits++;
        }
        
        if (fV0COfflineTrigger[i]) fNV0C++;
    }
    
    for (Int_t i = 0; i < 64; i++) {
        fV0BG[i] = dataVZERO->GetBGFlag(i);
        fV0AvgPhi[i] = dataVZERO->GetVZEROAvgPhi(i);
        fV0EtaMin[i] = dataVZERO->GetVZEROEtaMin(i);
        fV0EtaMax[i] = dataVZERO->GetVZEROEtaMax(i);
        if (fV0BG[i]) {
            if (i<32) fV0CBGNHits++;
            else fV0ABGNHits++;
        }
    }
    
    // Number of matched tracks
    fV0CNMatched = 0;
    Int_t index1 = TagCellVZEROC(fTrkEta1, fTrkPhi1);
    if (index1 > 0) {if (fV0COfflineTrigger[index1]) fV0CNMatched++;}
    Int_t index2 = TagCellVZEROC(fTrkEta1, fTrkPhi1);
    if (index2 > 0) {if (fV0COfflineTrigger[index2]) fV0CNMatched++;}
    
    // check AD
    AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
    if(!dataAD){
        PostData(1, fOutputList);
        PostData(2, fAnaTree);
        return;
    }
    if (unlikeSign) fCounterH->Fill(iSelectionCounter); //  AD info
    iSelectionCounter++;
    fADADecision = dataAD->GetADADecision();
    fADCDecision = dataAD->GetADCDecision();
    fADABBNHits = 0;
    fADCBBNHits = 0;
    fADABGNHits = 0;
    fADCBGNHits = 0;
    for (Int_t i = 0; i < 16; i++) {
        fADBB[i] = dataAD->GetBBFlag(i);
        fADBG[i] = dataAD->GetBGFlag(i);
        if (fADBB[i]) {
            if (i<8) fADCBBNHits++;
            else fADABBNHits++;
        }
        if (fADBG[i]) {
            if (i<8) fADCBGNHits++;
            else fADABGNHits++;
        }
    }
    
    //Past-future protection maps
    fIR1Map = fAOD->GetHeader()->GetIRInt1InteractionMap();
    fIR2Map = fAOD->GetHeader()->GetIRInt2InteractionMap();
    
    // fill the tree
    fAnaTree->Fill();
    
    // post the data
    PostData(1, fOutputList);
    PostData(2, fAnaTree);
    
    // clean up
    // delete [] idxPosTrk;
    // delete [] idxNegTrk;
    
}
//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016::Terminate(Option_t *)
{
    cout << "\n\ndone\n\n" << endl;
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________

//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016::ProcessMCParticles()
{
    // process MC particles
    /*
     TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
     if (AODMCTrackArray == NULL) {std::cout << "\n\nNo MC object" << std::endl; return;}
     */
    
    //std::cout << "Process MC particles\n" << std::endl;
    AliMCEvent *AODMCData = dynamic_cast<AliMCEvent*>(MCEvent());
    if (!AODMCData) {std::cout << "\n\nNo MC object" << std::endl; return;}
    
    // This piece of code is taken from JGC
    // get MC particles
    Double_t MuonMass = 0.105658; // GeV/c^2
    Double_t JpsiMass = 3.096916; // GeV/c^2
    fMCTrkQ1 = fMCTrkQ2 = 0;
    Int_t n_mcp = AODMCData->GetNumberOfTracks();
    for(Int_t i_mcp = 0; i_mcp < n_mcp; i_mcp++) {
        AliMCParticle *mcp = static_cast<AliMCParticle*>(AODMCData->GetTrack(i_mcp));
        if(!mcp) {continue;}
        
        // assuming that j/psi was not stored, so muons have no mother
        if (mcp->GetMother() == -1) {
    
            if (mcp->PdgCode() == 13) {
                fMCTrkPt1 = mcp->Pt();
                fMCTrkEta1 = mcp->Eta();
                fMCTrkPhi1 = mcp->Phi();
                fMCTrkQ1--;
            }
            if (mcp->PdgCode() == -13) {
                fMCTrkPt2 = mcp->Pt();
                fMCTrkEta2 = mcp->Eta();
                fMCTrkPhi2 = mcp->Phi();
                fMCTrkQ2++;
            }
        }
    }   // end of for(Int_t i_mcp = 0; i_mcp < n_mcp; i_mcp++)
    if (fMCTrkQ1*fMCTrkQ2==-1) {
        TLorentzVector PosLV, NegLV;
        if (fMCTrkQ1==-1 && fMCTrkQ2==1) {
            PosLV.SetPtEtaPhiM(fMCTrkPt1, fMCTrkEta1, fMCTrkPhi1, MuonMass);
            NegLV.SetPtEtaPhiM(fMCTrkPt2, fMCTrkEta2, fMCTrkPhi2, MuonMass);
        }
        else {
            PosLV.SetPtEtaPhiM(fMCTrkPt2, fMCTrkEta2, fMCTrkPhi2, MuonMass);
            NegLV.SetPtEtaPhiM(fMCTrkPt1, fMCTrkEta1, fMCTrkPhi1, MuonMass);
        }
        TLorentzVector TrkTrk = NegLV+PosLV;
        fMCTrkTrkPt = TrkTrk.Pt();
        fMCTrkTrkPhi = TrkTrk.Phi();
        fMCTrkTrkY = TrkTrk.Rapidity();
        fMCTrkTrkM = TrkTrk.M();
        fMCTrkTrkEnergy = TrkTrk.Energy();
        Double_t rap;
        if (fPeriod == 0 ) rap = -fMCTrkTrkY;
        else rap = fMCTrkTrkY;
        fMCTrkTrkW = TMath::Sqrt(2*6500*JpsiMass*TMath::Exp(-rap));
        fMCTrkTrkZ = 6500*fMCTrkTrkEnergy * (1- TMath::Sqrt(1- (Square(fMCTrkTrkPt) + Square(JpsiMass))/Square(fMCTrkTrkEnergy)) )/ (Square(fMCTrkTrkW)/2);
    } else {
        return;
        fMCTrkTrkPt = -1;
        fMCTrkTrkPhi = -999;
        fMCTrkTrkY = -999;
        fMCTrkTrkM = -1;
        fMCTrkTrkEnergy = -1;
        fMCTrkTrkZ = -1;
    }
    
    fGenTree->Fill();
    PostData(3, fGenTree);
}
//_____________________________________________________________________________
