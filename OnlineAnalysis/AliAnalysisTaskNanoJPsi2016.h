/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskNanoJPsi2016_H
#define AliAnalysisTaskNanoJPsi2016_H

#include "AliAnalysisTaskSE.h"

class AliMuonTrackCuts; 	// Include class for standard muon tack cuts
class AliPIDResponse;
class AliMCEvent;
class AliAODMCParticle;


class AliAnalysisTaskNanoJPsi2016 : public AliAnalysisTaskSE
{
public:
	AliAnalysisTaskNanoJPsi2016();
	AliAnalysisTaskNanoJPsi2016(const char *name);
	virtual                 ~AliAnalysisTaskNanoJPsi2016();
	
	virtual void            UserCreateOutputObjects();
	virtual void            UserExec(Option_t* option);
	virtual void            Terminate(Option_t* option);
	virtual void   			NotifyRun();// Implement the Notify run to search for the new parameters at each new runs
	void TrkTrkKine(Int_t pos, Int_t neg, Double_t mass);
	void CheckTrigger(Bool_t &trig);
	void SetPeriod(Int_t period);
	void SetMC(Bool_t flag);
    Int_t TagCellVZEROC(Double_t ItsEta, Double_t ItsPhi);
	Bool_t GoodMUONTrack(Int_t iTrack);
	Bool_t GoodCentralTrack(Int_t iTrack);
	Int_t GetMassHypothesis(Int_t *idxPosTrk, Int_t *idxNegTrk);
	
	virtual void ProcessMCParticles();
	
	AliMuonTrackCuts* 		fMuonTrackCuts; 					// Use the class as a data member
	
private:
	Int_t fPeriod;
	Bool_t fIsMC;
	
	AliPIDResponse *fPIDResponse;
	
	AliAODEvent*            fAOD;       //! input event
	AliMCEvent*             fMCEvent;       //! corresponding MC event
	
	TList*                  fOutputList; //! output list
	TH1F*                   fCounterH; //! counter for events passing each cut
	TH1F*                   fMuonTrackCounterH; //! counter for tracks passing each cut
	TH1F*                   fCentralTrackCounterH; //! counter for tracks passing each cut
	TH1F*                   fTriggerCounterFwdH; //! counter for triggersper run
	TH1F*                   fTriggerCounterCentH; //! counter for triggersper run
	TH1F*                   fTriggerCounterSemiFwdH; //! counter for triggersper run
	
	TTree *fAnaTree; //! analysis tree
	Int_t fRunNum;
	Int_t fOrbitNum;
	Int_t fBunchCrossNum;
	UInt_t fL0inputs;
	Int_t fTracklets;
	Int_t fAnaType;
	Bool_t fZNAfired;
	Bool_t fZNCfired;
	Double_t fZNATDC[4];
	Double_t fZNCTDC[4];
	Bool_t fZNAgoodTiming;
	Bool_t fZNCgoodTiming;
	Bool_t fZNAgoodTiming2;
	Bool_t fZNCgoodTiming2;
	Double_t fZNCEnergy;
	Double_t fZNAEnergy;
	Double_t fZPCEnergy;
	Double_t fZPAEnergy;
	Double_t fZPATDC[4];
	Double_t fZPCTDC[4];
	Bool_t fZPAgoodTiming;
	Bool_t fZPCgoodTiming;
	Bool_t fZPAgoodTiming2;
	Bool_t fZPCgoodTiming2;
	Int_t fV0ADecision;
	Int_t fV0CDecision;
	Int_t fV0ACounts;
	Int_t fV0CCounts;
    Int_t fNV0C;
	Double_t fV0AMultiplicity;
	Double_t fV0CMultiplicity;
	Double_t fV0ARingMultiplicity[4];
	Double_t fV0CRingMultiplicity[4];
	Double_t fV0ATime;
	Double_t fV0CTime;
	Int_t fV0ABBNHits;
	Int_t fV0CBBNHits;
	Int_t fV0ABGNHits;
	Int_t fV0CBGNHits;
	//Bool_t fV0BB[64];
    Bool_t fV0AOnlineTrigger[32];
    Bool_t fV0COnlineTrigger[32];
    Bool_t fV0AOfflineTrigger[32];
    Bool_t fV0COfflineTrigger[32];
	Bool_t fV0BG[64];
    Int_t fV0CNMatched;
	Double_t fV0AvgPhi[64];
	Double_t fV0EtaMin[64];
	Double_t fV0EtaMax[64];
	Int_t fADADecision;
	Int_t fADCDecision;
	Int_t fADABBNHits;
	Int_t fADCBBNHits;
	Int_t fADABGNHits;
	Int_t fADCBGNHits;
	Bool_t fADBB[16];
	Bool_t fADBG[16];
	TBits fIR1Map;
	TBits fIR2Map;
    TLorentzVector TrkTrk;
	Double_t fTrkTrkPt;
	Double_t fTrkTrkPhi;
	Double_t fTrkTrkY;
	Double_t fTrkTrkM;
    Double_t fTrkTrkEnergy;
    Double_t fTrkTrkZ;
    Double_t fTrkTrkZBis;
    Double_t fTrkTrkW;
	Double_t fTrkPt1;
	Double_t fTrkPt2;
	Double_t fTrkEta1;
	Double_t fTrkEta2;
	Double_t fTrkPhi1;
	Double_t fTrkPhi2;
	Double_t fTrkQ1;
	Double_t fTrkQ2;
	Double_t fTrkRabs1;
	Double_t fTrkRabs2;
    std::vector<Double_t> fTrkProtonEta;
    std::vector<Double_t> fTrkProtonPhi;
    Int_t fEventLabel;
	
	TTree *fGenTree; //! MC tree
	Double_t fMCTrkTrkPt;
	Double_t fMCTrkTrkPhi;
	Double_t fMCTrkTrkY;
	Double_t fMCTrkTrkM;
    Double_t fMCTrkTrkEnergy;
    Double_t fMCTrkTrkZ;
    Double_t fMCTrkTrkW;
	Double_t fMCTrkPt1;
	Double_t fMCTrkPt2;
	Double_t fMCTrkEta1;
	Double_t fMCTrkEta2;
	Double_t fMCTrkPhi1;
	Double_t fMCTrkPhi2;
	Int_t fMCTrkQ1;
	Int_t fMCTrkQ2;
    std::vector<Double_t> fMCTrkProtonEta;
    std::vector<Double_t> fMCTrkProtonPhi;
    std::vector<Double_t> fMCTrkProtonP;
    std::vector<Bool_t> fMCIsPrimary;
    std::vector<Int_t> fMCPdg;
    Double_t fMCTrkProtonRemnantsM;
    Double_t fMCProtonEta;
    Int_t fMCnProton;
	
	AliAnalysisTaskNanoJPsi2016(const AliAnalysisTaskNanoJPsi2016&); // not implemented
	AliAnalysisTaskNanoJPsi2016& operator=(const AliAnalysisTaskNanoJPsi2016&); // not implemented
	
	ClassDef(AliAnalysisTaskNanoJPsi2016, 1);
};

#endif
