/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskMu_H
#define AliAnalysisTaskMu_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTriggerScalers;
class AliLHCData;

class AliAnalysisTaskMu : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskMu();
    AliAnalysisTaskMu(const char *name);
    virtual                 ~AliAnalysisTaskMu();
    
    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    
    void CheckTrigger(Bool_t &trig);
    void SetPeriod(Int_t period);
    
private:
    AliAODEvent*            fAOD;           //! input event
    TTree*                  fOutputTree;    //! output list
    
    Int_t fPeriod;
    Int_t fRunNum;
    UInt_t fL0inputs;
    Bool_t f0VBA;
    Bool_t f0VGA;
    Bool_t f0UBA;
    Bool_t f0VBC;
    Bool_t f0UGC;
    Bool_t f0UBC;
    Bool_t f0VC5;
    Bool_t f0SH2;
    Double_t fZNCEnergy;
    Double_t fZNAEnergy;
    Double_t fZNATDC[4];
    Double_t fZNCTDC[4];
    Bool_t fZNAgoodTiming;
    Bool_t fZNCgoodTiming;
    Bool_t fZNAgoodTiming2;
    Bool_t fZNCgoodTiming2;
    Int_t fV0ADecision;
    Int_t fV0CDecision;
    Int_t fADADecision;
    Int_t fADCDecision;
    /*
     Bool_t fCtrue;
     Bool_t fCut14on;
     Bool_t fCut14offDAA;
     Bool_t fCut14offDVA;
     Bool_t fCut14offDVC;
     Bool_t fCut14offZN;
     Bool_t fCut14off;
     Bool_t fCut14;
     Double_t fMu;
     */
    Int_t fTrdNTracks;
    
    AliAnalysisTaskMu(const AliAnalysisTaskMu&); // not implemented
    AliAnalysisTaskMu& operator=(const AliAnalysisTaskMu&); // not implemented
    
    ClassDef(AliAnalysisTaskMu, 1);
};

#endif
