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

/* AliAnaysisTaskMu
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

// c++ headers
#include <iostream>
#include <fstream>

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTriggerScalers.h"
#include "AliLHCData.h"

// my headers
#include "AliAnalysisTaskMu.h"

class AliAnalysisTaskMu;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskMu) // classimp: necessary for root

//______________________________________________________________________________
Bool_t IsMainSatelliteCollision(const char* triggerClassName)
{
    TString tcn(triggerClassName);
    tcn.ToUpper();
    
    return (tcn.Contains("-S-") ||  tcn.Contains("-SC-") || tcn.Contains("-SA-"));
}

AliAnalysisTaskMu::AliAnalysisTaskMu() : AliAnalysisTaskSE(),
fAOD(0), fOutputTree(0), fRunNum(0),
fPeriod(0),
fZNCEnergy(0),
fZNAEnergy(0),
fZNAgoodTiming(0),
fZNCgoodTiming(0),
fZNAgoodTiming2(0),
fZNCgoodTiming2(0),
fV0ADecision(-10),
fV0CDecision(-10),
fADADecision(-10),
fADCDecision(-10),
f0VBA(0),
f0VGA(0),
f0UBA(0),
f0VBC(0),
f0UGC(0),
f0UBC(0),
f0VC5(0),
f0SH2(0),
/*
fCtrue(0),
 fCut14on(0),
 fCut14offDAA(0),
 fCut14offDVA(0),
 fCut14offDVC(0),
 fCut14offZN(0),
 fCut14off(0),
 fCut14(0),
 fMu(-1),
 */
fTrdNTracks(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskMu::AliAnalysisTaskMu(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputTree(0), fRunNum(0),
fPeriod(0),
fZNCEnergy(0),
fZNAEnergy(0),
fZNAgoodTiming(0),
fZNCgoodTiming(0),
fZNAgoodTiming2(0),
fZNCgoodTiming2(0),
fV0ADecision(-10),
fV0CDecision(-10),
fADADecision(-10),
fADCDecision(-10),
f0VBA(0),
f0VGA(0),
f0UBA(0),
f0VBC(0),
f0UGC(0),
f0UBC(0),
f0VC5(0),
f0SH2(0),
/*
 fCtrue(0),
 fCut14on(0),
 fCut14offDAA(0),
 fCut14offDVA(0),
 fCut14offDVC(0),
 fCut14offZN(0),
 fCut14off(0),
 fCut14(0),
 fMu(-1),
 */
fTrdNTracks(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
    // this chain is created by the analysis manager, so no need to worry about it,
    // it does its work automatically
    DefineOutput(1, TTree::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
    // you can add more output objects by calling DefineOutput(2, classname::Class())
    // if you add more output objects, make sure to call PostData for all of them, and to
    // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskMu::~AliAnalysisTaskMu()
{
    // destructor
    if(fOutputTree) {
        delete fOutputTree;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskMu::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputTree = new TTree("fAnaTree", "fAnaTree");          // this is a list which will contain all of your histograms
    fOutputTree->Branch("fRunNum", &fRunNum, "fRunNum/I");
    // at the end of the analysis, the contents of this list are written
    // to the output file
    fOutputTree->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
    // Booleans here
    fOutputTree->Branch("f0VBA", &f0VBA,"f0VBA/O");
    fOutputTree->Branch("f0VGA", &f0VGA,"f0VGA/O");
    fOutputTree->Branch("f0UBA", &f0UBA,"f0UBA/O");
    fOutputTree->Branch("f0VBC", &f0VBC,"f0VBC/O");
    fOutputTree->Branch("f0UGC", &f0UGC,"f0UGC/O");
    fOutputTree->Branch("f0UBC", &f0UBC,"f0UBC/O");
    fOutputTree->Branch("f0VC5", &f0VC5,"f0VC5/O");
    fOutputTree->Branch("f0SH2", &f0SH2,"f0SH2/O");
    
    fOutputTree->Branch("fTrdNTracks", &fTrdNTracks, "fTrdNTracks/I");
    fOutputTree->Branch("fZNCEnergy", &fZNCEnergy, "fZNCEnergy/D");
    fOutputTree->Branch("fZNAEnergy", &fZNAEnergy, "fZNAEnergy/D");
    fOutputTree->Branch("fZNATDC", &fZNATDC[0], "fZNATDC[4]/D");
    fOutputTree->Branch("fZNCTDC", &fZNCTDC[0], "fZNCTDC[4]/D");
    fOutputTree->Branch("fZNAgoodTiming", &fZNAgoodTiming, "fZNAgoodTiming/O");
    fOutputTree->Branch("fZNCgoodTiming", &fZNCgoodTiming, "fZNCgoodTiming/O");
    fOutputTree->Branch("fZNAgoodTiming2", &fZNAgoodTiming2, "fZNAgoodTiming2/O");
    fOutputTree->Branch("fZNCgoodTiming2", &fZNCgoodTiming2, "fZNCgoodTiming2/O");
    fOutputTree->Branch("fV0ADecision", &fV0ADecision, "fV0ADecision/I");
    fOutputTree->Branch("fV0CDecision", &fV0CDecision, "fV0CDecision/I");
    fOutputTree->Branch("fADADecision", &fADADecision, "fADADecision/I");
    fOutputTree->Branch("fADCDecision", &fADCDecision, "fADCDecision/I");
    /*
     fOutputTree ->Branch("fCtrue",       &fCtrue,       "fCtrue/O"      );
     fOutputTree->Branch("fCut14on", &fCut14on, "fCut14on/O");
     fOutputTree->Branch("fCut14offDAA", &fCut14offDAA, "fCut14offDAA/O");
     fOutputTree->Branch("fCut14offDVA", &fCut14offDVA, "fCut14offDVA/O");
     fOutputTree->Branch("fCut14offDVC", &fCut14offDVC, "fCut14offDVC/O");
     fOutputTree->Branch("fCut14offZN", &fCut14offZN, "fCut14offZN/O");
     fOutputTree->Branch("fCut14off", &fCut14off, "fCut14off/O");
     fOutputTree->Branch("fCut14", &fCut14, "fCut14/O");
    fOutputTree->Branch("fMu", &fMu, "fMu/D");
     */
    
    //fOutputTree->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
    // if requested (dont worry about this now)
    
    // example of a histogram
    //fHistPt = new TH1F("fHistPt", "fHistPt", 100, 0, 10);       // create your histogra
    //fOutputTree->Add(fHistPt);          // don't forget to add it to the list! the list will be written to file, so if you want
    // your histogram in the output file, add it to the list!
    
    PostData(1, fOutputTree);           // postdata will notify the analysis manager of changes / updates to the
    // fOutputTree object. the manager will in the end take care of writing your output to file
    // so it needs to know what's in the output
}

//_____________________________________________________________________________
void AliAnalysisTaskMu::SetPeriod(Int_t period)
{
    // period = 0 => 2016 s, = 1 => 2016 r
    fPeriod = period;
}

void AliAnalysisTaskMu::CheckTrigger(Bool_t &isTriggered)
// checks if event is triggered according to period and analysis type
{
    // Initialise: 0 = fwd, 1 = cent, 2 = semi-fwd
    
    // read trigger info
    
    TString trigger = fAOD->GetFiredTriggerClasses();
    if (trigger.Contains("CTRUE-B-NOPF-CENT")) {isTriggered = kTRUE;}
}

//_____________________________________________________________________________
void AliAnalysisTaskMu::UserExec(Option_t *)
{
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you
    // have access to the current event.
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file
    // there's another event format (ESD) which works in a similar wya
    // but is more cpu/memory unfriendly. for now, we'll stick with aod's
    if(!fAOD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event
    // example part: i'll show how to loop over the tracks in an event
    fRunNum = fAOD->GetRunNumber();
    fTrdNTracks = fAOD->GetNumberOfTrdTracks();
    Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
    if (iTracks<1) {
        PostData(1, fOutputTree);
        return;
    }
    
    // read trigger info
    Bool_t isTriggered = kFALSE;
    CheckTrigger(isTriggered);
    
    if (!isTriggered) {
        PostData(1, fOutputTree);
        return;
    }
    //std::cout << trigger << std::endl;
    
    
    //  ZDC
    AliAODZDC *dataZDC = dynamic_cast<AliAODZDC*>(fAOD->GetZDCData());
    if(!dataZDC) {
        PostData(1, fOutputTree);
        return;
    }
    fZNAEnergy = dataZDC->GetZNATowerEnergy()[0];
    fZNCEnergy = dataZDC->GetZNCTowerEnergy()[0];
    
    fZNAgoodTiming = kFALSE; fZNCgoodTiming = kFALSE; fZNAgoodTiming2 = kFALSE; fZNCgoodTiming2 = kFALSE;
    for (Int_t i=0;i<4;i++) {
        fZNATDC[i] = dataZDC->GetZNATDCm(i);
        fZNCTDC[i] = dataZDC->GetZNCTDCm(i);
        
        if (fZNATDC[i] > -2. && fZNATDC[i] < 2.) fZNAgoodTiming = kTRUE;
        if (fZNCTDC[i] > -2. && fZNCTDC[i] < 2.) fZNCgoodTiming = kTRUE;
        if (fZNATDC[i] > -8. && fZNATDC[i] < 8.) fZNAgoodTiming2 = kTRUE;
        if (fZNCTDC[i] > -8. && fZNCTDC[i] < 8.) fZNCgoodTiming2 = kTRUE;
    }
    
    // V0
    AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
    if(!dataVZERO) {
        PostData(1, fOutputTree);
        return;
    }
    fV0ADecision = dataVZERO->GetV0ADecision();
    fV0CDecision = dataVZERO->GetV0CDecision();

    // check AD
    AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
    if(!dataAD) {
        PostData(1, fOutputTree);
        return;
    }
    fADADecision = dataAD->GetADADecision();
    fADCDecision = dataAD->GetADCDecision();
    
    
    // Check L0 input
    Int_t fL0inputs  = fAOD->GetHeader()->GetL0TriggerInputs();
    
    AliAnalysisTriggerScalers ts(fRunNum,"alien://folder=/alice/data/2016/OCDB");
    const char* triggerClassName("CTRUE-B-NOPF-CENT");
    //const char* triggerClassName("CINT7-B-NOPF-CENT");
    
    /*
    //AliAnalysisTriggerScalerItem* b = ts.GetTriggerScaler(fRunNum,"L0B",triggerClassName);
    AliAnalysisTriggerScalerItem* b = ts.GetTriggerScaler(fRunNum,"LMB",triggerClassName);
    
    if (!b) {
        AliWarning(Form("Could not get L0B for trigger %s for run %d",triggerClassName,fRunNum));
        PostData(1, fOutputTree);
        return;
    }
     */
    
    // Numbering of trigger inputs
    // Checked for one run in LHC16r and one in LHC16s, assumed to be valid for all other runs in LHC16r and LHC16s
    Int_t inputId_0VBA = 1;
    Int_t inputId_0VGA = 4;
    Int_t inputId_0UBA = 6;
    Int_t inputId_0VBC = 2;
    Int_t inputId_0UGC = 12;
    Int_t inputId_0UBC = 7;
    Int_t inputId_0VC5 = 8;
    Int_t inputId_0SH2(0);
    if (fPeriod==0) inputId_0SH2 = 9;
    else if (fPeriod==1) inputId_0SH2 = 17;
    
    // get L0 trigger flags
    f0VBA = fL0inputs & 1 << (inputId_0VBA-1);
    f0VGA = fL0inputs & 1 << (inputId_0VGA-1);
    f0UBA = fL0inputs & 1 << (inputId_0UBA-1);
    f0VBC = fL0inputs & 1 << (inputId_0VBC-1);
    f0UGC = fL0inputs & 1 << (inputId_0UGC-1);
    f0UBC = fL0inputs & 1 << (inputId_0UBC-1);
    f0VC5 = fL0inputs & 1 << (inputId_0VC5-1);
    f0SH2 = fL0inputs & 1 << (inputId_0SH2-1);
    
    /*
    // study the behaviour of the online tigger
    fCut14on = f0UBA || f0VBA;
    // study the ADA decision
    fCut14offDAA = !(f0UBA || f0VBA) && ((fADADecision != 0));
    // study the V0A decision
    fCut14offDVA = !(f0UBA || f0VBA || (fADADecision != 0) ) && (fV0ADecision != 0);
    // study the V0C decision
    fCut14offDVC = !(f0UBA || f0VBA || (fADADecision != 0) || (fV0ADecision != 0)) && (fADCDecision != 0);
    // study the ZN decision
    fCut14offZN = !(f0UBA || f0VBA || (fADADecision != 0) || (fV0ADecision != 0) || (fADCDecision != 0) ) && ( (fZNAEnergy>400) || (fZNCEnergy>18) );
    // study all offline vetoes
    fCut14off = !(f0UBA || f0VBA) && ((fADADecision != 0) || (fV0ADecision != 0) || (fADCDecision != 0) || (fZNAEnergy>400) || (fZNCEnergy>18));
    // all vetoes at the same time
    // The cut to be actually used in the analysis is Cut14. Cuts Cut14on and Cut14off are used for a cross check. Other cuts are used only for illustration.
    fCut14 = f0UBA || f0VBA || (fADADecision != 0) || (fV0ADecision != 0) || (fADCDecision != 0) || (fZNAEnergy>400) || (fZNCEnergy>18);
     */
    
    
    /*
    Bool_t mainSat = IsMainSatelliteCollision(triggerClassName);
    AliLHCData* lhc = static_cast<AliLHCData*>(ts.GetOCDBObject("GRP/GRP/LHCData",fRunNum));
    Int_t numberOfInteractingBunches = ts.NumberOfInteractingBunches(*lhc,fRunNum,mainSat);

     std::cout << "b->ValueCorrectedForDownscale() " << b->ValueCorrectedForDownscale() << std::endl;
     std::cout << "b->Duration() " << b->Duration() << std::endl;
     std::cout << "numberOfInteractingBunches " << numberOfInteractingBunches << std::endl;
    
    //fMu = ts.Mu((double)b->ValueCorrectedForDownscale()/(double)b->Duration(),(double)numberOfInteractingBunches);
    //fMu = ts.Mu(b->ValueCorrectedForDownscale()/b->Duration(),numberOfInteractingBunches);
    //std::cout << "mu = " << fMu << std::endl;
    //if (fMu < 0.001) std::cout << "\n\n\nrunNumber = " << fRunNum << std::endl << std::endl << std::endl << std::endl << std::endl;
    */
    
    
    fOutputTree->Fill();
    PostData(1, fOutputTree);                           // stream the results the analysis of this event to
    // the output manager which will take care of writing
    // it to a file
    
}
//_____________________________________________________________________________
void AliAnalysisTaskMu::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
