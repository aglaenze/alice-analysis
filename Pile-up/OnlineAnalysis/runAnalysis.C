// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly
#include "AliAnalysisTaskMu.h"

// LHC16r
const int LHC16rRuns[] = {265589, 265594, 265596, 265607, 265669, 265691, 265694, 265697, 265698, 265700, 265701, 265709, 265713, 265714, 265740, 265741, 265742, 265744, 265746, 265754, 265756, 265785, 265787, 265788, 265789, 265792, 265795, 265797, 265840, 266022, 266023, 266025, 266034, 266074, 266076, 266081, 266084, 266085, 266086, 266117, 266187, 266189, 266190, 266193, 266196, 266197, 266208, 266234, 266235, 266296, 266299, 266300, 266304, 266305, 266312, 266316, 266318};
const int n16r = sizeof(LHC16rRuns) / sizeof(int);

// LHC16s
const int LHC16sRuns[]  = {266439, 266441, 266472, 266479, 266480, 266487, 266514, 266516, 266518, 266520, 266522, 266523, 266525, 266533, 266534, 266539, 266543, 266549, 266584, 266587, 266588, 266591, 266593, 266595, 266613, 266614, 266618, 266621, 266630, 266657, 266658, 266659, 266665, 266668, 266669, 266674, 266676, 266702, 266703, 266706, 266708, 266775, 266776, 266800, 266805, 266807, 266857, 266878, 266880, 266882, 266883, 266885, 266886, 266912, 266915, 266940, 266942, 266943, 266944, 266988, 266993, 266994, 266997, 266998, 267020, 267022, 267062, 267063, 267067, 267070, 267072, 267077, 267109, 267110, 267130, 267131};
const int n16s = sizeof(LHC16sRuns) / sizeof(int);

void runAnalysis(Int_t period = 0)
// period = 0 => 2016 r,
//        = 1 => 2016 s
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t local = kTRUE;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t gridTest = kFALSE;
    
    // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif
     
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskMu");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);



    // compile the class and load the add task macro
    // here we have to differentiate between using the just-in-time compiler
    // from root6, or the interpreter of root5
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AliAnalysisTaskMu.cxx++g");
    char txt_cmd[120];
    sprintf(txt_cmd,"AddMyTask.C(%d)", period);
    AliAnalysisTaskMu *task = reinterpret_cast<AliAnalysisTaskMu*>(gInterpreter->ExecuteMacro(txt_cmd));
#else
    gROOT->LoadMacro("AliAnalysisTaskMu.cxx++g");
    gROOT->LoadMacro("AddMyTask.C");
    AliAnalysisTaskMu *task = AddMyTask();
#endif


    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);
    
    //cout << "\n\n\nStarting\n\n\n" << endl;
    
    char tmpstr[120];

    if(local) {
        if (period == 1) {
            cout << "Period has to be 0 in local" << endl;
            return;
        }
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        // add a few files to the chain (change this so that your local files are added)
        //chain->Add("AliAOD.root");
        
        chain->Add("AliAOD_pass1_UD_266076_1.root");
        //chain->Add("AliAOD_pass1_UD_266076_2.root");
        //chain->Add("AliAOD_pass1_UD_266076_3.root");
         
        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);
        //cout << "\n\n\nStarting\n\n\n" << endl;
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskMu.cxx AliAnalysisTaskMu.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskMu.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20191015_ROOT6-1");
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        // select the input data
        if (period == 0) {
            alienHandler->SetGridDataDir("/alice/data/2016/LHC16r");
            alienHandler->SetDataPattern("*pass1_UD/PWGUD/UD_pPb_AOD/143_20191104-1006/*AliAOD.UPCNano.root"); // with ITS tracks
              // MC has no prefix, data has prefix 000
            alienHandler->SetRunPrefix("000");
              // runnumber
            for (int k = 0; k<n16r; k++) alienHandler->AddRunNumber(LHC16rRuns[k]);
            // working dir
            alienHandler->SetGridWorkingDir("LHC16r_ana_mu");
            alienHandler->SetExecutable("LHC16r_Task_mu.sh");
            alienHandler->SetJDLName("LHC16r_Task_mu.jdl");
        }
        else if (period == 1) {
            alienHandler->SetGridDataDir("/alice/data/2016/LHC16s");
            alienHandler->SetDataPattern("*pass1_UD/PWGUD/UD_pPb_AOD/144_20191104-1006/*AliAOD.UPCNano.root"); // with ITS tracks
            // MC has no prefix, data has prefix 000
            alienHandler->SetRunPrefix("000");
              // runnumber
            for (int k = 0; k<n16s; k++) alienHandler->AddRunNumber(LHC16sRuns[k]);
            // working dir
            alienHandler->SetGridWorkingDir("LHC16s_ana_mu");
            alienHandler->SetExecutable("LHC16s_Task_mu.sh");
            alienHandler->SetJDLName("LHC16s_Task_mu.jdl");
        }
        else {
            cout << " not a valid option ... bye!" << endl;
            return;
        }

        
        // number of files per subjob
        //alienHandler->SetSplitMaxInputFileNumber(40);`
        alienHandler->SetSplitMaxInputFileNumber(20);
        // specify how many seconds your job may take
        alienHandler->SetTTL(10000);
        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate") 
        // (see below) mode, set SetMergeViaJDL(kFALSE) 
        // to collect final results
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(kTRUE);

        // define the output folders
        //alienHandler->SetGridWorkingDir("myWorkingDir");
        alienHandler->SetGridOutputDir("myOutputDir");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if (gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            alienHandler->SetRunMode("full");
            //alienHandler->SetRunMode("terminate");
            mgr->StartAnalysis("grid");
        }
    }
}
