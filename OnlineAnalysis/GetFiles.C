
// header files
// c++ headers
#include <iostream>

// root headers
#include "TGrid.h"
#include "TGridResult.h"
#include "TFileMerger.h"

// main program
void GetFiles(Int_t opt = 1, Bool_t MC = kFALSE)
// get a specifc set of files from grid
{
	TString MCtype = "kCohJpsiToMu";
	//TString MCtype = "kIncohJpsiToMu";
	//TString MCtype = "kIncohPsi2sToMu";
	//TString MCtype = "kIncohPsi2sToMuPi";
	//TString MCtype = "kTwoGammaToMuLow";
	//TString MCtype = "kTwoGammaToMuMedium";
	
	//connect to the GRID
	TGrid::Connect("alien://");
	
	TGridResult* result = NULL;
	if (opt == 0) {
		if (!MC) result = gGrid->Query("/alice/cern.ch/user/a/aglaenze/LHC16r_ana/myOutputDir/","AnalysisResults.root");
		else result = gGrid->Query("/alice/cern.ch/user/a/aglaenze/LHC16r_MC_"+MCtype+"/myOutputDir/","AnalysisResults.root");
	}
	if (opt == 1) {
		if(!MC) result = gGrid->Query("/alice/cern.ch/user/a/aglaenze/LHC16s_ana/myOutputDir/","AnalysisResults.root");
		else result = gGrid->Query("/alice/cern.ch/user/a/aglaenze/LHC16s_MC_"+MCtype+"/myOutputDir/","AnalysisResults.root");
	}
	TFileMerger m;
	if (opt == 0) {
		if (!MC) m.OutputFile("AnalysisResults_LHC16r.root");
		else m.OutputFile("AnalysisResults_LHC16r_MC_"+MCtype+".root");
	}
	if (opt == 1) {
		if (!MC) m.OutputFile("AnalysisResults_LHC16s.root");
		else m.OutputFile("AnalysisResults_LHC16s_MC_"+MCtype+".root");
	}
	Int_t i = 0;
	//Loop over the TGridResult entries and add them to the TFileMerger
	while(result->GetKey(i,"turl")) {
		cout << " adding " << result->GetKey(i,"turl") << endl;
		m.AddFile(result->GetKey(i,"turl"));
		i++;
	}
	//Merge
	if(i) m.Merge();
	
}




