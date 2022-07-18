
// header files
// c++ headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <dirent.h>
#include <string.h>
#include <stdio.h>

// root headers
#include "TFileMerger.h"

int GetNumberOfFiles(TString path, TString name) {
    Int_t num = 0;
    struct dirent **namelist;
    Int_t n = scandir(path, &namelist, 0, alphasort);
    if (n < 1) {cout << "empty folder" << endl; return 0;}
    else {
        while (n--) {
            //if (strstr(namelist[n]->d_name, name) != NULL) num++;
            string filename = namelist[n]->d_name;
            if ( filename.rfind(name, 0) == 0) num++;
        }
        return num;
    }
}

// main program
int MergeAODs()
// get a specifc set of files from grid
{
    TFileMerger m;
    m.OutputFile("AliAOD_MC_UD_266025_rapgap.root", "RECREATE");
    
    int num = GetNumberOfFiles(".", "AliAOD-");
    //int num = 21;
    for (int i = 0; i < num; i++) {
        TFile* f1 = new TFile(Form("AliAOD-%d.root", i), "READ");
        //TFile* f1 = new TFile(Form("test%d/AliAOD.root", i), "READ");
        m.AddFile(f1);
        //cout << 
    }
    bool test = m.Merge();
    if (!test) return -1;
    
    return 0;
	
}




