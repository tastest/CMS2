
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include "TString.h"
#include "TROOT.h"
#include "../../Tools/goodrun.cc"

void skimLeptonTree (TString runlist, TString inputFile, TString outputFile)
{

    set_goodrun_file(runlist.Data());

    TFile* fin = new TFile(inputFile);
    TTree* ch=(TTree*)fin->Get("leptons"); 
    if (ch==0x0) return; 

    TFile *newfile= new TFile(outputFile,"recreate");
    TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");

    unsigned int run_ = 0;
    unsigned int lumi_ = 0;
    ch->SetBranchAddress( "run"           , &run_     );     
    ch->SetBranchAddress( "lumi"          , &lumi_     );     

    for(int ievt = 0; ievt < ch->GetEntries() ;ievt++) {
        ch->GetEntry(ievt); 

                if (!goodrun_json(run_, lumi_)) continue;


        evt_tree->Fill();
    }

newfile->cd(); 
evt_tree->Write(); 
newfile->Close();
}  

