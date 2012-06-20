
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include "TString.h"
#include "TROOT.h"
#include "../../Tools/goodrun.cc"

void skimLeptonTree (TString runlist, TString inputDir, TString outputFile, bool tp)
{

    set_goodrun_file(runlist.Data());

    TFile* fin = new TFile(inputDir+"/merged.root");
    TTree* ch=(TTree*)fin->Get("leptons"); 
    if (ch==0x0) return; 

    TFile *newfile= new TFile(inputDir+"/"+outputFile,"recreate");
    TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");

    unsigned int run_ = 0;
    unsigned int lumi_ = 0;
    ch->SetBranchAddress( "run"           , &run_     );     
    ch->SetBranchAddress( "lumi"          , &lumi_     );     

    unsigned int eventSelection_ = 0;
    ch->SetBranchAddress( "eventSelection"          , &eventSelection_     );

    for(int ievt = 0; ievt < ch->GetEntries() ;ievt++) {
        ch->GetEntry(ievt); 

            if (tp) {
                if (!((eventSelection_ & (1<<0)) == (1<<0)  ||
                    (eventSelection_ & (1<<1)) == (1<<1))) continue;
            } else {
                if (!((eventSelection_ & (1<<4)) == (1<<4)  ||
                    (eventSelection_ & (1<<5)) == (1<<5))) continue;
            }
            
                if (!goodrun_json(run_, lumi_)) continue;


        evt_tree->Fill();
    }

newfile->cd(); 
evt_tree->Write(); 
newfile->Close();
}  

