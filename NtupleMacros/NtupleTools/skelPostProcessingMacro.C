#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TList.h"
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"
#include "TROOT.h"
#include <iostream>
#include <vector>
#include "TChain.h"
#include "TString.h"
#include "TObjArray.h"


using namespace std;

void doPostProcessing(TString infname, TString outfile, Int_t events, 
		      Float_t xsec, Float_t kfactor,
		      Float_t filt_eff=1) {
  
  cout << "Processing File " << infname << endl;
  
  TFile *f = TFile::Open(infname.Data(), "READ");
  if (! f || f->IsZombie()) {
    cout << "File does not exist!" << endl;
    return;
  }
  
  TTree* t = (TTree*)f->Get("Events");
  if (! t || t->IsZombie()) {
    cout << "Tree does not exist!" << endl;
    return;
  }
        
  //-------------------------------------------------------------
  // Removes the branches (if they exist) that we want to replace
  //-------------------------------------------------------------`
  //evt_xsec_excl
  TString bName = t->GetAlias("evt_xsec_excl");
  //cout << "evt_xsec_excl " << bName << endl;
  if(bName != "") {
    bName.ReplaceAll(".obj", "*");
    t->SetBranchStatus(bName.Data(), 0); 
  }

  //evt_xsec_incl
  bName = t->GetAlias("evt_xsec_incl");
  //cout << "evt_xsec_incl " << bName << endl;
  if(bName != "") {
    bName.ReplaceAll(".obj", "*");
    t->SetBranchStatus(bName.Data(), 0);   
  }
  
  //evt_kfactor
  bName = t->GetAlias("evt_kfactor");
  //cout << "evt_kfactor " << bName << endl;
  if(bName != "") {
    bName.ReplaceAll(".obj", "*");
    t->SetBranchStatus(bName.Data(), 0); 
  }


  //-------------------------------------------------------------

 //Calculate scaling factor and put variables into tree 
  Float_t scale1fb = xsec*kfactor*1000*filt_eff/(Float_t)events;
  cout << "scale1fb: " << scale1fb << endl; 

  TBranch* b1 = t->Branch("evtscale1fb", &scale1fb, "evt_scale1fb/F");
  TBranch* b2 = t->Branch("evtxsecexcl", &xsec, "evt_xsec_excl/F");
  TBranch* b3 = t->Branch("evtxsecincl", &xsec, "evt_xsec_incl/F");
  TBranch* b4 = t->Branch("evtkfactor", &kfactor, "evt_kfactor/F");
  TBranch* b5 = t->Branch("evtnEvts", &events, "evt_nEvts/I");
  TBranch* b6 = t->Branch("evtfilteff", &filt_eff, "evt_filt_eff/F");
   
  t->SetAlias("evt_scale1fb",  "evtscale1fb");
  t->SetAlias("evt_xsec_excl", "evtxsecexcl");
  t->SetAlias("evt_xsec_incl",  "evtxsecincl");
  t->SetAlias("evt_kfactor",   "evtkfactor");
  t->SetAlias("evt_nEvts",     "evtnEvts");
  t->SetAlias("evt_filt_eff",     "evtfilteff");

  Int_t nentries = t->GetEntries();
  for(Int_t i = 0; i < nentries; i++) {
    b1->Fill();
    b2->Fill();
    b3->Fill();
    b4->Fill();
    b5->Fill();
    b6->Fill();
  }
  //-------------------------------------------------------------


 TFile *out = TFile::Open(outfile.Data(), "RECREATE");
  TTree *clone = t->CloneTree(-1, "fast");
  clone->Write(); 
  out->Close();
  f->Close();
  return;
  
}

