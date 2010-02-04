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
  // Removes all non *_CMS2.* branches
  //-------------------------------------------------------------`
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("*_CMS2.*", 1);

 TFile *out = TFile::Open(outfile.Data(), "RECREATE");
  TTree *clone = t->CloneTree(-1, "fast");

  //-------------------------------------------------------------

 //Calculate scaling factor and put variables into tree 
  Float_t scale1fb = xsec*kfactor*1000*filt_eff/(Float_t)events;
  cout << "scale1fb: " << scale1fb << endl; 

  TBranch* b1 = clone->Branch("evtscale1fb", &scale1fb, "evt_scale1fb/F");
  TBranch* b2 = clone->Branch("evtxsecexcl", &xsec, "evt_xsec_excl/F");
  TBranch* b3 = clone->Branch("evtxsecincl", &xsec, "evt_xsec_incl/F");
  TBranch* b4 = clone->Branch("evtkfactor", &kfactor, "evt_kfactor/F");
  TBranch* b5 = clone->Branch("evtnEvts", &events, "evt_nEvts/I");
  TBranch* b6 = clone->Branch("evtfilteff", &filt_eff, "evt_filt_eff/F");
   
  clone->SetAlias("evt_scale1fb",  "evtscale1fb");
  clone->SetAlias("evt_xsec_excl", "evtxsecexcl");
  clone->SetAlias("evt_xsec_incl",  "evtxsecincl");
  clone->SetAlias("evt_kfactor",   "evtkfactor");
  clone->SetAlias("evt_nEvts",     "evtnEvts");
  clone->SetAlias("evt_filt_eff",     "evtfilteff");

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

  clone->Write(); 
  out->Close();
  f->Close();
  return;
  
}

