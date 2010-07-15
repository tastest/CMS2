/* Usage:
   root [0] .L ScanChain.C++
   root [1] TFile *_file0 = TFile::Open("merged_ntuple.root")
   root [2] TChain *chain = new TChain("Events")
   root [3] chain->Add("merged_ntuple.root")

   There are several places where one may create CMS2 cms2
   It can be done here (in a doAll.C script), i.e.:

   root [4] CMS2 cms2 

   It can be done in the source as is done below, or it can be
   ascertained by including CORE/CMS2.cc as is commented out
   below.  They are all the same, and everything will work so
   long as it is created somewhere globally.

   root [5] ScanChain(chain)
*/
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "CMS2.cc"

#include "tauify.C"

using namespace tas;

int ScanChain( TChain* chain, int nEvents = -1, std::string skimFilePrefix="") {

  Tauify *t = new Tauify("decay.txt");

  //
  TFile *f = new TFile("test.root","RECREATE");
  
  float pi = TMath::Pi();

  //
  TH1F *h_zmass       = new TH1F("zmass", "zmass", 150, 0, 150);
  TH1F *h_p_cm        = new TH1F("p_cm", "p_cm", 100, 0, 2);
  TH1F *h_costheta_cm = new TH1F("costheta_cm", "costheta_cm", 100, -1, 1);

  TH1F *h_lep_pt      = new TH1F("lep_pt", "lep_pt", 200, 0, 200);
  TH1F *h_lep_eta     = new TH1F("lep_eta", "lep_eta", 200, -2.5, 2.5);
  TH1F *h_lep_phi     = new TH1F("lep_phi", "lep_phi", 200, -1*pi, pi);

  TH1F *h_taulep_pt   = new TH1F("taulep_pt", "taulep_pt", 500, 0, 5000);
  TH1F *h_taulep_eta  = new TH1F("taulep_eta", "taulep_eta", 200, -2.5, 2.5);
  TH1F *h_taulep_phi  = new TH1F("taulep_phi", "taulep_phi", 200, -1*pi, pi);


  TObjArray *listOfFiles = chain->GetListOfFiles();
  unsigned int nEventsChain=0;
  if(nEvents==-1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;
      // Progress feedback to the user
      if(nEventsTotal%1000 == 0) {
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)) {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
          "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }//if(nEventsTotal%20000 == 0) {

      //
      for(unsigned int iHyp = 0; iHyp < hyp_p4().size(); iHyp++) {

        // muons
        if( abs(hyp_ll_id().at(iHyp) ) != 13 ) continue;
        if( abs(hyp_lt_id().at(iHyp) ) != 13 ) continue;

        // eta
        if( hyp_ll_p4().at(iHyp).eta() > 2.5 ) continue;
        if( hyp_lt_p4().at(iHyp).eta() > 2.5 ) continue;

        // pt
        //if( hyp_ll_p4().at(iHyp).pt() < 10 ) continue;
        //if( hyp_lt_p4().at(iHyp).pt() < 10 ) continue;

        //
        LorentzVector p4_z = hyp_ll_p4().at(iHyp) + hyp_lt_p4().at(iHyp);
        h_zmass->Fill( p4_z.mass() );

        //
        t->SetLepton( hyp_lt_p4().at(iHyp), 0, 0, 0 );

        //LorentzVector id_test       = t->ParticleId();
        //LorentzVector p_test        = t->MomentumCM();
        //LorentzVector costheta_test = t->CosThetaCM();

        LorentzVector p4_test  = t->TauP4();
  
        int   id_test       = t->ParticleId();

        float met_test      = t->TauMET();
        float iso_test      = t->TauIso();
        float ip_test       = t->TauIP();


        h_p_cm->Fill( t->MomentumCM() );
        h_costheta_cm->Fill( t->CosThetaCM() );

        h_lep_pt->Fill( hyp_lt_p4().at(iHyp).pt() );
        h_lep_eta->Fill( hyp_lt_p4().at(iHyp).eta() );
        h_lep_phi->Fill( hyp_lt_p4().at(iHyp).phi() );

        h_taulep_pt->Fill( p4_test.pt() );
        h_taulep_eta->Fill( p4_test.eta() );
        h_taulep_phi->Fill( p4_test.phi() );


      }


    }
    delete tree;
    f.Close();
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  
  //
  f->Write();

  return 0;
}
