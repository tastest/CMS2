//
#include "ScanChain.h"

//
#include <iostream>
#include <vector>

//
#include "tauify.C"

using namespace tas;

void myBabyMaker::ScanChain( TChain* chain) {

  int nEvents = -1;
 
  // Make a baby ntuple
  MakeBabyNtuple("tautest.root");

  int itest = 0;

  Tauify *t = new Tauify("decay.txt");

  //
  TFile *f = new TFile("test.root","RECREATE");
  
  float pi = TMath::Pi();

  //
  TH1F *h_zmass       = new TH1F("zmass", "zmass", 150, 0, 150);
  TH1F *h_p_cm        = new TH1F("p_cm", "p_cm", 100, 0, 2);
  TH1F *h_costheta_cm = new TH1F("costheta_cm", "costheta_cm", 100, -1, 1);

  TH1F *h_lep_pt      = new TH1F("lep_pt", "lep_pt", 1000, 0, 1000);
  TH1F *h_lep_eta     = new TH1F("lep_eta", "lep_eta", 200, -2.5, 2.5);
  TH1F *h_lep_phi     = new TH1F("lep_phi", "lep_phi", 200, -1*pi, pi);
  TH1F *h_lep_mass      = new TH1F("lep_mass", "lep_mass", 100, 0, .2);

  TH1F *h_taulep_pt   = new TH1F("taulep_pt", "taulep_pt", 1000, 0, 1000);
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

        // Initialize baby ntuple
        InitBabyNtuple();

        // muons
        if( abs(hyp_ll_id().at(iHyp) ) != 13 ) continue;
        if( abs(hyp_lt_id().at(iHyp) ) != 13 ) continue;

        // eta
        if( hyp_ll_p4().at(iHyp).eta() > 2.5 ) continue;
        if( hyp_lt_p4().at(iHyp).eta() > 2.5 ) continue;

        // pt
        if( hyp_ll_p4().at(iHyp).pt() < 20 ) continue;
        if( hyp_lt_p4().at(iHyp).pt() < 20 ) continue;

        // Z mass
        LorentzVector p4_z = hyp_ll_p4().at(iHyp) + hyp_lt_p4().at(iHyp);
        h_zmass->Fill( p4_z.mass() );


/* test it out on the tight lepton */

        // Set
        t->SetLepton( hyp_lt_p4().at(iHyp), 0, 0, 0 );

        // Get
        LorentzVector p4_test  = t->TauP4();
        int   id_test       = t->ParticleId();
        float met_test      = t->TauMET();
        float iso_test      = t->TauIso();
        float ip_test       = t->TauIP();

        if( abs(id_test) != 13 ) continue;

        // don't test this on M <=0 leptons
        if( hyp_lt_p4().at(iHyp).mass() <= 0 ) continue;
      
        // need to come back to understand this
        if( p4_test.pt() == 0 ) continue;

        // CM quantities
        h_p_cm->Fill( t->MomentumCM() );
        h_costheta_cm->Fill( t->CosThetaCM() );

        // lepton in the lab
        h_lep_pt->Fill( hyp_lt_p4().at(iHyp).pt() );
        h_lep_eta->Fill( hyp_lt_p4().at(iHyp).eta() );
        h_lep_phi->Fill( hyp_lt_p4().at(iHyp).phi() );
        h_lep_mass->Fill( hyp_lt_p4().at(iHyp).mass() );

        // lepton made to be a tau
        h_taulep_pt->Fill( p4_test.pt() );
        h_taulep_eta->Fill( p4_test.eta() );
        h_taulep_phi->Fill( p4_test.phi() );

//        // sanity
//        double epsilon = .001;
//        if( 
//            fabs( hyp_lt_p4().at(iHyp).pt() - p4_test.pt() ) > epsilon
//        ){
//          cout << itest << endl;
//          itest++;
//          cout << "pt:\t" << hyp_lt_p4().at(iHyp).pt() << "\t->\t" << p4_test.pt() << endl;
//          cout << "eta:\t" << hyp_lt_p4().at(iHyp).eta() << "\t->\t" << p4_test.eta() << endl;
//          cout << "phi:\t" << hyp_lt_p4().at(iHyp).phi() << "\t->\t" << p4_test.phi() << endl;
//          cout << "mass:\t" << hyp_lt_p4().at(iHyp).mass() << "\t->\t" << p4_test.mass() << endl;
//          cout << p4_test.isTimelike() << "\t" << p4_test.isLightlike() << "\t" << p4_test.isSpacelike() << endl;
//          cout << hyp_lt_p4().at(iHyp).isTimelike() << "\t" << hyp_lt_p4().at(iHyp).isLightlike() << "\t" << hyp_lt_p4().at(iHyp).isSpacelike() << endl;
//        }


        // Fill
        FillBabyNtuple();

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

  //
  CloseBabyNtuple();

  return;
}


void myBabyMaker::InitBabyNtuple () {

  lep_id_       = -999;
  taulep_id_    = -999; 

  p_cm_         = -999;
  costheta_cm_  = -999;  
   
  //p4_lep.clear();
  //p4_taulep.clear();

}

void myBabyMaker::MakeBabyNtuple(const char *babyFilename)
{
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();
    babyFile_ = new TFile(Form("%s", babyFilename), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");

    babyTree_->Branch("lep_id",    &lep_id_,     "lep_id/I");
    babyTree_->Branch("taulep_id", &taulep_id_,  "taulep_id/I");

    babyTree_->Branch("p_cm",         &p_cm_,        "p_cm/F");
    babyTree_->Branch("costheta_cm",  &costheta_cm_, "costheta_cm/F");

    //babyTree_->Branch("p4_lep",     &p4_lep_);
    //babyTree_->Branch("p4_taulep",  &p4_taulep_);

}

//----------------------------------
// Fill the baby
//----------------------------------
void myBabyMaker::FillBabyNtuple()
{
    babyTree_->Fill();
}
//--------------------------------
// Close the baby
//--------------------------------
void myBabyMaker::CloseBabyNtuple()
{   
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();
}

