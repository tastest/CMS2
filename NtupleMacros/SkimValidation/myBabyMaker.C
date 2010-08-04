#include <iostream>

#include "TSystem.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TChainElement.h"
#include "Math/VectorUtil.h"

#include "CMS2.cc"
#include "myBabyMaker.h"
#include "../CORE/electronSelections.cc"
#include "../CORE/electronSelectionsParameters.cc"
#include "../CORE/muonSelections.cc"
#include "../Tools/goodrun.cc"
// for dorky event a la Jake (utilities.h)
#include "../CORE/utilities.h"

using namespace std;
using namespace tas;

//PUT TRIGGERS HERE

myBabyMaker::myBabyMaker() :
  mutrig_("HLT_Mu5"),
  eltrig_("HLT_Ele15_LW_L1R")
{
}

//-----------------------------------
// Looper code starts here
//-----------------------------------
void myBabyMaker::ScanChain( TChain* chain, const char *babyFilename, const char* GoodRunFile, const bool doels, const bool domus) {

  // Set the JSON file
  //set_goodrun_file("./jsonlist_132440_138751.txt");
  set_goodrun_file( GoodRunFile );

  // Make a baby ntuple
  MakeBabyNtuple(babyFilename);

  // duplicate removal
  DorkyEventIdentifier dei;

  int i_permilleOld = 0;
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = 0;

  nEventsChain = chain->GetEntries();
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  map<int, int> m_events;

  while( TChainElement *currentFile = (TChainElement*)fileIter.Next() ) {
    TString filename = currentFile->GetTitle();
    
    TFile f(filename.Data());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    unsigned int nEntries = tree->GetEntries();
    unsigned int nLoop = nEntries;

    for( unsigned int event = 0; event < nLoop; event++) {	// Event Loop
      cms2.GetEntry(event);

      // looper progress
      ++nEventsTotal;
      int i_permille = (int)floor(1000 * float(nEventsTotal) / float(nEventsChain));
      if (i_permille != i_permilleOld) {
        printf("  \015\033[32m ---> \033[1m\033[31m%4.1f%%" "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
        fflush(stdout);
        i_permilleOld = i_permille;
      }

      // Good  Runs
      if( !goodrun(evt_run(), evt_lumiBlock()) )
		continue;

      // remove duplicate events
      if ( dei.is_duplicate(DorkyEvent()) ) {
        //        std::cout<<".dr";
        continue;
      }


	  // Initialize baby ntuple--clear every event
	  InitBabyNtuple();

	  vector<int> goodels;
	  vector<int> goodmus;

    // Select events passing trigger w/ leptons
      if( domus && passHLTTrigger(mutrig_) ) {
		goodmus = doMuons();
	  }	  

       if( doels && passHLTTrigger(eltrig_) ) {
 		goodels = doElectrons();
 	  }


	  // Time to fill the baby
	  //if( goodels.size() > 0 || goodmus.size() > 0 ) {
	  if( nels_ > 0 || nmus_ > 0 ) {
		//cout << "I'm a good event" << endl;

		// event info
		run_    = evt_run();
		event_  = evt_event();
		lumi_   = evt_lumiBlock();

		// Jets
		jets_     = CleanJets( jets_p4()   , goodels, goodmus, true ); //true for cor for calo only
		pfjets_   = CleanJets( pfjets_p4() , goodels, goodmus );
		trkjets_  = CleanJets( trkjets_p4(), goodels, goodmus );

		doMet();
		FillBabyNtuple();
	  }

    } // closes loop over events
  } // closes loop over files
  cout << endl << endl;

  CloseBabyNtuple();
  return;
} // closes myLooper function  

//------------------------------------------
// Initialize baby ntuple variables
//------------------------------------------
void myBabyMaker::InitBabyNtuple () {
  //ints
  nels_  = 0;
  nmus_   = 0;
  nJets_ = 0;

  run_    = -999;
  event_  = -999;
  lumi_   = -999;

  // floats--mus
  musd0corr_.clear();
  muscharge_.clear();
  musp4_.clear();
  musreliso_.clear();
  mustrkiso_.clear();
  musecliso_.clear();
  mushcliso_.clear();
  musId_.clear();

  //els
  elsd0corr_.clear();
  elscharge_.clear();
  elsp4_.clear();
  elsreliso_.clear();
  elstrkiso_.clear();
  elsecliso_.clear();
  elshcliso_.clear();
  elsId_.clear();

  // Jets
  jets_   .clear();
  pfjets_ .clear();
  trkjets_.clear();
  //hypjets_.clear();

  //met
  clmet_ 		= -999.;
  tcmet_ 		= -999.;
  pfmet_ 		= -999.;
  clmetphi_ 	= -999.;
  tcmetphi_ 	= -999.;
  pfmetphi_ 	= -999.;
}

//-------------------------------------
// Book the baby ntuple
//-------------------------------------
void myBabyMaker::MakeBabyNtuple(const char *babyFilename)
{
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();
    babyFile_ = new TFile(Form("%s", babyFilename), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");

    //babyTree_->Branch("pt",   &pt_,   "pt/F"); //leftover
	//els
    babyTree_->Branch("nel"      , &nels_		);
    babyTree_->Branch("elsid"    , &elsId_		);	
    babyTree_->Branch("elsd0corr", &elsd0corr_	);
    babyTree_->Branch("elscharge", &elscharge_	);
    babyTree_->Branch("elsp4"    , &elsp4_ 		);
	babyTree_->Branch("elsiso"   , &elsreliso_	); //this is rel, but not max(20,pt)
    babyTree_->Branch("elstrkiso", &elstrkiso_	);
    babyTree_->Branch("elsecliso", &elsecliso_	);
    babyTree_->Branch("elshcliso", &elshcliso_	);

    babyTree_->Branch("nmu",       &nmus_		);
    babyTree_->Branch("musid",     &musId_		);
    babyTree_->Branch("musd0corr", &musd0corr_	);
    babyTree_->Branch("muscharge", &muscharge_	);
    babyTree_->Branch("musp4",     &musp4_		);
    babyTree_->Branch("musiso",    &musreliso_	);
    babyTree_->Branch("mustrkiso", &mustrkiso_	);
    babyTree_->Branch("musecliso", &musecliso_	);
    babyTree_->Branch("mushcliso", &mushcliso_	);

    babyTree_->Branch("clmet"   , &clmet_   );
    babyTree_->Branch("tcmet"   , &tcmet_   );
    babyTree_->Branch("pfmet"   , &pfmet_   );
    babyTree_->Branch("clmetphi", &clmetphi_);
    babyTree_->Branch("tcmetphi", &tcmetphi_);
    babyTree_->Branch("pfmetphi", &pfmetphi_);

    babyTree_->Branch("run"  , &run_  , "run/i"  );
    babyTree_->Branch("event", &event_, "event/i");
    babyTree_->Branch("lumi" , &lumi_ , "lumi/i" );
	//jets
    babyTree_->Branch("pfjets",     &pfjets_);
    babyTree_->Branch("trkjets",    &trkjets_);
    babyTree_->Branch("jets",       &jets_);

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

vector<int> myBabyMaker::doElectrons() {

  //take the ttbarV1 but turn off iso
  const cuts_t elselnoiso = electronSelection_ttbarV1 & ~(1ll<<ELEISO_REL015);
  const cuts_t elsel      = electronSelection_ttbarV1;  
  vector<int> goodels;
  //bool havegoodel = false;
  for( unsigned int iel = 0; iel < els_p4().size(); iel++) {
	if( !pass_electronSelection(iel, elselnoiso, false) )
	  continue;

	//write data members
	nels_++;
	elsp4_.push_back(els_p4()[iel]);
	elscharge_.push_back(els_charge()[iel]);
	elsd0corr_.push_back(els_d0corr()[iel]);
	elsreliso_.push_back( (els_ecalIso()[iel] + els_hcalIso()[iel] + els_tkIso()[iel]) /
						els_p4()[iel].pt() );
	elstrkiso_.push_back(els_tkIso()[iel] );
	elsecliso_.push_back(els_ecalIso()[iel]);
	elshcliso_.push_back(els_hcalIso()[iel]);

	bool elId = pass_electronSelection(iel, elsel, false);
	elsId_.push_back( elId );
	if( elId )
	  //havegoodel = true;
	  goodels.push_back(iel);
  }// closes loop over electrons

  //return havegoodel;
  return goodels;
}


vector<int> myBabyMaker::doMuons() {
  vector<int> goodmus;

  // Loop over muons
  int nmus = 0;
  for( unsigned int iMu = 0; iMu < mus_p4().size(); iMu++) {

	if( mus_p4().at(iMu).pt() < 10 ) continue;
	if( fabs( mus_p4().at(iMu).eta() ) > 2.5 ) continue; 
	if( !muonIdNotIsolated( iMu, NominalTTbar ) ) continue;

	nmus++; //num id'd, not necessarily isolated

	muscharge_.push_back(mus_charge().at(iMu));

	musp4_.push_back( mus_p4().at(iMu) ); // pt, eta, phi 

	musd0corr_.push_back( mus_d0corr().at(iMu) ); // d0 corrected

	musreliso_.push_back( muonIsoValue(iMu) ); // Relative Isolation
	mustrkiso_.push_back( mus_iso03_sumPt().at(iMu) );
	musecliso_.push_back( mus_iso03_emEt() .at(iMu) );
	mushcliso_.push_back( mus_iso03_hadEt().at(iMu) );

	// Muon Selection
	bool muId = muonId( iMu, NominalTTbar );
	musId_.push_back( muId );
	if( muId ) {
	  goodmus.push_back(iMu); //jet cleaning only
	}
  }

  // number of muons passing full selection
  nmus_ = nmus;

  return goodmus;
}

void myBabyMaker::doMet() {
  clmet_ 		= evt_met();
  tcmet_ 		= evt_tcmet();
  pfmet_ 		= evt_pfmet();
  clmetphi_ 	= evt_metPhi();
  tcmetphi_ 	= evt_tcmetPhi();
  pfmetphi_ 	= evt_pfmetPhi();
}

void myBabyMaker::doJets() {

}


// this is meant to be passed as the third argument, the predicate, of the
// standard library sort algorithm
bool sortByPt(const LorentzVector &vec1, 
			  const LorentzVector &vec2 )
{
    return vec1.pt() > vec2.pt();
}

// Function to apply Pt, Eta, and lep-jet dR selections on jets
vector<LorentzVector> myBabyMaker::CleanJets(
  const vector<LorentzVector> vect_p4_jets,
  const vector<int> elidxs,
  const vector<int> muidxs,
  const bool docor,
  const float jet_pt_threshold,
  const float jet_eta_threshold,
  const float jet_lepton_dR_veto_cone
  )
{
  // check arguments--these should be ok, but you never know
  assert( jet_pt_threshold > 0 && jet_pt_threshold < 55. );
  assert( jet_eta_threshold > 0 && jet_eta_threshold < 3. );
  assert( jet_lepton_dR_veto_cone > 0. && jet_lepton_dR_veto_cone < 2. );
  // Clean Jets
  vector<LorentzVector> cleaned_vect_p4_jets;
  for (unsigned int ijet = 0; ijet < vect_p4_jets.size(); ijet++) {         // loop on jets supplied by user
	float pt = vect_p4_jets.at(ijet).pt();
	if( docor )
	  pt *= jets_cor().at(ijet);
    if( pt < jet_pt_threshold ) // apply pt threshold
	  continue;
    if( fabs( vect_p4_jets.at(ijet).eta() ) > jet_eta_threshold ) // apply eta threshold
	  continue;
	
    // apply electron-jet dR vetos
    for (unsigned int iel = 0; iel < elidxs.size(); iel++) {
	  const int idx = elidxs.at(iel);
      if( ROOT::Math::VectorUtil::DeltaR( vect_p4_jets.at(ijet), els_p4().at(idx) ) < jet_lepton_dR_veto_cone )
		continue;  
    }
    // apply muon-jet dR vetos
    for (unsigned int imu = 0; imu < muidxs.size(); imu++) {
	  const int idx = muidxs.at(imu);
      if( ROOT::Math::VectorUtil::DeltaR( vect_p4_jets.at(ijet), mus_p4().at(idx) ) < jet_lepton_dR_veto_cone )
		continue;  
    }
	LorentzVector jet = vect_p4_jets.at(ijet);
	if( docor )
	  jet *= jets_cor().at(ijet);
    cleaned_vect_p4_jets.push_back( jet ); // keep jets within thresholds
  }
  // sort jets by pt
  sort( cleaned_vect_p4_jets.begin(), cleaned_vect_p4_jets.end(), sortByPt );  // sort jets by pt
  return cleaned_vect_p4_jets;
}
