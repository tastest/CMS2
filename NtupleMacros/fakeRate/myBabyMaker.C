#include <iostream>
#include <set>

// ROOT includes
#include "TSystem.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TChainElement.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Math/VectorUtil.h"

// TAS includes
#include "./CMS2.cc"
#include "../CORE/trackSelections.cc"
#include "../CORE/eventSelections.cc"

#include "../CORE/muonSelections.cc"
#include "../CORE/electronSelections.cc"
#include "../CORE/electronSelectionsParameters.cc"
#include "../CORE/metSelections.cc"


//#include "../CORE/fakerates.cc"
#include "../CORE/triggerUtils.cc"
#include "../Tools/goodrun.cc"
#include "./myBabyMaker.h"
using namespace std;
using namespace tas;

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

 
// transverse mass
float Mt( LorentzVector p4, float met, float met_phi ){
  return sqrt( 2*met*( p4.pt() - ( p4.Px()*cos(met_phi) + p4.Py()*sin(met_phi) ) ) );
}



//-----------------------------------
// Looper code starts here
// eormu=-1 do both e and mu
//      =11 do electrons
//      =13 do muons
//-----------------------------------
void myBabyMaker::ScanChain( TChain* chain, const char *babyFilename, bool isData, int eormu) {

  already_seen.clear();

  // Make a baby ntuple
  MakeBabyNtuple(babyFilename);

  // Set the JSON file
  if(isData){
    //set_goodrun_file("./jsonlist_133446_140387_254.4nb.txt");
    set_goodrun_file("jsonlist_132440_139239.txt");
  }

  // The deltaR requirement between objects and jets to remove the jet trigger dependence
  float deltaRCut   = 1.0;
  float deltaPhiCut = 2.5;

  //--------------------------
  // File and Event Loop
  //---------------------------
  int i_permilleOld = 0;
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = 0;
  int nEvents = -1;
  if (nEvents==-1){
    nEventsChain = chain->GetEntries();
  } else {
    nEventsChain = nEvents;
  }
  nEventsChain = chain->GetEntries();
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  map<int,int> m_events;
  while(TChainElement *currentFile = (TChainElement*)fileIter.Next() ) {
    TString filename = currentFile->GetTitle();
    
    TFile f(filename.Data());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    unsigned int nEntries = tree->GetEntries();
    unsigned int nLoop = nEntries;
    unsigned int z;
    for( z = 0; z < nLoop; z++) {	// Event Loop
      cms2.GetEntry(z);

      if(isData){
        // Good  Runs
        if(!goodrun( evt_run(), evt_lumiBlock() )) continue;

        // check for duplicated
        DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
        if (is_duplicate(id) ){ 
          cout << "\t! ERROR: found duplicate." << endl;
          continue;
        }
      }

      // looper progress
      ++nEventsTotal;
      int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
      if (i_permille != i_permilleOld) {
        printf("  \015\033[32m ---> \033[1m\033[31m%4.1f%%" "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
        fflush(stdout);
        i_permilleOld = i_permille;
      }
      
      // Event cleaning (careful, it requires technical bits)
      // if (!cleaning_standard(true)) continue;
      if (!cleaning_BPTX(isData))   continue;
      if (!cleaning_beamHalo())   continue;
      if (!cleaning_goodVertex()) continue;
      if (!cleaning_goodTracks()) continue;

      // Loop over jets and see what is btagged
      // Medium operating point from https://twiki.cern.ch/twiki/bin/view/CMS/BTagPerformanceOP
      int this_nbjet = 0;
      vector<unsigned int> bindex;
      for (unsigned int iJet = 0; iJet < jets_p4().size(); iJet++) {
	      if (jets_p4().at(iJet).pt() < 15.) continue;
	      if (jets_simpleSecondaryVertexHighEffBJetTag().at(iJet) < 1.74) continue;
	      this_nbjet++;
	      bindex.push_back(iJet);
      }
      
/* Electrons */

      if (eormu == -1 || eormu==11) {
	    for (unsigned int iEl = 0 ; iEl < els_p4().size(); iEl++) {

	      // ECAL spike cleaning
	      //float r19 = cms2.els_eMax()[iEl]/cms2.els_e5x5()[iEl];
	      //if (r19 >= 0.95) continue;

	      // Apply a pt cut (Changed it from 5 GeV to 10 GeV...Claudio 10 July 2010)
	      if ( els_p4().at(iEl).pt() < 10.) continue;


	      // Initialize baby ntuple
	      InitBabyNtuple();

	      // Add spike veto
	      num_ = pass_electronSelection( iEl, electronSelection_ttbarV1 ) && (!isSpikeElectron(iEl));
	      numv1_ = pass_electronSelection( iEl, electronSelection_ttbarV1 );
	      v1_  = pass_electronSelection( iEl, electronSelectionFO_el_ttbarV1_v1 );
	      v2_  = pass_electronSelection( iEl, electronSelectionFO_el_ttbarV1_v2 );
	      v3_  = pass_electronSelection( iEl, electronSelectionFO_el_ttbarV1_v3 );

	      numSS_ = pass_electronSelection(iEl, electronSelection_ss);
	      v1SS_  = pass_electronSelection(iEl, electronSelectionFO_ssVBTF80_v1);
	      v2SS_  = pass_electronSelection(iEl, electronSelectionFO_ssVBTF80_v2);
	      v3SS_  = pass_electronSelection(iEl, electronSelectionFO_ssVBTF80_v3);

        numAug9_ = pass_electronSelection( iEl, electronSelection_ttbarV1, true, true ) && (!isSpikeElectron(iEl));
	      v1Aug9_  = v1_;
	      v2Aug9_  = v2_;
	      v3Aug9_  = pass_electronSelection( iEl, electronSelectionFO_el_ttbarV1_v3, true, true );


	      // Sanity
	      if (num_ && (!v1_)) cout << "bad v1" << endl;
	      if (num_ && (!v2_)) cout << "bad v2" << endl;
	      if (num_ && (!v3_)) cout << "bad v3" << endl;

	      if (numSS_ && (!v1SS_)) cout << "bad v1" << endl;
	      if (numSS_ && (!v2SS_)) cout << "bad v2" << endl;
	      if (numSS_ && (!v3SS_)) cout << "bad v3" << endl;

	      if (numAug9_ && (!v1Aug9_)) cout << "bad v1" << endl;
	      if (numAug9_ && (!v2Aug9_)) cout << "bad v2" << endl;
	      if (numAug9_ && (!v3Aug9_)) cout << "bad v3" << endl;


	      // If there is no v1/v2/v3 lepton quit
	      if ( (!v1_) && (!v2_) && (!v3_) && (!v1SS_) && (!v2SS_) && (!v3SS_) ) continue;
	      
	      // If it is above 20 GeV see if we can make a 
	      // Z with another pt>20 FO.  Will use the v1 FO since 
	      // these are the loosest
	      bool isaZ = false;
	      if (els_p4().at(iEl).pt() > 20.) {
		for (unsigned int jEl = 0 ; jEl < els_p4().size(); jEl++) {
		  if (iEl == jEl)                             continue;
		  if (els_p4().at(jEl).pt() < 20.)            continue;
		  if ( ! pass_electronSelection( jEl, electronSelectionFO_el_ttbarV1_v1 ) ) continue;
		  if ( ! v1_ ) continue;
		  LorentzVector w = els_p4().at(iEl) + els_p4().at(jEl);
		  if (abs(w.mass()-91.) > 20.) continue;
		  isaZ = true;
		}
	      }
	      if (isaZ) continue;
		
	      // Load the electron and event quantities
	      run_   = evt_run();
	      ls_    = evt_lumiBlock();
	      evt_   = evt_event();
	      pt_    = els_p4().at(iEl).pt();
	      eta_   = els_p4().at(iEl).eta();
	      phi_   = els_p4().at(iEl).phi();
	      scet_  = els_eSC()[iEl] / cosh( els_etaSC()[iEl] );
	      id_    = 11*els_charge().at(iEl);
	      tcmet_ = evt_tcmet();
	      tcmetphi_ = evt_tcmetPhi();
		
	      // W transverse mass
	      mt_ = Mt( els_p4().at(iEl), tcmet_, tcmetphi_ );
	      
	      // The btag information
	      nbjet_ = this_nbjet;
	      dRbNear_ = 99.;
	      dRbFar_  = -99.;
	      for (int ii=0; ii<nbjet_; ii++) {
		unsigned int iJet = bindex[ii];
		float dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), jets_p4().at(iJet));
		if (dr < dRbNear_) dRbNear_ = dr;
		if (dr > dRbFar_)   dRbFar_  = dr;
	      }
			    
		
	      // Our jet trigger flags
	      hlt15u_ = min(2,nHLTObjects("HLT_Jet15U")); 
	      hlt30u_ = min(2,nHLTObjects("HLT_Jet30U")); 
	      hlt50u_ = min(2,nHLTObjects("HLT_Jet50U")); 
	      l16u_   = min(2,nHLTObjects("HLT_L1Jet6U"));
	      l110u_  = min(2,nHLTObjects("HLT_L1Jet10U"));

	      // If only one jet triggered, see if it is far enough away 
	      if (hlt15u_ == 1) {
		LorentzVector p4j = p4HLTObject("HLT_Jet15U",0);
		double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4j);
		if (dr > deltaRCut) hlt15u_ = 2;
	      }
	      if (hlt30u_ == 1) {
		LorentzVector p4j = p4HLTObject("HLT_Jet30U",0);
		double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4j);
		if (dr > deltaRCut) hlt30u_ = 2;
	      }
	      if (hlt50u_ == 1) {
		LorentzVector p4j = p4HLTObject("HLT_Jet50U",0);
		double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4j);
		if (dr > deltaRCut) hlt50u_ = 2;
	      }
	      if (l16u_ == 1) {
		LorentzVector p4j = p4HLTObject("HLT_L1Jet6U",0);
		double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4j);
		if (dr > deltaRCut) l16u_ = 2;
	      }
	      if (l110u_ == 1) {
		LorentzVector p4j = p4HLTObject("HLT_L1Jet10U",0);
		double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4j);
		if (dr > deltaRCut) l110u_ = 2;
	      }
	      
	      // Now fill the egamma trigger flags  (look at both cleaned and uncleaned photons_
	      int ph10cl = nHLTObjects("HLT_Photon10_Cleaned_L1R");
	      int ph10   = nHLTObjects("HLT_Photon10_L1R");
	      int ph15cl = nHLTObjects("HLT_Photon15_Cleaned_L1R");
	      int ph15   = nHLTObjects("HLT_Photon15_L1R");
	      el10_   = nHLTObjects("HLT_Ele10_LW_L1R");
	      el15_   = nHLTObjects("HLT_Ele15_LW_L1R");
	      eg5_    = nHLTObjects("HLT_L1SingleEG5");
	      eg8_    = nHLTObjects("HLT_L1SingleEG8");
		
	      // For Ele10 the HLT objects are saved for all data we have
	      if (el10_ > 0) {
		bool match = false;
		for (int itrg=0; itrg<el10_; itrg++) {
		  LorentzVector p4tr = p4HLTObject("HLT_Ele10_LW_L1R",itrg);
		  double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4tr);
		  if (dr < drel10_) drel10_ = dr;
		  if (dr < 0.4) match=true;
		}
		if (match) {
		  el10_ = 2;
		} else {
		  el10_ = 1;
		}
	      }
	      
	      // For Ele15 the HLT objects are saved for all data we have
	      if (el15_ > 0) {
		bool match = false;
		for (int itrg=0; itrg<el15_; itrg++) {
		  LorentzVector p4tr = p4HLTObject("HLT_Ele15_LW_L1R",itrg);
		  double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4tr);
		  if (dr < drel15_) drel15_ = dr;
		  if (dr < 0.4) match=true;
		}
		if (match) {
			      el15_ = 2;
		} else {
		  el15_ = 1;
		}
	      }
		
	      // Now for photon10 we look at cleaned and at uncleaned
	      if (ph10 == 0 && ph10cl == 0) ph10_=0;   // trigger failed
	      if (ph10 <  0 || ph10cl <  0) ph10_=-1;  // passed but no object
	      if (ph10cl > 0 || ph10 > 0) {
		bool match = false;
		for (int itrg=0; itrg<ph10; itrg++) {
		  LorentzVector p4tr = p4HLTObject("HLT_Photon10_L1R",itrg);
		  double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4tr);
		  if (dr < drph10_) drph10_ = dr;
		  if (dr < 0.4) match=true;
		}
		for (int itrg=0; itrg<ph10cl; itrg++) {
		  LorentzVector p4tr = p4HLTObject("HLT_Photon10_Cleaned_L1R",itrg);
		  double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4tr);
		  if (dr < drph10_) drph10_ = dr;
		  if (dr < 0.4) match=true;
		}
		if (match) {
		  ph10_ = 2;
		} else {
		  ph10_ = 1;
		}
	      }
	      
	      // Now for photon15 we look at cleaned and at uncleaned
	      if (ph15 == 0 && ph15cl == 0) ph15_=0;   // trigger failed
	      if (ph15 <  0 || ph15cl <  0) ph15_=-1;  // passed but no object
	      if (ph15cl > 0 || ph15 > 0) {
		bool match = false;
		for (int itrg=0; itrg<ph15; itrg++) {
		  LorentzVector p4tr = p4HLTObject("HLT_Photon15_L1R",itrg);
		  double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4tr);
		  if (dr < drph15_) drph15_ = dr;
		  if (dr < 0.4) match=true;
		}
		for (int itrg=0; itrg<ph15cl; itrg++) {
		  LorentzVector p4tr = p4HLTObject("HLT_Photon15_Cleaned_L1R",itrg);
		  double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4tr);
		  if (dr < drph15_) drph15_ = dr;
		  if (dr < 0.4) match=true;
		}
		if (match) {
		  ph15_ = 2;
		} else {
		  ph15_ = 1;
		}
	      }


	      // For EG5 and EG8 the HLT object is missing, so we'll try with the L1 info
	      // There appear to be two L1 EM objects: one with "iso" and one without.
	      // From some event dumps it looks like the iso one is the one we want
	      // Also: the Photon10 HLT object is sometimes missing, but I know that EG5 is its prerequisite
	      // so if it is missing we will match to EG5  (TAKE THIS OUT.... July 5th 2010)
	      if ( (eg5_ == -1 || eg8_ == -1 || ph10_ == -1) && l1_emiso_p4().size()>0) {
			    
		for (unsigned int ig = 0 ; ig < l1_emiso_p4().size(); ig++) {
		  
		  //if (ph10_ == -1 && l1_emiso_p4().at(ig).pt() > 5.) {
		  //	double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl),l1_emiso_p4().at(ig)); 
		  //	if (dr < drph10_) drph10_ = dr;
		  //}
		  if (eg5_ == -1 && l1_emiso_p4().at(ig).pt() > 5.) {
		    double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl),l1_emiso_p4().at(ig)); 
		    if (dr < dreg5_) dreg5_ = dr;
		  }
		  if (eg8_ == -1 && l1_emiso_p4().at(ig).pt() > 8.) {
		    double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl),l1_emiso_p4().at(ig)); 
		    if (dr < dreg8_) dreg8_ = dr;
		  }
		  
		} // closes loop over L1 EG objects
		
		//if (ph10_ == -1) {
		//  if (drph10_ < 100.) ph10_ = 1; // objects found, but no match
		//  if (drph10_ < 0.4)  ph10_ = 2; // objects found, and matched
		//}
		if (eg5_ == -1) {
		  if (dreg5_ < 100.) eg5_ = 1; // objects found, but no match
		  if (dreg5_ < 0.4)  eg5_ = 2; // objects found, and matched
		}
		if (eg8_ == -1) {
		  if (dreg8_ < 100.) eg8_ = 1; // objects found, but no match
		  if (dreg8_ < 0.4)  eg8_ = 2; // objects found, and matched
		}
		
	      } // closes if-block of EG5, EG8, PH10 objects missing
 
	      // Find the highest Pt jet separated by at least dRcut from this lepton and fill the jet Pt
	      ptj1_       = -999.0;
	      ptj1_b2b_   = -999.0;
	      dphij1_b2b_ = -999.0;
	      nj1_        = 0;
	      for (unsigned int iJet = 0; iJet < jets_p4().size(); iJet++) {
		      double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), jets_p4().at(iJet) );
		      if( dr > deltaRCut && jets_p4().at(iJet).pt() > 10 ) nj1_++;
		      if ( dr > deltaRCut && jets_p4().at(iJet).pt() > ptj1_ ){
		        ptj1_ = jets_p4().at(iJet).pt();
		  
		        // back to back in phi
		        float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iEl), jets_p4().at(iJet) ) );
		        if( dphi > deltaPhiCut && jets_p4().at(iJet).pt() > ptj1_b2b_ ){ 
              ptj1_b2b_   = jets_p4().at(iJet).pt();
		          dphij1_b2b_ = dphi;
            }
		      }
	      }

	      // Time to fill the baby for the electrons
	      FillBabyNtuple();
    
			}// closes loop over electrons
      } // closes if statements about whether we want to fill electrons

/* Muons */

      if (eormu == -1 || eormu==13) {
	for (unsigned int iMu = 0 ; iMu < mus_p4().size(); iMu++) {
	  
	  // Apply a pt cut --- moved the cut to 10 GeV (Claudio, 10 Jul 2010)
	  if ( mus_p4().at(iMu).pt() < 10.) continue;
	  
	  // If it is not a muon FO, quit
	  //if (!(isFakeableMuon(iMu, mu_ttbar))) continue;
	  if ( ! muonId(iMu, muonSelectionFO_mu_ttbar) ) continue;
		
	  // If it is above 20 GeV see if we can make a 
	  // Z with another pt>20 FO.  
	  bool isaZ = false;
	  if (mus_p4().at(iMu).pt() > 20.) {
	    for (unsigned int jMu = 0 ; jMu < mus_p4().size(); jMu++) {
	      if (iMu == jMu)                             continue;
	      if (mus_p4().at(jMu).pt() < 20.)            continue;
	      if ( ! muonId( jMu, muonSelectionFO_mu_ttbar) ) continue;
	      if ( ! muonId( iMu, muonSelectionFO_mu_ttbar) ) continue;
	      LorentzVector w = mus_p4().at(iMu) + mus_p4().at(jMu);
	      if (abs(w.mass()-91.) > 20.) continue;
	      isaZ = true;
	    }
	  }
	  if (isaZ) continue;
	  
	  // Initialize baby ntuple
	  InitBabyNtuple();
	  
	  // Load the muon and event quantities
	  run_  = evt_run();
	  ls_   = evt_lumiBlock();
	  evt_  = evt_event();
	  pt_   = mus_p4().at(iMu).pt();
	  eta_  = mus_p4().at(iMu).eta();
	  phi_  = mus_p4().at(iMu).phi();
	  id_   = 13*mus_charge().at(iMu);
	  num_  = muonId(iMu, NominalTTbarV2);
	  numv1_  = muonId(iMu, NominalTTbar);
	  numSS_  = muonId(iMu, Nominal);
	  tcmet_ = evt_tcmet();
	  tcmetphi_ = evt_tcmetPhi();
	  
    numAug9_ = num_;


	  // Now REALLY fix it (July 14, 2010)
	  if (pt_ > 10.) {
	    if (!wasMetCorrectedForThisMuon(iMu, usingTcMet)) {
	      float metX = tcmet_ * cos(evt_tcmetPhi());
	      float metY = tcmet_ * sin(evt_tcmetPhi());
	      fixMetForThisMuon(iMu, metX, metY, usingTcMet);
	      tcmet_ = sqrt(metX*metX + metY*metY);
	    }
	  }
	  
	  // W transverse mass
	  mt_ = Mt( mus_p4().at(iMu), tcmet_, tcmetphi_ );
	  
	  // The btag information
	  nbjet_ = this_nbjet;
	  dRbNear_ =  99.;
	  dRbFar_  = -99.;
	  for (int ii=0; ii<nbjet_; ii++) {
	    unsigned int iJet = bindex[ii];
	    float dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), jets_p4().at(iJet));
	    if (dr < dRbNear_) dRbNear_ = dr;
	    if (dr > dRbFar_)  dRbFar_  = dr;
	  }
	  
	  // Our jet trigger flags
	  hlt15u_ = min(2,nHLTObjects("HLT_Jet15U")); 
	  hlt30u_ = min(2,nHLTObjects("HLT_Jet30U")); 
	  hlt50u_ = min(2,nHLTObjects("HLT_Jet50U")); 
	  l16u_   = min(2,nHLTObjects("HLT_L1Jet6U"));
	  l110u_  = min(2,nHLTObjects("HLT_L1Jet10U"));
		
	  // If only one jet triggered, see if it is far enough away 
	  if (hlt15u_ == 1) {
	    LorentzVector p4j = p4HLTObject("HLT_Jet15U",0);
	    double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4j);
	    if (dr > deltaRCut) hlt15u_ = 2;
	  }
	  if (hlt30u_ == 1) {
	    LorentzVector p4j = p4HLTObject("HLT_Jet30U",0);
	    double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4j);
	    if (dr > deltaRCut) hlt30u_ = 2;
	  }
	  if (hlt50u_ == 1) {
	    LorentzVector p4j = p4HLTObject("HLT_Jet50U",0);
	    double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4j);
	    if (dr > deltaRCut) hlt50u_ = 2;
	  }
	  if (l16u_ == 1) {
	    LorentzVector p4j = p4HLTObject("HLT_L1Jet6U",0);
	    double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4j);
	    if (dr > deltaRCut) l16u_ = 2;
	  }
	  if (l110u_ == 1) {
	    LorentzVector p4j = p4HLTObject("HLT_L1Jet10U",0);
	    double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4j);
	    if (dr > deltaRCut) l110u_ = 2;
			  }
		
	  // Now fill the muon trigger flags
	  mu3_ = nHLTObjects("HLT_Mu3");
	  mu5_ = nHLTObjects("HLT_Mu5");
	  mu9_ = nHLTObjects("HLT_Mu9");
	  
	  // Explicit match with Mu3 trigger
	  if (mu3_ > 0) {
	    bool match = false;
	    for (int itrg=0; itrg<mu3_; itrg++) {
	      LorentzVector p4tr = p4HLTObject("HLT_Mu3",itrg);
	      double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4tr);
	      if (dr < drmu3_) drmu3_ = dr;
	      if (dr < 0.4) match=true;
	    }
	    if (match) {
	      mu3_ = 2;
	    } else {
	      mu3_ = 1;
	    }
	  }
	  
	  // Explicit match with Mu5 trigger
	  if (mu5_ > 0) {
	    bool match = false;
	    for (int itrg=0; itrg<mu5_; itrg++) {
	      LorentzVector p4tr = p4HLTObject("HLT_Mu5",itrg);
	      double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4tr);
	      if (dr < drmu5_) drmu5_ = dr;
	      if (dr < 0.4) match=true;
	    }
	    if (match) {
	      mu5_ = 2;
	    } else {
	      mu5_ = 1;
	    }
	  }
	  
	  // Explicit match with Mu9 trigger
	  if (mu9_ > 0) {
	    bool match = false;
	    for (int itrg=0; itrg<mu9_; itrg++) {
	      LorentzVector p4tr = p4HLTObject("HLT_Mu5",itrg);
	      double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4tr);
	      if (dr < drmu9_) drmu9_ = dr;
	      if (dr < 0.4) match=true;
	    }
	    if (match) {
	      mu9_ = 2;
	    } else {
	      mu9_ = 1;
	    }
	  }
			  
    // Find the highest Pt jet separated by at least dRcut from this lepton and fill the jet Pt
    ptj1_       = -999.0;
    ptj1_b2b_   = -999.0;
    dphij1_b2b_ = -999.0;
    nj1_        = 0;
    for (unsigned int iJet = 0; iJet < jets_p4().size(); iJet++) {
      double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), jets_p4().at(iJet) );
      if( dr > deltaRCut && jets_p4().at(iJet).pt() > 10 ) nj1_++;
      if ( dr > deltaRCut && jets_p4().at(iJet).pt() > ptj1_ ){
        ptj1_ = jets_p4().at(iJet).pt();

        // back to back in phi
        float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iMu), jets_p4().at(iJet) ) );
        if( dphi > deltaPhiCut && jets_p4().at(iJet).pt() > ptj1_b2b_ ){        
          ptj1_b2b_   = jets_p4().at(iJet).pt();
          dphij1_b2b_ = dphi;
        }
      }
    }

	  
  // Time to fill the baby for the muons
  FillBabyNtuple();

    }// closes loop over muons
  } // closes if statements about whether we want to fill muons

}// closes loop over events
}  // closes loop over files
cout << "   " <<endl;
CloseBabyNtuple();
return;
}    // closes myLooper function  

//------------------------------------------
// Initialize baby ntuple variables
//------------------------------------------
void myBabyMaker::InitBabyNtuple () {
  run_ = -1;
  ls_  = -1;
  evt_ = -1;
  id_  = -1;
  pt_  = -999.;
  eta_ = -999.;
  phi_ = -999.;
  scet_ = -999.;
  tcmet_ = -999.;
  tcmetphi_ = -999.;
  hlt15u_ = 0;
  hlt30u_ = 0;
  hlt50u_ = 0;
  l16u_   = 0;
  l110u_  = 0;
  v1_  = false;
  v2_  = false;
  v3_  = false;
  v1SS_  = false;
  v2SS_  = false;
  v3SS_  = false;
  num_ = false;
  numSS_ = false;
  numv1_ = false;
  numAug9_ = false;
  v1Aug9_  = false;
  v2Aug9_  = false;
  v3Aug9_  = false;
  ph10_ = 0;
  ph15_ = 0;
  el10_ = 0;
  el15_ = 0;
  eg5_  = 0;
  eg8_  = 0;
  mu5_  = 0;
  mu9_  = 0;
  mu3_  = 0;
  drph10_ = 99.;
  drph15_ = 99.;
  drel10_ = 99.;
  drel15_ = 99.;
  dreg5_  = 99.;
  dreg8_  = 99.;
  drmu5_  = 99.;
  drmu9_  = 99.;
  drmu3_  = 99.;
  nbjet_  = 0;
  dRbNear_ = 99.;
  dRbFar_ = -99.;
  ptj1_   = 0.;
  nj1_    = 0;
  ptj1_b2b_ = -999.;
  dphij1_b2b_ = -999.;
  mt_ = -999;
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

    babyTree_->Branch("run",          &run_,         "run/I"         );
    babyTree_->Branch("ls",           &ls_,          "ls/I"          );
    babyTree_->Branch("evt",          &evt_,         "evt/I"         );

    babyTree_->Branch("pt",           &pt_,          "pt/F"         );
    babyTree_->Branch("eta",          &eta_,         "eta/F"         );
    babyTree_->Branch("phi",          &phi_,         "phi/F"         );
    babyTree_->Branch("scet",          &scet_,         "scet/F"         );
    babyTree_->Branch("tcmet",          &tcmet_,         "tcmet/F"         );
    babyTree_->Branch("tcmetphi",          &tcmetphi_,         "tcmetphi/F"         );
    babyTree_->Branch("id",          &id_,         "id/I"         );

    babyTree_->Branch("hlt15u",       &hlt15u_,       "hlt15u/I"      );
    babyTree_->Branch("hlt30u",       &hlt30u_,       "hlt30u/I"      );
    babyTree_->Branch("hlt50u",       &hlt50u_,       "hlt50u/I"      );
    babyTree_->Branch("l16u",         &l16u_,         "l16uu/I"      );
    babyTree_->Branch("l110",         &l110u_,        "l110u/I"      );

    babyTree_->Branch("v1",         &v1_,        "v1/O"      );
    babyTree_->Branch("v2",         &v2_,        "v2/O"      );
    babyTree_->Branch("v3",         &v3_,        "v3/O"      );
    babyTree_->Branch("num",         &num_,        "num/O"      );
    babyTree_->Branch("numv1",         &numv1_,        "numv1/O"      );

    babyTree_->Branch("v1SS",         &v1SS_,        "v1SS/O"      );
    babyTree_->Branch("v2SS",         &v2SS_,        "v2SS/O"      );
    babyTree_->Branch("v3SS",         &v3SS_,        "v3SS/O"      );
    babyTree_->Branch("numSS",         &numSS_,        "numSS/O"      );

    babyTree_->Branch("v1Aug9",         &v1Aug9_,        "v1Aug9/O"      );
    babyTree_->Branch("v2Aug9",         &v2Aug9_,        "v2Aug9/O"      );
    babyTree_->Branch("v3Aug9",         &v3Aug9_,        "v3Aug9/O"      );
    babyTree_->Branch("numAug9",         &numAug9_,        "numAug9/O"      );

    babyTree_->Branch("ph10",       &ph10_,       "ph10/I"      );
    babyTree_->Branch("ph15",       &ph15_,       "ph15/I"      );
    babyTree_->Branch("el10",         &el10_,         "el10/I"      );
    babyTree_->Branch("el15",         &el15_,         "el15/I"      );
    babyTree_->Branch("eg5",         &eg5_,        "eg5/I"      );
    babyTree_->Branch("eg8",         &eg8_,        "eg8/I"      );

    babyTree_->Branch("drph10",       &drph10_,       "drph10/F"      );
    babyTree_->Branch("drph15",       &drph15_,       "drph15/F"      );
    babyTree_->Branch("drel10",         &drel10_,         "drel10/F"      );
    babyTree_->Branch("drel15",         &drel15_,         "drel15/F"      );
    babyTree_->Branch("dreg5",         &dreg5_,        "dreg5/F"      );
    babyTree_->Branch("dreg8",         &dreg8_,        "dreg8/F"      );

    babyTree_->Branch("mu9",       &mu9_,       "mu9/I"      );
    babyTree_->Branch("mu5",       &mu5_,       "mu5/I"      );
    babyTree_->Branch("mu3",       &mu3_,       "mu3/I"      );

    babyTree_->Branch("drmu9",       &drmu9_,       "drmu9/F"      );
    babyTree_->Branch("drmu5",       &drmu5_,       "drmu5/F"      );
    babyTree_->Branch("drmu3",       &drmu3_,       "drmu3/F"      );

    babyTree_->Branch("nbjet",       &nbjet_,       "nbjet/I"      );
    babyTree_->Branch("dRNear",       &dRbNear_,       "dRbNear/F"      );
    babyTree_->Branch("dRFar",       &dRbFar_,       "dRbFar/F"      );

    babyTree_->Branch("ptj1",       &ptj1_,       "ptj1/F"      );
    babyTree_->Branch("nj1",       &nj1_,       "nj1/I"      );
    babyTree_->Branch("ptj1_b2b",       &ptj1_b2b_,       "ptj1_b2b/F"      );
    babyTree_->Branch("dphij1_b2b",       &dphij1_b2b_,       "dphij1_b2b/F"      );

    babyTree_->Branch("mt",          &mt_,         "mt/F"         );

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




