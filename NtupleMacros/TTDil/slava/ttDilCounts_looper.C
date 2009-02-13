/* Usage:
   root[0] .L ttDilCounts_looper.C++
   root [1] TFile *_file0 = TFile::Open("ntuple_file.root")
   root [2] TChain *chain = new TChain("Events")
   root [3] chain->Add("ntuple_file.root")
   root [4] ttDilCounts_looper a 
   root [5] a.ScanChain(chain) // will give the same results
*/
#include <iostream>
#include <vector>
#include <map>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "TChain.h"
#include "TFile.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TROOT.h"

//from slava
#include "TH1F.h"
#include "TH2F.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include <algorithm>
#include "TRandom2.h"
#include <fstream>

#include "ttDilCounts_looper.h"
#include "selections.C"


typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >  VofP4;

// this is Jake's magic to sort jets by Pt
Bool_t comparePt(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lv1,
                 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lv2) {
  return lv1.pt() > lv2.pt();
}





int ttDilCounts_looper::ScanChain ( TChain* chain, char * prefix, float kFactor, int prescale, bool oldjets, unsigned int cutsMask){
  
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  using namespace std;
  
  //ostream& out = cout;
  ofstream outf;
  outf.open(Form("eventsCMS2_%d.txt",cutsMask));
  outf << "Starting " << prefix << " bitmask: " << cutsMask << endl;
  std::cout << "Starting " << prefix << " bitmask: " << cutsMask << std::endl;
  
  
  //book Histograms
  bookHistos(prefix);
  
  bool idcuts = false;
  bool isolationcuts = false;
  bool dilepMassVetoCut = false;
  bool METcut = false;
  bool nJets2 = false; //true if we want >=2 jets
  bool applyMuTag = false;
  bool METveto = false;
  bool applyMuTag5 = false;
  int  isoLooseMode = 0;
  bool looseDilSelectionTTDil08 = false;
  bool fillMultipleHypsOnly = false;
  bool applyZWindow = false;
  bool osSelection = false;
  bool fillMaxWeightDilOnly = false;
  bool leptonIsolationDilSelectionTTDil08 = false;
  bool looseDilSelectionNoIsoTTDil08 = false;
  bool lepton20Eta2p4DilSelection = false;
  bool metBaselineSelectionTTDil08 = false;
  bool dilepMassVetoCutTTDil08 = false;

  bool applyTriggersMu9orLisoE15 = false;
  bool applyTriggersTTDil08JanTrial = false;

  if( cutsMask & 1) {
    idcuts = true;
    cout << "Id cuts enabled" << endl;
  } else{
    idcuts = false;
    cout << "Id cuts disabled" << endl;
  }

  if((cutsMask>>1) & 1 ) {
    isolationcuts = true;
    cout << "Isolation cuts enabled" << endl;
  } else {
    isolationcuts = false;
    cout << "Isolation cuts disabled" << endl;
  }

  if( (cutsMask>>2) & 1) {
    dilepMassVetoCut = true;
    cout << "DiLeptonMassVetoCut enabled" << endl;
  } else {
    dilepMassVetoCut = false;
    cout << "DiLeptonMassVetoCut disabled" << endl;
  }

  if( (cutsMask>>3) & 1) {
    METcut = true;
    cout << "METCut enabled" << endl;
  } else {
    METcut = false;
    cout << "METCut disabled" << endl;
  }

  if( (cutsMask>>4) & 1 ) {
    nJets2 = true;
    cout << "NJets>=2 cut enabled" << endl;
  } else {
    nJets2 = false;
    cout << "NJets>=2 cut disabled" << endl;
  }
  if( (cutsMask>>5) & 1 ) {
    applyMuTag = true;
    cout << "Extra Muon tag cut enabled" << endl;
  } else {
    applyMuTag = false;
    cout << "Extra muon tag cut disabled" << endl;
  }
  

  if ( (cutsMask>>6) & 1 ){
    METveto = true;
    cout << "MET veto is enabled" << endl;
  } else {
    METveto = false;
    cout << "MET veto is disabled" << endl;
  }

  if( (cutsMask>>7) & 1) {
    applyMuTag5 = true;
    cout << "Extra Muon 5GeV tag cut enabled" << endl;
  } else {
    applyMuTag5 = false;
    cout << "Extra muon 5GeV tag cut disabled" << endl;
  }

  
  isoLooseMode = ((cutsMask>>8) & 3);
  if (isoLooseMode == 0){
    cout << "Require both leptons to be isolated" << endl;
  } else if( isoLooseMode == 1) {
    cout << "Require one-only hyp lepton (electron in emu) with iso (0.6, 0.92)" << endl;
  } else if (isoLooseMode == 2 ) {
    cout << "Require one-only hyp lepton (muon in emu) with iso (0.6, 0.92)" << endl;
  } else if (isoLooseMode == 3){
    cout << "Require two hyp leptons  with iso (0.6, 0.92)" << endl;
  }

  looseDilSelectionTTDil08 = ((cutsMask>>10)&1);
  if (looseDilSelectionTTDil08 ){
    std::cout<< "Require loose dil selection for TTbar-->dilepton ana 2008/09"<<std::endl;
  }

  fillMultipleHypsOnly = ((cutsMask >> 11)&1);
  if (fillMultipleHypsOnly){
    std::cout<< "Fill only multiple hypotheses"<<std::endl;
  }

  applyZWindow = ((cutsMask >> 12) & 1);
  if (applyZWindow){
    std::cout<<"Events from Z-window only"<<std::endl;
  }
  
  osSelection = ((cutsMask >> 13 ) & 1);
  if (osSelection){
    std::cout<<"Require OS"<<std::endl;
  }

  fillMaxWeightDilOnly = ((cutsMask >> 14) & 1);
  if (fillMaxWeightDilOnly){
    std::cout<<"Fill only the dilepton with the max dilepton weight"<<std::endl;
  }

  leptonIsolationDilSelectionTTDil08 = ((cutsMask>> 15) & 1);
  if (leptonIsolationDilSelectionTTDil08){
    std::cout<<" Apply isolation cuts on leptons for TTbar-->dilepton ana 2008/09"<<std::endl;
  }

  looseDilSelectionNoIsoTTDil08 = ((cutsMask>>16)&1);
  if (looseDilSelectionNoIsoTTDil08 ){
    std::cout<< "Require loose dil selection for TTbar-->dilepton ana 2008/09; drop loose iso cuts"<<std::endl;
  }

  lepton20Eta2p4DilSelection = ((cutsMask>>17)&1);
  if (lepton20Eta2p4DilSelection){
    std::cout<<"Two leptons pt>20 and |eta|<2.4 are selected -- bare minimum"<<std::endl;
  }

  metBaselineSelectionTTDil08 = ((cutsMask>>18)&1);
  if (metBaselineSelectionTTDil08){
    std::cout<<"Apply TTDil08 baseline MET selection: use corrected pat-met emu met >20, mm,em met>30"<<std::endl;
  }

  dilepMassVetoCutTTDil08 = ((cutsMask>>19)&1);
  if (dilepMassVetoCutTTDil08){
    std::cout<<"Apply Z mass veto on same flavor dils, use TTDil08 selections"<<std::endl;
  }

  applyTriggersMu9orLisoE15 = ((cutsMask>>20)&1);
  if (applyTriggersMu9orLisoE15){
    std::cout<<"HLT bits: 47 (HLT_LooseIsoEle15_LW_L1R), 82 (HLT_Mu9): ee -- 47, em -- 47 OR 82, mm -- 82"<<std::endl;
  }
  applyTriggersTTDil08JanTrial = ((cutsMask>>21)&1);
  if (applyTriggersTTDil08JanTrial){
    std::cout<<"HLT bits 45 (IsoEle18_L1R), 54 (DoubleIsoEle12_L1R), 86 (Mu15_L1Mu7), 90 (DoubleMu3), 126 (IsoEle10_Mu10_L1R):\n"
	     <<"\t ee -- 45 OR 54, mm -- 86 or 90, em -- 45 OR 86 OR 126"<<std::endl;
  }

  // Check that prescale is OK
  if (prescale < 1) {
    cout << "Illegal Prescale = " << prescale << endl;
    return 0;
  }
  float probOfKeeping = 1./prescale;

  // Initialize the random number generator for the prescale
  TRandom2 * r = new TRandom2();

  //Event Loop
  int nAfterPrescale = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=chain->GetEntries();
  unsigned int nEventsTotal = 0;
  
  // file loop
  TIter fileIter(listOfFiles);
  int nAllEvents = 0;
  map<int,int> m_events;
  while(TChainElement *currentFile = (TChainElement*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    Init(tree);
    unsigned int nEntries = tree->GetEntries();
    unsigned int nLoop = nEntries;

    unsigned int z;

    for( z = 0; z < nLoop; z++) {
      nAllEvents++;
      // Progress feedback to the user
      int iz = nAllEvents/10000;
      if (nAllEvents-10000*iz == 0) cout << "Processing event " << nAllEvents+1 
				       << " of sample " << prefix << endl;
       
      // Prescale
      if (prescale != 1) {
	if (r->Uniform(1) > probOfKeeping) continue;
      }
      nAfterPrescale++;
      GetEntry(z);
      ++nEventsTotal;


      //check if it's a correct genp-event
      std::string prefixStr(prefix);
      if (prefixStr == "ttdil" && genpCountPDGId(11,13,15) != 2) continue;
      if (prefixStr == "ttotr" && genpCountPDGId(11,13,15) == 2) continue;
      if (prefixStr == "DYee" && genpCountPDGId(11) != 2) continue;
      if (prefixStr == "DYmm" && genpCountPDGId(13) != 2) continue;
      if (prefixStr == "DYtautau" && genpCountPDGId(15) != 2) continue;

      //decide weather or not the event passed
      bool eventPassed = false;
      

      std::vector<unsigned int> goodHyps(0);

      bool hltIsoEle18_L1R = ((evt_HLT2() & (1<<(45-32))) != 0);
      bool hltLooseIsoEle15_LW_L1R = ((evt_HLT2() & (1<<(47-32))) != 0);
      bool hltDoubleIsoEle12_L1R = ((evt_HLT2() & (1<<(54-32))) != 0); 
      bool hltMu9 = ((evt_HLT3() & (1<<(82-64))) != 0);
      bool hltMu15_L1Mu7 = ((evt_HLT3() & (1<<(86-64))) != 0); 
      bool hltDoubleMu3 = ((evt_HLT3() & (1<<(90-64))) != 0);
      bool hltIsoEle10_Mu10_L1R = ((evt_HLT4() & (1<<(126-96))) != 0);

      for(unsigned int hypIdx = 0; hypIdx < hyp_p4().size(); hypIdx++) {
       
	if (applyZWindow && fabs(hyp_p4().at(hypIdx).mass()-91)> 15) continue;
      
              
// 	// Triggers (comment out the unused ones)
// 	bool HLT1ElectronRelaxed     = ( (evt_HLT2() & 0x2) != 0);
// 	// bool HLT2ElectronRelaxed     = ( (evt_HLT2() & 0x8) != 0);
// 	bool HLT1MuonNonIso          = ( (evt_HLT2() & 0x8000) != 0);
// 	// bool HLT2MuonNonIso          = ( (evt_HLT2() & 0x20000) != 0);
// 	// bool HLTXElectronMuonRelaxed = ( (evt_HLT3() & 0x80000) != 0); 

       

	// if (hyp_type().at(hypIdx) == 3 && !(HLT1ElectronRelaxed)) continue; // mm
	// 	if (hyp_type().at(hypIdx) == 0 && !(HLT1MuonNonIso))      continue; // ee
	// 	if (  (hyp_type().at(hypIdx) == 1 || hyp_type().at(hypIdx) == 2) && 
	// 	      !(HLT1ElectronRelaxed || HLT1MuonNonIso)) continue; // emu
	//if(hyp_type().at(hypIdx) !=1 && hyp_type().at(hypIdx) !=2 ) continue;
       
       
	unsigned int i_lt = hyp_lt_index().at(hypIdx);
	unsigned int i_ll = hyp_ll_index().at(hypIdx);
	
	// Cut on lepton Pt
// 	if (hyp_lt_p4().at(hypIdx).pt() < 20.0) continue;
// 	if (hyp_ll_p4().at(hypIdx).pt() < 20.0) continue;

	if(dilepMassVetoCut) {
	  // Z mass veto using hyp_leptons for ee and mumu final states
	  if (hyp_type().at(hypIdx) == 0 || hyp_type().at(hypIdx) == 3) {
	    if (inZmassWindow(hyp_p4().at(hypIdx).mass())) continue;
	  }
      
	  // Z veto using additional leptons in the event
	  if (additionalZveto()) continue;
	}
        if(dilepMassVetoCutTTDil08) {
          // Z mass veto using hyp_leptons for ee and mumu final states
          if (hyp_type().at(hypIdx) == 0 || hyp_type().at(hypIdx) == 3) {
            if (inZmassWindow(hyp_p4().at(hypIdx).mass())) continue;
          }
    
          // Z veto using additional leptons in the event
          if (additionalZveto(true)) continue; //"true" to use TTDil lepton selections                                                
        }

      
	// Dima's MET requirement
	bool pass2MetPassed = pass2Met(hypIdx);
	if (! pass2MetPassed && METcut) continue;
	if ( pass2MetPassed && METveto ) continue;

	// ! for TTDil analysis this should be made for the event-qualifying hyp only
	if (!fillMaxWeightDilOnly && metBaselineSelectionTTDil08){
	  if (! passPatMet_OF20_SF30(hypIdx)) continue;
	}

	// this is for per-hypothesis choice
	if (! fillMaxWeightDilOnly && applyTriggersMu9orLisoE15){
	  if (hyp_type().at(hypIdx) == 0 && ! (hltMu9) ) continue;
	  if ((hyp_type().at(hypIdx) == 1 || hyp_type().at(hypIdx) == 2) && ! (hltMu9 || hltLooseIsoEle15_LW_L1R)) continue;
	  if (hyp_type().at(hypIdx) == 3 && ! hltLooseIsoEle15_LW_L1R) continue; 
	}
	if (! fillMaxWeightDilOnly && applyTriggersTTDil08JanTrial){
	  if (hyp_type().at(hypIdx) == 0 && ! (hltMu15_L1Mu7 || hltDoubleMu3) ) continue;
	  if ((hyp_type().at(hypIdx) == 1 || hyp_type().at(hypIdx) == 2) 
	      && ! (hltIsoEle18_L1R || hltMu15_L1Mu7 || hltIsoEle10_Mu10_L1R)) continue;
	  if (hyp_type().at(hypIdx) == 3 && ! (hltIsoEle18_L1R || hltDoubleIsoEle12_L1R) ) continue; 
	}

	if(idcuts) {
	  // Muon quality cuts, excluding isolation
	  if (abs(hyp_lt_id().at(hypIdx)) == 13 && !goodMuonWithoutIsolation(i_lt) ) continue;
	  if (abs(hyp_ll_id().at(hypIdx)) == 13 && !goodMuonWithoutIsolation(i_ll) ) continue;
      
	  // Electron quality cuts, excluding isolation
	  if (abs(hyp_lt_id().at(hypIdx)) == 11 && !goodElectronWithoutIsolation(i_lt) ) continue;
	  if (abs(hyp_ll_id().at(hypIdx)) == 11 && !goodElectronWithoutIsolation(i_ll) ) continue;
	}

	if(isolationcuts) {
	  if (!passDilAntiIsolation(isoLooseMode, hypIdx)) continue; 
	}
    
  
	if (applyMuTag && ! haveExtraMuon(hypIdx)) continue;
	if (applyMuTag5 && ! haveExtraMuon5(hypIdx)) continue;


	if (lepton20Eta2p4DilSelection){
	  // Muon pt eta cuts
	  if (abs(hyp_lt_id().at(hypIdx)) == 13 && !muon20Eta2p4(i_lt) ) continue;
	  if (abs(hyp_ll_id().at(hypIdx)) == 13 && !muon20Eta2p4(i_ll) ) continue;

	  // Electron pt, eta cuts
	  if (abs(hyp_lt_id().at(hypIdx)) == 11 && !electron20Eta2p4(i_lt) ) continue;
	  if (abs(hyp_ll_id().at(hypIdx)) == 11 && !electron20Eta2p4(i_ll) ) continue;	  
	}

	if (looseDilSelectionNoIsoTTDil08){
	  // Muon quality cuts, no isolation
	  if (abs(hyp_lt_id().at(hypIdx)) == 13 && !looseMuonSelectionNoIsoTTDil08(i_lt) ) continue;
	  if (abs(hyp_ll_id().at(hypIdx)) == 13 && !looseMuonSelectionNoIsoTTDil08(i_ll) ) continue;

	  // Electron quality cuts, no isolation
	  if (abs(hyp_lt_id().at(hypIdx)) == 11 && !looseElectronSelectionNoIsoTTDil08(i_lt) ) continue;
	  if (abs(hyp_ll_id().at(hypIdx)) == 11 && !looseElectronSelectionNoIsoTTDil08(i_ll) ) continue;
	}

	if (looseDilSelectionTTDil08){
	  // Muon quality cuts, loose isolation
	  if (abs(hyp_lt_id().at(hypIdx)) == 13 && !looseMuonSelectionTTDil08(i_lt) ) continue;
	  if (abs(hyp_ll_id().at(hypIdx)) == 13 && !looseMuonSelectionTTDil08(i_ll) ) continue;

	  // Electron quality cuts, loose isolation
	  if (abs(hyp_lt_id().at(hypIdx)) == 11 && !looseElectronSelectionTTDil08(i_lt) ) continue;
	  if (abs(hyp_ll_id().at(hypIdx)) == 11 && !looseElectronSelectionTTDil08(i_ll) ) continue;
	}

	if (leptonIsolationDilSelectionTTDil08){
	  // muon passes ttdil08 iso cuts
	  if (abs(hyp_lt_id().at(hypIdx)) == 13 && !passMuonIsolationTTDil08(i_lt) ) continue;
	  if (abs(hyp_ll_id().at(hypIdx)) == 13 && !passMuonIsolationTTDil08(i_ll) ) continue;

	  // Electron passes ttdil08 iso cuts
	  if (abs(hyp_lt_id().at(hypIdx)) == 11 && !passElectronIsolationTTDil08(i_lt) ) continue;
	  if (abs(hyp_ll_id().at(hypIdx)) == 11 && !passElectronIsolationTTDil08(i_ll) ) continue;

	}

	if (osSelection){
	  if ( hyp_lt_id().at(hypIdx) * hyp_ll_id().at(hypIdx) > 0 ) continue;
	}

	goodHyps.push_back(hypIdx);
	//done with cuts on hyps		   

	eventPassed = true;
      }
      
      unsigned int nGoodHyps = goodHyps.size();
      std::vector<float> hypWeights(0);
      unsigned int maxWeightIndex = 0;
      int strasbourgDilType = -1;
      if (nGoodHyps > 0){
	//some debug printouts
	float maxWeight = -1;
	
	for (unsigned int hypIdxL=0; hypIdxL < nGoodHyps; ++hypIdxL){
	  unsigned int hypIdx = goodHyps[hypIdxL];
	  float hypWeight_lt = 0;
	  float hypWeight_ll = 0;
	  float hypWeight_iso = 0;
	  float hypWeight = 0;
	  unsigned int i_lt = hyp_lt_index().at(hypIdx);
	  unsigned int i_ll = hyp_ll_index().at(hypIdx);
	  //ad-hoc selection of weights
	  if (abs(hyp_lt_id().at(hypIdx)) == 11){
	    //I want to select "trk & cal"-isolated ones
	    //shift by 0.25 to be positive-definite
	    hypWeight_iso += (electronTrkIsolation(i_lt)*electronCalIsolation(i_lt)-0.25);
	    if (els_tightId().at(i_lt)) hypWeight_lt += 0.2;
	  }
	  if (abs(hyp_lt_id().at(hypIdx)) == 13){
	    //I want to select "trk & cal"-isolated ones
	    //shift by 0.25 to be positive-definite
	    hypWeight_iso += (muonTrkIsolation(i_lt)*muonCalIsolation(i_lt)-0.25);
	    hypWeight_lt += 0.4;
	  }
	  if (abs(hyp_ll_id().at(hypIdx)) == 11){
	    //I want to select "trk & cal"-isolated ones
	    //shift by 0.25 to be positive-definite
	    hypWeight_iso *= (electronTrkIsolation(i_ll)*electronCalIsolation(i_ll)-0.25);
	    if (els_tightId().at(i_ll)) hypWeight_ll += 0.2;
	  }
	  if (abs(hyp_ll_id().at(hypIdx)) == 13){
	    //I want to select "trk & cal"-isolated ones
	    //shift by 0.25 to be positive-definite
	    hypWeight_iso *= (muonTrkIsolation(i_ll)*muonCalIsolation(i_ll)-0.25);
	    hypWeight_ll += 0.4;
	  }
	  float pt_lt = hyp_lt_p4().at(hypIdx).pt();
	  float pt_ll = hyp_ll_p4().at(hypIdx).pt();
	  hypWeight_lt += (1. - 20./pt_lt*20./pt_lt);
	  hypWeight_ll += (1. - 20./pt_ll*20./pt_ll);

	  hypWeight = hypWeight_ll*hypWeight_lt*hypWeight_iso; //again, desire to have both good

	  if (hypWeight > maxWeight){
	    maxWeight = hypWeight;
	    maxWeightIndex = hypIdx;
	  }

	  //Now let's implement the Strasbourg-type disambiguation/dispatch
	  //ee
	  {
	    std::vector<unsigned int> looseEls(0);
	    std::vector<unsigned int> looseMus(0);
	    for (unsigned int iEl =0; iEl < els_p4().size(); ++iEl){
	      if (looseElectronSelectionTTDil08(iEl)){
		looseEls.push_back(iEl);
	      }
	    }
	    for (unsigned int iMu =0; iMu < mus_p4().size(); ++iMu){
	      if (looseMuonSelectionTTDil08(iMu)){
		looseMus.push_back(iMu);
	      }
	    }
	    
	    bool pass_elec = false;
	    if (looseEls.size()>1){
	      if (els_charge().at(looseEls[0]) != els_charge().at(looseEls[1])){
		pass_elec = true;
	      }
	      if (looseMus.size()>0 && mus_p4().at(looseMus[0]).pt() > els_p4().at(looseEls[1]).pt()) pass_elec = false;
	      if (looseMus.size()>0 && 
		  ( ( muonTrkIsolation(looseMus[0]) > electronTrkIsolation(looseEls[0]) 
		      && mus_charge().at(looseMus[0]) != els_charge().at(looseEls[0]) )
		    || ( muonTrkIsolation(looseMus[0]) > electronTrkIsolation(looseEls[1])
			 && mus_charge().at(looseMus[0]) != els_charge().at(looseEls[0]))
		    )
		  ) pass_elec = false; 
	    }
	    bool pass_muon = false;
	    if (looseMus.size()>1){
	      for (unsigned int iMu=0; iMu < looseMus.size(); ++iMu){
		for (unsigned int jMu=iMu+1; jMu < looseMus.size(); ++jMu){
		  if (mus_charge().at(looseMus[iMu]) != mus_charge().at(looseMus[jMu])) pass_muon = true;
		}
	      }
	      if (looseEls.size()>0 && els_p4().at(looseEls[0]).pt() > mus_p4().at(looseMus[1]).pt()
		  && mus_charge().at(looseMus[1]) != els_charge().at(looseEls[0])) pass_muon = false;
	    }
	    bool pass_elecmuon = false;
	    if (looseMus.size() > 0 && looseEls.size() > 0){
	      if (! pass_elec && ! pass_muon ){
		if (mus_charge().at(looseMus[0]) != els_charge().at(looseEls[0])) pass_elecmuon = true;
		if (! pass_elecmuon && looseEls.size()>1){
		  if (mus_charge().at(looseMus[0]) != els_charge().at(looseEls[0])) pass_elecmuon = true;
		}
	      }
	    }

	    unsigned int passStatus = 0;
	    if (pass_muon) passStatus++;
	    if (pass_elecmuon) passStatus++;
	    if (pass_elec) passStatus++;
	    if (passStatus > 1) std::cout<<"ERROR: inconsistent assignment"<<std::endl;
	    if (passStatus == 1){
	      if (pass_muon) strasbourgDilType = 0;
	      if (pass_elecmuon) strasbourgDilType = 1;
	      if (pass_elec) strasbourgDilType = 2;
	    }
	  }
	}
	int genpDilType = genpDileptonType();
	if (genpDilType>=0 && 1>2){ std::cout<<"Dil type "<<genpDilType<<std::endl;
	  if (nGoodHyps > 1){
	    int maxWeightType = hyp_type().at(maxWeightIndex);
	    if ((maxWeightType == 0 && genpDilType == 0)
		|| ( (maxWeightType == 1 || maxWeightType == 2) && genpDilType == 1)
		|| (maxWeightType == 3 && genpDilType == 2)){
	      std::cout<<"Dil type "<<genpDilType<<" ; Strasbourg dil type "<<strasbourgDilType 
		       <<" assigned correctly by maxWeight method";
	      std::cout<<" out of"; for(unsigned int iih=0;iih<nGoodHyps;++iih)std::cout<<" "<<hyp_type().at(goodHyps[iih]);
	      std::cout<<std::endl;
	    } else {
	      std::cout<<"Dil type "<<genpDilType<<" ; Strasbourg dil type "<<strasbourgDilType 
		       <<" assigned incorrectly by maxWeight method";
	      std::cout<<" out of"; for(unsigned int iih=0;iih<nGoodHyps;++iih)std::cout<<" "<<hyp_type().at(goodHyps[iih]);
	      std::cout<<std::endl;	    
	    }
	  }
	}
      }

      // ! event level cut here, can reset the eventPassed to false
      if (fillMaxWeightDilOnly && metBaselineSelectionTTDil08 && nGoodHyps > 0){
	if (! passPatMet_OF20_SF30(maxWeightIndex)){
	  eventPassed = false;
	  continue;
	}
      }
      // this is for per-hypothesis choice
      if ( fillMaxWeightDilOnly && applyTriggersMu9orLisoE15 && nGoodHyps>0){
	if (hyp_type().at(maxWeightIndex) == 0 && ! (hltMu9) ) continue;
	if ((hyp_type().at(maxWeightIndex) == 1 || hyp_type().at(maxWeightIndex) == 2) 
	    && ! (hltMu9 || hltLooseIsoEle15_LW_L1R)) continue;
	if (hyp_type().at(maxWeightIndex) == 3 && ! hltLooseIsoEle15_LW_L1R) continue; 
      }
      if ( fillMaxWeightDilOnly && applyTriggersTTDil08JanTrial && nGoodHyps >0){
	if (hyp_type().at(maxWeightIndex) == 0 && ! (hltMu15_L1Mu7 || hltDoubleMu3) ) continue;
	if ((hyp_type().at(maxWeightIndex) == 1 || hyp_type().at(maxWeightIndex) == 2) 
	    && ! (hltIsoEle18_L1R || hltMu15_L1Mu7 || hltIsoEle10_Mu10_L1R)) continue;
	if (hyp_type().at(maxWeightIndex) == 3 && ! (hltIsoEle18_L1R || hltDoubleIsoEle12_L1R) ) continue; 
      }

      for(unsigned int hypIdxL=0; hypIdxL< nGoodHyps; ++ hypIdxL){
	unsigned int hypIdx = goodHyps[hypIdxL];
	if (fillMaxWeightDilOnly && hypIdx != maxWeightIndex) continue;
	//count the number of tight leptons:
	// This is an FKW variable which I turned off since I do not understand it 
	// and also uses simpleIdPlus, so it should be fixed up before being turned back on.
	//    int inumTightLep = numTightLeptons();
	int inumTightLep = 0;    

	// The event weight including the kFactor (scaled to 1 fb-1) and the prescale
	//float weight = evt_scale1fb * kFactor * prescale;
	//float weight = CalculateWeight(evt_CSA07Process(), evt_scale1fb(), kFactor, prescale);
	float weight = kFactor*evt_scale1fb(); // /100; //10pb^-1

	if ( (prefixStr == "ppMuX" || prefixStr == "EM") 
	     && (hyp_type().at(hypIdx) == 1 || hyp_type().at(hypIdx) == 2)) weight *= 0.5;//this isn't quite right :(
	//and works only if both em and ppmux are in play
      
	// If we made it to here, we passed all cuts and we are ready to fill
	m_events.insert(pair<int,int>(evt_event(), 1));

	int myType = 99;
	if (hyp_type().at(hypIdx) == 3) myType = 0;  // ee
	if (hyp_type().at(hypIdx) == 0) myType = 1;  // mm
	if (hyp_type().at(hypIdx) == 1 || hyp_type().at(hypIdx) == 2) myType=2; // em
	if (myType == 99) {
	  cout << "YUK:  unknown dilepton type = " << hyp_type().at(hypIdx) << endl;
	  continue;
	}

	// Now we have to manipulate the jets.
	// For the old jet selection (odjets=true) we use the hyp_jets
	// However: the old ntuples had uncorrected jets. Now the jets
	// are corrected, so we have to undo the correction before filling
	// the jet histograms.
	// For the new jet selection, we will need to re-count the jets.
	// General solution: we will make a new vector, new_hyp_jets_p4
	// with the hyp_jets for the selection that we have chosen.
	// In the case of the old jet selection, the p4 will be uncorrected.
	//
	int new_hyp_njets=0;  // jet count
	VofP4 jp4;            // vector of jets 
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > blah; // temp variable
	// First the case where we take the default hyp_jets
	if (oldjets) {
	  new_hyp_njets = hyp_njets().at(hypIdx);
	  for (unsigned int ijet=0; 
	       ijet<(unsigned int)(hyp_njets().at(hypIdx)); 
	       ijet++) {
	    blah = hyp_jets_pat_noCorrF().at(hypIdx).at(ijet) * hyp_jets_p4().at(hypIdx).at(ijet);
	    jp4.push_back(blah);
	  }
	} else {
	  // Look among the hyp_jets
	  for (unsigned int ijet=0; ijet<(unsigned int)(hyp_njets().at(hypIdx)); ijet++) {
	    blah = hyp_jets_p4().at(hypIdx).at(ijet);
	    if (blah.pt() > 30 && abs(blah.eta()) < 2.4) {
	      jp4.push_back(blah);
	      new_hyp_njets++;
	    }
	  }
	  // Now look among the other jets
	  for (unsigned int ijet=0; ijet<hyp_other_jets_p4().at(hypIdx).size(); ijet++) {
	    blah = hyp_other_jets_p4().at(hypIdx).at(ijet);
	    if (blah.pt() > 30 && abs(blah.eta()) < 2.4) {
	      jp4.push_back(blah);
	      new_hyp_njets++;
	    }
	  }
	}
	VofP4* new_hyp_jets_p4 = &jp4;
			   
	//     // Last chance to reject...
	int arrNjets = min(new_hyp_njets, 2);


	// jet count
	hnJet[myType]->Fill(new_hyp_njets, weight);
	hnJet[3]->Fill(new_hyp_njets, weight);
	if( inumTightLep < 3) {
	  hnJetLepVeto[myType]->Fill(new_hyp_njets, weight);
	  hnJetLepVeto[3]->Fill(new_hyp_njets, weight);
	}
	numTightLep[myType][arrNjets]->Fill(inumTightLep,weight);
	numTightLep[3][arrNjets]->Fill(inumTightLep,weight);
      
      
	// lepton Pt
	if (abs(hyp_lt_id().at(hypIdx)) == 11) helePt[myType][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).pt(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 11) helePt[myType][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).pt(), weight);
	if (abs(hyp_lt_id().at(hypIdx)) == 13) hmuPt[myType][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).pt(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 13) hmuPt[myType][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).pt(), weight);
	if (abs(hyp_lt_id().at(hypIdx)) == 13) 
	  hmuPtFromSilicon[myType][arrNjets]->Fill(mus_trk_p4().at(hyp_lt_index().at(hypIdx)).pt(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 13)
	  hmuPtFromSilicon[myType][arrNjets]->Fill(mus_trk_p4().at(hyp_ll_index().at(hypIdx)).pt(), weight);
	hminLepPt[myType][arrNjets]->Fill(min(hyp_ll_p4().at(hypIdx).pt(), hyp_lt_p4().at(hypIdx).pt()), weight);
	hmaxLepPt[myType][arrNjets]->Fill(max(hyp_ll_p4().at(hypIdx).pt(), hyp_lt_p4().at(hypIdx).pt()), weight );
    
	if (abs(hyp_lt_id().at(hypIdx)) == 11) helePt[3][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).pt(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 11) helePt[3][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).pt(), weight);
	if (abs(hyp_lt_id().at(hypIdx)) == 13) hmuPt[3][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).pt(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 13) hmuPt[3][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).pt(), weight);
	if (abs(hyp_lt_id().at(hypIdx)) == 13) 
	  hmuPtFromSilicon[3][arrNjets]->Fill(mus_trk_p4().at(hyp_lt_index().at(hypIdx)).pt(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 13) 
	  hmuPtFromSilicon[3][arrNjets]->Fill(mus_trk_p4().at(hyp_ll_index().at(hypIdx)).pt(), weight);
	hminLepPt[3][arrNjets]->Fill(min(hyp_ll_p4().at(hypIdx).pt(), hyp_lt_p4().at(hypIdx).pt()), weight);
	hmaxLepPt[3][arrNjets]->Fill(max(hyp_ll_p4().at(hypIdx).pt(), hyp_lt_p4().at(hypIdx).pt()), weight );


	// lepton Phi
	if (abs(hyp_lt_id().at(hypIdx)) == 11) helePhi[myType][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).phi(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 11) helePhi[myType][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).phi(), weight);
	if (abs(hyp_lt_id().at(hypIdx)) == 13) hmuPhi[myType][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).phi(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 13) hmuPhi[myType][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).phi(), weight);
    
	if (abs(hyp_lt_id().at(hypIdx)) == 11) helePhi[3][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).phi(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 11) helePhi[3][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).phi(), weight);
	if (abs(hyp_lt_id().at(hypIdx)) == 13) hmuPhi[3][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).phi(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 13) hmuPhi[3][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).phi(), weight);
    
	// dilepton mass
	hdilMass[myType][arrNjets]->Fill(hyp_p4().at(hypIdx).mass(), weight);
	hdilMassTightWindow[myType][arrNjets]->Fill(hyp_p4().at(hypIdx).mass(), weight);
	hdilMass[3][arrNjets]->Fill(hyp_p4().at(hypIdx).mass(), weight);
	hdilMassTightWindow[3][arrNjets]->Fill(hyp_p4().at(hypIdx).mass(), weight);
    
	// delta phi btw leptons
	double dphi = fabs(hyp_lt_p4().at(hypIdx).phi() - hyp_ll_p4().at(hypIdx).phi());
	if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
	hdphiLep[myType][arrNjets]->Fill(dphi, weight);
	hdphiLep[3][arrNjets]->Fill(dphi, weight);
			   
	// dphill vs mll, i.e. the 2d correlation between the previous two variables
	hdphillvsmll[myType][arrNjets]->Fill(hyp_p4().at(hypIdx).mass(), dphi, weight);
	hdphillvsmll[3][arrNjets]->Fill(hyp_p4().at(hypIdx).mass(), dphi, weight);
 
	// lepton Eta
	if (abs(hyp_lt_id().at(hypIdx)) == 11) heleEta[myType][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).eta(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 11) heleEta[myType][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).eta(), weight);
	if (abs(hyp_lt_id().at(hypIdx)) == 13) hmuEta[myType][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).eta(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 13) hmuEta[myType][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).eta(), weight);
    
	if (abs(hyp_lt_id().at(hypIdx)) == 11) heleEta[3][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).eta(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 11) heleEta[3][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).eta(), weight);
	if (abs(hyp_lt_id().at(hypIdx)) == 13) hmuEta[3][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).eta(), weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 13) hmuEta[3][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).eta(), weight);
       
      
	// electron trk isolation 
	double temp_lt_iso = hyp_lt_iso().at(hypIdx);  // so that min works
	double temp_ll_iso = hyp_ll_iso().at(hypIdx);  // so that min works
	if (abs(hyp_lt_id().at(hypIdx)) == 11) heleSumPt[myType][arrNjets]->Fill(min(temp_lt_iso,24.99),weight);
	if (abs(hyp_lt_id().at(hypIdx)) == 11) heleSumPt[3][arrNjets]->Fill(min(temp_lt_iso,24.99),weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 11) heleSumPt[myType][arrNjets]->Fill(min(temp_ll_iso,24.99),weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 11) heleSumPt[3][arrNjets]->Fill(min(temp_ll_iso,24.99),weight);

	// muon trk isolation
	if (abs(hyp_lt_id().at(hypIdx)) == 13) hmuSumPt[myType][arrNjets]->Fill(min(temp_lt_iso,24.99),weight);
	if (abs(hyp_lt_id().at(hypIdx)) == 13) hmuSumPt[3][arrNjets]->Fill(min(temp_lt_iso,24.99),weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 13) hmuSumPt[myType][arrNjets]->Fill(min(temp_ll_iso,24.99),weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 13) hmuSumPt[3][arrNjets]->Fill(min(temp_ll_iso,24.99),weight);
      
      
	// muon trk+calo isolation
	double combIso_lt = -1.;
	double combIso_ll = -1.;
	if (abs(hyp_lt_id().at(hypIdx)) == 13)
	  combIso_lt = mus_iso03_sumPt().at(hyp_lt_index().at(hypIdx))
	    +mus_iso03_emEt().at(hyp_lt_index().at(hypIdx))
	    +mus_iso03_hadEt().at(hyp_lt_index().at(hypIdx));
	if (abs(hyp_ll_id().at(hypIdx)) == 13)
	  combIso_ll = mus_iso03_sumPt().at(hyp_ll_index().at(hypIdx))
	    +mus_iso03_emEt().at(hyp_ll_index().at(hypIdx))
	    +mus_iso03_hadEt().at(hyp_ll_index().at(hypIdx));
	if (abs(hyp_lt_id().at(hypIdx)) == 13) hmuSumIso[myType][arrNjets]->Fill(min(combIso_lt,24.99),weight);
	if (abs(hyp_lt_id().at(hypIdx)) == 13) hmuSumIso[3][arrNjets]->Fill(min(combIso_lt,24.99),weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 13) hmuSumIso[myType][arrNjets]->Fill(min(combIso_ll,24.99),weight);
	if (abs(hyp_ll_id().at(hypIdx)) == 13) hmuSumIso[3][arrNjets]->Fill(min(combIso_ll,24.99),weight);


	// Relative isolation... muons
	if (abs(hyp_lt_id().at(hypIdx)) == 13) {
	  double thisSum =  mus_iso03_sumPt().at(hyp_lt_index().at(hypIdx)) +  
	    mus_iso03_emEt().at(hyp_lt_index().at(hypIdx))  +
	    mus_iso03_hadEt().at(hyp_lt_index().at(hypIdx));
	  double thisPt  = mus_p4().at(hyp_lt_index().at(hypIdx)).pt();
	  double temp    = thisPt / (thisPt+thisSum);
	  hmuRelIso[myType][arrNjets]->Fill(temp, weight);
	  hmuRelIso[3][arrNjets]->Fill(temp, weight);
	}
	if (abs(hyp_ll_id().at(hypIdx)) == 13) {
	  double thisSum =  mus_iso03_sumPt().at(hyp_ll_index().at(hypIdx)) +  
	    mus_iso03_emEt().at(hyp_ll_index().at(hypIdx))  +
	    mus_iso03_hadEt().at(hyp_ll_index().at(hypIdx));
	  double thisPt  = mus_p4().at(hyp_ll_index().at(hypIdx)).pt();
	  double temp    = thisPt / (thisPt+thisSum);
	  hmuRelIso[myType][arrNjets]->Fill(temp, weight);
	  hmuRelIso[3][arrNjets]->Fill(temp, weight);
	}


	// Relative isolation... electrons
	if (abs(hyp_lt_id().at(hypIdx)) == 11) {
	  double thisSum =  hyp_lt_iso().at(hypIdx);
	  double thisPt  = hyp_lt_p4().at(hypIdx).pt();
	  double temp    = thisPt / (thisPt+thisSum);
	  heleRelIso[myType][arrNjets]->Fill(temp, weight);
	  heleRelIso[3][arrNjets]->Fill(temp, weight);
	}
	if (abs(hyp_ll_id().at(hypIdx)) == 11) {
	  double thisSum =  hyp_ll_iso().at(hypIdx);
	  double thisPt  = hyp_ll_p4().at(hypIdx).pt();
	  double temp    = thisPt / (thisPt+thisSum);
	  heleRelIso[myType][arrNjets]->Fill(temp, weight);
	  heleRelIso[3][arrNjets]->Fill(temp, weight);
	}
      
	// dilepton pt
	hdilPt[myType][arrNjets]->Fill(hyp_p4().at(hypIdx).pt(), weight);
	hdilPt[3][arrNjets]->Fill(hyp_p4().at(hypIdx).pt(), weight);
    
	// Met and Met phi
	hmet[myType][arrNjets]->Fill(hyp_met().at(hypIdx), weight);      
	hmetPhi[myType][arrNjets]->Fill(hyp_metPhi().at(hypIdx), weight);      
	hmet[3][arrNjets]->Fill(hyp_met().at(hypIdx), weight);      
	hmetPhi[3][arrNjets]->Fill(hyp_metPhi().at(hypIdx), weight);      
    
	// Met vs dilepton Pt
	hmetVsDilepPt[myType][arrNjets]->Fill(hyp_met().at(hypIdx), hyp_p4().at(hypIdx).pt(), weight);
	hmetVsDilepPt[3][arrNjets]->Fill(hyp_met().at(hypIdx), hyp_p4().at(hypIdx).pt(), weight);
    
	// Met over dilepton Pt vs deltaphi btw the two
	double dphi2 = fabs(hyp_p4().at(hypIdx).phi() - hyp_metPhi().at(hypIdx));
	if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
	dphi2 = TMath::Pi() - dphi2;  // changed the definition CC 28 March 08
	hmetOverPtVsDphi[myType][arrNjets]->Fill(hyp_met().at(hypIdx)/hyp_p4().at(hypIdx).pt(), dphi2, weight);
	hmetOverPtVsDphi[3][arrNjets]->Fill(hyp_met().at(hypIdx)/hyp_p4().at(hypIdx).pt(), dphi2, weight);
    
      

	// Make a vector of sorted jets, fill jet histograms
	if (new_hyp_njets > 0) {
	  vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > my_hyp_jets_p4(*new_hyp_jets_p4);
	  sort(my_hyp_jets_p4.begin(), my_hyp_jets_p4.end(), comparePt);   // sort them by Pt
	  hptJet1[myType][arrNjets]->Fill(my_hyp_jets_p4[0].Pt(), weight);
	  hptJet1[3][arrNjets]->Fill(my_hyp_jets_p4[0].Pt(), weight);
	  hetaJet1[myType][arrNjets]->Fill(my_hyp_jets_p4[0].Eta(), weight);
	  hetaJet1[3][arrNjets]->Fill(my_hyp_jets_p4[0].Eta(), weight);
	  if (new_hyp_njets > 1) {
	    hptJet2[myType][arrNjets]->Fill(my_hyp_jets_p4[1].Pt(), weight);
	    hptJet2[3][arrNjets]->Fill(my_hyp_jets_p4[1].Pt(), weight);
	    hetaJet2[myType][arrNjets]->Fill(my_hyp_jets_p4[1].Eta(), weight);
	    hetaJet2[3][arrNjets]->Fill(my_hyp_jets_p4[1].Eta(), weight);
	  }
	  if (new_hyp_njets > 2) {
	    hptJet3[myType][arrNjets]->Fill(my_hyp_jets_p4[2].Pt(), weight);
	    hptJet3[3][arrNjets]->Fill(my_hyp_jets_p4[2].Pt(), weight);
	    hetaJet3[myType][arrNjets]->Fill(my_hyp_jets_p4[2].Eta(), weight);
	    hetaJet3[3][arrNjets]->Fill(my_hyp_jets_p4[2].Eta(), weight);
	  }
	  if (new_hyp_njets > 3) {
	    hptJet4[myType][arrNjets]->Fill(my_hyp_jets_p4[3].Pt(), weight);
	    hptJet4[3][arrNjets]->Fill(my_hyp_jets_p4[3].Pt(), weight);
	    hetaJet4[myType][arrNjets]->Fill(my_hyp_jets_p4[3].Eta(), weight);
	    hetaJet4[3][arrNjets]->Fill(my_hyp_jets_p4[3].Eta(), weight);
	  }
	}
      }//hypothesis loop
	  
    }//event loop
  }//file loop

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  std::cout<<"Done with "<<prefix<<std::endl;
  outf.close();
  std::cout<<"Closed outf "<<std::endl;
  rootdir = gDirectory->GetDirectory("Rint:"); 
  rootdir->cd(); 

  return 0;
}




void ttDilCounts_looper::bookHistos(char *prefix) {

  //  Book histograms...
  //  Naming Convention:
  //  Prefix comes from the sample and it is passed to the scanning function
  //  Suffix is "ee" "em" "em" "all" which depends on the final state
  //  For example: histogram named tt_hnJet_ee would be the Njet distribution
  //  for the ee final state in the ttbar sample.
  
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  
  char *jetbins[5] = {"0", "1", "2", "3", "#geq 4"};
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  for (int i=0; i<4; i++) {
    for (int j=0; j<3; j++) {
      char *suffixall[4];
      suffixall[0] = "ee";
      suffixall[1] = "mm";
      suffixall[2] = "em";
      suffixall[3] = "all";
  
      if (j == 0){
	hnJet[i] = new TH1F(Form("%s_hnJet_%s",prefix,suffixall[i]),Form("%s_nJet_%s",prefix,suffixall[i]),
			    5,0.,5.);	
	hnJet[i]->SetDirectory(rootdir);
	hnJet[i]->GetXaxis()->SetTitle("nJets");
	
	hnJetLepVeto[i] = new TH1F(Form("%s_hnJetLepVeto_%s",prefix,suffixall[i]),Form("%s_nJetLepVeto_%s",prefix,suffixall[i]),
				   5,0.,5.);	
	hnJetLepVeto[i]->SetDirectory(rootdir);
	hnJetLepVeto[i]->GetXaxis()->SetTitle("nJets");
	
	
	for(int k = 0; k<5; k++) {
	  hnJet[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
	  hnJet[i]->GetXaxis()->SetLabelSize(0.07);
	  
	  hnJetLepVeto[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
	  hnJetLepVeto[i]->GetXaxis()->SetLabelSize(0.07);
	}
      }
    
      char *suffix[4];
      char *njetCh[3] = { "0j", "1j", "2j" };
      suffix[0] = Form("%s_ee", njetCh[j]);
      suffix[1] = Form("%s_mm", njetCh[j]);
      suffix[2] = Form("%s_em", njetCh[j]);
      suffix[3] = Form("%s_all", njetCh[j]);

      numTightLep[i][j] = new TH1F(Form("%s_numTightLep_%s",prefix,suffix[i]),Form("%s_numTightLep_%s",prefix,suffix[i]),
				   10,0.,10.);	
      numTightLep[i][j]->SetDirectory(rootdir);
      
      helePt[i][j] = new TH1F(Form("%s_helePt_%s",prefix,suffix[i]),Form("%s_elePt_%s",prefix,suffix[i]),
			      150,0.,150.);
      helePt[i][j]->SetDirectory(rootdir);
      helePt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
      
    
      hmuPt[i][j]  = new TH1F(Form("%s_hmuPt_%s",prefix,suffix[i]),Form("%s_muPt_%s",prefix,suffix[i]),
			      150,0.,150.);
      hmuPt[i][j]->SetDirectory(rootdir);
      hmuPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
      
    
      hmuPtFromSilicon[i][j]  = new TH1F(Form("%s_hmuPtFromSilicon_%s",prefix,suffix[i]),
					 Form("%s_muPtFromSilicon_%s",prefix,suffix[i]),150,0.,150.);
      hmuPtFromSilicon[i][j]->SetDirectory(rootdir);
      hmuPtFromSilicon[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
      
    
      hminLepPt[i][j]  = new TH1F(Form("%s_hminLepPt_%s",prefix,suffix[i]),
				  Form("%s_minLepPt_%s",prefix,suffix[i]),150,0.,150.);
      hminLepPt[i][j]->SetDirectory(rootdir);
      hminLepPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
      
    
      hmaxLepPt[i][j]  = new TH1F(Form("%s_hmaxLepPt_%s",prefix,suffix[i]),
				  Form("%s_maxLepPt_%s",prefix,suffix[i]),150,0.,150.);
      hmaxLepPt[i][j]->SetDirectory(rootdir);
      hmaxLepPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
      
    
      helePhi[i][j] = new TH1F(Form("%s_helePhi_%s",prefix,suffix[i]),Form("%s_elePhi_%s",prefix,suffix[i]),
			       50,-1*TMath::Pi(), TMath::Pi());
      helePhi[i][j]->SetDirectory(rootdir);
      helePhi[i][j]->GetXaxis()->SetTitle("#phi");
      

      hmuPhi[i][j]  = new TH1F(Form("%s_hmuPhi_%s",prefix,suffix[i]),Form("%s_muPhi_%s",prefix,suffix[i]),
			       50,-1*TMath::Pi(), TMath::Pi());
      hmuPhi[i][j]->SetDirectory(rootdir);
      hmuPhi[i][j]->GetXaxis()->SetTitle("#phi");
      
    
      hdphiLep[i][j]  = new TH1F(Form("%s_hdphiLep_%s",prefix,suffix[i]),Form("%s_dphiLep_%s",prefix,suffix[i]),
				 50,0., TMath::Pi());
      hdphiLep[i][j]->SetDirectory(rootdir);
      hdphiLep[i][j]->GetXaxis()->SetTitle("#delta#phi_{ll}");
      
      
      heleEta[i][j] = new TH1F(Form("%s_heleEta_%s",prefix,suffix[i]),Form("%s_eleEta_%s",prefix,suffix[i]),
			       60, -3., 3.);
      heleEta[i][j]->SetDirectory(rootdir);
      heleEta[i][j]->GetXaxis()->SetTitle("#eta");
      
	
      hmuEta[i][j]  = new TH1F(Form("%s_hmuEta_%s",prefix,suffix[i]),Form("%s_muEta_%s",prefix,suffix[i]),
			       60, -3., 3.);
      hmuEta[i][j]->SetDirectory(rootdir);
      hmuEta[i][j]->GetXaxis()->SetTitle("#eta");
      
 
      hdilMass[i][j] = new TH1F(Form("%s_hdilMass_%s",prefix,suffix[i]),Form("%s_dilMass_%s",prefix,suffix[i]),
				100, 0., 300.);
      hdilMass[i][j]->SetDirectory(rootdir);
      hdilMass[i][j]->GetXaxis()->SetTitle("Mass_{ll} (GeV)");
      

      hdilMassTightWindow[i][j] = new TH1F(Form("%s_hdilMassTightWindow_%s",prefix,suffix[i]),
					   Form("%s_dilMassTightWindow_%s",prefix,suffix[i]),
					   120, 60., 120.);
      hdilMassTightWindow[i][j]->SetDirectory(rootdir);
      hdilMassTightWindow[i][j]->GetXaxis()->SetTitle("Mass_{ll} (GeV)");
      
    
      hdilPt[i][j] = new TH1F(Form("%s_hdilPt_%s",prefix,suffix[i]),Form("%s_dilPt_%s",prefix,suffix[i]),
			      100, 0., 300.);
      hdilPt[i][j]->SetDirectory(rootdir);
      hdilPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      //changed binning from 2 GeV to 10 GeV
      hmet[i][j] = new TH1F(Form("%s_hmet_%s",prefix,suffix[i]),Form("%s_met_%s",prefix,suffix[i]),20,0.,200.);
      hmet[i][j]->SetDirectory(rootdir);
      hmet[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      hmetPhi[i][j] = new TH1F(Form("%s_hmetPhi_%s",prefix,suffix[i]),Form("%s_metPhi_%s",prefix,suffix[i]),
			       50,-1*TMath::Pi(), TMath::Pi());
      hmetPhi[i][j]->SetDirectory(rootdir);
      hmetPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hmetVsDilepPt[i][j] = new TH2F(Form("%s_hmetVsDilepPt_%s",prefix,suffix[i]),
				     Form("%s_metVsDilepPt_%s",prefix,suffix[i]),
				     100,0.,200.,100,0.,200.);
      hmetVsDilepPt[i][j]->SetDirectory(rootdir);
      hmetVsDilepPt[i][j]->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
      hmetVsDilepPt[i][j]->GetYaxis()->SetTitle("Met (GeV)");
    
      hmetOverPtVsDphi[i][j] = new TH2F(Form("%s_hmetOverPtVsDphi_%s",prefix,suffix[i]),
					Form("%s_metOverPtVsDphi_%s",prefix,suffix[i]),
					30,0.,3.,25,0.,TMath::Pi());
      hmetOverPtVsDphi[i][j]->SetDirectory(rootdir);
      hmetVsDilepPt[i][j]->GetXaxis()->SetTitle("#Delta#Phi");
      hmetVsDilepPt[i][j]->GetYaxis()->SetTitle("MET/Pt_{ll}");
    
      hdphillvsmll[i][j] = new TH2F(Form("%s_dphillvsmll_%s",prefix,suffix[i]),
				    Form("%s_dphillvsmll_%s",prefix,suffix[i]),
				    100,10.,210.,50,0., TMath::Pi());
      hdphillvsmll[i][j]->SetDirectory(rootdir);
      hdphillvsmll[i][j]->GetXaxis()->SetTitle("Mass_{ll} (GeV)");
      hdphillvsmll[i][j]->GetYaxis()->SetTitle("#delta#phi_{ll}");

      hptJet1[i][j] = new TH1F(Form("%s_hptJet1_%s",prefix,suffix[i]),Form("%s_ptJet1_%s",prefix,suffix[i]),
			       100, 0., 300.);
      hptJet1[i][j]->SetDirectory(rootdir);
      hptJet1[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJet2[i][j] = new TH1F(Form("%s_hptJet2_%s",prefix,suffix[i]),Form("%s_ptJet2_%s",prefix,suffix[i]),
			       100, 0., 300.);
      hptJet2[i][j]->SetDirectory(rootdir);
      hptJet2[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
  
      hptJet3[i][j] = new TH1F(Form("%s_hptJet3_%s",prefix,suffix[i]),Form("%s_ptJet3_%s",prefix,suffix[i]),
			       100, 0., 300.);
      hptJet3[i][j]->SetDirectory(rootdir);
      hptJet3[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
    
      hptJet4[i][j] = new TH1F(Form("%s_hptJet4_%s",prefix,suffix[i]),Form("%s_ptJet4_%s",prefix,suffix[i]),
			       100, 0., 300.);
      hptJet4[i][j]->SetDirectory(rootdir);
      hptJet4[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
    
      hetaJet1[i][j] = new TH1F(Form("%s_hetaJet1_%s",prefix,suffix[i]),Form("%s_etaJet1_%s",prefix,suffix[i]),
				50, -4., 4.);
      hetaJet1[i][j]->SetDirectory(rootdir);
      hetaJet1[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJet2[i][j] = new TH1F(Form("%s_hetaJet2_%s",prefix,suffix[i]),Form("%s_etaJet2_%s",prefix,suffix[i]),
				50, -4., 4.);
      hetaJet2[i][j]->SetDirectory(rootdir);
      hetaJet2[i][j]->GetXaxis()->SetTitle("#eta");
 
      hetaJet3[i][j] = new TH1F(Form("%s_hetaJet3_%s",prefix,suffix[i]),Form("%s_etaJet3_%s",prefix,suffix[i]),
				50, -4., 4.);
      hetaJet3[i][j]->SetDirectory(rootdir);
      hetaJet3[i][j]->GetXaxis()->SetTitle("#eta");
    
      hetaJet4[i][j] = new TH1F(Form("%s_hetaJet4_%s",prefix,suffix[i]),Form("%s_etaJet4_%s",prefix,suffix[i]),
				50, -4., 4.);
      hetaJet4[i][j]->SetDirectory(rootdir);
      hetaJet4[i][j]->GetXaxis()->SetTitle("#eta");
    
      heleSumPt[i][j] = new TH1F(Form("%s_heleSumPt_%s",prefix,suffix[i]),Form("%s_heleSumPt_%s",prefix,suffix[i]),
				 100, 0., 25.);
      heleSumPt[i][j]->SetDirectory(rootdir);
      heleSumPt[i][j]->GetXaxis()->SetTitle("#SigmaPt");
    
      hmuSumPt[i][j] = new TH1F(Form("%s_hmuSumPt_%s",prefix,suffix[i]),Form("%s_hmuSumPt_%s",prefix,suffix[i]),
				100, 0., 25.);
      hmuSumPt[i][j]->SetDirectory(rootdir);
      hmuSumPt[i][j]->GetXaxis()->SetTitle("#SigmaPt");
    
      hmuSumIso[i][j] = new TH1F(Form("%s_hmuIsoSum_%s",prefix,suffix[i]),Form("%s_hmuIsoSum_%s",prefix,suffix[i]),
				 100, 0., 25.);
      hmuSumIso[i][j]->SetDirectory(rootdir);
      hmuSumIso[i][j]->GetXaxis()->SetTitle("#SigmaPt");
    
      heleRelIso[i][j] = new TH1F(Form("%s_heleRelIso_%s",prefix,suffix[i]),Form("%s_heleRelIso_%s",prefix,suffix[i]),
				  100, 0., 1.0001);
      heleRelIso[i][j]->SetDirectory(rootdir);
      hmuRelIso[i][j] = new TH1F(Form("%s_hmuRelIso_%s",prefix,suffix[i]),Form("%s_hmuRelIso_%s",prefix,suffix[i]),
				 100, 0., 1.0001);

      if (j==0){
	hnJet[i]->Sumw2();
	hnJetLepVeto[i]->Sumw2();
      }
      numTightLep[i][j]->Sumw2();
      helePt[i][j]->Sumw2();
      hmuPt[i][j]->Sumw2();
      hmuPtFromSilicon[i][j]->Sumw2();
      hminLepPt[i][j]->Sumw2();
      hmaxLepPt[i][j]->Sumw2();
      helePhi[i][j]->Sumw2();
      hmuPhi[i][j]->Sumw2();
      hdphiLep[i][j]->Sumw2();
      heleEta[i][j]->Sumw2();
      hmuEta[i][j]->Sumw2();
      hdilMass[i][j]->Sumw2();
      hdilMassTightWindow[i][j]->Sumw2();
      hdilPt[i][j]->Sumw2();
      hmet[i][j]->Sumw2();
      hmetPhi[i][j]->Sumw2();
      hptJet1[i][j]->Sumw2();
      hptJet2[i][j]->Sumw2();
      hptJet3[i][j]->Sumw2();
      hptJet4[i][j]->Sumw2();
      hetaJet1[i][j]->Sumw2();
      hetaJet2[i][j]->Sumw2();
      hetaJet3[i][j]->Sumw2();
      hetaJet4[i][j]->Sumw2();
      heleSumPt[i][j]->Sumw2();
      hmuSumPt[i][j]->Sumw2();
      hmuSumIso[i][j]->Sumw2();
      heleRelIso[i][j]->Sumw2(); 
      hmuRelIso[i][j]->Sumw2(); 
    }
  }//channel loop
}//CMS2::bookHistos()
