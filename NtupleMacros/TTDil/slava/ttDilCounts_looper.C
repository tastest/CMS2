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


unsigned int cutConf(unsigned int cutsMask, unsigned int shift, unsigned int mask, std::string* sConf, bool doPrint){
  unsigned int res = ((cutsMask>> shift) & mask);
  if ( res > 3 ) { std::cout<<"Config error "<<std::endl; exit (99);}
  if (doPrint ) if (sConf[res].size() ) std::cout << sConf[res].c_str() << std::endl;
  //  std::cout<< cutsMask<<" "<<shift<<" "<<mask<<" "<<res<<" "<<sConf[res].c_str()<<std::endl;
  return (res);
}

#define DEFINE_CUT(CUT, OFFSET, MASK, C1, C2, C3, C4, S1, S2, S3, S4)	\
  unsigned int CUT##_##shift = OFFSET; unsigned int CUT##_##mask = MASK;\
  std::string CUT##_##confS[4] = { C1, C2, C3, C4 };			\
  std::string CUT##_##shortS[4] = { S1, S2, S3, S4 }

#define SET_CUT(BITS, CUT, SHORTS, FLAG)					\
  CUT = cutConf(BITS, CUT##_##shift, CUT##_##mask, CUT##_##confS, FLAG);\
  SHORTS += CUT##_##shortS[CUT] == "" ? "" : "_"+CUT##_##shortS[CUT]

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

  std::string compactConfig = "";

  // this will likely change (playing with options for now)...
  // To define, set and print the new cut:
  // declare the flag "bool flag = false;"
  // then "DEFINE_CUT(flag, bitIndex, mask, descriptionIfFalse, descrIfTrue, "", "", shortDescrIfFalse, shortDescrIfTrue)
  // then actually set the flag based on input cutsMask "SET_CUT(cutsMask, flag, compactConfig, printDescription)"
  // all these macro do is hide a bunch of copy-paste looking code from you
  // Eventually the DEFINE_CUT piece will move outside ScanChain to smth like SetConfig
  // and compactConfig will be usable as a part of the output file name
  // so that instead of currrent, e.g.,  myHist_2122752.root you get myHist_preDil08_OS_noDupWt_hltTry08.root
  bool idcuts = false;
  DEFINE_CUT(idcuts, 0, 1, 
	     "Id cuts disabled", "Id cuts enabled", "", "", 
	     "", "idOld", "", "");
  SET_CUT(cutsMask, idcuts, compactConfig, true);

  bool isolationcuts = false;
  DEFINE_CUT( isolationcuts, 1, 1,
	      "Isolation cuts disabled", "Isolation cuts enabled", "", "",
	      "", "isoOld", "", "");
  SET_CUT(cutsMask, isolationcuts, compactConfig, true);

  bool dilepMassVetoCut = false;
  DEFINE_CUT(dilepMassVetoCut, 2, 1,
	     "DiLeptonMassVetoCut disabled", "DiLeptonMassVetoCut enabled", "", "",
	     "", "zVetOld", "", "");
  SET_CUT(cutsMask, dilepMassVetoCut, compactConfig, true);

  bool METcut = false;
  DEFINE_CUT(METcut, 3, 1, 
	     "METCut disabled", "METCut enabled", "", "",
	     "", "wMET", "", "");
  SET_CUT(cutsMask, METcut, compactConfig, true);

  bool nJets2 = false; 
  DEFINE_CUT(nJets2, 4, 1,
	     "NJets>=2 cut disabled", "NJets>=2 cut enabled", "", "",
	     "", "2J", "", "");
  SET_CUT(cutsMask, nJets2, compactConfig, true);

  bool applyMuTag = false;
  DEFINE_CUT(applyMuTag, 5, 1, 
	     "Extra muon tag cut disabled", "Extra Muon tag cut enabled", "", "",
	     "", "muTag", "", "");
  SET_CUT(cutsMask, applyMuTag, compactConfig, true);

  bool METveto = false;
  DEFINE_CUT(METveto, 6, 1,
	     "MET veto is disabled", "MET veto is enabled", "", "",
	     "", "vetoMET", "", "");
  SET_CUT(cutsMask, METveto, compactConfig, true);

  bool applyMuTag5 = false;
  DEFINE_CUT(applyMuTag5, 7, 1,
	     "Extra muon 5GeV tag cut disabled", "Extra Muon 5GeV tag cut enabled", "", "",
	     "", "muTag5", "", "");
  SET_CUT(cutsMask, applyMuTag5, compactConfig, true);

  int  isoLooseMode = 0;
  DEFINE_CUT(isoLooseMode, 8, 3,
	     "Require both leptons to be isolated", 
	     "Require one-only hyp lepton (electron in emu) with iso (0.6, 0.92)", 
	     "Require one-only hyp lepton (muon in emu) with iso (0.6, 0.92)",
	     "Require two hyp leptons  with iso (0.6, 0.92)",
	     "", "isoReg1", "isoReg2", "isoReg3");
  SET_CUT(cutsMask, isoLooseMode, compactConfig, true);

  bool looseDilSelectionTTDil08 = false;
  DEFINE_CUT(looseDilSelectionTTDil08, 10, 1,
	     "", "Require loose dil selection for TTbar-->dilepton ana 2008/09", "", "",
	     "", "preDil08", "", "");
  SET_CUT(cutsMask, looseDilSelectionTTDil08, compactConfig, true);

  bool fillMultipleHypsOnly = false;
  DEFINE_CUT(fillMultipleHypsOnly, 11, 1,
	     "", "Fill only multiple hypotheses", "", "",
	     "", "dupOnly", "", "");
  SET_CUT(cutsMask, fillMultipleHypsOnly, compactConfig, true);

  bool applyZWindow = false;
  DEFINE_CUT(applyZWindow, 12, 1, 
	     "", "Events from Z-window only", "", "",
	     "", "inZ", "", "");
  SET_CUT(cutsMask, applyZWindow, compactConfig, true);

  bool osSelection = false;
  DEFINE_CUT(osSelection, 13, 1,
	     "", "Require OS", "", "",
	     "", "OS", "", "");
  SET_CUT(cutsMask,osSelection , compactConfig, true);

  bool fillMaxWeightDilOnly = false;
  DEFINE_CUT(fillMaxWeightDilOnly, 14, 1,
	     "", "Fill only the dilepton with the max dilepton weight", "", "",
	     "", "noDupWt", "", "");
  SET_CUT(cutsMask, fillMaxWeightDilOnly, compactConfig, true);

  bool leptonIsolationDilSelectionTTDil08 = false;
  DEFINE_CUT(leptonIsolationDilSelectionTTDil08, 15, 1,
	     "", "Apply isolation cuts on leptons for TTbar-->dilepton ana 2008/09", "", "",
	     "", "isoDil08", "", "");
  SET_CUT(cutsMask, leptonIsolationDilSelectionTTDil08, compactConfig, true);

  bool looseDilSelectionNoIsoTTDil08 = false;
  DEFINE_CUT(looseDilSelectionNoIsoTTDil08, 16, 1,
	     "", "Require loose dil selection for TTbar-->dilepton ana 2008/09; drop loose iso cuts", "", "",
	     "", "preDil08noIso", "", "");
  SET_CUT(cutsMask, looseDilSelectionNoIsoTTDil08, compactConfig, true);

  bool lepton20Eta2p4DilSelection = false;
  DEFINE_CUT(lepton20Eta2p4DilSelection, 17, 1,
	     "", "Two leptons pt>20 and |eta|<2.4 are selected -- bare minimum", "", "",
	     "", "2pt20", "", "");
  SET_CUT(cutsMask, lepton20Eta2p4DilSelection, compactConfig, true);

  bool metBaselineSelectionTTDil08 = false;
  DEFINE_CUT(metBaselineSelectionTTDil08, 18, 1,
	     "", "Apply TTDil08 baseline MET selection: use corrected pat-met emu met >20, mm,em met>30", "", "",
	     "", "preMet08", "", "");
  SET_CUT(cutsMask, metBaselineSelectionTTDil08, compactConfig, true);

  bool dilepMassVetoCutTTDil08 = false;
  DEFINE_CUT(dilepMassVetoCutTTDil08, 19, 1,
	     "", "Apply Z mass veto on same flavor dils, use TTDil08 selections", "", "",
	     "", "outZ08", "", "");
  SET_CUT(cutsMask, dilepMassVetoCutTTDil08, compactConfig, true);

  bool applyTriggersMu9orLisoE15 = false;
  DEFINE_CUT(applyTriggersMu9orLisoE15, 20, 1,
	     "", "HLT bits: 47 (HLT_LooseIsoEle15_LW_L1R), 82 (HLT_Mu9): ee -- 47, em -- 47 OR 82, mm -- 82", "", "",
	     "", "hltMu9E15", "", "");
  SET_CUT(cutsMask, applyTriggersMu9orLisoE15, compactConfig, true);

  bool applyTriggersTTDil08JanTrial = false;
  DEFINE_CUT(applyTriggersTTDil08JanTrial, 21, 1,
	     "", 
	     "HLT bits 45 (IsoEle18_L1R), 54 (DoubleIsoEle12_L1R), 86 (Mu15_L1Mu7), 90 (DoubleMu3),\
 126 (IsoEle10_Mu10_L1R):\n\t ee -- 45 OR 54, mm -- 86 or 90, em -- 45 OR 86 OR 126", "", "",
	     "", "hltTry08", "", "");
  SET_CUT(cutsMask, applyTriggersTTDil08JanTrial, compactConfig, true);

  bool dilepAdditionalMassVetoCutTTDil08 = false;
  DEFINE_CUT(dilepAdditionalMassVetoCutTTDil08, 22, 1,
	     "", "Apply additional z-veto (reject if there is a pair of loose,\
 at least one isolated same flavor OS leptons with mass inside z-window", "", "",
	     "", "extraZv", "", "");
  SET_CUT(cutsMask, dilepAdditionalMassVetoCutTTDil08, compactConfig, true);

  std::cout<<"Compact config string is "<<compactConfig.c_str()<<std::endl;

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


      for(unsigned int hypIdx = 0; hypIdx < hyp_p4().size(); hypIdx++) {
       
	unsigned int i_lt = hyp_lt_index().at(hypIdx);
	unsigned int i_ll = hyp_ll_index().at(hypIdx);

	int id_lt = hyp_lt_id().at(hypIdx);
	int id_ll = hyp_ll_id().at(hypIdx);

	{// scope out old/legacy selections
	  if (applyZWindow && fabs(hyp_p4().at(hypIdx).mass()-91)> 15) continue;
	  
	  if(dilepMassVetoCut) {
	    // Z mass veto using hyp_leptons for ee and mumu final states
	    if (hyp_type().at(hypIdx) == 0 || hyp_type().at(hypIdx) == 3) {
	      if (inZmassWindow(hyp_p4().at(hypIdx).mass())) continue;
	    }
	    
	    // Z veto using additional leptons in the event
	    if (additionalZveto()) continue;
	  }
	  
	  // Dima's MET requirement
	  bool pass2MetPassed = pass2Met(hypIdx);
	  if (! pass2MetPassed && METcut) continue;
	  if ( pass2MetPassed && METveto ) continue;

	  if(idcuts) {
	    // Muon quality cuts, excluding isolation
	    if (abs(id_lt) == 13 && !goodMuonWithoutIsolation(i_lt) ) continue;
	    if (abs(id_ll) == 13 && !goodMuonWithoutIsolation(i_ll) ) continue;
      
	    // Electron quality cuts, excluding isolation
	    if (abs(id_lt) == 11 && !goodElectronWithoutIsolation(i_lt) ) continue;
	    if (abs(id_ll) == 11 && !goodElectronWithoutIsolation(i_ll) ) continue;
	  }

	  if(isolationcuts) {
	    if (!passDilAntiIsolation(isoLooseMode, hypIdx)) continue; 
	  }
    
  
	  if (applyMuTag && ! haveExtraMuon(hypIdx)) continue;
	  if (applyMuTag5 && ! haveExtraMuon5(hypIdx)) continue;

	}// end scope-out of old selections

	// this is for per-hypothesis choice
	if (! fillMaxWeightDilOnly && applyTriggersMu9orLisoE15){
	  if (! passTriggersMu9orLisoE15(hyp_type().at(hypIdx)) ) continue;
	}

	if (! fillMaxWeightDilOnly && applyTriggersTTDil08JanTrial){
	  if (! passTriggersTTDil08JanTrial(hyp_type().at(hypIdx)) ) continue;
	}

        if(dilepMassVetoCutTTDil08) {
          // Z mass veto using hyp_leptons for ee and mumu final states
          if (hyp_type().at(hypIdx) == 0 || hyp_type().at(hypIdx) == 3) {
            if (inZmassWindow(hyp_p4().at(hypIdx).mass())) continue;
          }    
        }
	if(dilepAdditionalMassVetoCutTTDil08){
          // Z veto using additional leptons in the event
          if (additionalZveto(true)) continue; //"true" to use TTDil lepton selections                                                
	}
      
	// ! for TTDil analysis this should be made for the event-qualifying hyp only
	if (!fillMaxWeightDilOnly && metBaselineSelectionTTDil08){
	  if (! passPatMet_OF20_SF30(hypIdx)) continue;
	}

	if (lepton20Eta2p4DilSelection){
	  //pt eta cuts
	  if (! lepton20Eta2p4(id_lt, i_lt) ) continue;
	  if (! lepton20Eta2p4(id_ll, i_ll) ) continue;
	}

	if (looseDilSelectionNoIsoTTDil08){
	  // Muon quality cuts, no isolation
	  if (! looseLeptonSelectionNoIsoTTDil08(id_lt, i_lt)) continue;
	  if (! looseLeptonSelectionNoIsoTTDil08(id_ll, i_ll)) continue;
	}

	if (looseDilSelectionTTDil08){
	  // Muon quality cuts, loose isolation
	  if (! looseLeptonSelectionTTDil08(id_lt, i_lt)) continue;
	  if (! looseLeptonSelectionTTDil08(id_ll, i_ll)) continue;
	}

	if (leptonIsolationDilSelectionTTDil08){
	  if (! passLeptonIsolationTTDil08(id_lt, i_lt)) continue;
	  if (! passLeptonIsolationTTDil08(id_ll, i_ll)) continue;
	}

	if (osSelection){
	  if ( id_lt * id_ll > 0 ) continue;
	}

	goodHyps.push_back(hypIdx);
	//done with cuts on hyps		   

	eventPassed = true;
      }
      
      unsigned int nGoodHyps = goodHyps.size();

      unsigned int maxWeightIndex = 0;
      int strasbourgDilType = -1;

      if (nGoodHyps > 0){
	maxWeightIndex = eventDilIndexByWeightTTDil08(goodHyps, strasbourgDilType);

	// ! event level cut here, can reset the eventPassed to false
	if (fillMaxWeightDilOnly && metBaselineSelectionTTDil08){
	  if (! passPatMet_OF20_SF30(maxWeightIndex)) continue;
	}

	// this is for per-hypothesis choice
	if ( fillMaxWeightDilOnly && applyTriggersMu9orLisoE15){
	  if (! passTriggersMu9orLisoE15(hyp_type().at(maxWeightIndex)) ) continue;
	}

	if ( fillMaxWeightDilOnly && applyTriggersTTDil08JanTrial){
	  if (! passTriggersTTDil08JanTrial(hyp_type().at(maxWeightIndex)) ) continue;
	}
      }

      //=============================================================================================

      //now fill the histograms

      //=============================================================================================

      for(unsigned int hypIdxL=0; hypIdxL< nGoodHyps; ++ hypIdxL){
	unsigned int hypIdx = goodHyps[hypIdxL];
	if (fillMaxWeightDilOnly && hypIdx != maxWeightIndex) continue;


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

	float pt_lt = hyp_lt_p4().at(hypIdx).pt();
	float pt_ll = hyp_ll_p4().at(hypIdx).pt();

	unsigned int i_lt = hyp_lt_index().at(hypIdx);
	unsigned int i_ll = hyp_ll_index().at(hypIdx);

	int id_lt = hyp_lt_id().at(hypIdx);
	int id_ll = hyp_ll_id().at(hypIdx);


	// jet count
	hnJet[myType]->Fill(new_hyp_njets, weight);
	hnJet[3]->Fill(new_hyp_njets, weight);

	// lepton Pt
	if (abs(id_lt) == 11) helePt[myType][arrNjets]->Fill(pt_lt, weight);
	if (abs(id_ll) == 11) helePt[myType][arrNjets]->Fill(pt_ll, weight);
	if (abs(id_lt) == 13) hmuPt[myType][arrNjets]->Fill(pt_lt, weight);
	if (abs(id_ll) == 13) hmuPt[myType][arrNjets]->Fill(pt_ll, weight);
	if (abs(id_lt) == 13) 
	  hmuPtFromSilicon[myType][arrNjets]->Fill(mus_trk_p4().at(i_lt).pt(), weight);
	if (abs(id_ll) == 13)
	  hmuPtFromSilicon[myType][arrNjets]->Fill(mus_trk_p4().at(i_ll).pt(), weight);
	hminLepPt[myType][arrNjets]->Fill(min(pt_ll, pt_lt), weight);
	hmaxLepPt[myType][arrNjets]->Fill(max(pt_ll, pt_lt), weight );
    
	if (abs(id_lt) == 11) helePt[3][arrNjets]->Fill(pt_lt, weight);
	if (abs(id_ll) == 11) helePt[3][arrNjets]->Fill(pt_ll, weight);
	if (abs(id_lt) == 13) hmuPt[3][arrNjets]->Fill(pt_lt, weight);
	if (abs(id_ll) == 13) hmuPt[3][arrNjets]->Fill(pt_ll, weight);
	if (abs(id_lt) == 13) 
	  hmuPtFromSilicon[3][arrNjets]->Fill(mus_trk_p4().at(i_lt).pt(), weight);
	if (abs(id_ll) == 13) 
	  hmuPtFromSilicon[3][arrNjets]->Fill(mus_trk_p4().at(i_ll).pt(), weight);
	hminLepPt[3][arrNjets]->Fill(min(pt_ll, pt_lt), weight);
	hmaxLepPt[3][arrNjets]->Fill(max(pt_ll, pt_lt), weight );


	// lepton Phi
	if (abs(id_lt) == 11) helePhi[myType][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).phi(), weight);
	if (abs(id_ll) == 11) helePhi[myType][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).phi(), weight);
	if (abs(id_lt) == 13) hmuPhi[myType][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).phi(), weight);
	if (abs(id_ll) == 13) hmuPhi[myType][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).phi(), weight);
    
	if (abs(id_lt) == 11) helePhi[3][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).phi(), weight);
	if (abs(id_ll) == 11) helePhi[3][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).phi(), weight);
	if (abs(id_lt) == 13) hmuPhi[3][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).phi(), weight);
	if (abs(id_ll) == 13) hmuPhi[3][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).phi(), weight);
    
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
	if (abs(id_lt) == 11) heleEta[myType][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).eta(), weight);
	if (abs(id_ll) == 11) heleEta[myType][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).eta(), weight);
	if (abs(id_lt) == 13) hmuEta[myType][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).eta(), weight);
	if (abs(id_ll) == 13) hmuEta[myType][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).eta(), weight);
    
	if (abs(id_lt) == 11) heleEta[3][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).eta(), weight);
	if (abs(id_ll) == 11) heleEta[3][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).eta(), weight);
	if (abs(id_lt) == 13) hmuEta[3][arrNjets]->Fill(hyp_lt_p4().at(hypIdx).eta(), weight);
	if (abs(id_ll) == 13) hmuEta[3][arrNjets]->Fill(hyp_ll_p4().at(hypIdx).eta(), weight);
       
      
	// electron trk isolation 
	double temp_lt_iso = hyp_lt_iso().at(hypIdx);  // so that min works
	double temp_ll_iso = hyp_ll_iso().at(hypIdx);  // so that min works
	if (abs(id_lt) == 11) heleSumPt[myType][arrNjets]->Fill(min(temp_lt_iso,24.99),weight);
	if (abs(id_lt) == 11) heleSumPt[3][arrNjets]->Fill(min(temp_lt_iso,24.99),weight);
	if (abs(id_ll) == 11) heleSumPt[myType][arrNjets]->Fill(min(temp_ll_iso,24.99),weight);
	if (abs(id_ll) == 11) heleSumPt[3][arrNjets]->Fill(min(temp_ll_iso,24.99),weight);

	// muon trk isolation
	if (abs(id_lt) == 13) hmuSumPt[myType][arrNjets]->Fill(min(temp_lt_iso,24.99),weight);
	if (abs(id_lt) == 13) hmuSumPt[3][arrNjets]->Fill(min(temp_lt_iso,24.99),weight);
	if (abs(id_ll) == 13) hmuSumPt[myType][arrNjets]->Fill(min(temp_ll_iso,24.99),weight);
	if (abs(id_ll) == 13) hmuSumPt[3][arrNjets]->Fill(min(temp_ll_iso,24.99),weight);
      
      
	// muon trk+calo isolation
	double combIso_lt = -1.;
	double combIso_ll = -1.;
	if (abs(id_lt) == 13)
	  combIso_lt = mus_iso03_sumPt().at(i_lt)
	    +mus_iso03_emEt().at(i_lt)
	    +mus_iso03_hadEt().at(i_lt);
	if (abs(id_ll) == 13)
	  combIso_ll = mus_iso03_sumPt().at(i_ll)
	    +mus_iso03_emEt().at(i_ll)
	    +mus_iso03_hadEt().at(i_ll);
	if (abs(id_lt) == 13) hmuSumIso[myType][arrNjets]->Fill(min(combIso_lt,24.99),weight);
	if (abs(id_lt) == 13) hmuSumIso[3][arrNjets]->Fill(min(combIso_lt,24.99),weight);
	if (abs(id_ll) == 13) hmuSumIso[myType][arrNjets]->Fill(min(combIso_ll,24.99),weight);
	if (abs(id_ll) == 13) hmuSumIso[3][arrNjets]->Fill(min(combIso_ll,24.99),weight);


	// Relative isolation... muons
	if (abs(id_lt) == 13) {
	  double thisSum =  mus_iso03_sumPt().at(i_lt) +  
	    mus_iso03_emEt().at(i_lt)  +
	    mus_iso03_hadEt().at(i_lt);
	  double thisPt  = mus_p4().at(i_lt).pt();
	  double temp    = thisPt / (thisPt+thisSum);
	  hmuRelIso[myType][arrNjets]->Fill(temp, weight);
	  hmuRelIso[3][arrNjets]->Fill(temp, weight);
	}
	if (abs(id_ll) == 13) {
	  double thisSum =  mus_iso03_sumPt().at(i_ll) +  
	    mus_iso03_emEt().at(i_ll)  +
	    mus_iso03_hadEt().at(i_ll);
	  double thisPt  = mus_p4().at(i_ll).pt();
	  double temp    = thisPt / (thisPt+thisSum);
	  hmuRelIso[myType][arrNjets]->Fill(temp, weight);
	  hmuRelIso[3][arrNjets]->Fill(temp, weight);
	}


	// Relative isolation... electrons
	if (abs(id_lt) == 11) {
	  double thisSum =  hyp_lt_iso().at(hypIdx);
	  double thisPt  = pt_lt;
	  double temp    = thisPt / (thisPt+thisSum);
	  heleRelIso[myType][arrNjets]->Fill(temp, weight);
	  heleRelIso[3][arrNjets]->Fill(temp, weight);
	}
	if (abs(id_ll) == 11) {
	  double thisSum =  hyp_ll_iso().at(hypIdx);
	  double thisPt  = pt_ll;
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
	
	for(int k = 0; k<5; k++) {
	  hnJet[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
	  hnJet[i]->GetXaxis()->SetLabelSize(0.07);
	  
	}
      }
    
      char *suffix[4];
      char *njetCh[3] = { "0j", "1j", "2j" };
      suffix[0] = Form("%s_ee", njetCh[j]);
      suffix[1] = Form("%s_mm", njetCh[j]);
      suffix[2] = Form("%s_em", njetCh[j]);
      suffix[3] = Form("%s_all", njetCh[j]);

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
      }
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
