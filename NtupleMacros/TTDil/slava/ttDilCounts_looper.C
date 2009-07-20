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

#include "CORE/CMS2.cc"
#include "ttDilCounts_looper.h"
#include "CORE/selections.cc"
#include "CORE/utilities.cc"

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

  compactConfig = "";

  // this will likely change (playing with options for now)...
  // To define, set and print the new cut:
  // declare the flag "bool flag = false;"
  // then "DEFINE_CUT(flag, bitIndex, mask, descriptionIfFalse, descrIfTrue, "", "", shortDescrIfFalse, shortDescrIfTrue)
  // then actually set the flag based on input cutsMask "SET_CUT(cutsMask, flag, compactConfig, printDescription)"
  // all these macro do is hide a bunch of copy-paste looking code from you
  // Eventually the DEFINE_CUT piece will move outside ScanChain to smth like SetConfig
  // and compactConfig will be usable as a part of the output file name
  // so that instead of currrent, e.g.,  myHist_2122752.root you get myHist_preDil08_OS_noDupWt_hltTry08.root
  bool idcuts =  cutsMask & 1;
  if( cutsMask & 1) {
    idcuts = true;
    cout << "Id cuts enabled" << endl;
    compactConfig = compactConfig+ "_idOld";
  }

  bool isolationcuts = (cutsMask>>1) & 1;
  if( isolationcuts ) {
    cout << "Isolation cuts enabled" << endl;
    compactConfig = compactConfig + "_isoOld";
  }

  bool dilepMassVetoCut =  (cutsMask>>2) & 1;
  if( dilepMassVetoCut ) {
    cout << "DiLeptonMassVetoCut enabled" << endl;
    compactConfig = compactConfig + "_zVetoOld";
  }
  
  bool METcut =  (cutsMask>>3) & 1;
   if( METcut ) {
     cout << "METCut enabled" << endl;
    compactConfig = compactConfig + "_wMET";
  } 

   bool nJets2 =  (cutsMask>>4) & 1; 
   if( nJets2 ) {
     cout << "NJets>=2 cut enabled" << endl;
     compactConfig = compactConfig + "_nJets2";
   }

   bool applyMuTag = (cutsMask>>5) & 1;
   if( applyMuTag ) {
     cout << "Extra Muon tag cut enabled" << endl;
     compactConfig = compactConfig + "_muTag";
   } 
   
   bool METveto = (cutsMask>>6) & 1;
   if ( METveto ){
     cout << "MET veto is enabled" << endl;
     compactConfig = compactConfig + "_vetoMET";
   } 

   bool applyMuTag5 =  (cutsMask>>7) & 1;
   if( applyMuTag5 ) {
     cout << "Extra Muon 5GeV tag cut enabled" << endl;
     compactConfig = compactConfig + "_muTag5";
   } 
   
   int  isoLooseMode = 0;
   isoLooseMode = ((cutsMask>>8) & 3);
   if (isoLooseMode == 0){
     cout << ".............." << endl;
   } else {
     std::cout<< " This is a dummy setting: don't turn it on"<<std::endl;
   }
   
   bool looseDilSelectionTTDil08 = (cutsMask>>10) & 1;
   if (looseDilSelectionTTDil08 ){
     cout << "Require loose dil selection for TTbar-->dilepton ana 2008/09"<< endl;
     compactConfig = compactConfig + "_looseDil08";
   }
   
   bool fillMultipleHypsOnly = ((cutsMask >> 11)&1);
   if (fillMultipleHypsOnly){
     std::cout<< "Fill only multiple hypotheses"<<std::endl;
     compactConfig = compactConfig + "_dupOnly";
   }

   bool applyZWindow = ((cutsMask >> 12) & 1);
   if (applyZWindow){
     std::cout<<"Events from Z-window only"<<std::endl;
     compactConfig = compactConfig + "_inZ";
   }
   

   bool osSelection = ((cutsMask >> 13 ) & 1);
   if (osSelection){
     std::cout<<"Require OS"<<std::endl;
     compactConfig = compactConfig + "_OS";
   }
   
   bool  fillMaxWeightDilOnly = ((cutsMask >> 14) & 1);
   if (fillMaxWeightDilOnly){
     std::cout<<"Fill only the dilepton with the max dilepton weight"<<std::endl;
     compactConfig = compactConfig + "_noDupWt";
  }
   
   bool leptonIsolationDilSelectionTTDil08 = ((cutsMask>> 15) & 1);
   if (leptonIsolationDilSelectionTTDil08){
     std::cout<<" Apply isolation cuts on leptons for TTbar-->dilepton ana 2008/09"<<std::endl;
     compactConfig = compactConfig + "_isoDil08";
   }
   
   bool looseDilSelectionNoIsoTTDil08 = ((cutsMask>>16)&1);
   if (looseDilSelectionNoIsoTTDil08 ){
     std::cout<< "Require loose dil selection for TTbar-->dilepton ana 2008/09; drop loose iso cuts"<<std::endl;
     compactConfig = compactConfig +"_preDil08noIso";
   }
   
   bool lepton20Eta2p4DilSelection = ((cutsMask>>17)&1);
   if (lepton20Eta2p4DilSelection){
     std::cout<<"Two leptons pt>20 and |eta|<2.4 are selected -- bare minimum"<<std::endl;
     compactConfig = compactConfig;
   }

   bool useTcMet = ((cutsMask>>25)&1);

   bool metBaselineSelectionTTDil08 = ((cutsMask>>18)&1);
   if (metBaselineSelectionTTDil08){
     if (useTcMet) {
       std::cout<<"Apply TTDil08 baseline MET selection: use tcMET emu met >20, mm,ee met>30"<<std::endl;
       compactConfig = compactConfig + "_preTcMet08";
     } else {
       std::cout<<"Apply TTDil08 baseline MET selection: use corrected pat-met emu met >20, mm,ee met>30"<<std::endl;
       compactConfig = compactConfig + "_preMet08";
     }
   }
   
   // careful: tcmet switch is only allowed for metBaselineSelectionTTDil08
   if (useTcMet) {
     if (!metBaselineSelectionTTDil08){
       std::cout<<" ***** tcmet is only allowed in conjunction with the TTDil08 baseline MET selection *** STOP"<<std::endl;
       return 1;
     }
   }

   bool dilepMassVetoCutTTDil08 = ((cutsMask>>19)&1);
   if (dilepMassVetoCutTTDil08){
     std::cout<<"Apply Z mass veto on same flavor dils, use TTDil08 selections"<<std::endl;
     compactConfig = compactConfig + "_outZ08";
   }
   
   bool applyTriggersMu9orLisoE15 = ((cutsMask>>20)&1);
   if (applyTriggersMu9orLisoE15){
     std::cout<<"HLT bits: 47 (HLT_LooseIsoEle15_LW_L1R), 82 (HLT_Mu9): ee -- 47, em -- 47 OR 82, mm -- 82"<<std::endl;
     compactConfig = compactConfig + "_hltMu9E15";
   }

   bool applyTriggersTTDil08JanTrial = ((cutsMask>>21)&1);
   if (applyTriggersTTDil08JanTrial){
     std::cout<<"HLT bits 45 (IsoEle18_L1R), 54 (DoubleIsoEle12_L1R), 86 (Mu15_L1Mu7), 90 (DoubleMu3), 126 (IsoEle10_Mu10_L1R):\n"
	      <<"\t ee -- 45 OR 54, mm -- 86 or 90, em -- 45 OR 86 OR 126"<<std::endl;
     compactConfig = compactConfig + "_hltTry08";
   }
   
   bool dilepAdditionalMassVetoCutTTDil08 = ((cutsMask>>22)&1);
   if( dilepAdditionalMassVetoCutTTDil08 ) {
     cout << "Apply additional z-veto. Reject event if there is a pair of opp. "
         << " sign same flavor leptons (with one loosely isolated) that has inv. mass " 
	  << "inside Z mass veto region" << endl;
     compactConfig = compactConfig + "_xtraZv";
   }

   bool corJES10ptUp = ((cutsMask>>23)&1);
   if(corJES10ptUp) {
     cout << "Jets are scaled 10% up" << endl;
     compactConfig = compactConfig + "_jets10Uphyp";
  }
   
   bool corJES10ptDn = ((cutsMask>>24)&1);
   if(corJES10ptDn) {
     cout << "Jets are scaled 10% down" << endl;
     compactConfig = compactConfig + "_jets10Dnhyp";
   }
   
   if (corJES10ptUp && corJES10ptDn){
     std::cout<<"Inconsistent config: JES up and down requested at the same time: bailing"<<std::endl;
     return 99;
   }

   if ((corJES10ptUp || corJES10ptDn) && useTcMet) {
     std::cout<<"*********************************************************************" <<std::endl;
     std::cout<<"CAUTION: You are rescaling the jets by 10%, and you want to use tcMet" <<std::endl;
     std::cout<<"         The met will NOT be rescaled in any way"                      <<std::endl;
     std::cout<<"*********************************************************************" <<std::endl;
   }

  float globalJESscaleRescale = 1.;
  if (corJES10ptUp)  globalJESscaleRescale = 1.1;
  if (corJES10ptDn)  globalJESscaleRescale = 0.9;

  bool useJPT = ((cutsMask>>26)&1);
  if (useJPT) {
    if (oldjets) {
      std::cout<<"Inconsistent config: JPT and oldjets requested: bailing"<<std::endl;
      return 99;
    }
    std::cout<<"Using JPT jets"<<std::endl;
    compactConfig = compactConfig + "_JPT";
  }

  bool muJetClean = ((cutsMask>>27)&1);
  if (muJetClean) {
    std::cout<<"Cleaning jets near muons"<<std::endl;
    compactConfig = compactConfig + "_JmuCl";
  }

  bool dilTruthMatch = ((cutsMask>>28)&1);
  if (dilTruthMatch){
    std::cout<<"Dileptons are required to match to W/Z mother"<<std::endl;
    compactConfig = compactConfig + "_MCtruth";
  }

  bool dilWeightMaxMass = ((cutsMask>>29)&1);
  if (dilWeightMaxMass){
    std::cout<<"Will order dileptons by highest mass and not by pt*iso weight"<<std::endl;
    compactConfig = compactConfig + "_maxMass";
  }
  bool dilWeightMaxPt = ((cutsMask>>30)&1);
  if (dilWeightMaxPt){
    std::cout<<"Will order dileptons by highest pt"<<std::endl;
    compactConfig = compactConfig + "_maxPt";
  }

  if (dilWeightMaxPt && dilWeightMaxMass){
    std::cout<<"Inconsistent configuration of dilWeightMaxPt && dilWeightMaxMass: can not both be true"<<std::endl;
    return 0;
  }

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
    cms2.Init(tree);
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
      cms2.GetEntry(z);
      ++nEventsTotal;


      //check if it's a correct genp-event
      std::string prefixStr(prefix);
      if (prefixStr == "ttdil" && genpCountPDGId(11,13,15) != 2) continue;
      if (prefixStr == "ttotr" && genpCountPDGId(11,13,15) == 2) continue;
      if (prefixStr == "DYeemm" && genpCountPDGId(11) != 2 && genpCountPDGId(13) != 2) continue;
      if (prefixStr == "DYee" && genpCountPDGId(11) != 2) continue;
      if (prefixStr == "DYmm" && genpCountPDGId(13) != 2) continue;
      if (prefixStr == "DYtautau" && genpCountPDGId(15) != 2) continue;

      //decide weather or not the event passed
      bool eventPassed = false;
      

      std::vector<unsigned int> goodHyps(0);


      for(unsigned int hypIdx = 0; hypIdx < cms2.hyp_p4().size(); hypIdx++) {
       
	unsigned int i_lt = cms2.hyp_lt_index()[hypIdx];
	unsigned int i_ll = cms2.hyp_ll_index()[hypIdx];

	int id_lt = cms2.hyp_lt_id()[hypIdx];
	int id_ll = cms2.hyp_ll_id()[hypIdx];

	{// scope out old/legacy selections
	  if (applyZWindow && fabs(cms2.hyp_p4()[hypIdx].mass()-91)> 15) continue;
	  
	  if(dilepMassVetoCut) {
	    // Z mass veto using hyp_leptons for ee and mumu final states
	    if (cms2.hyp_type()[hypIdx] == 0 || cms2.hyp_type()[hypIdx] == 3) {
	      if (inZmassWindow(cms2.hyp_p4()[hypIdx].mass())) continue;
	    }
	    
	    // Z veto using additional leptons in the event
	    if (additionalZveto()) continue;
	  }
	  
	  // Dima's MET requirement
	  TVector3 vecDummy;
	  bool pass2MetPassed = pass2Met(hypIdx, vecDummy);
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
	    if (!passLeptonIsolation(id_ll, i_ll, true)) continue; //use caloiso if it's ele
	    if (!passLeptonIsolation(id_lt, i_lt, true)) continue; //use caloiso it it's ele
	  }
    
  
	  if (applyMuTag && ! haveExtraMuon(hypIdx)) continue;
	  if (applyMuTag5 && ! haveExtraMuon5(hypIdx)) continue;

	}// end scope-out of old selections

	// this is for per-hypothesis choice
	if (! fillMaxWeightDilOnly && applyTriggersMu9orLisoE15){
	  if (! passTriggersMu9orLisoE15(cms2.hyp_type()[hypIdx]) ) continue;
	}

	if (! fillMaxWeightDilOnly && applyTriggersTTDil08JanTrial){
	  if (! passTriggersTTDil08JanTrial(cms2.hyp_type()[hypIdx]) ) continue;
	}

        if(dilepMassVetoCutTTDil08) {
          // Z mass veto using hyp_leptons for ee and mumu final states
          if (cms2.hyp_type()[hypIdx] == 0 || cms2.hyp_type()[hypIdx] == 3) {
            if (inZmassWindow(cms2.hyp_p4()[hypIdx].mass())) continue;
          }    
        }
	if(dilepAdditionalMassVetoCutTTDil08){
          // Z veto using additional leptons in the event
          if (additionalZvetoTTDil08()) continue; //"true" to use TTDil lepton selections                                                
	}
      
	// ! for TTDil analysis this should be made for the event-qualifying hyp only
	if (!fillMaxWeightDilOnly && metBaselineSelectionTTDil08){
	  if (globalJESscaleRescale == 1 && useTcMet && ! passMet_OF20_SF30(hypIdx,useTcMet)) continue;
	  if (globalJESscaleRescale != 1. && (!useTcMet)) {
	    float metx = cms2.met_pat_metCor()*cos(cms2.met_pat_metPhiCor());
	    float mety = cms2.met_pat_metCor()*sin(cms2.met_pat_metPhiCor());

	    unsigned int nJused = 0;
	    unsigned int nJ = cms2.hyp_jets_p4()[hypIdx].size();
	    for (unsigned int iJ = 0; iJ < nJ; ++iJ){
	      if (cms2.hyp_jets_p4()[hypIdx][iJ].pt() > 0 
		  && fabs( cms2.hyp_jets_p4()[hypIdx][iJ].eta()) < 12.4){
		metx -= cms2.hyp_jets_p4()[hypIdx][iJ].x()*(globalJESscaleRescale - 1.); 
		mety -= cms2.hyp_jets_p4()[hypIdx][iJ].y()*(globalJESscaleRescale - 1.); 
		nJused++;
	      }
	    }
	    unsigned int nOJ = cms2.hyp_other_jets_p4()[hypIdx].size();
	    for (unsigned int iJ = 0; iJ < nOJ; ++iJ){
	      if (cms2.hyp_other_jets_p4()[hypIdx][iJ].pt() > 0 
		  && fabs( cms2.hyp_other_jets_p4()[hypIdx][iJ].eta()) < 12.4){
		metx -= cms2.hyp_other_jets_p4()[hypIdx][iJ].x()*(globalJESscaleRescale - 1.); 
		mety -= cms2.hyp_other_jets_p4()[hypIdx][iJ].y()*(globalJESscaleRescale - 1.); 
		nJused++;
	      }
	    }
	    if (! passPatMet_OF20_SF30(metx, mety, hypIdx)) continue;
	  }
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

	if (dilTruthMatch){
	  //this better be in the selections.cc
	  bool isTrueLepton_ll = false;
	  bool isTrueLepton_lt = false;
	  isTrueLepton_ll = ( (abs(cms2.hyp_ll_id()[hypIdx]) == abs(cms2.hyp_ll_mc_id()[hypIdx]) &&
			       abs(cms2.hyp_ll_mc_motherid()[hypIdx]) < 50 //I wish I could match to W or Z explicitely, not in MGraph
			       )
			      || (cms2.hyp_ll_mc_id()[hypIdx]==22 && 
				  TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[hypIdx],cms2.hyp_ll_mc_p4()[hypIdx])) <0.05
				  && abs(cms2.hyp_ll_id()[hypIdx]) == abs(cms2.hyp_ll_mc_motherid()[hypIdx])
				  )
			      );
	  isTrueLepton_lt = ( (abs(cms2.hyp_lt_id()[hypIdx]) == abs(cms2.hyp_lt_mc_id()[hypIdx]) &&
			       abs(cms2.hyp_lt_mc_motherid()[hypIdx]) < 50 //I wish I could match to W or Z explicitely, not in MGraph
			       )
			      || (cms2.hyp_lt_mc_id()[hypIdx]==22 && 
				  TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[hypIdx],cms2.hyp_lt_mc_p4()[hypIdx])) <0.05
				  && abs(cms2.hyp_lt_id()[hypIdx]) == abs(cms2.hyp_lt_mc_motherid()[hypIdx])
				  )
			      );
	  if (!isTrueLepton_lt && !isTrueLepton_ll) continue;
	}

	goodHyps.push_back(hypIdx);
	//done with cuts on hyps		   

	eventPassed = true;
      }
      
      unsigned int nGoodHyps = goodHyps.size();

      unsigned int maxWeightIndex = 0;
      int strasbourgDilType = -1;

      if (nGoodHyps > 0){
	bool debugPrintDispatch = (prefixStr == "ttdil" || prefixStr == "ttotr");
	if (dilWeightMaxMass){
	  maxWeightIndex = eventDilIndexByMaxMass(goodHyps, debugPrintDispatch);
	} else if (dilWeightMaxPt) {
	  bool usePtOnlyForWeighting = true;
	  maxWeightIndex = eventDilIndexByWeightTTDil08(goodHyps, strasbourgDilType, debugPrintDispatch, usePtOnlyForWeighting);
	} else {
	  bool usePtOnlyForWeighting = false;
	  maxWeightIndex = eventDilIndexByWeightTTDil08(goodHyps, strasbourgDilType, debugPrintDispatch, usePtOnlyForWeighting);
	}

	// ! event level cut here, can reset the eventPassed to false
	if (fillMaxWeightDilOnly && metBaselineSelectionTTDil08){
	  if (globalJESscaleRescale == 1. && useTcMet && ! passMet_OF20_SF30(maxWeightIndex,useTcMet)) continue;
	  if ( (!useTcMet)){
	    float metx = cms2.met_pat_metCor()*cos(cms2.met_pat_metPhiCor());
	    float mety = cms2.met_pat_metCor()*sin(cms2.met_pat_metPhiCor());

	    unsigned int nJ = cms2.hyp_jets_p4()[maxWeightIndex].size();
	    for (unsigned int iJ = 0; iJ < nJ; ++iJ){
	      if (cms2.hyp_jets_p4()[maxWeightIndex][iJ].pt() > 0  
		  && fabs( cms2.hyp_jets_p4()[maxWeightIndex][iJ].eta()) < 12.4){ 
		metx -= cms2.hyp_jets_p4()[maxWeightIndex][iJ].x()*(globalJESscaleRescale - 1.); 
		mety -= cms2.hyp_jets_p4()[maxWeightIndex][iJ].y()*(globalJESscaleRescale - 1.); 
	      }
	    }
	    unsigned int nOJ = cms2.hyp_other_jets_p4()[maxWeightIndex].size();
	    for (unsigned int iJ = 0; iJ < nOJ; ++iJ){
	      if (cms2.hyp_other_jets_p4()[maxWeightIndex][iJ].pt() > 0  
		  && fabs( cms2.hyp_other_jets_p4()[maxWeightIndex][iJ].eta()) < 12.4){ 
		metx -= cms2.hyp_other_jets_p4()[maxWeightIndex][iJ].x()*(globalJESscaleRescale - 1.); 
		mety -= cms2.hyp_other_jets_p4()[maxWeightIndex][iJ].y()*(globalJESscaleRescale - 1.); 
	      }
	    }
	    if (! passPatMet_OF20_SF30(metx, mety, maxWeightIndex)) continue;
	  }
	}

	// this is for per-hypothesis choice
	if ( fillMaxWeightDilOnly && applyTriggersMu9orLisoE15){
	  if (! passTriggersMu9orLisoE15(cms2.hyp_type()[maxWeightIndex]) ) continue;
	}

	if ( fillMaxWeightDilOnly && applyTriggersTTDil08JanTrial){
	  if (! passTriggersTTDil08JanTrial(cms2.hyp_type()[maxWeightIndex]) ) continue;
	}
      }

      //=============================================================================================

      //now fill the histograms

      //=============================================================================================

      for(unsigned int hypIdxL=0; hypIdxL< nGoodHyps; ++ hypIdxL){
	unsigned int hypIdx = goodHyps[hypIdxL];
	if (fillMaxWeightDilOnly && hypIdx != maxWeightIndex) continue;


	// The event weight including the kFactor (scaled to 1 fb-1) and the prescale
	//float weight = cms2.evt_scale1fb * kFactor * prescale;
	//float weight = CalculateWeight(evt_CSA07Process(), cms2.evt_scale1fb(), kFactor, prescale);
	float weight = kFactor*cms2.evt_scale1fb()*0.01; // /100; //10pb^-1

	if ( (prefixStr == "ppMuX" || prefixStr == "EM" || prefixStr == "QCD") 
	     && (cms2.hyp_type()[hypIdx] == 1 || cms2.hyp_type()[hypIdx] == 2)) weight *= 0.5;//this isn't quite right :(
	//and works only if both em and ppmux are in play
      
	// If we made it to here, we passed all cuts and we are ready to fill
	m_events.insert(pair<int,int>(cms2.evt_event(), 1));

	int myType = 99;
	if (cms2.hyp_type()[hypIdx] == 3) myType = 0;  // ee
	if (cms2.hyp_type()[hypIdx] == 0) myType = 1;  // mm
	if (cms2.hyp_type()[hypIdx] == 1 || cms2.hyp_type()[hypIdx] == 2) myType=2; // em
	if (myType == 99) {
	  cout << "YUK:  unknown dilepton type = " << cms2.hyp_type()[hypIdx] << endl;
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
	  for (unsigned int ijet=0; 
	       ijet<(unsigned int)(cms2.hyp_njets()[hypIdx]); 
	       ijet++) {
	    // muon-jet cleaning
	    if (muJetClean) {
		LorentzVector vjet = cms2.hyp_jets_p4()[hypIdx][ijet];
		LorentzVector vlt  = cms2.hyp_lt_p4()[hypIdx];
		LorentzVector vll  = cms2.hyp_ll_p4()[hypIdx];
		if (abs(cms2.hyp_ll_id()[hypIdx]) == 13) {
		  if (dRBetweenVectors(vjet, vll) < 0.4) continue;
		}
		if (abs(cms2.hyp_lt_id()[hypIdx]) == 13) {
		  if (dRBetweenVectors(vjet, vlt) < 0.4) continue;
		}
	    }
	    float thisJetRescale =  globalJESscaleRescale;
	    blah = cms2.hyp_jets_pat_noCorrF()[hypIdx][ijet] * cms2.hyp_jets_p4()[hypIdx][ijet] * thisJetRescale;
	    //FIXME : this should be combined the same way as below
	    jp4.push_back(blah);
	    new_hyp_njets++;
	  }
	} else if (!useJPT) {
	  // Look among the hyp_jets
	  for (unsigned int ijet=0; ijet<(unsigned int)(cms2.hyp_njets()[hypIdx]); ijet++) {
	    // muon-jet cleaning
	    if (muJetClean) {
		LorentzVector vjet = cms2.hyp_jets_p4()[hypIdx][ijet];
		LorentzVector vlt  = cms2.hyp_lt_p4()[hypIdx];
		LorentzVector vll  = cms2.hyp_ll_p4()[hypIdx];
		if (abs(cms2.hyp_ll_id()[hypIdx]) == 13) {
		  if (dRBetweenVectors(vjet, vll) < 0.4) continue;
		}
		if (abs(cms2.hyp_lt_id()[hypIdx]) == 13) {
		  if (dRBetweenVectors(vjet, vlt) < 0.4) continue;
		}
	    }
	    float thisJetRescale =  globalJESscaleRescale;
	    blah = cms2.hyp_jets_p4()[hypIdx][ijet]* thisJetRescale;
	    if (blah.pt() > 30 && fabs(blah.eta()) < 2.4) {
	      jp4.push_back(blah);
	      new_hyp_njets++;
	    }
	  }
	  // Now look among the other jets
	  for (unsigned int ijet=0; ijet<cms2.hyp_other_jets_p4()[hypIdx].size(); ijet++) {
	    // muon-jet cleaning
	    if (muJetClean) {
		LorentzVector vjet = cms2.hyp_other_jets_p4()[hypIdx][ijet];
		LorentzVector vlt  = cms2.hyp_lt_p4()[hypIdx];
		LorentzVector vll  = cms2.hyp_ll_p4()[hypIdx];
		if (abs(cms2.hyp_ll_id()[hypIdx]) == 13) {
		  if (dRBetweenVectors(vjet, vll) < 0.4) continue;
		}
		if (abs(cms2.hyp_lt_id()[hypIdx]) == 13) {
		  if (dRBetweenVectors(vjet, vlt) < 0.4) continue;
		}
	    }
	    float thisJetRescale = globalJESscaleRescale;
	    blah = cms2.hyp_other_jets_p4()[hypIdx][ijet]* thisJetRescale;
	    if (blah.pt() > 30 && fabs(blah.eta()) < 2.4) {
	      jp4.push_back(blah);
	      new_hyp_njets++;
	    }
	  }
	} else {       
	  // This is with useJPT=true oldjets=false....
	  // We need to remove both electron and muon jets in this case
	  for (unsigned int ijet=0; ijet< cms2.evt_njpts(); ijet++) {
	    // jet cleaning (muons & electrons)
	    LorentzVector vjet = cms2.jpts_p4()[ijet];
	    LorentzVector vlt  = cms2.hyp_lt_p4()[hypIdx];
	    LorentzVector vll  = cms2.hyp_ll_p4()[hypIdx];
	    if (dRBetweenVectors(vjet, vll) < 0.4) continue;
	    if (dRBetweenVectors(vjet, vlt) < 0.4) continue;

	    float thisJetRescale = globalJESscaleRescale;
	    blah = cms2.jpts_p4()[ijet] * thisJetRescale;
	    if (blah.pt() > 30 && fabs(blah.eta()) < 2.4) {
	      jp4.push_back(blah);
	      new_hyp_njets++;
	    }
	  }
	}

	VofP4* new_hyp_jets_p4 = &jp4;
			   
	//     // Last chance to reject...
	int arrNjets = min(new_hyp_njets, 2);

	float pt_lt = cms2.hyp_lt_p4()[hypIdx].pt();
	float pt_ll = cms2.hyp_ll_p4()[hypIdx].pt();

	unsigned int i_lt = cms2.hyp_lt_index()[hypIdx];
	unsigned int i_ll = cms2.hyp_ll_index()[hypIdx];

	int id_lt = cms2.hyp_lt_id()[hypIdx];
	int id_ll = cms2.hyp_ll_id()[hypIdx];


	// jet count
	hnJet[myType]->Fill(min(new_hyp_njets,4), weight);
	hnJet[3]->Fill(min(new_hyp_njets,4), weight);
	if (inZmassWindow(cms2.hyp_p4().at(hypIdx).mass())) {
	  hnJetinZwindow[myType]->Fill(min(new_hyp_njets,4), weight);
	  hnJetinZwindow[3]->Fill(min(new_hyp_njets,4), weight);
	} else {
	  hnJetoutZwindow[myType]->Fill(min(new_hyp_njets,4), weight);
	  hnJetoutZwindow[3]->Fill(min(new_hyp_njets,4), weight);
	}

	// lepton Pt
	if (abs(id_lt) == 11) helePt[myType][arrNjets]->Fill(pt_lt, weight);
	if (abs(id_ll) == 11) helePt[myType][arrNjets]->Fill(pt_ll, weight);
	if (abs(id_lt) == 13) hmuPt[myType][arrNjets]->Fill(pt_lt, weight);
	if (abs(id_ll) == 13) hmuPt[myType][arrNjets]->Fill(pt_ll, weight);
	if (abs(id_lt) == 13) 
	  hmuPtFromSilicon[myType][arrNjets]->Fill(cms2.mus_trk_p4()[i_lt].pt(), weight);
	if (abs(id_ll) == 13)
	  hmuPtFromSilicon[myType][arrNjets]->Fill(cms2.mus_trk_p4()[i_ll].pt(), weight);
	hminLepPt[myType][arrNjets]->Fill(min(pt_ll, pt_lt), weight);
	hmaxLepPt[myType][arrNjets]->Fill(max(pt_ll, pt_lt), weight );
    
	if (abs(id_lt) == 11) helePt[3][arrNjets]->Fill(pt_lt, weight);
	if (abs(id_ll) == 11) helePt[3][arrNjets]->Fill(pt_ll, weight);
	if (abs(id_lt) == 13) hmuPt[3][arrNjets]->Fill(pt_lt, weight);
	if (abs(id_ll) == 13) hmuPt[3][arrNjets]->Fill(pt_ll, weight);
	if (abs(id_lt) == 13) 
	  hmuPtFromSilicon[3][arrNjets]->Fill(cms2.mus_trk_p4()[i_lt].pt(), weight);
	if (abs(id_ll) == 13) 
	  hmuPtFromSilicon[3][arrNjets]->Fill(cms2.mus_trk_p4()[i_ll].pt(), weight);
	hminLepPt[3][arrNjets]->Fill(min(pt_ll, pt_lt), weight);
	hmaxLepPt[3][arrNjets]->Fill(max(pt_ll, pt_lt), weight );


	// lepton Phi
	if (abs(id_lt) == 11) helePhi[myType][arrNjets]->Fill(cms2.hyp_lt_p4()[hypIdx].phi(), weight);
	if (abs(id_ll) == 11) helePhi[myType][arrNjets]->Fill(cms2.hyp_ll_p4()[hypIdx].phi(), weight);
	if (abs(id_lt) == 13) hmuPhi[myType][arrNjets]->Fill(cms2.hyp_lt_p4()[hypIdx].phi(), weight);
	if (abs(id_ll) == 13) hmuPhi[myType][arrNjets]->Fill(cms2.hyp_ll_p4()[hypIdx].phi(), weight);
    
	if (abs(id_lt) == 11) helePhi[3][arrNjets]->Fill(cms2.hyp_lt_p4()[hypIdx].phi(), weight);
	if (abs(id_ll) == 11) helePhi[3][arrNjets]->Fill(cms2.hyp_ll_p4()[hypIdx].phi(), weight);
	if (abs(id_lt) == 13) hmuPhi[3][arrNjets]->Fill(cms2.hyp_lt_p4()[hypIdx].phi(), weight);
	if (abs(id_ll) == 13) hmuPhi[3][arrNjets]->Fill(cms2.hyp_ll_p4()[hypIdx].phi(), weight);
    
	// dilepton mass
	hdilMass[myType][arrNjets]->Fill(cms2.hyp_p4()[hypIdx].mass(), weight);
	hdilMassTightWindow[myType][arrNjets]->Fill(cms2.hyp_p4()[hypIdx].mass(), weight);
	hdilMass[3][arrNjets]->Fill(cms2.hyp_p4()[hypIdx].mass(), weight);
	hdilMassTightWindow[3][arrNjets]->Fill(cms2.hyp_p4()[hypIdx].mass(), weight);
    
	// delta phi btw leptons
	double dphi = fabs(cms2.hyp_lt_p4()[hypIdx].phi() - cms2.hyp_ll_p4()[hypIdx].phi());
	if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
	hdphiLep[myType][arrNjets]->Fill(dphi, weight);
	hdphiLep[3][arrNjets]->Fill(dphi, weight);
			   
	// dphill vs mll, i.e. the 2d correlation between the previous two variables
	hdphillvsmll[myType][arrNjets]->Fill(cms2.hyp_p4()[hypIdx].mass(), dphi, weight);
	hdphillvsmll[3][arrNjets]->Fill(cms2.hyp_p4()[hypIdx].mass(), dphi, weight);
 
	// lepton Eta
	if (abs(id_lt) == 11) heleEta[myType][arrNjets]->Fill(cms2.hyp_lt_p4()[hypIdx].eta(), weight);
	if (abs(id_ll) == 11) heleEta[myType][arrNjets]->Fill(cms2.hyp_ll_p4()[hypIdx].eta(), weight);
	if (abs(id_lt) == 13) hmuEta[myType][arrNjets]->Fill(cms2.hyp_lt_p4()[hypIdx].eta(), weight);
	if (abs(id_ll) == 13) hmuEta[myType][arrNjets]->Fill(cms2.hyp_ll_p4()[hypIdx].eta(), weight);
    
	if (abs(id_lt) == 11) heleEta[3][arrNjets]->Fill(cms2.hyp_lt_p4()[hypIdx].eta(), weight);
	if (abs(id_ll) == 11) heleEta[3][arrNjets]->Fill(cms2.hyp_ll_p4()[hypIdx].eta(), weight);
	if (abs(id_lt) == 13) hmuEta[3][arrNjets]->Fill(cms2.hyp_lt_p4()[hypIdx].eta(), weight);
	if (abs(id_ll) == 13) hmuEta[3][arrNjets]->Fill(cms2.hyp_ll_p4()[hypIdx].eta(), weight);
       
      
	// electron trk isolation 
	double temp_lt_iso = cms2.hyp_lt_iso()[hypIdx];  // so that min works
	double temp_ll_iso = cms2.hyp_ll_iso()[hypIdx];  // so that min works
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
	double combIsoSum_lt = -1.;
	double combIsoSum_ll = -1.;
	double relIsoComb_lt = -1.;
	double relIsoComb_ll = -1.;
	double relIsoTrack_lt = -1.;
	double relIsoTrack_ll = -1.;
	double relIsoCalo_lt = -1.;
	double relIsoCalo_ll = -1.;
	if (abs(id_lt) == 13){
	  combIsoSum_lt = cms2.mus_pat_trackIso()[i_lt]
	    +cms2.mus_pat_ecalIso()[i_lt]
	    +cms2.mus_pat_hcalIso()[i_lt];
	  relIsoComb_lt = cms2.mus_p4()[i_lt].pt()/(cms2.mus_p4()[i_lt].pt() + combIsoSum_lt);
	  relIsoTrack_lt = cms2.mus_p4()[i_lt].pt()/(cms2.mus_p4()[i_lt].pt() + cms2.mus_pat_trackIso()[i_lt]);
	  relIsoCalo_lt = cms2.mus_p4()[i_lt].pt()/(cms2.mus_p4()[i_lt].pt() + cms2.mus_pat_ecalIso()[i_lt]+cms2.mus_pat_hcalIso()[i_lt]);
	}
	if (abs(id_ll) == 13){
	  combIsoSum_ll = cms2.mus_pat_trackIso()[i_ll]
	    +cms2.mus_pat_ecalIso()[i_ll]
	    +cms2.mus_pat_hcalIso()[i_ll];
	  relIsoComb_ll = cms2.mus_p4()[i_ll].pt()/(cms2.mus_p4()[i_ll].pt() + combIsoSum_ll);
	  relIsoTrack_ll = cms2.mus_p4()[i_ll].pt()/(cms2.mus_p4()[i_ll].pt() + cms2.mus_pat_trackIso()[i_ll]);
	  relIsoCalo_ll = cms2.mus_p4()[i_ll].pt()/(cms2.mus_p4()[i_ll].pt() + cms2.mus_pat_ecalIso()[i_ll]+cms2.mus_pat_hcalIso()[i_ll]);
	}
	if (abs(id_lt) == 13) hmuSumIso[myType][arrNjets]->Fill(min(combIsoSum_lt,24.99),weight);
	if (abs(id_lt) == 13) hmuSumIso[3][arrNjets]->Fill(min(combIsoSum_lt,24.99),weight);
	if (abs(id_ll) == 13) hmuSumIso[myType][arrNjets]->Fill(min(combIsoSum_ll,24.99),weight);
	if (abs(id_ll) == 13) hmuSumIso[3][arrNjets]->Fill(min(combIsoSum_ll,24.99),weight);
	//relative combined
	if (abs(id_lt) == 13) hmuRelIso[myType][arrNjets]->Fill(min(relIsoComb_lt,0.999),weight);
	if (abs(id_lt) == 13) hmuRelIso[3][arrNjets]->Fill(min(relIsoComb_lt,0.999),weight);
	if (abs(id_ll) == 13) hmuRelIso[myType][arrNjets]->Fill(min(relIsoComb_ll,0.999),weight);
	if (abs(id_ll) == 13) hmuRelIso[3][arrNjets]->Fill(min(relIsoComb_ll,0.999),weight);
	//relative tracker
	if (abs(id_lt) == 13) hmuRelIsoTrack[myType][arrNjets]->Fill(min(relIsoTrack_lt,0.999),weight);
	if (abs(id_lt) == 13) hmuRelIsoTrack[3][arrNjets]->Fill(min(relIsoTrack_lt,0.999),weight);
	if (abs(id_ll) == 13) hmuRelIsoTrack[myType][arrNjets]->Fill(min(relIsoTrack_ll,0.999),weight);
	if (abs(id_ll) == 13) hmuRelIsoTrack[3][arrNjets]->Fill(min(relIsoTrack_ll,0.999),weight);
	//relative calo
	if (abs(id_lt) == 13) hmuRelIsoCalo[myType][arrNjets]->Fill(min(relIsoCalo_lt,0.999),weight);
	if (abs(id_lt) == 13) hmuRelIsoCalo[3][arrNjets]->Fill(min(relIsoCalo_lt,0.999),weight);
	if (abs(id_ll) == 13) hmuRelIsoCalo[myType][arrNjets]->Fill(min(relIsoCalo_ll,0.999),weight);
	if (abs(id_ll) == 13) hmuRelIsoCalo[3][arrNjets]->Fill(min(relIsoCalo_ll,0.999),weight);

	//electrons now
	if (abs(id_lt) == 11){
	  combIsoSum_lt = cms2.els_pat_trackIso()[i_lt]
	    +cms2.els_pat_ecalIso()[i_lt]
	    +cms2.els_pat_hcalIso()[i_lt];
	  relIsoComb_lt = cms2.els_p4()[i_lt].pt()/(cms2.els_p4()[i_lt].pt() + combIsoSum_lt);
	  relIsoTrack_lt = cms2.els_p4()[i_lt].pt()/(cms2.els_p4()[i_lt].pt() + cms2.els_pat_trackIso()[i_lt]);
	  relIsoCalo_lt = cms2.els_p4()[i_lt].pt()/(cms2.els_p4()[i_lt].pt() + cms2.els_pat_ecalIso()[i_lt]+cms2.els_pat_hcalIso()[i_lt]);
	}
	if (abs(id_ll) == 11){
	  combIsoSum_ll = cms2.els_pat_trackIso()[i_ll]
	    +cms2.els_pat_ecalIso()[i_ll]
	    +cms2.els_pat_hcalIso()[i_ll];
	  relIsoComb_ll = cms2.els_p4()[i_ll].pt()/(cms2.els_p4()[i_ll].pt() + combIsoSum_ll);
	  relIsoTrack_ll = cms2.els_p4()[i_ll].pt()/(cms2.els_p4()[i_ll].pt() + cms2.els_pat_trackIso()[i_ll]);
	  relIsoCalo_ll = cms2.els_p4()[i_ll].pt()/(cms2.els_p4()[i_ll].pt() + cms2.els_pat_ecalIso()[i_ll]+cms2.els_pat_hcalIso()[i_ll]);
	}
	if (abs(id_lt) == 11) helSumIso[myType][arrNjets]->Fill(min(combIsoSum_lt,24.99),weight);
	if (abs(id_lt) == 11) helSumIso[3][arrNjets]->Fill(min(combIsoSum_lt,24.99),weight);
	if (abs(id_ll) == 11) helSumIso[myType][arrNjets]->Fill(min(combIsoSum_ll,24.99),weight);
	if (abs(id_ll) == 11) helSumIso[3][arrNjets]->Fill(min(combIsoSum_ll,24.99),weight);
	//relative combined
	if (abs(id_lt) == 11) helRelIso[myType][arrNjets]->Fill(min(relIsoComb_lt,0.999),weight);
	if (abs(id_lt) == 11) helRelIso[3][arrNjets]->Fill(min(relIsoComb_lt,0.999),weight);
	if (abs(id_ll) == 11) helRelIso[myType][arrNjets]->Fill(min(relIsoComb_ll,0.999),weight);
	if (abs(id_ll) == 11) helRelIso[3][arrNjets]->Fill(min(relIsoComb_ll,0.999),weight);
	//relative tracker
	if (abs(id_lt) == 11) helRelIsoTrack[myType][arrNjets]->Fill(min(relIsoTrack_lt,0.999),weight);
	if (abs(id_lt) == 11) helRelIsoTrack[3][arrNjets]->Fill(min(relIsoTrack_lt,0.999),weight);
	if (abs(id_ll) == 11) helRelIsoTrack[myType][arrNjets]->Fill(min(relIsoTrack_ll,0.999),weight);
	if (abs(id_ll) == 11) helRelIsoTrack[3][arrNjets]->Fill(min(relIsoTrack_ll,0.999),weight);
	//relative calo
	if (abs(id_lt) == 11) helRelIsoCalo[myType][arrNjets]->Fill(min(relIsoCalo_lt,0.999),weight);
	if (abs(id_lt) == 11) helRelIsoCalo[3][arrNjets]->Fill(min(relIsoCalo_lt,0.999),weight);
	if (abs(id_ll) == 11) helRelIsoCalo[myType][arrNjets]->Fill(min(relIsoCalo_ll,0.999),weight);
	if (abs(id_ll) == 11) helRelIsoCalo[3][arrNjets]->Fill(min(relIsoCalo_ll,0.999),weight);


      
	// dilepton pt
	hdilPt[myType][arrNjets]->Fill(cms2.hyp_p4()[hypIdx].pt(), weight);
	hdilPt[3][arrNjets]->Fill(cms2.hyp_p4()[hypIdx].pt(), weight);
    
	// Met and Met phi
	hmet[myType][arrNjets]->Fill(cms2.hyp_met()[hypIdx], weight);      
	hmetPhi[myType][arrNjets]->Fill(cms2.hyp_metPhi()[hypIdx], weight);      
	hmet[3][arrNjets]->Fill(cms2.hyp_met()[hypIdx], weight);      
	hmetPhi[3][arrNjets]->Fill(cms2.hyp_metPhi()[hypIdx], weight);      
	// pat Met and Met phi
	hpatmet[myType][arrNjets]->Fill(cms2.met_pat_metCor(), weight);      
	hpatmetPhi[myType][arrNjets]->Fill(cms2.met_pat_metPhiCor(), weight);      
	hpatmet[3][arrNjets]->Fill(cms2.met_pat_metCor(), weight);      
	hpatmetPhi[3][arrNjets]->Fill(cms2.met_pat_metPhiCor(), weight);      
	// tc Met and Met phi
	htcmet[myType][arrNjets]->Fill(cms2.evt_tcmet(), weight);      
	htcmetPhi[myType][arrNjets]->Fill(cms2.evt_tcmetPhi(), weight);      
	htcmet[3][arrNjets]->Fill(cms2.evt_tcmet(), weight);      
	htcmetPhi[3][arrNjets]->Fill(cms2.evt_tcmetPhi(), weight);      
    
	// Met vs dilepton Pt
	hmetVsDilepPt[myType][arrNjets]->Fill(cms2.hyp_met()[hypIdx], cms2.hyp_p4()[hypIdx].pt(), weight);
	hmetVsDilepPt[3][arrNjets]->Fill(cms2.hyp_met()[hypIdx], cms2.hyp_p4()[hypIdx].pt(), weight);
	//pat  Met vs dilepton Pt
	hpatmetVsDilepPt[myType][arrNjets]->Fill(cms2.met_pat_metCor(), cms2.hyp_p4()[hypIdx].pt(), weight);
	hpatmetVsDilepPt[3][arrNjets]->Fill(cms2.met_pat_metCor(), cms2.hyp_p4()[hypIdx].pt(), weight);
	//tc  Met vs dilepton Pt
	htcmetVsDilepPt[myType][arrNjets]->Fill(cms2.evt_tcmet(), cms2.hyp_p4()[hypIdx].pt(), weight);
	htcmetVsDilepPt[3][arrNjets]->Fill(cms2.evt_tcmet(), cms2.hyp_p4()[hypIdx].pt(), weight);
    
	// Met over dilepton Pt vs deltaphi btw the two
	double dphi2 = fabs(cms2.hyp_p4()[hypIdx].phi() - cms2.hyp_metPhi()[hypIdx]);
	if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
	dphi2 = TMath::Pi() - dphi2;  // changed the definition CC 28 March 08
	hmetOverPtVsDphi[myType][arrNjets]->Fill(cms2.hyp_met()[hypIdx]/cms2.hyp_p4()[hypIdx].pt(), dphi2, weight);
	hmetOverPtVsDphi[3][arrNjets]->Fill(cms2.hyp_met()[hypIdx]/cms2.hyp_p4()[hypIdx].pt(), dphi2, weight);
	//pat Met over dilepton Pt vs deltaphi btw the two
	dphi2 = fabs(cms2.hyp_p4()[hypIdx].phi() - cms2.met_pat_metPhiCor());
	if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
	dphi2 = TMath::Pi() - dphi2;  // changed the definition CC 28 March 08
	hpatmetOverPtVsDphi[myType][arrNjets]->Fill(cms2.met_pat_metCor()/cms2.hyp_p4()[hypIdx].pt(), dphi2, weight);
	hpatmetOverPtVsDphi[3][arrNjets]->Fill(cms2.met_pat_metCor()/cms2.hyp_p4()[hypIdx].pt(), dphi2, weight);
	//tc Met over dilepton Pt vs deltaphi btw the two
	dphi2 = fabs(cms2.hyp_p4()[hypIdx].phi() - cms2.evt_tcmetPhi());
	if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
	dphi2 = TMath::Pi() - dphi2;  // changed the definition CC 28 March 08
	htcmetOverPtVsDphi[myType][arrNjets]->Fill(cms2.evt_tcmet()/cms2.hyp_p4()[hypIdx].pt(), dphi2, weight);
	htcmetOverPtVsDphi[3][arrNjets]->Fill(cms2.evt_tcmet()/cms2.hyp_p4()[hypIdx].pt(), dphi2, weight);
    
      

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
	hnJetinZwindow[i] = new TH1F(Form("%s_hnJetinZwindow_%s",prefix,suffixall[i]),Form("%s_hnJetinZwindow_%s",prefix,suffixall[i]),
				     5,0.,5.);	
	hnJetoutZwindow[i] = new TH1F(Form("%s_hnJetoutZwindow_%s",prefix,suffixall[i]),Form("%s_hnJetoutZwindow_%s",prefix,suffixall[i]),
				     5,0.,5.);	
	
	hnJet[i]->SetDirectory(rootdir);
	hnJet[i]->GetXaxis()->SetTitle("nJets");

	hnJetinZwindow[i]->SetDirectory(rootdir);
	hnJetinZwindow[i]->GetXaxis()->SetTitle("nJets");

	hnJetoutZwindow[i]->SetDirectory(rootdir);
	hnJetoutZwindow[i]->GetXaxis()->SetTitle("nJets");

	for(int k = 0; k<5; k++) {
	  hnJet[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
	  hnJet[i]->GetXaxis()->SetLabelSize(0.07);
	  
	  hnJetinZwindow[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
	  hnJetinZwindow[i]->GetXaxis()->SetLabelSize(0.07);
	  
	  hnJetoutZwindow[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
	  hnJetoutZwindow[i]->GetXaxis()->SetLabelSize(0.07);
	  
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
				300, 0., 300.);
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
      //pat
      hpatmet[i][j] = new TH1F(Form("%s_hpatmet_%s",prefix,suffix[i]),Form("%s_patmet_%s",prefix,suffix[i]),20,0.,200.);
      hpatmet[i][j]->SetDirectory(rootdir);
      hpatmet[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      hpatmetPhi[i][j] = new TH1F(Form("%s_hpatmetPhi_%s",prefix,suffix[i]),Form("%s_patmetPhi_%s",prefix,suffix[i]),
			       50,-1*TMath::Pi(), TMath::Pi());
      hpatmetPhi[i][j]->SetDirectory(rootdir);
      hpatmetPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hpatmetVsDilepPt[i][j] = new TH2F(Form("%s_hpatmetVsDilepPt_%s",prefix,suffix[i]),
				     Form("%s_patmetVsDilepPt_%s",prefix,suffix[i]),
				     100,0.,200.,100,0.,200.);
      hpatmetVsDilepPt[i][j]->SetDirectory(rootdir);
      hpatmetVsDilepPt[i][j]->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
      hpatmetVsDilepPt[i][j]->GetYaxis()->SetTitle("Met (GeV)");
    
      hpatmetOverPtVsDphi[i][j] = new TH2F(Form("%s_hpatmetOverPtVsDphi_%s",prefix,suffix[i]),
					Form("%s_patmetOverPtVsDphi_%s",prefix,suffix[i]),
					30,0.,3.,25,0.,TMath::Pi());
      hpatmetOverPtVsDphi[i][j]->SetDirectory(rootdir);
      hpatmetVsDilepPt[i][j]->GetXaxis()->SetTitle("#Delta#Phi");
      hpatmetVsDilepPt[i][j]->GetYaxis()->SetTitle("MET/Pt_{ll}");
    
      //tc
      htcmet[i][j] = new TH1F(Form("%s_htcmet_%s",prefix,suffix[i]),Form("%s_tcmet_%s",prefix,suffix[i]),20,0.,200.);
      htcmet[i][j]->SetDirectory(rootdir);
      htcmet[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      htcmetPhi[i][j] = new TH1F(Form("%s_htcmetPhi_%s",prefix,suffix[i]),Form("%s_tcmetPhi_%s",prefix,suffix[i]),
			       50,-1*TMath::Pi(), TMath::Pi());
      htcmetPhi[i][j]->SetDirectory(rootdir);
      htcmetPhi[i][j]->GetXaxis()->SetTitle("#phi");

      htcmetVsDilepPt[i][j] = new TH2F(Form("%s_htcmetVsDilepPt_%s",prefix,suffix[i]),
				     Form("%s_tcmetVsDilepPt_%s",prefix,suffix[i]),
				     100,0.,200.,100,0.,200.);
      htcmetVsDilepPt[i][j]->SetDirectory(rootdir);
      htcmetVsDilepPt[i][j]->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
      htcmetVsDilepPt[i][j]->GetYaxis()->SetTitle("Met (GeV)");
    
      htcmetOverPtVsDphi[i][j] = new TH2F(Form("%s_htcmetOverPtVsDphi_%s",prefix,suffix[i]),
					Form("%s_tcmetOverPtVsDphi_%s",prefix,suffix[i]),
					30,0.,3.,25,0.,TMath::Pi());
      htcmetOverPtVsDphi[i][j]->SetDirectory(rootdir);
      htcmetVsDilepPt[i][j]->GetXaxis()->SetTitle("#Delta#Phi");
      htcmetVsDilepPt[i][j]->GetYaxis()->SetTitle("MET/Pt_{ll}");
    
    

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
    
      hmuSumIso[i][j] = new TH1F(Form("%s_hmuSumIso_%s",prefix,suffix[i]),Form("%s_hmuSumIso_%s",prefix,suffix[i]),
				 100, 0., 25.);
      hmuSumIso[i][j]->SetDirectory(rootdir);
      hmuSumIso[i][j]->GetXaxis()->SetTitle("#SigmaPt");
      helSumIso[i][j] = new TH1F(Form("%s_helSumIso_%s",prefix,suffix[i]),Form("%s_helSumIso_%s",prefix,suffix[i]),
				 100, 0., 25.);
      helSumIso[i][j]->SetDirectory(rootdir);
      helSumIso[i][j]->GetXaxis()->SetTitle("#SigmaPt");
    
      hmuRelIso[i][j] = new TH1F(Form("%s_hmuRelIso_%s",prefix,suffix[i]),Form("%s_hmuRelIso_%s",prefix,suffix[i]),
				  100, 0., 1.0001);
      hmuRelIso[i][j]->SetDirectory(rootdir);
      helRelIso[i][j] = new TH1F(Form("%s_helRelIso_%s",prefix,suffix[i]),Form("%s_helRelIso_%s",prefix,suffix[i]),
				  100, 0., 1.0001);
      helRelIso[i][j]->SetDirectory(rootdir);
      // tracker
      hmuRelIsoTrack[i][j] = new TH1F(Form("%s_hmuRelIsoTrack_%s",prefix,suffix[i]),Form("%s_hmuRelIsoTrack_%s",prefix,suffix[i]),
				  100, 0., 1.0001);
      hmuRelIsoTrack[i][j]->SetDirectory(rootdir);
      helRelIsoTrack[i][j] = new TH1F(Form("%s_helRelIsoTrack_%s",prefix,suffix[i]),Form("%s_helRelIsoTrack_%s",prefix,suffix[i]),
				  100, 0., 1.0001);
      helRelIsoTrack[i][j]->SetDirectory(rootdir);
      // calorimeter
      hmuRelIsoCalo[i][j] = new TH1F(Form("%s_hmuRelIsoCalo_%s",prefix,suffix[i]),Form("%s_hmuRelIsoCalo_%s",prefix,suffix[i]),
				  100, 0., 1.0001);
      hmuRelIsoCalo[i][j]->SetDirectory(rootdir);
      helRelIsoCalo[i][j] = new TH1F(Form("%s_helRelIsoCalo_%s",prefix,suffix[i]),Form("%s_helRelIsoCalo_%s",prefix,suffix[i]),
				  100, 0., 1.0001);
      helRelIsoCalo[i][j]->SetDirectory(rootdir);


      if (j==0){
	hnJet[i]->Sumw2();
	hnJetinZwindow[i]->Sumw2();
	hnJetoutZwindow[i]->Sumw2();
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
      hpatmet[i][j]->Sumw2();
      hpatmetPhi[i][j]->Sumw2();
      htcmet[i][j]->Sumw2();
      htcmetPhi[i][j]->Sumw2();
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
      helSumIso[i][j]->Sumw2();
      hmuRelIso[i][j]->Sumw2();
      helRelIso[i][j]->Sumw2();
      hmuRelIsoTrack[i][j]->Sumw2();
      helRelIsoTrack[i][j]->Sumw2();
      hmuRelIsoCalo[i][j]->Sumw2();
      helRelIsoCalo[i][j]->Sumw2();

    }
  }//channel loop
}//CMS2::bookHistos()
