//now make the source file
#include <iostream>
#include <vector>
#include <set>

#include "TChain.h"
#include "TFile.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace std;

#ifndef __CINT__
#include "CMS2.h"
CMS2 cms2;
#endif
#include "../Tools/selections.C"
#include "../Tools/utilities.C"

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC                                                                                   
     unsigned long int run, event;
     float trks_d0;
     float hyp_lt_pt, hyp_lt_eta, hyp_lt_phi;
     bool operator < (const DorkyEventIdentifier &) const;
     bool operator == (const DorkyEventIdentifier &) const;
};

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
     if (run != other.run)
          return run < other.run;
     if (event != other.event)
          return event < other.event;
     // the floating point numbers are not easy, because we're                                                                                        
     // comapring ones that are truncated (because they were written                                                                                  
     // to file and read back in) with ones that are not truncated.                                                                                   
     if (fabs(trks_d0 - other.trks_d0) > 1e-6 * trks_d0)
       return trks_d0 < other.trks_d0;
     if (fabs(hyp_lt_pt - other.hyp_lt_pt) > 1e-6 * hyp_lt_pt)
       return hyp_lt_pt < other.hyp_lt_pt;
     if (fabs(hyp_lt_eta - other.hyp_lt_eta) > 1e-6 * hyp_lt_eta)
       return hyp_lt_eta < other.hyp_lt_eta;
     if (fabs(hyp_lt_phi - other.hyp_lt_phi) > 1e-6 * hyp_lt_phi)
       return hyp_lt_phi < other.hyp_lt_phi;
     // if the records are exactly the same, then r1 is not less than                                                                                 
     // r2.  Duh!                                                                                                                                     
     return false;
}

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
     if (run != other.run)
          return false;
     if (event != other.event)
          return false;
     // the floating point numbers are not easy, because we're                                                                                        
     // comapring ones that are truncated (because they were written                                                                                  
     // to file and read back in) with ones that are not truncated.                                                                                   
     if (fabs(trks_d0 - other.trks_d0) > 1e-6 * trks_d0)
          return false;
     if (fabs(hyp_lt_pt - other.hyp_lt_pt) > 1e-6 * hyp_lt_pt)
          return false;
     if (fabs(hyp_lt_eta - other.hyp_lt_eta) > 1e-6 * hyp_lt_eta)
          return false;
     if (fabs(hyp_lt_phi - other.hyp_lt_phi) > 1e-6 * hyp_lt_phi)
          return false;
     return true;
}

static std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id)
{
     std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
          already_seen.insert(id);
     return !ret.second;
}

int ScanChain( TChain* chain, char * prefix="", int specDY=-1, float kFactor=1.0) {

  // Make sure the specDY flags is kosher                                                                                                             
  if (specDY < -1 || specDY > 2) {
    std::cout << "specDY flag is not allowed...quit" << std::endl;
    return 1;
  }

  // clear list of duplicates                                                                                                                         
  already_seen.clear();
  int duplicates_total_n = 0;
  double duplicates_total_weight = 0;

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;

  const unsigned int allBuckets = 20;
  char *suffix[allBuckets+1];
  suffix[0]  = "mpmpmp";
  suffix[1]  = "mpmpmm";
  suffix[2]  = "mpmmmm";
  suffix[3]  = "mmmmmm";
  suffix[4]  = "mpmpep";
  suffix[5]  = "mpmpem";
  suffix[6]  = "mpmmep";
  suffix[7]  = "mpmmem";
  suffix[8]  = "mmmmep";
  suffix[9]  = "mmmmem";
  suffix[10] = "mpepep";
  suffix[11] = "mpepem";
  suffix[12] = "mpemem";
  suffix[13] = "mmepep";
  suffix[14] = "mmepem";
  suffix[15] = "mmemem";
  suffix[16] = "epepep";
  suffix[17] = "epepem";
  suffix[18] = "epemem";
  suffix[19] = "ememem";
  suffix[20] = "all";

  // declare histograms
  TH1F* hLowestPtGood[allBuckets];

  TH1F* hMETFinal[allBuckets];
  TH1F* hMETAllFinal[allBuckets];

  TH1F* hZmassFinal[allBuckets];
  TH1F* h2ndZmass[allBuckets];

  TH1F* hNjetsFinal[allBuckets];
  TH1F* hNjetsNoMetCut[allBuckets];
  TH1F* hNjetsGoodLeptonsVeto[allBuckets];
  TH1F* hNjets2ndZVeto[allBuckets];
  TH1F* hNjetsBothLeptonsVeto[allBuckets];
  TH1F* hNjetsBucket[allBuckets];

  TH1F* hNGoodLeptons[allBuckets];
  TH1F* hNAllLeptons[allBuckets];
    
  for (unsigned int i=0; i<=allBuckets; ++i) {

    // plots of pt of lowest pt lepton in trilepton candidate, only requires goodIsolated<lepton> 20/X/10
    const int nbins = 4;
    float bins[nbins+1] = {0.,5.,10.,20.,30.};
    hLowestPtGood[i] = book1DVarHist(Form("%s_hLowestPtGood_%s",prefix,suffix[i]),Form("%s_nLowestPtGood_%s",prefix,suffix[i]),nbins,bins,"p_{T} [GeV]","trilepton cand.");   

    // MET plots, require goodIsolated<lepton> 20/X/10 and Z candidate in Z window cut
    const int nBinsMET = 50;
    const float lowBinMET = 0.;
    const float highBinMET = 500;
    hMETFinal[i] = book1DHist(Form("%s_hMETFinal_%s",prefix,suffix[i]),Form("%s_nMETFinal_%s",prefix,suffix[i]),nBinsMET,lowBinMET,highBinMET,"MET [GeV]","trilepton cand.");   
    hMETAllFinal[i] = book1DHist(Form("%s_hMETAllFinal_%s",prefix,suffix[i]),Form("%s_nMETAllFinal_%s",prefix,suffix[i]),nBinsMET,lowBinMET,highBinMET,"MET [GeV]","trilepton cand.");   

    // Z mass plot, requires MET cut and goodIsolated<lepton> 20/X/10
    const int nBinsZmass = 100;
    const float lowBinZmass = 0.;
    const float highBinZmass = 200.;
    hZmassFinal[i] = book1DHist(Form("%s_hZmassFinal_%s",prefix,suffix[i]),Form("%s_nZmassFinal_%s",prefix,suffix[i]),nBinsZmass,lowBinZmass,highBinZmass,"Mass [GeV]","trilepton cand.");   

    // Z mass plot of 2nd Z candidate, requires Z candidate in Z windows, goodIsolated<lepton> 20/X/10 and MET cut 
    h2ndZmass[i] = book1DHist(Form("%s_h2ndZmass_%s",prefix,suffix[i]),Form("%s_n2ndZmass_%s",prefix,suffix[i]),nBinsZmass,lowBinZmass,highBinZmass,"Mass [GeV]","trilepton cand.");   

    // njets plot, requires Z candidate in Z windows, goodIsolated<lepton> 20/X/10 and MET cut
    const int nBinsNjets = 10;
    const float lowBinNjets = 0.;
    const float highBinNjets = 10.;
    hNjetsFinal[i] = book1DHist(Form("%s_hNjetsFinal_%s",prefix,suffix[i]),Form("%s_nNjetsFinal_%s",prefix,suffix[i]),nBinsNjets,lowBinNjets,highBinNjets,"Njets","trilepton cand.");   
    hNjetsNoMetCut[i] = book1DHist(Form("%s_hNjetsNoMetCut_%s",prefix,suffix[i]),Form("%s_nNjetsNoMetCut_%s",prefix,suffix[i]),nBinsNjets,lowBinNjets,highBinNjets,"Njets","trilepton cand.");   

    // njets plot, requires Z candidate in Z windows, goodIsolated<lepton> 20/X/10 and MET cut and goodLeptonsVeto
    hNjetsGoodLeptonsVeto[i] = book1DHist(Form("%s_hNjetsGoodLeptonsVeto_%s",prefix,suffix[i]),Form("%s_nNjetsGoodLeptonsVeto_%s",prefix,suffix[i]),nBinsNjets,lowBinNjets,highBinNjets,"Njets","trilepton cand.");   

    // njets plot, requires Z candidate in Z windows, goodIsolated<lepton> 20/X/10 and MET cut and 2nd Z candidate veto
    hNjets2ndZVeto[i] = book1DHist(Form("%s_hNjets2ndZVeto_%s",prefix,suffix[i]),Form("%s_nNjets2ndZVeto_%s",prefix,suffix[i]),nBinsNjets,lowBinNjets,highBinNjets,"Njets","trilepton cand.");   

    // njets plot, requires Z candidate in Z windows, goodIsolated<lepton> 20/X/10 and MET cut, goodLeptonsVeto and 2nd Z candidate veto
    hNjetsBothLeptonsVeto[i] = book1DHist(Form("%s_hNjetsBothLeptonsVeto_%s",prefix,suffix[i]),Form("%s_nNjetsBothLeptonsVeto_%s",prefix,suffix[i]),nBinsNjets,lowBinNjets,highBinNjets,"Njets","trilepton cand.");   

    // njets plot with the same cuts for all buckets: goodIsolated<lepton> 20/X/10 and MET cut
    hNjetsBucket[i] = book1DHist(Form("%s_hNjetsBucket_%s",prefix,suffix[i]),Form("%s_nNjetsBucket_%s",prefix,suffix[i]),nBinsNjets,lowBinNjets,highBinNjets,"Njets","trilepton cand.");   

    // nlepton plots, requires Z candidate in Z windows, goodIsolated<lepton> 20/X/10 and MET cut
    const int nBinsNGoodLeptons = 20;
    const float lowBinNGoodLeptons = 0.;
    const float highBinNGoodLeptons = 20.;
    hNGoodLeptons[i] = book1DHist(Form("%s_hNGoodLeptons_%s",prefix,suffix[i]),Form("%s_nNGoodLeptons_%s",prefix,suffix[i]),nBinsNGoodLeptons,lowBinNGoodLeptons,highBinNGoodLeptons,"NLeptons","trilepton cand.");   
    hNAllLeptons[i] = book1DHist(Form("%s_hNAllLeptons_%s",prefix,suffix[i]),Form("%s_nNAllLeptons_%s",prefix,suffix[i]),nBinsNGoodLeptons,lowBinNGoodLeptons,highBinNGoodLeptons,"NLeptons","trilepton cand.");   

  }

  // trilepton candidate per bucket
  float trilepCounter[allBuckets];
  for ( unsigned int i = 0; i < allBuckets; ++i ) {
    trilepCounter[i] = 0.;
  }

  //CONSTANTS
  const float zmass = 91.19; //just making sure that I use the same Z mass everywhere!

  // CUTS
  const int   goodLeptonsCut        = 3;   // good lepton cut
  const float triggerLeptonMinPtCut = 20.; // one of the leptons of the trilepton canidate has to have at least this pt
  const float leptonMinPtCut        = 10.; // all leptons of the trilepton canidate have to have at least this pt
  const float electronMETAllCut     = 40.; // cut on METAll if lepton not belonging to primary Z is an electron
  const float muonMETAllCut         = 20.; // cut on METAll if lepton not belonging to primary Z is an muon

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

      if (cms2.trks_d0().size() == 0)
	continue;
      DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.trks_d0()[0],
				  cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
      if (is_duplicate(id)) {
	duplicates_total_n++;
	duplicates_total_weight += cms2.evt_scale1fb();
	continue;
      }

      // Progress feedback to the user
      if ((nEventsTotal)%1000 == 0) std::cout << "Processing event: " << nEventsTotal << std::endl;

      // The event weight including the kFactor (scaled to 1 fb-1)
      float weight = cms2.evt_scale1fb() * kFactor;

      // special handling for DY
      bool processEvent=true;
      if (specDY == 0) {
	if ( !isDYee() ) processEvent = false;
      } else if (specDY == 1) {
	if ( !isDYmm() ) processEvent = false;
      } else if (specDY == 2) {
	if ( !isDYtt() ) processEvent = false;
      }
      if (!processEvent) continue;
      
      // metAll correct, buggy in ntuples
      double metAll = cms2.evt_met();
      double metAllPhi = cms2.evt_metPhi();
      
      for ( unsigned int muon = 0;
	    muon < cms2.mus_p4().size();
            ++muon ) {
	correctMETmuons_crossedE(metAll, metAllPhi, 
				 cms2.mus_p4()[muon].pt(), cms2.mus_p4()[muon].phi(),
				 cms2.mus_trk_p4()[muon].theta(), cms2.mus_trk_p4()[muon].phi(),
				 cms2.mus_e_em()[muon], cms2.mus_e_had()[muon],cms2.mus_e_ho()[muon]);
      }

      // loop over trilepton candidates
      for ( unsigned int cand = 0; 
	    cand < cms2.hyp_trilep_bucket().size();
	    ++cand ) {
	
	unsigned int bucket = cms2.hyp_trilep_bucket()[cand],;
	int first = cms2.hyp_trilep_first_index()[cand];
	int second = cms2.hyp_trilep_second_index()[cand];
	int third = cms2.hyp_trilep_third_index()[cand];

	// count good leptons and all leptons
	// good lepton is goodIsolated<lepton>
	// lepton is good<lepton>
	int goodLeptons = 0;
	int allLeptons  = 0;
	for ( unsigned int i = 0; i < cms2.evt_nels(); ++i ) {
	  if ( goodElectronIsolated(i) ) ++goodLeptons;
	  if ( goodElectronWithoutIsolation(i) ) ++allLeptons;
	}
	for ( unsigned int i = 0; i < cms2.mus_p4().size(); ++i ) {
	  if ( goodMuonIsolated(i) ) ++goodLeptons;
	  if ( goodMuonWithoutIsolation(i) ) ++allLeptons;
	}

	// CUT: goodLeptonIsolated
	if ( !goodLeptonIsolated(bucket,
				 first,
				 second,
				 third) ) continue;
	
	// determine pt of lowest pt lepton
	float lowest = ptLowestPtLepton(bucket,first,second,third);
      
	hLowestPtGood[bucket]->Fill(lowest,weight);
	hLowestPtGood[allBuckets]->Fill(lowest,weight);
	
	// CUT: do nothing unless pt = 20/X/10
	if ( !passTriggerLeptonMinPtCut(bucket,first,second,third,triggerLeptonMinPtCut) ) continue;
	if ( lowest < leptonMinPtCut ) continue;

	// identify primary Z candidate
	unsigned int zArray[3] = {999999,999999,999999};
	calcPrimZ(bucket,first,second,third,zmass,zArray);
	// calculate mass of primary Z candiate
	float mZPrim = calcPrimZMass(bucket,zArray[0],zArray[1]);

 	std::vector<LorentzVector> correctedJets = correctJetsForElectrons(bucket,first,second,third);
	// distributions for Zmass
	if ( passMETAllCut(bucket,metAll,electronMETAllCut,muonMETAllCut) ){
	  hNjetsBucket[bucket]->Fill(correctedJets.size(),weight);
	  hNjetsBucket[allBuckets]->Fill(correctedJets.size(),weight);

	  hZmassFinal[bucket]->Fill(mZPrim, weight);
	  hZmassFinal[allBuckets]->Fill(mZPrim, weight);
	}
	// CUT: do nothing unless one of the opposite sign same flavor combos is in the zmass window:
	if( !inZmassWindow(mZPrim) ) continue;

	hMETFinal[bucket]->Fill(cms2.hyp_trilep_met()[cand],weight);
	hMETFinal[allBuckets]->Fill(cms2.hyp_trilep_met()[cand],weight);
	hMETAllFinal[bucket]->Fill(metAll,weight);
	hMETAllFinal[allBuckets]->Fill(metAll,weight);

	hNjetsNoMetCut[bucket]->Fill(correctedJets.size(),weight);
	hNjetsNoMetCut[allBuckets]->Fill(correctedJets.size(),weight);

	// CUT: MET cut
	if ( !passMETAllCut(bucket,metAll,electronMETAllCut,muonMETAllCut) ) continue;

	hNjetsFinal[bucket]->Fill(correctedJets.size(),weight);
	hNjetsFinal[allBuckets]->Fill(correctedJets.size(),weight);

	hNGoodLeptons[bucket]->Fill(goodLeptons,weight);
	hNGoodLeptons[allBuckets]->Fill(goodLeptons,weight);

	hNAllLeptons[bucket]->Fill(allLeptons,weight);
	hNAllLeptons[allBuckets]->Fill(allLeptons,weight);

	float mZSec = calcSecZMass(bucket,zArray,zmass);

	if ( mZSec == 999999. ) {
	  // if no 2nd Z candidate is found, fill 0
	  h2ndZmass[bucket]->Fill(0., weight);
	  h2ndZmass[allBuckets]->Fill(0., weight);
	} else {
	  h2ndZmass[bucket]->Fill(mZSec, weight);
	  h2ndZmass[allBuckets]->Fill(mZSec, weight);
	}
	  
	// njet plots with all cuts except goodLeptonsVeto
	if ( !inZmassWindow(mZSec) ) {
	  hNjets2ndZVeto[bucket]->Fill(correctedJets.size(),weight);
	  hNjets2ndZVeto[allBuckets]->Fill(correctedJets.size(),weight);
	}

	// CUT: goodLeptonsVeto
	if ( goodLeptons > goodLeptonsCut ) continue;

	hNjetsGoodLeptonsVeto[bucket]->Fill(correctedJets.size(),weight);
	hNjetsGoodLeptonsVeto[allBuckets]->Fill(correctedJets.size(),weight);

	// CUT: veto 2nd Z candidate
	if ( inZmassWindow(mZSec) ) continue;
	trilepCounter[bucket] += weight;

	hNjetsBothLeptonsVeto[bucket]->Fill(correctedJets.size(),weight);
	hNjetsBothLeptonsVeto[allBuckets]->Fill(correctedJets.size(),weight);
	
      }
    }
  }
  
  std::cout << std::endl;
  for ( unsigned int i = 0; i < allBuckets; ++i ) {
    cout << "Bucket: " << i << " entries: " << trilepCounter[i] << endl;
  }
  std::cout << std::endl;

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  
  return 0;
}
