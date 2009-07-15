//now make the source file
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TStopwatch.h"
#include <algorithm>
#include <set>
#include "TCanvas.h"
#include "TRegexp.h"
#include "TLorentzVector.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"


using namespace std;

#ifndef __CINT__
#include "../CORE/CMS2.cc"
#include "../CORE/utilities.cc"
#include "../CORE/selections.cc"
#endif

TH1F* hypos_total;
TH1F* hypos_total_weighted;

enum Sample {WW, WZ, ZZ, Wjets, DY, DYee, DYmm, DYtt, ttbar, tW, LM0x, LM1x, LM2x, LM3x, LM4x, LM5x, LM6x, LM7x, LM8x, LM9x}; // signal samples
enum Hypothesis {MM, EM, EE, ALL}; // hypothesis types (em and me counted as same) and all


// this is Jake's magic to sort jets by Pt
Bool_t comparePt(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lv1, 
                 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lv2) {
   return lv1.pt() > lv2.pt();
}

struct DorkyEventIdentifier {
     // this is a workaround for not having unique event id's in MC
     unsigned long int run, event, lumi;
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

// filter events by process
bool filterByProcess( enum Sample sample ) {
  switch (sample) {
  case DYee: 
    return isDYee();
  case DYmm:
    return isDYmm();
  case DYtt:
    return isDYtt();
  case WW:
    return isWW();
  case WZ:
    return isWZ();
  case ZZ:
    return isZZ();
  default:
    return true;
  }
}

bool isIdentified( enum Sample sample ) {
  switch (sample) {
  case DYee:
  case DYmm:
  case DYtt:
    return getDrellYanType()!=999;
  case WW:
  case WZ:
  case ZZ:
    return getVVType()!=999;
  default:
    return true;
  }
}

// filter candidates by hypothesis
Hypothesis filterByHypothesis( int candidate ) {
  switch (candidate) {
  case 0:
    return MM;
  case 1: case 2:
    return EM;
  case 3:
    return EE;
  }
  cout << "Unknown type: " << candidate << "Abort" << endl;
  assert(0);
  return MM;
}


// Meson Classification

bool idIsCharm(int id) {
  id = abs(id);
  if (
      id == 4       ||
      id == 411     ||
      id == 421     ||
      id == 10411   ||
      id == 10421   ||
      id == 413     ||
      id == 423     ||
      id == 10413   ||
      id == 10423   ||
      id == 20413   ||
      id == 20423   ||
      id == 415     ||
      id == 425     ||
      id == 431     ||
      id == 10431   ||
      id == 433     ||
      id == 10433   ||
      id == 20433   ||
      id == 435     ||
      id == 441     ||
      id == 10441   ||
      id == 100441  ||
      id == 443     ||
      id == 10443   ||
      id == 20443   ||
      id == 100443  ||
      id == 30443   ||
      id == 9000443 ||
      id == 9010443 ||
      id == 9020443 ||
      id == 445     ||
      id == 9000445 ||
      id == 4122    ||
      id == 4222    ||
      id == 4212    ||
      id == 4112    ||
      id == 4224    ||
      id == 4214    ||
      id == 4114    ||
      id == 4232    ||
      id == 4132    ||
      id == 4322    ||
      id == 4312    ||
      id == 4324    ||
      id == 4314    ||
      id == 4332    ||
      id == 4334    ||
      id == 4412    ||
      id == 4422    ||
      id == 4414    ||
      id == 4424    ||
      id == 4432    ||
      id == 4434    ||
      id == 4444
      ) {
    return true;
  }
  else return false;
}

bool idIsBeauty(int id) {
  id = abs(id);
  if (
      id == 5       ||
      id == 511     ||
      id == 521     ||
      id == 10511   ||
      id == 10521   ||
      id == 513     ||
      id == 523     ||
      id == 10513   ||
      id == 10523   ||
      id == 20513   ||
      id == 20523   ||
      id == 515     ||
      id == 525     ||
      id == 531     ||
      id == 10531   ||
      id == 533     ||
      id == 10533   ||
      id == 20533   ||
      id == 535     ||
      id == 541     ||
      id == 10541   ||
      id == 543     ||
      id == 10543   ||
      id == 20543   ||
      id == 545     ||
      id == 551     ||
      id == 10551   ||
      id == 100551  ||
      id == 110551  ||
      id == 200551  ||
      id == 210551  ||
      id == 553     ||
      id == 10553   ||
      id == 20553   ||
      id == 30553   ||
      id == 100553  ||
      id == 110553  ||
      id == 120553  ||
      id == 130553  ||
      id == 200553  ||
      id == 210553  ||
      id == 220553  ||
      id == 300553  ||
      id == 9000553 ||
      id == 9010553 ||
      id == 555     ||
      id == 10555   ||
      id == 20555   ||
      id == 100555  ||
      id == 110555  ||
      id == 120555  ||
      id == 200555  ||
      id == 557     ||
      id == 100557  ||
      id == 5122    || 
      id == 5112    ||
      id == 5212    ||
      id == 5222    ||
      id == 5114    ||
      id == 5214    ||
      id == 5224    ||
      id == 5132    ||
      id == 5232    ||
      id == 5312    ||
      id == 5322    ||
      id == 5314    ||
      id == 5324    ||
      id == 5332    ||
      id == 5334    ||
      id == 5142    ||
      id == 5242    ||
      id == 5412    ||
      id == 5422    ||
      id == 5414    ||
      id == 5424    ||
      id == 5342    ||
      id == 5432    ||
      id == 5434    ||
      id == 5442    ||
      id == 5444    ||
      id == 5512    ||
      id == 5522    ||
      id == 5514    ||
      id == 5524    ||
      id == 5532    ||
      id == 5534    ||
      id == 5542    ||
      id == 5544    ||
      id == 5554 
      ) {
    return true;
  }
  else return false;
}



//  Book histograms...
//  Naming Convention:
//  Prefix comes from the sample and it is passed to the scanning function
//  Suffix is "ee" "em" "em" "all" which depends on the final state
//  For example: histogram named tt_hnJet_ee would be the Njet distribution
//  for the ee final state in the ttbar sample.

// MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!

TH1F* hnJet[4];       // Njet distributions
TH1F* hnJetttbarother[4];
TH1F* hnJetttbarlep[4];
TH1F* hnJetfakeCharge[4];       // Njet distributions
TH1F* hnJetSemiTop[4];       // Njet distributions
TH1F* hnJetPlus[4];       // Njet distributions
TH1F* hnJetMinus[4];       // Njet distributions
TH1F* hnJetSemiWTop[4];       // Njet distributions
TH1F* hnJetSemitrueTop[4];
TH1F* hnJetSemiOtherTop[4];
TH1F* hnJetfakeLep[4];       // Njet distributions

TH1F* htcmetZveto[4];

// fkw September 2008 final hist used for muon tags estimate of top bkg


vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > calo_jetsp4;

struct hypo_monitor{
  std::vector<std::pair<std::string,unsigned int> > counters;
  void count(unsigned int index, const char* name){
    unsigned int current_size = counters.size();
    for ( unsigned int i=current_size; i<=index; ++i ) 
      counters.push_back( std::pair<std::string,unsigned int>("",0) );
    counters[index].first = name;
    counters[index].second++;
  }
  void print(){
    for ( unsigned int i=0; i<counters.size(); ++i ) 
      std::cout << counters[i].first << "\t" << counters[i].second << std::endl;
  }
};
    
hypo_monitor monitor;

void checkIsolation(int i_hyp, double weight){
  // LT
  if ( abs(cms2.hyp_lt_id()[i_hyp]) == 11 ) { 

  }
  // LL
  if ( abs(cms2.hyp_ll_id()[i_hyp]) == 11 ) { 

  }
}

void getIsolationSidebandsAfterSelections(int i_hyp, double weight, RooDataSet* dataset, bool passedAllLeptonRequirements){

}

void find_most_energetic_jets(int i_hyp, double weight)
{

}  

void hypo (int i_hyp, double kFactor, RooDataSet* dataset = 0) 
{
  int myType = 99;
  if (cms2.hyp_type()[i_hyp] == 3) myType = 0;  // ee
  if (cms2.hyp_type()[i_hyp] == 0) myType = 1;  // mm
  if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2) myType=2; // em
  if (myType == 99) {
    std::cout << "YUK:  unknown dilepton type = " << cms2.hyp_type()[i_hyp] << std::endl;
    return;
  }

  // The event weight including the kFactor (scaled to 1 fb-1)
     float weight = cms2.evt_scale1fb() * kFactor;

     unsigned int icounter(0);
     monitor.count(icounter++,"Total number of hypothesis: ");
     
     //     if ( ! passTriggersMu9orLisoE15( cms2.hyp_type()[i_hyp] ) ) return;
     if (! GoodSusyTrigger( cms2.hyp_type()[i_hyp] ) ) return;
     monitor.count(icounter++,"Total number of hypothesis after trigger requirements: ");
     
     // Cut on lepton Pt and eta
     if (cms2.hyp_lt_p4()[i_hyp].pt() < 10.0) return;
     if (cms2.hyp_ll_p4()[i_hyp].pt() < 10.0) return;
     if (max(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()) < 20) return;


     bool conversion = false;
     bool mischarge = false;
     
     if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
       int elIndex = cms2.hyp_ll_index()[i_hyp];
       if ( conversionElectron(elIndex)) conversion = true;
       if ((cms2.els_trkidx().at(elIndex) >= 0) && (cms2.els_charge().at(elIndex) != cms2.trks_charge().at(cms2.els_trkidx().at(elIndex)))) mischarge = true;
     }

     if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
       int elIndex = cms2.hyp_lt_index()[i_hyp];
       if ( conversionElectron(elIndex)) conversion = true;
       if ((cms2.els_trkidx().at(elIndex) >= 0) && (cms2.els_charge().at(elIndex) != cms2.trks_charge().at(cms2.els_trkidx().at(elIndex)))) mischarge = true;
     }

     if (conversion) return;
     if (mischarge) return;


     monitor.count(icounter++,"Total number of hypothesis after lepton pt and eta cut: ");

     // Require opposite sign
     //  if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) return;
     // Same Sign
     if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 ) return;

     // check electron isolation and id (no selection at this point)
     // checkIsolation(i_hyp, weight);
     
     monitor.count(icounter++,"Total number of hypothesis after lepton pt + z vetos: ");
     
     bool goodEvent = true;
     bool passedAllLeptonRequirements = true;

     // Lepton Quality cuts and isolation according to VJets09

     if (!GoodSusyLeptonID(cms2.hyp_lt_id()[i_hyp], cms2.hyp_lt_index()[i_hyp])) passedAllLeptonRequirements = false;
     if (!GoodSusyLeptonID(cms2.hyp_ll_id()[i_hyp], cms2.hyp_ll_index()[i_hyp])) passedAllLeptonRequirements = false;
     if (!PassSusyLeptonIsolation(cms2.hyp_ll_id()[i_hyp], cms2.hyp_ll_index()[i_hyp])) passedAllLeptonRequirements = false;
     if (!PassSusyLeptonIsolation(cms2.hyp_lt_id()[i_hyp], cms2.hyp_lt_index()[i_hyp])) passedAllLeptonRequirements = false;

     // Z mass veto using hyp_leptons for ee and mumu final states
     if (cms2.hyp_type()[i_hyp] == 0 || cms2.hyp_type()[i_hyp] == 3) {
       if (inZmassWindow(cms2.hyp_p4()[i_hyp].mass()) || additionalZveto()) {
	 htcmetZveto[myType]->Fill(cms2.evt_tcmet(),weight);
	 htcmetZveto[3]->Fill(cms2.evt_tcmet(),weight);
       }
     }

     bool useTcMet = true;
     if (!passMetVJets09(80., useTcMet)) return;

     monitor.count(icounter++,"Total number of hypothesis after lepton pt + z vetos + MET cuts: ");     

     if ( !passedAllLeptonRequirements ) return;
     if ( additionalZvetoSUSY09(i_hyp)) return;

     monitor.count(icounter++,"Total number of hypothesis after full lepton selection + z vetos + MET cuts: ");
     
     if ( ! goodEvent ) return;

     // -------------------------------------------------------------------//
     // If we made it to here, we passed all cuts and we are ready to fill //
     // -------------------------------------------------------------------//


     calo_jetsp4.clear();

     double etMax_calo = 0.0;
     double sumet_calo = 0.0;
     calo_jetsp4 = getCaloJets(i_hyp);
     //    calo_jetsp4 = getJPTJets(i_hyp);


     for (unsigned int jj=0; jj < calo_jetsp4.size(); ++jj) {
       if (calo_jetsp4[jj].Et() > etMax_calo) etMax_calo = calo_jetsp4[jj].Et();
       sumet_calo += calo_jetsp4[jj].Et();
     }

     int njets = 0;
     if (calo_jetsp4.size() > 0) njets = calo_jetsp4.size();
     // Final cuts on jets
     if (njets < 1) return;
     if (calo_jetsp4[0].Et() < 100) return; 
     //     if (sumet_calo < 200) return;
     
     // Semileptonic top - ala Claudio
     bool chargefake = false;
     bool semilep = false;
     bool ttbarlep = true;
     bool ttbarother = true;
     int sumCh = 99999;

     if ((abs(cms2.hyp_ll_id()[i_hyp]) == abs(cms2.hyp_ll_mc_id()[i_hyp])) || (abs(cms2.hyp_lt_id()[i_hyp]) == abs(cms2.hyp_lt_mc_id()[i_hyp]))) {
       if (cms2.hyp_lt_id()[i_hyp] * cms2.hyp_lt_mc_id()[i_hyp] < 0 ) chargefake = true; 
       if (cms2.hyp_ll_id()[i_hyp] * cms2.hyp_ll_mc_id()[i_hyp] < 0 ) chargefake = true; 
       sumCh = cms2.hyp_lt_id()[i_hyp] + cms2.hyp_ll_id()[i_hyp];
     } 

     if ((cms2.hyp_ll_mc_id()[i_hyp]==22 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.hyp_ll_mc_p4()[i_hyp])) <0.05 
	  && abs(cms2.hyp_ll_id()[i_hyp]) == abs(cms2.hyp_ll_mc_motherid()[i_hyp])) || 
	 (cms2.hyp_lt_mc_id()[i_hyp]==22 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.hyp_lt_mc_p4()[i_hyp])) <0.05
	  && abs(cms2.hyp_lt_id()[i_hyp]) == abs(cms2.hyp_lt_mc_motherid()[i_hyp]))) {
       if (cms2.hyp_lt_id()[i_hyp] * cms2.hyp_lt_mc_motherid()[i_hyp] < 0 ) chargefake = true;  
       if (cms2.hyp_ll_id()[i_hyp] * cms2.hyp_ll_mc_motherid()[i_hyp] < 0 ) chargefake = true;  
       sumCh = cms2.hyp_lt_mc_motherid()[i_hyp] + cms2.hyp_ll_mc_motherid()[i_hyp];
     }

     if (cms2.hyp_lt_id()[i_hyp] * cms2.hyp_lt_mc_id()[i_hyp] < 0 ) chargefake = true;
     if (cms2.hyp_ll_id()[i_hyp] * cms2.hyp_ll_mc_id()[i_hyp] < 0 ) chargefake = true;
     
     if (genpCountPDGId(11,13,15) != 2) ttbarlep = false;
     if (genpCountPDGId(11,13,15) == 2) ttbarother = false;

     { // ll hypothesis test the heavy flavor
       int els_mo = 0; 
       int mus_mo = 0; 
       int els_id = 0;
       int mus_id = 0;
       if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) els_mo = cms2.els_mc3_motherid()[cms2.hyp_ll_index()[i_hyp]]; 
       if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) mus_mo = cms2.mus_mc3_motherid()[cms2.hyp_ll_index()[i_hyp]]; 
       if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) els_id = cms2.els_mc3_id()[cms2.hyp_ll_index()[i_hyp]];
       if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) mus_id = cms2.mus_mc3_id()[cms2.hyp_ll_index()[i_hyp]];
       if (idIsCharm(cms2.hyp_ll_mc_motherid()[i_hyp]) || idIsBeauty(cms2.hyp_ll_mc_motherid()[i_hyp])) semilep = true;
       if (idIsCharm(els_mo) || idIsBeauty(els_mo)) semilep = true;
       if (idIsCharm(mus_mo) || idIsBeauty(mus_mo)) semilep = true;
       if (idIsCharm(els_id) || idIsBeauty(els_id)) semilep = true;
       if (idIsCharm(mus_id) || idIsBeauty(mus_id)) semilep = true;
     }

     { // lt hypothesis test the heavy flavor
       int els_mo = 0; 
       int mus_mo = 0; 
       int els_id = 0;
       int mus_id = 0;
       if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) els_mo = cms2.els_mc3_motherid()[cms2.hyp_lt_index()[i_hyp]]; 
       if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) mus_mo = cms2.mus_mc3_motherid()[cms2.hyp_lt_index()[i_hyp]]; 
       if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) els_id = cms2.els_mc3_id()[cms2.hyp_lt_index()[i_hyp]];
       if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) mus_id = cms2.mus_mc3_id()[cms2.hyp_lt_index()[i_hyp]];
       if (idIsCharm(cms2.hyp_lt_mc_motherid()[i_hyp]) || idIsBeauty(cms2.hyp_lt_mc_motherid()[i_hyp])) semilep = true;
       if (idIsCharm(els_mo) || idIsBeauty(els_mo)) semilep = true;
       if (idIsCharm(mus_mo) || idIsBeauty(mus_mo)) semilep = true;
       if (idIsCharm(els_id) || idIsBeauty(els_id)) semilep = true;
       if (idIsCharm(mus_id) || idIsBeauty(mus_id)) semilep = true;
     }

     int tttype = ttbarconstituents(i_hyp);

     //  Fill the distribution
     hnJet[myType]->Fill(min(njets,4), weight);
     hnJet[3]->Fill(min(njets,4), weight);
     
     if (ttbarlep) {
       hnJetttbarlep[myType]->Fill(min(njets,4), weight);
       hnJetttbarlep[3]->Fill(min(njets,4), weight);
     } else if (ttbarother) { 
       hnJetttbarother[myType]->Fill(min(njets,4), weight);
       hnJetttbarother[3]->Fill(min(njets,4), weight);
     }
     if (chargefake) {
       hnJetfakeCharge[myType]->Fill(min(njets,4), weight);
       hnJetfakeCharge[3]->Fill(min(njets,4), weight);
     }

     if (tttype == 1) {
       hnJetSemiTop[myType]->Fill(min(njets,4), weight);
       hnJetSemiTop[3]->Fill(min(njets,4), weight);
       if (sumCh > 0 && sumCh != 99999) {
         hnJetPlus[myType]->Fill(min(njets,4), weight);
         hnJetPlus[3]->Fill(min(njets,4), weight);
       } else if (sumCh < 0) {
         hnJetMinus[myType]->Fill(min(njets,4), weight);
         hnJetMinus[3]->Fill(min(njets,4), weight);
       }
     } else if (tttype == 2) {
       hnJetSemiWTop[myType]->Fill(min(njets,4), weight);
       hnJetSemiWTop[3]->Fill(min(njets,4), weight);
       if (semilep) {
	 hnJetSemitrueTop[myType]->Fill(min(njets,4), weight); 
	 hnJetSemitrueTop[3]->Fill(min(njets,4), weight); 
	 //	 cout << cms2.hyp_lt_mc_motherid()[i_hyp] << "   " << cms2.hyp_ll_mc_motherid()[i_hyp] << "  " << cms2.hyp_lt_mc_id()[i_hyp] << "  " << cms2.hyp_ll_mc_id()[i_hyp] << endl;
	 //	 if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) cout << " elec " <<  cms2.els_mc3_motherid()[cms2.hyp_ll_index()[i_hyp]] <<  "  " << cms2.els_mc3_id()[cms2.hyp_ll_index()[i_hyp]] << endl;
	 //	 if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) cout << " muon " << cms2.mus_mc3_motherid()[cms2.hyp_ll_index()[i_hyp]] << "  " << cms2.mus_mc3_id()[cms2.hyp_ll_index()[i_hyp]] << endl;
	 //	 if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) cout << " elec lt " <<  cms2.els_mc3_motherid()[cms2.hyp_lt_index()[i_hyp]] <<  "  " << cms2.els_mc3_id()[cms2.hyp_lt_index()[i_hyp]] << endl;
	 //	 if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) cout << " muon lt " << cms2.mus_mc3_motherid()[cms2.hyp_lt_index()[i_hyp]] << "  " << cms2.mus_mc3_id()[cms2.hyp_lt_index()[i_hyp]] << endl;
	 //	 cout << "ps gen " << genpCountPDGId(11,13,15) << endl;
       } else {
	 hnJetSemiOtherTop[myType]->Fill(min(njets,4), weight); 
	 hnJetSemiOtherTop[3]->Fill(min(njets,4), weight); 
       }
     } else {
       hnJetfakeLep[myType]->Fill(min(njets,4), weight);
       hnJetfakeLep[3]->Fill(min(njets,4), weight);
     }

     hypos_total->Fill(myType);
     hypos_total->Fill(3);
     hypos_total_weighted->Fill(myType,weight);
     hypos_total_weighted->Fill(3,weight);

}//end of void hypo

RooDataSet* MakeNewDataset(const char* name)
{
  RooRealVar set_iso("iso","iso",0.,1.);
  RooRealVar set_event("event","event",0);
  RooRealVar set_run("run","run",0);
  RooRealVar set_lumi("lumi","lumi",0);
  RooRealVar set_weight("weight","weight",0);
  RooCategory set_selected("selected","Passed final WW selection requirements");
  set_selected.defineType("true",1);
  set_selected.defineType("false",0);

  RooCategory set_hyp_type("hyp_type","Hypothesis type");
  set_hyp_type.defineType("ee",0);
  set_hyp_type.defineType("mm",1);
  set_hyp_type.defineType("em",2);
  
  RooCategory set_fake_type("fake_type","Define type of lepton for which isolation is extracted");
  set_fake_type.defineType("electron",0);
  set_fake_type.defineType("muon",1);

  RooCategory set_sample_type("sample_type","Sample type");
  set_sample_type.defineType("data_relaxed_iso",0);  // full sample with final selection 
  set_sample_type.defineType("control_sample_signal_iso",1);

  RooDataSet* dataset = new RooDataSet(name, name,
				       RooArgSet(set_event,set_run,set_lumi,
						 set_iso,set_selected,set_weight,
						 set_hyp_type,set_fake_type,set_sample_type) );
  dataset->setWeightVar(set_weight);
  return dataset;
}

void AddIsoSignalControlSample( int i_hyp, double kFactor, RooDataSet* dataset = 0 ) {
  if ( !dataset ) return;
  // The event weight including the kFactor (scaled to 1 fb-1)
  float weight = cms2.evt_scale1fb() * kFactor;
  // Cut on lepton Pt
  if (cms2.hyp_lt_p4()[i_hyp].pt() < 20.0) return;
  if (cms2.hyp_ll_p4()[i_hyp].pt() < 20.0) return;
  // Require opposite sign
  if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) return;
  // Z mass veto using hyp_leptons for ee and mumu final states
  if ( cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2 ) return;
  if (! inZmassWindow(cms2.hyp_p4()[i_hyp].mass())) return;
  RooArgSet set( *(dataset->get()) );
  set.setCatIndex("selected",0);
  set.setRealValue("event",cms2.evt_event());
  set.setRealValue("run",cms2.evt_run());
  set.setRealValue("lumi",cms2.evt_lumiBlock());
  set.setCatLabel("sample_type","control_sample_signal_iso");
	
  if ( cms2.hyp_type()[i_hyp] == 3 ){
    set.setCatLabel("hyp_type","ee");
    set.setCatLabel("fake_type","electron");
    if ( goodElectronIsolated(cms2.hyp_lt_index()[i_hyp],true) &&
	 goodElectronWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ){
      set.setRealValue("iso",el_rel_iso(cms2.hyp_ll_index()[i_hyp],true));
      dataset->add(set,weight);
    }
    if ( goodElectronIsolated(cms2.hyp_ll_index()[i_hyp],true) &&
	 goodElectronWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ){
      set.setRealValue("iso",el_rel_iso(cms2.hyp_lt_index()[i_hyp],true));
      dataset->add(set,weight);
    }
  } else {
    set.setCatLabel("hyp_type","mm");
    set.setCatLabel("fake_type","muon");
    if ( goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) &&
	 goodMuonWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ){
      set.setRealValue("iso",mu_rel_iso(cms2.hyp_ll_index()[i_hyp]));
      dataset->add(set,weight);
    }
    if ( goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) &&
	 goodMuonWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ){
      set.setRealValue("iso",mu_rel_iso(cms2.hyp_lt_index()[i_hyp]));
      dataset->add(set,weight);
    }
  }
}

RooDataSet* ScanChain( TChain* chain, enum Sample sample, bool identifyEvents ) {
  
  unsigned int nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
  unsigned int nEventsTotal = 0;

  // const unsigned int numHypTypes = 4;  // number of hypotheses: MM, EM, EE, ALL

 // declare and create array of histograms
  //  const char sample_names[][1024] = { "ww", "wz", "zz", "wjets", "dyee", "dymm", "dytt", "ttbar", "tw", "lm0x", "lm1x", "lm2x", "lm3x", "lm4x", "lm5x", "lm6x", "lm7x", "lm8x", "lm9x" };
  const char sample_names[][1024] = { "ww", "wz", "zz", "wjets", "dy", "dyee", "dymm", "dytt", "ttbar", "tw", "lm0x", "lm1x", "lm2x", "lm3x", "lm4x", "lm5x", "lm6x", "lm7x", "lm8x", "lm9x"};
  const char *prefix = sample_names[sample];
  RooDataSet* dataset = MakeNewDataset(sample_names[sample]);
  double kFactor = 1; // 1fb-1

  //  double kFactor = .1; // 1fb-1
  //   switch (sample) {
  //   case WW:
  //        evt_scale1fb = 0.1538;
  //        break;
  //   default:
  //        break;

  char *jetbins[5] = {"0", "1", "2", "3", "#geq 4"};
  char *suffix[3];
  suffix[0] = "ee";
  suffix[1] = "mm";
  suffix[2] = "em";
  suffix[3] = "all";
  
  hypos_total          = new TH1F(Form("%s_hypos_total",prefix),"Total number of hypothesis counts",4,0,4);
  hypos_total_weighted = new TH1F(Form("%s_hypos_total_weighted",prefix),"Total number of hypotheses (weighted)",4,0,4);
  hypos_total_weighted->Sumw2();

  for (unsigned int i=0; i<4; ++i){
    hypos_total->GetXaxis()->SetBinLabel(i+1,suffix[i]);
    hypos_total_weighted->GetXaxis()->SetBinLabel(i+1,suffix[i]);
  }
  
  // The statement below should work but does not work due to bug in root when TH2 are also used
  // Rene Brun promised a fix.
  //TH1::SetDefaultSumw2(kTRUE); // do errors properly based on weights
  

  for (int i=0; i<4; i++) {

    hnJet[i] = new TH1F(Form("%s_hnJet_%s",prefix,suffix[i]),Form("%s_nJet_%s",prefix,suffix[i]),
			5,0.,5.);	

    for(int k = 0; k<5; k++) {
      hnJet[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
      hnJet[i]->GetXaxis()->SetLabelSize(0.07);
    }
    hnJetttbarother[i] = new TH1F(Form("%s_hnJetttbarother_%s",prefix,suffix[i]),Form("%s_hnJetttbarother_%s",prefix,suffix[i]),
			5,0.,5.);	
    hnJetttbarlep[i] = new TH1F(Form("%s_hnJetttbarlep_%s",prefix,suffix[i]),Form("%s_hnJetttbarlep_%s",prefix,suffix[i]),
			5,0.,5.);	
    hnJetfakeCharge[i] = new TH1F(Form("%s_hnJetfakeCharge_%s",prefix,suffix[i]),Form("%s_hnJetfakeCharge_%s",prefix,suffix[i]),
			5,0.,5.);	
    hnJetSemiWTop[i] = new TH1F(Form("%s_hnJetSemiWTop_%s",prefix,suffix[i]),Form("%s_hnJetSemiWTop_%s",prefix,suffix[i]),
			5,0.,5.);	    
    hnJetPlus[i] = new TH1F(Form("%s_hnJetPlus_%s",prefix,suffix[i]),Form("%s_hnJetPlus_%s",prefix,suffix[i]),
			5,0.,5.);	
    hnJetMinus[i] = new TH1F(Form("%s_hnJetMinus_%s",prefix,suffix[i]),Form("%s_hnJetMinus_%s",prefix,suffix[i]),
			5,0.,5.);	
    hnJetSemiTop[i] = new TH1F(Form("%s_hnJetSemiTop_%s",prefix,suffix[i]),Form("%s_hnJetSemiTop_%s",prefix,suffix[i]),
			5,0.,5.);	
    hnJetSemitrueTop[i] = new TH1F(Form("%s_hnJetSemitrueTop_%s",prefix,suffix[i]),Form("%s_hnJetSemitrueTop_%s",prefix,suffix[i]),
			5,0.,5.);	
    hnJetSemiOtherTop[i] = new TH1F(Form("%s_hnJetSemiOtherTop_%s",prefix,suffix[i]),Form("%s_hnJetSemiOtherTop_%s",prefix,suffix[i]),
			5,0.,5.);	
    hnJetfakeLep[i] = new TH1F(Form("%s_hnJetfakeLep_%s",prefix,suffix[i]),Form("%s_hnJetfakeLep_%s",prefix,suffix[i]),
			5,0.,5.);	
    htcmetZveto[i] = new TH1F(Form("%s_htcmetZveto_%s",prefix,suffix[i]),Form("%s_htcmetZveto_%s",prefix,suffix[i]),100,0.,200.);


    hnJet[i]->Sumw2();

  }
  
  // clear list of duplicates
  already_seen.clear();
  int duplicates_total_n = 0;
  double duplicates_total_weight = 0;
  int nFailedIdentification = 0;
  int nFilteredOut = 0;
  int i_permille_old = 0;
  // file loop
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  monitor.counters.clear();

  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
    // need to call TFile::Open(), since the file is not
    // necessarily a plain TFile (TNetFile, TDcacheFile, etc)
    //        printf("current file: %s (%s), %s\n", currentFile->GetName(), 
    // 	      currentFile->GetTitle(), currentFile->IsA()->GetName());

       TFile *f = TFile::Open(currentFile->GetTitle()); 
       assert(f);
       TTree *tree = (TTree*)f->Get("Events");
       
       cms2.Init(tree);  // set branch addresses for TTree tree

       TStopwatch t;
       //Event Loop

       unsigned int nEvents = tree->GetEntries();
       for( unsigned int event = 0; event < nEvents; ++event) {
	    cms2.GetEntry(event);  // get entries for Event number event from branches of TTree tree
	    ++nEventsTotal;
//  Loops
	    if (cms2.trks_d0().size() == 0)
	      continue;
	    
	    DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.evt_lumiBlock(), cms2.trks_d0()[0], 
					cms2.trks_trk_p4()[0].pt(), cms2.trks_trk_p4()[0].eta(), cms2.trks_trk_p4()[0].phi() };


	    if (is_duplicate(id)) {
	      duplicates_total_n++;
	      duplicates_total_weight += cms2.evt_scale1fb();
		 // cout << "Duplicate event found. Run: " << cms2.evt_run() << ", Event:" << cms2.evt_event() << ", Lumi: " << cms2.evt_lumiBlock() << endl;
		 continue;
	    }

	    int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
	    if (i_permille != i_permille_old) {
		 // xterm magic from L. Vacavant and A. Cerri
		 printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
			"\033[0m\033[32m <---\033[0m\015", i_permille/10.);
		 fflush(stdout);
		 i_permille_old = i_permille;
	    }
	    
	    if ( identifyEvents ){
	      // check if we know what we are looking at
	      if ( ! isIdentified(sample) ) nFailedIdentification++;
	      
	      // filter by process
	      if ( ! filterByProcess(sample) ) {
		nFilteredOut++;
		continue;
	      }
	    }

	    // fkw, here go all histos that should be filled per event instead of per hyp:
	    // loop over generator particles:
	    //Note: top = +-6, W = +-24, b = +-5
	    //cout << " Event = " << event << endl;
	    // fkw, end of per event filling of histos.
	    // loop over hypothesis candidates
 
	    unsigned int nHyps = cms2.hyp_type().size();
	    for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
	      hypo(i_hyp, kFactor, dataset);
	      AddIsoSignalControlSample(i_hyp, kFactor, dataset);
	    }
       }
       t.Stop();
       printf("Finished processing file: %s\n",currentFile->GetTitle());
       printf("Real time: %u events / %f s = %e event/s\n", nEvents, 
	      t.RealTime(), nEvents / t.RealTime());
       printf("CPU time: %u events / %f s = %e event/s\n", nEvents, 
	      t.CpuTime(), nEvents / t.CpuTime());
       printf("Total duplicate count: %d.  Total weight %f\n",   
	      duplicates_total_n, duplicates_total_weight);
       delete f;
  }


  monitor.print();
  if ( nEventsChain != nEventsTotal ) {
       printf("ERROR: number of events from files (%d) is not equal to total number"
	      " of events (%d)\n", nEventsChain, nEventsTotal);
  }
  printf("Total number of skipped events due to bad identification: %d (%0.0f %%)\n",   
	   nFailedIdentification, nFailedIdentification*100.0/(nEventsChain+1e-5));
  printf("Total number of filtered out events: %d (%0.0f %%)\n",   
	   nFilteredOut, nFilteredOut*100.0/(nEventsChain+1e-5));
  printf("Total candidate count (ee mm em all): %.0f %.0f %.0f %0.f.\n",
	 hypos_total->GetBinContent(1), hypos_total->GetBinContent(2), 
	 hypos_total->GetBinContent(3), hypos_total->GetBinContent(4));
  printf("Total weighted candidate yeild (ee mm em all): %f %f %f %f\n",   
	 hypos_total_weighted->GetBinContent(1), hypos_total_weighted->GetBinContent(2), 
	 hypos_total_weighted->GetBinContent(3), hypos_total_weighted->GetBinContent(4));

//  ofstream outf;
//  outf.open("susy_eventcounts.txt", ios::app);
//  outf<<"|" << prefix << "|" << nEventsTotal  << endl;
//  outf.close();

  return dataset;
}

