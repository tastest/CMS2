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
#include "../Tools/tools.cc"
#endif

TH1F*hypos_total;
TH1F* hypos_total_weighted;

enum Sample {WW, WZ, ZZ, Wjets, DY, DYee, DYmm, DYtt, ttbar, tW, LM0x, LM1x, LM2x, LM3x, LM4x, LM5x, LM6x, LM7x, LM8x, LM9x}; // signal samples
enum Hypothesis {MM, EM, EE, ALL}; // hypothesis types (em and me counted as same) and all

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

//  Book histograms...
//  Naming Convention:
//  Prefix comes from the sample and it is passed to the scanning function
//  Suffix is "ee" "em" "em" "all" which depends on the final state
//  For example: histogram named tt_hnJet_ee would be the Njet distribution
//  for the ee final state in the ttbar sample.

// MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!

TH1F* hnJet[4];       // Njet distributions
TH1F* hnJetWW[4];
TH1F* hnJetWO[4];
TH1F* hnJetWOSemilep[4];       // Njet distributions
TH1F* hnJetWOOther[4];       // Njet distributions
TH1F* hnJetOO[4];       // Njet distributions

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

     monitor.count(icounter++,"Total number of hypothesis after adding lepton pt cut: ");

     bool conversion = false;
     bool mischarge = false;
     
     if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
       int elIndex = cms2.hyp_ll_index()[i_hyp];
       if ( conversionElectron(elIndex)) conversion = true;
       if ( isChargeFlip(elIndex)) mischarge = true;
     }

     if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
       int elIndex = cms2.hyp_lt_index()[i_hyp];
       if ( conversionElectron(elIndex)) conversion = true;
       if ( isChargeFlip(elIndex)) mischarge = true;
     }

     if (conversion) return;
     if (mischarge) return;

     monitor.count(icounter++,"Total number of hypothesis after adding conversion and chargeflip cuts: ");

     // Require opposite sign
     //  if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) return;
     // Same Sign
     if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 ) return;

     monitor.count(icounter++,"Total number of hypothesis after adding charge cut: ");
     
     bool goodEvent = true;
     bool passedAllLeptonRequirements = true;

     // Lepton Quality cuts and isolation according to VJets09

     if (!GoodSusyLeptonID(cms2.hyp_lt_id()[i_hyp], cms2.hyp_lt_index()[i_hyp])) passedAllLeptonRequirements = false;
     if (!GoodSusyLeptonID(cms2.hyp_ll_id()[i_hyp], cms2.hyp_ll_index()[i_hyp])) passedAllLeptonRequirements = false;
     if (!PassSusyLeptonIsolation(cms2.hyp_ll_id()[i_hyp], cms2.hyp_ll_index()[i_hyp])) passedAllLeptonRequirements = false;
     if (!PassSusyLeptonIsolation(cms2.hyp_lt_id()[i_hyp], cms2.hyp_lt_index()[i_hyp])) passedAllLeptonRequirements = false;

     if ( !passedAllLeptonRequirements ) return;
     monitor.count(icounter++,"Total number of hypothesis after adding lepton id and isolation, including eta cuts: ");     

     // Z mass veto using hyp_leptons for ee and mumu final states
     if (cms2.hyp_type()[i_hyp] == 0 || cms2.hyp_type()[i_hyp] == 3) {
       if (inZmassWindow(cms2.hyp_p4()[i_hyp].mass()) || additionalZveto()) {
	 htcmetZveto[myType]->Fill(cms2.evt_tcmet(),weight);
	 htcmetZveto[3]->Fill(cms2.evt_tcmet(),weight);
       }
     }

     bool useTcMet = true;
     if (!passMetVJets09(80., useTcMet)) return;

     monitor.count(icounter++,"Total number of hypothesis after adding tcmet cut: ");     

     calo_jetsp4.clear();

     double etMax_calo = 0.0;
     double sumet_calo = 0.0;
     calo_jetsp4 = getCaloJets(i_hyp);
     //    calo_jetsp4 = getJPTJets(i_hyp);

     for (unsigned int jj=0; jj < calo_jetsp4.size(); ++jj) {
       if (calo_jetsp4[jj].pt() > etMax_calo) etMax_calo = calo_jetsp4[jj].pt();
       sumet_calo += calo_jetsp4[jj].pt();
     }

     int njets = 0;
     if (calo_jetsp4.size() > 0) njets = calo_jetsp4.size();
     // Final cuts on jets
     if (njets < 3) return;
     //if (calo_jetsp4[0].pt() < 100) return; 
     if (sumet_calo < 200) return;
     
     monitor.count(icounter++,"Total number of hypothesis after adding jet cuts: ");

     if ( additionalZvetoSUSY09(i_hyp)) return;

     monitor.count(icounter++,"Total number of hypothesis after adding additionalZvetoSUSY09 cut: ");
     
     if ( ! goodEvent ) return;

     // -------------------------------------------------------------------//
     // If we made it to here, we passed all cuts and we are ready to fill //
     // -------------------------------------------------------------------//

     //  Fill the distribution
     // This hist is used in doTable.C to fill the results table for the twiki.
     hnJet[myType]->Fill(min(njets,4), weight);
     hnJet[3]->Fill(min(njets,4), weight);
     
     //The only thing left to do now is to work out the composition of the ttbar bkg:

     //To distinguish WW (=1), WO (=2), and OO (=3) 
     int tttype = ttbarconstituents(i_hyp);

     // Semileptonic
     int lttype = leptonIsFromW(cms2.hyp_lt_index()[i_hyp],cms2.hyp_lt_id()[i_hyp],cms2.hyp_lt_p4()[i_hyp] );
     int lltype = leptonIsFromW(cms2.hyp_ll_index()[i_hyp],cms2.hyp_ll_id()[i_hyp],cms2.hyp_ll_p4()[i_hyp] );

     if (tttype == 1) {
       hnJetWW[myType]->Fill(min(njets,4), weight);
       hnJetWW[3]->Fill(min(njets,4), weight);
     } else if (tttype == 2) {
       hnJetWO[myType]->Fill(min(njets,4), weight);
       hnJetWO[3]->Fill(min(njets,4), weight);
       if ( lttype == -1 || lttype == -2 || lltype == -1 || lltype == -2) {
	 hnJetWOSemilep[myType]->Fill(min(njets,4), weight); 
	 hnJetWOSemilep[3]->Fill(min(njets,4), weight); 
       } else {
	 hnJetWOOther[myType]->Fill(min(njets,4), weight); 
	 hnJetWOOther[3]->Fill(min(njets,4), weight); 
       }
     } else {
       hnJetOO[myType]->Fill(min(njets,4), weight);
       hnJetOO[3]->Fill(min(njets,4), weight);
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
    hnJetWW[i] = new TH1F(Form("%s_hnJetWW_%s",prefix,suffix[i]),Form("%s_hnJetWW_%s",prefix,suffix[i]),
			5,0.,5.);	
    hnJetWO[i] = new TH1F(Form("%s_hnJetWO_%s",prefix,suffix[i]),Form("%s_hnJetWO_%s",prefix,suffix[i]),
			5,0.,5.);	
    hnJetWOSemilep[i] = new TH1F(Form("%s_hnJetWOSemilep_%s",prefix,suffix[i]),Form("%s_hnJetWOSemilep_%s",prefix,suffix[i]),
			5,0.,5.);	
    hnJetWOOther[i] = new TH1F(Form("%s_hnJetWOOther_%s",prefix,suffix[i]),Form("%s_hnJetWOOther_%s",prefix,suffix[i]),
			5,0.,5.);	    
    hnJetOO[i] = new TH1F(Form("%s_hnJetOO_%s",prefix,suffix[i]),Form("%s_hnJetOO_%s",prefix,suffix[i]),
			5,0.,5.);	
    htcmetZveto[i] = new TH1F(Form("%s_htcmetZveto_%s",prefix,suffix[i]),Form("%s_htcmetZveto_%s",prefix,suffix[i]),100,0.,200.);


    hnJet[i]->Sumw2();
    hnJetWW[i]->Sumw2();
    hnJetWO[i]->Sumw2();
    hnJetWOSemilep[i]->Sumw2();
    hnJetWOOther[i]->Sumw2();
    hnJetOO[i]->Sumw2();
    htcmetZveto[i]->Sumw2();

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
	    
	    DorkyEventIdentifier id(cms2);

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

