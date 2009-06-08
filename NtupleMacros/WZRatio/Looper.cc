#include <math.h>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

Looper::Looper (Sample s, cuts_t c, const char *fname) 
  : LooperBase(s, c, fname)
{
  // zero out the candidate counters (don't comment this out)
  memset(cands_passing_	, 0, sizeof(cands_passing_       ));
  memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
  memset(cands_count_	, 0, sizeof(cands_count_         ));
}

void Looper::BookHistos ()
{
  //------------------------------------------------------------
  // Example histo booking; edit for your application
  //------------------------------------------------------------

  // book histograms the manual way:
  for (int i = 0; i < 4; ++i) {
    helPt_[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), "elPt", dilepton_hypo_names[i]), ";el pt", 100, 0, 100);
    hmuPt_[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), "muPt", dilepton_hypo_names[i]), ";mu pt", 100, 0, 100);
    hCaloEtaPt_[i] = new TH2F(Form("%s_%s_%s", SampleName().c_str(), "CaloEtaPt", dilepton_hypo_names[i]), ";pt;eta", 100, 0, 100, 10, -5, 5);

    // call Sumw2 on all histograms
    helPt_[i]->Sumw2();
    hmuPt_[i]->Sumw2();
    hCaloEtaPt_[i]->Sumw2();

    // set histogram color according to definitions in Tools/Samples.cc
    helPt_[i]->SetFillColor(sample_.histo_color);
    helPt_[i]->SetLineColor(sample_.histo_color);
    hmuPt_[i]->SetFillColor(sample_.histo_color);
    hmuPt_[i]->SetLineColor(sample_.histo_color);
    hCaloEtaPt_[i]->SetFillColor(sample_.histo_color);
    hCaloEtaPt_[i]->SetLineColor(sample_.histo_color);

  }
  // or use the N - 1 technology (see NMinus1Hist.h)
  // arguments are as follows: sample, name, binning, required cuts, cuts that are relaxed for the N - 1 plot
  // for the lt N - 1 plot, we relax the lt pt requirement
  hmt_		= new NMinus1Hist(sample_, "mt"   ,	 60, 0, 120, cuts_, CUT_BIT(CUT_MT));
  // same for ll
  hmet_		= new NMinus1Hist(sample_, "met"   ,	 60, 0, 120, cuts_, CUT_BIT(CUT_MET) | CUT_BIT(CUT_ANTI_MET));
  // for the dilepton mass plot, we relax any cut to do with the Z 
  hdilMass_	= new NMinus1Hist(sample_, "dilMass",	 100, 0, 300, cuts_, CUT_BIT(CUT_ZMASS) | CUT_BIT(CUT_ANTI_ZMASS));
  // for the dilepton mass plot, we relax any cut to do with the Z 
  hLepMetMass_	= new NMinus1Hist(sample_, "LepMetMass", 100, 0, 300, cuts_, CUT_BIT(CUT_ZMASS) | CUT_BIT(CUT_ANTI_ZMASS));
  // eta distribution of "lost" lepton from drell yan decays
  hGenLepEta_	= new NMinus1Hist(sample_, "GenLepEta", 100, -5, 5, cuts_, CUT_BIT(CUT_MET) | CUT_BIT(CUT_ANTI_MET));
  // pt distribution of "lost" lepton from drell yan decays
  hGenLepPt_	= new NMinus1Hist(sample_, "GenLepPt", 100, 0, 100, cuts_, CUT_BIT(CUT_MET) | CUT_BIT(CUT_ANTI_MET));
}


bool Looper::FilterEvent()
{ 

  //
  // duplicate filter, based on trk information and dilepton hyp
  //
  // comment in following lines
  // 

  //   if (cms2.trks_d0().size() == 0)
  //     return true;
  //   DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.trks_d0()[0], 
  // 			      cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
  //   if (is_duplicate(id)) {
  //     duplicates_total_n_++;
  //     duplicates_total_weight_ += cms2.evt_scale1fb();
  //     cout << "Filtered duplicate run: " << cms2.evt_run() << " event: " << cms2.evt_event() << endl;
  //     return true;
  //   }


  return false; 
}



cuts_t Looper::EventSelect ()
{
  //------------------------------------------------------------
  // In an event-based analysis, you would make your cuts here
  //------------------------------------------------------------

  cuts_t ret = 0;

  const double weight = Weight(0);

  wmt_   = 0;
  boson  = LorentzVector (0, 0, 0, 0);
  wp4_   = LorentzVector(0, 0, 0, 0);
  genp4_ = LorentzVector(0, 0, 0, 0);

  // enough tracks?
  if (cms2.trks_trk_p4().size() > 2)
    ret |= CUT_BIT(CUT_MORE_THAN_TWO_TRACKS);
  
  if( cms2.evt_tcmet() > 20 )
    ret |= CUT_BIT(CUT_MET);

  if( cms2.evt_tcmet() < 30 )
    ret |= CUT_BIT(CUT_ANTI_MET);

  // loop over leptons and find candidates

  std::vector<int> lep_idx;

  for (unsigned int mus = 0; mus < cms2.mus_p4().size(); mus++) {

    if( cms2.mus_p4()[mus].pt() < 20 ) continue;

    if( !goodMuonIsolated(mus) ) continue;

    lep_idx.push_back( mus + mu_shift);

  }

  for (unsigned int els = 0; els < cms2.els_p4().size(); els++) {

    if( cms2.els_p4()[els].pt() < 20 ) continue;

    if( !goodElectronIsolated(els, true) ) continue;

    lep_idx.push_back( els );

  }

  // skip event if no good leptons found
  if( lep_idx.size() == 0 )
    ret |= CUT_BIT(CUT_NL);

  // skip event if more than 2 good leptons found
  if( lep_idx.size() > 2 )
    ret |= CUT_BIT(CUT_ML);

  // classify as DY, check if passes cuts on OS, SF, ZMASS, MET VETO
  if( lep_idx.size() == 2) {

    if( (lep_idx[0] + lep_idx[1]) / mu_shift == 0 ) {
      ret |= CUT_BIT(CUT_EE);
    
      if( cms2.els_charge()[ lep_idx[0] ] * cms2.els_charge()[ lep_idx[1] ] < 0 )
	ret |= CUT_BIT(CUT_OS);

      lep1 = cms2.els_p4()[ lep_idx[0] ];
      lep2 = cms2.els_p4()[ lep_idx[1] ];
      boson = lep1 + lep2;

      if( boson.M() > 76 && boson.M() < 106 )
	ret |= CUT_BIT(CUT_ZMASS);
    }

    else if( (lep_idx[0] + lep_idx[1]) / mu_shift == 1 ) {
      ret |= CUT_BIT(CUT_EM);

      if( cms2.mus_charge()[ lep_idx[0] % mu_shift ] * cms2.els_charge()[ lep_idx[1] % mu_shift ] < 0 )
	ret |= CUT_BIT(CUT_OS);
    }

    else if( (lep_idx[0] + lep_idx[1]) / mu_shift == 2 ) {
      ret |= CUT_BIT(CUT_MM);
    
      if( cms2.mus_charge()[ lep_idx[0] % mu_shift ] * cms2.mus_charge()[ lep_idx[1] % mu_shift ] < 0 )
	ret |= CUT_BIT(CUT_OS);

      lep1 = cms2.mus_p4()[ lep_idx[0] % mu_shift ];
      lep2 = cms2.mus_p4()[ lep_idx[1] % mu_shift ];
      boson = lep1 + lep2;
      
      if( boson.M() > 76 && boson.M() < 106 )
	ret |= CUT_BIT(CUT_ZMASS);
    }    
  }

  if( lep_idx.size() == 1 ) {

    if( lep_idx[0] / mu_shift ) {
      ret |= CUT_BIT(CUT_M);

      lep1 = cms2.mus_p4()[ lep_idx[0] % mu_shift ];

      wmt_ = TMath::Sqrt( 2. * lep1.pt() * cms2.evt_tcmet() * (1 - cos(cms2.evt_tcmetPhi() - lep1.phi() )));
      
      if( wmt_ > 40 && wmt_ < 100 )
	ret |= CUT_BIT(CUT_MT);

      for( unsigned int trks = 0; trks < cms2.trks_trk_p4().size(); trks++ ) {

	LorentzVector trkp4_ = cms2.trks_trk_p4()[trks];

	wp4_ = lep1 + trkp4_;

	if( wp4_.M() < 76 || wp4_.M() > 106 )
	  ret |= CUT_BIT(CUT_ANTI_ZMASS);
      }

      for( unsigned int gens = 0; gens < cms2.genps_p4().size(); gens++ ) {
	if( cms2.genps_status()[gens] != 3) continue;

	if( abs( cms2.genps_id()[gens] ) != 13 ) continue;

	if( cms2.mus_mcidx()[ lep_idx[0] % mu_shift ] != gens ) genp4_ = cms2.genps_p4()[gens];
      }
    }

    else if( !(lep_idx[0] / mu_shift) ) {
      ret |= CUT_BIT(CUT_E);

      lep1 = cms2.els_p4()[ lep_idx[0] ];

      wmt_ = TMath::Sqrt( 2. * lep1.pt() * cms2.evt_tcmet() * (1 - cos(cms2.evt_tcmetPhi() - lep1.phi() )));
      
      if( wmt_ > 40 && wmt_ < 100 )
	ret |= CUT_BIT(CUT_MT);

      for( unsigned int trks = 0; trks < cms2.trks_trk_p4().size(); trks++ ) {

	LorentzVector trkp4_ = cms2.trks_trk_p4()[trks];

	wp4_ = lep1 + trkp4_;

	if( wp4_.M() < 76 || wp4_.M() > 106 )
	  ret |= CUT_BIT(CUT_ANTI_ZMASS);
      }

      for( unsigned int gens = 0; gens < cms2.genps_p4().size(); gens++ ) {
	if( cms2.genps_status()[gens] != 3) continue;

	if( abs( cms2.genps_id()[gens] ) != 11 ) continue;

	if( cms2.els_mcidx()[ lep_idx[0] % mu_shift ] != gens ) genp4_ = cms2.genps_p4()[gens];
      }
    }
  }

  hmt_->Fill(ret, DILEPTON_ALL, wmt_, weight);
  hmet_->Fill(ret, DILEPTON_ALL, cms2.evt_tcmet(), weight);
  hdilMass_->Fill(ret, DILEPTON_ALL, boson.M(), weight);
  hGenLepEta_->Fill(ret, DILEPTON_ALL, genp4_.eta(), weight);
  hGenLepPt_->Fill(ret, DILEPTON_ALL, genp4_.pt(), weight);

  return ret;
}

cuts_t Looper::DilepSelect (int i_hyp)
{
  //------------------------------------------------------------
  // Example dilepton cuts; edit for your application
  //------------------------------------------------------------

  // cuts are failed until proven otherwise
  cuts_t ret = 0;

  // the return value gets cached, too
  return ret;
}

cuts_t Looper::TrilepSelect (int i_hyp)
{
  //------------------------------------------------------------
  // In a trilepton analysis, you would make your cuts here
  //------------------------------------------------------------

  cuts_t ret = 0;
  return ret;
}

cuts_t Looper::QuadlepSelect (int i_hyp)
{
  //------------------------------------------------------------
  // In a quadlepton analysis, you would make your cuts here
  //------------------------------------------------------------

  cuts_t ret = 0;
  return ret;
}

void Looper::FillEventHistos ()
{
  //------------------------------------------------------------
  // In an event-based analysis, you would fill your histos here
  //------------------------------------------------------------

  cuts_t cuts_passed = EventSelect();

  if((cuts_passed & cuts_) == cuts_) {
    
    const double weight = Weight(0);

    cands_passing_[0] += weight;
    cands_passing_w2_[0] += weight * weight;
  }
}

void Looper::FillDilepHistos (int i_hyp)
{
}

void Looper::FillTrilepHistos (int i_hyp)
{
  //------------------------------------------------------------
  // In a trilepton analysis, you would fill your histos here
  //------------------------------------------------------------
}

void Looper::FillQuadlepHistos (int i_hyp)
{
  //------------------------------------------------------------
  // In a quadlepton analysis, you would fill your histos here
  //------------------------------------------------------------
}

void Looper::End ()
{
  //------------------------------------------------------------
  //Example status message at the end of a looper; edit for your
  //application
  //------------------------------------------------------------

  int ret = fprintf(logfile_, 
		    "Sample %10s: Total candidate count (ee em mm all): %8u %8u %8u %8u."
		    " Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n",   
		    sample_.name.c_str(),
		    CandsCount(DILEPTON_EE), CandsCount(DILEPTON_EMU), CandsCount(DILEPTON_MUMU), CandsCount(DILEPTON_ALL), 
		    CandsPassing(DILEPTON_EE)  , RMS(DILEPTON_EE),  
		    CandsPassing(DILEPTON_EMU) , RMS(DILEPTON_EMU),  
		    CandsPassing(DILEPTON_MUMU), RMS(DILEPTON_MUMU), 
		    CandsPassing(DILEPTON_ALL) , RMS(DILEPTON_ALL));
  if (ret < 0)
    perror("writing to log file");
}
