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
  printevt = false;
  sumjetpt = 0;
  sjpbin = 0;
  counttt = 0;

  // zero out the candidate counters (don't comment this out)
  memset(cands_passing_	, 0, sizeof(cands_passing_       ));
  memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
  memset(cands_count_		, 0, sizeof(cands_count_         ));
}

void Looper::NewHist(TH1F*& h, char* name, char* title, int bins, double min, double max) {
  h = new TH1F(name, title, bins, min, max);
  h->Sumw2();
  h->SetFillColor(sample_.histo_color);
  h->SetLineColor(sample_.histo_color);
  //h.Sumw2();
  //h.SetFillColor(sample_.histo_color);
  //h.SetLineColor(sample_.histo_color);
}

//void Looper::SetupHist(TH1F* h) { //, char* name, char* title, int bins, double min, double max) {
//  //h = new TH1F(name, title, bins, min, max);
//  h->Sumw2();
//  h->SetFillColor(sample_.histo_color);
//  h->SetLineColor(sample_.histo_color);
//}

void Looper::BookHistos ()
{
  double sjpmax = 400;
  int sjpbin = 80; //bins of 5
  
  double tcmcmax = 200; //will be +/- this
  int tcmcbin = 80; // width = nbins*binwidth
  int tcmbin = 40;
 
  // book histograms the manual way:
  for (int i = 0; i < 4; ++i) {
	//hsumjetpt[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), "sumJetPt", dilepton_hypo_names[i]), ";sum jet pt", sjpbin, 0, sjpmax);
	NewHist(htcmet[i],    Form("%s_%s_%s", SampleName().c_str(), "tcMet",      dilepton_hypo_names[i]), ";tcMET",      tcmbin, 0, tcmcmax);
	NewHist(htcmetouz[i], Form("%s_%s_%s", SampleName().c_str(), "tcMet_outz", dilepton_hypo_names[i]), ";tcMET_outz", tcmbin, 0, tcmcmax);
	NewHist(htcmetinz[i], Form("%s_%s_%s", SampleName().c_str(), "tcMet_inz",  dilepton_hypo_names[i]), ";tcMET_inz",  tcmbin, 0, tcmcmax);

	//combined x component + y component
	NewHist(htcmetxy[i],    Form("%s_%s_%s", SampleName().c_str(), "tcMet_xy", dilepton_hypo_names[i]), ";tcMET_x+y", tcmcbin, -tcmcmax, tcmcmax);
	NewHist(htcmetouzxy[i], Form("%s_%s_%s", SampleName().c_str(), "tcMet_xy_outz", dilepton_hypo_names[i]), ";tcMET_x+y_outz", tcmcbin, -tcmcmax, tcmcmax);
	NewHist(htcmetinzxy[i], Form("%s_%s_%s", SampleName().c_str(), "tcMet_xy_inz", dilepton_hypo_names[i]), ";tcMET_x+y_inz", tcmcbin, -tcmcmax, tcmcmax);  
	//individual components
	NewHist(htcmetx[i],    Form("%s_%s_%s", SampleName().c_str(), "tcMet_x", dilepton_hypo_names[i]), ";tcMET_x", tcmcbin, -tcmcmax, tcmcmax);
	NewHist(htcmetouzx[i], Form("%s_%s_%s", SampleName().c_str(), "tcMet_x_outz", dilepton_hypo_names[i]), ";tcMET_x_outz", tcmcbin, -tcmcmax, tcmcmax);
	NewHist(htcmetinzx[i], Form("%s_%s_%s", SampleName().c_str(), "tcMet_x_inz", dilepton_hypo_names[i]), ";tcMET_x_inz", tcmcbin, -tcmcmax, tcmcmax);

	NewHist(htcmety[i],    Form("%s_%s_%s", SampleName().c_str(), "tcMet_y", dilepton_hypo_names[i]), ";tcMET_y", tcmcbin, -tcmcmax, tcmcmax);
	NewHist(htcmetouzy[i], Form("%s_%s_%s", SampleName().c_str(), "tcMet_y_outz", dilepton_hypo_names[i]), ";tcMET_y_outz", tcmcbin, -tcmcmax, tcmcmax);
	NewHist(htcmetinzy[i], Form("%s_%s_%s", SampleName().c_str(), "tcMet_y_inz", dilepton_hypo_names[i]), ";tcMET_y_inz", tcmcbin, -tcmcmax, tcmcmax);

	for(int j=0;j<nsjpbins;j++) {
	  NewHist(htcmetxy_sjp[i][j],    Form("%s_%s%i_%s", SampleName().c_str(), "tcMet_xy_sjp", j, dilepton_hypo_names[i]),
			  Form("%s%i", ";tcMET_x+y_sjp", j), tcmcbin, -tcmcmax, tcmcmax);	
	  NewHist(htcmetouzxy_sjp[i][j], Form("%s_%s%i_%s", SampleName().c_str(), "tcMet_xy_outz_sjp", j, dilepton_hypo_names[i]),
			  Form("%s%i", ";tcMET_x+y_outz_sjp", j), tcmcbin, -tcmcmax, tcmcmax);;
	  NewHist(htcmetinzxy_sjp[i][j], Form("%s_%s%i_%s", SampleName().c_str(), "tcMet_xy_inz_sjp", j, dilepton_hypo_names[i]),
			  Form("%s%i", ";tcMET_x+y_inz_sjp", j), tcmcbin, -tcmcmax, tcmcmax);  ;
	}

	NewHist(hsumjetpt[i],    Form("%s_%s_%s", SampleName().c_str(), "sumJetPt", dilepton_hypo_names[i]), ";sum jet pt", sjpbin, 0, sjpmax);
	NewHist(hsumjetptinZ[i], Form("%s_%s_%s", SampleName().c_str(), "sumJetPt_inz", dilepton_hypo_names[i]), ";sum jet pt in Z", sjpbin, 0, sjpmax);
	NewHist(hsumjetptouZ[i], Form("%s_%s_%s", SampleName().c_str(), "sumJetPt_outz", dilepton_hypo_names[i]), ";sum jet pt out Z", sjpbin, 0, sjpmax);

	htcmetxvy[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "tcMet_xvy", dilepton_hypo_names[i]), ";tcMET_x;tcMET_y", tcmcbin, -tcmcmax, tcmcmax, tcmcbin, -tcmcmax, tcmcmax);

	helPt_[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), "elPt", dilepton_hypo_names[i]), ";el pt", 100, 0, 100);
	hmuPt_[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), "muPt", dilepton_hypo_names[i]), ";mu pt", 100, 0, 100);

	// call Sumw2 on all histograms
	helPt_[i]->Sumw2();
	hmuPt_[i]->Sumw2();
	helPt_[i]->SetFillColor(sample_.histo_color);
	helPt_[i]->SetLineColor(sample_.histo_color);
	hmuPt_[i]->SetFillColor(sample_.histo_color);
	hmuPt_[i]->SetLineColor(sample_.histo_color);

  }

  // or use the N - 1 technology (see NMinus1Hist.h)
  // arguments are as follows: sample, name, binning, required cuts, cuts that are relaxed for the N - 1 plot
  // for the lt N - 1 plot, we relax the lt pt requirement
  //hltPt_		= new NMinus1Hist(sample_, "ltPt"   ,	 150, 0, 150, cuts_, CUT_BIT(CUT_LT_PT));
  // same for ll
  //hllPt_		= new NMinus1Hist(sample_, "llPt"   ,	 150, 0, 150, cuts_, CUT_BIT(CUT_LL_PT));
  // for the dilepton mass plot, we relax any cut to do with the Z 
  //hdilMass_		= new NMinus1Hist(sample_, "dilMass",	 100, 0, 300, cuts_, CUT_BIT(CUT_PASS_ZVETO) | CUT_BIT(CUT_PASS_ADDZVETO) | CUT_BIT(CUT_IN_Z_WINDOW));
  //cout << "kfactor: " << sample_.kFactor << "  weight: " << Weight(0) << endl;
  //cout << "  weight: " << Weight(1) << endl;

}


bool Looper::FilterEvent() {  //if  (FilterEvent()) { continue;
  //printevt = false;
  if( ( cms2.evt_event() == 43 && cms2.evt_lumiBlock() == 1023 )
	  || ( cms2.evt_event() == 45 && cms2.evt_lumiBlock() == 101185 )
	  || ( cms2.evt_event() == 152 && cms2.evt_lumiBlock() == 101793 ) )
	{
	//printevt = true;
  }

  return false; 
}

cuts_t Looper::EventSelect() {
  //this function is not called(?)
  cuts_t ret = 0;
  return ret;
}

cuts_t Looper::DilepSelect( int i_hyp ) {
  cuts_t ret = 0;

  //double mass = (cms2.hyp_lt_p4()[i_hyp] + cms2.hyp_ll_p4()[i_hyp]).mass();
  double mass = cms2.hyp_p4()[i_hyp].mass();
  if( mass < 106 && mass > 76 )
	ret |= CUT_BIT(CUT_IN_Z_WINDOW);

  if( (mass > 61 && mass < 76) || (mass > 106 && mass < 121) ) //15 above and below normal window
  	ret |= CUT_BIT(CUT_IN_OUTER_Z);

  //remove this one
  //this cut passes if opp flav or same flav and outside z window
  if( ((mass > 106 || mass < 76) && abs(cms2.hyp_lt_id()[i_hyp]) == abs(cms2.hyp_ll_id()[i_hyp]))
	  || abs(cms2.hyp_lt_id()[i_hyp]) != abs(cms2.hyp_ll_id()[i_hyp]) )
	ret |= CUT_BIT(CUT_Z_SF);

  //keep below
  if( cms2.hyp_lt_p4()[i_hyp].pt() > 20 && cms2.hyp_ll_p4()[i_hyp].pt() > 10 )
	ret |= CUT_BIT(CUT_HYP_PT);
  else
	cout << "hyps don't make 20/10\n";

  //our hypdilepmaker doesn't sort by pt if both are greater than the tight cut
  if( cms2.hyp_lt_p4()[i_hyp].pt() < cms2.hyp_ll_p4()[i_hyp].pt() && cms2.hyp_ll_p4()[i_hyp].pt() < 20)
	cout << "ll = " << cms2.hyp_ll_p4()[i_hyp].pt() << " > lt = " << cms2.hyp_lt_p4()[i_hyp].pt() << "\n"; //stupid check

  if( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 )
	ret |= CUT_BIT(CUT_OPP_SIGN);

  if( GoodSusyTrigger( cms2.hyp_type()[i_hyp] ) )
	ret |= CUT_BIT(CUT_SUSY_TRIGGER);
  //if (! GoodSusyTrigger( cms2.hyp_type()[i_hyp] ) ) return; //from fkw
  
  if( GoodSusyLeptonID( cms2.hyp_lt_id()[i_hyp], cms2.hyp_lt_index()[i_hyp]) )
	ret |= CUT_BIT(CUT_LT_GOOD);

  if( GoodSusyLeptonID( cms2.hyp_ll_id()[i_hyp], cms2.hyp_ll_index()[i_hyp]) )
	ret |= CUT_BIT(CUT_LL_GOOD);

  if( PassSusyLeptonIsolation( cms2.hyp_lt_id()[i_hyp], cms2.hyp_lt_index()[i_hyp]) )
	ret |= CUT_BIT(CUT_LT_ISO);

  if( PassSusyLeptonIsolation( cms2.hyp_ll_id()[i_hyp], cms2.hyp_ll_index()[i_hyp]) )
	ret |= CUT_BIT(CUT_LL_ISO);

  //remove
  if( passMetVJets09(100.0, true) ) //use only for reference--100gev
	ret |= CUT_BIT(CUT_TCMET);

  sumjetpt = 0; //reset from last hyp
  sjpbin = 0;
  vector<LorentzVector> calojets = getCaloJets( i_hyp );
	
  for( unsigned int i=0; i<calojets.size(); i++ )
	sumjetpt += calojets[i].pt(); //made this a private member so can use in filldilephistos

  //get right sjpbin--default is zero bin = 0 sumjetpt (all jets below threshold)
  //right now nsjpbins is 4
  if( sumjetpt >= 30. && sumjetpt < 100. ) //30.0 will be kept
	sjpbin = 1;
  else if( sumjetpt >= 100. && sumjetpt < 200. )
	sjpbin = 2;
  else if( sumjetpt >= 200. )
	sjpbin = 3;

  //remove both below
  if( sumjetpt > 200 )
	ret |= CUT_BIT(CUT_SUMJETPT200);

  if( calojets.size() >= 2 )
	ret |= CUT_BIT(CUT_NJETS2);

  //all of these should yield 616.2 events for ttbar with same flavor zveto

  return ret;
}

cuts_t Looper::TrilepSelect (int i_hyp) {
     cuts_t ret = 0;
     return ret;
}

cuts_t Looper::QuadlepSelect (int i_hyp) {
     cuts_t ret = 0;
     return ret;
}

void Looper::FillEventHistos () {
}

//THIS IS CALLED FOR EVERY HYP IN LOOPERBASE
void Looper::FillDilepHistos (int i_hyp) {
  // everybody histogram needs to know what hypothesis he is 
  const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
  // and what the event weight is 
  const double weight = Weight(i_hyp);

  // these are the cuts that the candidate passes:
  cuts_t cuts_passed = DilepSelect(i_hyp);

  double tcmetx = cms2.evt_tcmet()*cos( cms2.evt_tcmetPhi() );
  double tcmety = cms2.evt_tcmet()*sin( cms2.evt_tcmetPhi() );

  if( (cuts_passed & inz_metres) == inz_metres ) {
	//cout << "tcmet " << cms2.evt_tcmet() << endl;
	hsumjetptinZ[DILEPTON_ALL]->Fill( sumjetpt, weight );
	hsumjetptinZ[myType]->Fill( sumjetpt, weight );

	htcmetinz[DILEPTON_ALL]->Fill( cms2.evt_tcmet(), weight );
	htcmetinz[myType]->Fill( cms2.evt_tcmet(), weight );

	htcmetinzx[DILEPTON_ALL]->Fill( tcmetx, weight );
	htcmetinzx[myType]->Fill( tcmetx, weight );

	htcmetinzy[DILEPTON_ALL]->Fill( tcmety, weight );
	htcmetinzy[myType]->Fill( tcmety, weight );

	//the xy histogram is for each x and y components, so fill twice
	htcmetinzxy[DILEPTON_ALL]->Fill( tcmetx, weight );
	htcmetinzxy[myType]->Fill( tcmetx, weight );
	htcmetinzxy[DILEPTON_ALL]->Fill( tcmety, weight );
	htcmetinzxy[myType]->Fill( tcmety, weight );

	htcmetinzxy_sjp[DILEPTON_ALL][sjpbin]->Fill( tcmetx, weight );
	htcmetinzxy_sjp[myType][sjpbin]->Fill( tcmetx, weight );
	htcmetinzxy_sjp[DILEPTON_ALL][sjpbin]->Fill( tcmety, weight );
	htcmetinzxy_sjp[myType][sjpbin]->Fill( tcmety, weight );
  }
  else if( (cuts_passed & outerz_metres) == outerz_metres ) {
	hsumjetptouZ[DILEPTON_ALL]->Fill( sumjetpt, weight );
	hsumjetptouZ[myType]->Fill( sumjetpt, weight );

	htcmetouz[DILEPTON_ALL]->Fill( cms2.evt_tcmet(), weight );
	htcmetouz[myType]->Fill( cms2.evt_tcmet(), weight );

	htcmetouzx[DILEPTON_ALL]->Fill( tcmetx, weight );
	htcmetouzx[myType]->Fill( tcmetx, weight );

	htcmetouzy[DILEPTON_ALL]->Fill( tcmety, weight );
	htcmetouzy[myType]->Fill( tcmety, weight );

	//the xy histogram is for each x and y components, so fill twice
	htcmetouzxy[DILEPTON_ALL]->Fill( tcmetx, weight );
	htcmetouzxy[myType]->Fill( tcmetx, weight );
	htcmetouzxy[DILEPTON_ALL]->Fill( tcmety, weight );
	htcmetouzxy[myType]->Fill( tcmety, weight );

	htcmetouzxy_sjp[DILEPTON_ALL][sjpbin]->Fill( tcmetx, weight );
	htcmetouzxy_sjp[myType][sjpbin]->Fill( tcmetx, weight );
	htcmetouzxy_sjp[DILEPTON_ALL][sjpbin]->Fill( tcmety, weight );
	htcmetouzxy_sjp[myType][sjpbin]->Fill( tcmety, weight );
  }
  
  // for TH1/TH2, we have to check explicitly whether the candidate passes
  if( (cuts_passed & cuts_) == cuts_ ) {
	hsumjetpt[DILEPTON_ALL]->Fill( sumjetpt, weight );
	hsumjetpt[myType]->Fill( sumjetpt, weight );

	htcmet[DILEPTON_ALL]->Fill( cms2.evt_tcmet(), weight );
	htcmet[myType]->Fill( cms2.evt_tcmet(), weight );
	
	htcmetx[DILEPTON_ALL]->Fill( tcmetx, weight );
	htcmetx[myType]->Fill( tcmetx, weight );
	
	htcmety[DILEPTON_ALL]->Fill( tcmety, weight );
	htcmety[myType]->Fill( tcmety, weight );
	
	//the xy histogram is for each x and y components, so fill twice
	htcmetxy[DILEPTON_ALL]->Fill( tcmetx, weight );
	htcmetxy[myType]->Fill( tcmetx, weight );
	htcmetxy[DILEPTON_ALL]->Fill( tcmety, weight );
	htcmetxy[myType]->Fill( tcmety, weight );

	htcmetxy_sjp[DILEPTON_ALL][sjpbin]->Fill( tcmetx, weight );
	htcmetxy_sjp[myType][sjpbin]->Fill( tcmetx, weight );
	htcmetxy_sjp[DILEPTON_ALL][sjpbin]->Fill( tcmety, weight );
	htcmetxy_sjp[myType][sjpbin]->Fill( tcmety, weight );
	
	htcmetxvy[DILEPTON_ALL]->Fill( tcmetx, tcmety, weight ); //2d
	htcmetxvy[myType]->Fill( tcmetx, tcmety, weight ); //2d
	
	//if( sumjetpt > 0 ) {
	
	// and then fill defaults
	if( abs(cms2.hyp_lt_id()[i_hyp]) == 11 ) {
	  helPt_[DILEPTON_ALL]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
	  helPt_[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
	} else {
	  hmuPt_[DILEPTON_ALL]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
	  hmuPt_[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
	}
	if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  helPt_[DILEPTON_ALL]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
	  helPt_[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
	} else {
	  hmuPt_[DILEPTON_ALL]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
	  hmuPt_[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
	}
  }

  //if( (cuts_passed & cuts_) == cuts_ ) {
  if( printevt ) {
	cout << "run " << cms2.evt_run() << " evt " << cms2.evt_event() << " lumi " << cms2.evt_lumiBlock()
		 <<  " hyp " << i_hyp + 1 << " of " << cms2.hyp_lt_id().size() << endl; //+1 is just bc 0 is first
	cout << "cut bit: " << cuts_passed << endl;
	//dumpDocLines();
  }

  if( (cuts_passed & cuts_) == cuts_ ) {
	counttt += weight;
	//cout << "run " << cms2.evt_run() << " evt " << cms2.evt_event() << " lumi " << cms2.evt_lumiBlock() << endl;
	//dumpDocLines();
	// if the candidate passed, we count it
	cands_passing_[myType] += weight;
	cands_passing_w2_[myType] += weight * weight;
	cands_count_[myType]++;
	cands_passing_[DILEPTON_ALL] += weight;
	cands_passing_w2_[DILEPTON_ALL] += weight * weight;
	cands_count_[DILEPTON_ALL]++;
  }
  // for the NMinus1Hist, the histogram checks the cuts for us
  //hltPt_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
  //hllPt_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
  //hdilMass_->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].mass(), weight);
}

void Looper::FillTrilepHistos (int i_hyp) {
}

void Looper::FillQuadlepHistos (int i_hyp) {
}

void Looper::End () {

  cout << "kfactor: " << sample_.kFactor << "  weight: " << Weight(0) << endl;
  cout << "hyps passed: " << counttt << endl << endl;
  //cout << "  weight: " << Weight(1) << endl;  
  //cout << endl << baseline_susy << endl << cuts_ << endl << endl;

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
