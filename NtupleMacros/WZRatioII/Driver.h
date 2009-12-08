// -*- C++ -*-

#include <unistd.h>
#include <string>
#include "ABCDLooper.h"
#include "RLooper.h"
#include "Tools/Sample.h"
//#include "LocalSample.h"

#include "Tools/tools.h"

//#include "Cuts.h"

using std::string;

enum {
  // 2_1_X
  LOOP_QCD,
  LOOP_QCD_BCTOE,
  LOOP_QCD_MU15,
  LOOP_QCD_EM,
  LOOP_QCD_PHOTON,
  LOOP_WJET_ALP,
  LOOP_ZJET_ALP,
  LOOP_TTBAR,
};

const uint32 default_samples =
//1 << LOOP_QCD		| 
1 << LOOP_QCD_MU15		| 
1 << LOOP_QCD_EM		| 
1 << LOOP_QCD_PHOTON	| 
1 << LOOP_QCD_BCTOE     |
1 << LOOP_WJET_ALP		| 
1 << LOOP_ZJET_ALP 	|
1 << LOOP_TTBAR;

const uint32 sig_samples =
1 << LOOP_WJET_ALP		| 
1 << LOOP_ZJET_ALP ;

// helper function used to print yield tables
template <class Looper> void printTable (const Looper **hists, int n, const char *fname, 
				 uint32 which_ones)
{
  
  FILE *f = 0;
  if (fname == 0 || strlen(fname) == 0)
	f = stdin;
  else f = fopen(fname, "w");
  if (f == 0) {
	perror("printing table");
	return;
  }
  fprintf(f, "| %6s", "");
  for (int j = 0; j < n; ++j) {
	fprintf(f, "|*%16s*", hists[j]->SampleName().c_str());
  }
  fprintf(f, "|%12s  |\n", "total");
  //dilep     : static const char dilepton_hypo_names[][128] = { "all", "mm", "em", "ee" };
  //enum DileptonHypType {     DILEPTON_ALL,     DILEPTON_MUMU,     DILEPTON_EMU,     DILEPTON_EE  };

  //for (int i = 0; i < 4; ++i) {
  for (int i = 3; i >= 0; --i) {
	fprintf(f, "|%6s ", dilepton_hypo_names[i]);
	double cands = 0;
	double w2 = 0;
	for (int j = 0; j < n; ++j) {
	  fprintf(f, "|  %6.1f +- %6.1f", 
			  hists[j]->DCandsPassing(DileptonHypType(i)), //this is a type case of an int "i" to a DileptonHypType
			  hists[j]->DRMS(DileptonHypType(i)));
 	  cands += hists[j]->DCandsPassing(DileptonHypType(i));
	  w2 += hists[j]->DRMS(DileptonHypType(i)) * hists[j]->DRMS(DileptonHypType(i));
	}
	fprintf(f, "|  %6.1f +- %6.1f|\n", cands, sqrt(w2));
  }
  //single lep
  char* slepnames[3] = {"e", "m", "all"};
  for (int i = 0; i < 3; ++i) {
	fprintf(f, "|%6s ", slepnames[i]);
	double cands = 0;
	double w2 = 0;
	for (int j = 0; j < n; ++j) {
	  fprintf(f, "|  %6.1f +- %6.1f", 
			  hists[j]->SCandsPassing(DileptonHypType(i)),
			  hists[j]->SRMS(DileptonHypType(i)));
	  cands += hists[j]->SCandsPassing(DileptonHypType(i));
	  w2 += hists[j]->SRMS(DileptonHypType(i)) * hists[j]->SRMS(DileptonHypType(i));
	}
	fprintf(f, "|  %6.1f +- %6.1f|\n", cands, sqrt(w2));
  }

  //0=e, 1=m, 2=all
  fprintf(f, "\n\nElectrons\n\n");
  printR( hists, n, f, 0 );
  fprintf(f, "\n\nMuons\n\n");
  printR( hists, n, f, 1 );

  if (f != stdin) 
	fclose(f);
}


template <class Looper> void printR( const Looper **loopers, int n, FILE* f, int flv = -1 ) {
  double Nw = 0, Nz = 0; //yields
  //double Nwrms = 0, Nzrms = 0;
  double Fw = 0, Fz = 0; //fake yield
  //double Fwrms = 0, Fzrms = 0;
  double et = 0, siget = 0.; //tight efficiencies
  double el = 0, sigel = 0.; //loose efficiencies
  double enu = 0.99, sigenu = 0.01; //neutrino efficiencies
  double accz = 0.4, sigaccz = 0.02; //acceptance z
  double accw = 0.5, sigaccw = 0.01; //acceptance w

  DileptonHypType sidx, didx;
  if( flv == 0 ) { //electrons
	sidx = DileptonHypType(0);
	didx = DILEPTON_EE;
	et = 0.85;
	siget = 0.02;
	el = 0.90;
	sigel = 0.02;
  }
  else if( flv == 1 ) { //muons
	sidx = DileptonHypType(1);
	didx = DILEPTON_MUMU;
	et = 0.87;
	siget = 0.02;
	el = 0.92;
	sigel = 0.02;
  }
  else {
	cout << "printR: only e, mu individually supported now\n";
	return;
  }

  //get the quantities out of the loopers
  for( int i=0;i < n; i++ ) {
	if( loopers[i]->isSsignal() ) {
	  Nw += loopers[i]->SCandsPassing(sidx); //the all
	  //Nwrms += loopers[i]->SRMS(sidx)*loopers[i]->SRMS(sidx);
	  Fz += loopers[i]->DCandsPassing(didx); //real w->fake z
	  //Fzrms += loopers[i]->DRMS(didx)*loopers[i]->DRMS(didx);
	}
	else if( loopers[i]->isDsignal() ) {
	  Nz += loopers[i]->DCandsPassing(didx); //the all
	  //Nzrms += loopers[i]->DRMS(didx)*loopers[i]->DRMS(didx);
	  Fw += loopers[i]->SCandsPassing(sidx); //real z->fake w
	  //Fwrms += loopers[i]->SRMS(sidx)*loopers[i]->SRMS(sidx);
	}
	else { //bkg for both s, d
	  Fz += loopers[i]->DCandsPassing(didx); //fake z
	  //Fzrms += loopers[i]->DRMS(didx)*loopers[i]->DRMS(didx);
	  Fw += loopers[i]->SCandsPassing(sidx); //fake w
	  //Fwrms += loopers[i]->SRMS(sidx)*loopers[i]->SRMS(sidx);
	}
  }

  //calculate R
  double R = Nw*accz*et/(Nz*accw*enu);

  //calculate sigma^2(R)
  double signwnz2 = (1/(Nz*Nz))*(Nw + 2*Fw + Fw*Fw/4) + ((Nw*Nw)/(Nz*Nz*Nz*Nz))*(Nz + 2*Fz + Fz*Fz/4);
  // sigma^2(ez/ew)
  double sigezew2 = siget*siget*(1-2*el)*(1-2*el)/(enu*enu) + sigel*sigel*4*(1-et)*(1-et)/(enu*enu) + sigenu*sigenu*(et+2*el*(1-et))*(et+2*el*(1-et))/(enu*enu*enu*enu);
  //double sigezew2 = siget*siget/(enu*enu) + sigenu*sigenu*et*et/(enu*enu*enu*enu); // tight-tight formula
  // sigma^2(accz/accw)
  double sigazaw2 = sigaccz*sigaccz/(accw*accw) + sigaccw*sigaccw*accz*accz/(accw*accw*accw*accw);

  double t1 = signwnz2*accz*accz*et*et/(accw*accw*enu*enu);
  double t2 = sigazaw2*Nw*Nw*et*et/(Nz*Nz*enu*enu);
  double t3 = sigezew2*Nw*Nw*accz*accz/(Nz*Nz*accw*accw);
  double sigR = sqrt( t1 + t2 + t3 );

  fprintf(f, "R = %f +- %f\n", R, sigR);
  //fprintf(f, "\nsignwnz = %f, sigezew = %f, sigazaw = %f\n", sqrt(signwnz2), sqrt(sigezew2), sqrt(sigazaw2));
  fprintf(f, "\nt1 = %f, t2 = %f, t3 = %f\n", t1, t2, t3);
  cout << "R = " << R << " +- " << sigR << endl;
}


template <class Looper> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff, cuts_t cuts2=0)
{
  const string hist = name + ".root";
  const string tbl = name + ".tbl";
  const string log = name + ".log";

  Sample qcd_em = fQCDEMenrichedPt20to30() + fQCDEMenrichedPt30to80() + fQCDEMenrichedPt80to170();
  qcd_em.name = "QCDEMEnriched"; 
  Sample qcd_photon = fPhotonJetPt15to20() + fPhotonJetPt20to25() + fPhotonJetPt25to30() + 
	fPhotonJetPt30to35() + fPhotonJetPt35 ();
  qcd_photon.name = "PhotonJet";
  Sample qcd_bctoe = fQCDBCtoEPt20to30() + fQCDBCtoEPt30to80() + fQCDBCtoEPt80to170();
  qcd_bctoe.name = "QCDBCtoE";

  Looper looper_qcd30(fQCDpt30()					, cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_QCD)) 	looper_qcd30.Loop();
  Looper looper_qcd80(fQCDpt80()					, cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_QCD)) 	looper_qcd80.Loop();
  Looper looper_wejet_alp(fWejetsAlpgenSingle()		, cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_WJET_ALP)) 	looper_wejet_alp.Loop();
  Looper looper_wmjet_alp(fWmjetsAlpgenSingle()		, cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_WJET_ALP)) 	looper_wmjet_alp.Loop();
  Looper looper_wtjet_alp(fWtjetsAlpgenSingle()		, cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_WJET_ALP)) 	looper_wtjet_alp.Loop();
  Looper looper_qcd_mu15(fInclusiveMuPt15Single()	, cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_QCD_MU15)) 	looper_qcd_mu15.Loop();
  Looper looper_qcd_em(qcd_em				        , cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_QCD_EM)) 	    looper_qcd_em.Loop();
  Looper looper_qcd_photon(qcd_photon 			    , cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_QCD_PHOTON)) 	looper_qcd_photon.Loop();
  Looper looper_qcd_bctoe(qcd_bctoe                 , cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_QCD_BCTOE))   looper_qcd_bctoe.Loop();
  Looper looper_zeejet_alp(fZeejetsAlpgenSingle()	, cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_ZJET_ALP)) 	looper_zeejet_alp.Loop();
  Looper looper_zmmjet_alp(fZmmjetsAlpgenSingle()	, cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_ZJET_ALP)) 	looper_zmmjet_alp.Loop();
  Looper looper_zttjet_alp(fZttjetsAlpgenSingle()	, cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_ZJET_ALP)) 	looper_zttjet_alp.Loop();
  Looper looper_ttbar(fttbarSingle()			    , cuts, log.c_str(), cuts2); if (which_ones & (1 << LOOP_TTBAR))       looper_ttbar.Loop();

  // when all the loopers are done, we save the histograms to file
  saveHist(hist.c_str());

  // then we collect them all and print a table
  const Looper *loopers[] = { 
	//&looper_qcd30,
	//&looper_qcd80,
	&looper_qcd_mu15,
	&looper_qcd_em,
	&looper_qcd_photon,
	&looper_qcd_bctoe,
	&looper_ttbar,
	&looper_wejet_alp,
	&looper_wmjet_alp,
	&looper_wtjet_alp,
	&looper_zeejet_alp,
	&looper_zmmjet_alp,
	&looper_zttjet_alp,
  };

  looper_wejet_alp.setSsignal( true );
  looper_wmjet_alp.setSsignal( true );
  looper_wtjet_alp.setSsignal( true ); //dont use these b'c then need tau efficiency???
  looper_zeejet_alp.setDsignal( true );
  looper_zmmjet_alp.setDsignal( true );
  looper_zttjet_alp.setDsignal( true );

  cout << name << endl;
  printTable(loopers, sizeof(loopers) / sizeof(Looper *), tbl.c_str(), which_ones);
  return 0;
}

// default yield table
int Results(cuts_t cut1 = 0, cuts_t cut2 = 0)
{
  //return run<Looper>(0, "Results", default_samples );
  //return run<Looper>(0, "Results_sig", sig_samples );
  //return run<Looper>(0, "Results", samples );

}

int ABCDResults(cuts_t cut1 = 0, cuts_t cut2 = 0)
{
  return run<ABCDLooper>(cut1, "ABCDResults", default_samples, cut2 );
  //return run<ABCDLooper>(0, "ABCDResults", (1<<LOOP_WJET_ALP) );
  //return run<ABCDLooper>(0, "ABCDResults_sig", sig_samples );
  //return run<ABCDLooper>(0, "ABCDResults", samples );

}

int RResults(cuts_t cut1 = 0, cuts_t cut2 = 0, string name="RResults")
{
  return run<RLooper>(cut1, name, default_samples, cut2 );
  //return run<RLooper>(cut1, name, (1<<LOOP_ZJET_ALP), cut2 );
}
