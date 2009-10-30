// -*- C++ -*-

#include <unistd.h>
#include <string>
#include "ABCDLooper.h"
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
1 << LOOP_QCD		| 
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
void printTable (const ABCDLooper **hists, int n, const char *fname, 
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
	fprintf(f, "|*%14s*", hists[j]->SampleName().c_str());
  }
  fprintf(f, "|%12s  |\n", "total");
  //dilep
  //for (int i = 0; i < 4; ++i) {
  for (int i = 3; i >= 0; --i) {
	fprintf(f, "|%6s ", dilepton_hypo_names[i]);
	double cands = 0;
	double w2 = 0;
	for (int j = 0; j < n; ++j) {
	  fprintf(f, "|  %6.1f +- %6.1f", 
			  hists[j]->DCandsPassing(DileptonHypType(i)),
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

  if (f != stdin) 
	fclose(f);

}

template <class Looper> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff)
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

  Looper looper_qcd30(fQCDpt30()			, cuts, log.c_str()); if (which_ones & (1 << LOOP_QCD)) 	looper_qcd30.Loop();
  Looper looper_qcd80(fQCDpt80()			, cuts, log.c_str()); if (which_ones & (1 << LOOP_QCD)) 	looper_qcd80.Loop();
  Looper looper_wejet_alp(fWejetsAlpgenSingle()		, cuts, log.c_str()); if (which_ones & (1 << LOOP_WJET_ALP)) 	looper_wejet_alp.Loop();
  Looper looper_wmjet_alp(fWmjetsAlpgenSingle()		, cuts, log.c_str()); if (which_ones & (1 << LOOP_WJET_ALP)) 	looper_wmjet_alp.Loop();
  Looper looper_wtjet_alp(fWtjetsAlpgenSingle()		, cuts, log.c_str()); if (which_ones & (1 << LOOP_WJET_ALP)) 	looper_wtjet_alp.Loop();
  Looper looper_qcd_mu15(fInclusiveMuPt15Single()	, cuts, log.c_str()); if (which_ones & (1 << LOOP_QCD_MU15)) 	looper_qcd_mu15.Loop();
  Looper looper_qcd_em(qcd_em				        , cuts, log.c_str()); if (which_ones & (1 << LOOP_QCD_EM)) 	    looper_qcd_em.Loop();
  Looper looper_qcd_photon(qcd_photon 			    , cuts, log.c_str()); if (which_ones & (1 << LOOP_QCD_PHOTON)) 	looper_qcd_photon.Loop();
  Looper looper_qcd_bctoe(qcd_bctoe                 , cuts, log.c_str()); if (which_ones & (1 << LOOP_QCD_BCTOE))   looper_qcd_bctoe.Loop();
  Looper looper_zeejet_alp(fZeejetsAlpgenSingle()	, cuts, log.c_str()); if (which_ones & (1 << LOOP_ZJET_ALP)) 	looper_zeejet_alp.Loop();
  Looper looper_zmmjet_alp(fZmmjetsAlpgenSingle()	, cuts, log.c_str()); if (which_ones & (1 << LOOP_ZJET_ALP)) 	looper_zmmjet_alp.Loop();
  Looper looper_zttjet_alp(fZttjetsAlpgenSingle()	, cuts, log.c_str()); if (which_ones & (1 << LOOP_ZJET_ALP)) 	looper_zttjet_alp.Loop();
  Looper looper_ttbar(fttbarSingle()			    , cuts, log.c_str()); if (which_ones & (1 << LOOP_TTBAR))       looper_ttbar.Loop();

  // when all the loopers are done, we save the histograms to file
  saveHist(hist.c_str());

  // then we collect them all and print a table
  const Looper *loopers[] = { 
	&looper_qcd30,
	&looper_qcd80,
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

  printTable(loopers, sizeof(loopers) / sizeof(Looper *), tbl.c_str(), which_ones);
  return 0;
}

// default yield table
int Results ()
{
  //return run<Looper>(0, "Results", default_samples );
  //return run<Looper>(0, "Results_sig", sig_samples );
  //return run<Looper>(0, "Results", samples );

}

int ABCDResults ()
{
  return run<ABCDLooper>(0, "ABCDResults", default_samples );
  //return run<ABCDLooper>(0, "ABCDResults_sig", sig_samples );
  //return run<ABCDLooper>(0, "ABCDResults", samples );

}

