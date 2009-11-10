#include <unistd.h>
#include <string>
#include "Looper.h"
#include "Tools/Sample.h"
#include "LocalSample.h"

#include "Tools/tools.h"

#include "Cuts.h"

using std::string;

enum {

	// V02-00-08
	LOOP_WENU_7TeV,
	LOOP_QCD30_7TeV,
	LOOP_PHOTONJET_7TeV,

	// V02-00-12
     LOOP_WW    ,
     LOOP_WZ    ,
     LOOP_ZZ    ,
     LOOP_WE ,
     LOOP_WM ,
     LOOP_WT ,
     LOOP_DYEE  ,
     LOOP_DYMM  ,
     LOOP_DYTT  ,
     LOOP_TTBAR ,
     LOOP_QCD30,
     LOOP_MU30,
     LOOP_PHOTONJET20_170.
     LOOP_PHOTONJET20_30,
     LOOP_PHOTONJET30_50,
     LOOP_PHOTONJET50_80,
     LOOP_PHOTONJET80_120,
     LOOP_PHOTONJET120_170,
     LOOP_EM20_170,
     LOOP_EM20_30,
     LOOP_EM30_80,
     LOOP_EM80_170,
     LOOP_BCE20_170,
     LOOP_BC20_30,
     LOOP_BC30_80,
     LOOP_BC80_170,
     LOOP_BC20_170,

};

// helper function used to print yield tables
void printTable (const Looper **hists, int n, const char *fname, 
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
     fprintf(f, "| %10s", "");
     for (int j = 0; j < n; ++j) {
	  fprintf(f, "|  *%30s*  ", hists[j]->SampleName().c_str());
     }

/*
     fprintf(f, "|%30s  |\n", "total");
     for (int i = 0; i < 4; ++i) {
	  fprintf(f, "|%10s  ", dilepton_hypo_names[i]);
	  double cands = 0;
	  double w2 = 0;
	  for (int j = 0; j < n; ++j) {
	       fprintf(f, "|  %10.1f &plusmn; %10.1f", 
		       hists[j]->CandsPassing(DileptonHypType(i)),
		       hists[j]->RMS(DileptonHypType(i)));
	       cands += hists[j]->CandsPassing(DileptonHypType(i));
	       w2 += hists[j]->RMS(DileptonHypType(i)) * 
		    hists[j]->RMS(DileptonHypType(i));
	  }
	  fprintf(f, "|  %10.1f &plusmn; %10.1f|\n", cands, sqrt(w2));
     }
*/

     const char detectorNames[][128] = {"EB", "EE", "ALL"};

     fprintf(f, "WEfficiency Results");
     fprintf(f, "|%30s  |\n", "total");
     for (int i = 0; i < 3; ++i) {
          fprintf(f, "|%10s  ", detectorNames[i]);
          double cands = 0;
          double w2 = 0;
          for (int j = 0; j < n; ++j) {
               fprintf(f, "|  %10.1f &plusmn; %10.1f", 
                       hists[j]->CandsPassingW(i),
                       hists[j]->RMSW(i));
               cands += hists[j]->CandsPassingW(i);
               w2 += hists[j]->RMSW(i) * 
                    hists[j]->RMSW(i);
          }
          fprintf(f, "|  %10.1f &plusmn; %10.1f|\n", cands, sqrt(w2));
     }


     fprintf(f, "AN2009-98 Results");
     fprintf(f, "|%30s  |\n", "total");
     for (int i = 0; i < 3; ++i) {
          fprintf(f, "|%10s  ", detectorNames[i]);
          double cands = 0;
          double w2 = 0;
          for (int j = 0; j < n; ++j) {
               fprintf(f, "|  %10.1f &plusmn; %10.1f",
                       hists[j]->CandsPassingAN2009_98(i),
                       hists[j]->RMSAN2009_98(i));
               cands += hists[j]->CandsPassingAN2009_98(i);
               w2 += hists[j]->RMSAN2009_98(i) *
                    hists[j]->RMSAN2009_98(i);
          }
          fprintf(f, "|  %10.1f &plusmn; %10.1f|\n", cands, sqrt(w2));
     }
	


     if (f != stdin) 
	  fclose(f);
}

template <class Looper> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff)
{
     const string hist = name + ".root";
     const string tbl = name + ".tbl";
     const string log = name + ".log";

	// V02-00-08
     Looper looper_wenu_7TeV(fWenu_7TeV(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_WENU_7TeV)) looper_wenu_7TeV.Loop();
     Looper looper_qcd30_7TeV(fQCD_Pt30_7TeV(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_QCD30_7TeV)) looper_qcd30_7TeV.Loop();
     Looper looper_photonjet_7TeV(fPhotonJet_7TeV(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_PHOTONJET_7TeV)) looper_photonjet_7TeV.Loop();

        // V02-00-12
     Looper looper_ww(fWW(), cuts, log.c_str()); 
	if (which_ones & (1 << LOOP_WW)) looper_ww.Loop();
     Looper looper_wz(fWZ(), cuts, log.c_str()); 
	if (which_ones & (1 << LOOP_WZ)) looper_wz.Loop();
     Looper looper_zz(fZZ(), cuts, log.c_str()); 
	if (which_ones & (1 << LOOP_ZZ)) looper_zz.Loop();
     Looper looper_ttbar(fttbar(), cuts, log.c_str()); 
	if (which_ones & (1 << LOOP_TTBAR)) looper_ttbar.Loop();

     Looper looper_bc20_170(fQCDBCtoEPt20to170(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_BC20_170)) looper_bc20_170.Loop();
     Looper looper_bc20_30(fQCDBCtoEPt20to30(), cuts, log.c_str()); 
	if (which_ones & (1 << LOOP_BC20_30)) looper_bc20_30.Loop();
     Looper looper_bc30_80(fQCDBCtoEPt30to80(), cuts, log.c_str()); 
	if (which_ones & (1 << LOOP_BC30_80)) looper_bc30_80.Loop();
     Looper looper_bc80_170(fQCDBCtoEPt80to170(), cuts, log.c_str()); 
	if (which_ones & (1 << LOOP_BC80_170)) looper_bc80_170.Loop();

     Looper looper_em20_170(fQCDEMenrichedPt20to170(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_EM20_170)) looper_em20_170.Loop();
     Looper looper_em20_30(fQCDEMenrichedPt20to30(), cuts, log.c_str()); 
	if (which_ones & (1 << LOOP_EM20_30)) looper_em20_30.Loop();
     Looper looper_em30_80(fQCDEMenrichedPt30to80(), cuts, log.c_str()); 
	if (which_ones & (1 << LOOP_EM30_80)) looper_em30_80.Loop();
     Looper looper_em80_170(fQCDEMenrichedPt80to170(), cuts, log.c_str()); 
	if (which_ones & (1 << LOOP_EM80_170)) looper_em80_170.Loop();

     Looper looper_mu30(fInclusiveMu5Pt30(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_MU30)) looper_mu30.Loop();
     Looper looper_qcd30(fQCDpt30(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_QCD30)) looper_qcd30.Loop();

     Looper looper_we(fWe(), cuts, log.c_str()); 
        if (which_ones & (1 << LOOP_WE)) looper_we.Loop();
     Looper looper_wm(fWm(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_WM)) looper_wm.Loop();
     Looper looper_wt(fWt(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_WT)) looper_wt.Loop();

     Looper looper_dyee(fDYee(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_DYEE)) looper_dyee.Loop();
     Looper looper_dymm(fDYmm(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_DYMM)) looper_dymm.Loop();
     Looper looper_dytt(fDYtt(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_DYTT)) looper_dytt.Loop();

     Looper looper_photonjet20_170(fPhotonJetPt20to170(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_PHOTONJET20_170)) looper_photonjet20_170.Loop();
     Looper looper_photonjet20_30(fPhotonJetPt20to30(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_PHOTONJET20_30)) looper_photonjet20_30.Loop();
     Looper looper_photonjet30_50(fPhotonJetPt30to50(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_PHOTONJET30_50)) looper_photonjet30_50.Loop();
     Looper looper_photonjet50_80(fPhotonJetPt50to80(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_PHOTONJET50_80)) looper_photonjet50_80.Loop();
     Looper looper_photonjet80_120(fPhotonJetPt80to120(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_PHOTONJET80_120)) looper_photonjet80_120.Loop();
     Looper looper_photonjet120_170(fPhotonJetPt120to170(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_PHOTONJET120_170)) looper_photonjet120_170.Loop();

     // when all the loopers are done, we save the histograms to file
     saveHist(hist.c_str());

     // then we collect them all and print a table
     const Looper *loopers[] = { 

	// V02-00-08
	&looper_wenu_7TeV,
	&looper_qcd30_7TeV,
	&looper_photonjet_7TeV,

	// V02-00-12
	&looper_ww,
	&looper_wz,
	&looper_zz,
	&looper_we,
	&looper_wm,
	&looper_wt,
	&looper_dyee,
	&looper_dymm,
	&looper_dytt,
	&looper_ttbar,
	&looper_qcd30,
	&looper_mu30,
        &looper_photonjet20_170,
	&looper_photonjet20_30,
        &looper_photonjet30_50,
        &looper_photonjet50_80,
        &looper_photonjet80_120,
        &looper_photonjet120_170,
        &looper_em20_170,
	&looper_em20_30,
        &looper_em30_80,
        &looper_em80_170,
        &looper_bc20_170,
        &looper_bc20_30,
        &looper_bc30_80,
        &looper_bc80_170,

     };

     printTable(loopers, sizeof(loopers) / sizeof(Looper *), tbl.c_str(), which_ones);

     return 0;
}

// default yield table
int Results_tcmet30 ()
{
     return run<Looper>( (CUT_BIT(ELE_PT_20)) | (CUT_BIT(ELE_ISO_10) | CUT_BIT(EVT_JPT_25) | CUT_BIT(EVT_TCMET_30)),
	"Results_iso10_jpt25_tcmet30", 
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}


//
// Run all
//

uint32 all_samples_w = 
     (1<<LOOP_WE ) |
     (1<<LOOP_QCD30) |
     (1<<LOOP_PHOTONJET20_170) |
     (1<<LOOP_EM20_170) |
     (1<<LOOP_BC20_170) |

uint32 all_samples = (1<<LOOP_WW) |
     (1<<LOOP_WZ    ) |
     (1<<LOOP_ZZ    ) |
     (1<<LOOP_WE ) |
     (1<<LOOP_WM ) |
     (1<<LOOP_WT ) |
     (1<<LOOP_DYEE  ) |
     (1<<LOOP_DYMM  ) |
     (1<<LOOP_DYTT  ) |
     (1<<LOOP_TTBAR ) |
     (1<<LOOP_QCD30) |
     (1<<LOOP_MU30) |
     (1<<LOOP_PHOTONJET20_30) |
     (1<<LOOP_PHOTONJET30_50) |
     (1<<LOOP_PHOTONJET50_80) |
     (1<<LOOP_PHOTONJET80_120) |
     (1<<LOOP_PHOTONJET120_170) |
     (1<<LOOP_EM20_30) |
     (1<<LOOP_EM30_80) |
     (1<<LOOP_EM80_170) |
     (1<<LOOP_BC20_30) |
     (1<<LOOP_BC30_80) |
     (1<<LOOP_BC80_170);

int Results312_pt20_isoV0_tcmet30 ()
{
     return run<Looper>(
        // control
        (CUT_BIT(CONTROL_STUDYW)) |
        // cuts
        (CUT_BIT(ELE_PT_20)) | (CUT_BIT(ELE_ISO_V0) | CUT_BIT(EVT_TCMET_30)),
        "Results_ntupletest", all_samples_w);
}

int Results_ntupletest ()
{
     return run<Looper>(  
	// control
	(CUT_BIT(CONTROL_STUDYW)) |
	// cuts
	(CUT_BIT(ELE_PT_20)) | (CUT_BIT(ELE_ISO_V0) | CUT_BIT(EVT_TCMET_30)),
        "Results_ntupletest", all_samples);
}


//
// electron id (note that here cut bits don't really do anything
// ... probably should add one that is "do electron id study" or something
// ... like that)
//

int Results_isoV0_studies ()
{
     return run<Looper>(  (CUT_BIT(ELE_PT_20)) | (CUT_BIT(ELE_ISO_V0) | CUT_BIT(EVT_TCMET_30)),
        "Results_isoV0_studies", 
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}

//
// studies related to jet veto
//

// iso v1 is ecal, hcal, tkJura
int Results_pt30_isoV1_tcmet30 ()
{    
     return run<Looper>( (CUT_BIT(ELE_PT_30)) | (CUT_BIT(ELE_ISO_V1) | CUT_BIT(EVT_TCMET_30)),
        "Results_pt30_isoV1_tcmet30", 
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}

int Results_pt20_isoV1_tcmet30 ()
{
     return run<Looper>( (CUT_BIT(ELE_PT_20)) | (CUT_BIT(ELE_ISO_V1) | CUT_BIT(EVT_TCMET_30)),
        "Results_pt20_isoV1_tcmet30",
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}

int Results_pt20_isoV1_phimax130_tcmet30 ()
{
     return run<Looper>( (CUT_BIT(ELE_PT_20)) | (CUT_BIT(ELE_ISO_V1) | CUT_BIT(EVT_JPT_PHIMAX_130) | CUT_BIT(EVT_TCMET_30)),
        "Results_pt20_isoV1_phimax130_tcmet30",
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}

int Results_pt20_isoV1_phimax130_tasv1_tcmet30 ()
{    
     return run<Looper>( (CUT_BIT(ELE_PT_20)) | (CUT_BIT(ELE_ISO_V1) | CUT_BIT(EVT_JPT_PHIMAX_130) | CUT_BIT(ELE_TAS_V1) | CUT_BIT(EVT_TCMET_30)),
        "Results_pt20_isoV1_phimax130_tasv1_tcmet30",
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}


//
//
//

int Results_AN2009_098_studies_30v0 ()
{    
     return run<Looper>( (CUT_BIT(ELE_PT_30)) | (CUT_BIT(ELE_ISO_V0) | CUT_BIT(EVT_TCMET_30)),
        "Results_AN2009_098_studies_30v0", 
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}

int Results_AN2009_098_studies_30v1 ()
{ 
     return run<Looper>( (CUT_BIT(ELE_PT_30)) | (CUT_BIT(ELE_ISO_V1) | CUT_BIT(EVT_TCMET_30)),
        "Results_AN2009_098_studies_30v1", 
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}


// note: v2 is combined Ecal-Hcal iso
//
int Results_AN2009_098_studies_30v2 ()
{
     return run<Looper>( (CUT_BIT(ELE_PT_30)) | (CUT_BIT(ELE_ISO_V2) | CUT_BIT(EVT_TCMET_30)),
        "Results_AN2009_098_studies_30v2",
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}

int Results_AN2009_098_studies_20v1 ()
{
     return run<Looper>( (CUT_BIT(ELE_PT_20)) | (CUT_BIT(ELE_ISO_V1) | CUT_BIT(EVT_TCMET_30)),
        "Results_AN2009_098_studies_20v1",
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}

// default yield table
int Results_tcmet30_conv ()
{
     return run<Looper>(  (CUT_BIT(ELE_PT_20)) | (CUT_BIT(ELE_ISO_10) | CUT_BIT(EVT_JPT_25) | CUT_BIT(EVT_TCMET_30) | CUT_BIT(ELE_NOCONV)),
        "Results_iso10_jpt25_tcmet30_conv",
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}

int Results_tcmet30_phimax110 ()
{
     return run<Looper>(  (CUT_BIT(ELE_PT_20)) | (CUT_BIT(ELE_ISO_10) | CUT_BIT(EVT_JPT_PHIMAX_110) | CUT_BIT(EVT_TCMET_30)),
        "Results_iso10_jptphimax110_tcmet30",
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}

int Results_tcmet30_phimax100_conv ()
{
     return run<Looper>(  (CUT_BIT(ELE_PT_20)) | (CUT_BIT(ELE_ISO_10) | CUT_BIT(EVT_JPT_PHIMAX_100) | CUT_BIT(EVT_TCMET_30) | CUT_BIT(ELE_NOCONV)),
        "Results_iso10_jptphimax100_tcmet30_conv",
1 << LOOP_WENU_7TeV | 1 << LOOP_QCD30_7TeV | 1 << LOOP_PHOTONJET_7TeV);
}



