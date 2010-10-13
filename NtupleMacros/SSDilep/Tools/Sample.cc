#include "TChain.h"
#include "Sample.h"
#include "tools.h"
#include "../CORE/selections.h"
#include "../CORE/CMS2.h"
#include <cstdlib>
#include <string>
#include <iostream>

bool filterByProcess (enum Process sample)
{
     switch (sample) {
     case DYee: 
          return isDYee();
     case DYmm:
          return isDYmm();
     case DYtt:
          return isDYtt();
     case We: 
          return isWe();
     case Wm:
          return isWm();
     case Wt:
          return isWt();
     default:
	  return true;
     }
     return true;
}

Sample operator + (const Sample &a, const Sample &b)
{
     Sample ret = a;
     a.chain->GetEntries();
     ret.chain = dynamic_cast<TChain *>(a.chain->Clone());
     ret.chain->GetEntries();
//      printf("cloning chain: %llu (%d files) to %llu entries (%d files)\n", 
// 	    a.chain->GetEntries(), a.chain->GetListOfFiles()->GetEntries(),
// 	    ret.chain->GetEntries(), ret.chain->GetListOfFiles()->GetEntries());
     b.chain->GetEntries();      // ha ha, if you don't do this, the
				 // combined chain will have a random
				 // number of entries.
     ret.chain->Add(dynamic_cast<TChain *>(b.chain->Clone()));
//      printf("adding %llu (%d files), returned chain now has %llu entries (%d files)\n", 
// 	    b.chain->GetEntries(), b.chain->GetListOfFiles()->GetEntries(),
// 	    a.chain->GetEntries(), a.chain->GetListOfFiles()->GetEntries());
     return ret;
}

Sample fFreeForm (const char *sample_glob, enum Process p, int histo_color, 
		  double kFactor, std::string name, bool sm, double pthat)
{
     TChain *c = makeChain(sample_glob);
     Sample ret = { c, p, histo_color, kFactor, name, sm, pthat };
     return ret;
     
}

static const std::string prefix = (getenv("CMS2_NTUPLE_LOCATION") != 0) ? 
     std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" : "/data/tmp/";

//WW file
Sample fWW ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/WW_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*root";
     c->Add(sample.c_str());
     Sample ret = { c, WW, kRed, 1, "ww", true, 0. };
     return ret;
}

//WZ file
Sample fWZ ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/WZ_Summer09-MC_31X_V3_7TeV-v1_v2/V03-00-35/merge*root";
     c->Add(sample.c_str());
     Sample ret = { c, WZ, kBlue, 1, "wz", true, 0. };
     return ret;
}

//ZZ file
Sample fZZ ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/ZZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*root";
     c->Add(sample.c_str());
     Sample ret = { c, ZZ, kGreen, 1, "zz", true, 0. };
     return ret;
}

Sample fWe ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/Wenu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*root";
     c->Add(sample.c_str());

     Sample ret = { c, Wjets, 33, 1, "we", true, 0. };
     return ret;
}

Sample fWm ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/Wmunu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*root";
     c->Add(sample.c_str());

     Sample ret = { c, Wjets, 36, 1, "wm", true, 0. };
     return ret;
}

Sample fWt ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/Wtaunu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*root";
     c->Add(sample.c_str());

     Sample ret = { c, Wjets, 39, 1, "wt", true, 0. };
     return ret;
}

Sample fDYee ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/Zee_Summer09-MC_31X_V3_7TeV_TrackingParticles-v1/V03-00-35/merge*root";
     c->Add(sample.c_str());

     Sample ret = { c, DYee, 42, 1, "dyee", true, 0. };
     return ret;
}

Sample fDYmm ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/Zmumu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*root";
     c->Add(sample.c_str());

     Sample ret = { c, DYmm, 46, 1, "dymm", true, 0. };
     return ret;
}

Sample fDYtt ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/Ztautau_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*root";
     c->Add(sample.c_str());

     Sample ret = { c, DYtt, 49, 1, "dytt", true, 0. };
     return ret;
}

Sample fttbar ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/TTBar_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*root";
     c->Add(sample.c_str());
     Sample ret = { c, ttbar, kYellow, 1, "ttbar", true, 0. };
     return ret;
}

Sample fQCDpt30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/QCD_Pt30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt30, 29, 1, "QCDpt30", true, 999999999 };
     return ret;
}

Sample fInclusiveMu5Pt30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/InclusiveMu5_Pt30_Summer09-MC_31X_V3_7TeV-v1_2xv41m/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, InclusiveMuPt15, 30, 1, "InclusiveMuPt15", true, 0. };
     return ret;
}

Sample fQCDBCtoEPt30to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/QCD_BCtoE_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "cms2/QCD_BCtoE_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt30to170, 28, 1, "QCDBCtoEPt30to170", true, 0. };
     return ret;
}

Sample fQCDBCtoEPt20to170 ()
{    
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/QCD_BCtoE_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "cms2/QCD_BCtoE_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "cms2/QCD_BCtoE_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt20to170, 28, 1, "QCDBCtoEPt20to170", true, 0. };
     return ret;
}    


Sample fQCDBCtoEPt20to30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/QCD_BCtoE_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt20to30, 28, 1, "QCDBCtoEPt20to30", true, 0. };
     return ret;
}

Sample fQCDBCtoEPt30to80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/QCD_BCtoE_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt30to80, 28, 1, "QCDBCtoEPt30to80", true, 0. };
     return ret;
}

Sample fQCDBCtoEPt80to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/QCD_BCtoE_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt80to170, 28, 1, "QCDBCtoEPt80to170", true, 0. };
     return ret;
}

Sample fQCDEMenrichedPt30to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/QCD_EMEnriched_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "cms2/QCD_EMEnriched_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt30to170, 28, 1, "QCDEMenrichedPt30to170", true, 0. };
     return ret;
}

Sample fQCDEMenrichedPt20to170 ()
{     
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/QCD_EMEnriched_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "cms2/QCD_EMEnriched_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "cms2/QCD_EMEnriched_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt20to170, 28, 1, "QCDEMenrichedPt20to170", true, 0. };
     return ret;
}

Sample fQCDEMenrichedPt20to30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/QCD_EMEnriched_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt20to30, 28, 1, "QCDEMenrichedPt20to30", true, 0. };
     return ret;
}

Sample fQCDEMenrichedPt30to80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/QCD_EMEnriched_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt30to80, 28, 1, "QCDEMenrichedPt30to80", true, 0. };
     return ret;
}

Sample fQCDEMenrichedPt80to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/QCD_EMEnriched_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt80to170, 28, 1, "QCDEMenrichedPt80to170", true, 0. };
     return ret;
}

// photon + jet samples

Sample fPhotonJetPt20to170 ()
{    
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/PhotonJet_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "cms2/PhotonJet_Pt30to50_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "cms2/PhotonJet_Pt50to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "cms2/PhotonJet_Pt80to120_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "cms2/PhotonJet_Pt120to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, PhotonJet, 28, 1, "PhotonJetPt20to170", true, 0. };
     return ret;
}    

Sample fPhotonJetPt20to30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/PhotonJet_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, PhotonJet, 28, 1, "PhotonJetPt20to30", true, 0. };
     return ret;
}

Sample fPhotonJetPt30to50 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/PhotonJet_Pt30to50_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, PhotonJet, 28, 1, "PhotonJetPt30to50", true, 0. };
     return ret;
}

Sample fPhotonJetPt50to80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/PhotonJet_Pt50to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, PhotonJet, 28, 1, "PhotonJetPt50to80", true, 0. };
     return ret;
}

Sample fPhotonJetPt80to120 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/PhotonJet_Pt80to120_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, PhotonJet, 28, 1, "PhotonJetPt80to120", true, 0. };
     return ret;
}

Sample fPhotonJetPt120to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2/PhotonJet_Pt120to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, PhotonJet, 28, 1, "PhotonJetPt120to170", true, 0. };
     return ret;
}

