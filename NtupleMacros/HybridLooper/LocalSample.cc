#include "TChain.h"
#include "LocalSample.h"
#include "../Tools/tools.h"
#include "../CORE/selections.h"
#include "../CORE/CMS2.h"
#include <cstdlib>
#include <string>
#include <iostream>

static const std::string prefix = (getenv("CMS2_NTUPLE_LOCATION") != 0) ?
     std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" : "/data/tmp/";


Sample fPhotonJet_7TeV()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "/cms2-V02-00-08/PhotonJet_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "/cms2-V02-00-08/PhotonJet_Pt30to50_Summer09-MC_31X_V3_7TeV-v1/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "/cms2-V02-00-08/PhotonJet_Pt50to80_Summer09-MC_31X_V3_7TeV-v1/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "/cms2-V02-00-08/PhotonJet_Pt80to120_Summer09-MC_31X_V3_7TeV-v1/merged*.root";
     c->Add(sample.c_str());
     sample = prefix + "/cms2-V02-00-08/PhotonJet_Pt120to170_Summer09-MC_31X_V3_7TeV-v1/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, OTHER, kMagenta, 1, "photonjet", true, 0. };
     return ret;
}

Sample fQCD_Pt30_7TeV()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "/cms2-V02-00-08/QCD_Pt30_Summer09-MC_31X_V3_7TeV-v1/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, OTHER, kGreen, 1, "qcd_pt30", true, 0. };
     return ret;
}

Sample fWenu_7TeV()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "/cms2-V02-00-08/Wenu_Summer09-MC_31X_V3_7TeV-v1/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, OTHER, kBlue, 1, "wenu", true, 0. };
     return ret;
}

Sample fZZ_7TeV()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "/cms2-V02-00-08/ZZ_Summer09-MC_31X_V3_7TeV-v1/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, OTHER, kBlue, 1, "zz", true, 0. };
     return ret;
}


// validation sample
Sample fValidation()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "/dlevans/wenu_postprocessed.root";
     c->Add(sample.c_str());
     Sample ret = { c, OTHER, kBlue, 1, "wenu", true, 0. };
     return ret;
}

// Wenu
Sample fWenu ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/ESTest_v0/Wenu_Summer09-MC_31X_V2_preproduction_311-v1/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, OTHER, 40, 1, "wenu", true, 0. };
     return ret;
}

// QCD EM Enriched 30-80
Sample fEM30_80()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/ESTest_v0/QCD_EMEnriched_Pt30to80_Summer09-MC_31X_V2_preproduction_311-v1/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, OTHER, 28, 1, "em30_80", true, 0. };
     return ret;
}

// BC->E 30-80
Sample fBC30_80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/ESTest_v0/QCD_BCtoE_Pt30to80_Summer09-MC_31X_V2_preproduction_311-v1/merged*.root";
     c->Add(sample.c_str());
     Sample ret = { c, OTHER, 29, 1, "bc30_80", true, 0. };
     return ret;
}


