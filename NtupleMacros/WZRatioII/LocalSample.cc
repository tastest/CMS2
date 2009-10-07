//#include "TChain.h"
#include "LocalSample.h"
//#include "../Tools/tools.h"
//#include "../CORE/selections.h"
//#include "../CORE/CMS2.h"
//#include <cstdlib>
#include <string>
#include <iostream>

static const std::string prefix = (getenv("CMS2_NTUPLE_LOCATION") != 0) ?
     std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" : "/data/tmp/";

/*
// Z+0 jet
Sample fZ_0Jet ()
{    
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/Z_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged*.root";  
     c->Add(sample.c_str());
     Sample ret = { c, OTHER, 42, 1, "z_0jet", true, 0. };
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


*/
