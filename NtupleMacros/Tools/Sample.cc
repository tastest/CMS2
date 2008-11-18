#include "TChain.h"
#include "Sample.h"
#include "utilities.h"
#include "../CORE/selections.h"
#include "../CORE/CMS2.h"
#include <cstdlib>
#include <string>
#include <iostream>

bool filterByProcess (enum Process sample)
{
     switch (sample) {
     case WW: case WZ: case ZZ: case tW: case LM1: case LM2: case LM4:
	  return true;
     case Wjets:
	  return cms2.evt_CSA07Process() < 11;
     case DYee: 
	  return (cms2.evt_CSA07Process() > 10 && cms2.evt_CSA07Process() < 22 && isDYee() );
     case DYmm:
	  return (cms2.evt_CSA07Process() > 10 && cms2.evt_CSA07Process() < 22 && isDYmm() );
     case DYtt:
	  return (cms2.evt_CSA07Process() > 10 && cms2.evt_CSA07Process() < 22 && isDYtt() );
     case ttbar:
	  return (cms2.evt_CSA07Process() > 21 && cms2.evt_CSA07Process() < 27);
     }
     return false;
}

//WW file
Sample fWW ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-04-00/merge_WW.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-04-00/merge_WW.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, WW, kRed, 1, "ww", true };
     return ret;
}

//WZ file
Sample fWZ ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-05-00/merge_WZ.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-05-00/merge_WZ.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, WZ, kBlue, 1, "wz", true };
     return ret;
}

//ZZ file
Sample fZZ ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-05-00/merge_ZZ.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-05-00/merge_ZZ.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, ZZ, kGreen, 1, "zz", true };
     return ret;
}

//Wjets file
Sample fWjets ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-04-01/merge_Wjet.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-04-01/merge_Wjet.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, Wjets, 40, 1, "wjets", true };
     return ret;
}

//DYee file
Sample fDYee ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-04-01/merge_DY.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-04-01/merge_DY.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kMagenta, 1.12, "dyee", true };
     return ret;
}

//DYmm file
Sample fDYmm ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-04-01/merge_DY.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-04-01/merge_DY.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kCyan, 1.12, "dymm", true };
     return ret;
}

//DYtt file
Sample fDYtt ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-04-01/merge_DY.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-04-01/merge_DY.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, DYtt, kBlack, 1.12, "dytt", true };
     return ret;
}

//ttbar file
Sample fttbar ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-04-01/merge_ttbar.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-04-01/merge_ttbar.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, ttbar, kYellow, 1.85, "ttbar", true };
     return ret;
}

Sample ftW ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-04-00/merge_tW.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-04-00/merge_tW.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, tW, 63, 1, "tw", true };
     return ret;
}

Sample fLM1 ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-05-00/LM1/*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-05-00/LM1/*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, LM1, 37, 1, "LM1", false };
     return ret;
}

Sample fLM2 ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-05-00/LM2/*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-05-00/LM2/*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, LM2, 38, 1, "LM2", false };
     return ret;
}

Sample fLM4 ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-05-00/LM4/*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-05-00/LM4/*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, LM4, 28, 1, "LM4", false };
     return ret;
}
