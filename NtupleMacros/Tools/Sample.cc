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
     return true;
}

//WW file
Sample fWW ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-01/WW_2l-Pythia/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "/data/tmp/jmuelmen/2_2_ntuple_test/ntuple_*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, WW, kRed, 1, "ww", true };
     return ret;
}

//WZ file
Sample fWZ ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-01/WZ_3l-Pythia/merged_ntuple*.root";
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
     std::string sample = "/data/tmp/cms2-V01-02-01/ZZ_2l2n-Pythia/merged_ntuple*.root";
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
     std::string sample = "/data/tmp/cms2-V01-02-01/WJets-madgraph/merged_ntuple*.root";
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
     std::string sample = "/data/tmp/cms2-V01-02-01/ZJets-madgraph/merged_ntuple*.root";
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
     std::string sample = "/data/tmp/cms2-V01-02-01/ZJets-madgraph/merged_ntuple*.root";
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
     std::string sample = "/data/tmp/cms2-V01-02-01/ZJets-madgraph/merged_ntuple*.root";
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
     std::string sample = "/data/tmp/cms2-V01-02-01/TTJets-madgraph/merged_ntuple*.root";
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
     std::string sample = "/data/tmp/cms2-V01-02-01/SingleTop_tWChannel-madgraph-LHE/merged_ntuple*.root";
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
