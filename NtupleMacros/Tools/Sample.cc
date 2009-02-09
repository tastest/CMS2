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
     default:
	  return true;
     }
     return true;
}

//WW file
Sample fWW ()
{
     TChain *c = new TChain("Events");
//     std::string sample = "/data/tmp/cms2-V01-02-01/WW_2l-Pythia/merged_ntuple*.root";
     std::string sample = "/data/tmp/cms2-V01-02-06/WW_2l_Summer08_IDEAL_V9_v2/postprocessing/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
//       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/WW_2l-Pythia/merged_ntuple*.root";
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/WW_2l_Summer08_IDEAL_V9_v2/postprocessing/merged_ntuple*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, WW, kRed, 1, "ww", true };
     return ret;
}

//WZ file
Sample fWZ ()
{
     TChain *c = new TChain("Events");
//     std::string sample = "/data/tmp/cms2-V01-02-01/WZ_3l-Pythia/merged_ntuple*.root";
     std::string sample = "/data/tmp/cms2-V01-02-06/ZZ_2l2n_Summer08_IDEAL_V9_v2/postprocessing/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
//       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/WZ_3l-Pythia/merged_ntuple*.root";
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/ZZ_2l2n_Summer08_IDEAL_V9_v2/postprocessing/merged_ntuple*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, WZ, kBlue, 1, "wz", true };
     return ret;
}

//WZ file
Sample fWZ_incl ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-01/WZ_incl-Pythia/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/WZ_incl-Pythia/merged_ntuple*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, WZ, kBlue, 1, "wzincl", true };
     return ret;
}

//ZZ file
Sample fZZ ()
{
     TChain *c = new TChain("Events");
//     std::string sample = "/data/tmp/cms2-V01-02-01/ZZ_2l2n-Pythia/merged_ntuple*.root";
     std::string sample = "/data/tmp/cms2-V01-02-06/ZZ_2l2n_Summer08_IDEAL_V9_v2/postprocessing/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
//       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/ZZ_2l2n-Pythia/merged_ntuple*.root";
        sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/ZZ_2l2n_Summer08_IDEAL_V9_v2/postprocessing/merged_ntuple*.root";
    }
     c->Add(sample.c_str());
     Sample ret = { c, ZZ, kGreen, 1, "zz", true };
     return ret;
}

//Wjets file
Sample fWjets ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-01/WJets-madgraph/postprocessing/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/WJets-madgraph/postprocessing/merged_ntuple*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, Wjets, 40, 1, "wjets", true };
     return ret;
}

//DYee file
Sample fDYee ()
{
     TChain *c = new TChain("Events");
//     std::string sample = "/data/tmp/cms2-V01-02-01/ZJets-madgraph/merged_ntuple*.root";
     std::string sample = "/data/tmp/cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/postprocessing/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
//       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/ZJets-madgraph/merged_ntuple*.root";
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/postprocessing/merged_ntuple*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kMagenta, 1.12, "dyee", true };
     return ret;
}

//DYmm file
Sample fDYmm ()
{
     TChain *c = new TChain("Events");
//     std::string sample = "/data/tmp/cms2-V01-02-01/ZJets-madgraph/merged_ntuple*.root";
     std::string sample = "/data/tmp/cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/postprocessing/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
//       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/ZJets-madgraph/merged_ntuple*.root";
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/postprocessing/merged_ntuple*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kCyan, 1.12, "dymm", true };
     return ret;
}

//DYtt file
Sample fDYtt ()
{
     TChain *c = new TChain("Events");
//     std::string sample = "/data/tmp/cms2-V01-02-01/ZJets-madgraph/merged_ntuple*.root";
     std::string sample = "/data/tmp/cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/postprocessing/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
//       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/ZJets-madgraph/merged_ntuple*.root";
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/postprocessing/merged_ntuple*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, DYtt, kBlack, 1.12, "dytt", true };
     return ret;
}

//ttbar file
Sample fttbar ()
{
     TChain *c = new TChain("Events");
//     std::string sample = "/data/tmp/cms2-V01-02-01/TTJets-madgraph/merged_ntuple*.root";
     std::string sample = "/data/tmp/cms2-V01-02-06/TTJets-madgraph_Fall08_IDEAL_V9_v1/postprocessing/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
//       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/TTJets-madgraph/merged_ntuple*.root";
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/TTJets-madgraph_Fall08_IDEAL_V9_v1/postprocessing/merged_ntuple*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, ttbar, kYellow, 1.85, "ttbar", true };
     return ret;
}

//ttbar file
Sample fttbar_taula ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-01/TauolaTTbar-Pythia/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/TauolaTTbar-Pythia/merged_ntuple*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, ttbar, kYellow, 1.85, "ttbartauola", true };
     return ret;
}

Sample ftW ()
{
     TChain *c = new TChain("Events");
//     std::string sample = "/data/tmp/cms2-V01-02-01/SingleTop_tWChannel-madgraph-LHE/merged_ntuple*.root";
     std::string sample = "/data/tmp/cms2-V01-02-06/SingleTop_tWChannel-madgraph-LHE/postprocessing/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
//       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/SingleTop_tWChannel-madgraph-LHE/merged_ntuple*.root";
	sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "/cms2-V01-02-06/SingleTop_tWChannel-madgraph-LHE/postprocessing/merged_ntuple*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, tW, 63, 1, "tw", true };
     return ret;
}

Sample fSingleTop_tChannel ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-01/SingleTop_tChannel-madgraph-LHE/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/SingleTop_tChannel-madgraph-LHE/merged_ntuple*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, tW, 63, 1, "singletopt", true };
     return ret;
}

Sample fSingleTop_sChannel ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-01/SingleTop_sChannel-madgraph-LHE/merged_ntuple*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-01/SingleTop_sChannel-madgraph-LHE/merged_ntuple*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, tW, 63, 1, "singletops", true };
     return ret;
}


// still old, has to be changed
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

// still old, has to be changed
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

// still old, has to be changed
Sample fLM4 ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V00-05-00/LM4/*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V00-05-00/LM4/*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, LM4, 38, 1, "LM4", false };
     return ret;
}

// QCD samples
Sample fInclusiveMu5Pt50 ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-06/InclusiveMu5Pt50/*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/InclusiveMu5Pt50/*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, InclusiveMu5Pt50, 28, 1, "InclusiveMu5Pt50", true };
     return ret;
}

Sample fInclusiveMuPt15 ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-06/InclusiveMuPt15/*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/InclusiveMuPt15/*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, InclusiveMuPt15, 28, 1, "InclusiveMuPt15", true };
     return ret;
}

Sample fQCDBCtoEPt20to30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-06/QCD_BCtoE_Pt20to30/*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/QCD_BCtoE_Pt20to30/*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt20to30, 28, 1, "QCDBCtoEPt20to30", true };
     return ret;
}

Sample fQCDBCtoEPt30to80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-06/QCD_BCtoE_Pt30to80/*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/QCD_BCtoE_Pt30to80/*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt30to80, 28, 1, "QCDBCtoEPt30to80", true };
     return ret;
}

Sample fQCDBCtoEPt80to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-06/QCD_BCtoE_Pt80to170/*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/QCD_BCtoE_Pt80to170/*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt80to170, 28, 1, "QCDBCtoEPt80to170", true };
     return ret;
}

Sample fQCDEMenrichedPt20to30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-06/QCD_EMenriched_Pt20to30/*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/QCD_EMenriched_Pt20to30/*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt20to30, 28, 1, "QCDEMenrichedPt20to30", true };
     return ret;
}

Sample fQCDEMenrichedPt30to80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-06/QCD_EMenriched_Pt30to80/*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/QCD_EMenriched_Pt30to80/*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt30to80, 28, 1, "QCDEMenrichedPt30to80", true };
     return ret;
}

Sample fQCDEMenrichedPt80to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/data/tmp/cms2-V01-02-06/QCD_EMenriched_Pt80to170/*.root";
     char *location = getenv("CMS2_NTUPLE_LOCATION");
     if ( location != 0 ) {
       sample = std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" + "cms2-V01-02-06/QCD_EMenriched_Pt80to170/*.root";
     }
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt80to170, 28, 1, "QCDEMenrichedPt80to170", true };
     return ret;
}

