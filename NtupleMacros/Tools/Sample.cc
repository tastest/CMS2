#include "TChain.h"
#include "Sample.h"
#include "tools.h"
#include "CORE/mcSelections.h"
#include <string>

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
     std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" : "/nfs-3/userdata/cms2/";

struct file {
     char name[1024];
     char location[4096];
     Sample sample;
};

static const file files[] = {
     { "WW",		"WW_Spring10-START3X_V26_S09-v1_DiLep/V03-04-08/merged_ntuple*.root", 		{ 0, WW, kRed, 		1, "ww", 	true, 	0. } },
     { "WZ",		"WZ_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*.root", 		{ 0, WZ, kBlue, 	1, "wz", 	true, 	0. } },
     { "ZZ",		"ZZ_Spring10-START3X_V26_S09-v1_DiLep/V03-04-08/merged_ntuple*.root", 		{ 0, ZZ, kGreen,	1, "zz", 	true, 	0. } },
     { "We", 		"WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-13-07/diLepPt2010Skim/*.root", 		{ 0, We, 33, 		31314./28049., "we", 	true, 	0. } },
     { "Wmu",		"WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-13-07/diLepPt2010Skim/*.root", 		{ 0, Wm, 36, 		31314./28049., "wmu", 	true, 	0. } },
     { "Wtau",		"WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-13-07/diLepPt2010Skim/*.root", 		{ 0, Wt, 39, 		31314./28049., "wtau", 	true, 	0. } },
     { "Zee",		"ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/diLepPt2010Skim/*.root",		{ 0, DYee, 42, 		1, "zee", 	true, 	0. } },
     { "Zmm",		"ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/diLepPt2010Skim/*.root",		{ 0, DYmm, 46, 		1, "zmm", 	true, 	0. } },
     { "Ztt",		"ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/diLepPt2010Skim/*.root",		{ 0, DYtt, 49, 		1, "ztt", 	true, 	0. } },
     { "VV",		"VVJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/*.root",		{ 0, VV, kRed, 		1, "vv", 	true, 	0. } },
     { "ttbar",		"TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/diLepPt2010Skim/*.root",		{ 0, ttbar, kYellow, 	157.5/165., "ttbar", 	true, 	0. } },
     { "tW",		"SingleTop_tWChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/*.root",	{ 0, tW, 63, 		1, "tw", 	true, 	0. } },
     { "LM8",		"LM8_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged_ntuple*.root", 		{ 0, LM8, 37, 		1, "LM8", 	false, 	0. } },
};
static const int files_size = sizeof(files) / sizeof(file);

const file *find_file (const char *name)
{
     for (int i = 0; i < files_size; ++i) {
	  if (strcmp(files[i].name, name) == 0)
	       return &files[i];
     }
     return 0;
}

inline Sample f_ (const char *x) 
{ 
     const file *f = find_file(x);
     assert(f != 0);
     Sample ret = f->sample;
     ret.chain = makeChain((prefix + f->location).c_str());
     return ret;
}

// New test file
Sample fTest ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/home/dlevans/eleID/CMSSW_2_2_10/src/CMS2/NtupleMaker/test/postprocessed_ntuple.root";
     c->Add(sample.c_str());
     Sample ret = { c, TEST, kRed, 1, "test", true, 0. };
     return ret;
}

#define f__(x) Sample f##x () { return f_(#x); }

f__(WW)
f__(WZ)
f__(ZZ)
f__(VV)

f__(We)
f__(Wmu)
f__(Wtau)

Sample fdata ()
{
     std::string prefix_ = std::string(" ") + prefix; // prepend prefix with space so we can stack paths together 
     std::string data = 
	  prefix_ + "MinimumBias_Commissioning10-SD_EG-Jun14thSkim_v1_RECO/V03-04-26-02/singleLepPt10Skim/*.root" +
	  prefix_ + "MinimumBias_Commissioning10-SD_Mu-Jun14thSkim_v1_RECO/V03-04-26-02/merged*.root" +
	  
	  prefix_ + "EG_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/singleLepPt5Skim/*.root" +
	  prefix_ + "Mu_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/singleLepPt5Skim/*.root" +      
	  
	  //ADD THE JULY 16th RERECO HERE BEFORE THE PROMPT RECO
	  " /nfs-3/userdata/slava77/EG_Run2010A-Jul16thReReco-v2_RECO/V03-04-26-03/dilep1010/skimlep1010*.root" +
	  prefix_ + "Mu_Run2010A-Jul16thReReco-v1/V03-04-26-02/merged_ntuple*.root" +
	  
	  prefix_ + "EG_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/*.root" +
	  prefix_ + "EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/*.root" +            
	  prefix_ + "Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/*.root" +
	  prefix_ + "Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/*.root";
     printf("%s\n", data.c_str());
//      TChain *c = makeChain(data.c_str());
     TChain *c = makeChain("/home/users/jmuelmen/CMSSW_3_6_1_patch4/src/CMS2/NtupleMacros/NtupleTools/dilep_skim_2.root");
//     TChain *c = makeChain("/nfs-3/userdata/fgolf/SSskims/data/*.root");
//	 c->Add("/tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-12/singleLepPt10Skim/skimmed_ntuple_1425*.root");
//	 c->Add("/tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-12/singleLepPt10Skim/skimmed_ntuple_1426*.root");			   
//	 c->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-12/singleLepPt10Skim/skimmed_ntuple_1425*.root");
//	 c->Add("/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-12/singleLepPt10Skim/skimmed_ntuple_1426*.root");
     Sample ret = { c, DATA, kBlack, 1, "data", true, 0. };
     return ret;
}

//WeJets file
Sample fWeJets ()
{
     TChain *c = makeChain((prefix + "WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/dilep-skim/*skim*.root").c_str());
//     std::string sample = prefix + "WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root";
     Sample ret = { c, We, 33, 1, "wejets", true, 0. };
     return ret;
}

//WmJets file
Sample fWmJets ()
{
     TChain *c = makeChain((prefix + "WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/dilep-skim/*skim*.root").c_str());
//     std::string sample = prefix + "WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root";
     Sample ret = { c, Wm, 36, 1, "wmjets", true, 0. };
     return ret;
}

//WtJets file
Sample fWtJets ()
{
     TChain *c = makeChain((prefix + "WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/dilep-skim/*skim*.root").c_str());
//     std::string sample = prefix + "WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root";
     Sample ret = { c, Wt, 39, 1, "wtjets", true, 0. };
     return ret;
}

f__(Zee)
f__(Zmm)
f__(Ztt)

//Zee Madgraph file
Sample fZeejets ()
{
     std::string prefix_ = std::string(" ") + prefix; // prepend prefix with space so we can stack paths together 
     std::string str = 
	  prefix_ + "EarlyDataSamples/Zee_Spring10-START3X_V26_S09-v1/*.root" + 
	  prefix_ + "EarlyDataSamples/DYee_M10to20_Spring10-START3X_V26_S09-v1/V03-04-08-01/*.root" +
	  prefix_ + "ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08/diLepPt2020Skim/skimmed_ntuple.root";
     TChain *c = makeChain(str.c_str());
     Sample ret = { c, DYee, 42, 1, "zeejets", true, 0. };
     return ret;
}

//Zmm Madgraph file
Sample fZmmjets ()
{
     std::string prefix_ = std::string(" ") + prefix; // prepend prefix with space so we can stack paths together 
     std::string str = 
	  prefix_ + "EarlyDataSamples/Zmumu_Spring10-START3X_V26_S09-v1/*.root" +
	  prefix_ + "EarlyDataSamples/DYmumu_M10to20_Spring10-START3X_V26_S09-v1/*.root" + 
	  prefix_ + "ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08/diLepPt2020Skim/skimmed_ntuple.root";
     TChain *c = makeChain(str.c_str());
     Sample ret = { c, DYmm, 46, 1, "zmmjets", true, 0. };
     return ret;
}

//Ztt Madgraph file
Sample fZttjets ()
{
     std::string prefix_ = std::string(" ") + prefix; // prepend prefix with space so we can stack paths together 
     std::string str = 
	  prefix_ + "EarlyDataSamples/Ztautau_Spring10-START3X_V26_S09-v1/*.root" +
	  prefix_ + "ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08/diLepPt2020Skim/skimmed_ntuple.root";
     TChain *c = makeChain(str.c_str());
     Sample ret = { c, DYtt, 49, 1, "zttjets", true, 0. };
     return ret;
}

f__(ttbar)
f__(tW)

//ttbar madgraph file
Sample fttMadgraph ()
{
     TChain *c = makeChain((prefix + "TTbarJets-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, ttMadgraph, kYellow, 1, "ttjets", true, 0. };
     return ret;
}

//single top t channel file
Sample fSingleTop_tChannel ()
{
     TChain *c = makeChain((prefix + "SingleTop_tChannel-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, singleTop_tChannel, 63, 1, "singletopt", true, 0. };
     return ret;
}

//single top s channel file
Sample fSingleTop_sChannel ()
{
     TChain *c = makeChain((prefix + "SingleTop_sChannel-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, singleTop_tChannel, 63, 1, "singletops", true, 0. };
     return ret;
}

//LM0 file
Sample fLM0 ()
{
     TChain *c = makeChain((prefix + "LM0_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM0, 37, 1, "LM0", false, 0. };
     return ret;
}

//LM1 file
Sample fLM1 ()
{
     TChain *c = makeChain((prefix + "LM1_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, LM1, 37, 1, "LM1", false, 0. };
     return ret;
}

//LM2 file
Sample fLM2 ()
{
     TChain *c = makeChain((prefix + "LM2_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM2, 37, 1, "LM2", false, 0. };
     return ret;
}

//LM2mhfeq360 file
Sample fLM2mhfeq360 ()
{
     TChain *c = makeChain((prefix + "LM2mhfeq360_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM2mhfeq360, 37, 1, "LM2mhfeq360", false, 0. };
     return ret;
}

//LM3 file
Sample fLM3 ()
{
     TChain *c = makeChain((prefix + "LM3_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM3, 37, 1, "LM3", false, 0. };
     return ret;
}

//LM4 file
Sample fLM4 ()
{
     TChain *c = makeChain((prefix + "LM4_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM4, 37, 1, "LM4", false, 0. };
     return ret;
}

//LM5 file
Sample fLM5 ()
{
     TChain *c = makeChain((prefix + "LM5_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM5, 37, 1, "LM5", false, 0. };
     return ret;
}

//LM6 file
Sample fLM6 ()
{
     TChain *c = makeChain((prefix + "LM6_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM6, 37, 1, "LM6", false, 0. };
     return ret;
}

//LM7 file
Sample fLM7 ()
{
     TChain *c = makeChain((prefix + "LM7_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM7, 37, 1, "LM7", false, 0. };
     return ret;
}

f__(LM8)

//LM9 file
Sample fLM9 ()
{
     TChain *c = makeChain((prefix + "LM9_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM9, 37, 1, "LM9", false, 0. };
     return ret;
}

//LM9p file
Sample fLM9p ()
{
     TChain *c = makeChain((prefix + "LM9p_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM9p, 37, 1, "LM9p", false, 0. };
     return ret;
}

//LM9t175 file
Sample fLM9t175 ()
{
     TChain *c = makeChain((prefix + "LM9t175_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM9t175, 37, 1, "LM9t175", false, 0. };
     return ret;
}

//LM10 file
Sample fLM10 ()
{
     TChain *c = makeChain((prefix + "LM10_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM10, 37, 1, "LM10", false, 0. };
     return ret;
}

//LM11 file
Sample fLM11 ()
{
     TChain *c = makeChain((prefix + "LM11_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM11, 37, 1, "LM11", false, 0. };
     return ret;
}

//LM12 file
Sample fLM12 ()
{
     TChain *c = makeChain((prefix + "LM12_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM12, 37, 1, "LM12", false, 0. };
     return ret;
}

//LM13 file
Sample fLM13 ()
{
     TChain *c = makeChain((prefix + "LM13_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merge*.root").c_str());
     Sample ret = { c, LM13, 37, 1, "LM13", false, 0. };
     return ret;
}

//InclusiveMuPt15 file
Sample fInclusiveMuPt15 ()
{
     TChain *c = makeChain((prefix + "InclusiveMu15_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, InclusiveMuPt15, 30, 1, "InclusiveMuPt15", true, 0. };
     return ret;
}

//InclusiveMuPt15 dilepton filtered file
Sample fInclusiveMuPt15dilep ()
{
     TChain *c = makeChain((prefix + "InclusiveMu15_Summer09-MC_31X_V3_7TeV-v1_dilepfilt/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, InclusiveMuPt15dilep, 30, 1, "InclusiveMuPt15", true, 0. };
     return ret;
}

//QCDBCtoEPt20to30 file
Sample fQCDBCtoEPt20to30 ()
{
     TChain *c = makeChain((prefix + "QCD_BCtoE_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDBCtoEPt20to30, 28, 1, "QCDBCtoEPt20to30", true, 0. };
     return ret;
}

//QCDBCtoEPt30to80 file
Sample fQCDBCtoEPt30to80 ()
{
     TChain *c = makeChain((prefix + "QCD_BCtoE_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDBCtoEPt30to80, 28, 1, "QCDBCtoEPt30to80", true, 0. };
     return ret;
}

//QCDBCtoEPt80to170 file
Sample fQCDBCtoEPt80to170 ()
{
     TChain *c = makeChain((prefix + "QCD_BCtoE_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDBCtoEPt80to170, 28, 1, "QCDBCtoEPt80to170", true, 0. };
     return ret;
}

// combination of QCDBCtoE samples
Sample fQCDBCtoEPt20to170 ()
{    
     TChain *c = makeChain((prefix + "QCD_BCtoE_Pt{20to30,30to80,80to170}_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDBCtoEPt20to170, 28, 1, "QCDBCtoEPt20to170", true, 0. };
     return ret;
}    

//QCDEMenrichedPt20to30 file
Sample fQCDEMenrichedPt20to30 ()
{
     TChain *c = makeChain((prefix + "QCD_EMEnriched_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDEMenrichedPt20to30, 28, 1, "QCDEMenrichedPt20to30", true, 0. };
     return ret;
}

//QCDEMenrichedPt30to80 file
Sample fQCDEMenrichedPt30to80 ()
{
     TChain *c = makeChain((prefix + "QCD_EMEnriched_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDEMenrichedPt30to80, 28, 1, "QCDEMenrichedPt30to80", true, 0. };
     return ret;
}

//QCDEMenrichedPt80to170 file
Sample fQCDEMenrichedPt80to170 ()
{
     TChain *c = makeChain((prefix + "QCD_EMEnriched_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDEMenrichedPt80to170, 28, 1, "QCDEMenrichedPt80to170", true, 0. };
     return ret;
}

// combination of QCDEMentriched samples
Sample fQCDEMenrichedPt20to170 ()
{     
     TChain *c = makeChain((prefix + "QCD_EMEnriched_Pt{20to30,30to80,80to170}_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDEMenrichedPt20to170, 28, 1, "QCDEMenrichedPt20to170", true, 0. };
     return ret;
}

//QCDpt15 file
Sample fQCDpt15 ()
{
     TChain *c = makeChain((prefix + "QCD_Pt15_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDpt15, 28, 1, "QCDpt15", true, 30 };
     return ret;
}

//QCDpt30 file
Sample fQCDpt30 ()
{
     TChain *c = makeChain((prefix + "QCD_Pt30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDpt30, 28, 1, "QCDpt30", true, 80 };
     return ret;
}

//QCDpt80 file
Sample fQCDpt80 ()
{
     TChain *c = makeChain((prefix + "QCD_Pt80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDpt80, 28, 1, "QCDpt80", true, 170 };
     return ret;
}

//QCDpt170 file
Sample fQCDpt170 ()
{
     TChain *c = makeChain((prefix + "QCD_Pt170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDpt170, 28, 1, "QCDpt170", true, 300 };
     return ret;
}

//QCDpt1400 file
Sample fQCDpt1400 ()
{
     TChain *c = makeChain((prefix + "QCD_Pt1400_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, QCDpt1400, 28, 1, "QCDpt1400", true, 0. };
     return ret;
}

// photonJetPt20to30 sample
Sample fphotonJetPt20to30 ()
{
     TChain *c = makeChain((prefix + "PhotonJet_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, photonJetPt20to30, 26, 1, "photonJetPt20to30", true, 0. };
     return ret;
}

// photonJetPt20to30 sample
Sample fphotonJetPt30to50 ()
{
     TChain *c = makeChain((prefix + "PhotonJet_Pt30to50_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, photonJetPt30to50, 26, 1, "photonJetPt30to50", true, 0. };
     return ret;
}

// photonJetPt20to30 sample
Sample fphotonJetPt50to80 ()
{
     TChain *c = makeChain((prefix + "PhotonJet_Pt50to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, photonJetPt50to80, 26, 1, "photonJetPt50to80", true, 0. };
     return ret;
}

// photonJetPt20to30 sample
Sample fphotonJetPt80to120 ()
{
     TChain *c = makeChain((prefix + "PhotonJet_Pt80to120_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, photonJetPt80to120, 26, 1, "photonJetPt80to120", true, 0. };
     return ret;
}

// photonJetPt20to30 sample
Sample fphotonJetPt120to170 ()
{
     TChain *c = makeChain((prefix + "PhotonJet_Pt120to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, photonJetPt120to170, 26, 1, "photonJetPt120to170", true, 0. };
     return ret;
}

// photonJetPt20to30 sample
Sample fphotonJetPt170to300 ()
{
     TChain *c = makeChain((prefix + "PhotonJet_Pt170to300_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, photonJetPt170to300, 26, 1, "photonJetPt170to300", true, 0. };
     return ret;
}

// combination of all photonJet samples
Sample fphotonJetPt20to300 ()
{    
     TChain *c = makeChain((prefix + "PhotonJet_Pt{20to30,30to50,50to80,80to120,120to170,170to300}_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, photonJetPt20to300, 26, 1, "photonJetPt20to300", true, 0. };
     return ret;
}    

// single electron sample
Sample fsingleElectron ()
{
     TChain *c = makeChain((prefix + "SingleElectronPt5to100_336patch4MC31XV9/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, singleElectron, kBlue, 1, "singleElectron", true, 0. };
     return ret;
}

// single photon sample
Sample fsingleGamma ()
{
     TChain *c = makeChain((prefix + "SingleGammaPt10to100_31X/V03-00-35/merged_ntuple*.root").c_str());
     Sample ret = { c, singleGamma, kGreen, 1, "singleGamma", true, 0. };
     return ret;
}

