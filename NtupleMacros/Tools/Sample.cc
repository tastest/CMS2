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
     ret.chain->Add(b.chain);
//      printf("adding %llu (%d files), returned chain now has %llu entries (%d files)\n", 
// 	    b.chain->GetEntries(), b.chain->GetListOfFiles()->GetEntries(),
// 	    a.chain->GetEntries(), a.chain->GetListOfFiles()->GetEntries());
     return ret;
}

static const std::string prefix = (getenv("CMS2_NTUPLE_LOCATION") != 0) ? 
     std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" : "/data/tmp/";

// New test file
Sample fTest ()
{
     TChain *c = new TChain("Events");
     std::string sample = "/home/dlevans/eleID/CMSSW_2_2_10/src/CMS2/NtupleMaker/test/postprocessed_ntuple.root";
     c->Add(sample.c_str());
     Sample ret = { c, TEST, kRed, 1, "test", true, 0. };
     return ret;
}

     
//WW file
Sample fWW ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/WW_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, WW, kRed, 1, "ww", true, 0. };
     return ret;
}

//WW file
/*
Sample fWW_excl ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/WW_2l_Summer08_IDEAL_V9_v2/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, WW, kRed, 1, "ww", true, 0. };
     return ret;
}
*/
//WZ file
Sample fWZ ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/WZ_incl_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, WZ, kBlue, 1, "wz", true, 0. };
     return ret;
}

//ZZ file
Sample fZZ ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/ZZ_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, ZZ, kGreen, 1, "zz", true, 0. };
     return ret;
}

//Wjets file
Sample fWjets ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/WJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, Wjets, 40, 1, "wjets", true, 0. };
     return ret;
}

//WjetsAlpgen file
Sample fWjetsAlpgenSingle ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_1jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_1jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_1jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_1jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_1jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_2jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_2jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_2jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_2jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_2jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_3jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_3jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_3jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_3jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_3jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_4jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_4jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_4jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_4jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_4jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_5jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_5jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_5jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_5jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/W_5jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());

     Sample ret = { c, Wjets, 40, 1, "wjetsAlpgen", true, 0. };
     return ret;
}


Sample fZeejetsAlpgenSingle ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/Z_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());

     Sample ret = { c, DYee, 42, 1, "dyeeAlpgen", true, 0. };
     return ret;
}

Sample fZmmjetsAlpgenSingle ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/Z_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());

     Sample ret = { c, DYmm, 44, 1, "dymmAlpgen", true, 0. };
     return ret;
}

Sample fZttjetsAlpgenSingle ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/Z_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_1jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_2jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_3jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_4jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());
     sample = prefix + "cms2-V01-03-01/Z_5jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merge*root";
     c->Add(sample.c_str());

     Sample ret = { c, DYtt, 46, 1, "dyttAlpgen", true, 0. };
     return ret;
}



//Wjets file (single lepton filter)
/*
Sample fWjetsSingle ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/WJets-madgraph_Fall08_IDEAL_V11_redigi_v1_Single-lepton/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, Wjets, 40, 1, "wjets", true, 0. };
     return ret;
}
*/
//Wjets file
Sample fWc ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/Wc-madgraph_Fall08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, Wjets, 40, 1, "wc", true };
     return ret;
}

// "vlqq" sample
Sample fVlqq ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/VQQ-madgraph_Summer08_IDEAL_V11_redigi_v2/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DY, kBlack, 1, "vlqq", true };
     return ret;
}

//DYee file
Sample fDYee ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kMagenta, 1, "dyee", true, 0. };
     return ret;
}

//DYmm file
Sample fDYmm ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kCyan, 1, "dymm", true, 0. };
     return ret;
}

//DYtt file
Sample fDYtt ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYtt, kBlack, 1, "dytt", true, 0. };
     return ret;
}

// Pythia DY ntuples
//
// DYee
Sample fDY20ee    ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/Zee_M20_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kMagenta, 1, "dy20ee", true, 0. };
     return ret;
}
// DYmm
Sample fDY20mm    ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/Zmumu_M20_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kCyan, 1, "dy20mm", true, 0. };
     return ret;
}
// DYtt
Sample fDY20tt    ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/Ztautau_M20_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYtt, kBlack, 1, "dy20tt", true, 0. };
     return ret;
}
//
//

// low-mass DY sample
Sample fAstar ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/AstarJets-madgraph_Fall08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DY, kBlack, 1, "astar", true };
     return ret;
}
/*
// Wgamma
Sample fWgamma ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/Wgamma_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, Wgamma, kBlack, 1, "wgamma", true };
     return ret;
}

// Zgamma
Sample fZgamma ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/Zgamma_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, Zgamma, kBlack, 1, "zgamma", true };
     return ret;
}
*/
//ttbar file
Sample fttbar ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, ttbar, kYellow, 1, "ttbar", true, 0. };
     return ret;
}

//ttbar file
Sample fttbarSingle ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10-SingleLepton/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, ttbar, kYellow, 1, "ttbar", true, 0. };
     return ret;
}

//ttbar file
Sample fttbar_taula ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/TauolaTTbar_Summer08_IDEAL_V11_redigi_v2/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, ttbar, kYellow, 1, "ttbartauola", true, 0. };
     return ret;
}

Sample ftW ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "/cms2-V01-03-01/SingleTop_tWChannel_Summer08_IDEAL_V11_redigi_v3/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, tW, 63, 1, "tw", true, 0. };
     return ret;
}

Sample fSingleTop_tChannel ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SingleTop_tChannel_Summer08_IDEAL_V11_redigi_v3/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, tW, 63, 1, "singletopt", true, 0. };
     return ret;
}

Sample fSingleTop_sChannel ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SingleTop_sChannel_Summer08_IDEAL_V11_redigi_v3/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, tW, 63, 1, "singletops", true, 0. };
     return ret;
}

Sample fLM0 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM0-sftsht_Summer08_IDEAL_V11_v1/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM0, 37, 1, "LM0", false, 0. };
     return ret;
}

Sample fLM1 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM1-sftsht_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM1, 37, 1, "LM1", false, 0. };
     return ret;
}

Sample fLM2 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM2-sftsht_Summer08_IDEAL_V11_redigi_v1/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM2, 37, 1, "LM2", false, 0. };
     return ret;
}

Sample fLM3 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM3-sftsht_Summer08_IDEAL_V11_redigi_v1/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM3, 37, 1, "LM3", false, 0. };
     return ret;
}

Sample fLM4 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM4-sftsht_Summer08_IDEAL_V11_redigi_v1/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM4, 37, 1, "LM4", false, 0. };
     return ret;
}

Sample fLM5 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM5-sftsht_Summer08_IDEAL_V11_redigi_v1/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM5, 37, 1, "LM5", false, 0. };
     return ret;
}

Sample fLM6 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM6-sftsht_Summer08_IDEAL_V11_redigi_v1/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM6, 37, 1, "LM6", false, 0. };
     return ret;
}

Sample fLM7 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM7-sftsht_Summer08_IDEAL_V11_redigi_v1/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM7, 37, 1, "LM7", false, 0. };
     return ret;
}

Sample fLM8 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM8-sftsht_Summer08_IDEAL_V11_redigi_v1/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM8, 37, 1, "LM8", false, 0. };
     return ret;
}

Sample fLM9 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM9-sftsht_Summer08_IDEAL_V11_redigi_v1/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM9, 37, 1, "LM9", false, 0. };
     return ret;
}

Sample fLM9p ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM9p-sftsht_Summer08_IDEAL_V11_redigi_v1/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM9p, 37, 1, "LM9p", false, 0. };
     return ret;
}

Sample fLM10 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM10-sftsht_Summer08_IDEAL_V11_redigi_v1/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM10, 37, 1, "LM10", false, 0. };
     return ret;
}

Sample fLM11 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/SUSY_LM11-sftsht_Summer08_IDEAL_V11_redigi_v1/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM11, 37, 1, "LM11", false, 0. };
     return ret;
}

Sample fQCDpt30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCDpt30_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt30, 28, 1, "QCDpt30", true, 999999999 };
     return ret;
}

Sample fQCDpt80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCDpt80_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt80, 28, 1, "QCDpt80", true, 0. };
     return ret;
}

Sample fInclusiveMuPt15Single ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/InclusiveMuPt15_Summer08_IDEAL_V11_redigi_v1-SingleLepton/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, InclusiveMuPt15, 28, 1, "InclusiveMuPt15", true, 0. };
     return ret;
}


/*
// QCD samples
Sample fInclusiveMu5Pt50 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/InclusiveMu5Pt50/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, InclusiveMu5Pt50, 28, 1, "InclusiveMu5Pt50", true, 0. };
     return ret;
}

Sample fInclusiveMuPt15 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/InclusiveMuPt15/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, InclusiveMuPt15, 28, 1, "InclusiveMuPt15", true, 0. };
     return ret;
}

Sample fQCDBCtoEPt20to30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCD_BCtoE_Pt20to30/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt20to30, 28, 1, "QCDBCtoEPt20to30", true, 0. };
     return ret;
}

Sample fQCDBCtoEPt30to80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCD_BCtoE_Pt30to80/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt30to80, 28, 1, "QCDBCtoEPt30to80", true, 0. };
     return ret;
}

Sample fQCDBCtoEPt80to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCD_BCtoE_Pt80to170/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt80to170, 28, 1, "QCDBCtoEPt80to170", true, 0. };
     return ret;
}

Sample fQCDEMenrichedPt20to30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCD_EMenriched_Pt20to30/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt20to30, 28, 1, "QCDEMenrichedPt20to30", true, 0. };
     return ret;
}

Sample fQCDEMenrichedPt30to80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCD_EMenriched_Pt30to80/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt30to80, 28, 1, "QCDEMenrichedPt30to80", true, 0. };
     return ret;
}

Sample fQCDEMenrichedPt80to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCD_EMenriched_Pt80to170/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt80to170, 28, 1, "QCDEMenrichedPt80to170", true, 0. };
     return ret;
}


Sample fQCDpt30to80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCDpt30_v2/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt30to80, 28, 1, "QCDpt30to80", true, 80. };
     return ret;
}

Sample fQCDpt80to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCDpt80/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt80to170, 28, 1, "QCDpt80to170", true, 170 };
     return ret;
}

Sample fQCDpt170to300 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCDpt170/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt170to300, 28, 1, "QCDpt170to300", true, 300 };
     return ret;
}
Sample fQCDpt300to470 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCDpt300/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt300to470, 28, 1, "QCDpt300to470", true, 470 };
     return ret;
}
Sample fQCDpt470to800 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCDpt470/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt470to800, 28, 1, "QCDpt470to800", true, 800 };
     return ret;
}
Sample fQCDpt800toInf ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-03-01/QCDpt800/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt800toInf, 28, 1, "QCDpt800toInf", true, 999999999 };
     return ret;
}
*/
