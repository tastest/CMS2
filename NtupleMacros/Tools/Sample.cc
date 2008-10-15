#include "TChain.h"
#include "Sample.h"
#include "utilities.h"
#include "selections.h"
#include "CMS2.h"

bool filterByProcess (enum Process sample)
{
     switch (sample) {
     case WW: case WZ: case ZZ: case tW:
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
     c->Add("/data/tmp/cms2-V00-04-00/merge_WW.root");
     Sample ret = { c, WW, kRed, 1, "ww" };
     return ret;
}

//WZ file
Sample fWZ ()
{
     TChain *c = new TChain("Events");
     c->Add("/data/tmp/cms2-V00-05-00/merge_WZ.root");
     Sample ret = { c, WZ, kBlue, 1, "wz" };
     return ret;
}

//ZZ file
Sample fZZ ()
{
     TChain *c = new TChain("Events");
     c->Add("/data/tmp/cms2-V00-05-00/merge_ZZ.root");
     Sample ret = { c, ZZ, kGreen, 1, "zz" };
     return ret;
}

//Wjets file
Sample fWjets ()
{
     TChain *c = new TChain("Events");
     c->Add("/data/tmp/cms2-V00-04-01/merge_Wjet.root");
     Sample ret = { c, Wjets, 40, 1, "wjets" };
     return ret;
}

//DYee file
Sample fDYee ()
{
     TChain *c = new TChain("Events");
     c->Add("/data/tmp/cms2-V00-04-01/merge_DY.root");
     Sample ret = { c, DYee, kMagenta, 1.12, "dyee" };
     return ret;
}

//DYmm file
Sample fDYmm ()
{
     TChain *c = new TChain("Events");
     c->Add("/data/tmp/cms2-V00-04-01/merge_DY.root");
     Sample ret = { c, DYmm, kCyan, 1.12, "dymm" };
     return ret;
}

//DYtt file
Sample fDYtt ()
{
     TChain *c = new TChain("Events");
     c->Add("/data/tmp/cms2-V00-04-01/merge_DY.root");
     Sample ret = { c, DYtt, kBlack, 1.12, "dytt" };
     return ret;
}

//ttbar file
Sample fttbar ()
{
     TChain *c = new TChain("Events");
     c->Add("/data/tmp/cms2-V00-04-01/merge_ttbar.root");
     Sample ret = { c, ttbar, kYellow, 1.85, "ttbar" };
     return ret;
}

Sample ftW ()
{
     TChain *c = new TChain("Events");
     c->Add("/data/tmp/cms2-V00-04-00/merge_tW.root");
     Sample ret = { c, tW, 63, 1, "tw" };
     return ret;
}
