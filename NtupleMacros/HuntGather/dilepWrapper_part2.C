{
#include "cuts.h"
#include <vector>
  gROOT->ProcessLine(".L goodrun.cc+");
  gROOT->ProcessLine(".L gather.C+");
  gROOT->ProcessLine(".L gather_several_doAll.C+");
  std::vector<TCut> sels;
  //   sels.push_back(inclusivez_dilep);
  //   sels.push_back(wplusjetsv1_emu);
  //   sels.push_back(wprime100_emu); 
  //   sels.push_back(semileptonictop_emu);
  //   sels.push_back(semileptonictopv2_emu);
  //   sels.push_back(semileptonictopv3_emu);
  //   sels.push_back(highptmuon_emu);
  //   sels.push_back(highptmuonv2_emu);
  //   sels.push_back(multimuon_emu);
  //   sels.push_back(zplusjetsv1_dilep);
  //   sels.push_back(zplusjetsv2_dilep);
  //   sels.push_back(zplusjetsv3_dilep);
  //   sels.push_back(zplusbjetsv1_dilep);
  //   sels.push_back(zplusbjetsv2_dilep);
  
   sels.push_back(zmet30_dilep);
   sels.push_back(zptz50_dilep);
   sels.push_back(zptz100_dilep);
   sels.push_back(inclusivenonz_dilep);
   sels.push_back(inclusivenonzof_dilep);
   sels.push_back(inclusivenonzhighmass_dilep);
   sels.push_back(inclusivenonzhighmassv2_dilep);
   sels.push_back(inclusivenonzjets_dilep);
   sels.push_back(inclusivenonzjetsv2_dilep);
   sels.push_back(inclusivenonzmet20_dilep);
   sels.push_back(inclusivenonzmet20nojets_dilep);
   sels.push_back(inclusivenonzmet40jets_dilep);
   sels.push_back(inclusivenonzptz50_dilep);
   sels.push_back(inclusivenonzptz100_dilep);
  
//   sels.push_back(dileptonictop_dilep);
//   sels.push_back(dileptonictopv2_dilep);
//   sels.push_back(dileptonictopv3_dilep);
//   sels.push_back(dileptonictopv4_dilep);
//   sels.push_back(dileptonictopv4_dilep_mnjets);
//   sels.push_back(samesigninclusive_dilep);
//   sels.push_back(zplusworzprime_dilep);
//   sels.push_back(meff_dilep); 
  
  std::cout<<"Running on "<<sels.size()<<" selections."<<std::endl;
  
  for(unsigned int sel = 0; sel < sels.size(); ++sel) {
    std::cout<<"Processing selection: "<<sels[sel].GetName()<<std::endl;
    gather_several_doAll(sels[sel]);
    std::cout<<"done with "<<sels[sel].GetName()<<std::endl;
  }
}
