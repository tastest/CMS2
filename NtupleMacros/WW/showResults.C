#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <iomanip>
#include "TString.h"
void showResults()
{
  using namespace std;
  
  //hist file:
  TFile *ftt = TFile::Open("processed_data_tag.root");
  assert(ftt);

    //do the emu row:
    TH1F *DYee  = dynamic_cast<TH1F*>(ftt->Get("dyee_hypos_total_weighted"));
    TH1F *DYmm  = dynamic_cast<TH1F*>(ftt->Get("dymm_hypos_total_weighted"));
    TH1F *DYtt  = dynamic_cast<TH1F*>(ftt->Get("dytt_hypos_total_weighted"));
    TH1F *tt    = dynamic_cast<TH1F*>(ftt->Get("ttbar_hypos_total_weighted"));
    TH1F *wjets = dynamic_cast<TH1F*>(ftt->Get("wjets_hypos_total_weighted"));
    TH1F *wz    = dynamic_cast<TH1F*>(ftt->Get("wz_hypos_total_weighted"));
    TH1F *zz    = dynamic_cast<TH1F*>(ftt->Get("zz_hypos_total_weighted"));
    TH1F *ww    = dynamic_cast<TH1F*>(ftt->Get("ww_hypos_total_weighted"));
    TH1F *tw    = dynamic_cast<TH1F*>(ftt->Get("tw_hypos_total_weighted"));
  
    char* finalState[4];
    finalState[0] = " ee   ";
    finalState[1] = " mumu ";
    finalState[2] = " em   ";
    finalState[3] = " total";
    
    cout << "|      |  *DY ee*  | *DY mumu* |*DY tautau*|  *ttbar*  |  *Wjets*  |    *WZ*   |    *ZZ*   |    *WW*    |    *TW*   |" << endl;
    
    for (int i=0; i<4; i++){
      
      cout << "|" << finalState[i] << "| ";
      cout.setf(ios::fixed,ios::floatfield);
      cout.precision(1);
      if (DYee) cout << DYee->GetBinContent(i+1) << "&plusmn;" << DYee->GetBinError(i+1);
      cout << " | ";
      if (DYmm) cout << DYmm->GetBinContent(i+1) << "&plusmn;" << DYmm->GetBinError(i+1);
      cout << " | ";
      if (DYtt) cout << DYtt->GetBinContent(i+1) << "&plusmn;" << DYtt->GetBinError(i+1);
      cout << " | ";
      if (tt) cout << tt->GetBinContent(i+1) << "&plusmn;" << tt->GetBinError(i+1);
      cout << " | ";
      if (wjets) cout << wjets->GetBinContent(i+1) << "&plusmn;" << wjets->GetBinError(i+1);
      cout << " | ";
      if (wz) cout << wz->GetBinContent(i+1) << "&plusmn;" << wz->GetBinError(i+1);
      cout << " | ";
      if (zz) cout << zz->GetBinContent(i+1) << "&plusmn;" << zz->GetBinError(i+1);
      cout << " | ";
      if (ww) cout << ww->GetBinContent(i+1) << "&plusmn;" << ww->GetBinError(i+1);
      cout << " | ";
      if (tw) cout << tw->GetBinContent(i+1) << "&plusmn;" << tw->GetBinError(i+1);
      cout << " | " <<endl;
    }
}
