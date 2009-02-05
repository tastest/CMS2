#include "TFile.h"
#include "TH1F.h"
#include <iostream>
void doTable()
{
  using namespace std;
  
  //hist file:
  TFile *ftt = TFile::Open("processed_data_tag.root");
  assert(ftt);

  TH1F *DYee[4], *DYmm[4], *DYtt[4], *tt[4], *wjets[4], *wz[4], *zz[4], *ww[4], *tw[4];

  for ( int njet=1; njet<4; ++njet) {// actual njet is one less than this counter!

    //do the emu row:
    DYee[0]  = dynamic_cast<TH1F*>(ftt->Get("dyee_hnJet_em"));
    DYmm[0]  = dynamic_cast<TH1F*>(ftt->Get("dymm_hnJet_em"));
    DYtt[0]  = dynamic_cast<TH1F*>(ftt->Get("dytt_hnJet_em"));
    tt[0]    = dynamic_cast<TH1F*>(ftt->Get("ttbar_hnJet_em"));
    wjets[0] = dynamic_cast<TH1F*>(ftt->Get("wjets_hnJet_em"));
    wz[0]    = dynamic_cast<TH1F*>(ftt->Get("wz_hnJet_em"));
    zz[0]    = dynamic_cast<TH1F*>(ftt->Get("zz_hnJet_em"));
    ww[0]    = dynamic_cast<TH1F*>(ftt->Get("ww_hnJet_em"));
    tw[0]    = dynamic_cast<TH1F*>(ftt->Get("tw_hnJet_em"));
  
    //do the mumu row:
    DYee[1]  = dynamic_cast<TH1F*>(ftt->Get("dyee_hnJet_mm"));
    DYmm[1]  = dynamic_cast<TH1F*>(ftt->Get("dymm_hnJet_mm"));
    DYtt[1]  = dynamic_cast<TH1F*>(ftt->Get("dytt_hnJet_mm"));
    tt[1]    = dynamic_cast<TH1F*>(ftt->Get("ttbar_hnJet_mm"));
    wjets[1] = dynamic_cast<TH1F*>(ftt->Get("wjets_hnJet_mm"));
    wz[1]    = dynamic_cast<TH1F*>(ftt->Get("wz_hnJet_mm"));
    zz[1]    = dynamic_cast<TH1F*>(ftt->Get("zz_hnJet_mm"));
    ww[1]    = dynamic_cast<TH1F*>(ftt->Get("ww_hnJet_mm"));
    tw[1]    = dynamic_cast<TH1F*>(ftt->Get("tw_hnJet_mm"));
  
    //do the ee row:
    DYee[2]  = dynamic_cast<TH1F*>(ftt->Get("dyee_hnJet_ee"));
    DYmm[2]  = dynamic_cast<TH1F*>(ftt->Get("dymm_hnJet_ee"));
    DYtt[2]  = dynamic_cast<TH1F*>(ftt->Get("dytt_hnJet_ee"));
    tt[2]    = dynamic_cast<TH1F*>(ftt->Get("ttbar_hnJet_ee"));
    wjets[2] = dynamic_cast<TH1F*>(ftt->Get("wjets_hnJet_ee"));
    wz[2]    = dynamic_cast<TH1F*>(ftt->Get("wz_hnJet_ee"));
    zz[2]    = dynamic_cast<TH1F*>(ftt->Get("zz_hnJet_ee"));
    ww[2]    = dynamic_cast<TH1F*>(ftt->Get("ww_hnJet_ee"));
    tw[2]    = dynamic_cast<TH1F*>(ftt->Get("tw_hnJet_ee"));
  
    //do the total row:
    if (DYee[0]) 
      DYee[3] = new TH1F(*(DYee[0]));
    else
      DYee[3] = 0;
    if (DYmm[0])
      DYmm[3] = new TH1F(*(DYmm[0]));
    else
      DYmm[3] = 0;
    if (DYtt[0])
      DYtt[3] = new TH1F(*(DYtt[0]));
    else
      DYtt[3] = 0;
    if (tt[0]) 
      tt[3] = new TH1F(*(tt[0]));
    else
      tt[3] = 0;
    if (wjets[0])
      wjets[3] = new TH1F(*(wjets[0]));
    else
      wjets[3] = 0;
    if (wz[0])
      wz[3] = new TH1F(*(wz[0]));
    else 
      wz[0] = 0;
    if (zz[0])
      zz[3] = new TH1F(*(zz[0]));
    else
      zz[3] = 0;
    if (ww[0])
      ww[3] = new TH1F(*(ww[0]));
    else
      ww[3] = 0;
    if (tw[0])
      tw[3] = new TH1F(*(tw[0]));
    else
      tw[3] = 0;
  
    char* finalState[4];
    finalState[0] = "emu";
    finalState[1] = "mumu";
    finalState[2] = "ee";
    finalState[3] = "total";
  
    cout << "  " << endl;
    cout << " Njet = " << njet-1 << endl;
    cout << "| |  *DY ee*  |  *DY mumu*  |  *DY tautau*  |  *ttbar*  |  *Wjets*  |  *WZ*  |  *ZZ*  |  *WW*  |  *TW*  | " << endl;
    
    for (int i=0; i<4; i++){
      
      cout << "| " << finalState[i] << " |  ";
      cout.setf(ios::fixed,ios::floatfield);
      cout.precision(1);
      if (DYee[i]) cout << DYee[i]->GetBinContent(njet) << " &plusmn; " << DYee[i]->GetBinError(njet);
      cout << "  |  ";
      if (DYmm[i]) cout << DYmm[i]->GetBinContent(njet) << " &plusmn; " << DYmm[i]->GetBinError(njet);
      cout << "  |  ";
      if (DYtt[i]) cout << DYtt[i]->GetBinContent(njet) << " &plusmn; " << DYtt[i]->GetBinError(njet);
      cout << "  |  ";
      if (tt[i]) cout << tt[i]->GetBinContent(njet) << " &plusmn; " << tt[i]->GetBinError(njet);
      cout << "  |  ";
      if (wjets[i]) cout << wjets[i]->GetBinContent(njet) << " &plusmn; " << wjets[i]->GetBinError(njet);
      cout << "  |  ";
      if (wz[i]) cout << wz[i]->GetBinContent(njet) << " &plusmn; " << wz[i]->GetBinError(njet);
      cout << "  |  ";
      if (zz[i]) cout << zz[i]->GetBinContent(njet) << " &plusmn; " << zz[i]->GetBinError(njet);
      cout << "  |  ";
      if (ww[i]) cout << ww[i]->GetBinContent(njet) << " &plusmn; " << ww[i]->GetBinError(njet);
      cout << "  |  ";
      if (tw[i]) cout << tw[i]->GetBinContent(njet) << " &plusmn; " << tw[i]->GetBinError(njet);
      cout << "  |  " <<endl;
      
      if ( i>0 && i<3 ) {
	if (DYee[3] && DYee[i]) DYee[3]->Add(DYee[i]);
	if (DYmm[3] && DYmm[i]) DYmm[3]->Add(DYmm[i]);
	if (DYtt[3] && DYtt[i]) DYtt[3]->Add(DYtt[i]);
	if (tt[3]   && tt[i])   tt[3]->Add(tt[i]);
	if (wjets[3]&&wjets[i]) wjets[3]->Add(wjets[i]);
	if (wz[3]   && wz[i])   wz[3]->Add(wz[i]);
	if (zz[3]   && zz[i])   zz[3]->Add(zz[i]);
	if (ww[3]   && ww[i])   ww[3]->Add(ww[i]);
	if (tw[3]   && tw[i])   tw[3]->Add(tw[i]);
    
      }
    }//end of for loop over njet
  }
}
