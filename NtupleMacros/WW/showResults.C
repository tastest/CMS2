#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <iomanip>
#include <string>
void showResults(const char* file = "processed_data_tag.root")
{
  using namespace std;
  
  //hist file:
  TFile *ftt = TFile::Open(file);
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
    string pm = "+/-";
    // string pm = "&plusmn;";
    for (int i=0; i<4; i++){
      
      cout << "|" << finalState[i] << "| ";
      cout.setf(ios::fixed,ios::floatfield);
      cout.precision(1);
      if (DYee) cout << DYee->GetBinContent(i+1) << pm << DYee->GetBinError(i+1);
      cout << " | ";
      if (DYmm) cout << DYmm->GetBinContent(i+1) << pm << DYmm->GetBinError(i+1);
      cout << " | ";
      if (DYtt) cout << DYtt->GetBinContent(i+1) << pm << DYtt->GetBinError(i+1);
      cout << " | ";
      if (tt) cout << tt->GetBinContent(i+1) << pm << tt->GetBinError(i+1);
      cout << " | ";
      if (wjets) cout << wjets->GetBinContent(i+1) << pm << wjets->GetBinError(i+1);
      cout << " | ";
      if (wz) cout << wz->GetBinContent(i+1) << pm << wz->GetBinError(i+1);
      cout << " | ";
      if (zz) cout << zz->GetBinContent(i+1) << pm << zz->GetBinError(i+1);
      cout << " | ";
      if (ww) cout << ww->GetBinContent(i+1) << pm << ww->GetBinError(i+1);
      cout << " | ";
      if (tw) cout << tw->GetBinContent(i+1) << pm << tw->GetBinError(i+1);
      cout << " | " <<endl;
    }
    
    string samples[9];
    samples[0] = "dyee";
    samples[1] = "dymm";
    samples[2] = "dytt";
    samples[3] = "ttbar";
    samples[4] = "wjets";
    samples[5] = "wz";
    samples[6] = "zz";
    samples[7] = "ww";
    samples[8] = "tw";
    
    string types[4];
    types[0] = "ee";
    types[1] = "mm";
    types[2] = "em";
    types[3] = "all";

    // top background estimate
    /*
    for ( unsigned int type = 0; type < 4; ++type ){
      double ttbar_yield_before(0.);
      double ttbar_yield_after(0.);
      double tw_yield_before(0.);
      double tw_yield_after(0.);
      double ttbar_yield_before_noveto(0.);
      double ttbar_yield_after_noveto(0.);
      double tw_yield_before_noveto(0.);
      double tw_yield_after_noveto(0.);
      double total_yield_before(0.);
      double total_yield_after(0.);

      cout << "Type: " << types[type] << endl;
      for ( unsigned int sample = 0; sample < 9; ++sample ){
	if ( TH2F* h = dynamic_cast<TH2F*>(ftt->Get( (samples[sample]+"_extramuonsvsnjet_"+types[type]).c_str() )) ){
	  unsigned int mbins = h->GetNbinsX(); // muon bins
	  unsigned int jbins = h->GetNbinsY(); // jet bins
	  if ( samples[sample] == "ttbar" ) {
	    ttbar_yield_before = h->Integral(1,mbins+1,1,1);
	    ttbar_yield_after  = h->Integral(1,1,1,1);
	    ttbar_yield_before_noveto = h->Integral(1,mbins+1,1,jbins+1);
	    ttbar_yield_after_noveto  = h->Integral(1,1,1,jbins+1);
	  }
	  if ( samples[sample] == "tw" ) {
	    tw_yield_before = h->Integral(1,mbins+1,1,1);
	    tw_yield_after  = h->Integral(1,1,1,1);
	    tw_yield_before_noveto = h->Integral(1,mbins+1,1,jbins+1);
	    tw_yield_after_noveto  = h->Integral(1,1,1,jbins+1);
	  }
	  total_yield_before += h->Integral(1,mbins+1,1,1);
	  total_yield_after  += h->Integral(1,1,1,1);
	  cout << "\t" << samples[sample] << " \tbefore: " << h->Integral(1,mbins+1,1,1) 
	       << " \tafter: " << h->Integral(1,1,1,1) <<endl;
	}
      }
      if ( (tw_yield_before+ttbar_yield_before)-(tw_yield_after+ttbar_yield_after) <= 0 ) {
	cout << "\tFailed" << endl;
	continue;
      }
      cout.setf(ios::fixed,ios::floatfield);
      cout.precision(3);
      double top_tag_eff = 1-(tw_yield_after+ttbar_yield_after)/(tw_yield_before+ttbar_yield_before);
      double ttbar_tag_eff = 1-ttbar_yield_after/ttbar_yield_before;
      cout << "\t ttbar tag eff (no jet veto): " << 1-ttbar_yield_after_noveto/ttbar_yield_before_noveto << "\n";
      cout << "\t tw tag eff (no jet veto): " << 1-tw_yield_after_noveto/tw_yield_before_noveto << "\n";
      cout << "\t ttbar tag eff: " << ttbar_tag_eff << "\n";
      cout << "\t tw tag eff: " << 1-tw_yield_after/tw_yield_before << "\n";
      cout << "\t top tag eff: " << top_tag_eff << "\n";
      cout.precision(1);
      cout << "\t tagged yield: " << total_yield_before-total_yield_after << "\n"; 
      cout << "\t top yield: " << (tw_yield_after+ttbar_yield_after)<< "\n"; 
      cout << "\t estimated top yield: " << (total_yield_before-total_yield_after)*(1/top_tag_eff-1) << endl;
      cout << "\t estimated top yield (ttbar eff): " << (total_yield_before-total_yield_after)*(1/ttbar_tag_eff-1) << endl;
    }

    // W+jets background estimation
    cout << "Wjets sideband background estimation" << endl;
    double wjets_el_estatimate(0);
    double wjets_mu_estatimate(0);
    for ( unsigned int sample = 0; sample < 9; ++sample ){
      cout << "\t" << samples[sample];
      TH1F* hel = dynamic_cast<TH1F*>(ftt->Get( (samples[sample]+"_hemElRelIso").c_str() ));
      TH1F* hmu = dynamic_cast<TH1F*>(ftt->Get( (samples[sample]+"_hemMuRelIso").c_str() ));
      if ( hel && hmu ){
	wjets_el_estatimate += hel->Integral(77,84);
	wjets_mu_estatimate += hmu->Integral(77,84);
	cout << "\tel-fake: " << hel->Integral(77,84) << "\tmu-fake: " <<  hmu->Integral(77,84) <<endl;
      } else {
	cout << "no info" << endl;
      }
    }
    cout << "Wjets emu background estimates:" << endl;
    cout << "\ttotal el-fake: " << wjets_el_estatimate << endl;
    cout << "\ttotal mu-fake: " << wjets_mu_estatimate << endl;
    cout << "\tTOTAL: " << wjets_el_estatimate+wjets_mu_estatimate << endl;
    */
}
