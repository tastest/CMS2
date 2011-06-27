#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "wwtypes.cc"
#include <assert.h>
#include <math.h>
void showResults(const char* file = "processed_data.root")
{
  using namespace std;
  
  //hist file:
  TFile *ftt = TFile::Open(file);
  assert(ftt);

  std::vector<std::pair<TH1F*,std::string> > bkgs;
  if ( TH1F *hist  = dynamic_cast<TH1F*>(ftt->Get("dyee_hypos_total_weighted")) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*DY ee*"));
  if ( TH1F *hist  = dynamic_cast<TH1F*>(ftt->Get("dymm_hypos_total_weighted")) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*DY mumu*"));
  if ( TH1F *hist  = dynamic_cast<TH1F*>(ftt->Get("dytt_hypos_total_weighted")) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*DY tautau*"));
  if ( TH1F *hist  = dynamic_cast<TH1F*>(ftt->Get("ttbar_hypos_total_weighted")) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*ttbar*"));
  if ( TH1F *hist  = dynamic_cast<TH1F*>(ftt->Get("tw_hypos_total_weighted")) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*TW*"));
  if ( TH1F *hist = dynamic_cast<TH1F*>(ftt->Get("wjets_hypos_total_weighted")) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*Wjets*"));
  if ( TH1F *hist    = dynamic_cast<TH1F*>(ftt->Get("wz_hypos_total_weighted")) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*WZ*"));
  if ( TH1F *hist    = dynamic_cast<TH1F*>(ftt->Get("zz_hypos_total_weighted")) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*ZZ*"));
  TH1F *ww    = dynamic_cast<TH1F*>(ftt->Get("ww_hypos_total_weighted"));
  TH1F *hww130 = dynamic_cast<TH1F*>(ftt->Get("hww130_hypos_total_weighted"));
  TH1F *data  = dynamic_cast<TH1F*>(ftt->Get("data_hypos_total_weighted"));
  if (!hww130 && !ww && !data && bkgs.empty()){
    cout << "no data is found." << endl;
    return;
  }
  const char* patternTitle = " %16s |";
  const char* patternData  = " %7.2f%s%6.2f |";

  cout << "\n" << Form("| %3s |","");
  for (unsigned int i=0; i<bkgs.size(); ++i) cout << Form(patternTitle,bkgs.at(i).second.c_str());
  if ( bkgs.size()>0 ) cout << Form(patternTitle,"*Total BKG*");
  if (ww) cout << Form(patternTitle,"*WW*");
  if (hww130) cout << Form(patternTitle,"*HWW130*");
  if (data) cout << Form(patternTitle,"*Data*");
  cout << endl;
  string pm = "+/-";
  // string pm = "&plusmn;";

  double bkg[4] = {0, 0, 0, 0};
  double bkgerr2[4] = {0, 0, 0, 0};
  for (int i=0; i<4; i++){
    cout << "|" << Form(" %3s ",HypothesisTypeName(i)) << "|";
    for (unsigned int j=0; j<bkgs.size(); ++j){
      TH1F* hist = bkgs.at(j).first;
      cout << Form(patternData,hist->GetBinContent(i+1),pm.c_str(),hist->GetBinError(i+1));
      bkg[i]     += hist->GetBinContent(i+1);
      bkgerr2[i] += pow(hist->GetBinError(i+1),2);
    }
    if ( bkgs.size()>0 ) cout << Form(patternData,bkg[i],pm.c_str(),sqrt(bkgerr2[i]));
    if (ww) 
      cout << Form(patternData,ww->GetBinContent(i+1),pm.c_str(),ww->GetBinError(i+1));
    if (hww130) 
      cout << Form(patternData,hww130->GetBinContent(i+1),pm.c_str(),hww130->GetBinError(i+1));
    if (data)
      cout << Form(patternData,data->GetBinContent(i+1),pm.c_str(),data->GetBinError(i+1));
    cout <<endl;
  }
  cout <<endl;


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
