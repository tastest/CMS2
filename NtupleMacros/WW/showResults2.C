#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "wwtypes.cc"
#include <assert.h>
#include <math.h>
#include "TCut.h"
#include "TDirectory.h"
#include "TROOT.h"

TH1F* getYields(const char* file, const char* cut = 0, const char* tree = "tree")
{
  TH1F* output(0);
  TFile *f = dynamic_cast<TFile*>(gROOT->GetListOfFiles()->FindObject(file));
  if (!f ) f = TFile::Open(file);
  if (!f) return output;
  TTree* tt = dynamic_cast<TTree*>(f->Get(tree));
  if (!tt) return output;
  if (cut)
    tt->Draw("type>>h(4,0,4)",Form("scale1fb*(%s)",cut),"e goff");
  else
    tt->Draw("type>>h(4,0,4)","scale1fb","e goff");
  output = dynamic_cast<TH1F*>(gDirectory->Get("h"));
  assert(output);
  output->SetDirectory(0);
  return output;
}

double getValue(const TH1F* hist, unsigned int type){
  switch (type){
  case 0:
    return hist->GetBinContent(1);
    break;
  case 1:
    return hist->GetBinContent(2);
    break;
  case 2:
    return hist->GetBinContent(3);
    break;
  case 3:
    return hist->GetBinContent(4);
    break;
  default:
    return hist->Integral(1,4);
  }
}
 
double getError(const TH1F* hist, unsigned int type){
  switch (type){
  case 0:
    return hist->GetBinError(1);
    break;
  case 1:
    return hist->GetBinError(2);
    break;
  case 2:
    return hist->GetBinError(3);
    break;
  case 3:
    return hist->GetBinError(4);
    break;
  default:
    return sqrt(pow(hist->GetBinError(1),2) + pow(hist->GetBinError(2),2)+
		pow(hist->GetBinError(3),2) + pow(hist->GetBinError(4),2));
  }
}

void showResults2(const char* cut = 0, bool showTotalBkg=true)
{
  using namespace std;
  
  std::vector<std::pair<TH1F*,std::string> > bkgs;
  if ( TH1F *hist  = getYields("smurf/dyee.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*DY ee*"));
  if ( TH1F *hist  = getYields("smurf/dymm.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*DY mumu*"));
  if ( TH1F *hist  = getYields("smurf/dytt.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*DY tautau*"));
  if ( TH1F *hist  = getYields("smurf/ttbar.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*ttbar*"));
  if ( TH1F *hist  = getYields("smurf/tw.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*TW*"));
  if ( TH1F *hist  = getYields("smurf/wjets.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*Wjets*"));
  if ( TH1F *hist  = getYields("smurf/wz.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*WZ*"));
  if ( TH1F *hist  = getYields("smurf/zz.root",cut) )
    bkgs.push_back(std::pair<TH1F*,std::string>(hist,"*ZZ*"));
  TH1F *ww     = getYields("smurf/ww.root",cut);
  TH1F *hww130 = getYields("smurf/hww130.root",cut);
  TH1F *hww160 = getYields("smurf/hww160.root",cut);
  TH1F *hww200 = getYields("smurf/hww200.root",cut);
  TH1F *hww250 = getYields("smurf/hww250.root",cut);
  TH1F *data   = getYields("smurf/data.root",cut);
  if (!hww130 && !hww160 && !hww200 && !hww250 && 
      !ww && !data && bkgs.empty()){
    cout << "no data is found." << endl;
    return;
  }
  const char* patternTitle = " %12s |";
  const char* patternData  = " %5.2f%s%4.2f |";

  cout << "\n" << Form("| %3s |","");
  for (unsigned int i=0; i<bkgs.size(); ++i) cout << Form(patternTitle,bkgs.at(i).second.c_str());
  if ( showTotalBkg && bkgs.size()>0 ) cout << Form(patternTitle,"*Total BKG*");
  cout << endl;
  // string pm = "+/-";
  string pm = "&plusmn;";

  double bkg[5] = {0, 0, 0, 0, 0};
  double bkgerr2[5] = {0, 0, 0, 0, 0};
  const char* names[5] = {"mm","me","em","ee","all"};
  for (int i=0; i<5; i++){
    cout << "|" << Form(" %3s ",names[i]) << "|";
    for (unsigned int j=0; j<bkgs.size(); ++j){
      TH1F* hist = bkgs.at(j).first;
      cout << Form(patternData,getValue(hist,i),pm.c_str(),getError(hist,i));
      bkg[i]     += getValue(hist,i);
      bkgerr2[i] += pow(getError(hist,i),2);
    }
    if ( showTotalBkg && bkgs.size()>0 ) cout << Form(patternData,bkg[i],pm.c_str(),sqrt(bkgerr2[i]));
    cout <<endl;
  }
  cout <<endl;

  cout << "\n" << Form("| %3s |","");
  if (ww) cout << Form(patternTitle,"*WW*");
  if (hww130) cout << Form(patternTitle,"*HWW130*");
  if (hww160) cout << Form(patternTitle,"*HWW160*");
  if (hww200) cout << Form(patternTitle,"*HWW200*");
  if (hww250) cout << Form(patternTitle,"*HWW250*");
  if (data) cout << Form(patternTitle,"*Data*");
  cout << endl;

  for (int i=0; i<5; i++){
    cout << "|" << Form(" %3s ",names[i]) << "|";
    if (ww) 
      cout << Form(patternData,getValue(ww,i),pm.c_str(),getError(ww,i));
    if (hww130) 
      cout << Form(patternData,getValue(hww130,i),pm.c_str(),getError(hww130,i));
    if (hww160) 
      cout << Form(patternData,getValue(hww160,i),pm.c_str(),getError(hww160,i));
    if (hww200) 
      cout << Form(patternData,getValue(hww200,i),pm.c_str(),getError(hww200,i));
    if (hww250) 
      cout << Form(patternData,getValue(hww250,i),pm.c_str(),getError(hww250,i));
    if (data)
      cout << Form(patternData,getValue(data,i),pm.c_str(),getError(data,i));
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
