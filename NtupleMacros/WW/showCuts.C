#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "wwtypes.cc"
#include <assert.h>
#include <math.h>
struct Process{
  TH1* htot;
  TH1* hcuts_ee;
  TH1* hcuts_em;
  TH1* hcuts_mm;
  TH1* hcuts_all;
  const char* name;
  Process():htot(0),hcuts_ee(0),hcuts_em(0),hcuts_mm(0),hcuts_all(0){}
  Process(const char* _name, TFile* ftt, const char* prefix){
    htot      = dynamic_cast<TH1F*>(ftt->Get(Form("%s_hypos_total_weighted",prefix)));
    hcuts_ee  = dynamic_cast<TH1F*>(ftt->Get(Form("%s_hcuts_ee",prefix)));
    hcuts_em  = dynamic_cast<TH1F*>(ftt->Get(Form("%s_hcuts_em",prefix)));
    hcuts_mm  = dynamic_cast<TH1F*>(ftt->Get(Form("%s_hcuts_mm",prefix)));
    hcuts_all = dynamic_cast<TH1F*>(ftt->Get(Form("%s_hcuts_all",prefix)));
    name = _name;
  }
  bool isGood() const { return htot && hcuts_ee && hcuts_em && hcuts_mm && hcuts_all; }
};
void showCuts(const char* file = "processed_data.root")
{
  using namespace std;
  
  //hist file:
  TFile *ftt = TFile::Open(file);
  assert(ftt);

  std::vector<Process> mc;
  mc.push_back(Process("*DY ee*",     ftt, "dyee"));
  mc.push_back(Process("*DY mumu*",   ftt, "dymm"));
  mc.push_back(Process("*DY tautau*", ftt, "dytt"));
  mc.push_back(Process("*ttbar*",     ftt, "ttbar"));
  mc.push_back(Process("*tW*",        ftt, "tw"));
  mc.push_back(Process("*Wjets*",     ftt, "wjets"));
  mc.push_back(Process("*WZ*",        ftt, "wz"));
  mc.push_back(Process("*ZZ*",        ftt, "zz"));
  mc.push_back(Process("*WW*",        ftt, "ww"));
  mc.push_back(Process("*qqWW*",      ftt, "qqww"));
  // mc.push_back(Process("*QCD*",       ftt, "qcd"));
  
  std::vector<Process> data;
  data.push_back(Process("*Data*",        ftt, "data"));
  int nCuts(0);
  printf("MC samples: ");
  for(std::vector<Process>::const_iterator p=mc.begin(); p!=mc.end(); ++p){
    if (! p->isGood() ) continue;
    if (nCuts==0) nCuts=p->hcuts_ee->GetNbinsX();
    nCuts = std::max(nCuts,p->hcuts_ee->GetNbinsX());
    nCuts = std::max(nCuts,p->hcuts_em->GetNbinsX());
    nCuts = std::max(nCuts,p->hcuts_mm->GetNbinsX());
    nCuts = std::max(nCuts,p->hcuts_all->GetNbinsX());
    printf("%s ",p->name);
  }
  printf("\nData samples: ");
  for(std::vector<Process>::const_iterator p=data.begin(); p!=data.end(); ++p){
    if (! p->isGood() ) continue;
    if (nCuts==0) nCuts=p->hcuts_ee->GetNbinsX();
    nCuts = std::max(nCuts,p->hcuts_ee->GetNbinsX());
    nCuts = std::max(nCuts,p->hcuts_em->GetNbinsX());
    nCuts = std::max(nCuts,p->hcuts_mm->GetNbinsX());
    nCuts = std::max(nCuts,p->hcuts_all->GetNbinsX());
    printf("%s ",p->name);
  }
  printf("\nNumber of cuts: %d\n",nCuts);
  printf("Selection & MC ee & MC em & MC mm & MC Total & Data ee & Data em & Data mmu & Data Total\n");
  for (int i=0; i<nCuts; ++i){
    double mc_ee(0), mc_em(0), mc_mm(0), mc_all(0);
    double data_ee(0), data_em(0), data_mm(0), data_all(0);
    const char* name = "";
    for(std::vector<Process>::const_iterator p=mc.begin(); p!=mc.end(); ++p){
      if (! p->isGood() ) continue;
      double scale = p->hcuts_all->GetBinContent(nCuts)>0?p->htot->GetBinContent(4)/p->hcuts_all->GetBinContent(nCuts):0;
      name = p->hcuts_ee->GetXaxis()->GetBinLabel(i+1);
      mc_ee  += scale*p->hcuts_ee->GetBinContent(i+1);
      mc_em  += scale*p->hcuts_em->GetBinContent(i+1);
      mc_mm  += scale*p->hcuts_mm->GetBinContent(i+1);
      mc_all += scale*p->hcuts_all->GetBinContent(i+1);
    }
    for(std::vector<Process>::const_iterator p=data.begin(); p!=data.end(); ++p){
      if (! p->isGood() ) continue;
      data_ee  += p->hcuts_ee->GetBinContent(i+1);
      data_em  += p->hcuts_em->GetBinContent(i+1);
      data_mm  += p->hcuts_mm->GetBinContent(i+1);
      data_all += p->hcuts_all->GetBinContent(i+1);
    }
    printf("%2d & %9.1f & %9.1f & %9.1f & %9.1f & %9.1f & %9.1f & %9.1f & %9.1f \t%s\n",i+1,
	   mc_ee,mc_em,mc_mm,mc_all,
	   data_ee,data_em,data_mm,data_all,
	   name);
  } 
}
/*  if (!ww && !data && bkgs.empty()){
    cout << "no data is found." << endl;
    return;
  }
  const char* patternTitle = " %12s |";
  const char* patternData  = " %5.2f%s%4.2f |";

  cout << "\n" << Form("| %3s |","");
  for (unsigned int i=0; i<bkgs.size(); ++i) cout << Form(patternTitle,bkgs.at(i).second.c_str());
  if ( bkgs.size()>0 ) cout << Form(patternTitle,"*Total BKG*");
  if (ww) cout << Form(patternTitle,"*WW*");
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
    if (data)
      cout << Form(patternData,data->GetBinContent(i+1),pm.c_str(),data->GetBinError(i+1));
    cout <<endl;
  }
  cout <<endl;


    */
