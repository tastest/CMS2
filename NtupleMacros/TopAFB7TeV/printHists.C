#include "THStack.h"
#include <iostream>
#include <utility>
#include "TMath.h"
#include "histtools.h"
#include "TDirectory.h"
#include "TROOT.h"
#include <vector>
#include "TString.h"
#include "CommonFunctions.C"
using namespace std;


void printNJets( bool latex=false, const char* formatS = "%6.1f", const char* signalS= "ttprime", 
		 bool combineVVsamples = true, bool combineDYsamples = true, bool combineJetBins = false, bool printProbs=false, bool printErrorsForData = false){


  const char* suffix[4];
  suffix[0] = (const char*)"ee";
  suffix[1] = (const char*)"mm";
  suffix[2] = (const char*)"em";
  suffix[3] = (const char*)"all";


  std::string pmSign  = latex ? " $\\pm$ " : " &plusmn; ";
  std::string colSep  = latex ? " & " : " | ";
  std::string beginL  = latex ? ""   : "| ";
  std::string endL    = latex ? " \\\\ " : " | ";
  std::string mathSep = latex ? "$" : "";


 

  bool haveVVsamples = false;
  if(gROOT->FindObjectAny("ww_hnBtagJet_allj_ee") != NULL || 
     gROOT->FindObjectAny("zz_hnBtagJet_allj_ee") != NULL || 
     gROOT->FindObjectAny("wz_hnBtagJet_allj_ee") != NULL)
    haveVVsamples = true;

  bool haveData = false;
  if(gROOT->FindObjectAny("data_hnBtagJet_allj_ee") != NULL)
    haveData = true;


  
  //if(combineJetBins) { 
  vector<TString> v_prefixes;
  TList *list = gDirectory->GetList();
  for(int i = 0 ; i < list->GetSize(); i++) {
    TString name = list->At(i)->GetName();
    if(!name.EndsWith("all"))
      continue;
    TString prefix = name.Tokenize("_")->At(0)->GetName();
    cout <<prefix <<endl;
    bool old_hist=std::find(v_prefixes.begin(), v_prefixes.end(), prefix) != v_prefixes.end();
    if(!old_hist) {
      v_prefixes.push_back(prefix);
    }
  }
    
  vector<pair<int, int> > v_binranges;
  if(combineJetBins) 
    v_binranges.push_back(make_pair(1,7));
  else {
    v_binranges.push_back(make_pair(1,1)); //0 jet bin 
    v_binranges.push_back(make_pair(2,2)); //1 jet bin 
    v_binranges.push_back(make_pair(3,3)); //2 jet bin
    v_binranges.push_back(make_pair(4,4)); //3 jet bin
    v_binranges.push_back(make_pair(3,7)); //>=2 jet bin
    v_binranges.push_back(make_pair(1,7)); //all jet bins
  }
  
  for(unsigned int i_bins = 0; i_bins <  v_binranges.size(); i_bins++) {

    int lowBin  = v_binranges.at(i_bins).first;
    int highBin = v_binranges.at(i_bins).second;

    if(latex) {
      if(combineJetBins) 
	cout << "\\textrm{Combined:} " << endL << endl;
      else {
	if(i_bins <4)
	cout << "\\textrm{" << i_bins << " b-tagged Jet Bin:}" << endL << endl;
	else if(i_bins == 4) 
	  cout << "\\textrm{$>=$ 2 b-tagged Jets:}" << endL << endl;
	else if(i_bins == 5)
	  cout << "\\textrm{Total ($>=$ 0 b-tagged Jets):} " << endL << endl;
	cout << "\\textrm{ } " << endL << endl;
	cout << "\\begin{tabular}{l |  c  c  c  c}" << endl;
	cout << "\\hline" << endl;
	cout << "Sample " << colSep  << "ee" << colSep << mathSep << "\\mu\\mu" << mathSep
	     << colSep << "e" << mathSep << "\\mu" << mathSep << colSep << "all" << endL << "\\hline" << endl;
      }
    } else {
      cout << "|  *sample* " << colSep << " *ee* " << colSep << " *mumu* " << colSep 
	   << " *emu* " << colSep << " *all* " <<  endL << endl;
    }

    
    //combine DY samples
    if(combineDYsamples) {
      if(latex) 
	cout << Form("%9s ","DY") << colSep;
      else 
	cout << "|  *DY* " << colSep;

      float n_ee	= 0;
      float n_mm 	= 0;
      float n_em 	= 0;
      float n_all 	= 0;
      float nE_ee 	= 0;
      float nE_mm 	= 0;
      float nE_em 	= 0;
      float nE_all 	= 0;


      for(unsigned int i = 0; i < v_prefixes.size(); i++) {    

	if(!(v_prefixes.at(i).Contains("DY")))
	  continue;

	TH1F *hee = (TH1F*)gDirectory->Get(v_prefixes.at(i) + "_hnBtagJet_allj_ee");
	TH1F *hmm = (TH1F*)gDirectory->Get(v_prefixes.at(i) + "_hnBtagJet_allj_mm");
	TH1F *hem = (TH1F*)gDirectory->Get(v_prefixes.at(i) + "_hnBtagJet_allj_em");
	
	n_ee = n_ee + GetEntries(hee, lowBin, highBin);
	n_mm = n_mm + GetEntries(hmm, lowBin, highBin);
	n_em = n_em + GetEntries(hem, lowBin, highBin);
	nE_ee = sqrt(pow(nE_ee,2) + pow(GetTotalError(hee, lowBin, highBin),2));
	nE_mm = sqrt(pow(nE_mm,2) + pow(GetTotalError(hmm, lowBin, highBin),2));
	nE_em = sqrt(pow(nE_em,2) + pow(GetTotalError(hem, lowBin, highBin),2));

	n_all = n_ee + n_mm + n_em;
	nE_all = sqrt(pow(nE_ee, 2) + pow(nE_mm, 2) + pow(nE_em, 2) );
      }
      cout	<< formatFloat(n_ee,formatS)	<< pmSign	<< formatFloat(nE_ee, formatS)	<< colSep;
      cout	<< formatFloat(n_mm, formatS)	<< pmSign	<< formatFloat(nE_mm, formatS)	<< colSep;
      cout	<< formatFloat(n_em, formatS)	<< pmSign	<< formatFloat(nE_em, formatS)	<< colSep;
      cout	<< formatFloat(n_all, formatS)	<< pmSign	<< formatFloat(nE_all, formatS) << 	endL	<< endl;
    }//combineDYsamples
   

    if(haveVVsamples && combineVVsamples) {
      if(latex) 
	cout << Form("%9s","VV") << colSep;
      else 
	cout << beginL << " *VV* " << colSep;

      float n_ee	= 0;
      float n_mm 	= 0;
      float n_em 	= 0;
      float n_all 	= 0;
      float nE_ee 	= 0;
      float nE_mm 	= 0;
      float nE_em 	= 0;
      float nE_all 	= 0;
      
      for(unsigned int i = 0; i < v_prefixes.size(); i++) {    

      
	if(!(v_prefixes.at(i).Contains("WW") || v_prefixes.at(i).Contains("WZ") || v_prefixes.at(i).Contains("ZZ")))
	  continue;
	TH1F *hee = (TH1F*)gDirectory->Get(v_prefixes.at(i) + "_hnBtagJet_allj_ee");
	TH1F *hmm = (TH1F*)gDirectory->Get(v_prefixes.at(i) + "_hnBtagJet_allj_mm");
	TH1F *hem = (TH1F*)gDirectory->Get(v_prefixes.at(i) + "_hnBtagJet_allj_em");


	n_ee = n_ee + GetEntries(hee, lowBin, highBin);
	n_mm = n_mm + GetEntries(hmm, lowBin, highBin);
	n_em = n_em + GetEntries(hem, lowBin, highBin);
	nE_ee = sqrt(pow(nE_ee,2) + pow(GetTotalError(hee, lowBin, highBin),2));
	nE_mm = sqrt(pow(nE_mm,2) + pow(GetTotalError(hmm, lowBin, highBin),2));
	nE_em = sqrt(pow(nE_em,2) + pow(GetTotalError(hem, lowBin, highBin),2));

	n_all = n_ee + n_mm + n_em;
	nE_all = sqrt(pow(nE_ee, 2) + pow(nE_mm, 2) + pow(nE_em, 2) );
	cout << formatFloat(n_ee, formatS) << pmSign << formatFloat(nE_ee, formatS) << colSep;
	cout << formatFloat(n_mm, formatS) << pmSign << formatFloat(nE_mm, formatS) << colSep;
	cout << formatFloat(n_em, formatS) << pmSign << formatFloat(nE_em, formatS) << colSep;
	cout << formatFloat(n_all, formatS) << pmSign << formatFloat(nE_all, formatS) << endL  << endl;
      }
    } //combineVVsamples


    //now do the rest              
    for(unsigned int i = 0; i < v_prefixes.size(); i++) {    
    
      float n_ee	= 0;
      float n_mm 	= 0;
      float n_em 	= 0;
      float n_all 	= 0;
      float nE_ee 	= 0;
      float nE_mm 	= 0;
      float nE_em 	= 0;
      float nE_all 	= 0;
      
      if(v_prefixes.at(i).Contains("data"))
	continue;
      if(haveVVsamples && combineVVsamples && (v_prefixes.at(i).Contains("WW") || v_prefixes.at(i).Contains("WZ") || v_prefixes.at(i).Contains("ZZ")))
	continue;
      if(combineDYsamples && (v_prefixes.at(i).Contains("DY")))
	continue;

      //print the names
      if(latex)
	cout << beginL << Form("%9s ",v_prefixes.at(i).Data()) <<  colSep;
      else
	cout << beginL << " *" << v_prefixes.at(i) << "* " << colSep;
      TH1F *hee = (TH1F*)gDirectory->Get(v_prefixes.at(i) + "_hnBtagJet_allj_ee");
      TH1F *hmm = (TH1F*)gDirectory->Get(v_prefixes.at(i) + "_hnBtagJet_allj_mm");
      TH1F *hem = (TH1F*)gDirectory->Get(v_prefixes.at(i) + "_hnBtagJet_allj_em");
      
      n_ee =  GetEntries(hee, lowBin, highBin);
      n_mm =  GetEntries(hmm, lowBin, highBin);
      n_em =  GetEntries(hem, lowBin, highBin);
      nE_ee = GetTotalError(hee, lowBin, highBin);
      nE_mm = GetTotalError(hmm, lowBin, highBin);
      nE_em = GetTotalError(hem, lowBin, highBin);      

      n_all = n_ee + n_mm + n_em;
      nE_all = sqrt(pow(nE_ee, 2) + pow(nE_mm, 2) + pow(nE_em, 2) );
      cout << formatFloat(n_ee, formatS) << pmSign << formatFloat(nE_ee, formatS) << colSep;
      cout << formatFloat(n_mm, formatS) << pmSign << formatFloat(nE_mm, formatS) << colSep;
      cout << formatFloat(n_em, formatS) << pmSign << formatFloat(nE_em, formatS) << colSep;
      cout << formatFloat(n_all, formatS) << pmSign << formatFloat(nE_all, formatS) << endL;
      if((v_prefixes.at(i).Contains("tprime") || i == v_prefixes.size() -1) && latex)
	cout << " \\hline " << endl;
      else if ((v_prefixes.at(i).Contains("wprime") || i == v_prefixes.size() -1) && latex)
        cout << " \\hline " << endl;
      else if ((v_prefixes.at(i).Contains("axigluon") || i == v_prefixes.size() -1) && latex)
	cout << " \\hline " << endl;
      else
	cout << endl;
    }
    

    //now do it again, but for the total MC!
    if(latex)
      cout << Form("%9s ","Total MC") << colSep;
    else
      cout <<"|  *Total MC*  " << colSep;
    float n_ee	= 0;
    float n_mm 	= 0;
    float n_em 	= 0;
    float n_all 	= 0;
    float nE_ee 	= 0;
    float nE_mm 	= 0;
    float nE_em 	= 0;
    float nE_all = 0;

    for(unsigned int i = 0; i < v_prefixes.size(); i++) {    
     
     
      if(v_prefixes.at(i).Contains("data"))
	continue;
      if(v_prefixes.at(i).Contains("tprime"))
        continue;
      if(v_prefixes.at(i).Contains("wprime"))
        continue;
      if(v_prefixes.at(i).Contains("axigluon"))
        continue;
     
      TH1F *hee = (TH1F*)gDirectory->Get(v_prefixes.at(i) + "_hnBtagJet_allj_ee");
      TH1F *hmm = (TH1F*)gDirectory->Get(v_prefixes.at(i) + "_hnBtagJet_allj_mm");
      TH1F *hem = (TH1F*)gDirectory->Get(v_prefixes.at(i) + "_hnBtagJet_allj_em");
     
      n_ee = n_ee + GetEntries(hee, lowBin, highBin);
      n_mm = n_mm + GetEntries(hmm, lowBin, highBin);
      n_em = n_em + GetEntries(hem, lowBin, highBin);
      nE_ee = sqrt(pow(nE_ee,2) + pow(GetTotalError(hee, lowBin, highBin),2));
      nE_mm = sqrt(pow(nE_mm,2) + pow(GetTotalError(hmm, lowBin, highBin),2));
      nE_em = sqrt(pow(nE_em,2) + pow(GetTotalError(hem, lowBin, highBin),2));
     
      n_all = n_ee + n_mm + n_em;
      nE_all = sqrt(pow(nE_ee, 2) + pow(nE_mm, 2) + pow(nE_em, 2) );
    }
    cout << formatFloat(n_ee, formatS) << pmSign << formatFloat(nE_ee, formatS) << colSep;
    cout << formatFloat(n_mm, formatS) << pmSign << formatFloat(nE_mm, formatS) << colSep;
    cout << formatFloat(n_em, formatS) << pmSign << formatFloat(nE_em, formatS) << colSep;
    cout << formatFloat(n_all, formatS) << pmSign << formatFloat(nE_all, formatS) << endL;
    if(latex)
      cout << " \\hline" << endl;
    else 
      cout << endl;


    if(haveData) {
      
      if(latex)
	cout << beginL << Form("%9s ","Data") << colSep; 
      else
	cout << beginL << " *Data* " << colSep;
      for(unsigned int i = 0 ; i < 4; i++) {
	TString name = "data_hnBtagJet_allj_" + TString(suffix[i]);
	TH1F *h = (TH1F*)gDirectory->Get(name.Data());
	float n = 0;
	float nE = 0;

	n  = GetEntries(h, lowBin, highBin);
	nE = GetTotalError(h, lowBin, highBin);

	if (!printErrorsForData)
	  cout <<  (int)n;
	else
	  cout << formatFloat(n,formatS) << pmSign << formatFloat(nE, formatS);
	if(i!=3)
	  cout << colSep;
	else {
	  if(latex)
	    cout << endL << " \\hline" << endl;
	  else
	    cout << endL << endl;
	}
      }
    }//have data

    if(latex) {
      cout << "\\end{tabular}" << endl;
      cout << "\\vspace{2em} " << endL << endl;
    }

  }//bin ranges
    
}

