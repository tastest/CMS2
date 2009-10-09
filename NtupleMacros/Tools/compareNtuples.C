/*
  Macro to validate and compare ntuples. The macro loops through the ntuples and 
  plots common branches using TTree's Draw function. The macro does NOT plot 
  vectors of vectors or TStrings. To run, do:
  
  [1934] kalavase@uaf-6 ~/temp$ root -b 
  root [0] .L compareNtuples.C++
  root [1] compareNtuples("ntuple_1.root", "ntuple_2.root", true, true)  

  The first two arguments are self-explanatory - the files that you want to compare.
  The last 2 arguments are true by default. The third argument, if set false, will save 
  every canvas, even if the two histograms that are the result of the Tree Draw are 
  identical. If the
*/
  

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TClass.h"
#include "TKey.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include <fstream>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

vector<TString> getAliasNames(TTree *t) {

  vector<TString> v_aliasnames;
  
  TList *t_list =  t->GetListOfAliases();
  for(unsigned int i = 0; i < t_list->GetSize(); i++) {
    TString aliasname(t_list->At(i)->GetName());
    TBranch *branch = t->GetBranch(t->GetAlias(aliasname.Data()));
    TString branchname(branch->GetName());
    TString branchtitle(branch->GetTitle());
    if(branchname.BeginsWith("intss") ||
       branchname.BeginsWith("floatss") ||
       branchname.BeginsWith("doubless") ||
       branchname.Contains("LorentzVectorss") ) {
      
      cout << "Sorry, I dont know about vector of vectors of objects" 
	   << "Will be skipping " << aliasname << endl;
      
      continue;
    }

    if(branchname.Contains("TString") ) {
      cout << "Sorry, I don't know how to graphically represent TStrings in only 3 dimensions" 
	   << " Put in a feature request. Will be skipping " << aliasname << endl;
      continue;
    }
      v_aliasnames.push_back(aliasname);
  }

  return v_aliasnames;
}

//----------------------------------------------------------------------
bool areHistosTheSame(TH1F* h1, TH1F* h2) {
  
  if(h1->GetNbinsX() != h2->GetNbinsX()) 
    return false;
  
  //make sure that the bin range is the same
  float range1 = h1->GetBinCenter(1) - h1->GetBinCenter(h1->GetNbinsX());
  float range2 = h2->GetBinCenter(1) - h2->GetBinCenter(h2->GetNbinsX());

  if(TMath::Abs(range1 - range2) > 0.000001) 
    return false;
  
  for(unsigned int i = 1; i < h1->GetNbinsX()+1; i++) {
    if(TMath::Abs(h1->GetBinContent(i) - h2->GetBinContent(i)) > 0.000001) 
      return false;
  }
  
  
  return true;
}
//-----------------------------------------------------------------------

void compareNtuples(TString file1, TString file2, bool doNotSaveSameHistos="true", bool drawWithErrors="true") {
  
  //TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  
  TFile *f1 = TFile::Open(file1.Data(), "READ");
  if(f1 == NULL) {
    cout << "Exiting" << endl;
    return;
  }
  TTree *tree1 = (TTree*)f1->Get("Events");
  if(tree1 == NULL) {
    cout << "Can't find the tree \"Events\" in " << file1 << " Exiting " << endl;
    return;
  }
    
  

  TFile *f2 = TFile::Open(file2.Data(), "READ");
  if(f2 == NULL) {
    cout << "Exiting" << endl;
    return;
  }
  TTree *tree2 = (TTree*)f2->Get("Events");
  if(tree2 == NULL) {
    cout << "Can't find the tree \"Events\" in " << file2 << " Exiting " << endl;
    return;
  }


  
  vector<TString> t1_aliasnames = getAliasNames(tree1);
  vector<TString> t2_aliasnames = getAliasNames(tree2);
  vector<TString> v_commonBranches;
  
  
  for(vector<TString>::iterator it = t1_aliasnames.begin(); 
      it != t1_aliasnames.end(); it++) {
    
    if(find(t2_aliasnames.begin(), t2_aliasnames.end(), *it) != t2_aliasnames.end())
      v_commonBranches.push_back(*it);
  }



  
  for(unsigned int i =  0; i < v_commonBranches.size(); i++) {
        
    cout << "Comparing Branch: " << v_commonBranches.at(i) << endl;
    TString hist1name = "h1_"+(v_commonBranches.at(i));
    TString hist2name = "h2_"+(v_commonBranches.at(i));
    TString command1 = (v_commonBranches.at(i))+">>"+hist1name;
    TString command2 = (v_commonBranches.at(i))+">>"+hist2name;
    
    if((v_commonBranches.at(i)).Contains("p4")) {
      hist1name = "h1_"+(v_commonBranches.at(i)) + "_Pt";
      hist2name = "h2_"+(v_commonBranches.at(i)) + "_Pt";
      command1 = (v_commonBranches.at(i)) + ".Pt()>>"+hist1name;
      command2 = (v_commonBranches.at(i)) + ".Pt()>>"+hist2name;
    }
    TCanvas *c1 = new TCanvas(); 
    tree1->Draw(command1.Data());
    TH1F *h1 = (TH1F*)c1->GetPrimitive(hist1name.Data());
    if(h1==NULL) {
      cout << "********** Branch " << v_commonBranches.at(i) 
	   << " in file " << file1 << "exists, but is undrawable for some reason. " 
	   << "Skipping this branch" << endl;
      c1->Clear();
      continue; 
    }
    //set the Overflow at the last bin
    tree2->Draw(command2.Data());
    TH1F *h2 = (TH1F*)c1->GetPrimitive(hist2name.Data());
    if(h2==NULL) {
      cout << "********** Branch " << v_commonBranches.at(i) 
	   << " in file " << file2 << "exists, but is undrawable for some reason. " 
	   << "Skipping this branch" << endl;
      c1->Clear();
      continue;  
    }
    c1->Clear();

    if(areHistosTheSame(h1, h2) && doNotSaveSameHistos)
      continue;
    
    if(drawWithErrors) {
      h1->TH1F::Sumw2();
      h2->TH1F::Sumw2();
    }


    h1->Scale(1./h1->GetEntries());
    h2->Scale(1./h2->GetEntries());


    if(h1->GetNbinsX() != h2->GetNbinsX() ) {
      cout << "Branch " << v_commonBranches.at(i) << " not the same between the 2 files" 
	   << ". They will be drawn side by side" << endl;

      if(!drawWithErrors) {
	h1->SetLineColor(0);
	h1->SetMarkerSize(1.1);
	h1->SetMarkerStyle(3);
	h2->SetLineColor(kRed);
	h2->SetLineStyle(7);
      } else {
	h1->SetMarkerSize(1.3);
	h1->SetMarkerStyle(3);
	h2->SetMarkerSize(1.1);
	h2->SetMarkerStyle(8);
	h2->SetMarkerColor(kRed);
      }

      c1->Divide(2,1);
      c1->cd(1);
      h1->Draw();
      c1->cd(2);
      h2->Draw();
      if(i < v_commonBranches.size() - 1)
	c1->SaveAs("diff.eps(");
      else {
	cout << "done" << endl;
	c1->SaveAs("diff.eps)");
      }
      continue;
    }        
            
    if(h1->GetMaximum() >= h2->GetMaximum()) {
      
      double max = 1.1*h1->GetMaximum();
      h1->SetMaximum(max);
      h2->SetMaximum(max);
      
      if(!drawWithErrors) {
	h1->SetLineColor(0);
	h1->SetMarkerSize(1.1);
	h1->SetMarkerStyle(3);
	h2->SetLineColor(kRed);
	h2->SetLineStyle(7);
	h2->Draw();
	h1->Draw("samesh*");
      } else {
	h1->SetMarkerSize(1.3);
	h1->SetMarkerStyle(3);
	h2->SetMarkerSize(1.1);
	h2->SetMarkerStyle(8);
	h2->SetMarkerColor(kRed);
	h2->Draw("e");
	h1->Draw("samese");
      }
      
    } else {
      double max = 1.1*h2->GetMaximum();
      
      h1->SetMaximum(max);
      h2->SetMaximum(max);
      
       if(!drawWithErrors) {
	h1->SetLineColor(0);
	h1->SetMarkerSize(1.1);
	h1->SetMarkerStyle(3);
	h2->SetLineColor(kRed);
	h2->SetLineStyle(7);
	h1->Draw();
	h2->Draw("samesh*");
      } else {
	h1->SetMarkerSize(1.3);
	h1->SetMarkerStyle(3);
	h2->SetMarkerSize(1.1);
	h2->SetMarkerStyle(8);
	h2->SetMarkerColor(kRed);
	h1->Draw("e");
	h2->Draw("samese");
      }


    }
        
    if(i < v_commonBranches.size() - 1)
      c1->SaveAs("diff.eps(");
    else {
      cout << "done" << endl;
      c1->Clear();
      c1->SaveAs("diff.eps)");
    }
    
  }//for loop
    
    
}
