/* 
   This macro takes in the locations of 2 CMS2 ntuple files 
   and ouputs a macro, called compareFiles.C. The outputed macro will do a tree Draw
   for each branch in common to the two trees and will then subtract the 2 distributions.
   If there are distributions that are not the same, the subtracted histogram will have 
   entries that are not 0, and the macro will report those branches to you.

   To run do in a root session:
   .L makeValidationMacro.C++
   makeValidationMacro()
   
   Then, preferably in a new root session running in batch mode (root -b), do:
   .L compareFiles.C
   compareFiles()

   In your directory, you will see an eps file called diff.eps which will contain the diffed 
   histograms.
*/




#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TClass.h"
#include "TKey.h"
#include "TTree.h"
#include "TCanvas.h"
#include <fstream>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

//if you just have one file, and you want to draw every branch
void makeValidationMacro(TString f1) {

  ofstream myfile;
  myfile.open ("drawBranches.C");
  myfile << "#include \"TH1F.h\"" << endl;
  myfile << "#include \"TTree.h\"" << endl;
  myfile << "#include \"TFile.h\"" << endl;
  myfile << "#include <vector>"  << endl;
  myfile << "using namespace std;" << endl;
  myfile << "TCanvas *c = new TCanvas();" << endl;
  myfile << "\nvoid drawBranches() {" << endl;
  myfile << "   TFile *f1 = TFile::Open(\"" << f1 << "\", \"READ\");\n";
  myfile << "   TTree *Events1 = dynamic_cast<TTree*>(f1->Get(\"Events\"));\n\n";
  

  TFile *file1 = new TFile(f1.Data(), "READ");
  TTree *Events1 = dynamic_cast<TTree*>(file1->Get("Events"));
  TList *list1   = Events1->GetListOfAliases();
  TIterator *iter1 = list1->MakeIterator();
  vector<TString> bNames1;
  TKey *key;
  while(key=(TKey*)iter1->Next()) {
    TString name(key->GetName());
    if(name.Contains("hyp_jets") ||
       name.Contains("hyp_other_jets") 
       || name.Contains("genps_lepdaughter_p4") 
       || name.Contains("hlt8e29_trigObjs_p4") 
       || name.Contains("hlt_trigObjs_p4") 
       || name.Contains("taus_pf_isochargecand_p4") 
       || name.Contains("taus_pf_isogammacand_p4") 
       || name.Contains("taus_pf_isoneutrcand_p4") 
       || name.Contains("taus_pf_sigchargecand_p4") 
       || name.Contains("taus_pf_siggammacand_p4") 
       || name.Contains("taus_pf_signeutrcand_p4") 
       )
       continue;
    bNames1.push_back(key->GetName());
  }
   //loop over the  branches and make histos
  for(vector<TString>::iterator v_it = bNames1.begin();
      v_it != bNames1.end(); v_it++) {
    
    TString branch(*v_it);
    
    TString suffix1 = f1.ReplaceAll(".root", "");
    suffix1 = suffix1.Tokenize("/")->Last()->GetName();
    
    TString histoname1 = "h_" + branch + "_" + suffix1;
    
    if(branch.Contains("hyp_jets"))
      continue;
    if(branch.Contains("hyp_other_jets"))
      continue;
    if(branch.Contains("trigNames"))
      continue;
    if(branch.Contains("evt_dataset"))
      continue;
    if(branch.Contains("vtxs_position"))
      continue;


    if(branch.Contains("p4") || branch.Contains("P4") ) {
      myfile << "   Events1->Draw(\"" << branch << ".Pt()>>" << histoname1 << "Pt\");" << endl;
    } else {
      myfile << "   Events1->Draw(\"" << branch << ">>" << histoname1 << "\");" << endl;
    }

    if(v_it + 1 != bNames1.end())
      myfile << "   c->SaveAs(\"" << f1.ReplaceAll(".root","") << "_histos.eps(\");" << endl;
  }
  
  myfile << "   c->SaveAs(\"" << f1.ReplaceAll(".root","") << "_histos.eps)\");" << endl;
  myfile << "}";
  myfile.close();
  file1->Close();
  
    
}

//if you have 2 files and want to compare the common branches 
//in the two files
void makeValidationMacro(TString f1, TString f2) {


  ofstream myfile;
  myfile.open ("compareFiles.C");
  myfile << "#include \"TH1F.h\"" << endl;
  myfile << "#include \"TTree.h\"" << endl;
  myfile << "#include \"TFile.h\"" << endl;
  myfile << "#include <vector>"  << endl;
  myfile << "using namespace std;" << endl;
  
  myfile << "vector<TString> v_badhistos;" << endl; 
  myfile << "TCanvas *c = new TCanvas();" << endl;
  myfile << "bool compareHistos(TH1F *h1, TH1F *h2, TString diffHistoName) {" << endl;
  myfile << "   bool res = true;" << endl;
  myfile << "   int nbins1 = h1->GetNbinsX();" << endl;
  myfile << "   int nbins2 = h2->GetNbinsX();" << endl;
  myfile << "   if(nbins1 != nbins2) {" << endl;
  myfile << "      v_badhistos.push_back(h1->GetTitle());" << endl;
  myfile << "      res = false; " << endl;
  myfile << "   }" << endl;
  myfile << "   if (res){" << endl;
  myfile << "      for(int i = 0; i < nbins1; i++) { " << endl;
  myfile << "         if(h1->GetBinContent(i) - h2->GetBinContent(i) != 0) {" << endl;
  myfile << "            v_badhistos.push_back(h1->GetTitle());" << endl;
  myfile << "            res =  false; " << endl;
  myfile << "         }" << endl;
  myfile << "      }" << endl;
  myfile << "   }" << endl;
  myfile << "   if (!res){" << endl;
  myfile << "      //TH1F *h = dynamic_cast<TH1F*>(h1->Clone());" << endl;
  myfile << "      //h->SetName((diffHistoName+\"_diff\").Data());" << endl;
  myfile << "      //h->Add(h1, h2, 1, -1);" << endl;
  myfile << "      h1->Draw();" << endl;
  myfile << "      h2->Draw(\"sames\");" << endl;
  myfile << "      c->SaveAs(\"diff.eps(\");";
  myfile << "      cout << \"  ERRORDIFF: \"<< h1->GetName() << \"    \"<< h2->GetName() <<endl;" << endl;
  myfile << "   }" << endl;
  myfile << "   return res;" << endl;
  myfile << "} " << endl;
  
   
  TFile *file1 = new TFile(f1.Data(), "READ");
  TTree *Events1 = dynamic_cast<TTree*>(file1->Get("Events"));
  TList *list1   = Events1->GetListOfAliases();
  TIterator *iter1 = list1->MakeIterator();
  vector<TString> bNames1;
  
  TFile *file2 = new TFile(f2.Data(), "READ");
  TTree *Events2 = dynamic_cast<TTree*>(file2->Get("Events"));
  TList *list2   = Events2->GetListOfAliases();
  TIterator *iter2 = list2->MakeIterator();
  vector<TString> bNames2;
  
  vector<TString> commonBranches;
  
  TKey *key;
  while(key=(TKey*)iter1->Next()) {
    TString name(key->GetName());
    if(name.Contains("hyp_jets") ||
       name.Contains("hyp_other_jets") 
       || name.Contains("genps_lepdaughter_p4") 
       || name.Contains("hlt8e29_trigObjs_p4") 
       || name.Contains("hlt_trigObjs_p4") 
       || name.Contains("taus_pf_isochargecand_p4") 
       || name.Contains("taus_pf_isogammacand_p4") 
       || name.Contains("taus_pf_isoneutrcand_p4") 
       || name.Contains("taus_pf_sigchargecand_p4") 
       || name.Contains("taus_pf_siggammacand_p4") 
       || name.Contains("taus_pf_signeutrcand_p4") 
       )
      continue;
    bNames1.push_back(key->GetName());
  }

  while(key=(TKey*)iter2->Next()) {
    TString name(key->GetName());
    if(name.Contains("hyp_jets") ||
       name.Contains("hyp_other_jets") 
       || name.Contains("genps_lepdaughter_p4") 
       || name.Contains("hlt8e29_trigObjs_p4") 
       || name.Contains("hlt_trigObjs_p4") 
       || name.Contains("taus_pf_isochargecand_p4") 
       || name.Contains("taus_pf_isogammacand_p4") 
       || name.Contains("taus_pf_isoneutrcand_p4") 
       || name.Contains("taus_pf_sigchargecand_p4") 
       || name.Contains("taus_pf_siggammacand_p4") 
       || name.Contains("taus_pf_signeutrcand_p4") 
       )
      continue;
    
    bNames2.push_back(key->GetName());
  }

  
  if(bNames2.size() > bNames1.size()) {
    for(vector<TString>::iterator v_it = bNames2.begin();
	v_it != bNames2.end(); v_it++) {
      if(find(bNames1.begin(), bNames1.end(), *v_it) != bNames1.end()) 
	commonBranches.push_back(*v_it);
    }
  } else {
    for(vector<TString>::iterator v_it = bNames1.begin();
	v_it != bNames1.end(); v_it++) {
      if(find(bNames2.begin(), bNames2.end(), *v_it) != bNames2.end()) 
	commonBranches.push_back(*v_it);
    }
  }

  myfile << "\nvoid compareFiles() {" << endl;
  myfile << "   TFile *f1 = TFile::Open(\"" << f1 << "\", \"READ\");\n";
  myfile << "   TTree *Events1 = dynamic_cast<TTree*>(f1->Get(\"Events\"));\n\n";
  myfile << "   TFile *f2 = TFile::Open(\"" << f2 << "\", \"READ\");\n";
  myfile << "   TTree *Events2 = dynamic_cast<TTree*>(f2->Get(\"Events\"));\n\n";
  
  
  //loop over the common branches and make histos
  for(vector<TString>::iterator v_it = commonBranches.begin();
      v_it != commonBranches.end(); v_it++) {
    
    TString branch(*v_it);
    
    TString suffix1 = f1.ReplaceAll(".root", "");
    TString suffix2 = f2.ReplaceAll(".root", "");
    suffix1 = suffix1.Tokenize("/")->Last()->GetName();
    suffix2 = suffix2.Tokenize("/")->Last()->GetName();
    
    TString histoname1 = "h1_" + branch + "_" + suffix1;
    TString histoname2 = "h2_" + branch + "_" + suffix2;

    if(branch.Contains("hyp_jets"))
      continue;
    if(branch.Contains("hyp_other_jets"))
      continue;
    if(branch.Contains("trigNames"))
      continue;
    if(branch.Contains("evt_dataset"))
      continue;
    if(branch.Contains("vtxs_position"))
      continue;

    myfile << "   cout << \"comparing Branches: " << histoname1 
	   << " and " << histoname2 << "\" << endl;" << endl;
    if(branch.Contains("p4") || branch.Contains("P4") ) {
      myfile << "   Events1->Draw(\"" << branch << ".Pt()>>" << histoname1 << "Pt\");" << endl;
      myfile << "   Events2->Draw(\"" << branch << ".Pt()>>" << histoname2 << "Pt\");" << endl;
      myfile << "   compareHistos(" << histoname1 << "Pt," << histoname2 << "Pt,\"h_" 
	   << branch << "Pt\");" << endl;
    } else {
      myfile << "   Events1->Draw(\"" << branch << ">>" << histoname1 << "\");" << endl;
      myfile << "   Events2->Draw(\"" << branch << ">>" << histoname2 << "\");" << endl;
      myfile << "   compareHistos(" << histoname1 << "," << histoname2 << ",\"h_" 
	   << branch << "\");" << endl;
    }
        
  }
  
//   myfile << "   TCanvas *c1 = new TCanvas();" << endl;
//   myfile << "   c1->cd(); " << endl;
//   myfile << "   unsigned int index = 0;" << endl;
//   myfile << "   for(vector<TH1F> v_it = v_diffhistos.begin(); v_it != v_diffhistos.end(); v_it++, index++) {" 
// 	 << endl;
//   myfile << "      v_it->Draw();" << endl;
//   myfile << "      if(index + 1 != v_diffhistos.size())" << endl;
//   myfile << "         c1->SaveAs(\"diff.eps(\");" << endl;
//   myfile << "      else " << endl;
//   myfile << "         c1->SaveAs(\"diff.eps)\")" << endl;
//   myfile << "   }" << endl;

  myfile << "   for(unsigned int i = 0; i < v_badhistos.size(); i++) {" << endl;
  myfile << "       cout << \"Branch \" << v_badhistos.at(i) << \" does not agree between the two files\" << endl;" << endl;
  myfile <<    "}" << endl;
  

  myfile << "   c->SaveAs(\"diff.eps)\");";
  myfile << "}";
  myfile.close();
  file1->Close();
  file2->Close();

  

}
