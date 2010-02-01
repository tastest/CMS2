// Original author: Puneeth Kalavase (UCSB)
// 
/*ROOT macro to make CMS2.h and ScanChain.C files for basic analysis
  of CMS2 ntuples. Usage:

  [kalavase@stau ~/rootmacros]$ root
  root [0] .L makeCMS2ClassFiles.C++
  // 
  // the second and third arguments are optional
  //
  // if we are paranoid, each float, vector<float> and
  // vector<vector<float> > is checked for NaNs and infs
  //
  // The classname is CMS2 by default
  //
  root [1] makeCMS2ClassFiles("merged_ntuple.root",true,"classname")

*/


#include "TBranch.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>
#include <vector>
#include "Math/LorentzVector.h"

using namespace std;
ofstream headerf;
ofstream codef;
ofstream implf;
ofstream branchfile;

void makeHeaderFile(TFile *f, bool paranoid, string Classname) {
	
  
  
  headerf << "// -*- C++ -*-" << endl;
  headerf << "#ifndef " << Classname << "_H" << endl;
  headerf << "#define " << Classname << "_H" << endl;
  headerf << "#include \"Math/LorentzVector.h\"" << endl;
  headerf << "#include \"Math/Point3D.h\"" << endl;
  headerf << "#include \"TMath.h\"" << endl;
  headerf << "#include \"TBranch.h\"" << endl;
  headerf << "#include \"TTree.h\"" << endl;
  headerf << "#include \"TH1F.h\""  << endl;
  headerf << "#include \"TFile.h\"" << endl;
  headerf << "#include <vector> " << endl << endl;
  if (paranoid)
    headerf << "#define PARANOIA" << endl << endl;
  headerf << "using namespace std; " << endl;
  headerf << "class " << Classname << " {" << endl;
  headerf << "private: " << endl;
  headerf << "protected: " << endl;
  headerf << "\tunsigned int index;" << endl;
  TTree *ev = (TTree*)f->Get("Events");
  TList *fullarray =  ev->GetListOfAliases();
  TList *aliasarray = new TList();
    
  //for(Int_t i = 0; i < aliasarray->GetSize(); ++i) {
  for(Int_t i = 0; i < fullarray->GetSize(); ++i) {
    TString aliasname(fullarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    TString branchname(branch->GetName());
    TString branchtitle(branch->GetTitle());
    if(!branchname.BeginsWith("int") && 
       !branchname.BeginsWith("uint") && 
       !branchname.BeginsWith("float") &&
       !branchname.BeginsWith("double") &&
       !branchtitle.EndsWith("/F") && 
       !branchtitle.EndsWith("/I") &&
       !branchtitle.EndsWith("/i") &&
       !branchtitle.BeginsWith("TString"))
      continue;
    aliasarray->Add(fullarray->At(i));
  }
  
  
  for(Int_t i = 0; i< aliasarray->GetSize(); ++i) {
    
    //Class name is blank for a int of float
    TString aliasname(aliasarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    
    TString classname = branch->GetClassName();
    TString title     = branch->GetTitle();
    if ( classname.Contains("vector") ) {
      if(classname.Contains("edm::Wrapper<") ) {
	classname = classname(0,classname.Length()-2);
	classname.ReplaceAll("edm::Wrapper<","");
	headerf << "\t" << classname << " " << aliasname << "_;" << endl;
      } else {
	headerf << "\t" << classname << " *" << aliasname << "_;" << endl;
      }
    } else {
      
      if(classname != "" ) { //LorentzVector
	if(classname.Contains("edm::Wrapper<") ) {
	  classname = classname(0,classname.Length()-1);
	  classname.ReplaceAll("edm::Wrapper<","");
	  headerf << "\t" << classname << " " << aliasname << "_;" << endl;
	} else {
	  headerf << "\t" << classname << " *" << aliasname << "_;" << endl;
	}
      } else {
	if(title.EndsWith("/i"))
	  headerf << "\tunsigned int" << "\t" << aliasname << "_;" << endl;
	if(title.EndsWith("/F"))
	  headerf << "\tfloat" << "\t" << aliasname << "_;" << endl;
	if(title.EndsWith("/I"))
	  headerf << "\tint" << "\t" << aliasname << "_;" << endl;
      }
    }
    headerf << "\tTBranch *" << Form("%s_branch",aliasname.Data()) << ";" << endl;
    headerf << "\tbool " << Form("%s_isLoaded",aliasname.Data()) << ";" << endl;
  }
  
  
  headerf << "public: " << endl;
  headerf << "int ScanChain(class TChain* chain, int nEvents=-1, std::string skimFilePrefix=\"\");" << endl;
  headerf << "void Init(TTree *tree) {" << endl;
    

  // SetBranchAddresses for LorentzVectors
  for(Int_t i = 0; i< aliasarray->GetSize(); i++) {
    TString aliasname(aliasarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    TString classname = branch->GetClassName();
    if ( !classname.Contains("vector<vector") ) {
      if ( classname.Contains("Lorentz") || classname.Contains("PositionVector") ) {
	headerf << "\t" << Form("%s_branch",aliasname.Data()) << " = 0;" << endl;
	headerf << "\t" << "if (tree->GetAlias(\"" << aliasname << "\") != 0) {" << endl;
	headerf << "\t\t" << Form("%s_branch",aliasname.Data()) << " = tree->GetBranch(tree->GetAlias(\"" << aliasname << "\"));" << endl;
	headerf << "\t\t" << Form("%s_branch",aliasname.Data()) << "->SetAddress(&" << aliasname << "_);" << endl << "\t}" << endl;
	headerf << "\t" << "if(" << Form("%s_branch",aliasname.Data()) << " == 0 ) {" << endl;
	headerf << "\t" << "cout << \"Branch " << aliasname.Data() << " does not exist.\" << endl;" << endl;
	headerf << "\t" << "}" << endl;
      }
    }
  }


  // SetBranchAddresses for everything else
  headerf << "  tree->SetMakeClass(1);" << endl;
  for(Int_t i = 0; i< aliasarray->GetSize(); i++) {
    TString aliasname(aliasarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    TString classname = branch->GetClassName();
    if ( ! (classname.Contains("Lorentz") || classname.Contains("PositionVector")) || classname.Contains("vector<vector") ) {
      headerf << "\t" << Form("%s_branch",aliasname.Data()) << " = 0;" << endl;
      headerf << "\t" << "if (tree->GetAlias(\"" << aliasname << "\") != 0) {" << endl;
      headerf << "\t\t" << Form("%s_branch",aliasname.Data()) << " = tree->GetBranch(tree->GetAlias(\"" << aliasname << "\"));" << endl;
      headerf << "\t\t" << Form("%s_branch",aliasname.Data()) << "->SetAddress(&" << aliasname << "_);" << endl << "\t}" << endl;
      headerf << "\t" << "if(" << Form("%s_branch",aliasname.Data()) << " == 0 ) {" << endl;
      headerf << "\t" << "cout << \"Branch " << aliasname.Data() << " does not exist.\" << endl;" << endl;
      headerf << "\t" << "}" << endl;
    }
  }

  headerf << "  tree->SetMakeClass(0);" << endl;
  headerf << "}" << endl;

  // GetEntry
  headerf << "void GetEntry(unsigned int idx) " << endl;
  headerf << "\t// this only marks branches as not loaded, saving a lot of time" << endl << "\t{" << endl;
  headerf << "\t\tindex = idx;" << endl;
  for(Int_t i = 0; i< aliasarray->GetSize(); i++) {
    TString aliasname(aliasarray->At(i)->GetName());
    headerf << "\t\t" << Form("%s_isLoaded",aliasname.Data()) << " = false;" << endl;
  }
  headerf << "\t}" << endl << endl;

  // LoadAllBranches
  headerf << "void LoadAllBranches() " << endl;
  headerf << "\t// load all branches" << endl << "{" << endl;
  for(Int_t i = 0; i< aliasarray->GetSize(); i++) {
    TString aliasname(aliasarray->At(i)->GetName());
    headerf << "\t" << "if (" << aliasname.Data() <<  "_branch != 0) " << Form("%s();",aliasname.Data()) << endl;
  }
  headerf << "}" << endl << endl;

  // accessor functions
  for (Int_t i = 0; i< aliasarray->GetSize(); i++) {
    TString aliasname(aliasarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    TString classname = branch->GetClassName();
    TString title = branch->GetTitle();
    bool isSkimmedNtuple = false;
    if(!classname.Contains("edm::Wrapper<") &&
       (classname.Contains("vector") || classname.Contains("LorentzVector") ) )
      isSkimmedNtuple = true;
    if ( classname.Contains("vector") ) {
      if(classname.Contains("edm::Wrapper<") ) {
	classname = classname(0,classname.Length()-2);
	classname.ReplaceAll("edm::Wrapper<","");
      }
      headerf << "\t" << classname << " &" << aliasname << "()" << endl;
    } else {
      if(classname.Contains("edm::Wrapper<") ) {
	classname = classname(0,classname.Length()-1);
	classname.ReplaceAll("edm::Wrapper<","");
      }
      if(classname != "" ) {
	headerf << "\t" << classname << " &" << aliasname << "()" << endl;
      } else {
	if(title.EndsWith("/i"))
	  headerf << "\tunsigned int &" << aliasname << "()" << endl;
	if(title.EndsWith("/F"))
	  headerf << "\tfloat &" << aliasname << "()" << endl;
	if(title.EndsWith("/I"))
	  headerf << "\tint &" << aliasname << "()" << endl;
      }
    }
    aliasname = aliasarray->At(i)->GetName();
    headerf << "\t{" << endl;
    headerf << "\t\t" << "if (not " << Form("%s_isLoaded) {",aliasname.Data()) << endl;
    headerf << "\t\t\t" << "if (" << Form("%s_branch",aliasname.Data()) << " != 0) {" << endl;
    headerf << "\t\t\t\t" << Form("%s_branch",aliasname.Data()) << "->GetEntry(index);" << endl;
    if (paranoid) {
      headerf << "\t\t\t\t#ifdef PARANOIA" << endl;
      if (classname == "vector<vector<float> >") {
	if(isSkimmedNtuple) {
	  headerf << "\t\t\t\t" << "for (vector<vector<float> >::const_iterator i = " 
		  << aliasname << "_->begin(); i != "<< aliasname << "_->end(); ++i) {" << endl;
	} else {
	  headerf << "\t\t\t\t" << "for (vector<vector<float> >::const_iterator i = " 
		  << aliasname << "_.begin(); i != "<< aliasname << "_.end(); ++i) {" << endl;
	}
	headerf << "\t\t\t\t\t" << "for (vector<float>::const_iterator j = i->begin(); " 
	  "j != i->end(); ++j) {" << endl;
	headerf << "\t\t\t\t\t\t" << "if (not isfinite(*j)) {" << endl;
	headerf << "\t\t\t\t\t\t\t" << "printf(\"branch " << Form("%s_branch",aliasname.Data()) 
		<< " contains a bad float: %f\\n\", *j);" << endl << "\t\t\t\t\t\t\t" << "exit(1);"
		<< endl;
	headerf << "\t\t\t\t\t\t}\n\t\t\t\t\t}\n\t\t\t\t}" << endl;
      } else if (classname == "vector<float>") {
	if(isSkimmedNtuple) {
	  headerf << "\t\t\t\t" << "for (vector<float>::const_iterator i = " 
		  << aliasname << "_->begin(); i != "<< aliasname << "_->end(); ++i) {" << endl;
	} else {
	  headerf << "\t\t\t\t" << "for (vector<float>::const_iterator i = " 
		  << aliasname << "_.begin(); i != "<< aliasname << "_.end(); ++i) {" << endl;
	}
	headerf << "\t\t\t\t\t" << "if (not isfinite(*i)) {" << endl;
	headerf << "\t\t\t\t\t\t" << "printf(\"branch " << Form("%s_branch",aliasname.Data()) 
		<< " contains a bad float: %f\\n\", *i);" << endl << "\t\t\t\t\t\t" << "exit(1);"
		<< endl;
	headerf << "\t\t\t\t\t}\n\t\t\t\t}" << endl;
      } else if (classname == "float") {
	headerf << "\t\t\t\t" << "if (not isfinite(" << aliasname << "_)) {" << endl;
	headerf << "\t\t\t\t\t" << "printf(\"branch " << Form("%s_branch",aliasname.Data()) 
		<< " contains a bad float: %f\\n\", " << aliasname << "_);" << endl 
		<< "\t\t\t\t\t" << "exit(1);"
		<< endl;
	headerf << "\t\t\t\t}" << endl;
      } else if (classname.BeginsWith("vector<vector<ROOT::Math::LorentzVector")) {
	if(isSkimmedNtuple) {
	  headerf << "\t\t\t\t" << "for (" << classname.Data() <<"::const_iterator i = " 
		  << aliasname << "_->begin(); i != "<< aliasname << "_->end(); ++i) {" << endl;
	} else {
	  headerf << "\t\t\t\t" << "for (" << classname.Data() <<"::const_iterator i = " 
		  << aliasname << "_.begin(); i != "<< aliasname << "_.end(); ++i) {" << endl;
	}
	// this is a slightly hacky way to get rid of the outer vector< > ...
	std::string str = classname.Data() + 7;
	str[str.length() - 2] = 0;
	headerf << "\t\t\t\t\t" << "for (" << str.c_str() << "::const_iterator j = i->begin(); " 
	  "j != i->end(); ++j) {" << endl;
	headerf << "\t\t\t\t\t\t" << "int e;" << endl;
	headerf << "\t\t\t\t\t\t" << "frexp(j->pt(), &e);" << endl;
	headerf << "\t\t\t\t\t\t" << "if (not isfinite(j->pt()) || e > 30) {" << endl;
	headerf << "\t\t\t\t\t\t\t" << "printf(\"branch " << Form("%s_branch",aliasname.Data()) 
		<< " contains a bad float: %f\\n\", j->pt());" << endl << "\t\t\t\t\t\t\t" << "exit(1);"
		<< endl;
	headerf << "\t\t\t\t\t\t}\n\t\t\t\t\t}\n\t\t\t\t}" << endl;
      } else if (classname.BeginsWith("vector<ROOT::Math::LorentzVector")) {
	if(isSkimmedNtuple) {
	  headerf << "\t\t\t\t" << "for (" << classname.Data() << "::const_iterator i = " 
		  << aliasname << "_->begin(); i != "<< aliasname << "_->end(); ++i) {" << endl;
	} else {
	  headerf << "\t\t\t\t" << "for (" << classname.Data() << "::const_iterator i = " 
		  << aliasname << "_.begin(); i != "<< aliasname << "_.end(); ++i) {" << endl;
	}
	headerf << "\t\t\t\t\t" << "int e;" << endl;
	headerf << "\t\t\t\t\t" << "frexp(i->pt(), &e);" << endl;
	headerf << "\t\t\t\t\t" << "if (not isfinite(i->pt()) || e > 30) {" << endl;
	headerf << "\t\t\t\t\t\t" << "printf(\"branch " << Form("%s_branch",aliasname.Data()) 
		<< " contains a bad float: %f\\n\", i->pt());" << endl << "\t\t\t\t\t\t" << "exit(1);"
		<< endl;
	headerf << "\t\t\t\t\t}\n\t\t\t\t}" << endl;
      } else if (classname.BeginsWith("ROOT::Math::LorentzVector")) {
	headerf << "\t\t\t\t" << "int e;" << endl;
	if(isSkimmedNtuple) {
	  headerf << "\t\t\t\t" << "frexp(" << aliasname << "_->pt(), &e);" << endl;
	  headerf << "\t\t\t\t" << "if (not isfinite(" << aliasname << "_->pt()) || e > 30) {" << endl;
	  headerf << "\t\t\t\t\t" << "printf(\"branch " << Form("%s_branch",aliasname.Data()) 
		  << " contains a bad float: %f\\n\", " << aliasname << "_->pt());" << endl 
		  << "\t\t\t\t\t" << "exit(1);"
		  << endl;
	} else {
	  headerf << "\t\t\t\t" << "frexp(" << aliasname << "_.pt(), &e);" << endl;
	  headerf << "\t\t\t\t" << "if (not isfinite(" << aliasname << "_.pt()) || e > 30) {" << endl;
	  headerf << "\t\t\t\t\t" << "printf(\"branch " << Form("%s_branch",aliasname.Data()) 
		  << " contains a bad float: %f\\n\", " << aliasname << "_.pt());" << endl 
		  << "\t\t\t\t\t" << "exit(1);"
		  << endl;
	}
	headerf << "\t\t\t\t}" << endl;
      }
      headerf << "\t\t\t\t#endif // #ifdef PARANOIA" << endl;
    }
    headerf << "\t\t\t" << "} else { " << endl;
    headerf << "\t\t\t\t" << "printf(\"branch " << Form("%s_branch",aliasname.Data()) 
	    << " does not exist!\\n\");" << endl;
    headerf << "\t\t\t\t" << "exit(1);" << endl << "\t\t\t}" << endl;
    headerf << "\t\t\t" << Form("%s_isLoaded",aliasname.Data()) << " = true;" << endl;
    headerf << "\t\t" << "}" << endl;
    if(isSkimmedNtuple) {
      headerf << "\t\t" << "return *" << aliasname << "_;" << endl << "\t}" << endl;
    }
    else {
      headerf << "\t\t" << "return " << aliasname << "_;" << endl << "\t}" << endl;
    }
  }

  bool haveHLTInfo = false;
  bool haveL1Info  = false;
  bool haveHLT8E29Info = false;
  for(int i = 0; i < aliasarray->GetSize(); i++) {
    TString aliasname(aliasarray->At(i)->GetName());
    if(aliasname=="hlt_trigNames") 
      haveHLTInfo = true;
    if(aliasname=="l1_trigNames") 
      haveL1Info = true;
    if(aliasname=="hlt8e29_trigNames") 
      haveHLT8E29Info = true;
  }
   
  if(haveHLTInfo) {
    //functions to return whether or not trigger fired - HLT
    headerf << "\t" << "bool passHLTTrigger(TString trigName) {" << endl;
    headerf << "\t\t" << "int trigIndx;" << endl;
    headerf << "\t\t" << "vector<TString>::const_iterator begin_it = hlt_trigNames().begin();" << endl;
    headerf << "\t\t" << "vector<TString>::const_iterator end_it = hlt_trigNames().end();" << endl;
    headerf << "\t\t" << "vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);" << endl;
    headerf << "\t\t" << "if(found_it != end_it)" << endl;
    headerf << "\t\t\t" << "trigIndx = found_it - begin_it;" << endl;
    headerf << "\t\t" << "else {" << endl;
    headerf << "\t\t\t" << "cout << \"Cannot find Trigger \" << trigName << endl; " << endl;
    headerf << "\t\t\t" << "return 0;" << endl;
    headerf << "\t\t"   << "}" << endl << endl;
    //get the list of branches that hold the HLT bitmasks
    //store in a set 'cause its automatically sorted
    set<TString> s_HLTbitmasks;
    set<TString> s_L1bitmasks;
    for(int j = 0; j < aliasarray->GetSize(); j++) {
      TString aliasname(aliasarray->At(j)->GetName());
      TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
      TString classname = branch->GetClassName();
      if(aliasname.Contains("hlt_bits") && classname.Contains("int")) {
	s_HLTbitmasks.insert(aliasname);
      }
     
    }
    int i = 0;
    for(set<TString>::const_iterator s_it = s_HLTbitmasks.begin();
	s_it != s_HLTbitmasks.end(); s_it++, i++) {
      
      if(i==0) {
	headerf << "\t\t" << "if(trigIndx <= 31) {" << endl;
	headerf << "\t\t\t" << "unsigned int bitmask = 1;" << endl;
	headerf << "\t\t\t" << "bitmask <<= trigIndx;" << endl;	
	headerf << "\t\t\t" << "return " << *s_it << "() & bitmask;" << endl;
	headerf << "\t\t" << "}" << endl;
	continue;
      }
      headerf << "\t\t" << "if(trigIndx >= " << Form("%d && trigIndx <= %d", 32*i, 32*i+31) << ") {" << endl;
      headerf << "\t\t\t" << "unsigned int bitmask = 1;" << endl;
      headerf << "\t\t\t" << "bitmask <<= (trigIndx - " << Form("%d",32*i) << "); " << endl;	
      headerf << "\t\t\t" << "return " << *s_it << "() & bitmask;" << endl;
      headerf << "\t\t" << "}" << endl;
    }
    headerf << "\t" << "return 0;" << endl;
    headerf << "\t" << "}" << endl;
  }//if(haveHLTInfo) 

  if(haveHLT8E29Info) {
    //functions to return whether or not trigger fired - HLT
    headerf << "\t" << "bool passHLT8E29Trigger(TString trigName) {" << endl;
    headerf << "\t\t" << "int trigIndx;" << endl;
    headerf << "\t\t" << "vector<TString>::const_iterator begin_it = hlt8e29_trigNames().begin();" << endl;
    headerf << "\t\t" << "vector<TString>::const_iterator end_it = hlt8e29_trigNames().end();" << endl;
    headerf << "\t\t" << "vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);" << endl;
    headerf << "\t\t" << "if(found_it != end_it)" << endl;
    headerf << "\t\t\t" << "trigIndx = found_it - begin_it;" << endl;
    headerf << "\t\t" << "else {" << endl;
    headerf << "\t\t\t" << "cout << \"Cannot find Trigger \" << trigName << endl; " << endl;
    headerf << "\t\t\t" << "return 0;" << endl;
    headerf << "\t\t"   << "}" << endl << endl;
    //get the list of branches that hold the HLT bitmasks
    //store in a set 'cause its automatically sorted
    set<TString> s_HLTbitmasks;
    set<TString> s_L1bitmasks;
    for(int j = 0; j < aliasarray->GetSize(); j++) {
      TString aliasname(aliasarray->At(j)->GetName());
      TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
      TString classname = branch->GetClassName();
      if(aliasname.Contains("hlt8e29_bits") && classname.Contains("int")) {
	s_HLTbitmasks.insert(aliasname);
      }
     
    }
    int i = 0;
    for(set<TString>::const_iterator s_it = s_HLTbitmasks.begin();
	s_it != s_HLTbitmasks.end(); s_it++, i++) {
      
      if(i==0) {
	headerf << "\t\t" << "if(trigIndx <= 31) {" << endl;
	headerf << "\t\t\t" << "unsigned int bitmask = 1;" << endl;
	headerf << "\t\t\t" << "bitmask <<= trigIndx;" << endl;	
	headerf << "\t\t\t" << "return " << *s_it << "() & bitmask;" << endl;
	headerf << "\t\t" << "}" << endl;
	continue;
      }
      headerf << "\t\t" << "if(trigIndx >= " << Form("%d && trigIndx <= %d", 32*i, 32*i+31) << ") {" << endl;
      headerf << "\t\t\t" << "unsigned int bitmask = 1;" << endl;
      headerf << "\t\t\t" << "bitmask <<= (trigIndx - " << Form("%d",32*i) << "); " << endl;	
      headerf << "\t\t\t" << "return " << *s_it << "() & bitmask;" << endl;
      headerf << "\t\t" << "}" << endl;
    }
    headerf << "\t" << "return 0;" << endl;
    headerf << "\t" << "}" << endl;
  }//if(haveHLT8E29Info) 


  if(haveL1Info) {
    //functions to return whether or not trigger fired - L1
    headerf << "\t" << "bool passL1Trigger(TString trigName) {" << endl;
    headerf << "\t\t" << "int trigIndx;" << endl;
    headerf << "\t\t" << "vector<TString>::const_iterator begin_it = l1_trigNames().begin();" << endl;
    headerf << "\t\t" << "vector<TString>::const_iterator end_it = l1_trigNames().end();" << endl;
    headerf << "\t\t" << "vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);" << endl;
    headerf << "\t\t" << "if(found_it != end_it)" << endl;
    headerf << "\t\t\t" << "trigIndx = found_it - begin_it;" << endl;
    headerf << "\t\t" << "else {" << endl;
    headerf << "\t\t\t" << "cout << \"Cannot find Trigger \" << trigName << endl; " << endl;
    headerf << "\t\t\t" << "return 0;" << endl;
    headerf << "\t\t"   << "}" << endl << endl;
    //get the list of branches that hold the L1 bitmasks
    //store in a set 'cause its automatically sorted
    set<TString> s_L1bitmasks;
    for(int j = 0; j < aliasarray->GetSize(); j++) {
      TString aliasname(aliasarray->At(j)->GetName());
      TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
      TString classname = branch->GetClassName();
      if(aliasname.Contains("l1_bits") && classname.Contains("int")) {
	s_L1bitmasks.insert(aliasname);
      }
     
    }
    int i = 0;
    for(set<TString>::const_iterator s_it = s_L1bitmasks.begin();
	s_it != s_L1bitmasks.end(); s_it++, i++) {
      
      if(i==0) {
	headerf << "\t\t" << "if(trigIndx <= 31) {" << endl;
	headerf << "\t\t\t" << "unsigned int bitmask = 1;" << endl;
	headerf << "\t\t\t" << "bitmask <<= trigIndx;" << endl;	
	headerf << "\t\t\t" << "return " << *s_it << "() & bitmask;" << endl;
	headerf << "\t\t" << "}" << endl;
	continue;
      }
      headerf << "\t\t" << "if(trigIndx >= " << Form("%d && trigIndx <= %d", 32*i, 32*i+31) << ") {" << endl;
      headerf << "\t\t\t" << "unsigned int bitmask = 1;" << endl;
      headerf << "\t\t\t" << "bitmask <<= (trigIndx - " << Form("%d",32*i) << "); " << endl;	
      headerf << "\t\t\t" << "return " << *s_it << "() & bitmask;" << endl;
      headerf << "\t\t" << "}" << endl;
    }
    headerf << "\t" << "return 0;" << endl;
    headerf << "\t" << "}" << endl;
  }//if(haveL1Info)
    
  headerf << "};" << endl << endl;
    
  headerf << "#ifndef __CINT__" << endl;
  headerf << "extern " << Classname << " cms2;" << endl;
  headerf << "#endif" << endl << endl;

  // Create namespace that can be used to access the extern'd cms2
  // object methods without having to type cms2. everywhere.
  // Does not include cms2.Init and cms2.GetEntry because I think
  // it is healthy to leave those methods as they are
  headerf << "namespace tas {" << endl;
  implf   << "namespace tas {" << endl;
  for (Int_t i = 0; i< aliasarray->GetSize(); i++) {
    TString aliasname(aliasarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    TString classname = branch->GetClassName();
    TString title = branch->GetTitle();
    if ( classname.Contains("vector") ) {
      if(classname.Contains("edm::Wrapper") ) {
	classname = classname(0,classname.Length()-2);
	classname.ReplaceAll("edm::Wrapper<","");
      }
      headerf << "\t" << classname << " &" << aliasname << "()";
      implf   << "\t" << classname << " &" << aliasname << "()";
    } else {
      if(classname.Contains("edm::Wrapper<") ) {
	classname = classname(0,classname.Length()-1);
	classname.ReplaceAll("edm::Wrapper<","");
      }
      if(classname != "" ) {
	headerf << "\t" << classname << " &" << aliasname << "()";
	implf   << "\t" << classname << " &" << aliasname << "()";
      } else {
	if(title.EndsWith("/i")){
	  headerf << "\tunsigned int &" << aliasname << "()";
	  implf   << "\tunsigned int &" << aliasname << "()";
	}
	if(title.EndsWith("/F")){
	  headerf << "\tfloat &" << aliasname << "()";
	  implf   << "\tfloat &" << aliasname << "()";
	}
	if(title.EndsWith("/I")){
	  headerf << "\tint &" << aliasname << "()";
	  implf   << "\tint &" << aliasname << "()";
	}
      }
    }
    headerf << ";" << endl;
    implf   << " { return cms2." << aliasname << "(); }" << endl;
  }
  if(haveHLTInfo) {
    //functions to return whether or not trigger fired - HLT
    headerf << "\t" << "bool passHLTTrigger(TString trigName);" << endl;
    implf   << "\t" << "bool passHLTTrigger(TString trigName) { return cms2.passHLTTrigger(trigName); }" << endl;
  }//if(haveHLTInfo) 
  if(haveHLT8E29Info) {
    //functions to return whether or not trigger fired - HLT
    headerf << "\t" << "bool passHLT8E29Trigger(TString trigName);" << endl;
    implf   << "\t" << "bool passHLT8E29Trigger(TString trigName) { return cms2.passHLT8E29Trigger(trigName); }" << endl;
  }//if(haveHLT8E29Info) 
  if(haveL1Info) {
    //functions to return whether or not trigger fired - L1
    headerf << "\t" << "bool passL1Trigger(TString trigName);" << endl;
    implf   << "\t" << "bool passL1Trigger(TString trigName) { return cms2.passL1Trigger(trigName); }" << endl;
  }//if(haveL1Info)
    
}
  
  
// void makeSrcFile(TFile *f, bool paranoid, std::string Classname, std::string branchNamesFile) {
void makeSrcFile(std::string Classname, std::string branchNamesFile) {

  // TTree *ev = (TTree*)f->Get("Events");
  
  codef << "/* Usage:" << endl;
  codef << "   root [0] .L ScanChain.C++" << endl;
  codef << "   root [1] TFile *_file0 = TFile::Open(\"merged_ntuple.root\")" << endl;
  codef << "   root [2] TChain *chain = new TChain(\"Events\")" << endl;
  codef << "   root [3] chain->Add(\"merged_ntuple.root\")" << endl;
  codef << endl;
  codef << "   There are several places where one may create " << Classname << " cms2" << endl;
  codef << "   It can be done here (in a doAll.C script), i.e.:" << endl;
  codef << endl;
  codef << "   root [4] " << Classname << " cms2 " << endl;
  codef << endl;
  codef << "   It can be done in the source as is done below, or it can be" << endl;
  codef << "   ascertained by including CORE/CMS2.cc as is commented out" << endl;
  codef << "   below.  They are all the same, and everything will work so" << endl;
  codef << "   long as it is created somewhere globally." << endl;
  codef << endl;
  codef << "   root [5] ScanChain(chain)" << endl;
  codef << "*/" << endl;
  codef << "#include <iostream>" << endl;
  codef << "#include <vector>" << endl;
  codef << "" << endl;
  codef << "#include \"TChain.h\"" << endl;
  codef << "#include \"TFile.h\"" << endl;
  codef << "#include \"TDirectory.h\"" << endl;
  codef << "#include \"TROOT.h\"" << endl;
  codef << "" << endl;
  codef << "#include \"" + Classname+".h\"" << endl;
  if(branchNamesFile!="")
    codef << "#include \"branches.h\"" << endl;
  codef << Classname << " cms2;" << endl;
  codef << "/*" << endl;
  
  codef << "#include \"CORE/CMS2.cc\"" << endl;
  codef << "#include \"CORE/selections.cc\"" << endl;
  codef << "#include \"CORE/utilities.cc\"" << endl;
  codef << "*/" << endl;
  codef << endl;
  codef << "using namespace tas;" << endl;
  codef << endl;
  
  codef << "int ScanChain( TChain* chain, int nEvents = -1, std::string skimFilePrefix=\"\") {" << endl;
  codef << "" << endl;
  codef << "  TObjArray *listOfFiles = chain->GetListOfFiles();" << endl;
  codef << "" << endl;
  codef << "  unsigned int nEventsChain=0;" << endl;
  codef << "  if(nEvents==-1) " << endl << "    nEvents = chain->GetEntries();" << endl;
  codef << "  nEventsChain = nEvents;" << endl;
  
  codef << "  unsigned int nEventsTotal = 0;" << endl;
  if(branchNamesFile!="")
    codef << "  InitSkimmedTree(skimFilePrefix);" << endl;
  codef << "  TDirectory *rootdir = gDirectory->GetDirectory(\"Rint:\");" << endl << endl;
  codef << "  TH1F *samplehisto = new TH1F(\"samplehisto\", \"Example histogram\", 200,0,200);" << endl;
  codef << "  samplehisto->SetDirectory(rootdir);" << endl;

  codef << "  // file loop" << endl;
  codef << "  TIter fileIter(listOfFiles);" << endl;
  codef << "  TFile *currentFile = 0;" << endl;
  codef << "  while ( currentFile = (TFile*)fileIter.Next() ) {" << endl;
  codef << "    TFile f(currentFile->GetTitle());" << endl;
  codef << "    TTree *tree = (TTree*)f.Get(\"Events\");" << endl;
  codef << "    cms2.Init(tree);" << endl;
  codef << "    " << endl;
  codef << "    //Event Loop" << endl;
  codef << "    unsigned int nEvents = tree->GetEntries();" << endl;
  codef << "    for( unsigned int event = 0; event < nEvents; ++event) {" << endl;
  codef << "      cms2.GetEntry(event);" << endl;
  codef << "      ++nEventsTotal;" << endl;
  //codef << "      std::cout << \"els size: \" << els_p4().size() << \" \";" << endl;
  //codef << "      std::cout << \"mus size: \" << mus_p4().size() << std::endl;" << endl << endl;
  //codef << "      for (unsigned int mus = 0; " << endl;
  //codef << "           mus < mus_p4().size(); mus++) " << endl << endl;
  //codef << "         samplehisto->Fill(mus_p4().at(mus).Pt());" << endl << endl;
  //codef << "      for (unsigned int hyp = 0;" << endl;
  //codef << "           hyp < hyp_jets_p4().size();" << endl;
  //codef << "           ++hyp) {" << endl;
  //codef << "        std::cout << \"hyp: \" << hyp << \"jet corrections:\";" << endl;
  //codef << "        for ( unsigned int jet = 0;" << endl;
  //codef << "              jet < hyp_jets_p4()[hyp].size();" << endl;
  //codef << "              ++jet ) {" << endl;
  //codef << "          std::cout << \" \" << hyp_jets_p4()[hyp][jet].pt();" << endl;
  //codef << "        }" << endl;
  //codef << "        std::cout << endl;" << endl;
  //codef << "      }" << endl;
  //    codef << "      if ( hyp_jets_p4().size() == 0 ) {" << endl;
  //codef << "        std::cout << \"no hypothesis!\" << std::endl;" << endl;
  //codef << "      }" << endl;
  codef << "    }" << endl;
  codef << "  }" << endl;
  codef << "" << endl;
  codef << "  if ( nEventsChain != nEventsTotal ) {" << endl;
  codef << "    std::cout << \"ERROR: number of events from files is not equal to total number of events\" << std::endl;" << endl;
  codef << "  }" << endl;
  codef << "" << endl;
  if(branchNamesFile!="") {
    codef << "  outFile_->cd();" << endl;
    codef << "  outTree_->Write();" << endl;
    codef << "  outFile_->Close();" << endl;
  }
  codef << "  samplehisto->Draw();" << endl;
  codef << "  return 0;" << endl;
  codef << "}" << endl;
  
}

  
      
void makeBranchFile(std::string branchNamesFile) {
  
  ifstream branchesF(branchNamesFile.c_str());
  vector<TString> v_datatypes;
  vector<TString> v_varNames;
  while(!branchesF.eof()) {
    string temp;
    getline(branchesF, temp);
    TString line(temp);
    // do not use commented lines
    if(line.BeginsWith("//"))
      continue;
    vector<TString> v_line;
    TIter objIt((TObjArray*)line.Tokenize(" "));
    while(TObject *obj = (TObject*)objIt.Next()) {
      if(obj==NULL) 
	continue;
      v_line.push_back(obj->GetName());
    }
    
    if(v_line.size() == 0)
      continue;
    TString varName;
    //loop over v_line until you get to the first element thats not a space
    for(vector<TString>::iterator it = v_line.begin();
	it != v_line.end(); it++) {
      if( *it != " " ) {
	varName = *it;
      }
    }
    
    v_varNames.push_back(varName.Strip());  //paranoid....strips trailing spaces
    TString datatype("");
    for(unsigned int i = 0; i < v_line.size()-1; i++) {
      TString temp = v_line[i];
      if(temp.Contains("vector") && !temp.Contains("std::")) 
	temp.ReplaceAll("vector", "std::vector");
      if(temp.Contains("LorentzVector") && !temp.Contains("ROOT::Math::LorentzVector"))
	temp.ReplaceAll("LorentzVector", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >");
      temp.ReplaceAll(">>", "> >");
      temp.ReplaceAll(">>>", "> > >");
      if(i!=0)
	datatype = datatype+" " + temp;
      else
	datatype = datatype+temp;
    }
    v_datatypes.push_back(datatype);
    
  } 
  branchfile.open("branches.h");
  branchfile << "#ifndef BRANCHES_H" << endl << "#define BRANCHES_H" << endl;
  branchfile << "#include <vector>" << endl;
  branchfile << "#include \"TFile.h\"" << endl;
  branchfile << "#include \"TTree.h\"" << endl << endl << endl << endl;

  for(unsigned int i = 0; i < v_datatypes.size(); i++) {
    TString temp(v_varNames.at(i));
    if(v_datatypes.at(i).Contains("vector") || v_datatypes.at(i).Contains("LorentzVector")) {
      branchfile << v_datatypes.at(i) << " *" 
		 << temp.ReplaceAll("_","")+"_;" << endl;
      continue;
    }
    branchfile << v_datatypes.at(i) << " " 
	       << temp.ReplaceAll("_","")+"_;" << endl;
  }


  branchfile << "TFile *outFile_;" << endl;
  branchfile << "TTree *outTree_;" << endl;

  //now declare the branches and set aliases in the InitSkimmedTree function
  branchfile << "void InitSkimmedTree(std::string skimFilePrefix=\"\") {\n\n";
  branchfile << "   if(skimFilePrefix != \"\")" << endl;
  branchfile << "      outFile_ = TFile::Open(string(skimFilePrefix + \"_skimmednTuple.root\").c_str(),\"RECREATE\");\n";
  branchfile << "   else outFile_ = TFile::Open(\"skimmednTuple.root\", \"RECREATE\");\n";
  branchfile << "   outFile_->cd();" << endl;
  branchfile << "   outTree_ = new TTree(\"Events\", \"\");\n\n";
  branchfile << "   //book the branches\n";
  for(unsigned int i = 0; i < v_datatypes.size(); i++) {
    TString varName = v_varNames[i];
    varName = varName.ReplaceAll("_", "") + "_";
    TString varType = v_datatypes[i];
    if(varType.BeginsWith("std::vector") 
       || varType.BeginsWith("ROOT::Math") ) {
      TString prefix = "";
      if(varType.BeginsWith("std::vector<float>") )
        prefix = "floats";
      if(varType.BeginsWith("std::vector<int>") ) 
        prefix = "ints";
      if(varType.BeginsWith("std::vector<unsigned int>") ) 
        prefix = "uints";
      if(varType.BeginsWith("ROOT::Math::LorentzVector<") )
        prefix = "floatROOTMathPxPyPzE4DROOTMathLorentzVector";
      if(varType.BeginsWith("std::vector<ROOT::Math::LorentzVector<") )
        prefix = "floatROOTMathPxPyPzE4DROOTMathLorentzVectors";
      if(varType.Contains("std::vector<std::vector<ROOT::Math::LorentzVector<") )
        prefix = "floatROOTMathPxPyPzE4DROOTMathLorentzVectorss";
      branchfile << "   outTree_->Branch(\"" << prefix + "_" +varName << "\",   \"" 
		 << varType << "\",   &" << varName << ");" << endl;
      branchfile << "   outTree_->SetAlias(\"" << v_varNames[i] << "\",   " 
		 << "\"" << prefix+"_"+varName << "\");" << endl;
      continue;
    }
    if(varType=="float" || varType == "Float_t") {
      branchfile << "   outTree_->Branch(\"" << varName << "\",   &" << varName;
      branchfile << ",   \"" << varName + "/F\");" << endl;
      branchfile << "   outTree_->SetAlias(\"" << v_varNames[i] << "\",   " 
		 << "\"" << varName << "\");" << endl;
      continue;
    }
    if(varType=="unsigned int" || varType == "UInt_t") {
      branchfile << "   outTree_->Branch(\"" << varName << "\",   &" << varName;
      branchfile << ",   \"" << varName + "/i\");" << endl;
      branchfile << "   outTree_->SetAlias(\"" << v_varNames[i] << "\",   " 
		 << "\"" << varName << "\");" << endl;
      continue;
    }
    if(varType=="int" || varType == "Int_t") {
      branchfile << "   outTree_->Branch(\"" << varName << "\",   &" << varName;
      branchfile << ",   \"" << varName + "/I\");" << endl;
      branchfile << "   outTree_->SetAlias(\"" << v_varNames[i] << "\",   " 
		 << "\"" << varName << "\");" << endl;
      continue;
    }
  }
  branchfile << "} " <<  endl;
  branchfile << "#endif" << endl;
  
  branchfile.close();
}


void makeCMS2ClassFiles (std::string fname, bool paranoid = true, 
			 std::string branchNamesFile="", std::string className = "CMS2") {

  using namespace std;
  
  TFile *f = TFile::Open( fname.c_str() );
  if(f->IsZombie()) { 
    cout << "File is not a valid root file, or root file is corruped" << endl;
    cout << "Exiting..." << endl;
  }
  //class is CMS2 by default
  //std::string Classname = className=="" ? "CMS2" : className;
  
  headerf.open((className+".h").c_str());
  implf.open((className+".cc").c_str());
  codef.open("ScanChain.C");
  
  implf << "#include \"" << className+".h" << "\"\nCMS2 cms2;" << endl;
  makeHeaderFile(f, paranoid, className);
  makeSrcFile(className, branchNamesFile);
  if(branchNamesFile!="")
    makeBranchFile(branchNamesFile);
  
  implf << "}" << endl;
  implf.close();
  headerf << "}" << endl;
  headerf << "#endif" << endl;
  headerf.close();

  
  codef.close();
  f->Close();
}
  
