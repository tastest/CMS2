//example usage:
// > root
// root [0] .L NtupleMacros/Tools/PrintBranchSizes.C++
// ...
// root [1] PrintBranches("ntuple_nofilt_100evts_qcd30.root")



#include "TTree.h"
#include "TObjArray.h"
#include "TString.h"
#include "TBranch.h"
#include "TFile.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <iostream>
#include <iomanip>


using namespace std;

void PrintBranchSizes(TTree* tree);

//call function below with name of ntuple file as arg
void PrintSizes(string input) {
  TFile *infile = TFile::Open(input.c_str());
  TTree *tree = (TTree*)infile->Get("Events");
  PrintBranchSizes(tree);
}


void PrintBranchSizes(TTree* tree) {

  TList *t_list =  tree->GetListOfAliases();

  map<string, float> m_sortByMakerName;
  map<float, string> m_sortBySize;
  
  //float treeSizeZipBytes = tree->GetZipBytes("*");
  float nEvents  = tree->GetEntries();
  int maxAliasNameSize = 0;
  float sumBranchSize = 0.;
  
  for( int i = 0; i < t_list->GetSize(); i++) {

    string aliasname = t_list->At(i)->GetName();
    if(aliasname.size() > maxAliasNameSize)
      maxAliasNameSize = aliasname.size();
    TBranch *b = tree->GetBranch(tree->GetAlias(aliasname.c_str()));
    //only want branches made by us
    float size = b->GetZipBytes("*");
    sumBranchSize = sumBranchSize + size;
    m_sortBySize[-1*size]           = string(aliasname);

    TObjArray *temp = TString(b->GetName()).Tokenize("_");
    TString strppedbName = TString(temp->At(1)->GetName()) + "_" + TString(temp->At(2)->GetName());
    m_sortByMakerName[string(strppedbName.Data())] = size; 
  }

  string heading = "Branch Name";
  int padding = maxAliasNameSize - 11;
  for(unsigned int i = 0; i < padding; i++) 
    heading = heading + " ";
  cout << heading << "     Size      Percentage " << endl;

  
  for(map<float, string>::iterator it = m_sortBySize.begin();
      it != m_sortBySize.end(); it++) {

    padding = maxAliasNameSize - string(it->second).size();
    string aliasname = it->second;
    cout << setiosflags(ios::fixed);
    for(unsigned int i = 0; i < padding; i++) 
      aliasname = aliasname + " ";
    cout << aliasname << " " << -1*it->first/1000. << " kB (" 
	 <<-100*(it->first)/sumBranchSize << "%)" << endl;
        
  }
  
  cout << "Total size of branches: " << sumBranchSize/(1000.*1000.) << " MB" << endl;
 
  //now report by maker
  TString prevMakerName(m_sortByMakerName.begin()->first);
  prevMakerName = TString((prevMakerName.Tokenize("_"))->At(0)->GetName());
  float totMakerSize = 0.;
  map<float, string> m_sizeByMaker;
  float maxMakerSize = 0;
  for(map<string, float>::iterator it = m_sortByMakerName.begin();
      it != m_sortByMakerName.end(); it++) {

    TString makerName(it->first);
    makerName = TString(makerName.Tokenize("_")->At(0)->GetName());
    if(makerName.Length() > maxMakerSize)
      maxMakerSize = makerName.Length();
    if(makerName==prevMakerName) 
      totMakerSize = totMakerSize + it->second;
    else {
      m_sizeByMaker[-1*totMakerSize] = prevMakerName;
      totMakerSize = it->second;
    }
    prevMakerName = makerName;
  }

  
  float sumMakerSize = 0.0;
  cout << "Organized according to Producer:" << endl;
  for(map<float, string>::iterator it = m_sizeByMaker.begin(); 
      it != m_sizeByMaker.end(); it++) {

    float padding = maxMakerSize - it->second.size();
    string temp = it->second;
    sumMakerSize = sumMakerSize + it->first;
    for(unsigned int i = 0; i < padding; i++) 
      temp = temp + " ";
    cout << temp << " " << -1*it->first/1000. << " kB (" 
	 << -100.*it->first/(sumBranchSize) << "%)" <<endl;
  }

  cout << "Total size of Makers: " << -1*sumMakerSize/(1000.*1000.) << " MB" << endl;

}

    
  
  
