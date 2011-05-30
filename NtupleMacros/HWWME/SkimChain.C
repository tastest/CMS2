/* 
=== Description
   root [0] .L SkimChain.C++
   root [1] TChain *chain = new TChain("Events")
   root [2] chain->Add("merged_ntuple.root")
   root [3] SkimChain(chain)
*/

// C++
#include "TPRegexp.h"
#include "TRegexp.h"
#include "TSystem.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TDatabasePDG.h"



// CMS2
#include "CORE/CMS2.cc"

using namespace tas;

struct FileEntry{
  std::string inputDir;
  std::string outputDir;
  std::string fileName;
  unsigned int nIn;
  unsigned int nOut;
  bool operator<(const FileEntry& rhs) const { return inputDir < rhs.inputDir; }
};


// select only the events with colored parton pT > 10 GeV and |eta| < 2.5
bool passedSkimSelection()
{
  // skipping the first 6 genparticles which are the incoming partons
  for (unsigned int i = 6;  i < cms2.genps_id().size(); ++i) {
    // selecting status = 3 colored objects
    if (cms2.genps_status().at(i) != 3 || (cms2.genps_id().at(i) != 21 && TMath::Abs(cms2.genps_id().at(i)) > 8 )) continue;
    // apply acceptance cuts, which must be looser than the FO definition
    if (cms2.genps_p4().at(i).pt() >= 10 && TMath::Abs(cms2.genps_p4().at(i).eta()) <= 2.5 ) return true;
  }
  return false;
}


void SkimChain(TChain* chain, const TString DataDir, TString SkimDir, TString skimname){

  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TPRegexp fileNameMatchPattern("(.*?)/([^/]*)$");
  TPRegexp dataDirMatchPattern(DataDir);
  
  unsigned int nEventsTotal = 0;
  unsigned int nEventsSelectedTotal = 0;

  std::vector<FileEntry> files;
  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
    FileEntry file;
    TObjArray* matches = fileNameMatchPattern.MatchS(currentFile->GetTitle());
    assert(matches);
    assert(matches && matches->GetLast()==2);
    TString inputFileName(currentFile->GetTitle());
    TString fileName(((TObjString*)matches->At(2))->GetString());
    file.fileName = fileName.Data();
    TString outputDirectory(((TObjString*)matches->At(1))->GetString());
    file.inputDir = outputDirectory.Data();
    dataDirMatchPattern.Substitute(outputDirectory,SkimDir);
    outputDirectory += "/" + skimname;
    file.outputDir = outputDirectory.Data();
    // make output directory if it doesn't exist yet
    if ( gSystem->AccessPathName(outputDirectory.Data()) ){
      gSystem->mkdir(outputDirectory.Data(),true);
      assert( !gSystem->AccessPathName(outputDirectory.Data()) );
    }
    TString outputFileName(outputDirectory+"/"+fileName);
    cout << "Skimming " << inputFileName << " -> " << outputFileName << endl;

    TFile *output = TFile::Open(outputFileName.Data(), "RECREATE");
    assert(output);
    TFile *input = TFile::Open(inputFileName.Data());
    assert(input);
    TTree *tree = (TTree*)input->Get("Events");
    TTree *newtree = tree->CloneTree(0);
    newtree->SetDirectory(output);
    
    cms2.Init(newtree);
    cms2.Init(tree);
    
    // Event Loop
    const unsigned int nEvents = tree->GetEntries();
    unsigned int nEventsSelected = 0;
    int i_permille_old = 0;
    for (unsigned int event = 0; event < nEvents; ++event, ++nEventsTotal) {
      int i_permille = (int)floor(10000 * event / float(nEvents));
      if (i_permille != i_permille_old) {
	// xterm magic from L. Vacavant and A. Cerri
	if (isatty(1)) {
	  printf("\015\033[32m ---> \033[1m\033[31m%5.2f%%"
		 "\033[0m\033[32m <---\033[0m\015", i_permille/100.);
	  fflush(stdout);
	}
	i_permille_old = i_permille;
      }
      cms2.GetEntry(event);
      //set condition to skip event
      if ( not passedSkimSelection() ) continue;
            
      ++nEventsSelected;
      cms2.LoadAllBranches();
      // fill the new tree
      newtree->Fill();
    }
    file.nIn = nEvents;
    file.nOut = nEventsSelected;
    files.push_back(file);
    nEventsSelectedTotal += nEventsSelected;
    output->cd();
    newtree->Write();
    output->Close();
    input->Close();
  }
  cout << Form("Processed events: %u, \tselected: %u\n",nEventsTotal,nEventsSelectedTotal) << endl;
}





void ProcessSample(std::vector<std::string> file_patterns, TString DataDir, TString SkimDir, TString skimname) {
  TChain *tchain = new TChain("Events");
  for ( std::vector<std::string>::const_iterator pattern = file_patterns.begin();
	pattern != file_patterns.end(); ++pattern )
    tchain->Add(pattern->c_str());
  
  SkimChain(tchain, DataDir, SkimDir, skimname);

}

void ProcessSample(std::string file_pattern, TString DataDir, TString SkimDir, TString skimname ) {
  std::vector<std::string> vec;
  vec.push_back(file_pattern);
  ProcessSample(vec, DataDir, SkimDir, skimname);
}
