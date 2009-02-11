/* Usage:
   root[0] .L CMS2V010201_looper.C++
   root [1] TFile *_file0 = TFile::Open("ntuple_file.root")
   root [2] TChain *chain = new TChain("Events")
   root [3] chain->Add("ntuple_file.root")
   root [4] CMS2V010201_looper a 
   root [5] a.ScanChain(chain) // will give the same results
*/
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "CMS2V010201_looper.h"


int CMS2V010201_looper::ScanChain( TChain* chain, int nEvents) {

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=0;
  if(nEvents==-1) 
     nEvents = chain->GetEntries();
  else nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  TH1F *samplehisto = new TH1F("samplehisto", "Example histogram", 200,0,200);
  samplehisto->SetDirectory(rootdir);
  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents; ++event) {
      GetEntry(event);
      ++nEventsTotal;
      std::cout << "els size: " << els_p4().size() << " ";
      std::cout << "mus size: " << mus_p4().size() << std::endl;
      for (unsigned int mus = 0; 
           mus < mus_p4().size(); mus++) 

         samplehisto->Fill(mus_p4().at(mus).Pt());

      for (unsigned int hyp = 0;
           hyp < hyp_jets_p4().size();
           ++hyp) {
        std::cout << "hyp: " << hyp << "jet corrections:";
        for ( unsigned int jet = 0;
              jet < hyp_jets_p4()[hyp].size();
              ++jet ) {
          std::cout << " " << hyp_jets_p4()[hyp][jet].pt();
        }
        std::cout << endl;
      }
      if ( hyp_jets_p4().size() == 0 ) {
        std::cout << "no hypothesis!" << std::endl;
      }
    }
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  samplehisto->Draw();
  return 0;
}
