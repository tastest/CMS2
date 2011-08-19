/*
  Skim ntuples by requiring 2 isolated leptons above 20 GeV
  Usage:
     .L SkimNtuples.C+
     chain = new TChain("Events")
     chain->Add("~/scratch/TnS/CompHepZbbll_1_6_7/ntuple_*")
     SkimNtuples(chain,true,"/tmp/skim.root")
*/
#include <iostream>
#include <vector>
#include <set>
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TEventList.h"
#include "../../CORE/CMS2.cc"
#include "../../CORE/utilities.cc"
#include "../../CORE/selections.cc"
#include "TDatabasePDG.h"
#include "TProfile.h"

void SkimNtuples( TChain* chain, bool filterEvents = true, const char* outputFileName = "skim.root" )
{
   long long max_tree_size = 200000000000LL;
   TTree::SetMaxTreeSize(max_tree_size);
   
   chain->GetEntry(0);
   
   unsigned int nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
   unsigned int nEventsTotal = 0;
   
   // define the clone
   TFile* outputFile = new TFile(outputFileName,"recreate");
   TChain* outputChain = (TChain*)chain->CloneTree(0);
   TTree* outputTree = outputChain->GetTree();
   
   int i_permille_old = 0;
   TObjArray *listOfFiles = chain->GetListOfFiles();
   TIter fileIter(listOfFiles);
   while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
      printf("Processing file: %s\n",currentFile->GetTitle());
      TFile *f = TFile::Open(currentFile->GetTitle()); 
      assert(f);
      TTree *tree = (TTree*)f->Get("Events");
      cms2.Init(tree);
      
      //Event Loop
      unsigned int nEvents = tree->GetEntries();
      for( unsigned int event = 0; event < nEvents; ++event) {
	 ++nEventsTotal;
	 cms2.GetEntry(event);  // get entries for Event number event from branches of TTree tree
	 
	 int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
	 if (i_permille != i_permille_old) {
	    // xterm magic from L. Vacavant and A. Cerri
	    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
		   "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	    fflush(stdout);
	    i_permille_old = i_permille;
	 }
	 
	 // loop over hypothesis candidates
	 unsigned int nHyps = cms2.hyp_type().size();
	 bool goodEvent = false;
	 if ( ! filterEvents )
	   goodEvent = true;
	 else {
	    for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
	       // Cut on lepton Pt
	       if (cms2.hyp_lt_p4()[i_hyp].pt() < 20.0) continue;
	       if (cms2.hyp_ll_p4()[i_hyp].pt() < 20.0) continue;
	       if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) continue;
	       if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) continue;
	       
	       // Electron quality cuts, including isolation (trk only)
	       if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_lt_index()[i_hyp],false) ) continue;
	       if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !goodElectronIsolated(cms2.hyp_ll_index()[i_hyp],false) ) continue;
	    
	       goodEvent = true;
	       break;
	    }
	 }
	 
	 if ( goodEvent ) {
	    chain->GetEntry(chain->GetChainEntryNumber(event));
	    outputTree->Fill();
	 }
      }
      f->Close();
   }
   printf("Done\n");
   outputFile->cd();
   outputTree->Write();
   outputFile->Close();
}
