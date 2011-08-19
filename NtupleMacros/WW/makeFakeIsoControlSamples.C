//now make the source file
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TStopwatch.h"
#include <algorithm>
#include <set>
#include "TCanvas.h"
#include "TRegexp.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

using namespace std;

#ifndef __CINT__
#include "../CORE/CMS2.cc"
#include "../CORE/utilities.cc"
#include "../CORE/selections.cc"
#endif


// fake type: 
//  - electron
//  - muon 
TH1F* ScanChain( TChain* chain, const char* sample, const char* type ) {
  
  unsigned int nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
  unsigned int nEventsTotal = 0;

  TH1F* h = new TH1F(Form("h_%s_%s",type,sample),Form("Fake %s isolation control sample (%s)"),100,0.,1.);
  h->Sumw2();
  h->SetDirectory(0);

  int i_permille_old = 0;
  // file loop
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
    TFile *f = TFile::Open(currentFile->GetTitle()); 
    assert(f);
    TTree *tree = (TTree*)f->Get("Events");
    
    cms2.Init(tree);  // set branch addresses for TTree tree

    TStopwatch t;
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);  // get entries for Event number event from branches of TTree tree
      ++nEventsTotal;
      
      int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
      if (i_permille != i_permille_old) {
	// xterm magic from L. Vacavant and A. Cerri
	printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
	       "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	fflush(stdout);
	i_permille_old = i_permille;
      }
      
      float weight = cms2.evt_scale1fb();

      if ( strcmp("electron",type)==0 ){
	// loop over electrons
	for ( unsigned int i=0; i<cms2.els_p4().size(); ++i ){
	  if (cms2.els_p4()[i].pt() < 20.0) continue;
	  if ( goodElectronWithoutIsolation(i) )
	    h->Fill(el_rel_iso(i,true),weight);
	}
      } else {
	// loop over muons
	for ( unsigned int i=0; i<cms2.mus_p4().size(); ++i ){
	  if (cms2.mus_p4()[i].pt() < 20.0) continue;
	  if ( goodMuonWithoutIsolation(i) )
	    h->Fill(mu_rel_iso(i),weight);
	}
      }
    }
    t.Stop();
    printf("Finished processing file: %s\n",currentFile->GetTitle());
    printf("Real time: %u events / %f s = %e event/s\n", nEvents, 
	   t.RealTime(), nEvents / t.RealTime());
    printf("CPU time: %u events / %f s = %e event/s\n", nEvents, 
	   t.CpuTime(), nEvents / t.CpuTime());
    f->Close();
  }
  if ( nEventsChain != nEventsTotal ) {
    printf("ERROR: number of events from files (%d) is not equal to total number"
	   " of events (%d)\n", nEventsChain, nEventsTotal);
  }
  return h;
}

TH1F* makeFakeIsoControlSample(const char* files, const char* name, const char* type)
{
  // read dataset prefix
  string dataset;
  if ( ! gSystem->Getenv("CMS2_NTUPLE_LOCATION") ){
    cout << "ERROR: Dataset location is not set. Please set CMS2_NTUPLE_LOCATION." <<endl;
    return 0;
  }
  dataset = gSystem->Getenv("CMS2_NTUPLE_LOCATION");
  
  TChain *chain = new TChain("Events");
  chain->Add((dataset+files).c_str());
  if (chain->GetEntries() == 0 ){
    std::cout << "ERROR: chain is empty for sample: " << name << std::endl;
    return 0;
  }
  return ScanChain( chain, name, type );
}

void makeFakeIsoControlSamples()
{
  TH1F* h_qcd30_e = 
    makeFakeIsoControlSample("/cms2-V01-02-06/QCDpt30_v2/merged_ntuple*root","qcd30","electron");
  TH1F* h_qcd30_m =
    makeFakeIsoControlSample("/cms2-V01-02-06/QCDpt30_v2/merged_ntuple*root","qcd30","muon");

  TH1F* h_qcd80_e =
    makeFakeIsoControlSample("/cms2-V01-02-06/QCDpt80/merged_ntuple*root","qcd80","electron");
  TH1F* h_qcd80_m =
    makeFakeIsoControlSample("/cms2-V01-02-06/QCDpt80/merged_ntuple*root","qcd80","muon");

  TH1F* h_qcd170_e =
    makeFakeIsoControlSample("/cms2-V01-02-06/QCDpt170/merged_ntuple*root","qcd170","electron");
  TH1F* h_qcd170_m =
    makeFakeIsoControlSample("/cms2-V01-02-06/QCDpt170/merged_ntuple*root","qcd170","muon");
  
  TH1F* h_qcd300_e =
    makeFakeIsoControlSample("/cms2-V01-02-06/QCDpt300/merged_ntuple*root","qcd300","electron");
  TH1F* h_qcd300_m =
    makeFakeIsoControlSample("/cms2-V01-02-06/QCDpt300/merged_ntuple*root","qcd300","muon");

  if ( h_qcd30_e && h_qcd30_m && h_qcd80_e && h_qcd80_m && 
       h_qcd170_e && h_qcd170_m && h_qcd300_e && h_qcd300_m ) {
    
    TFile* f = TFile::Open("fakeIsoControlSamples.root","RECREATE");

    h_qcd30_e->Write();
    h_qcd30_m->Write();
    h_qcd80_e->Write();
    h_qcd80_m->Write();
    h_qcd170_e->Write();
    h_qcd170_m->Write();
    h_qcd300_e->Write();
    h_qcd300_m->Write();

    f->Close();
  } else {
    std::cout << "Failed to make fake control sample histograms. Abort" << std::endl;
    exit(1);
  }
}
