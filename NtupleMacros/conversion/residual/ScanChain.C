/* Usage:
   root [0] .L ScanChain.C++
   root [1] TFile *_file0 = TFile::Open("merged_ntuple.root")
   root [2] TChain *chain = new TChain("Events")
   root [3] chain->Add("merged_ntuple.root")

   There are several places where one may create CMS2 cms2
   It can be done here (in a doAll.C script), i.e.:

   root [4] CMS2 cms2 

   It can be done in the source as is done below, or it can be
   ascertained by including CORE/CMS2.cc as is commented out
   below.  They are all the same, and everything will work so
   long as it is created somewhere globally.

   root [5] ScanChain(chain)
*/
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "CMS2.h"
//#include "branches.h"

CMS2 cms2;

//#include "CORE/CMS2.cc"
// #include "CORE/selections.cc"
// #include "CORE/utilities.cc"
// #include "Tools/tools.cc"

using namespace tas;


int ScanChain( TChain* chain, int nEvents = -1, std::string skimFilePrefix="") {

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=0;
  if(nEvents==-1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  // InitSkimmedTree(skimFilePrefix);
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();



 
  TH1F *trks_residualX_nminus1 = new TH1F ("residualX_nminus1",
				    "residualX_nminus1",
				    200,
				    -10.0,
				    10.0
				    );
				    
 
 TH1F *trks_residualX_nplus1 = new TH1F("residualX_nplus1",
				    "residualX_nplus1",
				    200,
				    -10.0,
				    10.0
				    );
				    
  

 

  
  // samplehisto->SetDirectory(rootdir);
  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
   
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);
     
   
      int nels =  cms2.els_trk_p4().size();
      int ntrks =  cms2.trks_trk_p4().size();
      
   //    for(int i=0; i<nels; i++){


// 	int trk_index = els_trkidx().at(i);


// 	if(trk_index >= 0 ){
// 	  int nhits = cms2.trks_residualX()[trk_index].size();
// 	  for (int i_hit = 0; i_hit <nhits; i_hit++){
// 	    if(trks_hit_type()[trk_index].at(i_hit) ==0){
	    
// 	      if(i_hit+1 <nhits && trks_hit_type()[trk_index].at(i_hit+1) !=0){
// 		trks_residualX_nplus1->Fill( trks_residualX()[trk_index].at(i_hit+1));
		
// 	      }
	    
// 	      if(i_hit-1>=0 && trks_hit_type()[trk_index].at(i_hit-1) !=0 ){
// 		trks_residualX_nminus1 ->Fill(trks_residualX()[trk_index].at(i_hit-1));
		
// 	      }
// 	    }
// 	  }
// 	}

	
//       }




      for(int trk_index =0; trk_index <ntrks;trk_index++){

	int nhits = cms2.trks_residualX()[trk_index].size();
	for (int i_hit = 0; i_hit <nhits; i_hit++){
	  if(trks_hit_type()[trk_index].at(i_hit) ==0){
	    
	    if(i_hit+1 <nhits && trks_hit_type()[trk_index].at(i_hit+1) !=0){
	      trks_residualX_nplus1->Fill( trks_residualX()[trk_index].at(i_hit+1));
	      
	    }
	    
	    if(i_hit-1>=0 && trks_hit_type()[trk_index].at(i_hit-1) !=0 ){
	      trks_residualX_nminus1 ->Fill(trks_residualX()[trk_index].at(i_hit-1));
	      
	    }
	  }
	}
	
	
	
      }


   
      // outTree_->Fill();

      ++nEventsTotal;
      if(nEventsTotal%10000 ==0)std::cout << "number of events processed " << nEventsTotal<<std::endl;
    }
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    std::cout << "number of events processed " << nEventsTotal<<std::endl;
  }
 
  return 0;
}
