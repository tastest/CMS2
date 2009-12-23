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
#include <algorithm>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH2F.h"
#include "CMS2.h"
#include "TProfile.h"

#define MAXMET 19.8


CMS2 cms2;

using namespace tas;
bool isGoodTrk( unsigned int i );
//bool isGoodEvent();
bool passesTrigger(bool runningonGEN);
//void getMETQuantities(const float metin, const float metPhiin, 
//		      float& met, float& metPhi, float& metx, float& mety);
bool passesTrackCuts();


TString ScanChain( TChain* chain, bool runningonGEN, bool requireTrackCuts = true, std::vector<unsigned int> v_goodRuns = std::vector<unsigned int>(), int nEvents = -1) {


  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=0;
  if(nEvents==-1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  //TDirectory *rootdir = gDirectory->GetDirectory("Rint:"); //wasn't used

  // file loop
  TIter fileIter(listOfFiles);

  vector<string> v_prefix;
  vector<string> v_title;
  if(runningonGEN) {
    v_prefix.push_back("mcpt_");
    v_prefix.push_back("mcft_"); 
    v_prefix.push_back("mcall_");
    v_title.push_back(" for MC events passing the trigger Requirements");
    v_title.push_back(" for MC events failing the trigger Requirements");
    v_title.push_back(" for all MC events");
  } else {
    v_prefix.push_back("pt_");
    v_prefix.push_back("ft_"); //pt = passed trigger, ft = failed trigger
    v_prefix.push_back("all_");
    v_title.push_back(" for events passing the trigger Requirements");
    v_title.push_back(" for events failing the trigger Requirements");
    v_title.push_back(" for all events");
  }

  if(v_prefix.size() != v_title.size() ) {
    cout << "The vector of prefixes and the vector of title are not the same size!!! Exiting!" << endl;
    return 0;
  }

  unsigned int aSize = v_prefix.size();

  //fkw's new histos for hit patterns and conversion study
  TH1F *h_trks_innerlayers[aSize];
  TH1F *h_trks_outerlayers[aSize];
  TH1F *h_els_innerlayers[aSize];
  TH1F *h_els_outerlayers[aSize];
  TH1F *h_els_ecalEnergy[aSize];
  TH1F *h_els_type[aSize];

  TH1F *h_els_ecalEnergy_EcalDriven[aSize];
  TH1F *h_els_ecalEnergy_TrkDriven[aSize];
  TH1F *h_els_ecalEnergy_TrkEcalDriven[aSize];

  TH1F *h_els_innerlayers_EcalDriven[aSize];
  TH1F *h_els_innerlayers_TrkDriven[aSize];
  TH1F *h_els_innerlayers_TrkEcalDriven[aSize];

  TH2F *h_nvtxs[aSize];
  TH1F *h_nvtxsgood[aSize];
  TH1F *h_nvtxsbad[aSize];
  TH1F *h_ratioGoodTracks[aSize];
  TH1F *h_numGoodTracks[aSize];
  TH1F *h_numTracksVtx[aSize];
  TH1F *h_zvtx[aSize];
  TH1F *h_rvtx[aSize];

  
  //book histos
  //fkw's new histos come first
  for(unsigned int i = 0; i < v_prefix.size(); i++) {
    h_trks_innerlayers[i] = new TH1F((v_prefix.at(i)+"trks_innerlayers").c_str(), 
			    ("trks_innerlayers" + v_title.at(i)).c_str(), 10, 0, 10);
    h_trks_outerlayers[i] = new TH1F((v_prefix.at(i)+"trks_outerlayers").c_str(), 
			    ("trks_outerlayers" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_innerlayers[i] = new TH1F((v_prefix.at(i)+"els_innerlayers").c_str(), 
			    ("els_innerlayers" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_innerlayers_EcalDriven[i] = new TH1F((v_prefix.at(i)+"els_innerlayers_EcalDriven").c_str(), 
			    ("els_innerlayers_EcalDriven" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_innerlayers_TrkDriven[i] = new TH1F((v_prefix.at(i)+"els_innerlayers_TrkDriven").c_str(), 
			    ("els_innerlayers_TrkDriven" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_innerlayers_TrkEcalDriven[i] = new TH1F((v_prefix.at(i)+"els_innerlayers_TrkEcalDriven").c_str(), 
			    ("els_innerlayers_TrkEcalDriven" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_outerlayers[i] = new TH1F((v_prefix.at(i)+"els_outerlayers").c_str(), 
			    ("els_outerlayers" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_ecalEnergy[i] = new TH1F((v_prefix.at(i)+"els_ecalEnergy").c_str(), 
			    ("els_ecalEnergy" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_ecalEnergy_EcalDriven[i] = new TH1F((v_prefix.at(i)+"els_ecalEnergy_EcalDriven").c_str(), 
			    ("els_ecalEnergy_EcalDriven" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_ecalEnergy_TrkDriven[i] = new TH1F((v_prefix.at(i)+"els_ecalEnergy_TrkDriven").c_str(), 
			    ("els_ecalEnergy_TrkDriven" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_ecalEnergy_TrkEcalDriven[i] = new TH1F((v_prefix.at(i)+"els_ecalEnergy_TrkEcalDriven").c_str(), 
			    ("els_ecalEnergy_TrkEcalDriven" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_type[i] = new TH1F((v_prefix.at(i)+"els_type").c_str(), 
			    ("els_type" + v_title.at(i)).c_str(), 20, 0, 20);
    h_nvtxs[i] = new TH2F((v_prefix.at(i)+"nvtx").c_str(), 
			    ("nvtx" + v_title.at(i)).c_str(), 5, 0, 5, 5,0,5);
    h_nvtxsgood[i] = new TH1F((v_prefix.at(i)+"nvtxgood").c_str(), 
			    ("nvtxgood" + v_title.at(i)).c_str(), 15, 0, 5.0);
    h_nvtxsbad[i] = new TH1F((v_prefix.at(i)+"nvtxbad").c_str(), 
			    ("nvtxbad" + v_title.at(i)).c_str(), 15, 0, 5.0);
    h_ratioGoodTracks[i] = new TH1F((v_prefix.at(i)+"ratioGoodTracks").c_str(), 
			    ("ratioGoodTracks" + v_title.at(i)).c_str(), 101, 0, 1.01);
    h_numGoodTracks[i] = new TH1F((v_prefix.at(i)+"numGoodTracks").c_str(), 
			    ("numGoodTracks" + v_title.at(i)).c_str(), 100, 0, 100);
    h_numTracksVtx[i] = new TH1F((v_prefix.at(i)+"numTracksVtx").c_str(), 
			    ("numTracksVtx" + v_title.at(i)).c_str(), 50, 0, 50);
    h_zvtx[i] = new TH1F((v_prefix.at(i)+"zvtx").c_str(), 
			    ("zvtx" + v_title.at(i)).c_str(), 100, -50, 50);
    h_rvtx[i] = new TH1F((v_prefix.at(i)+"rvtx").c_str(), 
			    ("rvtx" + v_title.at(i)).c_str(), 400, 0, 4.0);

    h_trks_innerlayers[i]-> TH1F::Sumw2();
    h_trks_outerlayers[i]-> TH1F::Sumw2();
    h_els_innerlayers[i]-> TH1F::Sumw2();
    h_els_innerlayers_EcalDriven[i]-> TH1F::Sumw2();
    h_els_innerlayers_TrkDriven[i]-> TH1F::Sumw2();
    h_els_innerlayers_TrkEcalDriven[i]-> TH1F::Sumw2();
    h_els_outerlayers[i]-> TH1F::Sumw2();
    h_els_ecalEnergy[i]-> TH1F::Sumw2();
    h_els_ecalEnergy_EcalDriven[i]-> TH1F::Sumw2();
    h_els_ecalEnergy_TrkDriven[i]-> TH1F::Sumw2();
    h_els_ecalEnergy_TrkEcalDriven[i]-> TH1F::Sumw2();
    h_els_type[i]-> TH1F::Sumw2();
    h_nvtxs[i]->TH2F::Sumw2();
    h_nvtxsgood[i]->TH1F::Sumw2();
    h_nvtxsbad[i]->TH1F::Sumw2();
    h_ratioGoodTracks[i]->TH1F::Sumw2();
    h_numGoodTracks[i]->TH1F::Sumw2();
    h_numTracksVtx[i]->TH1F::Sumw2();
    h_zvtx[i]->TH1F::Sumw2();
    h_rvtx[i]->TH1F::Sumw2();

  }
    
  TFile *currentFile = 0;

  //pass fail counters
  int nGoodEvents = 0;
  int nPassTriggers = 0;
  int nPassTrackingCuts = 0;

  

  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    
    
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;
      
      //fkw: Here's what the next line does. As input argument of ScanChain we have a list of run numbers.
      //fkw: If the present run number is not in this set then find returns the end of the iterator, i.e. one
      //fkw: beyond the last index. In that case we give up right here.
      if(!runningonGEN && find(v_goodRuns.begin(), v_goodRuns.end(), evt_run()) == v_goodRuns.end())
	continue;
      nGoodEvents++;
      
      if(!passesTrackCuts() && requireTrackCuts) 
	continue;
      nPassTrackingCuts++;
      
      int index = (int)(!passesTrigger(runningonGEN));//index is used for histo booking
      if(passesTrigger(runningonGEN))
	nPassTriggers++;
      
      //first comes fkw's new stuff

      for (unsigned int i=0; i< trks_exp_innerlayers().size(); i++){//loop over tracks

	if ( isGoodTrk( i ) ) {
	  h_trks_innerlayers[index]->Fill(min((double)trks_exp_innerlayers().at(i),9.9));
	  h_trks_innerlayers[2]->Fill(min((double)trks_exp_innerlayers().at(i),9.9));
	  h_trks_outerlayers[index]->Fill(min((double)trks_exp_outerlayers().at(i),9.9));
	  h_trks_outerlayers[2]->Fill(min((double)trks_exp_outerlayers().at(i),9.9));
	}

      }

      for (unsigned int i=0; i< els_exp_innerlayers().size(); i++){//loop over electrons

	  h_els_ecalEnergy[index]->Fill(min((double)els_ecalEnergy().at(i),9.9));
	  h_els_ecalEnergy[2]->Fill(min((double)els_ecalEnergy().at(i),9.9));
	  h_els_innerlayers[index]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	  h_els_innerlayers[2]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	  h_els_outerlayers[index]->Fill(min((double)els_exp_outerlayers().at(i),9.9));
	  h_els_outerlayers[2]->Fill(min((double)els_exp_outerlayers().at(i),9.9));
	  h_els_type[index]->Fill(min((double)els_type().at(i),19.9));
	  h_els_type[2]->Fill(min((double)els_type().at(i),19.9));

	  if( els_type().at(i) & (1 << 2) ) {
	    h_els_ecalEnergy_EcalDriven[index]->Fill(min((double)els_ecalEnergy().at(i),9.9));
	    h_els_ecalEnergy_EcalDriven[2]->Fill(min((double)els_ecalEnergy().at(i),9.9));
	    h_els_innerlayers_EcalDriven[index]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	    h_els_innerlayers_EcalDriven[2]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	  }
	  if( els_type().at(i) & (1 << 3) ) {
	    h_els_ecalEnergy_TrkDriven[index]->Fill(min((double)els_ecalEnergy().at(i),9.9));
	    h_els_ecalEnergy_TrkDriven[2]->Fill(min((double)els_ecalEnergy().at(i),9.9));
	    h_els_innerlayers_TrkDriven[index]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	    h_els_innerlayers_TrkDriven[2]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	  }
	  if( (els_type().at(i) & (1 << 3)) && (els_type().at(i) & (1 << 2)) ) {
	    h_els_ecalEnergy_TrkEcalDriven[index]->Fill(min((double)els_ecalEnergy().at(i),9.9));
	    h_els_ecalEnergy_TrkEcalDriven[2]->Fill(min((double)els_ecalEnergy().at(i),9.9));
	    h_els_innerlayers_TrkEcalDriven[index]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	    h_els_innerlayers_TrkEcalDriven[2]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	  }
      }

      //loop over primary vertices
      int nGoodVtxs = 0;
      for(unsigned int i = 0; i < vtxs_isFake().size(); i++) {
	if(vtxs_isFake().at(i))
	  continue;
	h_numTracksVtx[index]->Fill(min((double)vtxs_tracksSize().at(i),49.9));
	h_numTracksVtx[2]->Fill(min((double)vtxs_tracksSize().at(i),49.9));
	h_zvtx[index]->Fill(vtxs_position().at(i).z());
	h_zvtx[2]->Fill(vtxs_position().at(i).z());
	h_rvtx[index]->Fill(vtxs_position().at(i).pt());
	h_rvtx[2]->Fill(vtxs_position().at(i).pt());
	if(vtxs_tracksSize().at(i) < 4 )
	  continue;
	if( fabs( vtxs_position().at(i).z() ) > 15 )
	  continue;
	if( vtxs_position().at(i).pt() > 2 )
	  continue;
	nGoodVtxs++;
      }
      h_nvtxs[index]->Fill(min((double)nGoodVtxs,4.9),min((double)(vtxs_isFake().size()-nGoodVtxs),4.9));  
      h_nvtxs[2]->Fill(min((double)nGoodVtxs,4.9),min((double)(vtxs_isFake().size()-nGoodVtxs),4.9));  
      h_nvtxsgood[index]->Fill(min((double)nGoodVtxs,4.9));  
      h_nvtxsgood[2]->Fill(min((double)nGoodVtxs,4.9));  
      h_nvtxsbad[index]->Fill(min((double)(vtxs_isFake().size()-nGoodVtxs),4.9));  
      h_nvtxsbad[2]->Fill(min((double)(vtxs_isFake().size()-nGoodVtxs),4.9));  

      //make some simple tracking plots
      int nGoodTrks = 0;
      for(unsigned int i = 0; i < trks_trk_p4().size(); i++) {
	if(trks_qualityMask().at(i) & 4)
	  nGoodTrks++;
      }
      double rat = 0;
      if (trks_trk_p4().size() == 0) rat = 0.0;
      else rat = ((double)nGoodTrks)/((double)trks_trk_p4().size());
      h_ratioGoodTracks[index]->Fill(rat);
      h_numGoodTracks[index]->Fill(min((double)nGoodTrks,99.9));
      h_ratioGoodTracks[2]->Fill(rat);
      h_numGoodTracks[2]->Fill(min((double)nGoodTrks,99.9));

      //now comes all of what was there before      
      
    }//event loop
  }//file loop

  cout << "********************SUMMARY********************" << endl;
  cout << "Total number of events: " << nEventsTotal << endl;
  cout << "Total number of events from selected runs : " << nGoodEvents << " (" << 100*(double)nGoodEvents/nEventsTotal << "%)" << endl;
  cout << "Total number of events from selected runs that pass the triggers: " << nPassTriggers << " (" << 100*(double)nPassTriggers/nGoodEvents << "%)" << endl;
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  
  
  TString cutDescription = "";
  if(requireTrackCuts)
    cutDescription = cutDescription + "Require Event level track quality cuts\n";
  if(v_goodRuns.size() != 0)
    cutDescription = cutDescription + "for selected Runs: ";
  for(unsigned int i = 0; i < v_goodRuns.size(); i++) {
    if(i == v_goodRuns.size() -1)
      cutDescription = cutDescription + Form("%d", v_goodRuns.at(i));
    else
      cutDescription = cutDescription + Form("%d", v_goodRuns.at(i)) + ",";
  }
  cutDescription = cutDescription + "\n";
  cout << cutDescription << endl;
  return cutDescription;
}

bool isGoodTrk(unsigned int idx ){
  //will fill this later
  if(trks_qualityMask().at(idx) & 4){
    return true;
  } else {
    return false;
  }

}

bool passesTrigger(bool runningonGEN) {
  
  //if(passL1Trigger("L1_SingleHfBitCountsRing1_1"))
  //return true;
  
  //if(passL1Trigger("L1_SingleHfBitCountsRing2_1"))
  //return true;
  
  //Beam Halo triggers
  if(l1_techbits2() & (1<<4) || l1_techbits2() & (1<<5) || l1_techbits2() & (1<<6) || l1_techbits2() & (1<<7) )
    return false;

  //BPTX triggers
  if(!(l1_techbits1() & (1<<0)) && !runningonGEN)
    return false;

  //BSC triggers
  if(l1_techbits2() & (1<<8) || l1_techbits2() & (1 << 9))
    return true;
  
  return false;
}


bool passesTrackCuts() {
    
  int nGoodVtxs = 0;
  for(unsigned int i = 0; i < vtxs_isFake().size(); i++) {
    if(vtxs_isFake().at(i))
      continue;
    if(vtxs_tracksSize().at(i) < 4 )
      continue;
    if( fabs( vtxs_position().at(i).z() ) > 15 )
      continue;
    if( vtxs_position().at(i).pt() > 2 )
      continue;
    nGoodVtxs++;
  }
  
  //require that there be at least one good vertex
  if(nGoodVtxs==0)
    return false;
  
  //require less than 100 tracks
  if(trks_trk_p4().size() < 10){
    return true;
  } else {

    //require that the fraction of highPurity tracks be > 50%
    int nGoodTrks = 0;
    for(unsigned int i = 0; i < trks_trk_p4().size(); i++) {
      if(trks_qualityMask().at(i) & 4)
	nGoodTrks++;
    }
    if((float)nGoodTrks/trks_trk_p4().size() < 0.2)
      return false;
  }

  return true;
}


