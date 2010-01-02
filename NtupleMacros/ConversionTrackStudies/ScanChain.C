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
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"

#define MAXMET 19.8


CMS2 cms2;

using namespace tas;
bool isGoodTrk( unsigned int i );
//bool isGoodEvent();
bool passesTrigger(bool runningonGEN);
//void getMETQuantities(const float metin, const float metPhiin, 
//		      float& met, float& metPhi, float& metx, float& mety);
bool passesTrackCuts();
std::pair<float, float> getConversionInfo(int idx1, int idx2, float bfiled); 

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
  TH1F *h_els_innerlayers_TrkDrivenLowEt[aSize];
  TH1F *h_els_innerlayers_TrkDrivenHiEt[aSize];
  TH1F *h_els_innerlayers_TrkEcalDriven[aSize];

  TH2F *h_nvtxs[aSize];
  TH1F *h_nvtxsgood[aSize];
  TH1F *h_nvtxsbad[aSize];
  TH1F *h_ratioGoodTracks[aSize];
  TH1F *h_numGoodTracks[aSize];
  TH1F *h_numTracksVtx[aSize];
  TH1F *h_zvtx[aSize];
  TH1F *h_rvtx[aSize];

  TH2F *h_dcotvsdistElectrons[aSize];//looking for conversion partner to an electron
  TH2F *h_dcotvsdistEcalElectrons[aSize];//looking for conversion partner to an electron
  TH2F *h_dcotvsdistTrackElectrons[aSize];//looking for conversion partner to an electron
  TH2F *h_dcotvsdistTracks[aSize];//looking for conversion partner to a track

  TH1F *h_dcotElectrons[aSize];//looking for conversion partner to an electron
  TH1F *h_dcotEcalElectrons[aSize];//looking for conversion partner to an electron
  TH1F *h_dcotTrackElectrons[aSize];//looking for conversion partner to an electron
  TH1F *h_dcotTracks[aSize];//looking for conversion partner to a track
  TH1F *h_dcotTracksMissInn[aSize];//looking for conversion partner to a track
  TH1F *h_dcotTracksNoMissInn[aSize];//looking for conversion partner to a track

  TH2F *h_scs_emaxvsediff[aSize];//eMax vs eMax-e3x3 for superclusters
  TH2F *h_scs_emaxvserat[aSize];//eMax vs (e3x3-eMax)/eMax for superclusters
  TH2F *h_scs_emaxvseratZoom[aSize];//eMax vs (e3x3-eMax)/eMax for superclusters
  TH1F *h_scs_erat[aSize];//(e3x3-eMax)/eMax for superclusters
  TH1F *h_scs_eratCut[aSize];//(e3x3-eMax)/eMax for superclusters
  TH1F *h_scs_emax[aSize];//eMax for superclusters
  TH1F *h_scs_ediff[aSize];//eMax-e3x3 for superclusters
 
  TH1F *h_scs_eta[aSize];//eta of spikes above 10GeV in energy
  TH2F *h_scs_etavsphi[aSize];//eta vs phi of spikes above 10GeV in energy
  TH1F *h_cmetCut[aSize];//met for events with cut on erat and emax
  TH1F *h_cmet[aSize];//met 
  
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
    h_els_innerlayers_TrkDrivenLowEt[i] = new TH1F((v_prefix.at(i)+"els_innerlayers_TrkDrivenLowEt").c_str(), 
			    ("els_innerlayers_TrkDrivenLowEt" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_innerlayers_TrkDrivenHiEt[i] = new TH1F((v_prefix.at(i)+"els_innerlayers_TrkDrivenHiEt").c_str(), 
			    ("els_innerlayers_TrkDrivenHiEt" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_innerlayers_TrkEcalDriven[i] = new TH1F((v_prefix.at(i)+"els_innerlayers_TrkEcalDriven").c_str(), 
			    ("els_innerlayers_TrkEcalDriven" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_outerlayers[i] = new TH1F((v_prefix.at(i)+"els_outerlayers").c_str(), 
			    ("els_outerlayers" + v_title.at(i)).c_str(), 10, 0, 10);
    h_els_ecalEnergy[i] = new TH1F((v_prefix.at(i)+"els_ecalEnergy").c_str(), 
			    ("els_ecalEnergy" + v_title.at(i)).c_str(), 100, 0, 40);
    h_els_ecalEnergy_EcalDriven[i] = new TH1F((v_prefix.at(i)+"els_ecalEnergy_EcalDriven").c_str(), 
			    ("els_ecalEnergy_EcalDriven" + v_title.at(i)).c_str(), 100, 0, 40);
    h_els_ecalEnergy_TrkDriven[i] = new TH1F((v_prefix.at(i)+"els_ecalEnergy_TrkDriven").c_str(), 
			    ("els_ecalEnergy_TrkDriven" + v_title.at(i)).c_str(), 100, 0, 40);
    h_els_ecalEnergy_TrkEcalDriven[i] = new TH1F((v_prefix.at(i)+"els_ecalEnergy_TrkEcalDriven").c_str(), 
			    ("els_ecalEnergy_TrkEcalDriven" + v_title.at(i)).c_str(), 100, 0, 40);
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

    h_dcotvsdistElectrons[i] = new TH2F((v_prefix.at(i)+"dcotvsdistElectron").c_str(), 
			    ("dcotvsdistElectron" + v_title.at(i)).c_str(), 40, -0.2, 0.2, 40, -0.2, 0.2 );
    h_dcotvsdistEcalElectrons[i] = new TH2F((v_prefix.at(i)+"dcotvsdistEcalElectron").c_str(), 
			    ("dcotvsdistEcalElectron" + v_title.at(i)).c_str(), 40, -0.2, 0.2, 40, -0.2, 0.2 );
    h_dcotvsdistTrackElectrons[i] = new TH2F((v_prefix.at(i)+"dcotvsdistTrackElectron").c_str(), 
			    ("dcotvsdistTrackElectron" + v_title.at(i)).c_str(), 40, -0.2, 0.2, 40, -0.2, 0.2 );
    h_dcotvsdistTracks[i] = new TH2F((v_prefix.at(i)+"dcotvsdistTrack").c_str(), 
			    ("dcotvsdistTrack" + v_title.at(i)).c_str(), 40, -0.2, 0.2, 40, -0.2, 0.2 );

    h_dcotElectrons[i] = new TH1F((v_prefix.at(i)+"dcotElectron").c_str(), 
			    ("dcotElectron" + v_title.at(i)).c_str(), 40, -0.2, 0.2 );
    h_dcotEcalElectrons[i] = new TH1F((v_prefix.at(i)+"dcotEcalElectron").c_str(), 
			    ("dcotEcalElectron" + v_title.at(i)).c_str(), 40, -0.2, 0.2 );
    h_dcotTrackElectrons[i] = new TH1F((v_prefix.at(i)+"dcotTrackElectron").c_str(), 
			    ("dcotTrackElectron" + v_title.at(i)).c_str(), 40, -0.2, 0.2 );
    h_dcotTracks[i] = new TH1F((v_prefix.at(i)+"dcotTrack").c_str(), 
			    ("dcotTrack" + v_title.at(i)).c_str(), 40, -0.2, 0.2 );
    h_dcotTracksMissInn[i] = new TH1F((v_prefix.at(i)+"dcotTrackMissInn").c_str(), 
			    ("dcotTrackMissInn" + v_title.at(i)).c_str(), 40, -0.2, 0.2 );
    h_dcotTracksNoMissInn[i] = new TH1F((v_prefix.at(i)+"dcotTrackNoMissInn").c_str(), 
			    ("dcotTrackNoMissInn" + v_title.at(i)).c_str(), 40, -0.2, 0.2 );

    h_scs_emaxvsediff[i] = new TH2F((v_prefix.at(i)+"scs_emaxvsediff").c_str(), 
			    ("scs_emaxvsediff" + v_title.at(i)).c_str(), 40, 0.0, 4.0, 40, 0.0, 10.0 );
    h_scs_emaxvserat[i] = new TH2F((v_prefix.at(i)+"scs_emaxvserat").c_str(), 
			    ("scs_emaxvserat" + v_title.at(i)).c_str(), 41, -0.025, 1.0, 40, 0.0, 40.0 );
    h_scs_emaxvseratZoom[i] = new TH2F((v_prefix.at(i)+"scs_emaxvseratZoom").c_str(), 
			    ("scs_emaxvseratZoom" + v_title.at(i)).c_str(), 40, -0.05, 0.05, 40, 0.0, 40.0 );
    h_scs_emax[i] = new TH1F((v_prefix.at(i)+"scs_emax").c_str(), 
			     ("scs_emax" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_scs_ediff[i] = new TH1F((v_prefix.at(i)+"scs_ediff").c_str(), 
			    ("scs_ediff" + v_title.at(i)).c_str(), 110, -1.0, 10.0);
    h_scs_erat[i] = new TH1F((v_prefix.at(i)+"scs_erat").c_str(), 
			    ("scs_erat" + v_title.at(i)).c_str(), 220, -0.2, 2.0);
    h_scs_eratCut[i] = new TH1F((v_prefix.at(i)+"scs_eratCut").c_str(), 
			    ("scs_eratCut" + v_title.at(i)).c_str(), 220, -0.2, 2.0);

    h_scs_eta[i] = new TH1F((v_prefix.at(i)+"scs_eta").c_str(), 
			    ("scs_eta" + v_title.at(i)).c_str(), 200, -3.0, 3.0);
    h_scs_etavsphi[i] = new TH2F((v_prefix.at(i)+"scs_etavsphi").c_str(), 
			    ("scs_etavsphi" + v_title.at(i)).c_str(), 30, -TMath::Pi(), TMath::Pi(), 60, -3.0, 3.0);

    h_cmetCut[i] = new TH1F((v_prefix.at(i)+"cmetCut").c_str(), 
			    ("cmetCut" + v_title.at(i)).c_str(), 100, 0.0, 50.0);

    h_cmet[i] = new TH1F((v_prefix.at(i)+"cmet").c_str(), 
			    ("cmet" + v_title.at(i)).c_str(), 100, 0.0, 50.0);


    h_trks_innerlayers[i]-> TH1F::Sumw2();
    h_trks_outerlayers[i]-> TH1F::Sumw2();
    h_els_innerlayers[i]-> TH1F::Sumw2();
    h_els_innerlayers_EcalDriven[i]-> TH1F::Sumw2();
    h_els_innerlayers_TrkDriven[i]-> TH1F::Sumw2();
    h_els_innerlayers_TrkDrivenLowEt[i]-> TH1F::Sumw2();
    h_els_innerlayers_TrkDrivenHiEt[i]-> TH1F::Sumw2();
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

    h_dcotvsdistElectrons[i]->TH2F::Sumw2();
    h_dcotvsdistEcalElectrons[i]->TH2F::Sumw2();
    h_dcotvsdistTrackElectrons[i]->TH2F::Sumw2();
    h_dcotvsdistTracks[i]->TH2F::Sumw2();

    h_dcotElectrons[i]->TH1F::Sumw2();
    h_dcotEcalElectrons[i]->TH1F::Sumw2();
    h_dcotTrackElectrons[i]->TH1F::Sumw2();
    h_dcotTracks[i]->TH1F::Sumw2();
    h_dcotTracksMissInn[i]->TH1F::Sumw2();
    h_dcotTracksNoMissInn[i]->TH1F::Sumw2();

    h_scs_emaxvsediff[i]->TH2F::Sumw2();
    h_scs_emaxvserat[i]->TH2F::Sumw2();
    h_scs_emaxvseratZoom[i]->TH2F::Sumw2();
    h_scs_erat[i]->TH1F::Sumw2();
    h_scs_eratCut[i]->TH1F::Sumw2();
    h_scs_emax[i]->TH1F::Sumw2();
    h_scs_ediff[i]->TH1F::Sumw2();
    h_scs_eta[i]->TH1F::Sumw2();
    h_scs_etavsphi[i]->TH2F::Sumw2();
    h_cmet[i]->TH1F::Sumw2();
    h_cmetCut[i]->TH1F::Sumw2();

  }
    
  TFile *currentFile = 0;

  //pass fail counters
  int nGoodEvents = 0;
  int nUnmatchedElectrons = 0;
  int nPassTriggers = 0;
  int nPassTrackingCuts = 0;

  std::vector<float> v_erat;  
  std::vector<float> v_emax;  
  std::vector<float> v_eta;  
  std::vector<float> v_phi;  
  std::vector<float> v_cmet;  
  std::vector<int> v_run;

  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    
    int thisRun = 0;
    
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;
      
      //fkw: Here's what the next line does. As input argument of ScanChain we have a list of run numbers.
      //fkw: If the present run number is not in this set then find returns the end of the iterator, i.e. one
      //fkw: beyond the last index. In that case we give up right here.
      if(!runningonGEN && find(v_goodRuns.begin(), v_goodRuns.end(), evt_run()) == v_goodRuns.end())
	continue;
      nGoodEvents++;
      
      if( thisRun != evt_run() ) {
	thisRun = evt_run();
	cout << "Reached run " << thisRun << endl; 
      }

      if(!passesTrackCuts() && requireTrackCuts) 
	continue;
      nPassTrackingCuts++;
      
      int index = (int)(!passesTrigger(runningonGEN));//index is used for histo booking
      if(passesTrigger(runningonGEN))
	nPassTriggers++;
      
      //first comes fkw's new stuff

      for (unsigned int i=0; i< trks_exp_innerlayers().size(); i++){//loop over tracks

	if( trks_trk_p4().at(i).pt() < 4.0 ) continue;
	if(!(trks_qualityMask().at(i) & 4)) continue;
	h_trks_innerlayers[index]->Fill(min((double)trks_exp_innerlayers().at(i),9.9));
	h_trks_innerlayers[2]->Fill(min((double)trks_exp_innerlayers().at(i),9.9));
	h_trks_outerlayers[index]->Fill(min((double)trks_exp_outerlayers().at(i),9.9));
	h_trks_outerlayers[2]->Fill(min((double)trks_exp_outerlayers().at(i),9.9));

      }

      for (unsigned int i=0; i< els_exp_innerlayers().size(); i++){//loop over electrons

	  h_els_ecalEnergy[index]->Fill(min((double)els_ecalEnergy().at(i),39.99));
	  h_els_ecalEnergy[2]->Fill(min((double)els_ecalEnergy().at(i),39.99));
	  h_els_innerlayers[index]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	  h_els_innerlayers[2]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	  h_els_outerlayers[index]->Fill(min((double)els_exp_outerlayers().at(i),9.9));
	  h_els_outerlayers[2]->Fill(min((double)els_exp_outerlayers().at(i),9.9));
	  h_els_type[index]->Fill(min((double)els_type().at(i),19.9));
	  h_els_type[2]->Fill(min((double)els_type().at(i),19.9));

	  if( els_type().at(i) & (1 << 2) ) {
	    h_els_ecalEnergy_EcalDriven[index]->Fill(min((double)els_ecalEnergy().at(i),39.99));
	    h_els_ecalEnergy_EcalDriven[2]->Fill(min((double)els_ecalEnergy().at(i),39.99));
	    h_els_innerlayers_EcalDriven[index]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	    h_els_innerlayers_EcalDriven[2]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	  }
	  if( els_type().at(i) & (1 << 3) ) {
	    if( els_p4().at(i).pt() > 4.0){
	      h_els_innerlayers_TrkDrivenHiEt[index]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	      h_els_innerlayers_TrkDrivenHiEt[2]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	    } else {
	      h_els_innerlayers_TrkDrivenLowEt[index]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	      h_els_innerlayers_TrkDrivenLowEt[2]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	    }
	  }
	  if( els_type().at(i) & (1 << 3) && !(els_type().at(i) & (1 << 2)) ) {
	    h_els_ecalEnergy_TrkDriven[index]->Fill(min((double)els_ecalEnergy().at(i),39.99));
	    h_els_ecalEnergy_TrkDriven[2]->Fill(min((double)els_ecalEnergy().at(i),39.99));
	    h_els_innerlayers_TrkDriven[index]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	    h_els_innerlayers_TrkDriven[2]->Fill(min((double)els_exp_innerlayers().at(i),9.9));
	  }
	  if( (els_type().at(i) & (1 << 3)) && (els_type().at(i) & (1 << 2)) ) {
	    h_els_ecalEnergy_TrkEcalDriven[index]->Fill(min((double)els_ecalEnergy().at(i),39.99));
	    h_els_ecalEnergy_TrkEcalDriven[2]->Fill(min((double)els_ecalEnergy().at(i),39.99));
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

      //if( thisRun == 123970) cout << "now starting conversions" << endl;

      //look for conversions

      //tracks on tracks first
      float dcot = 1.0;
      float dist = 1.0;
      for(unsigned int i1 = 0; i1 < trks_trk_p4().size(); i1++) {
	if( trks_trk_p4().at(i1).pt() < 1.0 ) continue;
	if(!(trks_qualityMask().at(i1) & 4)) continue;
	dcot = 1.0;
	dist = 1.0;
	for(unsigned int i2 = i1+1; i2 < trks_trk_p4().size(); i2++) {
	  if(!(trks_qualityMask().at(i2) & 4)) continue;
	  if( trks_charge().at(i1)*trks_charge().at(i2) > 0 ) continue;//skip same sign track pairs
	  std::pair<float, float> temp = getConversionInfo(i1, i2, evt_bField());
	  if( fabs(temp.second) < fabs(dcot) ){//find pair with smalles abs(dcot)
	    dcot = temp.second;
	    dist = temp.first;
	  }
	}
	h_dcotvsdistTracks[index]->Fill(dist,dcot);
	h_dcotvsdistTracks[2]->Fill(dist,dcot);
	if( dist > -0.07 && dist < 0.01) {
	  h_dcotTracks[index]->Fill(dcot);
	  h_dcotTracks[2]->Fill(dcot);
	  if ( trks_exp_innerlayers().at(i1) > 0) {
	    h_dcotTracksMissInn[index]->Fill(dcot);
	    h_dcotTracksMissInn[2]->Fill(dcot);
	  } else {
	    h_dcotTracksNoMissInn[index]->Fill(dcot);
	    h_dcotTracksNoMissInn[2]->Fill(dcot);
	  }
	}
      }

      //if( thisRun == 123970) cout << "passed the track block " << nGoodEvents << endl;

      //tracks on electrons next
      for(int i1 = 0; i1 < els_trk_p4().size(); i1++) {
	if( els_trkidx().at(i1) < 0 ) {
	  cout << " els not matched to any track " << endl;
	  nUnmatchedElectrons++;
	  continue;
	}
	if(!(trks_qualityMask().at(els_trkidx().at(i1)) & 4)) continue;
	dcot = 1.0;
	dist = 1.0;
	for(int i2 = 0; i2 < trks_trk_p4().size(); i2++) {
	  if(!(trks_qualityMask().at(i2) & 4)) continue;
	  if(els_trkidx().at(i1) == i2 && els_trkshFrac().at(i1) > 0.45 ) continue;//ignore trk if same as els trk
	  if( els_charge().at(i1)*trks_charge().at(i2) > 0 ) continue;//skip same sign track pairs
	  std::pair<float, float> temp = getConversionInfo(els_trkidx().at(i1), i2, evt_bField());
	  if( fabs(temp.second) < fabs(dcot) ){//find pair with smalles abs(dcot)
	    dcot = temp.second;
	    dist = temp.first;
	  }
	}
	h_dcotvsdistElectrons[index]->Fill(dist,dcot);
	h_dcotvsdistElectrons[2]->Fill(dist,dcot);
	h_dcotElectrons[index]->Fill(dcot);
	h_dcotElectrons[2]->Fill(dcot);
	if( (els_type().at(i1) & (1 << 3)) && !(els_type().at(i1) & (1 << 2)) ) {//track but not ecal
	  h_dcotvsdistTrackElectrons[index]->Fill(dist,dcot);
	  h_dcotvsdistTrackElectrons[2]->Fill(dist,dcot);
	  h_dcotTrackElectrons[index]->Fill(dcot);
	  h_dcotTrackElectrons[2]->Fill(dcot);
	} else if( (els_type().at(i1) & (1 << 2)) ) {//ecal, including those that are track as well
	  h_dcotvsdistEcalElectrons[index]->Fill(dist,dcot);
	  h_dcotvsdistEcalElectrons[2]->Fill(dist,dcot);
	  h_dcotEcalElectrons[index]->Fill(dcot);
	  h_dcotEcalElectrons[2]->Fill(dcot);
	}

      }

      //if( thisRun == 123970) cout << "passed the electron block" << endl;

      h_cmet[index]->Fill(evt_met());
      h_cmet[2]->Fill(evt_met());
      bool passed = false;
      bool failed = false;
      for( unsigned int i=0; i < scs_eMax().size(); i++ ) {
	if( (fabs(scs_eta().at(i) - 1.5) < 0.0500001 && fabs(scs_phi().at(i) - 1.57) < 0.11 ) ) {
	  failed = true;
	  continue;
	}
	h_scs_emax[index]->Fill(min((double)scs_eMax().at(i),49.9));
	h_scs_emax[2]->Fill(min((double)scs_eMax().at(i),49.9));
	h_scs_ediff[index]->Fill(min((double)(scs_e3x3().at(i)-scs_eMax().at(i)),9.99));
	h_scs_ediff[2]->Fill(min((double)(scs_e3x3().at(i)-scs_eMax().at(i)),9.99));
	h_scs_emaxvsediff[index]->Fill(scs_e3x3().at(i)-scs_eMax().at(i),scs_eMax().at(i));
	h_scs_emaxvsediff[2]->Fill(scs_e3x3().at(i)-scs_eMax().at(i),scs_eMax().at(i));
	double erat = (scs_e3x3().at(i)-scs_eMax().at(i))/scs_eMax().at(i) ;
	float emax = min((double)scs_eMax().at(i),39.9);
	if( fabs(erat) < 0.006 && emax > 10.0 ){
	  passed = true;
	  v_erat.push_back(erat);
	  v_emax.push_back(emax);
	  v_eta.push_back(scs_eta().at(i));
	  v_phi.push_back(scs_phi().at(i));
	  v_cmet.push_back(evt_met());
	  v_run.push_back(thisRun);
	  h_scs_eta[index]->Fill(scs_eta().at(i));
	  h_scs_eta[2]->Fill(scs_eta().at(i));
	  h_scs_etavsphi[index]->Fill(scs_phi().at(i),scs_eta().at(i));
	  h_scs_etavsphi[2]->Fill(scs_phi().at(i),scs_eta().at(i));
	} 
	h_scs_erat[index]->Fill(erat);
	h_scs_erat[2]->Fill(erat);
	if( emax > 10.0) {
	  h_scs_eratCut[index]->Fill(erat);
	  h_scs_eratCut[2]->Fill(erat);
	}
	h_scs_emaxvseratZoom[index]->Fill(erat, emax);
	h_scs_emaxvseratZoom[2]->Fill(erat, emax);
	erat = min( max( erat , -0.02499) , 0.99 );
	h_scs_emaxvserat[index]->Fill(erat, emax);
	h_scs_emaxvserat[2]->Fill(erat, emax);
      }
      if( passed && !failed ) {
	h_cmetCut[index]->Fill(evt_met());
	h_cmetCut[2]->Fill(evt_met());
      }

      //now comes all of what was there before      
      
    }//event loop
  }//file loop

  cout << "********************SUMMARY********************" << endl;
  cout << "Total number of events: " << nEventsTotal << endl;
  cout << "Total number of events from selected runs : " << nGoodEvents << " (" << 100*(double)nGoodEvents/nEventsTotal << "%)" << endl;
  cout << "Total number of events that pass tracking cuts: " << nPassTrackingCuts << " (" << 100*(double)nPassTrackingCuts/nGoodEvents << "%)" << endl;
  cout << "Total number of events from selected runs that pass the triggers and tracking cuts: " << nPassTriggers << " (" << 100*(double)nPassTriggers/nGoodEvents << "%)" << endl;
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  cout << "Total number of unmatched electrons = " << nUnmatchedElectrons  << endl; 
  
  //now print out the ecal spike info collected:
  for( unsigned int i=0; i < v_erat.size(); i++){
    cout << " Candidate spike  " << i << endl;
    cout << " erat = " << v_erat.at(i) << endl;
    cout << " emax = " << v_emax.at(i) << endl;
    cout << " eta = " << v_eta.at(i) << endl;
    cout << " phi = " << v_phi.at(i) << endl;
    cout << " calomet = " << v_cmet.at(i) << endl;
    cout << " run number = " << v_run.at(i) << endl;
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

//utility function to get the dist and delta cot theta
std::pair<float, float> getConversionInfo(int idx1, int idx2, float bField){
  
  int trk1_q = trks_charge().at(idx1);
  int trk2_q = trks_charge().at(idx2);
  double trk1_d0 =  trks_d0().at(idx1);
  double trk2_d0 =  trks_d0().at(idx2);
  double trk1_pt = trks_trk_p4().at(idx1).pt();
  double trk2_pt = trks_trk_p4().at(idx2).pt();
  double trk1_phi = trks_trk_p4().at(idx1).phi();
  double trk2_phi = trks_trk_p4().at(idx2).phi();
  double trk1_theta = trks_trk_p4().at(idx1).theta();
  double trk2_theta = trks_trk_p4().at(idx2).theta();

  if( trk1_pt == 0 || trk2_pt == 0 || bField == 0 || trk1_q == 0 || trk2_q == 0 || trk1_theta == 0 || trk2_theta == 0) 
    cout << "about to barf because of division by zero" << endl;
 
  double tk1Curvature = -0.3*bField*(trk1_q/trk1_pt)/100.;
  double rTk1 = fabs(1./tk1Curvature);
  double xTk1 = (1./tk1Curvature - trk1_d0)*cos(trk1_phi);
  double yTk1 = (1./tk1Curvature - trk1_d0)*sin(trk1_phi);
    
  double tk2Curvature = -0.3*bField*(trk2_q/trk2_pt)/100.;
  double rTk2 = fabs(1./tk2Curvature);
  double xTk2 = (1./tk2Curvature - trk2_d0)*cos(trk2_phi);
  double yTk2 = (1./tk2Curvature - trk2_d0)*sin(trk2_phi);
         
  double dist = sqrt(pow(xTk1-xTk2, 2) + pow(yTk1-yTk2 , 2));
  dist = dist - (rTk1 + rTk2);

  double dcot = 1/tan(trk1_theta) - 1/tan(trk2_theta);

  return make_pair(dist, dcot);
  
}
