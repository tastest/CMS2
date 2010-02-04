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
#include <math.h>

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
#include "/home/users/wandrews/CMSSW_3_3_6/src/CMS2/NtupleMacros/Tools/goodrun.cc"

#define MAXMET 19.8


CMS2 cms2;

using namespace tas;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

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

  TH1F *h_scs_eratmax[aSize];				//(e3x3-eMax)/eMax for superclusters
  TH1F *h_scs_eratmaxCut[aSize];			
  TH1F *h_scs_erat3x3[aSize];				
  TH1F *h_scs_erat3x3Cut[aSize]; 
  TH1F *h_scs_eratrat[aSize];				
  TH1F *h_scs_eratratCut[aSize];			
  TH1F *h_scs_eratratCutet[aSize];			
  TH1F *h_scs_eratratCutMet[aSize];			
  TH1F *h_scs_eratratp[aSize];              
  TH1F *h_scs_eratratpMet[aSize];
  TH1F *h_scs_eratratpCut[aSize];

  TH2F *h_scs_emaxvsediff[aSize];		//eMax vs eMax-e3x3 for superclusters
  TH2F *h_scs_e3x3vse2nd[aSize];
  TH2F *h_scs_r9vseratratp[aSize];

  TH2F *h_scs_emaxvseratmax[aSize];		//eMax vs (e3x3-eMax)/eMax for superclusters
  TH2F *h_scs_emaxvseratmaxZoom[aSize];	//eMax vs (e3x3-eMax)/eMax for superclusters
  TH2F *h_scs_emaxvserat3x3[aSize];		//eMax vs (e3x3-eMax)/eMax for superclusters
  TH2F *h_scs_emaxvserat3x3Zoom[aSize];	//eMax vs (e3x3-eMax)/eMax for superclusters
  TH2F *h_scs_emaxvseratrat[aSize];		//eMax vs (e3x3-eMax)/eMax for superclusters
  TH2F *h_scs_emaxvseratratZoom[aSize];	//eMax vs (e3x3-eMax)/eMax for superclusters

  TH2F *h_scs_etmaxvseratmax[aSize];		//eMax vs (e3x3-eMax)/eMax for superclusters
  TH2F *h_scs_etmaxvseratmaxZoom[aSize];	//eMax vs (e3x3-eMax)/eMax for superclusters
  TH2F *h_scs_etmaxvserat3x3[aSize];		//eMax vs (e3x3-eMax)/eMax for superclusters
  TH2F *h_scs_etmaxvserat3x3Zoom[aSize];	//eMax vs (e3x3-eMax)/eMax for superclusters
  TH2F *h_scs_etmaxvseratrat[aSize];		//eMax vs (e3x3-eMax)/eMax for superclusters
  TH2F *h_scs_etmaxvseratratZoom[aSize];	//eMax vs (e3x3-eMax)/eMax for superclusters
  TH2F *h_scs_etmaxvseratratpZoom[aSize];	//eMax vs (e3x3-eMax)/eMax for superclusters

  TH1F *h_scs_emax[aSize];				//eMax for superclusters
  TH1F *h_scs_etmax[aSize];				//eMax for superclusters
  TH1F *h_scs_etmaxnospike[aSize];		//eMax for superclusters
  TH1F *h_scs_emaxCut[aSize];			//eMax for superclusters with erat < eratcut (0.006 or so)
  TH1F *h_scs_etmaxCut[aSize];			//eMax for superclusters with erat < eratcut (0.006 or so)
  TH1F *h_scs_e2nd[aSize];				//eMax for superclusters
  TH1F *h_scs_et2nd[aSize];				//eMax for superclusters
  TH1F *h_scs_ediff[aSize];				//eMax-e3x3 for superclusters
 
  TH1F *h_scs_eta[aSize];            //eta of spikes above 10GeV in energy
  TH2F *h_scs_etavsphi[aSize];       //eta vs phi of spikes above 10GeV in energy
  TH2F *h_scs_etavsphiHot[aSize];    //eta vs phi of spikes above 10GeV in energy--include hot cell
  TH2F *h_scs_etavsphiHotp[aSize];    //eta vs phi of spikes above 10GeV in energy--include hot cell
  TH2F *h_scs_etavsphiNar[aSize];    //eta vs phi of spikes above 10GeV in energy
  TH2F *h_scs_etavsphiHotNar[aSize]; //eta vs phi of spikes above 10GeV in energy--include hot cell
  
  TH1F *h_cmetCut[aSize];//met for events with cut on erat and emax
  TH1F *h_cmet[aSize];//met 
  TH1F *h_cmetCutCorr[aSize];//met for events with cut on erat and emax
  //TH1F *h_cmetCorr[aSize];//met 
  TH1F *h_tcmet[aSize];//met 
  TH1F *h_tcmetCut[aSize];//met for events with cut on erat and emax
  TH1F *h_tcmetHF[aSize];
  TH1F *h_tcmetHFCorr[aSize];
  TH1F *h_tcmetCorr[aSize];
  TH1F *h_tcmetAllCorr[aSize];
  TH1F *h_pfmetCut[aSize];//met for events with cut on erat and emax
  TH1F *h_pfmet[aSize];//met 

  TH2F *h_scs_etvscmet[aSize]; //et of sc vs met
  TH2F *h_scs_etmaxvscmet[aSize]; //eMax_t of sc vs met
  TH2F *h_scs_petmaxvscmet[aSize]; //eMax_t of sc vs met
  TH2F *h_scs_phivscmetphi[aSize]; //sc_phi vs met_phi
  TH2F *h_scs_pphivscmetphi[aSize]; //sc_phi vs met_phi
  TH1F *h_scs_dphicmet[aSize]; //sc_phi - phi_met

  TH2F *h_scs_etvstcmet[aSize]; //et of sc vs met
  TH2F *h_scs_etmaxvstcmet[aSize]; //eMax_t of sc vs met
  TH2F *h_scs_phivstcmetphi[aSize]; //sc_phi vs met_phi
  TH1F *h_scs_dphitcmet[aSize]; //sc_phi - phi_met

  TH2F *h_scs_etvspfmet[aSize]; //et of sc vs met
  TH2F *h_scs_etmaxvspfmet[aSize]; //eMax_t of sc vs met
  TH2F *h_scs_phivspfmetphi[aSize]; //sc_phi vs met_phi
  TH1F *h_scs_dphipfmet[aSize]; //sc_phi - phi_met

  TH1F *h_scs_hoe[aSize];
  TH1F *h_scs_hoe_spike[aSize];
  TH1F *h_scs_hoe_spikep[aSize];

  TH1F *h_twrs_ass[aSize];

  TH1F *h_scs_spikeFlags[aSize];//reco flags for ECAL spikes
  TH1F *h_scs_notspikeFlags[aSize];//reco flags for ECAL !spikes
  TH1F *h_scs_allFlags[aSize];//reco flags for ECAL !spikes

  
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

	h_scs_eratmax[i] = new TH1F((v_prefix.at(i)+"scs_eratmax").c_str(), 
							 ("scs_eratmax" + v_title.at(i)).c_str(), 220, -0.2, 2.0);
    h_scs_eratmaxCut[i] = new TH1F((v_prefix.at(i)+"scs_eratmaxCut").c_str(), 
								("scs_eratmaxCut" + v_title.at(i)).c_str(), 220, -0.2, 2.0);
	h_scs_erat3x3[i] = new TH1F((v_prefix.at(i)+"scs_erat3x3").c_str(), 
							 ("scs_erat3x3" + v_title.at(i)).c_str(), 220, -0.2, 2.0);
    h_scs_erat3x3Cut[i] = new TH1F((v_prefix.at(i)+"scs_erat3x3Cut").c_str(), 
								("scs_erat3x3Cut" + v_title.at(i)).c_str(), 220, -0.2, 2.0);
	h_scs_eratrat[i] = new TH1F((v_prefix.at(i)+"scs_eratrat").c_str(), 
							 ("scs_eratrat" + v_title.at(i)).c_str(), 220, -0.2, 2.0);
    h_scs_eratratCut[i] = new TH1F((v_prefix.at(i)+"scs_eratratCut").c_str(), 
								("scs_eratratCut" + v_title.at(i)).c_str(), 220, -0.2, 2.0);
    h_scs_eratratCutMet[i] = new TH1F((v_prefix.at(i)+"scs_eratratMet").c_str(), 
								("scs_eratratMet" + v_title.at(i)).c_str(), 220, -0.2, 2.0);
    h_scs_eratratCutet[i] = new TH1F((v_prefix.at(i)+"scs_eratratCutet").c_str(), 
								("scs_eratratCutet" + v_title.at(i)).c_str(), 220, -0.2, 2.0);
	h_scs_eratratp[i] = new TH1F((v_prefix.at(i)+"scs_eratratp").c_str(), 
								 ("scs_eratratp" + v_title.at(i)).c_str(), 220, -0.2, 2.0); //prime
    h_scs_eratratpMet[i] = new TH1F((v_prefix.at(i)+"scs_eratratpMET").c_str(), 
								("scs_eratratpMET" + v_title.at(i)).c_str(), 220, -0.2, 2.0);
    h_scs_eratratpCut[i] = new TH1F((v_prefix.at(i)+"scs_eratratpCut").c_str(), 
									("scs_eratratpCut" + v_title.at(i)).c_str(), 220, -0.2, 2.0);


    h_scs_emaxvsediff[i] = new TH2F((v_prefix.at(i)+"scs_emaxvsediff").c_str(), 
			    ("scs_emaxvsediff" + v_title.at(i)).c_str(), 40, 0.0, 4.0, 40, 0.0, 10.0 );
	h_scs_e3x3vse2nd[aSize] = new TH2F((v_prefix.at(i)+"scs_e3x3vse2nd").c_str(), 
			    ("scs_e3x3vse2nd" + v_title.at(i)).c_str(), 40, 0.0, 4.0, 40, 0.0, 10.0 );
	h_scs_r9vseratratp[aSize] = new TH2F((v_prefix.at(i)+"scs_r9vseratratp").c_str(), 
			    ("scs_r9vseratratp" + v_title.at(i)).c_str(), 40, 0.0, 4.0, 40, 0.0, 10.0 );

    h_scs_emaxvseratmax[i] = new TH2F((v_prefix.at(i)+"scs_emaxvseratmax").c_str(), 
			    ("scs_emaxvseratmax" + v_title.at(i)).c_str(), 41, -0.025, 1.0, 40, 0.0, 40.0 );
    h_scs_emaxvseratmaxZoom[i] = new TH2F((v_prefix.at(i)+"scs_emaxvseratmaxZoom").c_str(), 
			    ("scs_emaxvseratmaxZoom" + v_title.at(i)).c_str(), 40, -0.05, 0.05, 40, 0.0, 40.0 );
    h_scs_emaxvserat3x3[i] = new TH2F((v_prefix.at(i)+"scs_emaxvserat3x3").c_str(), 
			    ("scs_emaxvserat3x3" + v_title.at(i)).c_str(), 41, -0.025, 1.0, 40, 0.0, 40.0 );
    h_scs_emaxvserat3x3Zoom[i] = new TH2F((v_prefix.at(i)+"scs_emaxvserat3x3Zoom").c_str(), 
			    ("scs_emaxvserat3x3Zoom" + v_title.at(i)).c_str(), 40, -0.05, 0.05, 40, 0.0, 40.0 );
    h_scs_emaxvseratrat[i] = new TH2F((v_prefix.at(i)+"scs_emaxvseratrat").c_str(), 
			    ("scs_emaxvseratrat" + v_title.at(i)).c_str(), 41, -0.025, 1.0, 40, 0.0, 40.0 );
    h_scs_emaxvseratratZoom[i] = new TH2F((v_prefix.at(i)+"scs_emaxvseratratZoom").c_str(), 
			    ("scs_emaxvseratratZoom" + v_title.at(i)).c_str(), 40, 0.95, 1.05, 40, 0.0, 40.0 );

    h_scs_etmaxvseratmax[i] = new TH2F((v_prefix.at(i)+"scs_etmaxvseratmax").c_str(), 
			    ("scs_etmaxvseratmax" + v_title.at(i)).c_str(), 41, -0.025, 1.0, 40, 0.0, 40.0 );
    h_scs_etmaxvseratmaxZoom[i] = new TH2F((v_prefix.at(i)+"scs_etmaxvseratmaxZoom").c_str(), 
			    ("scs_etmaxvseratmaxZoom" + v_title.at(i)).c_str(), 40, -0.05, 0.05, 40, 0.0, 40.0 );
    h_scs_etmaxvserat3x3[i] = new TH2F((v_prefix.at(i)+"scs_etmaxvserat3x3").c_str(), 
			    ("scs_etmaxvserat3x3" + v_title.at(i)).c_str(), 41, -0.025, 1.0, 40, 0.0, 40.0 );
    h_scs_etmaxvserat3x3Zoom[i] = new TH2F((v_prefix.at(i)+"scs_etmaxvserat3x3Zoom").c_str(), 
			    ("scs_etmaxvserat3x3Zoom" + v_title.at(i)).c_str(), 40, -0.05, 0.05, 40, 0.0, 40.0 );
    h_scs_etmaxvseratrat[i] = new TH2F((v_prefix.at(i)+"scs_etmaxvseratrat").c_str(), 
			    ("scs_etmaxvseratrat" + v_title.at(i)).c_str(), 41, -0.025, 1.0, 40, 0.0, 40.0 );
    h_scs_etmaxvseratratZoom[i] = new TH2F((v_prefix.at(i)+"scs_etmaxvseratratZoom").c_str(), 
			    ("scs_etmaxvseratratZoom" + v_title.at(i)).c_str(), 40, 0.95, 1.05, 40, 0.0, 40.0 );
    h_scs_etmaxvseratratpZoom[i] = new TH2F((v_prefix.at(i)+"scs_etmaxvseratratpZoom").c_str(), 
			    ("scs_etmaxvseratratpZoom" + v_title.at(i)).c_str(), 40, 0.95, 1.05, 40, 0.0, 40.0 );

    h_scs_emax[i] = new TH1F((v_prefix.at(i)+"scs_emax").c_str(), 
			     ("scs_emax" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_scs_etmax[i] = new TH1F((v_prefix.at(i)+"scs_etmax").c_str(), 
			     ("scs_etmax" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_scs_emaxCut[i] = new TH1F((v_prefix.at(i)+"scs_emaxCut").c_str(), 
			     ("scs_emaxCut" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_scs_etmaxCut[i] = new TH1F((v_prefix.at(i)+"scs_etmaxCut").c_str(), 
			     ("scs_etmaxCut" + v_title.at(i)).c_str(), 50, 0.0, 50.0);
    h_scs_e2nd[i] = new TH1F((v_prefix.at(i)+"scs_e2nd").c_str(), 
			     ("scs_e2nd" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_scs_et2nd[i] = new TH1F((v_prefix.at(i)+"scs_et2nd").c_str(), 
			     ("scs_et2nd" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_scs_etmaxnospike[i] = new TH1F((v_prefix.at(i)+"scs_etmaxnospike").c_str(), 
			     ("scs_etmaxnospike" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_scs_ediff[i] = new TH1F((v_prefix.at(i)+"scs_ediff").c_str(), 
			    ("scs_ediff" + v_title.at(i)).c_str(), 110, -1.0, 10.0);

    h_scs_eta[i] = new TH1F((v_prefix.at(i)+"scs_eta").c_str(), 
			    ("scs_eta" + v_title.at(i)).c_str(), 200, -3.0, 3.0);
    h_scs_etavsphi[i] = new TH2F((v_prefix.at(i)+"scs_etavsphi").c_str(), 
			    ("scs_etavsphi" + v_title.at(i)).c_str(), 30, -TMath::Pi(), TMath::Pi(), 60, -3.0, 3.0);
    h_scs_etavsphiHot[i] = new TH2F((v_prefix.at(i)+"scs_etavsphiHot").c_str(), 
			    ("scs_etavsphiHot" + v_title.at(i)).c_str(), 90, -TMath::Pi(), TMath::Pi(), 90, -3.0, 3.0);
    h_scs_etavsphiHotp[i] = new TH2F((v_prefix.at(i)+"scs_etavsphiHotp").c_str(), 
			    ("scs_etavsphiHotp" + v_title.at(i)).c_str(), 90, -TMath::Pi(), TMath::Pi(), 90, -3.0, 3.0);
    h_scs_etavsphiNar[i] = new TH2F((v_prefix.at(i)+"scs_etavsphiNar").c_str(), 
									("scs_etavsphiNar" + v_title.at(i)).c_str(), 420, -TMath::Pi(), TMath::Pi(), 400, -3.0, 3.0); //.015
    h_scs_etavsphiHotNar[i] = new TH2F((v_prefix.at(i)+"scs_etavsphiHotNar").c_str(), 
			    ("scs_etavsphiHotNar" + v_title.at(i)).c_str(), 420, -TMath::Pi(), TMath::Pi(), 400, -3.0, 3.0);

    h_cmetCut[i]  = new TH1F((v_prefix.at(i)+"cmetCut").c_str(), ("cmetCut" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_cmet[i]     = new TH1F((v_prefix.at(i)+"cmet").c_str(), ("cmet" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_tcmetCut[i] = new TH1F((v_prefix.at(i)+"tcmetCut").c_str(), ("tcmetCut" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_tcmetHF[i] = new TH1F((v_prefix.at(i)+"tcmetHF").c_str(), ("tcmetHF" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_tcmet[i]    = new TH1F((v_prefix.at(i)+"tcmet").c_str(), ("tcmet" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_pfmetCut[i] = new TH1F((v_prefix.at(i)+"pfmetCut").c_str(), ("pfmetCut" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_pfmet[i]    = new TH1F((v_prefix.at(i)+"pfmet").c_str(), ("pfmet" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_cmetCutCorr[i]  = new TH1F((v_prefix.at(i)+"cmetCutCorr").c_str(), ("cmetCutCorr" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    //h_cmetCorr[i]     = new TH1F((v_prefix.at(i)+"cmetCorr").c_str(), ("cmetCorr" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_tcmetCorr[i]    = new TH1F((v_prefix.at(i)+"tcmetCorr").c_str(), ("tcmetCorr" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_tcmetHFCorr[i]    = new TH1F((v_prefix.at(i)+"tcmetHFCorr").c_str(), ("tcmetHFCorr" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_tcmetAllCorr[i]    = new TH1F((v_prefix.at(i)+"tcmetAllCorr").c_str(), ("tcmetAllCorr" + v_title.at(i)).c_str(), 100, 0.0, 50.0);

	h_scs_etvscmet[i] = new TH2F((v_prefix.at(i)+"scet_vs_cmet").c_str(), ("scet_vs_cmet" + v_title.at(i)).c_str(), 100, 0.0, 50.0, 100, 0.0, 50.0);
	h_scs_etmaxvscmet[i] = new TH2F((v_prefix.at(i)+"scetmax_vs_cmet").c_str(), ("scetmax_vs_cmet" + v_title.at(i)).c_str(), 100, 0.0, 50.0, 100, 0.0, 50.0);
	h_scs_phivscmetphi[i] = new TH2F((v_prefix.at(i)+"scphi_vs_cmetphi").c_str(), ("scphi_vs_cmetphi" + v_title.at(i)).c_str(), 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi() );
	h_scs_dphicmet[i] = new TH1F((v_prefix.at(i)+"dphi_sc_cmet").c_str(), ("dphi_sc_cmet" + v_title.at(i)).c_str(), 100, -TMath::Pi(), TMath::Pi() );

	h_scs_petmaxvscmet[i] = new TH2F((v_prefix.at(i)+"pscetmax_vs_cmet").c_str(), ("pscetmax_vs_cmet" + v_title.at(i)).c_str(), 100, 0.0, 50.0, 100, 0.0, 50.0);
	h_scs_pphivscmetphi[i] = new TH2F((v_prefix.at(i)+"pscphi_vs_cmetphi").c_str(), ("pscphi_vs_cmetphi" + v_title.at(i)).c_str(), 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi() );

	h_scs_etvstcmet[i] = new TH2F((v_prefix.at(i)+"scet_vs_tcmet").c_str(), ("scet_vs_tcmet" + v_title.at(i)).c_str(), 100, 0.0, 50.0, 100, 0.0, 50.0);
	h_scs_etmaxvstcmet[i] = new TH2F((v_prefix.at(i)+"scetmax_vs_tcmet").c_str(), ("scetmax_vs_tcmet" + v_title.at(i)).c_str(), 100, 0.0, 50.0, 100, 0.0, 50.0);
	h_scs_phivstcmetphi[i] = new TH2F((v_prefix.at(i)+"scphi_vs_tcmetphi").c_str(), ("scphi_vs_tcmetphi" + v_title.at(i)).c_str(), 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi() );
	h_scs_dphitcmet[i] = new TH1F((v_prefix.at(i)+"dphi_sc_tcmet").c_str(), ("dphi_sc_tcmet" + v_title.at(i)).c_str(), 100, -TMath::Pi(), TMath::Pi() );

	h_scs_etvspfmet[i] = new TH2F((v_prefix.at(i)+"scet_vs_pfmet").c_str(), ("scet_vs_pfmet" + v_title.at(i)).c_str(), 100, 0.0, 50.0, 100, 0.0, 50.0);
	h_scs_etmaxvspfmet[i] = new TH2F((v_prefix.at(i)+"scetmax_vs_pfmet").c_str(), ("scetmax_vs_pfmet" + v_title.at(i)).c_str(), 100, 0.0, 50.0, 100, 0.0, 50.0);
	h_scs_phivspfmetphi[i] = new TH2F((v_prefix.at(i)+"scphi_vs_pfmetphi").c_str(), ("scphi_vs_pfmetphi" + v_title.at(i)).c_str(), 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi() );
	h_scs_dphipfmet[i] = new TH1F((v_prefix.at(i)+"dphi_sc_pfmet").c_str(), ("dphi_sc_pfmet" + v_title.at(i)).c_str(), 100, -TMath::Pi(), TMath::Pi() );

	h_scs_hoe[i]  		= new TH1F((v_prefix.at(i)+"scs_hoe").c_str(),        ("scs_hoe" + v_title.at(i)).c_str()       , 100, 0.0, 0.1);
	h_scs_hoe_spike[i]  = new TH1F((v_prefix.at(i)+"scs_hoe_spike").c_str(),  ("scs_hoe_spike" + v_title.at(i)).c_str() , 100, 0.0, 0.1);
	h_scs_hoe_spikep[i] = new TH1F((v_prefix.at(i)+"scs_hoe_spikep").c_str(), ("scs_hoe_spikep" + v_title.at(i)).c_str(), 100, 0.0, 0.1);

	h_twrs_ass[i] = new TH1F((v_prefix.at(i)+"twrs_ass").c_str(), ("twrs_ass" + v_title.at(i)).c_str(), 201, -1.0, 1.01); 

    h_scs_allFlags[i] = new TH1F((v_prefix.at(i)+"scs_allFlags").c_str(), 
								("scs_allFlags" + v_title.at(i)).c_str(), 20, 0.0, 20.0);
    h_scs_spikeFlags[i] = new TH1F((v_prefix.at(i)+"scs_spikeFlags").c_str(), 
								("scs_spikeFlags" + v_title.at(i)).c_str(), 20, 0.0, 20.0);
    h_scs_notspikeFlags[i] = new TH1F((v_prefix.at(i)+"scs_notspikeFlags").c_str(), 
								("scs_notspikeFlags" + v_title.at(i)).c_str(), 20, 0.0, 20.0);


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

    h_scs_emaxvseratmax[i]->TH2F::Sumw2();
    h_scs_emaxvseratmaxZoom[i]->TH2F::Sumw2();
    h_scs_emaxvserat3x3[i]->TH2F::Sumw2();
    h_scs_emaxvserat3x3Zoom[i]->TH2F::Sumw2();
    h_scs_emaxvseratrat[i]->TH2F::Sumw2();
    h_scs_emaxvseratratZoom[i]->TH2F::Sumw2();

    h_scs_eratmax[i]->TH1F::Sumw2();
    h_scs_eratmaxCut[i]->TH1F::Sumw2();
    h_scs_erat3x3[i]->TH1F::Sumw2();
    h_scs_erat3x3Cut[i]->TH1F::Sumw2();
    h_scs_eratrat[i]->TH1F::Sumw2();
    h_scs_eratratCut[i]->TH1F::Sumw2();
    h_scs_eratratCutet[i]->TH1F::Sumw2();
	h_scs_eratratCutMet[i]->TH1F::Sumw2();	
	h_scs_eratratp[i]->TH1F::Sumw2();
	h_scs_eratratpMet[i]->TH1F::Sumw2();
	h_scs_eratratpCut[i]->TH1F::Sumw2();

    h_scs_emax[i]->TH1F::Sumw2();
    h_scs_etmax[i]->TH1F::Sumw2();
    h_scs_ediff[i]->TH1F::Sumw2();
	h_scs_e2nd[i]->TH1F::Sumw2();
	h_scs_et2nd[i]->TH1F::Sumw2();
    h_scs_eta[i]->TH1F::Sumw2();
    h_scs_etavsphi[i]->TH2F::Sumw2();
    h_cmet[i]->TH1F::Sumw2();
    h_cmetCut[i]->TH1F::Sumw2();
    h_cmetCutCorr[i]->TH1F::Sumw2();
    h_tcmet[i]->TH1F::Sumw2();
    h_tcmetCut[i]->TH1F::Sumw2();
    h_tcmetHF[i]->TH1F::Sumw2();
    h_tcmetCorr[i]->TH1F::Sumw2();
    h_tcmetHFCorr[i]->TH1F::Sumw2();
    h_tcmetAllCorr[i]->TH1F::Sumw2();
    h_pfmet[i]->TH1F::Sumw2();
    h_pfmetCut[i]->TH1F::Sumw2();

	h_scs_dphicmet[i]->TH1F::Sumw2();
	h_scs_dphitcmet[i]->TH1F::Sumw2();
	h_scs_dphipfmet[i]->TH1F::Sumw2();

	h_scs_hoe[i]  	   ->TH1F::Sumw2();
	h_scs_hoe_spike[i] ->TH1F::Sumw2();
	h_scs_hoe_spikep[i]->TH1F::Sumw2();

	h_twrs_ass[i]->TH1F::Sumw2();

    h_scs_allFlags[i]->TH1F::Sumw2();
    h_scs_spikeFlags[i]->TH1F::Sumw2();
    h_scs_notspikeFlags[i]->TH1F::Sumw2();
  }
    
  TFile *currentFile = 0;

  //pass fail counters
  int nGoodEvents = 0;
  int nUnmatchedElectrons = 0;
  int nPassTriggers = 0;
  int nPassTrackingCuts = 0;
  vector<int> nGoodEventsPerRun;
  for( unsigned int i=0; i<v_goodRuns.size(); i++ ) {
	nGoodEventsPerRun.push_back(0);
  }

  std::vector<float> v_erat;  
  std::vector<float> v_emax;  
  std::vector<float> v_eta;  
  std::vector<float> v_phi;  
  std::vector<float> v_cmet;  
  std::vector<int> v_run;
  std::vector<int> v_trkmch;
  int nscs_goodrun = 0;           //total scs
  int nscs_badrun = 0;
  int nscs5gevet_goodrun = 0;     //n scs > 5gevet
  int nscs5gevet_badrun = 0;
  int ngoodscs5gevet_goodrun = 0; //n scs > 5gevet and |r9-1|>0.01
  int ngoodscs5gevet_badrun = 0;
  int nspikes_goodrun = 0;        //n scs > 5gevet and |r9-1|<0.01
  int nspikes_badrun = 0;

  int tothfspikes = 0;
    
  int thisRun = 0;
  int npassgoodrun = 0;
  bool passgoodrun = false;

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
      
      if( thisRun != (int)evt_run() ) {
		thisRun = evt_run();
		//cout << "Reached run " << thisRun << endl; 
      }

      if(!passesTrackCuts() && requireTrackCuts) 
		continue;
      nPassTrackingCuts++;
      
      int index = (int)(!passesTrigger(runningonGEN));//index is used for histo booking
      if( passesTrigger(runningonGEN) ) {
		nPassTriggers++;

		//count total good events per run
		for( unsigned int j=0; j<v_goodRuns.size(); j++ ) {
		  if( v_goodRuns.at(j) == evt_run() ) { //match run from runlist with spike's run
			nGoodEventsPerRun.at(j)++;
		  }
		}
	  }
	  else
		continue;

      //first comes fkw's new stuff
	  passgoodrun = false;
	  //insert jmu's goodrun fn + yj's list
	  if( goodrun( evt_run(), evt_lumiBlock() ) ){
	    npassgoodrun++;
		passgoodrun = true;
		//continue;
	  }

	  /*
	  //////////////////////////////////////////////////////////////////////////start block


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
      for(int i1 = 0; i1 < (int)els_trk_p4().size(); i1++) {
		if( els_trkidx().at(i1) < 0 ) {
		  //cout << " els not matched to any track " << endl;
		  nUnmatchedElectrons++;
		  continue;
		}
		if(!(trks_qualityMask().at(els_trkidx().at(i1)) & 4)) continue;
		dcot = 1.0;
		dist = 1.0;
		for(int i2 = 0; i2 < (int)trks_trk_p4().size(); i2++) {
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

*/
	  //////////////////////////////////////////////////////////////////////////end block

      //if( thisRun == 123970) cout << "passed the electron block" << endl;

	  //scs superclusters info
	  const float metmax = 49.9;
	  const float scdiff = 9.99;
	  const float emaxmax = 39.9;
	  const float etmaxcut = 5.0;
	  const float cmet  = min(evt_met(), metmax);
	  const float tcmet = min(evt_tcmet(), metmax);
	  const float pfmet = min(evt_pfmet(), metmax);
	  float tcmetallcorr = evt_tcmet();
      h_cmet[index]		->Fill(cmet);
      h_cmet[2]			->Fill(cmet);
      h_tcmet[index]	->Fill(tcmet);
      h_tcmet[2]		->Fill(tcmet);
      h_pfmet[index]	->Fill(pfmet);
      h_pfmet[2]		->Fill(pfmet);
      bool passed = false;
      bool failed = false;
	  int evt_nspikes = 0;
      for( unsigned int i=0; i < scs_eMax().size(); i++ ) {
		const float eratmax = (scs_e3x3().at(i)-scs_eMax().at(i))/scs_eMax().at(i);
		const float erat3x3 = (scs_e3x3().at(i)-scs_eMax().at(i))/scs_e3x3().at(i);
		const float eratrat = scs_eMax().at(i)/scs_e3x3().at(i);
		const float eratratp= scs_eMax().at(i)/(scs_e3x3().at(i) - scs_e2nd().at(i));
		const float emax    = min( scs_eMax().at(i), emaxmax );
		const float scet    = min( scs_energy().at(i)*sin( scs_pos_p4().at(i).Theta() ), metmax );
		const float scetmax = min( scs_eMax().at(i)  *sin( scs_pos_p4().at(i).Theta() ), metmax );
		const float scetmaxp= min( (scs_eMax().at(i)+scs_e2nd().at(i))*sin( scs_pos_p4().at(i).Theta() ), metmax );
		h_scs_etavsphiHotNar[index]->Fill(scs_phi().at(i),scs_eta().at(i));
		h_scs_etavsphiHotNar[2]    ->Fill(scs_phi().at(i),scs_eta().at(i));

		if( fabs(eratrat-1) < 0.01 && scetmax > etmaxcut ){ //spikes including hot cell
		  h_scs_etavsphiHot[index]->Fill(scs_phi().at(i),scs_eta().at(i));
		  h_scs_etavsphiHot[2]    ->Fill(scs_phi().at(i),scs_eta().at(i));
		}

		//deal with hot cell
		//if( fabs(scs_eta().at(i) - 1.5) < 0.0500001 && fabs(scs_phi().at(i) - 1.57) < 0.11 //eta phi of hot cell--fkw
		if( fabs(scs_eta().at(i) - 1.53) < 0.016 && fabs(scs_phi().at(i) - 1.66) < 0.016 ) { //eta phi of hot cell--warren
		  //&& fabs(eratrat-1) < 0.01 && scetmax > etmaxcut ) { //eratX, emax of hot cell
		  failed = true; //flag event with hot cell
		  continue; //skip hot cell
		}

		// counting of scs
		if( passgoodrun )
		  nscs_goodrun++; //total number
		else
		  nscs_badrun++;

		if( scetmax > etmaxcut ) { //nscs with et cut only
		  if( passgoodrun )
			nscs5gevet_goodrun++;
		  else
			nscs5gevet_badrun++;
		}

		if( fabs(eratrat-1) > 0.01 && scetmax > etmaxcut ) { //r9 cut inverted -- 'good scs' (aka 'no spike')
		  cout << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
		  if( passgoodrun )
			ngoodscs5gevet_goodrun++;
		  else
			ngoodscs5gevet_badrun++;
		}

		if( fabs(eratrat-1) < 0.01 && scetmax > etmaxcut ) { //spikes (5gevet)
		  //cout << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
		  if( passgoodrun )
			nspikes_goodrun++;
		  else
			nspikes_badrun++;
		}

		h_scs_hoe[index]->Fill( scs_hoe().at(i) );
		h_scs_hoe[2]    ->Fill( scs_hoe().at(i) );
		h_scs_allFlags[index]->Fill(scs_severitySeed().at(i));
		h_scs_allFlags[2]->Fill(scs_severitySeed().at(i));
		h_scs_emax[index]->Fill( min(scs_eMax().at(i), metmax) );
		h_scs_emax[2]    ->Fill( min(scs_eMax().at(i), metmax) );
		h_scs_etmax[index]->Fill( scetmax, metmax );
		h_scs_etmax[2]    ->Fill( scetmax, metmax );
		if( fabs(eratrat-1) < 0.01 ) {
		  h_scs_emaxCut[index]->Fill( min(scs_eMax().at(i), metmax) );
		  h_scs_emaxCut[2]    ->Fill( min(scs_eMax().at(i), metmax) );
		  h_scs_etmaxCut[index]->Fill( scetmax );
		  h_scs_etmaxCut[2]    ->Fill( scetmax );
		}
		h_scs_ediff[index]->Fill( min( scs_e3x3().at(i)-scs_eMax().at(i), scdiff ) );
		h_scs_ediff[2]    ->Fill( min( scs_e3x3().at(i)-scs_eMax().at(i), scdiff ) );
		h_scs_emaxvsediff[index]->Fill( scs_e3x3().at(i)-scs_eMax().at(i), scs_eMax().at(i));
		h_scs_emaxvsediff[2]    ->Fill( scs_e3x3().at(i)-scs_eMax().at(i), scs_eMax().at(i));

		h_scs_eratmax[index]	->Fill(eratmax);
		h_scs_eratmax[2]		->Fill(eratmax);
		h_scs_erat3x3[index]	->Fill(erat3x3);
		h_scs_erat3x3[2]		->Fill(erat3x3);
		h_scs_eratrat[index]	->Fill(eratrat);
		h_scs_eratrat[2]		->Fill(eratrat); 
		if( scetmax > 5. ) {
		  h_scs_eratratCutet[index]	->Fill(eratrat);
		  h_scs_eratratCutet[2]		->Fill(eratrat);
		}
		if( scetmax > etmaxcut ) {
		  h_scs_eratmaxCut[index]	->Fill(eratmax);
		  h_scs_eratmaxCut[2]		->Fill(eratmax);
		  h_scs_erat3x3Cut[index]	->Fill(erat3x3);
		  h_scs_erat3x3Cut[2]		->Fill(erat3x3);
		  h_scs_eratratCut[index]	->Fill(eratrat);
		  h_scs_eratratCut[2]		->Fill(eratrat);
		}
		//const float eratratmin = min( eratrat, (float)39.99 ); //this isn't necessary
		const float scetmaxmin = min( scetmax, (float)39.99 );
		h_scs_emaxvseratmaxZoom[index]		->Fill(eratmax, emax);
		h_scs_emaxvseratmaxZoom[2]			->Fill(eratmax, emax);
		h_scs_emaxvserat3x3Zoom[index]		->Fill(erat3x3, emax);
		h_scs_emaxvserat3x3Zoom[2]			->Fill(erat3x3, emax);
		h_scs_emaxvseratratZoom[index]		->Fill(eratrat, emax);
		h_scs_emaxvseratratZoom[2]			->Fill(eratrat, emax);
		h_scs_etmaxvseratmaxZoom[index]		->Fill(eratmax, scetmax);
		h_scs_etmaxvseratmaxZoom[2]			->Fill(eratmax, scetmax);
		h_scs_etmaxvserat3x3Zoom[index]		->Fill(erat3x3, scetmax);
		h_scs_etmaxvserat3x3Zoom[2]			->Fill(erat3x3, scetmax);
		h_scs_etmaxvseratratZoom[index]		->Fill(eratrat, scetmaxmin);
		h_scs_etmaxvseratratZoom[2]			->Fill(eratrat, scetmaxmin);
		const float eratmaxmod = min( max( eratmax, (float)-0.02499) , (float)0.99 );
		const float erat3x3mod = min( max( erat3x3, (float)-0.02499) , (float)0.99 );
		const float eratratmod = min( max( eratrat, (float)-0.02499) , (float)0.99 );
		h_scs_emaxvseratmax[index]		->Fill(eratmaxmod, emax);
		h_scs_emaxvseratmax[2]			->Fill(eratmaxmod, emax);
		h_scs_emaxvserat3x3[index]		->Fill(erat3x3mod, emax);
		h_scs_emaxvserat3x3[2]			->Fill(erat3x3mod, emax);
		h_scs_emaxvseratrat[index]		->Fill(eratratmod, emax);
		h_scs_emaxvseratrat[2]			->Fill(eratratmod, emax);
		h_scs_etmaxvseratmax[index]		->Fill(eratmaxmod, scetmax);
		h_scs_etmaxvseratmax[2]			->Fill(eratmaxmod, scetmax);
		h_scs_etmaxvserat3x3[index]		->Fill(erat3x3mod, scetmax);
		h_scs_etmaxvserat3x3[2]			->Fill(erat3x3mod, scetmax);
		h_scs_etmaxvseratrat[index]		->Fill(eratratmod, scetmaxmin);
		h_scs_etmaxvseratrat[2]			->Fill(eratratmod, scetmaxmin);

		const float cmetx  = evt_met()*cos( cms2.evt_metPhi() );
		const float cmety  = evt_met()*sin( cms2.evt_metPhi() );
		//select only spiking cells
		if( fabs(eratrat-1) < 0.01 && scetmax > etmaxcut ){
		  passed = true;
		  if( failed ) failed = false;
		  evt_nspikes++;
		  //totspikes++;
		  v_erat.push_back(eratmax);
		  v_emax.push_back(emax);
		  v_eta.push_back(scs_eta().at(i));
		  v_phi.push_back(scs_phi().at(i));
		  v_cmet.push_back(evt_met());
		  v_run.push_back(thisRun);
		  
		  h_scs_eta[index]->Fill(scs_eta().at(i));
		  h_scs_eta[2]    ->Fill(scs_eta().at(i));
		  h_scs_etavsphi[index]->Fill(scs_phi().at(i),scs_eta().at(i));
		  h_scs_etavsphi[2]    ->Fill(scs_phi().at(i),scs_eta().at(i));
		  h_scs_etavsphiNar[index]->Fill(scs_phi().at(i),scs_eta().at(i));
		  h_scs_etavsphiNar[2]    ->Fill(scs_phi().at(i),scs_eta().at(i));

		  h_scs_etvscmet[index]    ->Fill( scet,    cmet );
		  h_scs_etvscmet[2]        ->Fill( scet,    cmet );
		  h_scs_etmaxvscmet[index] ->Fill( scetmax, cmet );
		  h_scs_etmaxvscmet[2]     ->Fill( scetmax, cmet );
		  h_scs_etvstcmet[index]   ->Fill( scet,    tcmet );
		  h_scs_etvstcmet[2]       ->Fill( scet,    tcmet );
		  h_scs_etmaxvstcmet[index]->Fill( scetmax, tcmet );
		  h_scs_etmaxvstcmet[2]    ->Fill( scetmax, tcmet );
		  h_scs_etvspfmet[index]   ->Fill( scet,    pfmet );
		  h_scs_etvspfmet[2]       ->Fill( scet,    pfmet );
		  h_scs_etmaxvspfmet[index]->Fill( scetmax, pfmet );
		  h_scs_etmaxvspfmet[2]    ->Fill( scetmax, pfmet );
		  h_scs_hoe_spike[index]->Fill( scs_hoe().at(i) );
		  h_scs_hoe_spike[2]    ->Fill( scs_hoe().at(i) );
		  h_scs_spikeFlags[index]->Fill(scs_severitySeed().at(i));
		  h_scs_spikeFlags[2]->Fill(scs_severitySeed().at(i));
		  
		  const float tcmetx = evt_tcmet()*cos( cms2.evt_tcmetPhi() );
		  const float tcmety = evt_tcmet()*sin( cms2.evt_tcmetPhi() );
		  const float pfmetx = evt_pfmet()*cos( cms2.evt_pfmetPhi() );
		  const float pfmety = evt_pfmet()*sin( cms2.evt_pfmetPhi() );
		  const float dphicmet  = ROOT::Math::VectorUtil::DeltaPhi( LorentzVector(cmetx, cmety, 0, evt_met()), scs_pos_p4().at(i) );
		  const float dphitcmet = ROOT::Math::VectorUtil::DeltaPhi( LorentzVector(tcmetx, tcmety, 0, evt_met()), scs_pos_p4().at(i) );
		  const float dphipfmet = ROOT::Math::VectorUtil::DeltaPhi( LorentzVector(pfmetx, pfmety, 0, evt_met()), scs_pos_p4().at(i) );
		  const float scmaxetcorr = scs_eMax().at(i)*sin( scs_pos_p4().at(i).Theta() );
		  const float tcmetcorrx = scmaxetcorr*cos(scs_phi().at(i)) + tcmetx;
		  const float tcmetcorry = scmaxetcorr*sin(scs_phi().at(i)) + tcmety;
		  const float cmetcorrx  = scmaxetcorr*cos(scs_phi().at(i)) + cmetx;
		  const float cmetcorry  = scmaxetcorr*sin(scs_phi().at(i)) + cmety;
		  const float cmetcorr   = sqrt( cmetcorrx*cmetcorrx + cmetcorry*cmetcorry );
		  h_scs_dphicmet[index]			->Fill( dphicmet );
		  h_scs_dphicmet[2]				->Fill( dphicmet );
		  h_scs_phivscmetphi[index]		->Fill( scs_phi().at(i), atan2( -cmety, -cmetx ) );
		  h_scs_phivscmetphi[2]			->Fill( scs_phi().at(i), atan2( -cmety, -cmetx ) );
		  h_scs_dphitcmet[index]		->Fill( dphitcmet );
		  h_scs_dphitcmet[2]			->Fill( dphitcmet );
		  h_scs_phivstcmetphi[index]	->Fill( scs_phi().at(i), atan2( -tcmety, -tcmetx ) );
		  h_scs_phivstcmetphi[2]		->Fill( scs_phi().at(i), atan2( -tcmety, -tcmetx ) );
		  h_scs_dphipfmet[index]		->Fill( dphipfmet );
		  h_scs_dphipfmet[2]			->Fill( dphipfmet );
		  h_scs_phivspfmetphi[index]	->Fill( scs_phi().at(i), atan2( -pfmety, -pfmetx ) );
		  h_scs_phivspfmetphi[2]		->Fill( scs_phi().at(i), atan2( -pfmety, -pfmetx ) );
		  h_cmetCutCorr[index]          ->Fill( cmetcorr );
		  h_cmetCutCorr[2]              ->Fill( cmetcorr );
		  //track matching
		  int nmatch = 0;
		  for( unsigned int j=0; j<trks_outer_position().size(); j++ ) {
			float dr = ROOT::Math::VectorUtil::DeltaR(trks_outer_position().at(j), scs_pos_p4().at(i));
			if( dr < 0.1 ) {
			  //cout << "Track Match " << trks_trk_p4().at(j).pt() << "  " << scmaxetcorr << "  " << dr << "   evt " << evt_event() << endl;
			  nmatch++;
			}
		  }
		  v_trkmch.push_back(nmatch);
		  tcmetallcorr = sqrt( tcmetcorrx*tcmetcorrx + tcmetcorry*tcmetcorry );
		}
		else { //exclude spikes
		  h_scs_eratratp[index]	->Fill(eratratp);
		  h_scs_eratratp[2]   	->Fill(eratratp);
		  h_scs_etmaxvseratratpZoom[index]		->Fill(eratratp, scetmaxp);
		  h_scs_etmaxvseratratpZoom[2]			->Fill(eratratp, scetmaxp);

		  if( evt_met() > 10 ) {
			h_scs_eratratpMet[index]	->Fill(eratratp);
			h_scs_eratratpMet[2]   	    ->Fill(eratratp);
			h_scs_eratratCutMet[index]	->Fill(eratratp);
			h_scs_eratratCutMet[2]   	->Fill(eratratp);
		  }

		  if( scetmaxp > etmaxcut ) {
			h_scs_eratratpCut[index]	->Fill(eratratp);
			h_scs_eratratpCut[2]    	->Fill(eratratp);
		  }
		  if( fabs(eratratp-1) < 0.01 ){
			h_scs_etmaxnospike[index]->Fill( scetmaxp );
			h_scs_etmaxnospike[2]    ->Fill( scetmaxp );
		  }
		  if( fabs(eratratp-1) < 0.01 && scetmaxp > etmaxcut ){ //second type of spike
			h_scs_etavsphiHotp[index]		->Fill( scs_phi().at(i), scs_eta().at(i) );
			h_scs_etavsphiHotp[2]   		->Fill( scs_phi().at(i), scs_eta().at(i) );
			h_scs_pphivscmetphi[index]		->Fill( scs_phi().at(i), atan2( -cmety, -cmetx ) );
			h_scs_pphivscmetphi[2]			->Fill( scs_phi().at(i), atan2( -cmety, -cmetx ) );
			h_scs_petmaxvscmet[index] ->Fill( scetmaxp, cmet );
			h_scs_petmaxvscmet[2]     ->Fill( scetmaxp, cmet );
			h_scs_hoe_spikep[index]->Fill( scs_hoe().at(i) );
			h_scs_hoe_spikep[2]    ->Fill( scs_hoe().at(i) );
		  }

		}
		
      } //end loop on SCs
	  
      if( passed && !failed ) {
		h_cmetCut[index] ->Fill(cmet);
		h_cmetCut[2]     ->Fill(cmet);
		h_tcmetCut[index]->Fill(tcmet);
		h_tcmetCut[2]    ->Fill(tcmet);
		h_pfmetCut[index]->Fill(pfmet);
		h_pfmetCut[2]    ->Fill(pfmet);
      }
	  if( evt_nspikes > 1 ) //multi-spike events
		cout << "Multi Spike Event ****    " << evt_run() << evt_event() << endl;

	  //loop on calo towers
	  double tcmethf = evt_tcmet();
	  for( unsigned int i=0; i<twrs_emEnergy().size(); i++ ) {
		if( fabs(twrs_eta().at(i)) > 3. ) {
		  double ass = twrs_emEnergy().at(i)/(twrs_emEnergy().at(i) + twrs_hadEnergy().at(i));
		  if(  twrs_emEt().at(i) + twrs_hadEt().at(i) > 5. ) {
			h_twrs_ass[index]->Fill( ass );
			h_twrs_ass[2]->Fill( ass );
		  }
		  //double twret = (twrs_emEnergy().at(i) + twrs_hadEnergy().at(i))*sin( twrs_;
		  if( fabs(ass) > 0.98 && twrs_emEt().at(i) + twrs_hadEt().at(i) > 5. ) {
			tothfspikes++;
			const float tcmetx = evt_tcmet()*cos( cms2.evt_tcmetPhi() );
			const float tcmety = evt_tcmet()*sin( cms2.evt_tcmetPhi() );
			const float twretcorr = twrs_emEt().at(i) + twrs_hadEt().at(i); 
			const float cmetcorrx = twretcorr*cos(twrs_phi().at(i)) + tcmetx;
			const float cmetcorry = twretcorr*sin(twrs_phi().at(i)) + tcmety;
			const float cmetcorr  = sqrt( cmetcorrx*cmetcorrx + cmetcorry*cmetcorry );
			tcmethf = cmetcorr;
			tcmetallcorr = cmetcorr;
			
			h_tcmetCorr[index]->Fill( cmetcorr );
			h_tcmetCorr[2]->Fill( cmetcorr );
			h_tcmetHF[index]->Fill( evt_tcmet() );
			h_tcmetHF[2]    ->Fill( evt_tcmet() );
		  }
		}
	  }
	  h_tcmetHFCorr[index]->Fill( tcmethf );
	  h_tcmetHFCorr[2]    ->Fill( tcmethf );
	  h_tcmetAllCorr[index]->Fill( tcmetallcorr );
	  h_tcmetAllCorr[2]    ->Fill( tcmetallcorr );

      //now comes all of what was there before      
      
    }//event loop
  }//file loop

  cout << "\n\n********************SUMMARY********************" << endl;
  cout << "Total number of events: " << nEventsTotal << endl;
  cout << "Total number of events from selected runs : " << nGoodEvents << " (" << 100*(double)nGoodEvents/nEventsTotal << "%)" << endl;
  cout << "Total number of events that pass tracking cuts: " << nPassTrackingCuts << " (" << 100*(double)nPassTrackingCuts/nGoodEvents << "%)" << endl;
  cout << "Total number of events from selected runs that pass the triggers and tracking cuts: " << nPassTriggers << " (" << 100*(double)nPassTriggers/nGoodEvents << "%)" << endl;
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  cout << "Total number of unmatched electrons = " << nUnmatchedElectrons  << endl << endl; 
  
  //now print out the ecal spike info collected:
  //int crun = vrun.at(0);
  vector<int> nspikesrun;
  for( unsigned int i=0; i<v_goodRuns.size(); i++ ) {
	nspikesrun.push_back(0);
  }

  for( unsigned int i=0; i < v_erat.size(); i++){
  /*cout << " Candidate spike  " << i << endl;
    cout << " erat = " << v_erat.at(i) << endl;
    cout << " emax = " << v_emax.at(i) << endl;
    cout << " eta = " << v_eta.at(i) << endl;
    cout << " phi = " << v_phi.at(i) << endl;
    cout << " calomet = " << v_cmet.at(i) << endl;
    cout << " run number = " << v_run.at(i) << endl;
  */
	for( unsigned int j=0; j<v_goodRuns.size(); j++ ) {
	  if( (int)v_goodRuns.at(j) == v_run.at(i) ) { //match run from runlist with spike's run
		nspikesrun.at(j)++;
	  }
	}
  }

  bool printspikesummary = false;
  if( printspikesummary ) {
	int totgevts = 0;
	cout << "run & NSpikes & NEvents & Ratio \\\\\\hline" << endl;
	for( unsigned int i=0; i<v_goodRuns.size(); i++ ) {
	  cout << v_goodRuns.at(i) << " & " << nspikesrun.at(i) << " & " << nGoodEventsPerRun.at(i)
		   << " & " << double(nspikesrun.at(i))/double(nGoodEventsPerRun.at(i)) << " \\\\ " << endl;
	  totgevts += nGoodEventsPerRun.at(i);
	}
	cout << "tot spikes & " << nspikes_goodrun << "   " << v_erat.size() << "   " << totgevts << endl << endl;

	cout << "Total HF spikes " << tothfspikes << endl << endl;
  }

  cout << "good run selected nevts " << npassgoodrun << endl << endl;

  cout << "\t\tGood Run\tBad Run\n"
	   << "total scs\t\t" << nscs_goodrun << "\t" << nscs_badrun << endl
	   << "total scs 5gevet\t" << nscs5gevet_goodrun << "\t" << nscs5gevet_badrun << endl
	   << "good scs 5gevet \t" << ngoodscs5gevet_goodrun << "\t" << ngoodscs5gevet_badrun << endl
	   << "n spikes \t\t" << nspikes_goodrun << "\t" << nspikes_badrun << endl
	   << "ratio spikes to total scs  " << (float)nspikes_goodrun/(float)nscs_goodrun       << "\t" << (float)nspikes_badrun/(float)nscs_badrun << endl
	   << "ratio spikes to scs5gevet  " << (float)nspikes_goodrun/(float)nscs5gevet_goodrun << "\t" << (float)nspikes_badrun/(float)nscs5gevet_badrun << endl;

  cout << endl << endl;

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
