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
float getTwrHFSwiss( int seedidx, bool em );
float deltaPhi(float phi1,float phi2);
int CaloTwr_ieta( int detid );
int CaloTwr_iphi( int detid );


TString ScanChain( TChain* chain, bool runningonGEN, bool requireTrackCuts = true, std::vector<unsigned int> v_goodRuns = std::vector<unsigned int>(), bool is2tev=false, int nEvents = -1) {

  //cout << "starting" << endl;
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

  /* //get rid of these--track plots
  //remove old unused histsos
*/

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

  TH2F *h_twrs_etavsphiNar[aSize];
  TH2F *h_twrs_ietavsiphiNar[aSize];
  TH2F *h_twrs_twrr9vinvE[aSize];
  
  TH1F *h_cmetCut[aSize];       //met for events with cut on erat and emax
  TH1F *h_cmet[aSize];          //met 
  TH1F *h_tcmetx[aSize];
  TH1F *h_tcmety[aSize];
  TH1F *h_cmetHFCorr[aSize];          //met 
  TH1F *h_cmetAllCorr[aSize];          //met 
  TH1F *h_cmetCutCorr[aSize];   //met for events with cut on erat and emax
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
  //TH1F *h_scs_timeseed[aSize];
  TH1F *h_scs_timeseed_all[aSize];    
  TH1F *h_scs_timeseed_spike[aSize];  
  TH1F *h_scs_timeseed_goodscs[aSize];

  TH1F *h_scs_spikeFlags[aSize];//reco flags for ECAL spikes
  TH1F *h_scs_notspikeFlags[aSize];//reco flags for ECAL !spikes
  TH1F *h_scs_allFlags[aSize];//reco flags for ECAL !spikes
  TH1F *h_scs_er4[aSize]; //compare r4 from scs to r4 from twrs

  TH1F *h_twrs_er4[aSize];              //ecal swiss, aka r4
  TH1F *h_twrs_er4Met[aSize];
  TH1F *h_twrs_er4Cut[aSize];
  TH1F *h_twrs_er4CutL[aSize];
  TH2F *h_twrs_etmaxvser4Zoom[aSize];	    
  TH2F *h_twrsec_etmaxvstcmet[aSize];
  TH2F *h_twrsec_phivstcmetphi[aSize];
  TH2F *h_twrshf_etmaxvstcmet[aSize];
  TH2F *h_twrshf_phivstcmetphi[aSize];
  TH2F *h_twrsec_etmaxvsclmet[aSize];
  TH2F *h_twrsec_phivsclmetphi[aSize];
  TH2F *h_twrshf_etmaxvsclmet[aSize];
  TH2F *h_twrshf_phivsclmetphi[aSize];
 
  TH1F *h_twrs_ass[aSize];
  TH1F *h_twrs_timeseed_all  [aSize];    
  TH1F *h_twrs_timeseed_spike[aSize];  
  TH1F *h_twrs_timeseed_goodt[aSize];
  TH1F *h_twrs_adc_all  [aSize];
  TH1F *h_twrs_adc_spike[aSize];
  TH1F *h_twrs_adc_goodt[aSize];

  float metmax = 40;
  int metbins = 40;
  float metmaxN = 30;
  int metbinsN = 30;
  
  //book histos
  //fkw's new histos come first
  for(unsigned int i = 0; i < v_prefix.size(); i++) {
	/* //get rid of unused histos
	*/ //remove old

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

    h_twrs_etavsphiNar[i] = new TH2F((v_prefix.at(i)+"twrs_etavsphiNar").c_str(), 
									("twrs_etavsphiNar" + v_title.at(i)).c_str(), 420, -TMath::Pi(), TMath::Pi(), 400, -3.0, 3.0); //.015
	h_twrs_ietavsiphiNar[i] = new TH2F((v_prefix.at(i)+"twrs_ietavsiphiNar").c_str(), 
									   ("twrs_ietavsiphiNar" + v_title.at(i)).c_str(), 400, 0., 100., 400, -50.0, 50.0);

	h_twrs_twrr9vinvE[i]  = new TH2F((v_prefix.at(i)+"twrs_twrr9vinvE").c_str(), 
									("twrs_twrr9vinvE" + v_title.at(i)).c_str(), 200, 0., 0.05, 400, 0., 0.2); 

    h_pfmetCut[i] 		= new TH1F((v_prefix.at(i)+"pfmetCut").c_str(), ("pfmetCut" + v_title.at(i)).c_str(), 100, 0.0, 50.0);
    h_pfmet[i]    		= new TH1F((v_prefix.at(i)+"pfmet").c_str(), ("pfmet" + v_title.at(i)).c_str(), 100, 0.0, 50.0);

    h_cmetCut[i]  		= new TH1F((v_prefix.at(i)+"cmetCut").c_str(), ("cmetCut" + v_title.at(i)).c_str(), metbins, 0.0, metmax);
    h_cmet[i]     		= new TH1F((v_prefix.at(i)+"cmet").c_str(), ";calo MET", metbins, 0.0, metmax);
    h_cmetCutCorr[i]  	= new TH1F((v_prefix.at(i)+"cmetCutCorr").c_str(), ("cmetCutCorr" + v_title.at(i)).c_str(), metbins, 0.0, metmax);
	h_cmetHFCorr[i]     = new TH1F((v_prefix.at(i)+"cmetHFCorr").c_str(), ";calo MET", metbins, 0.0, metmax);
	h_cmetAllCorr[i]    = new TH1F((v_prefix.at(i)+"cmetAllCorr").c_str(), ";calo MET", metbins, 0.0, metmax);

    h_tcmetCut[i] 		= new TH1F((v_prefix.at(i)+"tcmetCut").c_str(), ("tcmetCut" + v_title.at(i)).c_str(), metbins, 0.0, metmax);
    h_tcmetHF[i] 		= new TH1F((v_prefix.at(i)+"tcmetHF").c_str(), ("tcmetHF" + v_title.at(i)).c_str(), metbins, 0.0, metmax);
    h_tcmet[i]    		= new TH1F((v_prefix.at(i)+"tcmet").c_str(), ";tcMET", metbins, 0.0, metmax);
    h_tcmetx[i]    		= new TH1F((v_prefix.at(i)+"tcmetx").c_str(), ";tcMET_{X}", 40, -20., 20.);
    h_tcmety[i]    		= new TH1F((v_prefix.at(i)+"tcmety").c_str(), ";tcMET_{Y}", 40, -20., 20.);
    h_tcmetCorr[i]    	= new TH1F((v_prefix.at(i)+"tcmetCorr").c_str(), ("tcmetCorr" + v_title.at(i)).c_str(), metbins, 0.0, metmax);
    h_tcmetHFCorr[i]    = new TH1F((v_prefix.at(i)+"tcmetHFCorr").c_str(), ";tcMET", metbins, 0.0, metmax);
    h_tcmetAllCorr[i]   = new TH1F((v_prefix.at(i)+"tcmetAllCorr").c_str(), ";tcMET", metbins, 0.0, metmax);

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
	h_scs_timeseed_all[i]     = new TH1F((v_prefix.at(i)+"scs_timeseed_all").c_str(),     ("scs_timeseed_all" + v_title.at(i)).c_str(), 50, -25, 25);
	h_scs_timeseed_spike[i]   = new TH1F((v_prefix.at(i)+"scs_timeseed_spike").c_str(),   ("scs_timeseed_spike" + v_title.at(i)).c_str(), 50, -25, 25);
	h_scs_timeseed_goodscs[i] = new TH1F((v_prefix.at(i)+"scs_timeseed_goodscs").c_str(), ("scs_timeseed_goodscs" + v_title.at(i)).c_str(), 50, -25, 25);

    h_scs_allFlags[i] = new TH1F((v_prefix.at(i)+"scs_allFlags").c_str(), 
								 ("scs_allFlags" + v_title.at(i)).c_str(), 20, 0.0, 20.0);
    h_scs_spikeFlags[i] = new TH1F((v_prefix.at(i)+"scs_spikeFlags").c_str(), 
								   ("scs_spikeFlags" + v_title.at(i)).c_str(), 20, 0.0, 20.0);
    h_scs_notspikeFlags[i] = new TH1F((v_prefix.at(i)+"scs_notspikeFlags").c_str(), 
									  ("scs_notspikeFlags" + v_title.at(i)).c_str(), 20, 0.0, 20.0);
	h_scs_er4[i] = new TH1F((v_prefix.at(i)+"scs_er4").c_str(), 
							("scs_er4" + v_title.at(i)).c_str(), 220, -0.2, 2.0);

	h_twrs_ass[i] = new TH1F((v_prefix.at(i)+"twrs_ass").c_str(), ";#alpha", 201, -1.0, 1.01); 
	h_twrs_timeseed_all[i]   = new TH1F((v_prefix.at(i)+"twrs_timeseed_all").c_str(),   ("twrs_timeseed_all" + v_title.at(i)).c_str()  , 50, -25.0, 25.); 
	h_twrs_timeseed_spike[i] = new TH1F((v_prefix.at(i)+"twrs_timeseed_spike").c_str(), ("twrs_timeseed_spike" + v_title.at(i)).c_str(), 50, -25.0, 25.); 
	h_twrs_timeseed_goodt[i] = new TH1F((v_prefix.at(i)+"twrs_timeseed_goodt").c_str(), ("twrs_timeseed_goodt" + v_title.at(i)).c_str(), 50, -25.0, 25.); 
	h_twrs_adc_all  [i] = new TH1F((v_prefix.at(i)+"twrs_adc_all").c_str(), ("twrs_adc_all" + v_title.at(i)).c_str(), 10, 0, 10); 
	h_twrs_adc_spike[i] = new TH1F((v_prefix.at(i)+"twrs_adc_spike").c_str(), ("twrs_adc_spike" + v_title.at(i)).c_str(), 10, 0, 10); 
	h_twrs_adc_goodt[i] = new TH1F((v_prefix.at(i)+"twrs_adc_goodt").c_str(), ("twrs_adc_goodt" + v_title.at(i)).c_str(), 10, 0, 10); 

	h_twrs_er4[i]    = new TH1F((v_prefix.at(i)+"twrs_er4").c_str(), 
								";R4", 220, -0.2, 2.0);
	h_twrs_er4Met[i] = new TH1F((v_prefix.at(i)+"twrs_er4Met").c_str(), 
								  ("twrs_er4Met" + v_title.at(i)).c_str(), 220, -0.2, 2.0);
	h_twrs_er4Cut[i] = new TH1F((v_prefix.at(i)+"twrs_er4Cut").c_str(), 
								";R4 (Emax_T > 5)", 220, -0.2, 2.0);
	h_twrs_er4CutL[i] = new TH1F((v_prefix.at(i)+"twrs_er4CutL").c_str(), 
								";R4 (E > 1)", 220, -0.2, 2.0);
	h_twrs_etmaxvser4Zoom[i] = new TH2F((v_prefix.at(i)+"twrs_etmaxvser4").c_str(), 
										//("twrs_etmaxvser4" + v_title.at(i)).c_str(), 100, -0.2, 0.2, 20, 0.0, 20.0 );
										";R4;Tower Emax_{T}", 120, -0.2, 1, 20, 0.0, 20.0 );
	h_twrsec_etmaxvstcmet[i] = new TH2F((v_prefix.at(i)+"twrsec_etmax_vs_tcmet").c_str(),
									  ";Tower Emax_{T};tcMET", metbinsN, 0.0, metmaxN, metbinsN, 0.0, metmaxN);
	h_twrsec_phivstcmetphi[i] = new TH2F((v_prefix.at(i)+"twrsec_phivstcmetphi").c_str(), 
									   ";Tower #phi;tcMET #phi",  100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
	h_twrshf_etmaxvstcmet[i] = new TH2F((v_prefix.at(i)+"twrshf_etmax_vs_tcmet").c_str(),
										";Tower E_{T};tcMET", metbinsN, 0.0, metmaxN, metbinsN, 0.0, metmaxN); //no max in title for hf bc i take emet + hadet
	h_twrshf_phivstcmetphi[i] = new TH2F((v_prefix.at(i)+"twrshf_phivstcmetphi").c_str(), 
									   ";Tower #phi;tcMET #phi",  100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
	h_twrsec_etmaxvsclmet[i] = new TH2F((v_prefix.at(i)+"twrsec_etmax_vs_clmet").c_str(),
									  ";Tower Emax_{T};caloMET", metbinsN, 0.0, metmaxN, metbinsN, 0.0, metmaxN);
	h_twrsec_phivsclmetphi[i] = new TH2F((v_prefix.at(i)+"twrsec_phivsclmetphi").c_str(), 
									   ";Tower #phi;-calo MET #phi",  100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
	h_twrshf_etmaxvsclmet[i] = new TH2F((v_prefix.at(i)+"twrshf_etmax_vs_clmet").c_str(),
										";Tower E_{T};caloMET", metbinsN, 0.0, metmaxN, metbinsN, 0.0, metmaxN); //no max in title for hf bc i take emet + hadet
	h_twrshf_phivsclmetphi[i] = new TH2F((v_prefix.at(i)+"twrshf_phivsclmetphi").c_str(), 
									   ";Tower #phi;-caloMet #phi",  100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());


	/* //old
	 */ //old

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
    h_tcmetx[i]->TH1F::Sumw2();
    h_tcmety[i]->TH1F::Sumw2();
    h_cmetCut[i]->TH1F::Sumw2();
    h_cmetHFCorr[i]->TH1F::Sumw2();
    h_cmetAllCorr[i]->TH1F::Sumw2();
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
	h_scs_timeseed_all[i]->TH1F::Sumw2();
	h_scs_timeseed_spike[i]->TH1F::Sumw2();
	h_scs_timeseed_goodscs[i]->TH1F::Sumw2();

    h_scs_allFlags[i]->TH1F::Sumw2();
    h_scs_spikeFlags[i]->TH1F::Sumw2();
    h_scs_notspikeFlags[i]->TH1F::Sumw2();
	h_scs_er4[i]->TH1F::Sumw2();

	h_twrs_ass[i]->TH1F::Sumw2();
	h_twrs_timeseed_all[i]->TH1F::Sumw2();
	h_twrs_timeseed_spike[i]->TH1F::Sumw2();
	h_twrs_timeseed_goodt[i]->TH1F::Sumw2();
	h_twrs_er4[i]   ->TH1F::Sumw2();
	h_twrs_er4Met[i]->TH1F::Sumw2();
	h_twrs_er4Cut[i]->TH1F::Sumw2();
	h_twrs_er4CutL[i]->TH1F::Sumw2();
  }
    
  //new tree
  TTree *outTree_ ;
  std::string fileName = "FlatTree.root";
  TFile *outFile_ = new TFile(fileName.c_str(), "RECREATE");
  outFile_->cd();
  outTree_ = new TTree("T1", "Tree");

  Float_t sumet_;
  Float_t tcsumet_;
  Float_t MET_;
  Float_t METPhi_;
  Float_t tcMET_;
  Float_t tcMETPhi_;

  outTree_->Branch("sumet", &sumet_, "sumet/F");
  outTree_->Branch("tcsumet", &tcsumet_, "tcsumet/F");

  outTree_->Branch("met", &MET_, "met/F");
  outTree_->Branch("metphi", &METPhi_, "metphi/F");

  outTree_->Branch("tcmet", &tcMET_, "tcmet/F");
  outTree_->Branch("tcmetphi", &tcMETPhi_, "tcmetphi/F");


  TFile *currentFile = 0;

  //pass fail counters
  int nGoodEvents = 0;
  //int nUnmatchedElectrons = 0;
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
  int ntwrs_eta3_5gev_all   = 0;
  int ntwrs_eta3_5gev_spike = 0;
  int ntwrs_eta3_5gev_goodt = 0;
    
  int thisRun = 0;
  int npassgoodrun = 0;
  bool passgoodrun = false;
  bool havetimeseed = true; //using this for adc counts as well as time

  while ( currentFile = (TFile*)fileIter.Next() ) {
	//cout << "starting file loop" << endl;
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
	  if( runningonGEN || goodrun( evt_run(), evt_lumiBlock() ) ){
	    npassgoodrun++;
		passgoodrun = true;
	  }
	  else if( !is2tev )
		continue;

	  //cout << "Done header of event loop" << endl;

	  /*
	  //////////////////////////////////////////////////////////////////////////start block
	  //////////////////////////////////////////////////////////////////////////end block
	  */

      //if( thisRun == 123970) cout << "passed the electron block" << endl;

	  //scs superclusters info
	  const float oldmetmax = 49.9;
	  const float scdiff = 9.99;
	  const float emaxmax = 39.9;
	  const float etmaxcut = 5.0;
	  const float cmet  = min(evt_met(), metmax);
	  const float tcmet = min(evt_tcmet(), metmax);
	  const float pfmet = min(evt_pfmet(), metmax);
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
		const float scet    = min( scs_energy().at(i)*sin( scs_pos_p4().at(i).Theta() ), oldmetmax );
		const float scetmax = min( scs_eMax().at(i)  *sin( scs_pos_p4().at(i).Theta() ), oldmetmax );
		const float scetmaxp= min( (scs_eMax().at(i)+scs_e2nd().at(i))*sin( scs_pos_p4().at(i).Theta() ), oldmetmax );
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
		  //cout << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
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

		//fill for all scs except hot cell
		h_scs_hoe[index]->Fill( scs_hoe().at(i) );
		h_scs_hoe[2]    ->Fill( scs_hoe().at(i) );
		h_scs_allFlags[index]->Fill(scs_severitySeed().at(i));
		h_scs_allFlags[2]->Fill(scs_severitySeed().at(i));
		h_scs_emax[index]->Fill( min(scs_eMax().at(i), oldmetmax) );
		h_scs_emax[2]    ->Fill( min(scs_eMax().at(i), oldmetmax) );
		h_scs_etmax[index]->Fill( min( scetmax, oldmetmax) );
		h_scs_etmax[2]    ->Fill( min( scetmax, oldmetmax) );
		if( fabs(eratrat-1) < 0.01 ) {
		  h_scs_emaxCut[index]->Fill( min(scs_eMax().at(i), oldmetmax) );
		  h_scs_emaxCut[2]    ->Fill( min(scs_eMax().at(i), oldmetmax) );
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
		//make latest r4 for scs to compare to twrs
		float r4 = (scs_e1x3().at(i) + scs_e3x1().at(i) - scs_eMax().at(i))/scs_eMax().at(i);
		h_scs_er4[index]->Fill(r4);
		h_scs_er4[2]->Fill(r4);
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
		if( havetimeseed ) {
		  h_scs_timeseed_all[index]           ->Fill( scs_timeSeed().at(i) );
		  h_scs_timeseed_all[2]               ->Fill( scs_timeSeed().at(i) );
		}
		
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
		  //const float tcmetcorrx = scmaxetcorr*cos(scs_phi().at(i)) + tcmetx;
		  //const float tcmetcorry = scmaxetcorr*sin(scs_phi().at(i)) + tcmety;
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
		  if( havetimeseed ) {
			h_scs_timeseed_spike[index]   ->Fill( scs_timeSeed().at(i) );
			h_scs_timeseed_spike[2]       ->Fill( scs_timeSeed().at(i) );
		  }
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
		  //tcmetallcorr = sqrt( tcmetcorrx*tcmetcorrx + tcmetcorry*tcmetcorry );
		}
		else if( fabs(eratrat-1) > 0.01 && scetmax > etmaxcut ) { //"good scs"
		  if( havetimeseed ) {
			h_scs_timeseed_goodscs[index]           ->Fill( scs_timeSeed().at(i) );
			h_scs_timeseed_goodscs[2]               ->Fill( scs_timeSeed().at(i) );
		  }
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

		  //for for two-crystal spikes
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

	  //cout << "begin loop on calo twrs size " << twrs_emEnergy().size() << "  "
	  //	   << twrs_emMax().size() << "  "
	  //	   << twrs_emMaxTime().size() << "  "
	  //   << endl;
	  //loop on calo towers
	  int nhfspikes = 0;
	  float tcmetallcorr    = evt_tcmet();
	  float tcmetallcorrphi = evt_tcmetPhi();
	  float tcmethfcorr     = evt_tcmet();
	  float tcmethfcorrphi  = evt_tcmetPhi();
	  float cmetallcorr    = evt_met();
	  float cmetallcorrphi = evt_metPhi();
	  float cmethfcorr     = evt_met();
	  float cmethfcorrphi  = evt_metPhi();
	  for( unsigned int i=0; i<twrs_emEnergy().size(); i++ ) {
		if( havetimeseed ) {
		  //h_twrs_etavsphiNar[index]->Fill( twrs_phi()[i], twrs_eta()[i] );
		  //h_twrs_etavsphiNar[2]    ->Fill( twrs_phi()[i], twrs_eta()[i] );
		  //h_twrs_ietavsiphiNar[index]->Fill( CaloTwr_iphi( twrs_detid()[i] ), CaloTwr_ieta( twrs_detid()[i] ) );
		  //h_twrs_ietavsiphiNar[2]    ->Fill( CaloTwr_iphi( twrs_detid()[i] ), CaloTwr_ieta( twrs_detid()[i] ) );
		}

		//if( fabs(scs_eta().at(i) - 1.5) < 0.0500001 && fabs(scs_phi().at(i) - 1.57) < 0.11 //eta phi of hot cell--fkw
		if( fabs(twrs_eta().at(i) - 1.53) < 0.016 && fabs(twrs_phi().at(i) - 1.66) < 0.016 ) { //eta phi of hot cell--warren
		  //&& fabs(eratrat-1) < 0.01 && scetmax > etmaxcut ) { //eratX, emax of hot cell
		  //failed = true; //flag event with hot cell
		  continue; //skip hot cell
		}
					   
		if( fabs(twrs_eta().at(i)) > 3. ) {
		  double ass = twrs_emEnergy().at(i)/(twrs_emEnergy().at(i) + twrs_hadEnergy().at(i));
		  if(  twrs_emEt().at(i) + twrs_hadEt().at(i) > 5. ) {
			h_twrs_ass[index]->Fill( ass );
			h_twrs_ass[2]->Fill( ass );
		  }
		  //double twret = (twrs_emEnergy().at(i) + twrs_hadEnergy().at(i))*sin( twrs_;
		  if( (ass <= -0.8 || ass >= 0.99) && twrs_emEt().at(i) + twrs_hadEt().at(i) > 5. ) {
			tothfspikes++;
			nhfspikes++;
			const float twretcorr  = twrs_emEt().at(i) + twrs_hadEt().at(i); 
			const float methfcorrx  = twretcorr*cos(twrs_phi().at(i)) + tcmethfcorr* cos(tcmethfcorrphi);
			const float methfcorry  = twretcorr*sin(twrs_phi().at(i)) + tcmethfcorr* sin(tcmethfcorrphi);
			const float metallcorrx = twretcorr*cos(twrs_phi().at(i)) + tcmetallcorr*cos(tcmetallcorrphi);
			const float metallcorry = twretcorr*sin(twrs_phi().at(i)) + tcmetallcorr*sin(tcmetallcorrphi);
			tcmethfcorr  = sqrt( methfcorrx*methfcorrx   + methfcorry* methfcorry );
			tcmetallcorr = sqrt( metallcorrx*metallcorrx + metallcorry*metallcorry );
			tcmetallcorrphi = atan2( metallcorry, metallcorrx ); 

			const float cmethfcorrx  = twretcorr*cos(twrs_phi().at(i)) + cmethfcorr* cos(cmethfcorrphi);
			const float cmethfcorry  = twretcorr*sin(twrs_phi().at(i)) + cmethfcorr* sin(cmethfcorrphi);
			const float cmetallcorrx = twretcorr*cos(twrs_phi().at(i)) + cmetallcorr*cos(cmetallcorrphi);
			const float cmetallcorry = twretcorr*sin(twrs_phi().at(i)) + cmetallcorr*sin(cmetallcorrphi);
			cmethfcorr  = sqrt( cmethfcorrx*cmethfcorrx   + cmethfcorry* cmethfcorry );
			cmetallcorr = sqrt( cmetallcorrx*cmetallcorrx + cmetallcorry*cmetallcorry );
			cmetallcorrphi = atan2( cmetallcorry, cmetallcorrx ); 

			float tmptwret = min( twretcorr, float(metmaxN-0.01) );
			h_twrshf_phivstcmetphi[index]->Fill( twrs_phi().at(i), atan2( -evt_tcmet()*sin(evt_tcmetPhi()), -evt_tcmet()*cos(evt_tcmetPhi()) ) );
			h_twrshf_phivstcmetphi[2]    ->Fill( twrs_phi().at(i), atan2( -evt_tcmet()*sin(evt_tcmetPhi()), -evt_tcmet()*cos(evt_tcmetPhi()) ) );
			h_twrshf_etmaxvstcmet[index]->Fill( tmptwret, min( evt_tcmet(), float(metmaxN-0.01) ) );
			h_twrshf_etmaxvstcmet[2]    ->Fill( tmptwret, min( evt_tcmet(), float(metmaxN-0.01) ) );
			h_twrshf_phivsclmetphi[index]->Fill( twrs_phi().at(i), atan2( -evt_met()*sin(evt_metPhi()), -evt_met()*cos(evt_metPhi()) ) );
			h_twrshf_phivsclmetphi[2]    ->Fill( twrs_phi().at(i), atan2( -evt_met()*sin(evt_metPhi()), -evt_met()*cos(evt_metPhi()) ) );
			h_twrshf_etmaxvsclmet[index]->Fill( tmptwret, min( evt_met(), float(metmaxN-0.01) ) );
			h_twrshf_etmaxvsclmet[2]    ->Fill( tmptwret, min( evt_met(), float(metmaxN-0.01) ) );
			h_tcmetCorr[index]->Fill( tcmethfcorr );
			h_tcmetCorr[2]->Fill( tcmethfcorr );
			h_tcmetHF[index]->Fill( evt_tcmet() );
			h_tcmetHF[2]    ->Fill( evt_tcmet() );
		  }

		  //PJ plot
		  float twrr9 = getTwrHFSwiss( i, true )/twrs_emEnergy()[i]; //"true" for em--r9 = s9/s1
		  if( twrr9 > 0. ) {
			h_twrs_twrr9vinvE[index]->Fill( 1/twrs_emEnergy()[i], twrr9 );
			h_twrs_twrr9vinvE[2]    ->Fill( 1/twrs_emEnergy()[i], twrr9 ); 
		  }


		} //end if eta > 3
		else {
		  //cout << "starting eta < 3 i = " << i << endl;

		  //r4 is (the four neighbors not including max) over max
		  //definition of spike is now r4<0.05 && et>5
		  const float emmaxet = twrs_emMax().at(i)/cosh(twrs_eta().at(i));
		  const float r4      = (twrs_emSwiss().at(i) - twrs_emMax().at(i))/twrs_emMax().at(i);

		  //temporary measure to check if border events are causing discontinuity
		  //if( fabs(twrs_eta().at(i)) > 1.3 ) continue; //barrel only //&& fabs(twrs_eta().at(i)) < 1.7

		  //overflow positive
		  h_twrs_etmaxvser4Zoom[index]->Fill(r4, min(emmaxet,(float)19.9));
		  h_twrs_etmaxvser4Zoom[2]    ->Fill(r4, min(emmaxet,(float)19.9));
		  h_twrs_er4[index]          ->Fill(r4);
		  h_twrs_er4[2]              ->Fill(r4);
		  if( twrs_emMax().at(i) > 1. ) { //loose cut on energy, not et
			h_twrs_er4CutL[index]   ->Fill(r4);
			h_twrs_er4CutL[2]       ->Fill(r4);
		  }
		  if( emmaxet > 5. ) {
			h_twrs_er4Cut[index]   ->Fill(r4);
			h_twrs_er4Cut[2]       ->Fill(r4);
		  }
		  if( evt_tcmet() > 15 ) {
			h_twrs_er4Met[index]     ->Fill(r4);
			h_twrs_er4Met[2]         ->Fill(r4);
		  }

		  //this if for spikes
		  if( emmaxet > 5. && r4 < 0.05 ) {
			const float metallcorrx = emmaxet*cos(twrs_phi().at(i)) + tcmetallcorr*cos(tcmetallcorrphi);
			const float metallcorry = emmaxet*sin(twrs_phi().at(i)) + tcmetallcorr*sin(tcmetallcorrphi);
			tcmetallcorr = sqrt( metallcorrx*metallcorrx + metallcorry*metallcorry );
			const float cmetallcorrx = emmaxet*cos(twrs_phi().at(i)) + cmetallcorr*cos(cmetallcorrphi);
			const float cmetallcorry = emmaxet*sin(twrs_phi().at(i)) + cmetallcorr*sin(cmetallcorrphi);
			cmetallcorr = sqrt( cmetallcorrx*cmetallcorrx + cmetallcorry*cmetallcorry );
			cmetallcorrphi = atan2( cmetallcorry, cmetallcorrx );
			
			float tmpemmaxet = min( emmaxet, float(metmaxN-0.01) );
			h_twrsec_phivstcmetphi[index]->Fill( twrs_phi().at(i), atan2( -evt_tcmet()*sin(evt_tcmetPhi()), -evt_tcmet()*cos(evt_tcmetPhi()) ) );
			h_twrsec_phivstcmetphi[2]    ->Fill( twrs_phi().at(i), atan2( -evt_tcmet()*sin(evt_tcmetPhi()), -evt_tcmet()*cos(evt_tcmetPhi()) ) );
			h_twrsec_etmaxvstcmet[index]->Fill( tmpemmaxet, min( evt_tcmet(), float(metmaxN-0.01) ) );
			h_twrsec_etmaxvstcmet[2]    ->Fill( tmpemmaxet, min( evt_tcmet(), float(metmaxN-0.01) ) );
			h_twrsec_phivsclmetphi[index]->Fill( twrs_phi().at(i), atan2( -evt_met()*sin(evt_metPhi()), -evt_met()*cos(evt_metPhi()) ) );
			h_twrsec_phivsclmetphi[2]    ->Fill( twrs_phi().at(i), atan2( -evt_met()*sin(evt_metPhi()), -evt_met()*cos(evt_metPhi()) ) );
			h_twrsec_etmaxvsclmet[index]->Fill( tmpemmaxet, min( evt_met(), float(metmaxN-0.01) ) );
			h_twrsec_etmaxvsclmet[2]    ->Fill( tmpemmaxet, min( evt_met(), float(metmaxN-0.01) ) );
		  }

		  //this if for fancy time/adc stuff
		  if( emmaxet > 5. && havetimeseed ) {
			//cout << "r l e: " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
			if( twrs_em3x3().at(i) == 0 || !isfinite(twrs_em3x3().at(i)) 
				|| !isfinite(twrs_emMaxTime().at(i)) ){
			  cout << "found BAD twr " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
			  continue;
			}

			//const float r9 = twrs_emMax().at(i)/twrs_em3x3().at(i);
			ntwrs_eta3_5gev_all++;
			h_twrs_timeseed_all[index]->Fill( twrs_emMaxTime().at(i) );
			h_twrs_timeseed_all[2]->Fill( twrs_emMaxTime().at(i) );

			//pulse shape
			for (unsigned int s = 0; s < 10; ++s) {
			  h_twrs_adc_all[index]->Fill(s, cms2.twrs_emMaxEcalMGPASampleADC()[i][s]);
			  h_twrs_adc_all[2]    ->Fill(s, cms2.twrs_emMaxEcalMGPASampleADC()[i][s]);
			}

			//if( fabs(r9-1) < 0.01 ) { //spike
			if( r4 < 0.05 ) {
			  h_twrs_timeseed_spike[index]->Fill( twrs_emMaxTime().at(i) );
			  h_twrs_timeseed_spike[2]->Fill( twrs_emMaxTime().at(i) );
			  ntwrs_eta3_5gev_spike++;
			  //pulse shape
			  for (unsigned int s = 0; s < 10; ++s) {
				h_twrs_adc_spike[index]->Fill(s, cms2.twrs_emMaxEcalMGPASampleADC()[i][s]);
				h_twrs_adc_spike[2]    ->Fill(s, cms2.twrs_emMaxEcalMGPASampleADC()[i][s]);
			  }
			}
			else {
			  h_twrs_timeseed_goodt[index]->Fill( twrs_emMaxTime().at(i) );
			  h_twrs_timeseed_goodt[2]->Fill( twrs_emMaxTime().at(i) );
			  ntwrs_eta3_5gev_goodt++;
			  //pulse shape
			  for (unsigned int s = 0; s < 10; ++s) {
				h_twrs_adc_goodt[index]->Fill(s, cms2.twrs_emMaxEcalMGPASampleADC()[i][s]);
				h_twrs_adc_goodt[2]    ->Fill(s, cms2.twrs_emMaxEcalMGPASampleADC()[i][s]);
			  }
			}
		  } //end emmaxet + time if
		  //cout << "end eta < 3" << endl;
		}//end eta if (else)

	  } //end loop on towers
	  if( nhfspikes > 1 )
		cout << "More one HF spike  " << nhfspikes << endl;

	  h_tcmetHFCorr[index]->Fill( tcmethfcorr );
	  h_tcmetHFCorr[2]    ->Fill( tcmethfcorr );
	  h_tcmetAllCorr[index]->Fill( tcmetallcorr );
	  h_tcmetAllCorr[2]    ->Fill( tcmetallcorr );
	  h_cmetHFCorr[index]->Fill( cmethfcorr );
	  h_cmetHFCorr[2]    ->Fill( cmethfcorr );
	  h_cmetAllCorr[index]->Fill( cmetallcorr );
	  h_cmetAllCorr[2]    ->Fill( cmetallcorr );
	  h_tcmetx[index]    ->Fill( tcmetallcorr*cos(tcmetallcorrphi) );
	  h_tcmetx[2]        ->Fill( tcmetallcorr*cos(tcmetallcorrphi) );
	  h_tcmety[index]    ->Fill( tcmetallcorr*sin(tcmetallcorrphi) );
	  h_tcmety[2]        ->Fill( tcmetallcorr*sin(tcmetallcorrphi) );

	  sumet_   = evt_sumet();
	  tcsumet_ = evt_tcsumet();
	  MET_     = cmetallcorr;
	  METPhi_  = cmetallcorrphi;
	  tcMET_   = tcmetallcorr;
	  tcMETPhi_= tcmetallcorrphi;

	  outFile_->cd();
      outTree_->Fill();


      //now comes all of what was there before      
      //cout << "end of event loop" << endl;
    }//event loop
  }//file loop

  //Avi -- Save tree
  outFile_->cd();
  outTree_->Write();
  outFile_->Close();
  delete outFile_; 


  if( havetimeseed ) {
	//scale pulse shape (adc) hist by 1/ntwrs_eta3_5gev to get avg pulse
	for( unsigned int i=0; i<aSize; i++) {
	  h_twrs_adc_all[i]	->Scale( 1./(float)ntwrs_eta3_5gev_all );
	  h_twrs_adc_spike[i]	->Scale( 1./(float)ntwrs_eta3_5gev_spike );
	  h_twrs_adc_goodt[i]	->Scale( 1./(float)ntwrs_eta3_5gev_goodt );
	}
  }

  cout << "\n\n********************SUMMARY********************" << endl;
  cout << "Total number of events: " << nEventsTotal << endl;
  cout << "Total number of events from selected runs : " << nGoodEvents << " (" << 100*(double)nGoodEvents/nEventsTotal << "%)" << endl;
  cout << "Total number of events that pass tracking cuts: " << nPassTrackingCuts << " (" << 100*(double)nPassTrackingCuts/nGoodEvents << "%)" << endl;
  cout << "Total number of events from selected runs that pass the triggers and tracking cuts: " << nPassTriggers << " (" << 100*(double)nPassTriggers/nGoodEvents << "%)" << endl;
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  //cout << "Total number of unmatched electrons = " << nUnmatchedElectrons  << endl << endl; 
  
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

  cout << endl << "n     twrs eta < 3 et > 5 : " << ntwrs_eta3_5gev_all << endl;
  cout << "n good twrs eta < 3 et > 5 : " << ntwrs_eta3_5gev_goodt << endl;
  cout << "n spike twrs eta < 3 et > 5 : " << ntwrs_eta3_5gev_spike << endl;
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


//fn to reproduce what PJ calls 'S9': the (em or had) of a tower plus the em+had of the four immediate neighbors
//input: index in twrs block of seed, bool for em or had
//this is ONLY for HF--to deal with borders:
//  do not calculate this for one tower on high and low eta border:
//   first tower is 2.853<|eta|<2.964, last tower is 4.889<|eta|<5.191
//algo: loop over all towers, for each, 
//  if abs delta tower phi is < 10 degrees (one tower width in phi), compare eta. 
//    If delta eta is < two tower sizes (0.175*2 = 0.35) then add energy (em+had) to result
//  and then reverse roles of eta/phi

// range of ieta is -41 <= ieta <= 41, EXCLUDING zero
// range of iphi is   1 <= iphi <= 72 in barrel (-20 <= ieta <= 20),
//                    1 <= iphi <= 71, ODD NUMBERS ONLY in endcap      : 21 <= abs(ieta) <= 39
//                    3 <= iphi <= 71, EVERY FOURTH ONLY in far forward: 40 <= abs(ieta) <= 41
//                    ie, 3, 7, 11, 15

float getTwrHFSwiss( int seedidx, bool em ) {

  if( twrs_eta()[seedidx] < 2.964 || twrs_eta()[seedidx] > 4.889 ) //see above
	return -1;

  float result = ( em ? twrs_emEnergy()[seedidx] : twrs_hadEnergy()[seedidx] );
  unsigned int added = 0;
  const int seedieta = CaloTwr_ieta( twrs_detid()[seedidx] );
  const int seediphi = CaloTwr_iphi( twrs_detid()[seedidx] );
  const int maxiphibarrel = 72;
  const int maxiphiendcap = 71;

  //cout << evt_run() << "  " << evt_lumiBlock() << "  " << evt_event() << endl;

  for( unsigned int i=0; i<twrs_eta().size(); i++ ) {
	if( int(i) == seedidx )
	  continue;

	bool phineighbors = false;
	bool etaneighbors = false;
	const int thisieta = CaloTwr_ieta( twrs_detid()[i] );
	const int thisiphi = CaloTwr_iphi( twrs_detid()[i] );
	const int maxiphi = (abs(thisieta) <= 20 ? maxiphibarrel : maxiphiendcap );
	int diffiphi = 1; //barrel
	if( abs(thisieta) > 20 && abs(thisieta) <= 39 )
	  diffiphi = 2;
	else if( abs(thisieta) > 39 )
	  diffiphi = 4;

	if( (seediphi == 1       && thisiphi == maxiphi) || //take account of periodicity of phi
		(seediphi == maxiphi && thisiphi == 1      ) ||
		abs(thisiphi-seediphi) == diffiphi ) { //abs bc doesn't matter if + or - diffiphi
	  phineighbors = true;
	}

	if( abs(thisieta-seedieta) == 1 ||
		(thisieta == 1  && seedieta == -1) ||
		(thisieta == -1 && seedieta ==  1) )
	  etaneighbors = true;
		
	if( ( seediphi == thisiphi && etaneighbors ) || //same phi, neighbors in eta
		( phineighbors && seedieta == thisieta ) ) { //same eta, neighbors in phi
	  result += twrs_emEnergy()[i] + twrs_hadEnergy()[i];
	  //if( added == 0 )
		//cout << "seed " << seedieta << "  " << seediphi << endl;
	  added++;
	  //cout << "added " << thisieta << "  " << thisiphi << endl;
	}

	/* ///////////////////////////////////////////
	//DON"T USE THIS---IT"S BROKEN
	//float dphi = ROOT::Math::VectorUtil::DeltaPhi( LorentzVector(, , , ),  );
	float dphi = deltaPhi( twrs_phi()[seedidx], twrs_phi()[i] );
	float deta = twrs_eta()[seedidx] - twrs_eta()[i];

	if( fabs(dphi) * 180./TMath::Pi() < 10. && fabs(deta) < 0.35 ) { //same phi, neighbors in eta
	  result += twrs_emEnergy()[i] + twrs_hadEnergy()[i]; 
	  added++;
	}
	else if( fabs(deta) < 0.175 && fabs(dphi) * 180./TMath::Pi() < 20. ) { //same eta, neighbors in phi
	  result += twrs_emEnergy()[i] + twrs_hadEnergy()[i];
	  added++;
	}
	*/
  }

  //cout << endl;
  if( added > 4 )
	cout << "added too many towers--fix" << endl;

  return result;

}


float deltaPhi(float phi1,float phi2){

  float deltaphi=phi1-phi2;
  if(deltaphi> acos(-1.)) deltaphi=2*acos(-1.)-deltaphi;
  if(deltaphi<-acos(-1.)) deltaphi=2*acos(-1.)+deltaphi;
  return deltaphi;

}

//port of code from CaloTowerDetId.h
int CaloTwr_ieta( int detid ) {
  int zside = (detid & 0x2000) ? 1 : -1 ;
  //warning: multiplication has higher precendnce than bitwise anding
  return zside * ((detid >> 7) & 0x3f);
}

int CaloTwr_iphi( int detid ) {
  return detid & 0x7F;
}
