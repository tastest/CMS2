#include <iostream>
#include <iomanip>


/* Usage:
   root [0] .L Ana_looper.C++
   root [1] TFile *_file0 = TFile::Open("merged_ntuple.root")
   root [2] TChain *chain = new TChain("Events")
   root [3] chain->Add("merged_ntuple.root")
   root [4] Ana_looper a 
   root [5] a.ScanChain(chain) // will give the same results
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



//#include "CORE/CMS2.h"
//#include "CORE/selections.cc"
//#include "CORE/utilities.cc"
#include "Tools/tools.cc"
#include "Ana_looper.h"

CMS2 cms2;
//using namespace tas;
using namespace std;

double Ana_looper::dRBetweenVectors(LorentzVector v1, LorentzVector v2) {
    double deta = v1.eta() - v2.eta();
    double dphi = fabs(v1.phi() - v2.phi());
    if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
    return sqrt(deta*deta+dphi*dphi);
}
std::pair<float, float> Ana_looper::getConversionInfo(LorentzVector trk1_p4, 
					  int trk1_q, float trk1_d0, 
					  LorentzVector trk2_p4,
					  int trk2_q, float trk2_d0,
					  float bField) {
  
  
  double tk1Curvature = -0.3*bField*(trk1_q/trk1_p4.pt())/100.;
  double rTk1 = fabs(1./tk1Curvature);
  double xTk1 = (1./tk1Curvature - trk1_d0)*cos(trk1_p4.phi());
  double yTk1 = (1./tk1Curvature - trk1_d0)*sin(trk1_p4.phi());
    
  double tk2Curvature = -0.3*bField*(trk2_q/trk2_p4.pt())/100.;
  double rTk2 = fabs(1./tk2Curvature);
  double xTk2 = (1./tk2Curvature - trk2_d0)*cos(trk2_p4.phi());
  double yTk2 = (1./tk2Curvature - trk2_d0)*sin(trk2_p4.phi());
	 
  double dist = sqrt(pow(xTk1-xTk2, 2) + pow(yTk1-yTk2 , 2));
  dist = dist - (rTk1 + rTk2);

  double dcot = 1/tan(trk1_p4.theta()) - 1/tan(trk2_p4.theta());

  return make_pair(dist, dcot);
  
}

bool Ana_looper::isconversionElectron09(int elIdx) {

  for(unsigned int tkIdx = 0; tkIdx < cms2.trks_trk_p4().size(); tkIdx++) {
    if(dRBetweenVectors(cms2.els_trk_p4()[elIdx], cms2.trks_trk_p4()[tkIdx]) > 0.5)
      continue;
    //skip the electron's track
    if(cms2.els_trkidx()[elIdx] == tkIdx && cms2.els_trkshFrac()[elIdx] > 0.45)
      continue;
    //ship non-opp sign candidates
    if(cms2.trks_charge()[tkIdx] + cms2.els_charge()[elIdx] != 0)
      continue;
    
    std::pair<float, float> temp = getConversionInfo(cms2.els_trk_p4()[elIdx], cms2.els_charge()[elIdx], cms2.els_d0()[elIdx], 
						     cms2.trks_trk_p4()[tkIdx], cms2.trks_charge()[tkIdx], cms2.trks_d0()[tkIdx],
						     cms2.evt_bField());
    
    if(fabs(temp.first) < 0.02 && fabs(temp.second) < 0.02)
      return true;
    
  }//track loop
  
  return false;
  
}
bool Ana_looper::isconversionElectron_PIXHIT(int ele_index) {
 // true if electron is a conversion electron

 if(cms2.els_p4()[ele_index].eta()>1.47){
   if ( cms2.els_valid_pixelhits()[ele_index] == 0) return true;
   else if ( cms2.els_valid_pixelhits()[ele_index] ==1 ||cms2.els_valid_pixelhits()[ele_index] ==2)
     {
	if((cms2.els_layer1_det()[ele_index] == 1) &&  cms2.els_layer1_layer()[ele_index]>1 ) return true;
	else if((cms2.els_layer1_det()[ele_index] == 2) &&  cms2.els_layer1_charge()[ele_index]>40000 ) return true;
     }
 }
 return false;
}


void Ana_looper::bookHistos(char* sample, int nchannels, int nhistsets){
 
  for (int i_ch=0; i_ch<nchannels; i_ch++) {
    for (int j_hist=0; j_hist<nhistsets; j_hist++) {  
      els_pt[i_ch][j_hist] = book1DHist(Form("%s_%s_%s%i%s%i",sample,"elsPt","Ch",i_ch,"H",j_hist),Form("%s_%s_%s%i%s%i",sample,"elsPt","Ch",i_ch,"H",j_hist),50,0,100,"Electron Pt", "Events", kBlue);  
      els_eta[i_ch][j_hist] = book1DHist(Form("%s_%s_%s%i%s%i",sample,"elsEta","Ch",i_ch,"H",j_hist),Form("%s_%s_%s%i%s%i",sample,"elsEta","Ch",i_ch,"H",j_hist),100,-3.2,3.2,"Electron Eta", "Events", kBlue);  

      njets[i_ch][j_hist] = book1DHist(Form("%s_%s_%s%i%s%i",sample,"nJets","Ch",i_ch,"H",j_hist),Form("%s_%s_%s%i%s%i",sample,"nJets","Ch",i_ch,"H",j_hist),10,0,10,"Number of Jets", "Events", kBlue);  
      
    }
  }
}




int Ana_looper::ScanChain( TChain* chain, int nEvents ,char* sample, float kFactor , int prescale) {

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=0;
  if(nEvents==-1) 
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  
  
  bookHistos(sample,NCHANNELS, NHISTS);
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
     
      //float weight = kFactor*cms2.evt_scale1fb();  //scale to 1 fb-1
      float weight =  1.0;
      int channel = -1;
      int hist  = 0;
      int nels =   cms2.els_p4().size();
      
      vector<int>  good_els_idx;
      int ngoodels = 0;
      channel = 0;
      ///loop over electrons
      for(int i_els=0; i_els<nels; i_els++){
        if(cms2.els_p4()[i_els].pt() <20) continue;
	if(TMath::Abs(cms2.els_p4()[i_els].eta()) > 2.4) continue;
	els_eta[channel][0]->Fill(cms2.els_p4()[i_els].eta(), weight);
	
	if(cms2.els_egamma_tightId()[i_els] !=1)continue;
	if(cms2.els_closestMuon()[i_els] != -1)continue;
	els_eta[channel][1]->Fill(cms2.els_p4()[i_els].eta(), weight);
 // 	if((cms2.els_tkIso03()[i_els] + cms2.els_hcalIso03()[i_els])/max(cms2.els_p4()[i_els].pt(),20.) >= 0.1) continue;  //ecal iso is buggy in release 3_1_0;
//  	els_eta[channel][2]->Fill(cms2.els_p4()[i_els].eta(), weight);
       	
	if(TMath::Abs(cms2.els_d0corr()[i_els]) >= 0.020) continue;
	els_eta[channel][3]->Fill(cms2.els_p4()[i_els].eta(), weight);


	if(cms2.els_n_inner_layers()[i_els] >1) continue;
	els_eta[channel][4]->Fill(cms2.els_p4()[i_els].eta(), weight);
	
	if(isconversionElectron09(i_els))continue;
	els_eta[channel][5]->Fill(cms2.els_p4()[i_els].eta(), weight);
	
	if(isconversionElectron_PIXHIT(i_els))continue;
	els_eta[channel][6]->Fill(cms2.els_p4()[i_els].eta(), weight);
        

	good_els_idx.push_back(i_els);

      }
      ngoodels = good_els_idx.size();
     
     

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
