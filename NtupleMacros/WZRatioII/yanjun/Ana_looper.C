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

//#include "CMS2.h"
//#include "branches.h"



#include "CORE/CMS2.h"
#include "CORE/selections.cc"
#include "CORE/utilities.cc"
#include "Tools/tools.cc"
#include "Ana_looper.h"

CMS2 cms2;
//using namespace tas;
using namespace std;

void Ana_looper::bookHistos(char* sample, int nchannels, int nhistsets){
 
  for (int i_ch=0; i_ch<nchannels; i_ch++) {
    for (int j_hist=0; j_hist<nhistsets; j_hist++) {  
      els_pt[i_ch][j_hist] = book1DHist(Form("%s_%s_%s%i%s%i",sample,"elsPt","Ch",i_ch,"H",j_hist),Form("%s_%s_%s%i%s%i",sample,"elsPt","Ch",i_ch,"H",j_hist),50,0,100,"Electron Pt", "Events", kBlue);  

      njets[i_ch][j_hist] = book1DHist(Form("%s_%s_%s%i%s%i",sample,"nJets","Ch",i_ch,"H",j_hist),Form("%s_%s_%s%i%s%i",sample,"nJets","Ch",i_ch,"H",j_hist),10,0,10,"Number of Jets", "Events", kBlue);  
      
    }
  }
}

double Ana_looper::Trans_W_Mass(TVector3& tcMET, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >&  p4){
  double mass = sqrt((tcMET.Pt()+p4.Et())*(tcMET.Pt()+p4.Et())
		     -(tcMET.Pt()*cos(tcMET.Phi())+ p4.Px())*(tcMET.Pt()*cos(tcMET.Phi())+ p4.Px())
		     -(tcMET.Pt()*sin(tcMET.Phi())+ p4.Py())*(tcMET.Pt()*cos(tcMET.Phi())+ p4.Py()));
  return mass;
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
     
      float weight = kFactor*cms2.evt_scale1fb();  //scale to 1 fb-1
      int channel = -1;
      int hist  = 0;
      int nels =   cms2.els_p4().size();
      int nmus =   cms2.mus_p4().size();
      int ntrks =  cms2.trks_trk_p4().size();
      int njpts =  cms2.jpts_p4().size();
      vector<int>  good_els_idx;
      vector<int>  good_mus_idx;
      vector<int>  good_trks_idx;
      vector<int>  good_jpts_idx;
      int ngoodels = 0;
      int ngoodmus = 0;
      int ngoodtrks = 0;
      int ngoodjpts = 0;
      ///loop over electrons
      for(int i_els=0; i_els<nels; i_els++){
        if(cms2.els_p4()[i_els].pt() <20) continue;
        if(!goodElectronIsolated(i_els,1)) continue;
	//if(!conversionElectron(i_els)) continue; 
	
        good_els_idx.push_back(i_els);

      }
      ngoodels = good_els_idx.size();

      
      ///loop over muons
      for(int i_mus=0; i_mus<nmus; i_mus++){
        if(cms2.mus_p4()[i_mus].pt()< 20) continue;
        if(!goodMuonIsolated(i_mus)) continue;
	good_mus_idx.push_back(i_mus);
      }
      ngoodmus = good_mus_idx.size();
      ///loop over tracks
      for(int i_trks=0; i_trks<ntrks; i_trks++){
        bool pass_trk = true;
        if(cms2.trks_trk_p4()[i_trks].pt()< 10) continue;
        if(!passTrackIsolation(i_trks)) continue;
        for(int i_els=0; i_els<ngoodels; i_els++){
          if ( (TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.els_p4()[i_els],cms2.trks_trk_p4()[i_trks])) < 0.15) )
          {
            pass_trk = false;
            break;
          }
        }
	for(int i_mus=0; i_mus<ngoodmus; i_mus++){
          if ( (TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.mus_p4()[i_mus],cms2.trks_trk_p4()[i_trks])) < 0.15) )
          {
            pass_trk = false;
            break;
          }
        }
        if(pass_trk){
	  good_trks_idx.push_back(i_trks);
        }
      }
      
      ngoodtrks = good_trks_idx.size();
      ///loop over jpt jets
      for ( int i_jpts=0; i_jpts < njpts; i_jpts++) {
        bool pass_jpt = true;
        if ( cms2.jpts_p4()[i_jpts].Et() < 20. ) continue;
        if ( TMath::Abs(cms2.jpts_p4()[i_jpts].eta()) > 2.4 ) continue;
        for(int i_els=0; i_els< ngoodels; i_els++){
          if ( (TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.els_p4()[i_els],cms2.jpts_p4()[i_jpts])) < 0.4) )
          {
            pass_jpt = false;
            break;
          }
        }

        for(int i_mus=0; i_mus<ngoodmus; i_mus++){
          if ( (TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.mus_p4()[i_mus],cms2.jpts_p4()[i_jpts])) < 0.4) )
          {
            pass_jpt = false;
            break;
          }
        }
        if(pass_jpt) good_jpts_idx.push_back(i_jpts);
      }
      ngoodjpts = good_jpts_idx.size();
      
      

      /////event selection; just a example, you can change it;
      int idword = 0;

      TVector3 tcmet;
      tcmet.SetPtEtaPhi(cms2.evt_tcmet(), 0, cms2.evt_tcmetPhi());
      if(tcmet.Pt() > 20) idword |= kPassWMETBit;
      if(tcmet.Pt() < 30) idword |= kPassZMETBit;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vec;
    
      double w_trans_mass;
      if(ngoodels==1 &&  ngoodmus==0){
	idword |= kOneEleBit;
	w_trans_mass =  Trans_W_Mass(tcmet, cms2.els_p4()[good_els_idx.at(0)]);
	
	if(w_trans_mass > 40 && w_trans_mass <100)  idword |= kPassWMassBit; 
	// else printf("%s, %f \n", "W trans mass  ",  w_trans_mass);
	

	bool z_hyp=false;
	for(int i=0; i<ngoodtrks; i++){
	 
	  vec = cms2.trks_trk_p4()[good_trks_idx.at(i)]+cms2.els_p4()[good_els_idx.at(0)];
	  if(inZmassWindow(vec.mass()) )  z_hyp = true; 
	  //       else printf("%s, %f \n", "Z mass  ",  vec.mass());
	}
	if(z_hyp ==false) idword |= kPassZVetoBit;
      }
      
      if(ngoodels==0 &&  ngoodmus==1){
	idword |= kOneMuonBit;
	w_trans_mass =  Trans_W_Mass(tcmet, cms2.mus_p4()[good_mus_idx.at(0)]);
	
	if(w_trans_mass > 40 && w_trans_mass <100)  idword |= kPassWMassBit; 
	// else printf("%s, %f \n", "W trans mass  ",  w_trans_mass);
	

	bool z_hyp=false;
	for(int i=0; i<ngoodtrks; i++){
	 
	  vec = cms2.trks_trk_p4()[good_trks_idx.at(i)]+cms2.mus_p4()[good_mus_idx.at(0)];
	  if(inZmassWindow(vec.mass()) )  z_hyp = true; 
	  //       else printf("%s, %f \n", "Z mass  ",  vec.mass());
	}
	if(z_hyp ==false) idword |= kPassZVetoBit;
      }
   
      if (ngoodels==2 &&  ngoodmus==0){
	idword |= kTwoEleBit;
	for(int i=0; i<ngoodels; i++){
	  for(int j=i+1; j<ngoodels; j++){
	    vec = cms2.els_p4()[good_els_idx.at(i)]+cms2.els_p4()[good_els_idx.at(j)];
	    if(inZmassWindow(vec.mass()) ) {
	      idword |= kPassZMassBit;
	    }
	    if((cms2.els_charge().at(good_els_idx.at(i))*cms2.els_charge().at(good_els_idx.at(j))) == -1){
	      idword |= kOppChargeBit;
	    }
	  }
	}
      } 
      
      if (ngoodels==0 &&  ngoodmus==2){
	idword |= kTwoMuonBit;
	for(int i=0; i<ngoodmus; i++){
	  for(int j=i+1; j<ngoodmus; j++){
	    vec = cms2.mus_p4()[good_mus_idx.at(i)]+cms2.mus_p4()[good_mus_idx.at(j)];
	    if(inZmassWindow(vec.mass()) ) {
	      idword |= kPassZMassBit;
	    }
	    if((cms2.mus_charge().at(good_mus_idx.at(i))*cms2.mus_charge().at(good_mus_idx.at(j))) == -1){
	      idword |= kOppChargeBit;
	    }
	  }
	}
      }

      int we_id = kOneEleBit | kPassWMETBit | kPassWMassBit | kPassZVetoBit ;
      int wm_id = kOneMuonBit | kPassWMETBit | kPassWMassBit | kPassZVetoBit ;
      int zee_id = kTwoEleBit  | kPassZMETBit | kPassZMassBit  | kOppChargeBit;
      int zmm_id = kTwoMuonBit  | kPassZMETBit | kPassZMassBit  | kOppChargeBit;

      
      if((idword & we_id)==we_id    )   channel =0; 
      if((idword & wm_id)==wm_id    )   channel =1; 
      if((idword & zee_id)==zee_id  )   channel =2; 
      if((idword & zmm_id)==zmm_id  )   channel =3; 
      
      if(channel >= 0){
	njets[channel][0]->Fill( ngoodjpts, weight);
	
	for(int i_els=0; i_els<ngoodels; i_els++){
	  if(channel < 0 || hist <0 || channel > NCHANNELS || hist > NHISTS)
	    std::cout << "ERROR: number of channels or histograms exceeds the maximum values" <<std::endl;
	  els_pt[channel][0]->Fill(cms2.els_p4().at(i_els).pt(), weight);
	}
      }

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
