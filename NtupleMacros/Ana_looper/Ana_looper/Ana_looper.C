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
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "CMS2.h"
//#include "branches.h"



//#include "CORE/CMS2.h"
// #include "CORE/selections.cc"
// #include "CORE/utilities.cc"
//#include "Tools/tools.cc"
#include "Ana_looper.h"

CMS2 cms2;
using namespace tas;
using namespace std;

void Ana_looper::bookHistos(char* sample, int nchannels, int nhistsets){
 
  for (int i_ch=0; i_ch<nchannels; i_ch++) {
    for (int j_hist=0; j_hist<nhistsets; j_hist++) {  
      //  els_pt[i_ch][j_hist] = book1DHist(Form("%s_%s_%s%i%s%i",sample,"elsPt","Ch",i_ch,"H",j_hist),Form("%s_%s_%s%i%s%i",sample,"elsPt","Ch",i_ch,"H",j_hist),50,0,100,"Electron Pt", "Events", kBlue);  

      // njets[i_ch][j_hist] = book1DHist(Form("%s_%s_%s%i%s%i",sample,"nJets","Ch",i_ch,"H",j_hist),Form("%s_%s_%s%i%s%i",sample,"nJets","Ch",i_ch,"H",j_hist),10,0,10,"Number of Jets", "Events", kBlue);  
      trkIso03[i_ch][j_hist] = new TH1F(Form("%s_%s_%s%i%s%i",sample,"ran_trkIso03_mu","",i_ch,"",j_hist), Form("%s_%s_%s%i%s%i",sample,"ran_trkIso03_mu","",i_ch,"",j_hist), 50,-0.5, 2);
      //      trkIso03[i_ch][j_hist] = new TH1F("ran_trkIso03_mu","ran_trkIso03_mu", 50,-0.5,2);
    }
  }
}

Ana_looper::Ana_looper(){
  ran_trksp4_=0;
  trks_trk_p4_=0;
  ran_isoTrk03_mu_=0;
 
}
Ana_looper::~Ana_looper(){
}


int Ana_looper::ScanChain( TChain* chain, int nEvents ,char* sample, float kFactor , int prescale,  std::string  skimFilePrefix) {
  if(skimFilePrefix != "")outFile_ = TFile::Open(string( skimFilePrefix+ "_skimmednTuple.root").c_str(),"RECREATE");
  else outFile_ = TFile::Open("skimmednTuple.root", "RECREATE");
  outFile_->cd();
  outTree_ = new TTree("Events", "");
  //book the branches
 
  outTree_->Branch("ran_trksp4",   "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >",   &ran_trksp4_);
  outTree_->Branch("trks_trk_p4",   "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >",   &trks_trk_p4_);
  outTree_->Branch("ran_trkIso03_mu",   "std::vector<float>",   &ran_isoTrk03_mu_);
  


  
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
     
      ran_isoTrk03_mu_->clear();

      ran_trksp4_= &(cms2.ran_trksp4());
      trks_trk_p4_= &(cms2.trks_trk_p4());
      int ntrks = 0;
      int npsuedo = 0;
      int channel =0;
      ntrks = cms2.trks_trk_p4().size();
      npsuedo = cms2.ran_trksp4().size();
      for(int i_psuedo=0; i_psuedo<npsuedo; i_psuedo++){
      ///loop over tracks
	double ptSum =0.0;
	
	for(int i_trks=0; i_trks<ntrks; i_trks++){
	
	  // note that pseudo directions are with respect to a 0, 0, 0 vertex

	  double this_pt  = cms2.trks_trk_p4().at(i_trks).pt();
	  if ( this_pt < 0 ) 
	    continue;
	  //upstream filter require at least 1 non-fake vertex. 
	  
	  // double this_dz = cms2.trks_z0().at(i_trks)- cms2.vtxs_position().at(0).z()  ;

	  double this_dz = cms2.trks_z0corr().at(i_trks)+ cms2.evt_bsp4().z()- cms2.vtxs_position().at(0).z()  ;
	  
	 
	  //std::cout <<cms2.vtxs_position().at(0).z() <<std::endl;
	  if (fabs( this_dz )> 0.2 )
	  // if (fabs( this_dz )> 1 )
	    continue;
	   if (fabs(cms2.trks_d0corr().at(i_trks) ) > 0.1   )
	 
	    continue;
	  // double dr = ROOT::Math::VectorUtil::DeltaR(cms2.trks_vertex_p4().at(i_trks),cms2.ran_trksp4().at(i_psuedo)) ;
	  double dr = ROOT::Math::VectorUtil::DeltaR(cms2.trks_trk_p4().at(i_trks),cms2.ran_trksp4().at(i_psuedo)) ;
	  if ( fabs(dr) < 0.3 && fabs(dr) >= 0.01  )
	  
	    {

	      ptSum += this_pt;
	    }
	  
	}
	trkIso03[0][0]->Fill(ptSum);
	ran_isoTrk03_mu_->push_back(ptSum);
      }
     

      outTree_->Fill();  
      ++nEventsTotal;
      if(nEventsTotal%10000 ==0)std::cout << "number of events processed " << nEventsTotal<<std::endl;
    }
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    std::cout << "number of events processed " << nEventsTotal<<std::endl;
  }
  //  delete trkIso03[0][0];
  outFile_->cd();
  outTree_->Write();
  outFile_->Close();
  return 0;
}
