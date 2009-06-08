#ifndef BRANCHES_H
#define BRANCHES_H
#include <vector>
#include "TFile.h"
#include "TTree.h"

float wmt_;
int lepid_;
int nJPTS_;
TFile *outFile_;
TTree *outTree_;
float iso_;
float zmass_;

TH1F* h_mupt1   = new TH1F("h_mupt1"  , "h_mupt1"  , 25, 0 , 100);
TH1F* h_mupt2   = new TH1F("h_mupt2"  , "h_mupt2"  , 25, 0 , 100);
TH1F* h_muett1  = new TH1F("h_muett1" , "h_muett1" , 50, -3, 3  );
TH1F* h_muett2  = new TH1F("h_muett2" , "h_muett2" , 50, -3, 3  );
TH1F* h_muiso   = new TH1F("h_muiso"  , "h_muiso"  , 50, 0 , 1  );
TH1F* h_mumet   = new TH1F("h_mumet"  , "h_mumet"  , 60, 0 , 120);
TH1F* h_mumass  = new TH1F("h_mumass" , "h_mumass" , 40, 0 , 200);

TH1F* h_elpt1   = new TH1F("h_elpt1"  , "h_elpt1"  , 25, 0 , 100);
TH1F* h_elpt2   = new TH1F("h_elpt2"  , "h_elpt2"  , 25, 0 , 100);
TH1F* h_elett1  = new TH1F("h_elett1" , "h_elett1" , 50, -3, 3  );
TH1F* h_elett2  = new TH1F("h_elett2" , "h_elett2" , 50, -3, 3  );
TH1F* h_eliso   = new TH1F("h_eliso"  , "h_eliso"  , 50, 0 , 1  );
TH1F* h_elmet   = new TH1F("h_elmet"  , "h_elmet"  , 60, 0 , 120);
TH1F* h_elmass  = new TH1F("h_elmass" , "h_elmass" , 40, 0 , 200);

void InitSkimmedTree() {

   outFile_ = TFile::Open("skimmednTuple.root","RECREATE");
   outFile_->cd();
   outTree_ = new TTree("Events", "");

   //book the branches
   outTree_->Branch("wmt_", &wmt_, "wmt_/F");
   outTree_->SetAlias("w_mt", "wmt_");
   outTree_->Branch("lepid_", &lepid_, "lepid_/I");
   outTree_->SetAlias("lep_id", "lepid_");
   outTree_->Branch("iso_", &iso_, "iso_/F");
   outTree_->SetAlias("iso", "iso_");
   outTree_->Branch("zmass_", &zmass_, "zmass_/F");
   outTree_->SetAlias("zmass", "zmass_");   
} 
#endif
