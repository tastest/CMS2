#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include "TChain.h"
#include "histtools.C"
#include "MyScanChain.cc"
#include <iostream>
void doData_long() {
  
  //     //
  //     // the looper
  //     //
  //     gSystem->Load("libTree.so");
  //     gSystem->Load("libPhysics.so");
  //     gSystem->Load("libEG.so");
  //     gSystem->Load("libMathCore.so");
  //     gSystem->Load("libCMS2NtupleMacrosCORE.so");
  //     gSystem->Load("libCMS2NtupleMacrosLooper.so");
  //     //    gROOT->ProcessLine(".L histtools.C+");
  
  //
  // output file for histograms
  //
  
  int NSTEPS = 20;
  
  TString infiles[NSTEPS];
  
  infiles[0]  = "/tas05/disk00/slava77/reltestdata/CMSSW_3_5_6-cms2-data/merged*1324*.root";
  infiles[1]  = "/tas05/disk00/slava77/reltestdata/CMSSW_3_5_6-cms2-data/merged*1325*.root";
  infiles[2]  = "/tas05/disk00/slava77/reltestdata/CMSSW_3_5_6-cms2-data/merged*1326*.root";
  infiles[3]  = "/tas05/disk00/slava77/reltestdata/CMSSW_3_5_6-cms2-data/merged*1327*.root";
  infiles[4]  = "/tas05/disk00/slava77/reltestdata/CMSSW_3_5_6-cms2-data/merged*1329*.root";
  infiles[5]  = "/tas05/disk00/slava77/reltestdata/CMSSW_3_5_6-cms2-data/merged*1331*.root";
  infiles[6]  = "/tas05/disk00/slava77/reltestdata/CMSSW_3_5_6-cms2-data/merged*1332*.root";
  infiles[7]  = "/tas05/disk00/slava77/reltestdata/CMSSW_3_5_6-cms2-data/merged*1333*.root";
  infiles[8]  = "/tas05/disk00/slava77/reltestdata/CMSSW_3_5_6-cms2-data/merged*1334*.root";
  infiles[9]  = "/tas05/disk00/slava77/reltestdata/CMSSW_3_5_6-cms2-data/merged*1335*.root";
  infiles[10] = "/tas03/disk03/slava77/reltestdata/CMSSW_3_5_7-cms2-data/merged_ntuple_1343*.root";
  infiles[11] = "/tas03/disk03/slava77/reltestdata/CMSSW_3_5_7-cms2-data/merged_ntuple_1344*.root";
  infiles[12] = "/tas03/disk03/slava77/reltestdata/CMSSW_3_5_7-cms2-data/merged_ntuple_1345*.root";
  infiles[13] = "/tas03/disk03/slava77/reltestdata/CMSSW_3_5_7-cms2-data/merged_ntuple_1346*.root";
  infiles[14] = "/tas03/disk03/slava77/reltestdata/CMSSW_3_5_7-cms2-data/merged_ntuple_1347*.root";
  infiles[15] = "/tas03/disk03/slava77/reltestdata/CMSSW_3_5_7-cms2-data/merged_ntuple_1335*.root";
  infiles[16] = "/tas03/disk03/slava77/reltestdata/CMSSW_3_5_7-cms2-data/merged_ntuple_1336*.root";
  infiles[17] = "/tas03/disk03/slava77/reltestdata/CMSSW_3_5_7-cms2-data/merged_ntuple_1337*.root";
  infiles[18] = "/tas03/disk03/slava77/reltestdata/CMSSW_3_5_7-cms2-data/merged_ntuple_1338*.root";
  infiles[19] = "/tas03/disk03/slava77/reltestdata/CMSSW_3_5_7-cms2-data/merged_ntuple_1339*.root";
  
  MyScanChain *looper = 0;
  TChain *chain_whunt_skim = 0;
  
  //  for(int splitNr = 0; splitNr < NSTEPS; ++splitNr) {
   std::cout<<"WARNING, temp restricted doData_long to subset of runs."<<std::endl;
   for(int splitNr = 18; splitNr < NSTEPS; ++splitNr) {
    
    TString filename = "histos_data_";
    filename+=splitNr;
    filename.Append(".root");
    
    std::cout<<"New output will be: "<<filename<<std::endl;
    
    looper = new MyScanChain();
    // DATA
    chain_whunt_skim = new TChain("Events");
    //    chain_whunt_skim->Add("/tas03/disk01/whunt/skim/emuskim_*.root");
    
    std::cout<<"Running on files: "<<infiles[splitNr]<<std::endl;
    
    chain_whunt_skim->Add(infiles[splitNr]);
    looper->ScanChain(true, "whunt", chain_whunt_skim);
    //
    // write histograms
    // 
    const char* outFile = filename.Data();
    hist::saveHist(outFile); 
    hist::deleteHistos();
    //
    // tidy up
    //
    looper = 0;
    chain_whunt_skim = 0;
  }
}

