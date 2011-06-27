#include <iostream>
#include <fstream>
#include "TChain.h"
#include "TTree.h"
#include "TChainElement.h"

using namespace std;

TChain* ChainFromText(const char* filename){

  // check files exists
  ifstream infile(filename);
  if( ! infile.is_open() ){
    cout << "Could not open file." << endl;
    //gSystem->Exit(1);
	exit(1);
  } 

  // chain to return   
  TChain* chain = new TChain("Events");

  // add files to chain
  cout << endl << "Adding Files to chain..." << endl << endl;
  while( ! infile.eof() ){
	string line;
	if( getline (infile,line) ){ 
	  cout << "\t" << Form("%s/*.root",line.c_str()) << endl;
	  int added = chain->Add( Form("%s/*.root",line.c_str()) );
	  // check there are root files in the path
	  if( added < 1 ){
		cout << "Error: No root files found... exiting." << endl;
		//gSystem->Exit(1);
		exit(1);
	  }
	  cout << "\t\t-> " << added << " files added." << endl;


	  TChain *tempChain = new TChain("Events");
	  tempChain->Add( Form("%s/*.root",line.c_str()) );
	  TString dset;
	  tempChain->SetBranchAddress("evt_dataset", &dset);
	  //Int_t nentries = (Int_t)tempChain->GetEntries();
	  //for(Int_t i=0; i < 1; i++){
	  //  tempChain->GetEntry(i); 
	  //  cout << dset.Data() << endl;
	  //}



	  // quick & Dirty
	  //TChain *tempChain = new TChain("Events");
	  //tempChain->Add( Form("%s/*.root",line) );
	  //tempChain->SetScanField(1);
	  //tempChain->Scan("evt_dataset","", "colsize=100") > dataset.temp;
	}
  }
  infile.close();

  cout << chain->GetEntries() << " total events in " << chain->GetListOfFiles()->GetEntries() << " files." << endl;
  cout << endl;
  return chain;

}
