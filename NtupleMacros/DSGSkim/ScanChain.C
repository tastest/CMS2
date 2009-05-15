#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TEventList.h"
#include "TBranch.h"

#include "CMS2.h"
CMS2 cms2;
#include "../CORE/selections.cc"
#include "../CORE/utilities.cc"

using namespace tas;

int ScanChain( TChain* chain, const char* outputname) {

  // output file and tree
  TFile *output = new TFile(outputname,"RECREATE");
  TTree *newtree = 0;

  // additional branches: declaration
  TBranch *anotherBranch = 0;
  TBranch *DSGBucketBranch = 0;
  
  // additional variables to be filled into an additional branch into the new tree
  int anotherVariable;
  std::vector<int> *anotherVector = new std::vector<int>;

  std::vector<int> *DSGBucketVector = new std::vector<int>;

  // list of files
  TObjArray *listOfFiles = chain->GetListOfFiles();

  // events to process
  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;

  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  bool first = true;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");

    // for the first file, clone the tree
    if ( first ) {
      newtree = chain->CloneTree(0);
      newtree->SetDirectory(output);
      first = false;

      // book branches
      anotherBranch = newtree->Branch("anotherbranch",&anotherVariable, "another_branch/I");
      newtree->SetAlias("another_branch","anotherbranch");

      anotherBranch = newtree->Branch("anothervector",&anotherVector);
      newtree->SetAlias("another_vector","anothervector");

      DSGBucketBranch = newtree->Branch("DSGBucketvector",&DSGBucketVector);
      newtree->SetAlias("DSGBucket_vector","DSGBucketvector");

    }

    // init
    cms2.Init(newtree);
    cms2.Init(tree);

    // Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents; ++event) {
      ++nEventsTotal;
      if ( nEventsTotal%10000 == 0 ) {
	cout << "Event: " << nEventsTotal << endl;
      }

      //      if(nEventsTotal > 1000) continue;

      // reset additional variables
      anotherVariable = 0;
      anotherVector->clear();

      DSGBucketVector->clear();

      // VERSION 0
      // inital DSGBucket definition for first strawman table
      //                     |njet<=1 | 2<=njet<=4 | njet>=5      |
      //=====================|========|============|==============|=
      //         ss + (of/sf)|   11   |     12     |      13      |
      //---------------------|--------|------------|--------------|-
      //         ss - (of/sf)|   21   |     22     |      23      |
      //---------------------|--------|------------|--------------|-
      //         os of       |   31   |     32     |      33      |
      //---------------------|--------|------------|--------------|-
      //         os sf       |   41   |     42     |      43      |
      //---------------------|--------|------------|--------------|-

      // VERSION 1- chosen
      // inital DSGBucket definition for first strawman table
      //                     |njet<=1 | 2<=njet<=4 | njet>=5      |
      //=====================|========|============|==============|=
      //         ss+ (of/sf) |   1    |      2     |       3      |
      //---------------------|--------|------------|--------------|-
      //         ss- (of/sf) |   4    |      5     |       6      |
      //---------------------|--------|------------|--------------|-
      //         os of       |   7    |      8     |       9      |
      //---------------------|--------|------------|--------------|-
      //         os sf       |   10   |     11     |      12      |
      //---------------------|--------|------------|--------------|-
      
      cms2.GetEntry(event);
      cms2.LoadAllBranches();

      // cut on njets >= 4 - placeholder for skim TCut here!
      if ( evt_njets() < 4 ) continue;

      // fill additional variables
      anotherVariable = evt_njets()-2;
      anotherVector->push_back(evt_njets()-4);
      anotherVector->push_back(evt_njets()-3);
      anotherVector->push_back(evt_njets()-2);
      anotherVector->push_back(evt_njets()-1);

      // loop the hypothesis, save a DSG bucket for each hyp
      for (unsigned int i_hyp = 0, nHyps = hyp_type().size(); i_hyp < nHyps; ++i_hyp ) {
        //dertermine DSGBucket
        bool sameFlavour     = false;
        bool sameSignPlus    = false;
        bool sameSignMinus   = false;
        bool oppSign         = false;

        if( hyp_lt_id()[i_hyp] > 0 && hyp_ll_id()[i_hyp] > 0 ) {
          sameSignPlus  = true;
          sameSignMinus = false;
          oppSign       = false;
        }
        else if( hyp_lt_id()[i_hyp] < 0 && hyp_ll_id()[i_hyp] < 0 ) {
          sameSignPlus  = false;
          sameSignMinus = true;
          oppSign       = false;
        }
        else {
          sameSignPlus  = false;
          sameSignMinus = false;
          oppSign       = true;
        }
        
       if( 
          ( TMath::Abs(hyp_lt_id()[i_hyp]) == 13 && TMath::Abs(hyp_ll_id()[i_hyp] ) == 13 )  ||
          ( TMath::Abs(hyp_lt_id()[i_hyp]) == 11 && TMath::Abs(hyp_ll_id()[i_hyp] ) == 11 )  
          ) {
         sameFlavour = true;
       }
       else {
         sameFlavour = false;
       }

       // these cuts need to be read from TCuts file!
        if( sameSignPlus && evt_njets()<=1 )                                             DSGBucketVector->push_back(1);
        else if( sameSignPlus && evt_njets()>=2 && evt_njets()<=4 )                      DSGBucketVector->push_back(2);
        else if( sameSignPlus && evt_njets()>=5 )                                        DSGBucketVector->push_back(3);
                                                                                         
        else if( sameSignMinus && evt_njets()<=1 )                                       DSGBucketVector->push_back(4);
        else if( sameSignMinus && evt_njets()>=2 && evt_njets()<=4 )                     DSGBucketVector->push_back(5);
        else if( sameSignMinus && evt_njets()>=5 )                                       DSGBucketVector->push_back(6);

        else if( oppSign && !sameFlavour && evt_njets()<=1 )                             DSGBucketVector->push_back(7);
        else if( oppSign && !sameFlavour && evt_njets()>=2 && evt_njets()<=4 )           DSGBucketVector->push_back(8);
        else if( oppSign && !sameFlavour && evt_njets()>=5 )                             DSGBucketVector->push_back(9);

        else if( oppSign && sameFlavour && evt_njets()<=1 )                              DSGBucketVector->push_back(10);
        else if( oppSign && sameFlavour && evt_njets()>=2 && evt_njets()<=4 )            DSGBucketVector->push_back(11);
        else if( oppSign && sameFlavour && evt_njets()>=5 )                              DSGBucketVector->push_back(12);

        else { // fallback, should never happen!
          std::cout<<"ALARM! unknown bucket in DSGSkim - check!"<<std::endl;
          DSGBucketVector->push_back(-999);
        }
      } // end hypothesis loop

      // fill the new tree
      newtree->Fill();

    }
  }

  output->cd();
  newtree->Write();
  delete output;

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  
  return 0;
}

