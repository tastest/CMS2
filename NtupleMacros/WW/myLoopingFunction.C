//now make the source file
#include <iostream>
#include <fstream>
#include <vector>

#include "TChain.h"
#include "TFile.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include <algorithm>
#include "TCanvas.h"
#include "TRegexp.h"

#include "CMS2.h"

void saveHist(const char* filename, const char* pat="*")
{
   TList* list = gDirectory->GetList() ;
   TIterator* iter = list->MakeIterator();

   //   TRegexp re(pat,kTRUE) ;
   
   Int_t test = 1;

   TRegexp re(pat, 1) ;

   TFile outf(filename,"RECREATE") ;
   while(TObject *obj=iter->Next()) {
      if (TString(obj->GetName()).Index(re)>=0) {
         obj->Write() ;
         cout << "." ;
         cout.flush() ;
      }
   }
   cout << endl ;
   outf.Close() ;

   delete iter ;
}


int ScanChain( TChain* chain) {

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
  unsigned int nEventsTotal = 0;
  unsigned int nHypsTotal   = 0;

  const unsigned int jetBuckets  = 4;  // count njet bins 0, 1, 2, 3+
  const unsigned int numHypTypes = 4;  // hypothesis type mumu, emu/mue, ee, all

  enum { MM, EM, ALL, EE }; // hypothesis types (em and me counted as same) and all

  // declare and create array of histograms
  TH1I* hist_njets[numHypTypes];
  //for(unsigned int h = 0; h < numHypTypes; ++h ) {
  //  hist_njets[h] = new TH1I(Form("h%I",h),Form("h%I",h), 0, 4, 4);
  //  hist_njets[h]->Sumw2();
  //}

  hist_njets[MM]  = new TH1I("mm"  ,"mm"  , 4, 0, 4);
  hist_njets[EM]  = new TH1I("em"  ,"em"  , 4, 0, 4);
  hist_njets[EE]  = new TH1I("ee"  ,"ee"  , 4, 0, 4);
  hist_njets[ALL] = new TH1I("all" ,"all" , 4, 0, 4);

  // jet count per bin
  unsigned int jetCounter[numHypTypes][jetBuckets];
  for( unsigned int k = 0; k < numHypTypes; ++k ) {
    for( unsigned int i = 0; i < jetBuckets; ++i ) {
      jetCounter[k][i] = 0;
    }
  }

  // CUTS
  const double jetMaxEtaCut = 3. ;  // absolute value eta must be less than this value for jet to be counted
  const double jetMinEtCut  = 15.;  // minimum transverse energy jet must have to be counted

  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    Init(tree);  // set branch addresses for TTree tree
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents; ++event) {
      GetEntry(event);  // get entries for Event number event from branches of TTree tree
      ++nEventsTotal;

      // Progress feedback to the user
      if ((nEventsTotal)%1000 == 0) std::cout << "Processing event: " << nEventsTotal << std::endl;

      unsigned int numJets = 0;

      // loop over jets
      for ( unsigned int jet = 0; jet < jets_p4.size(); ++jet ) {
	if( fabs( jets_p4[jet].Eta() ) < jetMaxEtaCut && jets_p4[jet].Et() > jetMinEtCut ) {
	  ++numJets;
	}
      }

      // loop over dilepton hypotheses and store numJets in appropriate bin for each hypothesis
      for( unsigned int hyps = 0; hyps < hyp_type.size(); ++hyps, ++nHypsTotal ) {
	if(hyp_type[hyps] == MM) {
	  if(numJets < jetBuckets) ++jetCounter[MM][numJets];
	  else ++jetCounter[MM][jetBuckets-1];
	  hist_njets[MM]->Fill( ( (numJets < 3) ? numJets : 3 ) );
	}
	else if(hyp_type[hyps] == EE) {
	  if(numJets < jetBuckets) ++jetCounter[EE][numJets];
	  else ++jetCounter[EE][jetBuckets-1];
	  hist_njets[EE]->Fill( ( (numJets < 3) ? numJets : 3 ) );
	}
	else {
	  if(numJets < jetBuckets) ++jetCounter[EM][numJets];
	  else ++jetCounter[EM][jetBuckets-1];
	  hist_njets[EM]->Fill( ( (numJets < 3) ? numJets : 3 ) );	  
	}
      }

      // store number of jets in ALL once for each event
      if(numJets < jetBuckets) ++jetCounter[ALL][numJets];
      else ++jetCounter[ALL][jetBuckets-1];
      hist_njets[ALL]->Fill( ( (numJets < 3) ? numJets : 3 ) );
    }
  }

  for( unsigned int l = 0; l < numHypTypes; ++l ) {
    for( unsigned int j = 0; j < jetBuckets; ++j ) {
      if( l == MM )      cout << "Njets, mm: "  << j << "entries: " << jetCounter[l][j] << endl;
      else if( l == EM ) cout << "Njets, em: "  << j << "entries: " << jetCounter[l][j] << endl;
      else if( l == EE ) cout << "Njets, ee: "  << j << "entries: " << jetCounter[l][j] << endl;
      else               cout << "Njets, all: " << j << "entries: " << jetCounter[l][j] << endl;
    }
    cout << endl;
  }

  cout << "Dilepton Hypotheses: " << nHypsTotal << endl << endl;

  //TCanvas* c1 = new TCanvas;
  //c1->Divide(2,2);
  //c1->cd(1);
  //hist_njets[EE]->Draw();
  //c1->cd(2);
  //hist_njets[MM]->Draw();
  //c1->cd(3);
  //hist_njets[EM]->Draw();
  //c1->cd(4);
  //hist_njets[ALL]->Draw();

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  return 0;
}
