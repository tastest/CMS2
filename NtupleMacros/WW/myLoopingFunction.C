//now make the source file
#include <iostream>
#include <fstream>
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include <algorithm>
#include "TCanvas.h"
#include "TRegexp.h"
#include "CMS2.h"
#include "ClaudioOriginal/selections.h"

using namespace std;

enum Sample {WW, WZ, ZZ, Wjets, DYee, DYmm, DYtt, ttbar}; // signal samples
enum Hypothesis {MM, EM, EE, ALL}; // hypothesis types (em and me counted as same) and all

TChain *TreePipe (const char *cmd, const char *treename) 
{
     FILE *f = popen(cmd, "r");
     if (!f) {
	  perror("Opening pipe");
	  return 0;
     }
     TChain *c = new TChain(treename);
     int s;
     do {
	  char fname[1024];
	  s = fscanf(f, " %1024s\n", fname);
	  if (s != 1) {
	       if (s != EOF)
		    perror("scanning file list");
	  } else {
	       printf("Adding %s\n", fname);
	       c->Add(fname);
	  }
     } while (s == 1);
     if (pclose(f) == -1) 
	  perror("Closing pipe");
     return c;
}

TChain *TreePipe_glob (const char *glob, const char *treename) 
{
     std::string cmd = "ls ";
     cmd += glob;
     return TreePipe(cmd.c_str(), treename);
}

// save histograms to outfile
void saveHist(const char* filename, const char* pat="*")
{
   TList* list = gDirectory->GetList() ;
   TIterator* iter = list->MakeIterator();

   //TRegexp re(pat,kTRUE) ;

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

//-------------------------------------------------
// Auxiliary function to scan the doc line and 
// identify DY-> ee vs mm vs tt
//-------------------------------------------------
int getDrellYanType() {
  bool foundZ;
  int size = genps_id.size();
  for (int jj=0; jj<size; jj++) {
    if (genps_id.at(jj) == 23) {
      foundZ = true;
      if (jj+3 > size) {
	std::cout << 
	  "Found Z but not enough room in doc lines for leptons?" << std::endl;
        return 999;
      }
      if (abs(genps_id.at(jj+1)) == 11) return 0;  //DY->ee
      if (abs(genps_id.at(jj+1)) == 13) return 1;  //DY->mm
      if (abs(genps_id.at(jj+1)) == 15) return 2;  //DY->tautau
    }
  }
  std::cout << "Does not look like a DY event" << std::endl;
  return 999;
}

//--------------------------------------------
// Booleans for DY
//------------------------------------------
bool isDYee() {
  if (getDrellYanType() == 0) return true;
  return false;
}
bool isDYmm() {
  if (getDrellYanType() == 1) return true;
  return false;
}
bool isDYtt() {
  if (getDrellYanType() == 2) return true;
  return false;
}

// filter events by process
bool filterByProcess( enum Sample sample ) {
  switch (sample) {
  case WW: case WZ: case ZZ: 
    return true;
  case Wjets:
    return evt_CSA07Process < 11;
  case DYee: 
    return (evt_CSA07Process > 10 && evt_CSA07Process < 22 && isDYee() );
  case DYmm:
    return (evt_CSA07Process > 10 && evt_CSA07Process < 22 && isDYmm() );
  case DYtt:
    return (evt_CSA07Process > 10 && evt_CSA07Process < 22 && isDYtt() );
  case ttbar:
    return (evt_CSA07Process > 21 && evt_CSA07Process < 27);
  }
}

// filter candidates by hypothesis
Hypothesis filterByHypothesis( int candidate ) {
  switch (candidate) {
  case 0:
    return MM;
  case 1: case 2:
    return EM;
  case 3:
    return EE;
  }
}

int ScanChain( TChain* chain, enum Sample sample ) {

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
  unsigned int nEventsTotal = 0;

  const unsigned int numHypTypes = 4;  // number of hypotheses: MM, EM, EE, ALL

 // declare and create array of histograms
  const char sample_names[][1024] = { "WW", "WZ", "ZZ", "Wjets", "DYee", "DYmm", "DYtt", "ttbar" };

  // declare and initialize array of histograms
  TH1D* hist_njets[numHypTypes];

  hist_njets[MM]  = new TH1D((string(sample_names[sample])+"_njets_"+"mm" ).c_str() ,"mm"  , 4, 0, 4);
  hist_njets[EM]  = new TH1D((string(sample_names[sample])+"_njets_"+"em" ).c_str() ,"em"  , 4, 0, 4);
  hist_njets[EE]  = new TH1D((string(sample_names[sample])+"_njets_"+"ee" ).c_str() ,"ee"  , 4, 0, 4);
  hist_njets[ALL] = new TH1D((string(sample_names[sample])+"_njets_"+"all").c_str() ,"all" , 4, 0, 4);

  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  int i_permille_old = 0;

  TTree *tree = chain;

  Init(tree);  // set branch addresses for TTree tree
  
  printf("%d events in chain\n", nEventsChain);
  //Event Loop
  unsigned int nEvents = tree->GetEntries();
  for( unsigned int event = 0; event < nEvents; ++event) {
       tree->GetEntry(event);  // get entries for Event number event from branches of TTree tree
      ++nEventsTotal;

      // Progress feedback to the user
//       if ( (nEventsTotal)%1000 == 0 ) std::cout << "Processing event: " << nEventsTotal << std::endl;
      int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
      if (i_permille != i_permille_old) {
	   // xterm magic from L. Vacanvant and A. Cerri
	   printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
		  "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	   fflush(stdout);
	   i_permille_old = i_permille;
      }

      // filter by process
      if ( !filterByProcess(sample) ) continue;

      // loop over hypothesis candidates
      unsigned int nHyps = hyp_type.size();
      for( unsigned int hyps = 0; hyps < nHyps; ++hyps ) {
	// fill njet bin for final state ALL
	unsigned int numJets = ( hyp_jets_p4[hyps].size() < 3 ) ? hyp_jets_p4[hyps].size() : 3;
	hist_njets[ALL]->Fill( numJets );

	printf("numJets = %d\n", numJets);

	// filter by hypothesis candidate
	enum Hypothesis final_state  = filterByHypothesis( hyp_type[hyps] );
	
	// fill final state njet bin
	hist_njets[final_state]->Fill( numJets );
      }
  }
  

  if ( nEventsChain != nEventsTotal ) {
       printf("ERROR: number of events from files (%d) is not equal to total number"
	      " of events (%d)", nEventsChain, nEventsTotal);
  }

  return 0;
}
