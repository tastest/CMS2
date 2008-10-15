#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "LooperBase.h"
#include "Sample.h"
#include "utilities.h"
#include "CMS2.h"

LooperBase::LooperBase (Sample s) : sample(s)
{

}

LooperBase::~LooperBase ()
{

}

void LooperBase::Begin ()
{
     printf("Processing %s\n", sample.name.c_str());
}

uint64 LooperBase::Loop ()
{
     Begin();
     TChain *chain = sample.chain;
     uint64 nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
     uint64 nEventsTotal = 0;
     memset(hypos_total_n, 0, sizeof(hypos_total_n));
     memset(hypos_total_weight, 0, sizeof(hypos_total_weight));

     // clear list of duplicates
     already_seen.clear();
     int duplicates_total_n = 0;
     double duplicates_total_weight = 0;

     int i_permille_old = 0;
     // file loop
     TObjArray *listOfFiles = chain->GetListOfFiles();
     TIter fileIter(listOfFiles);
     while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
	  // need to call TFile::Open(), since the file is not
	  // necessarily a plain TFile (TNetFile, TDcacheFile, etc)
	  TFile *f = TFile::Open(currentFile->GetTitle()); 
	  TTree *tree = (TTree*)f->Get("Events");
	  if (tree == 0) { 
	       fprintf(stderr, "***HELP ME!  Events TREE IS NOT ACCESSIBLE IN %s!***\n",
		       currentFile->GetTitle());
	       delete f;
	       continue;
	  }

	  cms2.Init(tree);  // set branch addresses for TTree tree

	  TStopwatch t;
	  //Event Loop
	  uint64 nEvents = tree->GetEntries();
// 	  nEvents = std::min(nEvents, 1000u);
	  for( uint64 event = 0; event < nEvents; ++event) {
	       cms2.GetEntry(event);  // get entries for Event number event from branches of TTree tree
	       ++nEventsTotal;
	       
	       if (cms2.trks_d0().size() == 0)
		    continue;
	       DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.trks_d0()[0], 
					   cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
	       if (is_duplicate(id)) {
		    duplicates_total_n++;
		    duplicates_total_weight += cms2.evt_scale1fb();
		    continue;
	       }
	       // Progress feedback to the user
//       if ( (nEventsTotal)%1000 == 0 ) std::cout << "Processing event: " << nEventsTotal << std::endl;
	       int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
	       if (i_permille != i_permille_old) {
		    // xterm magic from L. Vacavant and A. Cerri
		    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
			   "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
		    fflush(stdout);
		    i_permille_old = i_permille;
	       }
	       
	       // filter by process
	       if ( !filterByProcess(sample.process) ) continue;
	       
	       // call the event function
	       Event();

	       // call the dilepton candidate function
	       for (unsigned int i_hyp = 0, nHyps = cms2.hyp_type().size();
		    i_hyp < nHyps; ++i_hyp ) {
		    Dilep(i_hyp);
	       }
	       // trilepton
	       for (unsigned int i_hyp = 0, nHyps = cms2.hyp_trilep_bucket().size();
		    i_hyp < nHyps; ++i_hyp ) {
		    Trilep(i_hyp);
	       }
	       // quadlepton
	       for (unsigned int i_hyp = 0, nHyps = cms2.hyp_quadlep_bucket().size();
		    i_hyp < nHyps; ++i_hyp ) {
		    Quadlep(i_hyp);
	       }
	       
	  }
	  t.Stop();
	  printf("Real time: %llu events / %f s = %e event/s\n", nEvents, 
		 t.RealTime(), nEvents / t.RealTime());
	  printf("CPU time: %llu events / %f s = %e event/s\n", nEvents, 
		 t.CpuTime(), nEvents / t.CpuTime());
	  delete f;
     }
     if ( nEventsChain != nEventsTotal ) {
	  printf("ERROR: number of events from files (%llu) is not equal to total number"
		 " of events (%llu)\n", nEventsChain, nEventsTotal);
     }
     
//      printf("Total candidate count (ee em mm all): %llu %llu %llu %llu.  Total weight %f %f %f %f\n",   
// 	    hypos_total_n[0], hypos_total_n[1], hypos_total_n[2], hypos_total_n[3], 
// 	    hypos_total_weight[0], hypos_total_weight[1], hypos_total_weight[2], hypos_total_weight[3]);
//      printf("Total duplicate count: %d.  Total weight %f\n",   
// 	    duplicates_total_n, duplicates_total_weight);
     End();
     return nEventsTotal;
}
