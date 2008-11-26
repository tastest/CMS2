#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "LooperBase.h"
#include "Sample.h"
#include "tools.h"
#include "../CORE/CMS2.h"

LooperBase::LooperBase (Sample s, cuts_t c, const char *fname) : 
     sample(s), cuts(c)
{
     if (fname != 0 && strlen(fname) != 0) {
	  logfile = fopen(fname, "a");
	  if (logfile == 0)
	       perror("opening log file");
     } else {
	  logfile = stdout;
     }
}

LooperBase::~LooperBase ()
{

}

double LooperBase::Weight (int i_hyp)
{
     return cms2.evt_scale1fb() * sample.kFactor;
}

uint64 LooperBase::Loop ()
{
     printf("Processing %s\n", sample.name.c_str());
     // change to histogram directory
     TDirectory *old_gDirectory = gDirectory;
     gDirectory = histo_directory;
     // book histos
     BookHistos();
     // change back to current directory
     gDirectory = old_gDirectory;
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
	  if (f == 0) {
	       fprintf(stderr, "File %s could not be opened.  Aborting looper.\n", 
		       currentFile->GetTitle());
	       exit(1);
	  }
	  TTree *tree = (TTree*)f->Get("Events");
	  if (tree == 0) { 
	       fprintf(stderr, "Events tree is not accessible in file %s.  "
		       "Aborting looper.\n", currentFile->GetTitle());
	       exit(1);
	  }

	  cms2.Init(tree);  // set branch addresses for TTree tree

	  TStopwatch t;
	  //Event Loop
	  uint64 nEvents = tree->GetEntries();
	  for( uint64 event = 0; event < nEvents; ++event) {
	       cms2.GetEntry(event);  // get entries for Event number event from branches of TTree tree
	       ++nEventsTotal;
	       
	       // Progress feedback to the user
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
	       
	       // give the analysis a chance to filter out this event
	       // (for example because it's a duplicate)
	       if (FilterEvent())
		    continue;

	       // call the event function
	       FillEventHistos();
	       
	       // call the dilepton candidate function
	       for (unsigned int i_hyp = 0, nHyps = cms2.hyp_type().size();
		    i_hyp < nHyps; ++i_hyp ) {
		    FillDilepHistos(i_hyp);
	       }
	       // trilepton
	       for (unsigned int i_hyp = 0, nHyps = cms2.hyp_trilep_bucket().size();
		    i_hyp < nHyps; ++i_hyp ) {
		    FillTrilepHistos(i_hyp);
	       }
	       // quadlepton
	       for (unsigned int i_hyp = 0, nHyps = cms2.hyp_quadlep_bucket().size();
		    i_hyp < nHyps; ++i_hyp ) {
		    FillQuadlepHistos(i_hyp);
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
     if (logfile != stdout)
	  fclose(logfile);
     if (sample.chain != 0)
	  delete sample.chain;
     return nEventsTotal;
}
