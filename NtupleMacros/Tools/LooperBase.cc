#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "LooperBase.h"
#include "Sample.h"
#include "tools.h"
#include "../CORE/CMS2.h"

LooperBase::LooperBase (Sample s, cuts_t c, const char *fname) : 
     sample_(s), cuts_(c), hasRun_(false)
{
     if (fname != 0 && strlen(fname) != 0) {
	  logfile_ = fopen(fname, "a");
	  if (logfile_ == 0)
	       perror("opening log file");
     } else {
	  logfile_ = stdout;
     }
}

LooperBase::~LooperBase ()
{

}

double LooperBase::Weight (int i_hyp)
{
     return cms2.evt_scale1fb() * sample_.kFactor;
}

uint64 LooperBase::Loop ()
{
     printf("Processing %s\n", sample_.name.c_str());
     // change to histogram directory
     TDirectory *old_gDirectory = gDirectory;
     gDirectory = histo_directory;
     // book histos
     BookHistos();
     // change back to current directory
     gDirectory = old_gDirectory;
     Begin();
     TChain *chain = sample_.chain;
     uint64 nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
     uint64 nEventsTotal = 0;
     double weightEventsTotal = 0;
     memset(hypos_total_n_, 0, sizeof(hypos_total_n_));
     memset(hypos_total_weight_, 0, sizeof(hypos_total_weight_));

     // clear list of duplicates
     already_seen.clear();
     duplicates_total_n_ = 0;
     duplicates_total_weight_ = 0.;

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
// 	       fprintf(logfile_, "%f\n", LooperBase::Weight(0));
	       weightEventsTotal += LooperBase::Weight(0);
	       
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
 	       if (not filterByProcess(sample_.process)) 
		    continue;
	       
	       // give the analysis a chance to filter out this event
	       // (for example because it's a duplicate)
	       if (FilterEvent()) {
		    duplicates_total_n_++;
		    duplicates_total_weight_ += LooperBase::Weight(0);
		    continue;
	       }

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
	  exit(1);
     }
     
     int ret = fprintf(logfile_, 
		       "Sample %10s: Events before cut (weight): %10.1f.  Events: %8u\n",
		       sample_.name.c_str(),
		       weightEventsTotal, nEventsTotal);
     if (ret < 0)
	  perror("writing to log file");
     ret = fprintf(logfile_, 
		   "Sample %10s: Duplicate events: %8u.  Weight %10.1lf\n",
		   sample_.name.c_str(),
		   duplicates_total_n_, duplicates_total_weight_);
     if (ret < 0)
	  perror("writing to log file");
     End();
     if (logfile_ != stdout)
	  fclose(logfile_);
     if (sample_.chain != 0)
	  delete sample_.chain;
     hasRun_ = true;
     return nEventsTotal;
}
