
#include <assert.h>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TTreeCache.h"
//#include "../Tools/goodrun.cc"

#include "../CORE/CMS2.cc"
//CMS2 cms2;

#include "Rtypes.h"
typedef ULong64_t uint64;

using namespace tas;

// Used to filter an ntuple based on the 'select' function return
// infile is the path to the ntuple you want to filter (* allowed)
// outfile is the result
// cut is specified in the select() function

// eg:
// ntupleFilter("/data/tmp/cms2-V03-00-10/MinimumBias_BeamCommissioning09-rereco_FIRSTCOLL_v1/*.root","/data/tmp/wandrews/minbias/filtered_ntuple.root")
//infile and outfile must end in ".root"

// WARNING: all of the input files are put into 1 output file, so
// please be careful you don't create a file which is enormous (unless
// you can handle enormous files)

// ECAL Relative Isolation, Non-Truncated
float electronIsolation_ECAL_rel_v1(const unsigned int index, bool useEBps = true){
  float pt               = cms2.els_p4().at(index).pt();                                                                  // Electron Pt
  float ecal_sum_over_pt = 0.0;                                                                                           // ECAL Relative Isolation, NT
  if( fabs(cms2.els_etaSC().at(index)) > 1.479  ) ecal_sum_over_pt += cms2.els_ecalIso().at(index);                       // EE: Ecal Endcap  
  if( fabs(cms2.els_etaSC().at(index)) <= 1.479 ) {
      if (useEBps)
          ecal_sum_over_pt += max( 0.0, ( cms2.els_ecalIso().at(index) - 1.0 ) ); // EB: Ecal Barrel
      else
          ecal_sum_over_pt += cms2.els_ecalIso().at(index); // EB: Ecal Barrel
  }
  ecal_sum_over_pt /= pt;
  return ecal_sum_over_pt;
}

// HCAL Relative Isolation, Non-Truncated
float electronIsolation_HCAL_rel_v1(const unsigned int index){
  float pt               = cms2.els_p4().at(index).pt();      // Electron Pt
  float hcal_sum_over_pt = cms2.els_hcalIso().at(index) / pt; // HCAL Relative Isolation, NT
  return hcal_sum_over_pt;
}

// Relative Isolation, Non-Truncated
float electronIsolation_rel_v1(const unsigned int index, bool use_calo_iso){
    float pt               = cms2.els_p4().at(index).pt();          // Electron Pt
    float TRCK_sum_over_pt = cms2.els_tkIso().at(index) / pt;       // Tracker Relative Isolation, Non-Truncated
    float ECAL_sum_over_pt = electronIsolation_ECAL_rel_v1(index);  // ECAL    Relative Isolation, Non-Truncated
    float HCAL_sum_over_pt = electronIsolation_HCAL_rel_v1(index);  // HCAL    Relative Isolation, Non-Truncated

    float sum_over_pt      = TRCK_sum_over_pt;                      // Combined Subdetector Relative Isolation, Non-Truncated
    if(use_calo_iso){
      sum_over_pt += ECAL_sum_over_pt;
      sum_over_pt += HCAL_sum_over_pt;
    }
    return sum_over_pt;
}


bool select ()
{
	 //hyp filter
	 //if( hyp_p4().size() > 0 )
	 //return true;
/*
  if( !goodrun( evt_run(), evt_lumiBlock() ) )
  return false;
*/
//	 const float ptthresh = 10.;

/*
  for( unsigned int i=0; i<pfjets_p4().size(); i++ )
  if( pfjets_p4()[i].pt() > ptthresh )
  return true;
  for( unsigned int i=0; i<jets_p4().size(); i++ )
  if( jets_p4()[i].pt()*jets_cor()[i] > ptthresh )
  return true;

*/
/*
  for (size_t i = 0; i < cms2.mus_ndof().size(); ++i) 
  if (cms2.mus_p4()[i].Pt() > ptthresh) return true;
  for (size_t i = 0; i < cms2.evt_nels(); ++i) 
  if (cms2.els_p4()[i].Pt() > ptthresh) return true;
*/

	 for (unsigned int hypi = 0; hypi < cms2.hyp_p4().size(); ++hypi)
	 {

         int id_lt = hyp_lt_id()[hypi];
         int id_ll = hyp_ll_id()[hypi];
         int type = hyp_type()[hypi];


         if (type == 0) { // mm
           bool mumubool = true; 
           unsigned int iMut = hyp_lt_index()[hypi];
           unsigned int iMul = hyp_ll_index()[hypi];

           if (cms2.hyp_lt_p4()[hypi].pt() < 5.) mumubool = false; 
           if (cms2.hyp_ll_p4()[hypi].pt() < 5.) mumubool = false; 

           if (((cms2.mus_type().at(iMut)) & (1<<1)) == 0) mumubool = false;          
           if (((cms2.mus_type().at(iMul)) & (1<<1)) == 0) mumubool = false;          

           if (((cms2.mus_type().at(iMut)) & (1<<2)) == 0) mumubool = false;
           if (((cms2.mus_type().at(iMul)) & (1<<2)) == 0) mumubool = false;
           if (mumubool) return true;
         }

         if (type == 3) { // ee
           bool eebool = true; 
           unsigned int iElt = hyp_lt_index()[hypi];
           unsigned int iEll = hyp_ll_index()[hypi];

           if (cms2.hyp_lt_p4()[hypi].pt() < 10.) eebool = false;
           if (cms2.hyp_ll_p4()[hypi].pt() < 10.) eebool = false;
           if( electronIsolation_rel_v1(iElt, true ) > 1.0) eebool = false;
           if( electronIsolation_rel_v1(iEll, true ) > 1.0) eebool = false;

           if (eebool) return true;
         }

        
         if (type == 1 || type == 2) { // em
           bool embool = true;
           int iEl = 0;
           int iMu = 0;

           if (abs(id_lt) == 13) iMu = hyp_lt_index()[hypi];
           if (abs(id_ll) == 13) iMu = hyp_ll_index()[hypi];

           if (abs(id_lt) == 11) iEl = hyp_lt_index()[hypi];
           if (abs(id_ll) == 11) iEl = hyp_ll_index()[hypi];

           if (((cms2.mus_type().at(iMu)) & (1<<1)) == 0) embool = false;
           if (((cms2.mus_type().at(iMu)) & (1<<2)) == 0) embool = false;
           if (mus_p4().at(iMu).pt() < 5.) embool = false;
 
           if( electronIsolation_rel_v1(iEl, true ) > 1.0) embool = false;
           if (els_p4().at(iEl).pt() < 10.) embool = false;
           if (embool) return true;
         }
  
//		  if (cms2.hyp_lt_p4()[hypi].pt() > 20. && cms2.hyp_ll_p4()[hypi].pt() > 10.) return true;
//		  if (cms2.hyp_lt_p4()[hypi].pt() > 10. && cms2.hyp_ll_p4()[hypi].pt() > 20.) return true;
	 }

	 return false;

}

void ntupleFilterSSign (const std::string &infile, const std::string &outfile, bool printPass=false)  
{
     // output file and tree
     TFile *output =TFile::Open(outfile.c_str(), "RECREATE");
     assert(output != 0);
     TTree *newtree = 0;

     const long long max_tree_size = 20000000000000000LL;
     TTree::SetMaxTreeSize(max_tree_size);

	 FILE *log = 0; //for keeping any output desired on selection
	 if( printPass ) {
		  size_t pos = outfile.find(".root");
		  assert( pos != string::npos );
		  std::string outcpy = outfile;
		  log = fopen( outcpy.replace(pos, 5, "_run_lumi_event").c_str(), "w" );
	 }
	 
     TChain *chain = new TChain("Events");
     chain->Add(infile.c_str());
     TObjArray *listOfFiles = chain->GetListOfFiles();
     const uint64 nEventsChain = chain->GetEntries();
     uint64 nEventsTotal = 0;
     uint64 nEventsSelected = 0;

     // file loop
     TIter fileIter(listOfFiles);
     TFile *currentFile = 0;
     bool first = true;
     int i_permille_old = 0;
     while ( (currentFile = (TFile*)fileIter.Next()) ) {
		  TFile f(currentFile->GetTitle());
		  //const char *name = f.GetName();
		  TTree *tree = (TTree*)f.Get("Events");
//                  TTreeCache::SetLearnEntries(10);
//                  tree->SetCacheSize(128*1024*1024);

		  //chain->SetBranchStatus(chain->GetAlias("evt_scale1fb"		), 0);
		  //chain->SetBranchStatus(chain->GetAlias("evt_xsec_excl"	), 0);
		  //chain->SetBranchStatus(chain->GetAlias("evt_xsec_incl"	), 0);
		  //chain->SetBranchStatus(chain->GetAlias("evt_kfactor"		), 0);
		  //chain->SetBranchStatus(chain->GetAlias("evt_nEvts"		), 0);
		  //chain->SetBranchStatus(chain->GetAlias("evt_filt_eff"		), 0);

		  chain->SetBranchStatus("EventAuxiliary",0);
		  // for the first file, clone the tree
		  if ( first ) {
			   newtree = chain->CloneTree(0);
			   newtree->SetDirectory(output);
			   first = false;
	       
		  }
	  
		  // init
		  cms2.Init(newtree);
		  cms2.Init(tree);

		  // Event Loop
		  const unsigned int nEvents = tree->GetEntries();
		  for (unsigned int event = 0; event < nEvents; ++event, ++nEventsTotal) {

//                           tree->LoadTree(event); 

			   int i_permille = (int)floor(10000 * nEventsTotal / float(nEventsChain));
			   if (i_permille != i_permille_old) {
					// xterm magic from L. Vacavant and A. Cerri
					if (isatty(1)) {
						 printf("\015\033[32m ---> \033[1m\033[31m%5.2f%%"
								"\033[0m\033[32m <---\033[0m\015", i_permille/100.);
						 fflush(stdout);
					}
					i_permille_old = i_permille;
			   }

			   cms2.GetEntry(event);
			   //set condition to skip event
			   if (not select()) 
					continue;

			   ++nEventsSelected;
			   if( printPass ) {
					//cout << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl;
					fprintf(log, "%i %i %i\n", cms2.evt_run(), cms2.evt_lumiBlock(), cms2.evt_event() );
			   }

			   cms2.LoadAllBranches();
	       
			   // fill the new tree
			   newtree->Fill();
		  }
     }

	 if( printPass ) {
		  fprintf(log, "\nTotal events run on: %llu\n", nEventsTotal);
		  fprintf(log, "Num events selected: %llu\n", nEventsSelected ); //need two fprintf statements bc of some gcc bug
		  //cout << endl
		  //	  << "Total events run on: " << nEventsTotal << endl
		  //	  << "Num events selected: " << nEventsSelected << endl;
		  //<< "Copy finished. Closing Files" << endl;
	 }
	 
     output = newtree->GetCurrentFile();
	 output->cd();
     newtree->Write();
     delete output;
}
