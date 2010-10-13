#include <assert.h>
#include <map>
#include <utility>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TPRegexp.h"
#include "TTree.h"
#include <cstdio>

#include "CORE/CMS2.h"

#include "Rtypes.h"
typedef ULong64_t uint64;

bool select ()
{
    bool pass = false;
    for(unsigned int hypi = 0; hypi < cms2.hyp_p4().size(); ++hypi)
    {
        int index1 = cms2.hyp_lt_index()[hypi];
        int index2 = cms2.hyp_ll_index()[hypi];

        if(min(cms2.hyp_lt_p4()[hypi].pt(), cms2.hyp_ll_p4()[hypi].pt()) < 5.)
            continue;
        if(abs(cms2.hyp_lt_id()[hypi]) == 13 && !(cms2.mus_type()[index1] & 6))
            continue;
        if(abs(cms2.hyp_ll_id()[hypi]) == 13 && !(cms2.mus_type()[index2] & 6))
            continue;

        pass = true;
        break;
    }

    return pass;
}

void dilepskim (const std::string &infile, const std::string &outfile, bool printPass=false)  
{
     // output file and tree
     TFile *output =TFile::Open(outfile.c_str(), "RECREATE");
     assert(output != 0);
     TTree *newtree = 0;

     const long long max_tree_size = 20000000000000000LL;
     TTree::SetMaxTreeSize(max_tree_size);

     TChain *chain = new TChain("Events");
     chain->Add(infile.c_str());
     TObjArray *listOfFiles = chain->GetListOfFiles();
     const uint64 nEventsChain = chain->GetEntries();
     uint64 nEventsTotal = 0;
     uint64 nEventsSelected = 0;

     TPRegexp preg("\\S+/merged_ntuple_(\\d+_\\d+).root");
     // file loop
     TIter fileIter(listOfFiles);
     TFile *currentFile = 0;
     bool first = true;
     int i_permille_old = 0;
     while ( currentFile = (TFile*)fileIter.Next() ) {
         TFile f(currentFile->GetTitle());
         //const char *name = f.GetName();
         TTree *tree = (TTree*)f.Get("Events");

         FILE *log = 0; //for keeping any output desired on selection
         if( printPass ) {
             TString ident = ((TObjString*)preg.MatchS(TString(currentFile->GetTitle()))->At(1))->GetString();
             TString logFileName = "";
             logFileName.Append("logs/log_");
             logFileName.Append(ident);
             logFileName.Append(".txt");
             log = fopen( logFileName.Data(), "w" );
         }

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

         // Keep track of and record passed triggers
         // on file-by-file basis
         std::map<TString,int> trigCountMap;
         for(unsigned int trigi = 0; trigi < cms2.hlt_trigNames().size(); ++trigi)
             trigCountMap.insert(std::make_pair(cms2.hlt_trigNames()[trigi], 0));

         // Event Loop
         const unsigned int nEvents = tree->GetEntries();
         for (unsigned int event = 0; event < nEvents; ++event, ++nEventsTotal) {
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

             // triggers++
             std::map<TString,int>::iterator mapit;
             for(unsigned int trigi = 0; trigi < cms2.hlt_trigNames().size(); ++trigi)
                 if (cms2.passHLTTrigger(cms2.hlt_trigNames()[trigi]))
                     (trigCountMap[cms2.hlt_trigNames()[trigi]])++;

             //set condition to skip event
             if (not select()) 
                 continue;

             ++nEventsSelected;
             //if( printPass ) {
             //    fprintf(log, "%i %i %i\n", cms2.evt_run(), cms2.evt_lumiBlock(), cms2.evt_event() );
             //}

             cms2.LoadAllBranches();

             // fill the new tree
             newtree->Fill();
         }

         if( printPass ) {
             fprintf(log, "TotalEvents %i\n", nEventsTotal);
             fprintf(log, "SelectedEvents %i\n", nEventsSelected); //need two fprintf statements bc of some gcc bug

             // Dump trigCountMap
             std::map<TString,int>::const_iterator trigIter;
             for(trigIter = trigCountMap.begin(); trigIter != trigCountMap.end(); ++trigIter)
                 fprintf(log, "%s %i\n", trigIter->first.Data(), trigIter->second);
             //cout << endl
             //	  << "Total events run on: " << nEventsTotal << endl
             //	  << "Num events selected: " << nEventsSelected << endl;
             //<< "Copy finished. Closing Files" << endl;

             fclose(log);
         }
     }

     output->cd();
     newtree->Write();
     delete output;
}
