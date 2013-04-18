
#include "PlotUtilities.h"

#include "TROOT.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TRegexp.h"
#include "TClass.h"
#include "TKey.h"
#include <TIterator.h>

#include <iostream>

void deleteHistos() {
    // Delete all existing histograms in memory
    TObject* obj;
    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();
    while ( (obj = (iter->Next())) ) {
        if (obj->IsA()->InheritsFrom(TH1::Class()) ||
                obj->IsA()->InheritsFrom(TH2::Class()) ) {
                delete obj;
            }
    }
}

void saveHist(const char* filename, const char* pat)
{

    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();
    TRegexp re(pat,kTRUE) ;

    TFile outf(filename,"RECREATE") ;
    printf("[PlotUtilities::saveHist] Saving histograms to %s\n", filename);

    int counter = 0;
    while(TObject* obj=iter->Next()) {

        // don't save TH1Keys objects
        // this is a bug fudge
        if (TString(obj->GetName()).Contains("histokeys")) continue;

        // save other stuff
        if (TString(obj->GetName()).Index(re)>=0) {
            obj->Write() ;
            ++counter;
        }

    }

    printf("[PlotUtilities::saveHist] Saved %i histograms\n", counter);
    outf.Close() ;
    delete iter ;
}

