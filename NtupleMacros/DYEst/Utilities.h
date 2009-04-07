
#ifndef DLEUTILITIES_H
#define DLEUTILITIES_H

#include "TCanvas.h"
#include "TH1F.h"

#include <math.h>

namespace Utilities {
        Double_t getError2(const TH1F *hist, const Int_t binMin, const Int_t binMax);
        void saveCanvas(const TCanvas *c1, const TString name);
}

#endif

