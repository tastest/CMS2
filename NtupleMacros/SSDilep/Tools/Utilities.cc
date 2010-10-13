
#include "Utilities.h"

namespace Utilities {

        Double_t getError2(const TH1F *hist, const Int_t binMin, const Int_t binMax)
        {
                Double_t error2 = 0.0;
                for (Int_t i = binMin; i <= binMax; ++i)
                {
                        error2 += pow(hist->GetBinError(i), 2);
                }
                return error2;
        }

        void saveCanvas(const TCanvas *c1, const TString name)
        {
                c1->SaveAs(name + ".eps");
                c1->SaveAs(name + ".png");
                c1->SaveAs(name + ".root");
                c1->SaveAs(name + ".C");
        }

}

