#ifndef MACROS_histtools_H
#define MACROS_histtools_H
#include "TList.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TRegexp.h"
#include "TKey.h"
#include <iostream>

namespace hist {
  
  void add(const char* outHistName, const char* patORpfx);

   void add(const char* outHistName,   const char* patORpfx0,     const char* patORpfx1,     const char* patORpfx2 = 0, 
	    const char* patORpfx3 = 0, const char* patORpfx4 = 0, const char* patORpfx5 = 0, const char* patORpfx6 = 0, 
	    const char* patORpfx7 = 0, const char* patORpfx8 = 0, const char* patORpfx9 = 0);

   void color(const char* patORpfx, Color_t color);

   TLegend* legend(TCanvas* canvas, Option_t* option = "lpf", Bool_t addColor = kFALSE, Int_t token = -1,
                   Float_t xmin = 0.75, Float_t ymin = 0.75, Float_t xmax = 0.99, Float_t ymax = 0.99);

   TLegend* legend(THStack* stack, Option_t* option = "lpf", Bool_t addColor = kFALSE, Int_t token = -1,
                   Float_t xmin = 0.75, Float_t ymin = 0.65, Float_t xmax = 0.99, Float_t ymax = 0.99);

   void normalize(const char* patORpfx);
   void scale(const char* patORpfx, Double_t scale);
   void setrangey(TCanvas* canvas);

   void stack(const char* stackHistName, const char* patORpfx, Bool_t addColor = kFALSE, Option_t* drawOption = "");
   void xaxis(const char* patORpfx, const char* title);

   void yaxis(const char* patORpfx, const char* title);

   TH1D* eff(TH1D* h1, TH1D* h2, const char* name="eff");
   TH1D* eff(const char* name1, const char* name2, const char* name="eff");

   TH1D* eff_bg(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, const char* name="eff");

   void deleteHistos();

   void histio();
   void saveHist(const char* filename, const char* pat="*");
   void loadHist(const char* filename, const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE);
   void combineHists(const std::vector<TString> v_prfxsToCombine, const TString outPrefix);

}


#endif //MACROS_histtools_H
