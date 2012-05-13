#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TKey.h"
#include "TLegend.h"
#include "TRegexp.h"
#include <TList.h>
#include <TIterator.h>
#include <iostream>

using namespace std;

namespace hist {

     void add(const char* outHistName, const char* patORpfx);

   //Add all histograms whose names match one of ten possible regular expression
   //patterns or begin with one of ten possible given prefixes.  Feel free to
   //mix and match regular expression patterns and prefixes.  If the final hist-
   //ogram named outHistName does not exist it is created.

   void add(const char* outHistName, const char* patORpfx0, const char* patORpfx1, const char* patORpfx2, const char* patORpfx3, const char* patORpfx4, const char* patORpfx5, const char* patORpfx6, const char* patORpfx7, const char* patORpfx8, const char* patORpfx9);

   //Add all histograms whose names match the given regular expression pattern
   //or begin with the given prefix.  If the final histogram named outHistName
   //does not exist it is created.

   void add(const char* outHistName, const char* patORpfx);

   //For all histograms whose names match the given regular expression pattern
   //or begin with the given prefix, set the fill, line and marker colors to the
   //given value.

   void color(const char* patORpfx, Color_t color);

   //Return a pointer to a TLegend with an entry for each histogram drawn on a
   //given TCanvas.  Display either the line, point or fill values.  Optionally
   //apply colors to all histograms.  By default, entry labels are the names of
   //their respective histograms.  Optionally, if histogram names are of the
   //form XX_YY_ZZ_WW, entry labels can be XX (token=0), YY (token=1), etc.

     TLegend* legend(TCanvas* canvas, Option_t* option, Bool_t addColor, Int_t token,
             Float_t xmin, Float_t ymin, Float_t xmax, Float_t ymax);

   //Return a pointer to a TLegend with an entry for each histogram added to a
   //given THStack.  Display either the line, point or fill values.  Optionally
   //apply colors to all histograms.  By default, entry labels are the names of
   //their respective histograms.  Optionally, if histogram names are of the
   //form XX_YY_ZZ_WW, entry labels can be XX (token=0), YY (token=1), etc.

   TLegend* legend(THStack* stack, Option_t* option, Bool_t addColor, Int_t token,
                   Float_t xmin, Float_t ymin, Float_t xmax, Float_t ymax);

   //Normalize to one all histograms whose names match the given regular exp-
   //ression pattern or begin with the given prefix.

   void normalize(const char* patORpfx);

   //Scale by the given value all histograms whose names match the given regular
   //expression pattern or begin with the given prefix.

   void scale(const char* patORpfx, Double_t scale);

   //Don't you hate it when you draw multiple histograms on the same canvas only
   //to find that the bottom histogram's range does not encompass those of the
   //histograms drawn on top?  This method determines the maximum and minimum y
   //range of all the histograms drawn on a given TCanvas and appropriately re-
   //sizes the bottom histogram.

   void setrangey(TCanvas* canvas);

   //Create a stacked histogram consisting of all histograms whose names match
   //the given regular expression pattern or begin with the given prefix.  If
   //the THStack named stackHistName does not exist it is created.  Optionally
   //apply colors to all histograms.  Set drawOption to "nostack" if you do not
   //want to stack, to "hist" to display histograms without errors, to "histe"
   //to display histograms with errors, etc.

   void stack(const char* stackHistName, const char* patORpfx, Bool_t addColor, Option_t* drawOption);

   //Set the x-axis title of all histograms whose names match the given regular
   //expression pattern or begin with the given prefix.

   void xaxis(const char* patORpfx, const char* title);

   //Set the y-axis title of all histograms whose names match the given regular
   //expression pattern or begin with the given prefix.

   void yaxis(const char* patORpfx, const char* title);

}

// Input:  2 histogram
// Output: one histogram which is the efficiency:
// h1 :  TOTAL NUMBER OF EVENTS
// h2 :  NUMBER OF EVENTS THAT PASS

// Method by pointer
TH1F* eff(TH1F* h1, TH1F* h2, const char* name);

// Method by name
TH1F* eff(const char* name1, const char* name2, const char* name);

// Input:  4 histogram
// Output: one histogram which is the BG subtracted efficiency:
// h1 :  TOTAL NUMBER OF EVENTS, SIGNAL REGION
// h2 :  NUMBER OF EVENTS THAT PASS, SIGNAL REGION
// h3 :  TOTAL NUMBER OF EVENTS, SIDE BAND
// h4 :  NUMBER OF EVENTS THAT PASS, SIDE BAND
TH1F* eff_bg(TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4, const char* name);

void deleteHistos();
void histio();
void saveHist(const char* filename, const char* pat);
void loadHist(const char* filename, const char* pfx, const char* pat, Bool_t doAdd);
