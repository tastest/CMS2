
#ifndef PLOTRESULTS_H
#define PLOTRESULTS_H

#include "TString.h"
#include "TH2F.h"
#include "Tools/HistogramUtilities.h"
#include "Tools/DataSource.h"

//class HistogramUtilities;


void plotResults();

void plotResultsDilep(TString hyp);
void plotResultsLep(TString hyp);

void projectX(TH2F* h, const unsigned int n);
//void doabcd(TH2F* h, double x1=0., double x2=0., double x3=0., double x4=0., double y1=0., double y2=0., double y3=0., double y4=0.);
void oneabcd(TH2F* h, TString title, double x1=0., double x2=0., double x3=0., double x4=0., double y1=0., double y2=0., double y3=0., double y4=0.);
void Nabcd(TH2F** h, int n, double x1=0., double x2=0., double x3=0., double x4=0., double y1=0., double y2=0., double y3=0., double y4=0.);
double abcd(TH2F* h, double x1=0., double x2=0., double x3=0., double x4=0., double y1=0., double y2=0., double y3=0., double y4=0.);
double integrateTH2F(TH2F* h, double xlow=0., double xhgh=0., double ylow=0., double yhgh=0.);
void printCoords(double x1=0., double x2=0., double x3=0., double x4=0., double y1=0., double y2=0., double y3=0., double y4=0.);
void printTH2F(TH2F* h);
void compareTH2F(TH2F* h1, TH2F* h2);
//void aviCD( TH2F* single, TH2F* di, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4);
void aviCDtable( TH2F* data, TH2F* di, TH2F* ssig, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4);
void HShift( TH1D*& h, int shift );
TH2F* getTH2F( HistogramUtilities* h, sources_t theSources, TString var, TString var2, TString hyptyp, Int_t rebin, TString suffix, TString opt);
template <class TH> void printHistStats(TH* h, double xlow=0., double xhgh=-1.);
void makeStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suffix, bool cpylog=false, double ymin=1., double ymax=-999, TString name="");
template <class TH> TH2F* corrTH2F( TH2F* h, TH* met, TH* metacc);
template <class TH> TH1F* getAccCorr( TH* met, TH* metacc, TH* z ); //acceptance correctin--returns th1f even though args are templated

#endif

