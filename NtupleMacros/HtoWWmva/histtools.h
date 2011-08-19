#include "TH1F.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCanvas.h"

namespace hist {
    void add(const char* outHistName, const char* patORpfx0, const char* patORpfx1, const char* patORpfx2 = 0, const char* patORpfx3 = 0, const char* patORpfx4 = 0, const char* patORpfx5 = 0, const char* patORpfx6 = 0, const char* patORpfx7 = 0, const char* patORpfx8 = 0, const char* patORpfx9 = 0);
    void add(const char* outHistName, const char* patORpfx);
    void color(const char* patORpfx, Color_t color);
    void invertx(const char* patORpfx);
    void colors(TCanvas* canvas, Color_t color = 1);
    void colors(THStack* stack, Color_t color = 1);
    void drawsame(const char* canvasName, const char* patORpfx, Bool_t addColor = kFALSE, Bool_t addStyle = kFALSE, Option_t* drawOption = "");
    TLegend* legend(const char* canvasName, Option_t* option = "lpf", Bool_t addColor = kFALSE, Int_t token = -1, Float_t xmin = 0.75, Float_t ymin = 0.75, Float_t xmax = 0.99, Float_t ymax = 0.99);
    TLegend* legend(TCanvas* canvas, Option_t* option = "lpf", Bool_t addColor = kFALSE, Int_t token = -1, Float_t xmin = 0.75, Float_t ymin = 0.75, Float_t xmax = 0.99, Float_t ymax = 0.99);
    TLegend* legend(THStack* stack, Option_t* option = "lpf", Bool_t addColor = kFALSE, Int_t token = -1, Float_t xmin = 0.75, Float_t ymin = 0.75, Float_t xmax = 0.99, Float_t ymax = 0.99);
    void normalize(const char* patORpfx);
    void scale(const char* patORpfx, Double_t scale);
    void scalebins(const char* patORpfx, Double_t scale);
    void setrangey(const char* canvasName, Double_t maxy = -999999, Double_t miny = 999999);
    void setrangey(TCanvas* canvas, Double_t maxy = -999999, Double_t miny = 999999);
    void stack(const char* stackHistName, const char* patORpfx, Bool_t addColor = kFALSE, Option_t* drawOption = "");
    void style(const char* patORpfx, Style_t style = 20);
    void styles(TCanvas* canvas, Style_t style = 20);
    void styles(THStack* stack, Style_t style = 20);
    void xaxis(const char* patORpfx, const char* title);
    void yaxis(const char* patORpfx, const char* title);
}
TH1F* eff(TH1F* h1, TH1F* h2, const char* name="eff");
TH1F* eff(const char* name1, const char* name2, const char* name="eff");
TH1F* eff_bg(TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4, const char* name="eff");
void deleteHistos();
void saveHist(const char* filename, const char* pat="*");
void loadHist(const char* filename, const char* directory = 0, const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE);
