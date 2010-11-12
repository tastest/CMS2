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
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TRegexp.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TText.h"
#include "TString.h"
#include "TStyle.h"
#include "TExec.h"
#include <iostream>
#include <algorithm>
#include <set>

/////////////////////////////////////////////////////////
//
//  It's a real soup of god knows what :)
//
////////////////////////////////////////////////////////

//--------------------------------------------------------
// Our histogram are named XX_YY_ZZ
// where XX=tt, WW, WZ, etc
//       YY refers to what is actually plotted
//       ZZ=em, ee, mm, all
//
// It is useful to get a list of all the YY's
//
// We can get this list by looking at all the existing hostograms.
//
// This list is returned here as TObjArray*
//
// The passed prefix should be, eg, "tt" and the passed postfix
// should be, eg, "ee", so that the list is built only
// from the "tt_YY_ee" histograms
//
// Claudio 4 Sep 2007
//
// Added feature to skip 2D histograms
//--------------------------------------------------------

TObjArray* getMyHistosNames( const char* prefix, const char* postfix, bool keep2D) {

  TObjArray* array = new TObjArray(); 
  bool skip2D = !keep2D;

  // Get a list of object and their iterator 
  TList* list = gDirectory->GetList() ;
  TIterator* iter = list->MakeIterator();

  // Loop over objects
  TObject* obj;
  while((obj=iter->Next())) {

    // Only look at objects beginning with 'prefix' and ending with 'postfix'
    TString name = obj->GetName();
    if (name.BeginsWith(prefix) && name.EndsWith(postfix)) {
      
      if (skip2D && obj->IsA()->InheritsFrom(TH2::Class())) continue;

      // Only look at objects that inherit from TH1 or TH2
      if (obj->IsA()->InheritsFrom(TH1::Class()) ||
          obj->IsA()->InheritsFrom(TH2::Class())) {

	// Find the central key, ie, the YY
	TObjArray* t = name.Tokenize("_");
	TString z = TString(t->At(1)->GetName());
	if (t->GetEntries() == 4){
	  z = TString(Form("%s_%s",t->At(1)->GetName(),t->At(2)->GetName()));
	}
	// cout << t->At(1)->GetName() << endl;
	
	// add to the output array
	TObjString* zobj = new TObjString(z);
	array->Add(zobj);

      }
    }
  }

  return array;
}

float GetEntries(const TH1F *h, int lowBin, int highBin) {
  
  if(lowBin < 0 || highBin > h->GetNbinsX()+1) {
    cout << "Bins out of range. Setting the lowBin to the underflow and the highBin to the overflow" << endl;
    lowBin = 0;
    highBin = h->GetNbinsX();
  }
  
  float nentries = 0;
  for(int i = lowBin; i < highBin+1; i++) 
    nentries = nentries + h->GetBinContent(i);
  
  return nentries;
}

float GetTotalError(const TH1F *h, int lowBin, int highBin) {
  
  if(lowBin < 0 || highBin > h->GetNbinsX()+1) {
    cout << "Bins out of range. Setting the lowBin to the underflow and the highBin to the overflow" << endl;
    lowBin = 0;
    highBin = h->GetNbinsX();
  }
  
  float err2 = 0;
  for(int i = lowBin; i < highBin+1; i++) 
    err2 = err2 + pow(h->GetBinError(i),2);
  
  return sqrt(err2);
  
}


std::string formatFloat(double x, const char* formatS) {
  std::string xS = Form(Form("%s", formatS),x);
  double xB = atof(xS.c_str());
  if (x>0 && xB==0){
    xS = Form(" %6.1g",x);
  }
  return xS;
}



double GetMinimum(const vector<TH1F*> &v_hists) {

  TH1* h = NULL; //get the first non-empty histogram 
  for(unsigned int i = 0; v_hists.size(); i++) {
      h = v_hists.at(i);
      if(h->Integral() > 0)
	break;
  }
  if(h->GetMinimum()  < 5e-3)
    return 5e-3;
  else 
    return 
      0.2*(h->GetMinimum());

}
      
  

TLegend* makeLegend(const vector<TH1F*> &v_hists, vector<TString> v_legEntries, bool drawLogY, TString histName) {
  if (v_hists.empty()) return 0;
  //Prefer to draw the Legend on the right half of the 
  //canvas, so only look at the right half
  int nbins = v_hists.at(0)->GetNbinsX();
  TH1F *hdata = NULL;
  TH1F *hmax = NULL;

  for(vector<TH1F*>::const_reverse_iterator rit = v_hists.rbegin(); rit != v_hists.rend(); rit++) {
    if(TString((*rit)->GetName()).Contains("data")) {
      hdata = *rit;
      break;
    }
  }

  for(vector<TH1F*>::const_reverse_iterator rit = v_hists.rbegin(); rit != v_hists.rend(); rit++) {
    if(TString((*rit)->GetName()).Contains("data")) 
      continue;
    hmax = *rit;
    break;
  }
  if ( !hdata || !hmax ) return 0;

  //want the legend to be on the right
  //start with the last but 1 bin, and let this be the highX bin
  //skip the last bin 'cause it has the overflow
  float rangeX = hmax->GetXaxis()->GetXmax() - hmax->GetXaxis()->GetXmin() - hmax->GetBinWidth(nbins);
  int highBinX = nbins - 1;
  float highX = hmax->GetBinLowEdge(nbins) - 0.01*rangeX;
  float lowX = hmax->GetXaxis()->GetXmin() + 0.75*rangeX; 
  int lowBinX = hmax->FindBin(lowX);

  //now we need to figure out what the maxY is in the range
  //and set the Y's of the legend to accomodate
  float max = 1e-6;
  float lowY, highY;
  float rangeY = hmax->GetYaxis()->GetXmax() - hmax->GetYaxis()->GetXmin();
  for(int bin = lowBinX; bin < highBinX+1; bin++) {
    if(hmax->GetBinContent(bin) >= hdata->GetBinContent(bin)) {
      if(max < hmax->GetBinContent(bin))
	max = hmax->GetBinContent(bin);
    } else {
      if(max < hdata->GetBinContent(bin))
	 max = hdata->GetBinContent(bin);
    }
  }
      
  cout << "max: " << max << endl;
  rangeY = hdata->GetMinimum();
  if(hdata->GetMaximum() > hmax->GetMinimum())
    rangeY = hdata->GetMaximum() - rangeY;
  else 
    rangeY = hmax->GetMaximum() - rangeY;
  
  //if((histName.Contains("nJet") || histName.Contains("bTag")) && drawLogY) {

  if(drawLogY) {
    lowY = 1.5*max;
    if(histName == "hdilMass" || histName == "hdilMassNM1" || histName == "hmaxPFJetPtNM1" || histName == "hmetProj")
      highY = lowY + 1000*rangeY;
    else 
      highY = lowY + 30*rangeY;
  } else {
    lowY = 1.2*max;
    highY = lowY + 0.3*rangeY;
  }
  
  cout << "lowY, highY: "<< lowY << "," <<  highY << endl;
  TLegend *leg;
  if(drawLogY)
    leg = new TLegend(lowX,lowY,highX,highY, "", "br"); 
  else 
    leg = new TLegend(0.7, 0.55, 0.92, 0.90, "", "brNDC");
  if((histName.Contains("nJet") || histName.Contains("bTag")) && drawLogY) {
    int temp = hdata->GetNbinsX();
    leg = new TLegend(hdata->GetNbinsX()-1.5, lowY, temp - 0.45*hdata->GetBinWidth(temp), highY, "", "br");
          
  }
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(1001);    

  if(v_hists.size() != v_legEntries.size()) {
    cout << "the number of entries in the legend vector are not the same as the number"
	 << " of entries in the hists vector. Returning a null TLegend. " << endl;
    return NULL;
  }
  vector<TH1F*>::const_reverse_iterator ritH  = v_hists.rbegin();
  vector<TString>::const_reverse_iterator ritE  = v_legEntries.rbegin();
  for(; ritH != v_hists.rend(); ritH++, ritE++) {
    if(*ritE == "Data")
      leg->AddEntry(*ritH, *ritE, "P");
    else
      leg->AddEntry(*ritH, *ritE, "f");
  }




  return leg;
}





TPaveText *getPaveText(const vector<TH1F*> &v_hists, int i_channel, float lumi, bool drawFullErrors) {

  if(v_hists.size() == 0)
    return NULL;
  

  
  //for the MET plot, after the Zveto and the 2jet cut
  TPaveText *pt1 = new TPaveText(0.60, 0.77, 0.80, 0.90, "brNDC");
  if(i_channel == 3)
    pt1 = new TPaveText(0.22, 0.77, 0.45, 0.90, "brNDC");
  

  pt1->SetName("pt1name");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  
  
  TText *blah;
  if(drawFullErrors)
    blah = pt1->AddText("CMS");
  else
    blah = pt1->AddText("CMS Preliminary");
  blah->SetTextSize(0.036);
  blah->SetTextAlign(11);  
  TString temp = formatFloat(100*lumi,"%6.1f");
  temp.ReplaceAll(" " , "" );
  temp = temp + TString(" pb^{-1} at   #sqrt{s}=7 TeV");
  blah = pt1->AddText(temp.Data());
  blah->SetTextSize(0.036);       
  blah->SetTextAlign(11);

  if(i_channel==0)
    blah = pt1->AddText("Events with ee");
  if(i_channel == 1)
    blah = pt1->AddText("Events with #mu#mu");
  if(i_channel == 2)
    blah = pt1->AddText("Events with e#mu");
  if(i_channel == 3)
    blah = pt1->AddText("Events with ee/#mu#mu/e#mu");  
  blah->SetTextSize(0.036);
  blah->SetTextAlign(11);


  return pt1;
}

namespace hist {
  void add(const char* outHistName, const char* patORpfx);

  void add(const char* outHistName,   const char* patORpfx0,     const char* patORpfx1,     const char* patORpfx2 = 0, 
	   const char* patORpfx3 = 0, const char* patORpfx4 = 0, const char* patORpfx5 = 0, const char* patORpfx6 = 0, 
	    const char* patORpfx7 = 0, const char* patORpfx8 = 0, const char* patORpfx9 = 0);
  
  void color(const char* patORpfx, Color_t color);
  
  TLegend* legend(TCanvas* canvas, Option_t* option = "lpf", Bool_t addColor = kFALSE, Int_t token = -1,
		  Float_t xmin = 0.75, Float_t ymin = 0.75, Float_t xmax = 0.99, Float_t ymax = 0.99);
  
  TLegend* legend(THStack* stack, Option_t* option = "lpf", Bool_t addColor = kFALSE, Int_t token = -1,
		  Float_t xmin = 0.75, Float_t ymin = 0.65, Float_t xmax = 0.99, Float_t ymax = 0.99,
		  TH1F* hdata = NULL);
  
   void normalize(const char* patORpfx);
   void scale(const char* patORpfx, Double_t scale);
   void setrangey(TCanvas* canvas);

   void stack(const char* stackHistName, const char* patORpfx, Bool_t addColor = kFALSE, Option_t* drawOption = "", 
	      Int_t orderScheme = 0, const char* bsmName =0, bool doRefPats = false, std::string prefixToSkip = "");//, float yMinimum=1e-3);
   void xaxis(const char* patORpfx, const char* title);

   void yaxis(const char* patORpfx, const char* title);

   TH1F* eff(TH1F* h1, TH1F* h2, const char* name="eff");
   TH1F* eff(const char* name1, const char* name2, const char* name="eff");

   TH1F* eff_bg(TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4, const char* name="eff");

   void deleteHistos(const TString pat = "*");

   void histio();
   void saveHist(const char* filename, const char* pat="*");
   void loadHist(const char* filename, const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE);

   void combineHists(const std::vector<TString> v_prfxsToCombine, const TString outPrefix);
   
   double GetMinimum(const vector<TH1F*> &v_hists);
   TLegend* makeLegend(const vector<TH1F*> &v_hists, vector<TString> v_legEntries, bool drawLogY, TString histName);
   TPaveText *getPaveText(const vector<TH1F*> &v_hists, int i_channel, float lumi, bool drawFullErrors); //call this after makeLegend please

  //Add all histograms whose names match the given regular expression pattern
  //or begin with the given prefix.  If the final histogram named outHistName
  //does not exist it is created.

  void add(const char* outHistName, const char* patORpfx) {
    TRegexp reg(patORpfx, kFALSE);

    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();

    TObject* obj = 0;
    TObject* hist = 0;
    Bool_t makeOutHist = false;

    hist = gDirectory->Get(outHistName);
    //If out hist does not exist, remember to create it
    if (! hist) makeOutHist = true;

    while ((obj = iter->Next())) {
      if (! obj->InheritsFrom(TH1::Class())) continue;

      TString name = obj->GetName();
      //Don't add out hist
      if (name == TString(outHistName)) continue;

      if (TString(patORpfx).MaybeRegexp()) {
	if (TString(obj->GetName()).Index(reg) < 0 ) continue;
      } else if (! name.BeginsWith(patORpfx)) continue;

      if (makeOutHist) {
	hist = obj->Clone(outHistName);

	if (hist->InheritsFrom(TH2::Class()))
	  ((TH2*)hist)->Reset();
	else
	  ((TH1*)hist)->Reset();

	((TH1*)hist)->SetTitle(outHistName);
	if(((TH1*)hist)->GetSumw2()->GetSize() == 0)
	  ((TH1*)hist)->Sumw2();
	makeOutHist = false;
      }

      ((TH1*)hist)->Add((TH1*)obj);
    }
  }

  //Add all histograms whose names match one of ten possible regular expression
  //patterns or begin with one of ten possible given prefixes.  Feel free to
  //mix and match regular expression patterns and prefixes.  If the final hist-
  //ogram named outHistName does not exist it is created.

  void add(const char* outHistName,   const char* patORpfx0,     const char* patORpfx1,     const char* patORpfx2, 
	   const char* patORpfx3, const char* patORpfx4, const char* patORpfx5, const char* patORpfx6, 
	   const char* patORpfx7, const char* patORpfx8, const char* patORpfx9)
  {
    add(outHistName, patORpfx0);
    add(outHistName, patORpfx1);
    if (patORpfx2) add(outHistName, patORpfx2);
    if (patORpfx3) add(outHistName, patORpfx3);
    if (patORpfx4) add(outHistName, patORpfx4);
    if (patORpfx5) add(outHistName, patORpfx5);
    if (patORpfx6) add(outHistName, patORpfx6);
    if (patORpfx7) add(outHistName, patORpfx7);
    if (patORpfx8) add(outHistName, patORpfx8);
    if (patORpfx9) add(outHistName, patORpfx9);
  }


  //For all histograms whose names match the given regular expression pattern
  //or begin with the given prefix, set the fill, line and marker colors to the
  //given value.

  void color(const char* patORpfx, Color_t color) {
    TRegexp reg(patORpfx, kFALSE);

    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();

    TObject* obj = 0;

    while ((obj = iter->Next())) {
      if (! obj->InheritsFrom(TH1::Class())) continue;

      TString name = obj->GetName();

      if (TString(patORpfx).MaybeRegexp()) {
	if (TString(obj->GetName()).Index(reg) < 0 ) continue;
      } else if (! name.BeginsWith(patORpfx)) continue;

      ((TH1*)obj)->SetFillColor(color);
      ((TH1*)obj)->SetLineColor(color);
      ((TH1*)obj)->SetMarkerColor(color);
    }
  }


  void pattern(const char* patORpfx, Color_t pattern) {
    TRegexp reg(patORpfx, kFALSE);
    
    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();
    
    TObject* obj = 0;

    while ((obj = iter->Next())) {
      if (! obj->InheritsFrom(TH1::Class())) continue;

      TString name = obj->GetName();

      if (TString(patORpfx).MaybeRegexp()) {
	if (TString(obj->GetName()).Index(reg) < 0 ) continue;
      } else if (! name.BeginsWith(patORpfx)) continue;

      ((TH1*)obj)->SetFillStyle(pattern);
    }
  }




  //Return a pointer to a TLegend with an entry for each histogram drawn on a
  //given TCanvas.  Display either the line, point or fill values.  Optionally
  //apply colors to all histograms.  By default, entry labels are the names of
  //their respective histograms.  Optionally, if histogram names are of the
  //form XX_YY_ZZ_WW, entry labels can be XX (token=0), YY (token=1), etc.

  TLegend* legend(TCanvas* canvas, Option_t* option, Bool_t addColor, Int_t token,
		  Float_t xmin, Float_t ymin, Float_t xmax, Float_t ymax) {
    if(! canvas) return 0;

    TLegend* leg = new TLegend(xmin, ymin, xmax, ymax);
    TList* list = canvas->GetListOfPrimitives();
    //       TIterator* iter = list->MakeIterator();

    TObject* obj = 0;

    //Hist color iterator
    Int_t colorIt = 1;

    //while (obj = iter->Next()) {
    for(int i = list->GetSize() - 1 ; i >=0 ; i--) {
      obj = list->At(i);
      if (! obj->InheritsFrom(TH1::Class())) continue;

      if (addColor) {
	hist::color(obj->GetName(), colorIt);
	++colorIt;
      }
	 

      if (token == -1)
	leg->AddEntry(obj, obj->GetName(), option);
      else {
	TString name(obj->GetName());
	TObjArray* a = name.Tokenize("_");
	if (a->GetEntries() <= token)
	  leg->AddEntry(obj, obj->GetName(), option);
	else
	  leg->AddEntry(obj, a->At(token)->GetName(), option);
      }
    }

    return leg;
  }

  //Return a pointer to a TLegend with an entry for each histogram added to a
  //given THStack.  Display either the line, point or fill values.  Optionally
  //apply colors to all histograms.  By default, entry labels are the names of
  //their respective histograms.  Optionally, if histogram names are of the
  //form XX_YY_ZZ_WW, entry labels can be XX (token=0), YY (token=1), etc.

  TLegend* legend(THStack* stack, Option_t* option, Bool_t addColor, Int_t token,
		  Float_t xmin, Float_t ymin, Float_t xmax, Float_t ymax,
		  TH1F *hdata) {
    
    if(! stack) return 0;

    TLegend* leg = new TLegend(xmin, ymin, xmax, ymax);
    leg->SetFillColor(kWhite);
    TList* list = stack->GetHists();
    //       TIterator* iter = list->MakeIterator();

    TObject* obj = 0;


    if(hdata != NULL) 
      leg->AddEntry(hdata, "Data", "P");
      

    //Hist color iterator
    Int_t colorIt = 1;

    //while (obj = iter->Next()) {
    for(int i = list->GetSize() - 1 ; i >=0 ; i--) {
      obj = list->At(i);
      if (! obj->InheritsFrom(TH1::Class())) continue;
      TString name(obj->GetName());

      if(name.Contains("QCD"))
	continue;

      if (addColor) {
	hist::color(obj->GetName(), colorIt);
	++colorIt;
      }
	 
      if (token == -1)
	leg->AddEntry(obj, obj->GetName(), option);
      else {
	TString name(obj->GetName());
	TObjArray* a = name.Tokenize("_");
	TString entry(a->GetEntries() <= token ? obj->GetName() : a->At(token)->GetName());
	
	if (entry == "t"){
	  entry = "Single top";
	}else if (entry == "ttbar"){
	  entry = "#font[12]{t#bar{t}}";
	} else if (entry == "ttdil"){
	  entry = "#font[12]{t#bar{t}} signal";
	} else if (entry == "ttotr"){
	  entry = "#font[12]{t#bar{t}} other";
	} else if (entry == "VV"){
	  //entry = "#font[12]{ZZ}, #font[12]{WZ}, #font[12]{ZZ}";
	  entry = "#font[12]{VV}";
	} else if (entry == "DYtautau"){
	  //entry = "DY#rightarrow #tau#tau";
	  entry = "Z/#gamma*#rightarrow#tau^{+}#tau^{-}";
	} else if (entry == "DYeemm"){
	  //entry = "DY#rightarrow #font[12]{ee}, #mu#mu";
	  entry = "Z/#gamma*#rightarrowl^{+}l^{-}";
	} else if (entry == "wjets"){
	  //entry = "#font[12]{W}+jets";
	  entry = "#font[12]{W}#rightarrowl#nu";
	}else {
	  entry.ReplaceAll("tautau", "#tau#tau");
	  entry.ReplaceAll("mm", "#mu#mu");
	}
	leg->AddEntry(obj, entry, option);
      }
    }

    return leg;
  }

  //Normalize to one all histograms whose names match the given regular exp-
  //ression pattern or begin with the given prefix.

  void normalize(const char* patORpfx) {
    TRegexp reg(patORpfx, kFALSE);

    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();

    TObject* obj = 0;

    while ((obj = iter->Next())) {
      if (! obj->InheritsFrom(TH1::Class())) continue;

      TString name = obj->GetName();

      if (TString(patORpfx).MaybeRegexp()) {
	if (TString(obj->GetName()).Index(reg) < 0 ) continue;
      } else if (! name.BeginsWith(patORpfx)) continue;

      Double_t integral = 0;

      if (obj->InheritsFrom(TH2::Class()))
	integral = ((TH2*)obj)->Integral();
      else
	integral = ((TH1*)obj)->Integral();

      if (integral) {
	if( ((TH1*)obj)->GetSumw2()->GetSize() == 0)
	  ((TH1*)obj)->Sumw2();
	((TH1*)obj)->Scale(1./integral);
      }
    }
  }

  //Scale by the given value all histograms whose names match the given regular
  //expression pattern or begin with the given prefix.

  void scale(const char* patORpfx, Double_t scale) {
    TRegexp reg(patORpfx, kTRUE);

    TList* list = gDirectory->GetList()   ;
    TIterator* iter = list->MakeIterator();

    TObject* obj = 0;

    while ((obj = iter->Next())) {
      if (! obj->InheritsFrom(TH1::Class())) continue;

      TString name = obj->GetName();

      if (TString(patORpfx).MaybeRegexp()) {
	if (TString(obj->GetName()).Index(reg) < 0 ) continue;
      } else if (! name.BeginsWith(patORpfx)) continue;

      if( ((TH1*)obj)->GetSumw2()->GetSize() == 0)
	  ((TH1*)obj)->Sumw2();
      ((TH1*)obj)->Scale(scale);
    }
  }

  //Don't you hate it when you draw multiple histograms on the same canvas only
  //to find that the bottom histogram's range does not encompass those of the
  //histograms drawn on top?  This method determines the maximum and minimum y
  //range of all the histograms drawn on a given TCanvas and appropriately re-
  //sizes the bottom histogram.

  void setrangey(TCanvas* canvas) {
    if(! canvas) return ;

    TList* list = canvas->GetListOfPrimitives();
    TIterator* iter = list->MakeIterator();

    TObject* obj = 0;
    TObject* top = 0;

    //Extremes
    Double_t maxy = -999999;
    Double_t miny = 999999;

    while ((obj = iter->Next())) {
      if (! obj->InheritsFrom(TH1::Class())) continue;

      if (! top) top = obj;

      if (((TH1*)obj)->GetMaximum() > maxy) maxy = ((TH1*)obj)->GetMaximum();
      if (((TH1*)obj)->GetMinimum() < miny) miny = ((TH1*)obj)->GetMinimum();
    }

    ((TH1*)top)->SetMaximum(maxy*1.3);
    //Protect against log scale
    if (canvas->GetLogy() && ! miny)
      ((TH1*)top)->SetMinimum(1E-3);
    else
      ((TH1*)top)->SetMinimum(miny*0.7);
  }

  //Create a stacked histogram consisting of all histograms whose names match
  //the given regular expression pattern or begin with the given prefix.  If
  //the THStack named stackHistName does not exist it is created.  Optionally
  //apply colors to all histograms.  Set drawOption to "nostack" if you do not
  //want to stack, to "hist" to display histograms without errors, to "histe"
  //to display histograms with errors, etc.

  void stack(const char* stackHistName, const char* patORpfx, Bool_t addColor, Option_t* drawOption, Int_t orderScheme,
	     const char* bsmName, bool doRefPats, std::string prefixToSkip){//, float yMinimum) {
    TRegexp reg(patORpfx, kTRUE);
      
    TList* list = gDirectory->GetList() ;

    TObject* obj = 0;
    TObject* stack = 0;
    Bool_t makeStackHist = false;

    stack = gDirectory->Get(stackHistName);
    //If stack hist does not exist, remember to create it
    if (! stack) makeStackHist = true;

    //Hist color iterator
    Int_t colorIt = 1;

    std::vector<TString*>  pats(0);
    if (orderScheme == 0){
      pats.push_back(new TString(""));
      //pats.push_back(new TString(""));
      //pats.push_back(new TString("ttdil_"));
      //pats.push_back(new TString("ttotr_"));
      //pats.push_back(new TString("DYeemm_"));
      //pats.push_back(new TString("DYtautau_"));
      //pats.push_back(new TString("wjets_"));
      //pats.push_back(new TString("VV_"));
      //pats.push_back(new TString("tw_"));
      //pats.push_back(new TString("data_"));      
    } else if (orderScheme == 1) {
      //this requires some care: only the histograms beginning with what's in the pattern will be processed
      pats.push_back(new TString("ttdil_"));
      pats.push_back(new TString("ttotr_"));
      pats.push_back(new TString("t_"));
      pats.push_back(new TString("VV_"));
      pats.push_back(new TString("DYtautau_"));
      pats.push_back(new TString("DYeemm_"));
      pats.push_back(new TString("wjets_"));
      pats.push_back(new TString("QCD_"));
      //      pats.push_back(new TString("Vgamma_"));

      if (doRefPats){
	pats.push_back(new TString("ref_ttdil_"));
	pats.push_back(new TString("ref_ttotr_"));
	pats.push_back(new TString("ref_t_"));
	pats.push_back(new TString("ref_VV_"));
	pats.push_back(new TString("ref_DYtautau_"));
	pats.push_back(new TString("ref_DYeemm_"));
	pats.push_back(new TString("ref_wjets_"));
	pats.push_back(new TString("ref_QCD_"));
	//      pats.push_back(new TString("ref_Vgamma_"));
      }
    } else if (orderScheme == 2) {
      //this requires some care: only the histograms beginning with what's in the pattern will be processed
      pats.push_back(new TString("ttdil_"));
      pats.push_back(new TString("ttotr_"));
      pats.push_back(new TString("t_"));
      pats.push_back(new TString("VV_"));
      pats.push_back(new TString("DYtautau_"));
      pats.push_back(new TString("DYeemm_"));
      pats.push_back(new TString("wjets_"));
      pats.push_back(new TString("QCD_"));
      pats.push_back(new TString("Vgamma_"));
      if (bsmName!=0) pats.push_back(new TString(Form("%s_",bsmName)));

      if (doRefPats){
	pats.push_back(new TString("ref_ttdil_"));
	pats.push_back(new TString("ref_ttotr_"));
	pats.push_back(new TString("ref_t_"));
	pats.push_back(new TString("ref_VV_"));
	pats.push_back(new TString("ref_DYtautau_"));
	pats.push_back(new TString("ref_DYeemm_"));
	pats.push_back(new TString("ref_wjets_"));
	pats.push_back(new TString("ref_QCD_"));
	pats.push_back(new TString("ref_Vgamma_"));
	if(bsmName!=0) pats.push_back(new TString(Form("ref_%s_",bsmName)));
      }
    } else {
      pats.push_back(new TString(""));
    }

    for (unsigned int iPat = 0; iPat < pats.size(); ++iPat) {
      TIterator* iter = list->MakeIterator();
      while ((obj = iter->Next())) {
	if (! obj->InheritsFrom(TH1::Class())) continue;
	//((TH1F*)obj)->SetMinimum(yMinimum);
	TString name = obj->GetName();
	if(prefixToSkip != "" && name.BeginsWith(prefixToSkip.c_str())) continue;
	//if(name.Contains("data"))
	if (!( name.BeginsWith(pats[iPat]->Data()))) continue;
	
	if (TString(patORpfx).MaybeRegexp()) {
	  if (TString(obj->GetName()).Index(reg) < 0 ) continue;
	} else if (! name.BeginsWith(patORpfx)) continue;
	
	if (makeStackHist) {
	  stack = new THStack(stackHistName, stackHistName);
	  makeStackHist = false;
	}
	
	if (addColor) {
	  hist::color(obj->GetName(), colorIt);
	  ++colorIt;
	}
	((THStack*)stack)->Add((TH1*)obj, drawOption);
	
      }
      //      delete pats[iPat];//cleanup
    }
  

    // Currently breaks .ls
    //gDirectory->Append(stack);
  }
  
  //Set the x-axis title of all histograms whose names match the given regular
  //expression pattern or begin with the given prefix.
  
  void xaxis(const char* patORpfx, const char* title) {
    TRegexp reg(patORpfx, kFALSE);

    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();

    TObject* obj = 0;

    while ((obj = iter->Next())) {
      if (! (obj->InheritsFrom(TH1::Class()) || obj->InheritsFrom(THStack::Class()))) continue;

      TString name = obj->GetName();

      if (TString(patORpfx).MaybeRegexp()) {
	if (TString(obj->GetName()).Index(reg) < 0 ) continue;
      } else if (! name.BeginsWith(patORpfx)) continue;

      if (obj->InheritsFrom(TH1::Class()))
	((TH1*)obj)->GetXaxis()->SetTitle(title);
      if (obj->InheritsFrom(THStack::Class())) {
	((THStack*)obj)->Draw();
	((THStack*)obj)->GetXaxis()->SetTitle(title);
      }
    }
  }

  //Set the y-axis title of all histograms whose names match the given regular
  //expression pattern or begin with the given prefix.

  void yaxis(const char* patORpfx, const char* title) {
    TRegexp reg(patORpfx, kFALSE);

    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();

    TObject* obj = 0;

    while ((obj = iter->Next())) {
      if (! (obj->InheritsFrom(TH1::Class()) || obj->InheritsFrom(THStack::Class()))) continue;

      TString name = obj->GetName();

      if (TString(patORpfx).MaybeRegexp()) {
	if (TString(obj->GetName()).Index(reg) < 0 ) continue;
      } else if (! name.BeginsWith(patORpfx)) continue;

      if (obj->InheritsFrom(TH1::Class()))
	((TH1*)obj)->GetYaxis()->SetTitle(title);
      if (obj->InheritsFrom(THStack::Class())) {
	((THStack*)obj)->Draw();
	((THStack*)obj)->GetYaxis()->SetTitle(title);
      }
    }
  }

  // Input:  2 histogram
  // Output: one histogram which is the efficiency:
  // h1 :  TOTAL NUMBER OF EVENTS
  // h2 :  NUMBER OF EVENTS THAT PASS

#include "TH1.h"

  // Method by pointer
  TH1F* eff(TH1F* num, TH1F*denom, const char* name){

    // first, verify that all histograms have same binning
    // nx is the number of visible bins
    // nxtot = nx+2 includes underflow and overflow
    Int_t nx = num->GetNbinsX();
    if (denom->GetNbinsX() != nx) {
      cout << "Histograms must have same number of bins" << endl;
      return 0;
    }

    // get the new histogram
    TH1F* temp = (TH1F*) num->Clone(name);
    temp->SetTitle(name);
    temp->Reset();
    if(temp->GetSumw2()->GetSize() == 0 )
      temp->Sumw2();

    // Do the calculation
    temp->Divide(num,denom,1.,1.,"B");

    // Done
    return temp;
  }


  // Method by name
  TH1F* eff(const char* name1, const char* name2, const char* name){

    // Get a list of object and their iterator
    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();

    // Loop over objects, set the pointers
    TObject* obj;
    TH1F* h1=0;
    TH1F* h2=0;
    TString str1 = Form("%s",name1);
    TString str2 = Form("%s",name2);
    while((obj=iter->Next())) {
      TString objName = obj->GetName();
      if (objName == str1) h1 = (TH1F*) obj;
      if (objName == str2) h2 = (TH1F*) obj;
    }

    // quit if not found
    if (h1 == 0) {
      cout << "Histogram " << name1 << " not found" << endl;
      return 0;
    }
    if (h2 == 0) {
      cout << "Histogram " << name2 << " not found" << endl;
      return 0;
    }

    // Call the method by pointer
    TH1F* temp = eff(h1, h2, name);
    return temp;
  }
  // Input:  4 histogram
  // Output: one histogram which is the BG subtracted efficiency:
  // h1 :  TOTAL NUMBER OF EVENTS, SIGNAL REGION
  // h2 :  NUMBER OF EVENTS THAT PASS, SIGNAL REGION
  // h3 :  TOTAL NUMBER OF EVENTS, SIDE BAND
  // h4 :  NUMBER OF EVENTS THAT PASS, SIDE BAND

#include "TH1.h"


  TH1F* eff_bg(TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4, const char* name){

    // first, verify that all histograms have same binning
    // nx is the number of visible bins
    // nxtot = nx+2 includes underflow and overflow
    Int_t nx = h1->GetNbinsX();
    Int_t nxtot = nx + 2;
    if (h2->GetNbinsX() != nx) {
      cout << "Histograms must have same number of bins" << endl;
      return 0;
    }
    if (h3->GetNbinsX() != nx) {
      cout << "Histograms must have same number of bins" << endl;
      return 0;
    }
    if (h3->GetNbinsX() != nx) {
      cout << "Histograms must have same number of bins" << endl;
      return 0;
    }

    // get the new histogram
    TH1F* temp = (TH1F*) h1->Clone(name);
    temp->SetTitle(name);
    temp->Reset();
    if(temp->GetSumw2()->GetSize() == 0)
      temp->Sumw2();

    // Loop over bins, calculate efficiency and error, put it in histogram
    for (Int_t i=0; i<nxtot; i++) {
      Double_t x1 = h1->GetBinContent(i);
      Double_t x2 = h2->GetBinContent(i);
      Double_t x3 = h3->GetBinContent(i);
      Double_t x4 = h4->GetBinContent(i);
      Double_t denom = x1 - x3;
      Double_t eff;
      if (denom == 0.) {
	eff = 0;
      } else {
	eff   = (x2-x4)/denom;
      }
      Double_t failSig = x1 - x2;
      Double_t failBg  = x3 - x4;
      Double_t blah    = (1-eff)*(1-eff)*(x2+x4) + eff*eff*(failSig+failBg);
      if (blah <= 0.) blah=0.0;
      Double_t err;
      if (denom == 0) {
	err = 0.;
      } else {
	err = sqrt(blah)/denom;
      }
      temp->SetBinContent(i,eff);
      temp->SetBinError(i,err);
    }

    // Done
    return temp;
  }

#include <TList.h>
#include <TIterator.h>

  void deleteHistos(const TString pat) {
    // Delete all existing histograms in memory
    TObject* obj;
    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();
    

    while ((obj=iter->Next())) {
      if (obj->IsA()->InheritsFrom(TH1::Class()) ||
	  obj->IsA()->InheritsFrom(TH2::Class()) ) {
	if(pat.Contains("*")) {
	  delete obj;
	  continue;
	}
	if(TString(obj->GetName()).Contains(pat)) {
	  delete obj;
	}
      }
    }
  }
  void histio()
  {
  }

  void saveHist(const char* filename, const char* pat)
  {
    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();

    TRegexp re(pat,kTRUE) ;

    TFile outf(filename,"RECREATE") ;
    while(TObject* obj=iter->Next()) {
      if (TString(obj->GetName()).Index(re)>=0) {
	obj->Write() ;
	//cout << "." ;
	cout.flush() ;
      }
    }
    cout << endl ;
    outf.Close() ;

    delete iter ;
  }


  void loadHist(const char* filename, const char* pfx, const char* pat, Bool_t doAdd)
  {
    TFile inf(filename) ;
    //inf.ReadAll() ;
    TList* list = inf.GetListOfKeys() ;
    TIterator* iter = list->MakeIterator();

    TRegexp re(pat,kTRUE) ;
    if (pat!=0) cout << "pat = " << pat << endl ;
    else cout<<"no pattern: read all"<<endl;

    gDirectory->cd("Rint:") ;

    TObject* obj ;
    TKey* key ;
    cout << "doAdd = " << (doAdd?"T":"F") << endl ;
    cout << "loadHist: reading." ;
    while((key=(TKey*)iter->Next())) {

      //      Int_t ridx = TString(key->GetName()).Index(re) ;
      //      cout<<"reading "<<key->GetName()<<std::endl;
      if (pat!=0){
	//	cout<<"Check patt "<<pat<<endl;
	if (TString(key->GetName()).Index(re)==-1) {
	  continue ;
	}
      }

      obj = inf.Get(key->GetName()) ;
      TObject* clone ;
      if (pfx) {

	// Find existing TH1-derived objects
	TObject* oldObj = 0 ;
	if (doAdd){
	  oldObj = gDirectory->Get(Form("%s_%s",pfx,obj->GetName())) ;
	  if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	    oldObj = 0 ;
	  }
	}
	if (oldObj) {
	  clone = oldObj ;
	  ((TH1*)clone)->Add((TH1*)obj) ;
	} else {
	  clone = obj->Clone(Form("%s_%s",pfx,obj->GetName())) ;
	}


      } else {

	// Find existing TH1-derived objects
	TObject* oldObj = 0 ;
	if (doAdd){
	  oldObj = gDirectory->Get(key->GetName()) ;
	  if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	    oldObj = 0 ;
	  }
	}

	if (oldObj) {
	  clone = oldObj ;
	  ((TH1*)clone)->Add((TH1*)obj) ;
	} else {
	  clone = obj->Clone() ;
	}
      }
      if (!gDirectory->GetList()->FindObject(clone)) {
	gDirectory->Append(clone) ;
      }
      //cout << "." ;
      cout.flush() ;
    }
    cout << endl;
    inf.Close() ;
    delete iter ;
  }


  //takes in a vector of prefixes, and a output prefix. Combines all histograms 
  //whose prefixes matches those in the input vector of TStrings into a new histogram
  //with the the outPrefix
  void combineHists(const std::vector<TString> v_prfxsToCombine, const TString outPrefix) {

    
    if(v_prfxsToCombine.size() == 0) {
      cout << "Input vector size is 0" << endl;
    }
    //histogram of common suffixes
    vector<TString> v_suffixes;
    //get the list of histograms that match the first entry in the prefix
    TString temp = v_prfxsToCombine.at(0);
    TList *rootlist = gDirectory->GetList();
    
    for(int i = 0; i < rootlist->GetSize(); i++) {
      
      TObject *obj = rootlist->At(i);
      TString name = obj->GetName();
      if(!obj->InheritsFrom(TH1::Class()))
	continue;
      if(!name.BeginsWith(temp+"_"))
	continue;

      name.ReplaceAll(temp+"_", "_");

      
      if(obj->InheritsFrom(TH2::Class())) {
	TH2F *h = dynamic_cast<TH2F*>(obj->Clone());
	h->SetName((outPrefix+name).Data());
	h->SetTitle((outPrefix+name).Data());
	for(unsigned int j = 1; j < v_prfxsToCombine.size(); j++) {
	  TH2F *htemp = dynamic_cast<TH2F*>(gDirectory->Get((v_prfxsToCombine.at(j) + name).Data()));
	  h->Add(htemp);
	}	
      } else if(obj->InheritsFrom(TH1::Class())) {
	TH1F *h = dynamic_cast<TH1F*>(obj->Clone());
	h->SetName((outPrefix+name).Data());
	h->SetTitle((outPrefix+name).Data());
	for(unsigned int j = 1; j < v_prfxsToCombine.size(); j++) {
	  TH1F *htemp = dynamic_cast<TH1F*>(gDirectory->Get((v_prfxsToCombine.at(j) + name).Data()));
	  //cout << "TH1F: " << v_prfxsToCombine.at(j) + name << endl;
	  if(htemp == NULL)
	    cout << "htemp is NULL" << endl;
	  h->Add(htemp);
	}
      }
    }//rootlist loop


    //now delete the histos that match the prfxs to combine
    for(unsigned int i = 0; i < v_prfxsToCombine.size(); i++ ) {
      
      // Delete all existing histograms in memory
      TObject* obj;
      TList* list = gDirectory->GetList() ;
      TIterator* iter = list->MakeIterator();
      while ((obj=iter->Next())) {
	TString name = obj->GetName();
	if(name.BeginsWith(outPrefix+"_"))
	  continue;
	if(TString(obj->GetName()).BeginsWith(v_prfxsToCombine.at(i).Data())) {
	  if (obj->IsA()->InheritsFrom(TH1::Class()) ||
	      obj->IsA()->InheritsFrom(TH2::Class()) )
	    delete obj;
	}
      }
    }//loop over prefixes
   
  }//fnc declaration
  
}

// Make the stacks and then browse them interactive
// if (makePictures==true), saves all plots to out/stacks.ps
// Also: make out_stacks_xx.png where xx is the page number 
// on the stack.ps file
//
// keep2D=false to skip the annoying 2D histograms
//


void browseStacks(vector<TString> v_samples, vector<Color_t> v_colors, 
		  TString outfile, vector<TString> v_legEntries, bool drawLogY = false, 
		  vector<Style_t> v_style = vector<Style_t>(), bool drawFullErrors = false, float lumi = 1.0) {

  assert(!v_samples.empty());
  if(v_samples.size() != v_colors.size()) {
    cout << "Number of entries in the vector of samples is not the same as the number of entries in the vector of Color_t" << endl;
    return;
  }
  
  if(v_style.size()!=0 && v_style.size() != v_colors.size()) {
    cout << "Number of entries in the vector of styles is not the same as the number of entries in the vector of samples" << endl;
    return;
  }

  gStyle->SetOptTitle(0);
  
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  // Find out what the names of the existing histograms are
  // The histogram names are XX_YY_ZZ, where XX is the sample,
  // eg, "tt", YY is the actual name, ZZ is the final state, eg, "ee"
  //exclude data if it is found
  bool keep2D = false;
  TObjArray* myNames = getMyHistosNames(v_samples.at(0).Data(),"ee",keep2D);  
  if(myNames->GetSize() == 0) {
    cout << "could not find the sample at the first point in the vector of samples. Exiting." << endl;
    return;
  }
    
  // Now loop over histograms, and make stacks
  TCanvas *c = new TCanvas();  
  c->Divide(2,2);  
  vector<string> v_channel;
  v_channel.push_back("ee");
  v_channel.push_back("mm");
  v_channel.push_back("em");
  v_channel.push_back("all");
  for (int i=0; i<myNames->GetEntries(); i++) {
    for (int i_channel=0; i_channel<4; i_channel++) {
      vector<TH1F*> v_hists;
      TH1F *hdata = NULL;
            
      for(unsigned int i_prefix = 0; i_prefix < v_samples.size(); i_prefix++) {
	
	TString histoName = v_samples.at(i_prefix) + "_" + myNames->At(i)->GetName() + "_" + v_channel.at(i_channel);
	TObject *obj = gDirectory->Get(histoName.Data());
	// assert(obj);
	if(! obj->InheritsFrom(TH1::Class())) 
	  continue;
	
	TH1F *htemp = dynamic_cast<TH1F*>(obj);

	htemp->SetFillColor(v_colors.at(i_prefix));
	htemp->SetLineColor(kBlack);
	TString plot(myNames->At(i)->GetName());
	
	//cout << htemp->GetXaxis()->GetTitle() << endl;
	
	/*
	if(plot.Contains("hmetVal") ||
	   plot.Contains("hmetProjVal") ||
	   plot.Contains("dilMassVal") || 
	   plot.Contains("hmaxPFJetPtVal")) {
	     
	  string xtitle;
	  string ytitle;
	  if(plot.Contains("hmetVal")) {
	    xtitle = "tcMET [GeV]";
	    ytitle = "Events/(5 GeV)";
	  }  else if(plot.Contains("hmetProjVal")) {
	    xtitle = "Projected tcMET [GeV]";
	    ytitle = "Events/(5 GeV)";  
	  } else if(plot.Contains("hmaxPFJetPtVal")) {
	    xtitle = "Leading PF Jet Pt [GeV]";
	    ytitle = "Events/(5 GeV)";
	    //const char *jetbins[4] = {"0", "1", "2", "#geq 3"};
            //for(int k = 0; k<4; k++) 
	    //  htemp->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
	    //htemp->GetXaxis()->SetLabelSize(0.07);
	  } else if(plot.Contains("dilMass")) {
	    xtitle = "Dilepton mass [GeV/c^{2}]";
	    ytitle = "Events/(5 GeV)";
	  }
	  htemp->GetXaxis()->SetTitle(xtitle.c_str());
	  htemp->GetYaxis()->SetTitle(ytitle.c_str());
	*/
	  htemp->GetYaxis()->SetTitleOffset(1.5);
	  htemp->GetXaxis()->SetTitleSize(0.045);
	  htemp->GetYaxis()->SetTitleSize(0.045);
	    
	  //}
	  //cout << htemp->GetXaxis()->GetTitle() << endl;

	if(v_style.size() > 0) {
	  htemp->SetFillStyle(v_style.at(i_prefix));
	  if(v_samples.at(i_prefix) == "data") {
	    htemp->SetMarkerStyle(v_style.at(i_prefix));
	    htemp->SetMarkerColor(v_colors.at(i_prefix));
	  }
	}

	//don't add the data histogram to the stack
	if(v_samples.at(i_prefix) == "data") {
	  hdata = htemp;	  
	  continue;
	}

	if(i_prefix==0) {
	  v_hists.push_back(htemp);
	  continue;
	}
	
	htemp->Add(v_hists.back());
	v_hists.push_back(htemp);
      }//prefix loop			
      if(hdata != NULL)
	v_hists.push_back(hdata);

      c->cd(i_channel+1);

      //now set the Minimum if we want Log Scale
      float min = 0;
      if(drawLogY && i_channel != 2) 
	min = GetMinimum(v_hists);

      
      //set the minimum, before we pass the vector of hists to the legend function
      for(vector<TH1F*>::iterator it = v_hists.begin(); it != v_hists.end(); it++) 
	(*it)->SetMinimum(min);
	

      TLegend *leg = NULL;
      if(i_channel == 3) {
	leg = makeLegend(v_hists, v_legEntries, drawLogY, TString(myNames->At(i)->GetName()));

      } else {
	float max = 0;
	for(vector<TH1F*>::iterator it = v_hists.begin(); it != v_hists.end(); it++) {
	  if((*it)->GetMaximum() > max)
	    max = (*it)->GetMaximum();
	}
	
	

	for(vector<TH1F*>::iterator it = v_hists.begin(); it != v_hists.end(); it++) {
	  if(drawLogY && i_channel != 2 )
	    (*it)->SetMaximum(50*max);
	  else 
	    (*it)->SetMaximum(1.5*max);
	  if(!drawLogY && TString((*it)->GetName()).Contains("hnJet"))
	    (*it)->SetMaximum(6);
	}
      }
	  
	

      //now we gotta draw, in reverse order
      for(vector<TH1F*>::reverse_iterator it = v_hists.rbegin(); it != v_hists.rend(); it++) {
	
	//do not draw if the histogram is empty...1e-6 should be good enough
	if(drawLogY && (*it)->GetMaximum() < 1e-6)
	  continue;

	if(i_channel == 3 && leg) {  
	  float max = 1.1*leg->GetY2();
	  if(drawLogY)
	    max = 5*leg->GetY2();
	  if(v_hists.back()->GetMaximum() < max) 
	    (*it)->SetMaximum(max);
	}	

	if(TString((*it)->GetName()).Contains("data")) {
	  hdata = (*it);
	  if(it == v_hists.rbegin()) 
	    (*it)->Draw("Pe");	 
	  else 	  
	    (*it)->Draw("Pesame");
	} else {
	if(it == v_hists.rbegin()) 
	  (*it)->Draw("hist");	 
	else 	  
	  (*it)->Draw("histsame");
	}
      }
      
      if(hdata != NULL) {
	hdata->Draw("Pesame");
      }

      
      
      //need to draw the first histogram in the stack again, because of the tickmarks

      //now set the Log scale, if desired
      if(drawLogY && i_channel != 2) 
	gPad->SetLogy(1);
      

      //now get the top histogram (sum of all) and draw the errors
      if(drawFullErrors) {
	//if(true) {
	TH1F *h_sumBackGrounds = NULL;
	for(vector<TH1F*>::reverse_iterator it = v_hists.rbegin(); it != v_hists.rend(); it++) {
	  if(TString((*it)->GetName()).Contains("data"))
	    continue;
	  if(h_sumBackGrounds == NULL)
	    h_sumBackGrounds = (TH1F*)((*it)->Clone());
	}
	
	/*
	//temporary for btagging
	if(TString(h_sumBackGrounds->GetName()).Contains("bTag")) {
	  TH1F *h_temp = NULL;
	  for(vector<TH1F*>::reverse_iterator it = v_hists.rbegin(); it != v_hists.rend(); it++) {
	    if(TString((*it)->GetName()).Contains("ttdil"))
	    h_temp= (TH1F*)((*it)->Clone());
	  }
	  cout << __LINE__ << endl;
	  h_temp->SetName("blah");
	  cout << __LINE__ << endl;
	  h_sumBackGrounds->SetBinError(1, 0.49*h_temp->GetBinContent(1));
	  h_sumBackGrounds->SetBinError(2, 0.19*h_temp->GetBinContent(2));	
	  h_sumBackGrounds->SetBinError(3, 0.46*h_temp->GetBinContent(3));
	  //end temporary for btagging
	}
	*/
	h_sumBackGrounds->SetDirectory(rootdir);
	string name = "h_allButData_" + v_channel.at(i_channel);
	TExec *setex = new TExec("setex","gStyle->SetErrorX(0.5)");
	setex->Draw();
	h_sumBackGrounds->SetName(name.c_str());
	h_sumBackGrounds->SetFillColor(kBlack);
	h_sumBackGrounds->SetLineColor(kBlack);
	//h_sumBackGrounds->SetMarkerColor(kBlack);
	h_sumBackGrounds->SetFillStyle(3004);
	h_sumBackGrounds->Draw("samee2");
	
	//add histogram to the legend
	if(i_channel == 3 && leg) {
	  leg->SetX1(0.60);
	  leg->SetX2(0.96);
	  leg->SetY1(0.57);
	  leg->SetY2(0.92);
	  leg->AddEntry(h_sumBackGrounds, "Bckg. uncertainty" ,"f");
	  
	}	
	//go back to default
	TExec *setexDef = new TExec("setexDef","gStyle->SetErrorX(0.0)");
	setexDef->Draw();
	
	float totalErr2 = 0;
	for(int ibinx = 3 ; ibinx < h_sumBackGrounds->GetNbinsX() + 1; ibinx++) 
	  totalErr2 = totalErr2 + pow(h_sumBackGrounds->GetBinError(ibinx),2);
	//cout << "histname, bin,error " << h_sumBackGrounds->GetName() << " " << ibinx << " " << h_sumBackGrounds->GetBinError(ibinx) << endl;
	cout << "Total error in high jet bins in " << h_sumBackGrounds->GetName() << " is " << sqrt(totalErr2) << endl;

	
	gPad->RedrawAxis();
	
      }
      
      
      if(i_channel == 3 && leg) 
	leg->Draw();

      TPaveText *pt = getPaveText(v_hists, i_channel, lumi, drawFullErrors);
      pt->Draw();
      
      gPad->RedrawAxis();	
      
      //if(i_channel == 3)
      //gPad->SaveAs("nJets_final_allwithErrors.eps");

    }//channel loop
    c->Modified();
    c->Update();

    
    //cout << myNames->At(i)->GetName() << endl;
    TString tempstring = outfile.ReplaceAll(".root", "");
    if(drawFullErrors)
      tempstring = tempstring + "_withErrors";
    //c->SaveAs((tempstring + ".C").Data());

    if(i != myNames->GetEntries() - 1) 
      c->Print((tempstring + ".eps(").Data());
    else
    c->Print((tempstring + ".eps)").Data());

    


  }//histos loop



}

