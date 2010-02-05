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

#include "histtools.h"
namespace hist {

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

    while (obj = iter->Next()) {
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

    while (obj = iter->Next()) {
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
		  Float_t xmin, Float_t ymin, Float_t xmax, Float_t ymax) {
    if(! stack) return 0;

    TLegend* leg = new TLegend(xmin, ymin, xmax, ymax);
    leg->SetFillColor(kWhite);
    TList* list = stack->GetHists();
    //       TIterator* iter = list->MakeIterator();

    TObject* obj = 0;

    //Hist color iterator
    Int_t colorIt = 1;

    //while (obj = iter->Next()) {
    for(int i = list->GetSize() - 1 ; i >=0 ; i--) {
      obj = list->At(i);
      if (! obj->InheritsFrom(TH1::Class())) continue;
      TString name(obj->GetName());
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
	  entry = "t#font[12]{W}";
	} else if (entry == "totr"){
	  entry = "Top s/t-channel";
	} else if (entry == "ttdil"){
	  entry = "#font[12]{t#bar{t}}#rightarrow  #font[12]{ll}";
	} else if (entry == "ttotr"){
	  entry = "#font[12]{t#bar{t}} #rightarrow not #font[12]{ll}";
	} else if (entry == "VV"){
	  entry = "#font[12]{WW}, #font[12]{WZ}, #font[12]{ZZ}";
	} else if (entry == "DYtautau"){
	  entry = "DY#rightarrow #tau#tau";
	} else if (entry == "DYeemm"){
	  entry = "DY#rightarrow #font[12]{ee}, #mu#mu";
	} else if (entry == "wjets"){
	  entry = "#font[12]{W}+jets";
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

    while (obj = iter->Next()) {
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
	((TH1*)obj)->Sumw2();
	((TH1*)obj)->Scale(1./integral);
      }
    }
  }

  //Scale by the given value all histograms whose names match the given regular
  //expression pattern or begin with the given prefix.

  void scale(const char* patORpfx, Double_t scale) {
    TRegexp reg(patORpfx, kFALSE);

    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();

    TObject* obj = 0;

    while (obj = iter->Next()) {
      if (! obj->InheritsFrom(TH1::Class())) continue;

      TString name = obj->GetName();

      if (TString(patORpfx).MaybeRegexp()) {
	if (TString(obj->GetName()).Index(reg) < 0 ) continue;
      } else if (! name.BeginsWith(patORpfx)) continue;

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

    while (obj = iter->Next()) {
      if (! obj->InheritsFrom(TH1::Class())) continue;

      if (! top) top = obj;

      if (((TH1*)obj)->GetMaximum() > maxy) maxy = ((TH1*)obj)->GetMaximum();
      if (((TH1*)obj)->GetMinimum() < miny) miny = ((TH1*)obj)->GetMinimum();
    }

    ((TH1*)top)->SetMaximum(maxy*1.3);
    //Protect against log scale
    if (canvas->GetLogy() && ! miny)
      ((TH1*)top)->SetMinimum(1E-4);
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
	     const char* bsmName, bool doRefPats) {
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

    vector<TString*>  pats(0);
    if (orderScheme == 0){
      pats.push_back(new TString(""));
    } else if (orderScheme == 1) {
      //this requires some care: only the histograms beginning with what's in the pattern will be processed
      pats.push_back(new TString("ttdil_"));
      pats.push_back(new TString("ttotr_"));
      pats.push_back(new TString("t_"));
      pats.push_back(new TString("totr_"));
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
	pats.push_back(new TString("ref_totr_"));
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
      pats.push_back(new TString("totr_"));
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
	pats.push_back(new TString("ref_totr_"));
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

    for (unsigned int iPat = 0; iPat < pats.size(); ++iPat){
      TIterator* iter = list->MakeIterator();
      while (obj = iter->Next()) {
	if (! obj->InheritsFrom(TH1::Class())) continue;
	
	TString name = obj->GetName();
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

    while (obj = iter->Next()) {
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

    while (obj = iter->Next()) {
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
  TH1F* eff(TH1F* h1, TH1F* h2, const char* name){

    // first, verify that all histograms have same binning
    // nx is the number of visible bins
    // nxtot = nx+2 includes underflow and overflow
    Int_t nx = h1->GetNbinsX();
    if (h2->GetNbinsX() != nx) {
      cout << "Histograms must have same number of bins" << endl;
      return 0;
    }

    // get the new histogram
    TH1F* temp = (TH1F*) h1->Clone(name);
    temp->SetTitle(name);
    temp->Reset();
    temp->Sumw2();

    // Do the calculation
    temp->Divide(h2,h1,1.,1.,"B");

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
    while(obj=iter->Next()) {
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

  void deleteHistos() {
    // Delete all existing histograms in memory
    TObject* obj;
    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();
    while (obj=iter->Next()) {
      if (obj->IsA()->InheritsFrom(TH1::Class()) ||
	  obj->IsA()->InheritsFrom(TH2::Class()) ) {delete obj;}
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
	cout << "." ;
	cout.flush() ;
      }
    }
    cout << endl ;
    outf.Close() ;

    delete iter ;
  }


  void loadHist(const char* filename, const char* pfx, const char* pat, const char* antipat, Bool_t doAdd)
  {
    TFile inf(filename) ;
    //inf.ReadAll() ;
    TList* list = inf.GetListOfKeys() ;
    TIterator* iter = list->MakeIterator();

    TRegexp re(pat,kTRUE) ;
    if (pat!=0) cout << "pat = " << pat << endl ;
    else cout<<"no pattern: read all"<<endl;

    TRegexp antire(antipat,kTRUE) ;

    gDirectory->cd("Rint:") ;

    TObject* obj ;
    TKey* key ;
    cout << "doAdd = " << (doAdd?"T":"F") << endl ;
    cout << "loadHist: reading." ;
    while(key=(TKey*)iter->Next()) {

      //      Int_t ridx = TString(key->GetName()).Index(re) ;
      //      cout<<"reading "<<key->GetName()<<std::endl;
      if (pat!=0){
	//	cout<<"Check patt "<<pat<<endl;
	if (TString(key->GetName()).Index(re)==-1) {
	  continue ;
	}
      }
      if (antipat!=0){
	if (TString(key->GetName()).Index(antire)!=-1) {
	  continue;
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
      cout << "." ;
      cout.flush() ;
    }
    cout << endl;
    inf.Close() ;
    delete iter ;
  }
}
