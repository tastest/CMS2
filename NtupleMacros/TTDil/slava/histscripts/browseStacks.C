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
#include "TStyle.h"
#include <iostream>

#include "getMyHistosNames.h"
#include "histtools.h"

// Make the stacks and then browse them interactively
// if (makePictures==true), saves all plots to out/stacks.ps
// Also: make out_stacks_xx.png where xx is the page number 
// on the stack.ps file
//
// keep2D=false to skip the annoying 2D histograms
//

Int_t FavoriteHatches[10] = {3003,3004,3005,3006,3002,3244,3290,3315,3351,3025};

void hatch(const char* patORpfx,  Int_t hatch) {
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
	
    ((TH1*)obj)->SetFillStyle(hatch);
	
  }
}




void browseStacks( bool makePictures=false, bool wait=true , int addHistName = 1, Double_t maxYScaleF = 1., 
		   bool logScale = false, bool setMinZero = true, int colorScheme = 0, bool noLegend = false, int orderScheme = 0,
		   char* bsmName = "") {
  //addHistName =
  // 0 -- no name
  // 1 -- native name
  // 2 -- CMS Preliminary
  // 3 -- CMS Simulation
  // 4 -- CMS Preliminary, L = 10 pb-1
  // 5 -- CMS Simulation, L =  10 pb-1
  // 6 -- CMS Preliminary, L = 100 pb-1
  // 7 -- CMS Simulation, L = 100 pb-1

  // colorScheme = 
  // 0 -- do nothing, use defaults
  // 1 -- french style: 
  /*
    TTbar signal: 3
    TTbar background: 24
    Z+jets: 2
    W+jets: 9
    dibosons: 33
    Vqq: 42
    single top (tW, t-channel, s-channel): 46 
  */
  if (colorScheme ==1){
    hist::color("^ttdil_", 3);
    hist::color("^ttotr_", 24);
    hist::color("^DYeemm_", 2);
    hist::color("^DYtautau_", 1);
    hist::color("^VV_", 33);
    hist::color("^wjets_", 9);
    hist::color("^QCD_", 51); //use the same as in doAll
    hist::color("^t_", 46);      
    hist::color("^totr_", 49);      
  }
  if (colorScheme ==2){ //based on the style in the PAS
    hist::color("^ttdil_", kGreen);
    hist::color("^ttotr_", kYellow);
    hist::color("^DYeemm_", kRed);
    hist::color("^DYtautau_", kBlack);
    hist::color("^VV_", kGray);
    hist::color("^wjets_", kViolet);
    hist::color("^QCD_", kViolet+3); //use the same as in doAll
    hist::color("^t_", kRed-3);      
    hist::color("^totr_", kRed+3);      
    hist::color("^Vgamma_", kViolet-4);      
    hist::color(Form("^%s_",bsmName), kBlue);
  }

  gStyle->SetOptTitle(0);

  bool keep2D=false;

  //fix the hNJet histos
  TList *list = gDirectory->GetList();
  TIterator *iter = list->MakeIterator();
  TObject *obj = 0;
  while(obj = iter->Next()) {
    if (obj->InheritsFrom(TH1::Class())){
      TH1F* hobj = (TH1F*)obj;
      hobj->SetLineColor(1);
      hobj->SetLineWidth(1);
      if(TString(obj->GetName()).Contains("hnJet") ) {
	int nbins = hobj->GetNbinsX();
	float overflow = hobj->GetBinContent(nbins+1);
	float lastbinval = hobj->GetBinContent(nbins);
	hobj->SetBinContent(nbins, overflow+lastbinval);
	hobj->GetXaxis()->SetBinLabel(nbins, "#geq4");
      }
      if (TString(obj->GetName()).Contains("hpatmet_")){
	hobj->GetYaxis()->SetTitle("events/(10 GeV)");
      }
    }
  }
  
 
  // Find out what the names of the existing histograms are
  // The histogram names are XX_YY_ZZ, where XX is the sample,
  // eg, "tt", YY is the actual name, ZZ is the final state, eg, "ee"
  TObjArray* myNames = getMyHistosNames("ttdil","ee",keep2D);
    
  bool haveRefHistograms = gROOT->FindObjectAny(Form("ref_ttdil_%s_ee",myNames->At(0)->GetName())) != 0;

  // Now loop over histograms, and make stacks
  TCanvas *c = new TCanvas();
  c->Divide(2,2);
  char* suffix[4];
  suffix[0] = "ee";
  suffix[1] = "mm";
  suffix[2] = "em";
  suffix[3] = "all";
  if (makePictures) c->Print("out/stacks.ps[");
  for (int i=0; i<myNames->GetEntries(); i++) {
     
    for (int sample=0; sample<4; sample++) {
       
       
      hist::stack(Form("st_%s_%s",myNames->At(i)->GetName(),suffix[sample]),
		  Form("^[^_]+_%s_%s$",myNames->At(i)->GetName(), suffix[sample]), false, "",orderScheme, bsmName,haveRefHistograms);
      THStack* thisStack = (THStack*) gROOT->FindObjectAny(
							   Form("st_%s_%s", myNames->At(i)->GetName(), suffix[sample]));
      THStack* refStack = 0;
      if (haveRefHistograms){
	hist::stack(Form("st_ref_%s_%s",myNames->At(i)->GetName(),suffix[sample]),
		    Form("^ref_[^_]+_%s_%s$",myNames->At(i)->GetName(), suffix[sample]), false, "",orderScheme, bsmName,haveRefHistograms);
	refStack = (THStack*) gROOT->FindObjectAny(
						   Form("st_ref_%s_%s", myNames->At(i)->GetName(), suffix[sample]));
      }
      double thisMaximum = refStack ? refStack->GetMaximum()*maxYScaleF : thisStack->GetMaximum()*maxYScaleF;
      thisStack->SetMaximum(thisMaximum);

      if(TString(myNames->At(i)->GetName()).Contains("hnJet")) {
	//	TList* histolist = thisStack->GetHists();
	//	int hatchcount = 0;
	// 	for(int j = 0; j<histolist->GetSize();j++) {
	// 	  if(TString(histolist->At(j)->GetName()).Contains("tt") ||
	// 	     TString(histolist->At(j)->GetName()).Contains("tautau") ||
	// 	     TString(histolist->At(j)->GetName()).Contains("ww") ) continue;
	// 	  hatch(histolist->At(j)->GetName(), FavoriteHatches[hatchcount]);
	// 	  hatchcount++;
	// 	}
      }
	 
      TLegend* thisLeg = hist::legend(thisStack, "f", 0, 0, 0.725, 0.605, 0.99, 0.99);
      thisLeg->SetBorderSize(1);
      c->cd(sample+1);
      if (logScale) gPad->SetLogy(); else gPad->SetLogy(0);
      double stackMax = ((TH1*)thisStack->GetHists()->At(0))->GetMaximum();
      double stackMin = ((TH1*)thisStack->GetHists()->At(0))->GetMinimum();
      thisStack->SetMinimum(stackMin);
      if (setMinZero) thisStack->SetMinimum(0);
      if (logScale && stackMin <=0) thisStack->SetMinimum(1e-2*stackMax);
      if (logScale && stackMax == 0) thisStack->SetMinimum(1e-12); 
      if ( makePictures || wait ) {
	thisStack->Draw("hist");
	string xtitle( ((TH1*)gROOT->FindObjectAny(Form("ttdil_%s_%s", myNames->At(i)->GetName(), suffix[sample])))->GetXaxis()->GetTitle());
	string ytitle( ((TH1*)gROOT->FindObjectAny(Form("ttdil_%s_%s", myNames->At(i)->GetName(), suffix[sample])))->GetYaxis()->GetTitle());
	thisStack->GetXaxis()->SetTitle(xtitle.c_str());
	thisStack->GetYaxis()->SetTitle(ytitle.c_str());
	TString hname = thisStack->GetName();
	if(hname.Contains("hnJet")) {
	  thisStack->GetXaxis()->SetLabelSize(0.075);
	  thisStack->GetYaxis()->SetLabelSize(0.05);
	  thisStack->GetXaxis()->SetTitle("#font[12]{N}_{jets}");
	  thisStack->GetXaxis()->SetTitleOffset(1.05);
	}
        if (!noLegend && sample ==3) thisLeg->Draw();
	
	TPaveText *pt1 = new TPaveText(0.1, 0.95, 0.4, 0.999, "brNDC");
	pt1->SetName("pt1name");
	pt1->SetBorderSize(0);
	pt1->SetFillStyle(0);
	
	TText *blah=0;
	if (addHistName == 0){
	  blah = pt1->AddText("");
	} else if (addHistName == 1){
	  blah = pt1->AddText(hname);
	} else if (addHistName == 2){
	  blah = pt1->AddText("CMS Preliminary");
	} else if (addHistName == 3){
	  blah = pt1->AddText("CMS Simulation");
	} else if (addHistName == 4){
	  blah = pt1->AddText("CMS Preliminary  10 pb^{-1}");
	} else if (addHistName == 5){
	  blah = pt1->AddText("CMS Simulation   10 pb^{-1}");
	} else if (addHistName == 6){
	  blah = pt1->AddText("CMS Preliminary 100 pb^{-1}");
	} else if (addHistName == 7){
	  blah = pt1->AddText("CMS Simulation  100 pb^{-1}");
	}
	blah->SetTextSize(0.05);
        pt1->Draw();
        c->Modified();
	
	c->Update();
      }
    }
    if (makePictures) {
      std::cout<<myNames->At(i)->GetName()<<std::endl;
      c->Print("out/stacks.ps");
      //       c->Print(Form("out/stacks_%d.png",i+1));
      //c->Print(Form("out/stacks_%s.png",myNames->At(i)->GetName()));
      //      c->Print(Form("out/stacks_%s.eps",myNames->At(i)->GetName()));
    }
    if (wait) {
      cout << "Enter carriage return for the next set of plots....q to quit" << endl;
      char in = getchar();
      if (in == 'q') break;
    }
  }
  if (makePictures) c->Print("out/stacks.ps]");
}
