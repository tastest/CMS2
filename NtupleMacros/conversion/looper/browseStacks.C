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




void browseStacks( bool makePictures=false, bool wait=true , bool addHistName = false, Double_t maxYScaleF = 1., 
                  bool logScale = false, bool setMinZero = true) {


  gStyle->SetOptTitle(0);
  
  hist::color("^ttdil_", kGreen);
  hist::color("^ttotr_", kYellow);
  hist::color("^DYee_", kMagenta);
  hist::color("^DYmm_", kCyan);
  hist::color("^DYtautau_", kBlack);
  hist::color("^VV_", kGray);
  hist::color("^wjets_", kViolet);
  hist::color("^QCD_", kViolet+3); //use the same as in doAll
  hist::color("^t_", kRed-3);      
  hist::color("^Vgamma_", kViolet-4);      
  hist::color("^LM0_", kBlue);

  bool keep2D=false;

  //fix the hNJet histos
  TList *list = gDirectory->GetList();
  TIterator *iter = list->MakeIterator();
  TObject *obj = 0;
  while(obj = iter->Next()) {
  
    if(TString(obj->GetName()).Contains("nJets") && obj->InheritsFrom(TH1::Class())) {
      
      int nbins = ((TH1F*)obj)->GetNbinsX();
      float overflow = ((TH1F*)obj)->GetBinContent(nbins+1);
      float lastbinval = ((TH1F*)obj)->GetBinContent(nbins);
      ((TH1F*)obj)->SetBinContent(nbins, overflow+lastbinval);
      // ((TH1F*)obj)->GetXaxis()->SetBinLabel(nbins, "#geq4");
     
    }
  }
  
    
    
    
  // Find out what the names of the existing histograms are
  // The histogram names are XX_YY_ZZ, where XX is the sample,
  // eg, "tt", YY is the actual name, ZZ is the final state, eg, "ee"
  TObjArray* myNames = getMyHistosNames("DYee","Ch0H0",keep2D);
    

  // Now loop over histograms, and make stacks
  TCanvas *c = new TCanvas();
 
  
  c->Divide(1,2);
  char* suffix[2];
  suffix[0] = "Ch0H0";
  suffix[1] = "Ch1H0";
  // suffix[2] = "Ch0H1";
//   suffix[3] = "Ch1H1";
  if (makePictures) c->Print("out/stacks.ps[");
  for (int i=0; i<myNames->GetEntries(); i++) {
 
    for (int sample=0; sample<2; sample++) {
    
       
      hist::stack(Form("%s_%s",myNames->At(i)->GetName(),suffix[sample]),
		  Form("%s_%s$",myNames->At(i)->GetName(), suffix[sample]));
      // cout << myNames->At(i)->GetName() <<endl; 
      THStack* thisStack = (THStack*) gROOT->FindObjectAny(
							   Form("%s_%s", myNames->At(i)->GetName(), suffix[sample]));
       
      thisStack->SetMaximum(thisStack->GetMaximum()*maxYScaleF);
      if(TString(myNames->At(i)->GetName()).Contains("nJets")) {
	TList* histolist = thisStack->GetHists();
	// int hatchcount = 0;
// 		for(int j = 0; j<histolist->GetSize();j++) {
// 		  if(TString(histolist->At(j)->GetName()).Contains("tt") ||
// 		     TString(histolist->At(j)->GetName()).Contains("tautau") ||
// 		     TString(histolist->At(j)->GetName()).Contains("ww") ) continue;
// 		  hatch(histolist->At(j)->GetName(), FavoriteHatches[hatchcount]);
// 		  hatchcount++;
// 		}
      }
	 
     
      TLegend* thisLeg = hist::legend(thisStack, "lpf", 0, 0, 0.75, 0.65, 0.99, 0.99);
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
	string xtitle( ((TH1*)gROOT->FindObjectAny(Form("DYee_%s_%s", myNames->At(i)->GetName(), suffix[sample])))->GetXaxis()->GetTitle());
	string ytitle( ((TH1*)gROOT->FindObjectAny(Form("DYee_%s_%s", myNames->At(i)->GetName(), suffix[sample])))->GetYaxis()->GetTitle());
	thisStack->GetXaxis()->SetTitle(xtitle.c_str());
	thisStack->GetYaxis()->SetTitle(ytitle.c_str());
	TString hname = thisStack->GetName();

	if(hname.Contains("nJets")) {
	 
	  thisStack->GetXaxis()->SetLabelSize(0.075);
	  thisStack->GetYaxis()->SetLabelSize(0.05);
	  thisStack->GetXaxis()->SetTitle("N_{jets}");
	}
	
        thisLeg->Draw();

	TPaveText *pt1 = new TPaveText(0.1, 0.95, 0.4, 0.999, "brNDC");
	pt1->SetName("pt1name");
	pt1->SetBorderSize(0);
	pt1->SetFillStyle(0);
	
	TText *blah;
	if (addHistName) blah = pt1->AddText(hname);
	else blah = pt1->AddText("CMS Preliminary");
	blah->SetTextSize(0.05);
        pt1->Draw();
        c->Modified();
	
	c->Update();
      }
    }
    if (makePictures) {
      c->Print("out/stacks.ps");
      //c->Print(Form("out/stacks_%d.png",i+1));
      //c->Print(Form("out/stacks_%s.png",myNames->At(i)->GetName()));
      c->Print(Form("out/stacks_%s.eps",myNames->At(i)->GetName()));
    }
    if (wait) {
      cout << "Enter carriage return for the next set of plots....q to quit" << endl;
      char in = getchar();
      if (in == 'q') break;
    }
  }
  if (makePictures) c->Print("out/stacks.ps]");
}
