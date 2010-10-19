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
#include "CommonFunctions.C"
#include "getMyHistosNames.C"
#include "histtools.C"


using namespace std;

// Make the stacks and then browse them interactive
// if (makePictures==true), saves all plots to out/stacks.ps
// Also: make out_stacks_xx.png where xx is the page number 
// on the stack.ps file
//
// keep2D=false to skip the annoying 2D histograms
//


double GetMinimum(const vector<TH1F*> &v_hists);
TLegend* makeLegend(const vector<TH1F*> &v_hists, vector<TString> v_legEntries, bool drawLogY, TString histName);
TPaveText *getPaveText(const vector<TH1F*> &v_hists, int i_channel, float lumi, bool drawFullErrors); //call this after makeLegend please
void browseStacks(vector<TString> v_samples, vector<Color_t> v_colors, 
		  TString outfile, vector<TString> v_legEntries, bool drawLogY = false, 
		  vector<Style_t> v_style = vector<Style_t>(), bool drawFullErrors = false, float lumi = 1.0) {
		     

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
	if(! obj->InheritsFrom(TH1::Class())) 
	  continue;
	
	TH1F *htemp = dynamic_cast<TH1F*>(obj);

	htemp->SetFillColor(v_colors.at(i_prefix));
	htemp->SetLineColor(kBlack);
	TString plot(myNames->At(i)->GetName());
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
	  htemp->GetYaxis()->SetTitleOffset(1.5);
	  htemp->GetXaxis()->SetTitleSize(0.045);
	  htemp->GetYaxis()->SetTitleSize(0.045);
	    
	}
	
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

	if(i_channel == 3) {  
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
	if(i_channel == 3) {
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
      
      
      if(i_channel == 3) 
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
      0.5*(h->GetMinimum());

}
      
  

TLegend* makeLegend(const vector<TH1F*> &v_hists, vector<TString> v_legEntries, bool drawLogY, TString histName) {

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
      
  //cout << "max: " << max << endl;
  rangeY = hdata->GetMinimum();
  if(hdata->GetMaximum() > hmax->GetMinimum())
    rangeY = hdata->GetMaximum() - rangeY;
  else 
    rangeY = hmax->GetMaximum() - rangeY;

  if(drawLogY) {
    lowY = 5*max;
    highY = lowY + 10*rangeY;
  } else {
    lowY = 1.2*max;
    highY = lowY + 0.3*rangeY;
  }
  
  //cout << "lowY, highY: "<< lowY << "," <<  highY << endl;
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
