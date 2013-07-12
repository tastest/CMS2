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
#include "TSystem.h"
#include "TArrow.h"
#include "TLatex.h"
#include <iostream>
#include <algorithm>
#include <set>
#include "CommonFunctions.C"
#include "getMyHistosNames.h"
#include "histtools.h"

using namespace std;

// Make the stacks and then browse them interactive
// if (makePictures==true), saves all plots to out/stacks.ps
// Also: make out_stacks_xx.png where xx is the page number 
// on the stack.ps file
//
// keep2D=false to skip the annoying 2D histograms
//


double GetMinimum(const vector<TH1D*> &v_hists);
TLegend* makeLegend(const vector<TH1D*> &v_hists, vector<TString> v_legEntries, bool drawLogY, TString histName);
TPaveText *getPaveText(const vector<TH1D*> &v_hists, int i_channel, float lumi, bool drawFullErrors); //call this after makeLegend please
TH1D* getDiffHist(TH1D* h1, TH1D* h2);



void browseStacks(vector<TString> v_samples, vector<Color_t> v_colors, 
  TString outfile, vector<TString> v_legEntries, bool drawLogY = false, 
  vector<Style_t> v_style = vector<Style_t>(), bool drawFullErrors = false, 
bool drawDiffs = true, float lumi=1.0, int rebin=1, bool allhypOnly = false) {

  if(v_samples.size() != v_colors.size()) {
    cout << "Number of entries in the vector of samples is not the same as the number of entries in the vector of Color_t" << endl;
    return;
  }

  if(v_style.size()!=0 && v_style.size() != v_colors.size()) {
    cout << "Number of entries in the vector of styles is not the same as the number of entries in the vector of samples" << endl;
    return;
  }

  gStyle->SetOptTitle(0);
  gStyle->SetHatchesSpacing(2);
  gStyle->SetHatchesLineWidth(0.01);

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  // Find out what the names of the existing histograms are
  // The histogram names are XX_YY_ZZ, where XX is the sample,
  // eg, "tt", YY is the actual name, ZZ is the final state, eg, "ee"
  //exclude data if it is found
  bool keep2D = false;
  //TObjArray* myNames = getMyHistosNames(v_samples.at(0).Data(),"allj_ee",keep2D);  
  TObjArray* myNames = getMyHistosNames(v_samples.at(0).Data(),"ee",keep2D);  
  if(myNames->GetSize() == 0) {
    cout << "could not find the sample at the first point in the vector of samples. Exiting." << endl;
    return;
  }

  // Now loop over histograms, and make stacks
  vector<string> v_channel;
  v_channel.push_back("ee");
  v_channel.push_back("mm");
  v_channel.push_back("em");
  v_channel.push_back("all");

  TLine *l; //for drawing the diffs
  for (int i=0; i<myNames->GetEntries(); i++) {//histos loop

    TString plot(myNames->At(i)->GetName());
    TCanvas *c = new TCanvas();
    if(!allhypOnly) c->Divide(2,2);  

    vector<TH1D*> v_diffs;

    int first_channel_canvas = 0;
    if(allhypOnly) first_channel_canvas = 3;

    for (int i_channel=first_channel_canvas; i_channel<4; i_channel++) {//loop over the channels
      vector<TH1D*> v_hists;
      TH1D *hdata = NULL;
      TH1D *hwprime = NULL;

      for(unsigned int i_prefix = 0; i_prefix < v_samples.size(); i_prefix++) {

        TString histoName = v_samples.at(i_prefix) + "_" + myNames->At(i)->GetName() + "_" + v_channel.at(i_channel);
        TObject *obj = gDirectory->Get(histoName.Data());
        if(!obj->InheritsFrom(TH1::Class())) 
          continue;

        TH1D *htemp = dynamic_cast<TH1D*>(obj);

        htemp->SetFillColor(v_colors.at(i_prefix));
        htemp->SetLineColor(kBlack);
        if(plot.Contains("hnJet") ||
          plot.Contains("htcmet") ||
          plot.Contains("bTag") ||
          plot.Contains("dilMass") || 
          plot.Contains("hpfmet") || 
          plot.Contains("mt2")  ||
        plot.Contains("1Dmasscut") ) {

          string xtitle;
          string ytitle;
          if(plot.Contains("hnJet")) {
            xtitle = "Number of jets";
            ytitle = "Events";
          } else if(plot.Contains("htcmet")) {
      //xtitle = "Missing transverse energy [GeV]";
            xtitle = "tc Missing transverse energy [GeV]";
            ytitle = "Events/(10 GeV)";
          } else if (plot.Contains("hpfmet")) {
            xtitle = "PF Missing transverse energy [GeV]";
            ytitle = "Events/(10 GeV)";
          } else if(plot.Contains("hbTagsimpleSecVtxHighEff")) {
            xtitle = "Number of SSV b-tagged jets";
            ytitle = "Events";
            const char *jetbins[4] = {"0", "1", "2", "#geq 3"};
            for(int k = 0; k<4; k++) 
              htemp->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
            htemp->GetXaxis()->SetLabelSize(0.07);
          } else if(plot.Contains("hbTagtkCountHighEff")) {
            xtitle = "Number of TCHE b-tagged jets";
            ytitle = "Events";
            const char *jetbins[4] = {"0", "1", "2", "#geq 3"};
            for(int k = 0; k<4; k++) 
              htemp->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
            htemp->GetXaxis()->SetLabelSize(0.07);
          } else if(plot.Contains("dilMass")) {
            xtitle = "Dilepton mass [GeV/c^{2}]";
            ytitle = "Events/(5 GeV/c^{2})";
          } else if(plot.Contains("mt2_")) {
            xtitle = "MT2";
            ytitle = "";
          } else if(plot.Contains("mt2J_")) {
            xtitle = "MT2J";
            ytitle = "";
          } else if(plot.Contains("massltb1Dmasscut")) {
            xtitle = "M_{l1b1} (GeV/c^{2}) for M_{l2b2} > 170 GeV/c^{2}";
            ytitle = "Events";
          } else if(plot.Contains("massllb1Dmasscut")) {
            xtitle = "M_{l2b2} (GeV/c^{2}) for M_{l1b1} > 170 GeV/c^{2}";
            ytitle = "Events";
          }

          htemp->GetXaxis()->SetTitle(xtitle.c_str());
          htemp->GetYaxis()->SetTitle(ytitle.c_str());	  

        }
        string ytitle = "Events";
        htemp->GetYaxis()->SetTitle(ytitle.c_str());
        htemp->GetYaxis()->SetTitleOffset(1.5);
        htemp->GetXaxis()->SetTitleSize(0.045);
        htemp->GetYaxis()->SetTitleSize(0.045);


        if(v_style.size() > 0) {
          htemp->SetFillStyle(v_style.at(i_prefix));
          if(v_samples.at(i_prefix) == "data") {
            htemp->SetMarkerStyle(v_style.at(i_prefix));
            htemp->SetMarkerColor(v_colors.at(i_prefix));
          }
          if(v_samples.at(i_prefix).Contains("wprime")) {
      //htemp->SetMarkerStyle(v_style.at(i_prefix));
            htemp->SetLineColor(kBlue+2);
      //htemp->SetMarkerSize(0.5);
          }
        }

  //if(plot.Contains("_0j") || plot.Contains("_1j") || plot.Contains("_2j")|| plot.Contains("_allj") ) htemp->Rebin(rebin);
        if( htemp->GetNbinsX() % rebin == 0 && (plot.Contains("_0j") || plot.Contains("_1j") || plot.Contains("_2j")|| plot.Contains("_allj")) ) {
          htemp->Rebin(rebin);
          } else if ( plot.Contains("httpT") ) {
    // special rebinnig for httpT, choose numberOfBins = wanted number of bins for ttpT >= 0 and add 1
            const int numberOfBins = 21;
            Double_t rebinning[numberOfBins];
            rebinning[0] = -4.;
            rebinning[1] = 0.;
            for ( int index = 0; index < numberOfBins-1; ++index ) {
              rebinning[index+2] = 400./(numberOfBins-1)*(index+1);
            }
            htemp = dynamic_cast<TH1D*>(htemp->Rebin(numberOfBins,"",rebinning));    
          }
  //don't add the data histogram to the stack
          if(v_samples.at(i_prefix) == "data") {
            hdata = htemp;	  
            continue;
          }

  //don't add the wprime histogram to the stack
          if(v_samples.at(i_prefix).Contains("wprime")) {
            hwprime = htemp;	  
            continue;
          }

          if(i_prefix==0) {
            v_hists.push_back(htemp);
            continue;
          }

          htemp->Add(v_hists.back());
          v_hists.push_back(htemp);
        }//prefix loop			
        if(hwprime != NULL)
          v_hists.push_back(hwprime);	
        if(hdata != NULL)
          v_hists.push_back(hdata);

        c->cd(i_channel+1);

      //now set the Minimum if we want Log Scale
        float min = 0;
      //if(drawLogY && i_channel != 2) 
        if(drawLogY) 
          min = GetMinimum(v_hists);


        if(drawLogY && TString(myNames->At(i)->GetName()).Contains("1Dmasscut")) min = 0.01;     //for 1.09 fb-1


      //set the minimum, before we pass the vector of hists to the legend function
        for(vector<TH1D*>::iterator it = v_hists.begin(); it != v_hists.end(); it++) 
          (*it)->SetMinimum(min);


        TLegend *leg = NULL;
        if(i_channel == 3) {
          leg = makeLegend(v_hists, v_legEntries, drawLogY, TString(myNames->At(i)->GetName()));

        } //else {
          float max = 0;
          for(vector<TH1D*>::iterator it = v_hists.begin(); it != v_hists.end(); it++) {
            if((*it)->GetMaximum() > max)
              max = (*it)->GetMaximum();
          }



          for(vector<TH1D*>::iterator it = v_hists.begin(); it != v_hists.end(); it++) {
    //if(drawLogY && i_channel != 2 )
            if(drawLogY)
              (*it)->SetMaximum(50*max);
            else 
              (*it)->SetMaximum(1.45*max);
    //if(!drawLogY && TString((*it)->GetName()).Contains("hnJet")) 
    //(*it)->SetMaximum(6);
          }
     // end of commented out else  }


          TPad *p1, *p2;
          if(drawLogY && !drawDiffs) {
            p1 = new TPad("p1", "stack", 0.0, 0.0, 1, 1.);
            p1->Draw();
            p1->cd();

          }
          if(drawDiffs) {
            p1 = new TPad("p1", "stack", 0.0, 0.31, 1, 1.);
            p2 = new TPad("p2", "diff", 0.0, 0.0, 1, 0.3);
            p1->Draw();
            p2->Draw();
            p1->cd();
          }


          if(drawDiffs && hdata != NULL) {
            p2->cd();
            int tempnbins = hdata->GetNbinsX();
            TString s_hname = "diff_";
            TH1D *h_diff = new TH1D((s_hname +  myNames->At(i)->GetName() + "_" + v_channel.at(i_channel)).Data(), "", tempnbins,
              hdata->GetXaxis()->GetBinLowEdge(1),
              hdata->GetXaxis()->GetBinUpEdge(tempnbins));
            l = new TLine(hdata->GetXaxis()->GetBinLowEdge(1), 0.0, 
              hdata->GetXaxis()->GetBinUpEdge(tempnbins), 0.0);
            h_diff->TH1D::Sumw2();	
            for(int tempbin = 1; tempbin < hdata->GetNbinsX()+1; tempbin++) {
              double mc = (v_hists.at(v_hists.size()-2))->GetBinContent(tempbin);	  
              double data = hdata->GetBinContent(tempbin);
              double mcerr = (v_hists.at(v_hists.size()-2))->GetBinError(tempbin);
              double dataerr = hdata->GetBinError(tempbin);
              float diff = data - mc;
              float err2 = pow(mcerr*data/mc/mc,2) + pow(dataerr/mc,2);
              if(mc < 1e-6) 
                continue;
              h_diff->SetBinContent(tempbin, (data-mc)/mc);
              h_diff->SetBinError(tempbin, sqrt(err2));		    
              h_diff->GetXaxis()->SetBinLabel(tempbin, hdata->GetXaxis()->GetBinLabel(tempbin));
            }


            h_diff->GetXaxis()->SetTitle("");
            h_diff->GetXaxis()->SetLabelSize(0.);

            h_diff->GetYaxis()->SetTitle("(Data - MC)/MC");
            h_diff->GetYaxis()->SetTitleFont(hdata->GetYaxis()->GetTitleFont());
            h_diff->GetYaxis()->SetTitleOffset(0.5);
            h_diff->GetYaxis()->SetTitleSize(0.1);
            h_diff->GetYaxis()->SetLabelSize(0.105);
            h_diff->GetYaxis()->SetLabelFont(hdata->GetYaxis()->GetLabelFont());

            h_diff->SetMarkerSize(0.8);

            h_diff->Draw("Pe");

            h_diff->SetMinimum(-0.5);
            h_diff->SetMaximum(0.5);

            h_diff->Draw("Pesames");
            l->Draw();
            c->Modified();
            c->Update();
            p1->cd();	
          }

      //now we gotta draw, in reverse order
          for(vector<TH1D*>::reverse_iterator it = v_hists.rbegin(); it != v_hists.rend(); it++) {

  //do not draw if the histogram is empty...1e-6 should be good enough
            if(drawLogY && (*it)->GetMaximum() < 1e-6)
              continue;

            if(i_channel == 3) {  
              float max = 1.1*leg->GetY2();
              if(drawLogY)
                max = 5*leg->GetY2();
              if(drawLogY && TString(myNames->At(i)->GetName()).Contains("1Dmasscut")) max = 500.;  //for 1.09 fb-1
              if(v_hists.back()->GetMaximum() < max) 
                (*it)->SetMaximum(max);
            }	

            if(TString((*it)->GetName()).Contains("data")) {
              hdata = (*it);
              if(it == v_hists.rbegin()) 
                (*it)->Draw("Pe");	 
              else 	  
                (*it)->Draw("Pesame");

            }	
  //else if(TString((*it)->GetName()).Contains("wprime")) {
  //		hwprime = (*it);
  //		if(it == v_hists.rbegin()) 
  //    		(*it)->Draw("P");	 
  //		else 	  
  //    		(*it)->Draw("Psame");	
  //}	
            else {

              if(it == v_hists.rbegin()) 
                (*it)->Draw("hist");	 
              else 	  
                (*it)->Draw("histsame");
            }
          }//loop over MC

          if(hwprime != NULL) {
            hwprime->Draw("histsame");
          }

          if(hdata != NULL) {
            hdata->Draw("Pesame");
          }

      //now set the Log scale, if desired
      //if(drawLogY && i_channel != 2) 
          if(drawLogY) 
            p1->SetLogy(1);


      //now get the top histogram (sum of all) and draw the errors
          if(drawFullErrors) {

            TH1D *h_sumBackGrounds = NULL;
            for(vector<TH1D*>::reverse_iterator it = v_hists.rbegin(); it != v_hists.rend(); it++) {
              if(TString((*it)->GetName()).Contains("data"))
                continue;
              if(TString((*it)->GetName()).Contains("ttdil"))
                continue;
              if(h_sumBackGrounds == NULL)
                h_sumBackGrounds = (TH1D*)((*it)->Clone());
            }


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


          if(TString(myNames->At(i)->GetName()).Contains("1Dmasscut")) {
            TLine *line170 = new TLine();                                                                                                                                                                                                                                                                                                                                            
            line170->SetLineColor(kGreen+3);                                                                                                                                                                                                                                                                                                                                          
            line170->SetLineWidth(2);                                                                                                                                                                                                                                                                                                                                             
            if(drawLogY && hdata != NULL) line170->DrawLine(170.,0.0,170.,0.1*hdata->GetMaximum());                                                                                                                                                                                                                                                                                                     
            else if(hdata != NULL) line170->DrawLine(170.,0.0,170.,hdata->GetMaximum());

            TArrow *arrow170 = new TArrow();
            arrow170->SetLineWidth(2);
            arrow170->SetLineColor(kGreen+3);
            arrow170->SetFillColor(kGreen+3);        
            arrow170->DrawArrow(170.,3.5,340.,3.5, 0.03, "|>");  //for 1.09 fb-1

            TLatex *t170=new TLatex();                                                                                                                                                                                                                                                                                                                                                  
            t170->SetTextSize(0.034);                                                                                                                                                                                                                                                                                                                                                    
            t170->SetTextColor(kGreen+3);
            t170->DrawLatex(180,4.5,"Signal Region");        //for 1.09 fb-1

          }


        }//channel loop
        c->Modified();
        c->Update();

    //cout << myNames->At(i)->GetName() << endl;
        TString tempstring = outfile.ReplaceAll(".root", "");
        if(tempstring.Contains("_PU"))
          tempstring = tempstring.ReplaceAll("hist", "histwPU");
        if(drawFullErrors)
          tempstring = tempstring + "_withErrors";

        if(i != myNames->GetEntries() - 1) 
          c->Print((tempstring + ".pdf(").Data());
        else
          c->Print((tempstring + ".pdf)").Data());

        c->Print(("results/" + plot + ".png").Data());

      }//histos loop

    }




    double GetMinimum(const vector<TH1D*> &v_hists) {

  /*
  TH1* h = NULL; //get the first non-empty histogram 
  for(unsigned int i = 0; v_hists.size(); i++) {
  h = v_hists.at(i);
  if(h->Integral() > 0)
  break;
  }


if(h->GetMinimum()  < 1e-2)
return 1e-2;
else 
return 
0.5*(h->GetMinimum());
*/
if(v_hists.size() == 0 ) {
  cout << "There are no histograms in the vector of histograms" << endl;
  cout << "In function " << __FUNCTION__ << ". Exiting" << endl;
  gSystem->Exit(1);
}

  //the histogram thats drawn last
TH1D *h = v_hists.at(0);
float min = 99999999;
for(int i = 1; i < h->GetNbinsX()+1; i++) {

  if( h->GetBinContent(i) < min && h->GetBinContent(i) > 0 )
    min = h->GetBinContent(i);

}

return 0.5*min;

}



TLegend* makeLegend(const vector<TH1D*> &v_hists, vector<TString> v_legEntries, bool drawLogY, TString histName) {

  //Prefer to draw the Legend on the right half of the 
  //canvas, so only look at the right half
  int nbins = v_hists.at(0)->GetNbinsX();
  TH1D *hdata = NULL;
  TH1D *hmax = NULL;

  for(vector<TH1D*>::const_reverse_iterator rit = v_hists.rbegin(); rit != v_hists.rend(); rit++) {
    if(TString((*rit)->GetName()).Contains("data")) {
      hdata = *rit;
      break;
    }
  }

  for(vector<TH1D*>::const_reverse_iterator rit = v_hists.rbegin(); rit != v_hists.rend(); rit++) {
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
  if(hdata != NULL){
    for(int bin = lowBinX; bin < highBinX+1; bin++) {
      if(hmax->GetBinContent(bin) >= hdata->GetBinContent(bin)) {
        if(max < hmax->GetBinContent(bin))
          max = hmax->GetBinContent(bin);
      } else {
        if(max < hdata->GetBinContent(bin))
          max = hdata->GetBinContent(bin);
      }
    }
  }

  //cout << "max: " << max << endl;
  if(hdata != NULL){
    rangeY = hdata->GetMinimum();
    if(hdata->GetMaximum() > hmax->GetMinimum())
      rangeY = hdata->GetMaximum() - rangeY;
    else 
      rangeY = hmax->GetMaximum() - rangeY;
  }
  else rangeY = hmax->GetMaximum() - hmax->GetMinimum() ; 

  if(drawLogY) {  	
    lowY = 5*max;
    if(GetMinimum(v_hists)>max) lowY = 10*GetMinimum(v_hists);
    if(histName.Contains("1Dmasscut")) lowY = 1.; //for 1.09 fb-1
    highY = lowY + 10*rangeY;
    if(histName.Contains("1Dmasscut")) highY = 300.; //for 1.09 fb-1
  } else {
    lowY = 1.2*max;
    highY = lowY + 0.3*rangeY;
  }

  //cout << "lowY, highY: "<< lowY << "," <<  highY << endl;
  TLegend *leg;
  if(drawLogY)
    leg = new TLegend(lowX,lowY,highX,highY, "", "br"); 
  else if(histName.Contains("lepAzimAsym2_") || histName.Contains("hlepAngleBetween_"))
    leg = new TLegend(0.58, 0.65, 0.78, 0.93, "", "brNDC");
  else leg = new TLegend(0.72, 0.65, 0.92, 0.93, "", "brNDC");
  if((histName.Contains("nJet") || histName.Contains("bTag")) && drawLogY && hdata != NULL ) {
    int temp = hdata->GetNbinsX();
    leg = new TLegend(hdata->GetNbinsX()-1.5, lowY, temp - 0.45*hdata->GetBinWidth(temp), highY, "", "br");

  }

  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  //leg->SetFillStyle(1001);
  leg->SetFillStyle(0);     

  if(v_hists.size() != v_legEntries.size()) {
    cout << "the number of entries in the legend vector are not the same as the number"
      << " of entries in the hists vector. Returning a null TLegend. " << endl;
    return NULL;
  }
  vector<TH1D*>::const_reverse_iterator ritH  = v_hists.rbegin();
  vector<TString>::const_reverse_iterator ritE  = v_legEntries.rbegin();
  for(; ritH != v_hists.rend(); ritH++, ritE++) {
    if(*ritE == "Data")
      leg->AddEntry(*ritH, *ritE, "P");
  //else if( TString(*ritE).Contains("GeV") )
  //  leg->AddEntry(*ritH, *ritE, "P");
    else
      leg->AddEntry(*ritH, *ritE, "f");
  }




  return leg;
}





TPaveText *getPaveText(const vector<TH1D*> &v_hists, int i_channel, float lumi, bool drawFullErrors) {

  if(v_hists.size() == 0)
    return NULL;



  //for the MET plot, after the Zveto and the 2jet cut
  TPaveText *pt1 = new TPaveText(0.60, 0.77, 0.80, 0.90, "brNDC");
  if(i_channel == 3) {
    //if(TString(v_hists.at(0)->GetName()).Contains("hnJet")) 
      //pt1 = new TPaveText(0.33, 0.77, 0.56, 0.90, "brNDC");
    //else if(TString(v_hists.at(0)->GetName()).Contains("hdilMass"))
      //pt1 = new TPaveText(0.17, 0.77, 0.4, 0.90, "brNDC");
    //else 
    pt1 = new TPaveText(0.22, 0.77, 0.45, 0.90, "brNDC");
  }
  pt1->SetName("pt1name");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);


  TText *blah;
  if(drawFullErrors)
    blah = pt1->AddText("CMS");
  else
    blah = pt1->AddText("CMS Preliminary");
    //blah = pt1->AddText("");
  blah->SetTextSize(0.036);
  blah->SetTextAlign(11);  
  TString temp = formatFloat(lumi,"%6.1f");
  temp.ReplaceAll(" " , "" );
  temp = temp + TString(" fb^{-1} at   #sqrt{s}=7 TeV");
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


TH1D* getDiffHist(TH1D* h1, TH1D* h2){

  TH1D* hout = (TH1D*) h1->Clone(Form("%s_clonediff",h1->GetName()));

  for(int ibin = 1 ; ibin <= h1->GetNbinsX() ; ibin++){

    float val = h1->GetBinContent(ibin) - h2->GetBinContent(ibin);
    float err = sqrt(pow(h1->GetBinError(ibin),2)+pow(h2->GetBinError(ibin),2));
    hout -> SetBinContent( ibin, val );
    hout -> SetBinError(   ibin, err );
  }


  hout->SetMarkerStyle(24);
  hout->SetMarkerSize(1.3);
  return hout;
}
