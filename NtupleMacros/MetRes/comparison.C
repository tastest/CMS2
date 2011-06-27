
//#include <vector>
#include <iostream>
#include "TH1F.h"
//#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
//#include "TROOT.h"

using namespace std;

void over_save( TFile* Res223, TFile* Res2210, TString title, int color10 = 2 );
void over_save( TH1F* h3, TH1F* h10, bool log = false, bool logname = false, int color10 = 2 );
void over_save( TH1F* h1, TH1F* h2, TH1F* h3, TString savename = "", bool log = false, int color3 = 2 );
void over_save( TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4, TString savename = "", bool log = false, int color3 = 2 );

void comparison()
{
  TFile* Res223  = TFile::Open("Results223.root","READ");
  TFile* Res2210 = TFile::Open("Results2210.root");

  over_save( Res223, Res2210, "ww_ltPt-N-1_all" );
  over_save( Res223, Res2210, "ww_llPt-N-1_all" );
  over_save( Res223, Res2210, "ww_dilMass-N-1_all" );
  over_save( Res223, Res2210, "ww_JPTJetPt-N-1_all" );
  over_save( Res223, Res2210, "ttbar_JPTJetPt-N-1_all" );
  over_save( Res223, Res2210, "ww_metTrkCorr-N-1_all" );
  //el_rel_iso,mu_rel_iso, for ww+wjets,elPt,muPt
  over_save( Res223, Res2210, "ww_eleRelIso-N-1_all" );
  over_save( Res223, Res2210, "wjets_eleRelIso-N-1_all" );
  over_save( Res223, Res2210, "ww_muRelIso-N-1_all" );
  over_save( Res223, Res2210, "wjets_muRelIso-N-1_all" );
  over_save( Res223, Res2210, "ww_elPt-N-1_all" );
  over_save( Res223, Res2210, "wjets_elPt-N-1_all" );  

}

 //color10 = red by default
void over_save( TH1F* h3, TH1F* h10, bool log, bool logname, int color10 ) {

  //gStyle->SetTitleSize(0); //use legend instead
  //gPad->SetTitleSize(0);
  TString h3title = h3->GetTitle();
  TString h10title = h10->GetTitle();
  //h3->SetTitleSize(0); 
  //h10->SetTitleSize(0);
  h3-> SetTitle("");
  h10->SetTitle("");
  
  int blue = 4;

  TCanvas* c = new TCanvas(h3->GetName(), h3->GetName());
  //TLegend* leg = new TLegend(0.8, 0.6, 0.98, 0.75);
  //TLegend* leg = new TLegend(0.2, 0.85, 0.8, 0.98);
  TLegend* leg = new TLegend(0.4, 0.88, 0.9, 0.99);
  if( log )
	gPad->SetLogy();
  h10->SetStats(0);
  h3->SetStats(0);
  //h10->SetTitle(h10->GetName());
  //h3->SetTitle(h10->GetName());
  h10->SetMarkerColor( color10 );
  h10->SetLineColor( color10 );
  h3->SetMarkerColor( blue );
  h3->SetLineColor( blue );

  if( h3->GetMaximum() > h10->GetMaximum() ) {
	h3->Draw();
	h10->Draw("same");
  }
  else {
	h10->Draw();
	h3->Draw("same");
  }

  //h3->SetFillColor( blue );
  //h10->SetFillColor( color10 );
  leg->SetFillColor(0);
  leg->AddEntry(h3, h3->GetName());
  leg->AddEntry(h10, h10->GetName());
  leg->Draw();
  TString name;
  if( logname )
	name = (TString)h3->GetName()+"_log_comp.png";
  else
	name = (TString)h3->GetName()+"_comp.png";

  c->SaveAs(name);
  //c->SaveAs((TString)h3->GetName()+".eps");
  h3->SetTitle(h3title);
  h10->SetTitle(h10title);
}

void over_save( TFile* Res223, TFile* Res2210, TString title, int color10 ) {

  TH1F *h3  = (TH1F*)Res223->Get( title );
  TH1F *h10 = (TH1F*)Res2210->Get( title );

  over_save( h3, h10 );

}


//overloaded for three hists
void over_save( TH1F* h1, TH1F* h2, TH1F* h3, TString savename, bool log, int color3 ) {
  //cout << "starting the triple overlay" << endl;
  //gStyle->SetTitleSize(0); //use legend instead
  //gPad->SetTitleSize(0);
  TString h1title = h1->GetTitle();
  TString h2title = h2->GetTitle();
  TString h3title = h3->GetTitle();
  //h3->SetTitleSize(0); 
  //h10->SetTitleSize(0);
  //TString title = "Electron Component Isolation Efficiencies";
  TString title = "Muon Component Isolation Efficiencies";
  h1->SetTitle(title);
  h2->SetTitle(title);
  h3->SetTitle(title);
  h1->GetXaxis()->SetRangeUser(0,1.1);
  h2->GetXaxis()->SetRangeUser(0,1.1);
  h3->GetXaxis()->SetRangeUser(0,1.1);
  h1->GetYaxis()->SetRangeUser(0,1.0);
  h2->GetYaxis()->SetRangeUser(0,1.0);
  h3->GetYaxis()->SetRangeUser(0,1.0);
  h1->GetXaxis()->SetTitle("dR");
  h2->GetXaxis()->SetTitle("dR");
  h3->GetXaxis()->SetTitle("dR");
  h1->GetXaxis()->SetTitleOffset(0.4);
  h2->GetXaxis()->SetTitleOffset(0.4);
  h3->GetXaxis()->SetTitleOffset(0.4);
  
  int blue = 4;
  int green = 3;

  gStyle->SetPadTopMargin   ( 0.08);//was 0.08
  TCanvas* c = new TCanvas(h3->GetName(), h3->GetName());
  //TLegend* leg = new TLegend(0.8, 0.6, 0.98, 0.75);
  //TLegend* leg = new TLegend(0.2, 0.85, 0.8, 0.98);
  //TLegend* leg = new TLegend(0.4, 0.8, 0.9, 0.99);
  TLegend* leg = new TLegend(0.35, 0.17, 0.9, 0.34);
  if( log )
	gPad->SetLogy();
  h1->SetStats(0);
  h2->SetStats(0);
  h3->SetStats(0);
  //h1->SetTitle(h10->GetName());
  //h3->SetTitle(h10->GetName());
  h1->SetMarkerColor( color3 );
  h1->SetLineColor( color3 );
  h2->SetMarkerColor( blue );
  h2->SetLineColor( blue );
  h3->SetMarkerColor( green );
  h3->SetLineColor( green );

  if( h1->GetMaximum() > h2->GetMaximum() && h1->GetMaximum() > h3->GetMaximum() ) {
	h1->Draw();
	h2->Draw("same");
	h3->Draw("same");
  }
  else if( h2->GetMaximum() > h1->GetMaximum() && h2->GetMaximum() > h3->GetMaximum() ) {
	h2->Draw();
	h1->Draw("same");
	h3->Draw("same");
  }
  else {
	h3->Draw();
	h1->Draw("same");
	h2->Draw("same");
  }

  //h3->SetFillColor( blue );
  //h10->SetFillColor( color3 );
  leg->SetFillColor(0);
  //leg->AddEntry(h1, h1->GetName());
  //leg->AddEntry(h2, h2->GetName());
  //leg->AddEntry(h3, h3->GetName());
  leg->AddEntry(h1, "Track isolation efficiency");
  leg->AddEntry(h2, "Ecal isolation efficiency");
  leg->AddEntry(h3, "Hcal isolation efficiency");
  leg->Draw();
  //if( savename == "" ) 
  //if( savename.IsNull() ) 
  if( savename.Length() == 0 ) {
	//cout << "length is zero:" << savename << "::" << (TString)h1->GetName() << endl;
	c->SaveAs((TString)h1->GetName()+"_comp.png");
  }
  else {
	//cout << "length is not zero:" << savename << "::" << (TString)h1->GetName() << endl;
	c->SaveAs(savename+".png");
  }
  //c->SaveAs((TString)h1->GetName()+".eps");
  //need to reset these b'c they're pointers...
  h1->SetTitle(h1title);
  h2->SetTitle(h2title);
  h3->SetTitle(h3title);

  gStyle->SetPadTopMargin   ( 0.08);//was 0.08
}


//overloaded for four hists
void over_save( TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4, TString savename, bool log, int color3 ) {
  //cout << "starting the triple overlay" << endl;
  //gStyle->SetTitleSize(0); //use legend instead
  //gPad->SetTitleSize(0);
  TString h1title = h1->GetTitle();
  TString h2title = h2->GetTitle();
  TString h3title = h3->GetTitle();
  TString h4title = h4->GetTitle();
  //h3->SetTitleSize(0); 
  //h10->SetTitleSize(0);
  h1->SetTitle("");
  h2->SetTitle("");
  h3->SetTitle("");
  h4->SetTitle("");
  
  int green = 3;
  int blue = 4;
  int pink = 6;

  gStyle->SetPadTopMargin   ( 0.15);//was 0.08
  TCanvas* c = new TCanvas(h3->GetName(), h3->GetName());
  //TLegend* leg = new TLegend(0.8, 0.6, 0.98, 0.75);
  //TLegend* leg = new TLegend(0.2, 0.85, 0.8, 0.98);
  TLegend* leg = new TLegend(0.4, 0.8, 0.9, 0.99);
  if( log )
	gPad->SetLogy();
  h1->SetStats(0);
  h2->SetStats(0);
  h3->SetStats(0);
  h4->SetStats(0);
  //h1->SetTitle(h10->GetName());
  //h3->SetTitle(h10->GetName());
  h1->SetMarkerColor( color3 );
  h1->SetLineColor( color3 );
  h2->SetMarkerColor( blue );
  h2->SetLineColor( blue );
  h3->SetMarkerColor( green );
  h3->SetLineColor( green );
  h4->SetMarkerColor( pink );
  h4->SetLineColor( pink );

  if( h1->GetMaximum() > h2->GetMaximum() && h1->GetMaximum() > h3->GetMaximum() && h1->GetMaximum() > h4->GetMaximum() ) {
	h1->Draw();
	h2->Draw("same");
	h3->Draw("same");
	h4->Draw("same");
  }
  else if( h2->GetMaximum() > h1->GetMaximum() && h2->GetMaximum() > h3->GetMaximum() && h2->GetMaximum() > h4->GetMaximum() ) {
	h2->Draw();
	h1->Draw("same");
	h3->Draw("same");
	h4->Draw("same");
  }
  else if( h3->GetMaximum() > h1->GetMaximum() && h3->GetMaximum() > h2->GetMaximum() && h3->GetMaximum() > h4->GetMaximum() ) {
	h3->Draw();
	h1->Draw("same");
	h2->Draw("same");
	h4->Draw("same");
  }
  else {
	h4->Draw();
	h1->Draw("same");
	h2->Draw("same");
	h3->Draw("same");
  }

  //h3->SetFillColor( blue );
  //h10->SetFillColor( color3 );
  leg->SetFillColor(0);
  leg->AddEntry(h1, h1->GetName());
  leg->AddEntry(h2, h2->GetName());
  leg->AddEntry(h3, h3->GetName());
  leg->AddEntry(h4, h4->GetName());
  leg->Draw();
  //if( savename == "" ) 
  //if( savename.IsNull() ) 
  if( savename.Length() == 0 ) {
	//cout << "length is zero:" << savename << "::" << (TString)h1->GetName() << endl;
	c->SaveAs((TString)h1->GetName()+"_comp.png");
  }
  else {
	//cout << "length is not zero:" << savename << "::" << (TString)h1->GetName() << endl;
	c->SaveAs(savename+".png");
  }
  //c->SaveAs((TString)h1->GetName()+".eps");
  //need to reset these b'c they're pointers...
  h1->SetTitle(h1title);
  h2->SetTitle(h2title);
  h3->SetTitle(h3title);
  h4->SetTitle(h3title);

  gStyle->SetPadTopMargin   ( 0.08);//was 0.08
}
