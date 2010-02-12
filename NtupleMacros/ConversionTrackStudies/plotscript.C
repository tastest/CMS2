#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std; 

void plotscript(TString indata, TString inmc="", TString ind2="", TString inmc2="") {

  //plots i need to update:
  // tcmet, tcsumet, tcmet components, tcmet resolution
  // these 4 for 900 and 2360
  // this means figs 7, 8, 10.

  TH1::AddDirectory(false); 
  //gStyle->SetOptStat(0);
  //gStyle->SetPadRightMargin ( 0.1);

  TString f1a = "pt_twrs_ass";
  TString f1b = "mcpt_twrs_ass";
  TString f2a = "pt_twrshf_etmax_vs_clmet";
  TString f2b = "pt_twrshf_phivsclmetphi";
  TString f3a = "pt_cmet";
  TString f3b = "pt_cmetHFCorr";
  TString f4a = "pt_twrs_etmaxvser4";
  TString f4b = "mcpt_twrs_etmaxvser4";
  TString f5a = "pt_twrsec_etmax_vs_clmet";
  TString f5b = "pt_twrsec_phivsclmetphi";
  TString f6a = "pt_cmet";
  TString f6b = "pt_cmetHFCorr";
  TString f6c = "pt_cmetAllCorr";
  //TString f6b = "mcpt_cmet";
  //TString f6c = "";
  TString f7a = "pt_tcmetAllCorr";
  TString f7b = "mcpt_tcmet";
  TString f8a = "pt_tcmetxw";
  TString f8b = "pt_tcmetyw";
  TString f8c = "mcpt_tcmetxw";
  TString f8d = "mcpt_tcmetyw";
  TString f9a = "mcpt_cmet";  
  TString f10a= "pt_sumetAllCorr";
  TString f10b= "pt_tcsumetAllCorr";
  TString f10c= "mcpt_sumet";
  TString f10d= "mcpt_tcsumet";
  
  TFile* dafile_ = new TFile(indata.Data(), "READ");
  TFile* mcfile_ = new TFile(inmc.Data(), "READ");
  TFile* d2file_ = new TFile(ind2.Data(), "READ");
  TFile* m2file_ = new TFile(inmc2.Data(), "READ");

  string cms = "CMS Preliminary";
  string d9 = "#sqrt{s} = 900 GeV";
  string d2 = "#sqrt{s} = 2360 GeV";
  TLatex latex;
  latex.SetTextAlign(12);
  latex.SetTextSize(0.05);
  latex.SetTextFont(62);
  latex.SetNDC();


  TCanvas *c1 = new TCanvas();
  //cout << "getting " << f1a.Data() << endl;
  TH1F* h1a = (TH1F*)(dafile_->Get(f1a.Data())->Clone());
  TH1F* h1b = (TH1F*)(mcfile_->Get(f1b.Data())->Clone());
  h1a->GetXaxis()->SetTitleSize(0.06);
  h1a->SetMarkerStyle(20);
  h1a->SetMarkerSize(0.7);
  h1b->SetLineColor(2);
  h1b->Scale(184245./1988370.);
  gPad->SetLogy();
  // binning: 201, -1.0, 1.01);
  cout << "Figure 1: alpha: " << endl;
  cout << "data: total: " << h1a->Integral() << "  nspike: " << h1a->Integral(0,20)+h1a->Integral(200,202)
	   << " ratio: " << (h1a->Integral(0,20)+h1a->Integral(200,202))/h1a->Integral() << endl;
  cout << "mc: total: " << h1b->Integral() << "  nspike: " << h1b->Integral(0,20)+h1b->Integral(200,202)
	   << " ratio: " << (h1b->Integral(0,20)+h1b->Integral(200,202))/h1b->Integral() << endl << endl;
  h1a->Draw();
  h1b->Draw("same,hist");
  c1->SaveAs("pas_note/plots/alpha.eps");

  //fig 2
  TH2F* h2a = (TH2F*)(dafile_->Get(f2a.Data())->Clone());
  h2a->GetXaxis()->SetTitleSize(0.055);
  h2a->GetXaxis()->SetTitleOffset(0.95);
  h2a->GetYaxis()->SetTitle("caloMET");
  TH2F* h2b = (TH2F*)(dafile_->Get(f2b.Data())->Clone());
  h2b->GetXaxis()->SetTitleSize(0.06);
  //h2b->GetYaxis()->SetTitle("caloMET #phi"); //fixed in looper
  cout << "Fig 2: calo met vs twr et: " << endl;
  const int nbins = 30;
  float nevtmet29 = 0;
  float nevtmet29twr15 = 0;
  for( int i=0; i<=nbins+1; i++ ) {
	nevtmet29 += h2a->GetBinContent(i, nbins);
	if( i < 15 ) nevtmet29twr15 += h2a->GetBinContent(i, nbins);
  }
  cout << "num tot evts above met 29: " << nevtmet29 << endl;
  cout << "num tot evts above met 29, twr et < 15: " << nevtmet29twr15 << endl << endl;
  TCanvas *c2a = new TCanvas(); //2
  h2a->Draw("colz");
  c2a->SaveAs("pas_note/plots/cmet_v_hftwret.eps");
  TCanvas *c2b = new TCanvas(); //3
  h2b->Draw("colz");
  c2b->SaveAs("pas_note/plots/cmetphi_v_hftwrphi.eps");

  //fig 3
  TH1F* h3a = (TH1F*)(dafile_->Get(f3a.Data())->Clone());
  TH1F* h3b = (TH1F*)(dafile_->Get(f3b.Data())->Clone());
  cout << "Fig 3: num tot evts above met 29: " << h3a->Integral(30,41) << endl << endl;
  TCanvas *c3 = new TCanvas(); //4
  gPad->SetLogy();
  h3a->GetXaxis()->SetTitleSize(0.06);
  h3a->Draw("hist");
  h3b->SetLineColor(2);
  h3b->Draw("same");
  c3->SaveAs("pas_note/plots/clmet_hfclean.eps");

  //fig 4
  TH1F* h4a = (TH1F*)(dafile_->Get(f4a.Data())->Clone());
  h4a->GetYaxis()->SetTitleOffset(0.7);
  h4a->GetXaxis()->SetTitleSize(0.06);
  TH1F* h4b = (TH1F*)(mcfile_->Get(f4b.Data())->Clone());
  h4b->GetYaxis()->SetTitleOffset(0.7);
  h4b->GetXaxis()->SetTitleSize(0.06);
  TCanvas *c4a = new TCanvas(); //5
  gPad->SetLogz();
  h4a->Draw("colz");
  c4a->SaveAs("pas_note/plots/data_r4_v_twret.eps");
  TCanvas *c4b = new TCanvas(); //6
  gPad->SetLogz();
  h4b->Draw("colz");
  c4b->SaveAs("pas_note/plots/mc_r4_v_twret.eps");

  TH1F* h5a = (TH1F*)(dafile_->Get(f5a.Data())->Clone());
  //h5a->GetYaxis()->SetTitleOffset(0.8);
  h5a->GetXaxis()->SetTitleSize(0.06);
  TH1F* h5b = (TH1F*)(dafile_->Get(f5b.Data())->Clone());
  TCanvas *c5a = new TCanvas(); //7
  h5a->Draw("colz");
  TCanvas *c5b = new TCanvas(); //8
  h5b->GetXaxis()->SetTitleSize(0.06);
  h5b->Draw("colz");
  c5a->SaveAs("pas_note/plots/cmet_v_twret.eps");
  c5b->SaveAs("pas_note/plots/cmetphi_v_twrphi.eps");

  //fig 6
  TH1F* h6a = (TH1F*)(dafile_->Get(f6a.Data())->Clone());
  TH1F* h6b = (TH1F*)(dafile_->Get(f6b.Data())->Clone());
  TH1F* h6c = (TH1F*)(dafile_->Get(f6c.Data())->Clone());
  //TH1F* h6b = (TH1F*)(mcfile_->Get(f6b.Data())->Clone());
  TCanvas *c6a = new TCanvas(); //9
  //nmc evts=1988370, ndata evts=184245
  gPad->SetLogy();
  h6a->GetXaxis()->SetTitleSize(0.06);
  //h6b->Scale(184245./1988370.);
  h6b->SetLineColor(2); //red
  h6b->SetMarkerStyle(22);
  h6b->SetMarkerSize(1.2);
  h6b->SetMarkerColor(2);
  h6c->SetLineColor(4); //blue
  for(int i=0; i<=40; i++) h6b->SetBinError(i,0.001);
  h6a->Draw("hist");
  h6b->Draw("same");
  h6c->Draw("same");
  //c6a->SaveAs("pas_note/plots/calomet_datamc.eps");
  c6a->SaveAs("pas_note/plots/calomet_all.eps");

  //fig 7
  TH1F* h7a = (TH1F*)(dafile_->Get(f7a.Data())->Clone());
  TH1F* h7b = (TH1F*)(mcfile_->Get(f7b.Data())->Clone());
  TCanvas *c7a = new TCanvas(); //10
  //nmc evts=1988370, ndata evts=184245. ratio = 0.092661326
  //NEW: mc =1767180, data      =164703. ratio = 0.093201032
  gPad->SetLogy();
  h7a->SetMarkerStyle(20);
  h7a->SetMarkerSize(0.7);
  h7a->Draw();
  //h7b->Scale(184245./1988370.);
  h7b->Scale(164703./1767180.);
  h7b->SetLineColor(2);
  h7b->SetFillColor(2);
  h7b->GetXaxis()->SetTitle("tcMET [GeV]");
  h7b->GetXaxis()->SetTitleSize(0.06);
  h7b->GetXaxis()->SetTitleOffset(0.9);
  h7b->GetXaxis()->SetLabelOffset(0.005);
  h7b->GetYaxis()->SetTitle("Events/GeV");
  h7b->GetYaxis()->SetTitleOffset(1.1);
  h7b->GetYaxis()->SetLabelOffset(0.005);
  h7b->SetTickLength(0.03,"XYZ");
  h7b->Draw("hist");
  h7a->Draw("same");
  latex.DrawLatex(0.4, 0.8, (cms).c_str());
  latex.DrawLatex(0.4, 0.74, (d9).c_str());
  c7a->SaveAs("pas_note/plots/tcmet_datamc.eps");

  //fig 8
  TH1F* h8a = (TH1F*)(dafile_->Get(f8a.Data())->Clone()); //data
  TH1F* h8b = (TH1F*)(dafile_->Get(f8b.Data())->Clone());
  h8a->Add( h8b );
  h8a->SetMarkerStyle(20);
  h8a->SetMarkerSize(0.7);
  TH1F* h8c = (TH1F*)(mcfile_->Get(f8c.Data())->Clone()); //mc
  TH1F* h8d = (TH1F*)(mcfile_->Get(f8d.Data())->Clone());
  h8c->Add( h8d );
  //h8c->Scale(184245./1988370.);
  h8c->Scale(164703./1767180.);
  h8c->GetXaxis()->SetTitle("tcMET components [GeV]");
  h8c->GetXaxis()->SetTitleSize(0.06);
  h8c->GetXaxis()->SetTitleOffset(0.9);
  h8c->GetXaxis()->SetLabelOffset(0.005);
  h8c->GetYaxis()->SetTitle("Events/GeV");
  h8c->GetYaxis()->SetTitleOffset(1.1);
  h8c->GetYaxis()->SetLabelOffset(0.005);
  h8c->SetLineColor(2);
  h8c->SetFillColor(2);
  h8c->SetTickLength(0.03,"XYZ");
  TCanvas *c8a = new TCanvas(); //11
  h8c->Draw("hist");
  h8a->Draw("same");
  //gStyle->SetTickLength(-0.03,"XYZ");
  latex.DrawLatex(0.54, 0.9, (cms).c_str());
  latex.DrawLatex(0.6, 0.84, (d9).c_str());
  c8a->SaveAs("pas_note/plots/tcmet_comps_lin.eps");
  gPad->SetLogy();
  c8a->SaveAs("pas_note/plots/tcmet_comps.eps");

  //7 and 8 are repeated as 10 and 11 at 2360gev

  //fig 10 (fig 9 is actually from the other script)
  TH1F* h9a = (TH1F*)(d2file_->Get(f7a.Data())->Clone());
  TH1F* h9b = (TH1F*)(m2file_->Get(f7b.Data())->Clone());
  TCanvas *c9a = new TCanvas(); //12
  h9b->SetLineColor(2);
  h9b->SetFillColor(2);
  h9b->GetXaxis()->SetTitle("tcMET [GeV]");
  h9b->GetXaxis()->SetTitleSize(0.06);
  h9b->GetXaxis()->SetTitleOffset(0.9);
  h9b->GetXaxis()->SetLabelOffset(0.005);
  h9b->GetYaxis()->SetTitle("Events/GeV");
  h9b->GetYaxis()->SetTitleOffset(1.1);
  h9b->GetYaxis()->SetLabelOffset(0.005);
  h9b->SetTickLength(0.03,"XYZ");
  //num 2tev data evts = 10541, 2tev mc=402464. ratio = 0.026191162
  //NEW:          data = 9685,       mc=362673. ratio = 0.026704497
  //h9b->Scale(10541./402464.);
  h9b->Scale(9685./362673.);
  h9a->SetMarkerStyle(20);
  h9a->SetMarkerSize(0.7);
  gPad->SetLogy();
  h9b->Draw("hist");
  h9a->Draw("same"); //draw data on top of mc
  latex.DrawLatex(0.4, 0.8, (cms).c_str());
  latex.DrawLatex(0.4, 0.74, (d2).c_str());
  c9a->SaveAs("pas_note/plots/tcmet_datamc_2tev.eps");

  //fig 11
  TH1F* h10a = (TH1F*)(d2file_->Get(f8a.Data())->Clone());
  TH1F* h10b = (TH1F*)(d2file_->Get(f8b.Data())->Clone());
  h10a->Add( h10b );
  h10a->SetMarkerStyle(20);
  h10a->SetMarkerSize(0.7);
  //h10a->GetXaxis()->SetTitleOffset(1.0);
  TH1F* h10c = (TH1F*)(m2file_->Get(f8c.Data())->Clone()); //mc
  TH1F* h10d = (TH1F*)(m2file_->Get(f8d.Data())->Clone());
  h10c->Add( h10d );
  h10c->GetXaxis()->SetTitle("tcMET components [GeV]");
  h10c->GetXaxis()->SetTitleSize(0.06);
  h10c->GetXaxis()->SetTitleOffset(0.9);
  h10c->GetXaxis()->SetLabelOffset(0.005);
  h10c->GetYaxis()->SetTitle("Events/GeV");
  h10c->GetYaxis()->SetTitleOffset(1.1);
  h10c->GetYaxis()->SetLabelOffset(0.005);
  h10c->SetTickLength(0.03,"XYZ");
  //h10c->Scale(h10a->Integral()/h10c->Integral()); //ok to normalize by integral?
  //num 2tev data evts = 10541, 2tev mc=402464
  //h10c->Scale(10541./402464.);
  h10c->Scale(9685./362673.);
  h10c->SetLineColor(2);
  h10c->SetFillColor(2);
  TCanvas *c10a = new TCanvas(); //13
  h10c->Draw("hist");
  h10a->Draw("same");
  latex.DrawLatex(0.54, 0.9, (cms).c_str());
  latex.DrawLatex(0.6, 0.84, (d2).c_str());
  c10a->SaveAs("pas_note/plots/tcmet_comps_2tev_lin.eps");
  gPad->SetLogy();
  c10a->SaveAs("pas_note/plots/tcmet_comps_2tev.eps");

  //new plot of cmet data/mc
  TH1F* h11a = (TH1F*)(dafile_->Get(f3a.Data())->Clone()); //data
  TH1F* h11b = (TH1F*)(mcfile_->Get(f9a.Data())->Clone()); //mc
  TCanvas *c11 = new TCanvas(); //14
  gPad->SetLogy();
  h11a->SetMarkerStyle(20);
  h11a->SetMarkerSize(0.6);
  h11b->GetXaxis()->SetTitleSize(0.06);
  h11b->Scale(184245./1988370.);
  h11b->SetLineColor(2);
  h11b->Draw("hist");
  h11a->Draw("same");
  c11->SaveAs("pas_note/plots/clmet_mc_data.eps");

  //new plot 2: cmet allcorr with mc f9a and f6c
  TH1F* h12a = (TH1F*)(dafile_->Get(f6c.Data())->Clone()); //data
  TH1F* h12b = (TH1F*)(mcfile_->Get(f9a.Data())->Clone()); //mc
  TCanvas *c12 = new TCanvas(); //15
  gPad->SetLogy();
  h12a->SetMarkerStyle(20);
  h12a->SetMarkerSize(0.6);
  h12b->GetXaxis()->SetTitleSize(0.06);
  h12b->Scale(184245./1988370.);
  h12b->SetLineColor(2);
  h12b->Draw("hist");
  h12a->Draw("same");
  c12->SaveAs("pas_note/plots/clmet_mc_datacorr.eps");

  //new plot 3.1: sumet data corr mc 900 -- f10a, f10c
  TH1F* h13a = (TH1F*)(dafile_->Get(f10a.Data())->Clone()); //data
  TH1F* h13b = (TH1F*)(mcfile_->Get(f10c.Data())->Clone()); //mc
  TCanvas *c13 = new TCanvas(); //16
  gPad->SetLogy();
  h13a->SetMarkerStyle(20);
  h13a->SetMarkerSize(0.6);
  h13b->GetXaxis()->SetTitleSize(0.06);
  //h13b->Scale(184245./1988370.);
  h13b->Scale(164703./1767180.);
  h13b->SetLineColor(2);
  h13b->SetTitle("900 Gev Selected Runs");
  h13b->Draw("hist");
  h13a->Draw("same");
  c13->SaveAs("pas_note/plots/sumet_mc_datacorr.eps");

  //new plot 3.2: sumet data mc 2360 -- f10a, f10c
  TH1F* h14a = (TH1F*)(d2file_->Get(f10a.Data())->Clone()); //data
  TH1F* h14b = (TH1F*)(m2file_->Get(f10c.Data())->Clone()); //mc
  TCanvas *c14 = new TCanvas(); //17
  gPad->SetLogy();
  h14a->SetMarkerStyle(20);
  h14a->SetMarkerSize(0.6);
  h14b->GetXaxis()->SetTitleSize(0.06);
  //h14b->Scale(10541./402464.);
  h14b->Scale(9685./362673.);
  h14b->SetLineColor(2);
  //h14b->SetTitle("2360 Gev (Run 124120)");
  h14b->Draw("hist");
  h14b->GetYaxis()->SetRangeUser(0.04,1500);
  h14a->Draw("same");
  c14->SaveAs("pas_note/plots/sumet_2tev_mc_datacorr.eps");

  //new plot 3.3: tcsumet data mc 900 -- f10b, f10d
  TH1F* h15a = (TH1F*)(dafile_->Get(f10b.Data())->Clone()); //data
  TH1F* h15b = (TH1F*)(mcfile_->Get(f10d.Data())->Clone()); //mc
  TCanvas *c15 = new TCanvas(); //18
  gPad->SetLogy();
  h15a->SetMarkerStyle(20);
  h15a->SetMarkerSize(0.7);
  h15b->GetXaxis()->SetTitleSize(0.06);
  //h15b->Scale(184245./1988370.);
  h15b->Scale(164703./1767180.);
  h15b->SetLineColor(2);
  h15b->SetFillColor(2);
  h15b->GetXaxis()->SetTitle("tcSumET [GeV]");
  h15b->GetXaxis()->SetTitleSize(0.06);
  h15b->GetXaxis()->SetTitleOffset(0.9);
  h15b->GetXaxis()->SetLabelOffset(0.005);
  h15b->GetYaxis()->SetTitle("Events/GeV");
  h15b->GetYaxis()->SetTitleOffset(1.1);
  h15b->GetYaxis()->SetLabelOffset(0.005);
  h15b->SetTickLength(0.03,"XYZ");
  //h15b->SetTitle("900 Gev Selected Runs");
  h15b->Draw("hist");
  h15b->GetYaxis()->SetRangeUser(0.1,20000);
  h15a->Draw("same");
  latex.DrawLatex(0.4, 0.9, (cms).c_str());
  latex.DrawLatex(0.4, 0.84, (d9).c_str());
  c15->SaveAs("pas_note/plots/tcsumet_mc_datacorr.eps");

  //new plot 3.4: tcsumet data mc 2360 -- f10b, f10d
  TH1F* h16a = (TH1F*)(d2file_->Get(f10b.Data())->Clone()); //data
  TH1F* h16b = (TH1F*)(m2file_->Get(f10d.Data())->Clone()); //mc
  TCanvas *c16 = new TCanvas(); //19
  gPad->SetLogy();
  h16a->SetMarkerStyle(20);
  h16a->SetMarkerSize(0.7);
  h16b->GetXaxis()->SetTitleSize(0.06);
  //h16b->Scale(10541./402464.);
  h16b->Scale(9685./362673.);
  h16b->SetLineColor(2);
  h16b->SetFillColor(2);
  h16b->GetXaxis()->SetTitle("tcSumET [GeV]");
  h16b->GetXaxis()->SetTitleSize(0.06);
  h16b->GetXaxis()->SetTitleOffset(0.9);
  h16b->GetXaxis()->SetLabelOffset(0.005);
  h16b->GetYaxis()->SetTitle("Events/GeV");
  h16b->GetYaxis()->SetTitleOffset(1.1);
  h16b->GetYaxis()->SetLabelOffset(0.005);
  h16b->SetTickLength(0.03,"XYZ");
  //h16b->SetTitle("2360 Gev (Run 124120)");
  h16b->Draw("hist");
  h16b->GetYaxis()->SetRangeUser(0.1,800);
  h16a->Draw("same");
  latex.DrawLatex(0.4, 0.9, (cms).c_str());
  latex.DrawLatex(0.4, 0.84, (d2).c_str());
  c16->SaveAs("pas_note/plots/tcsumet_2tev_mc_datacorr.eps");



}
