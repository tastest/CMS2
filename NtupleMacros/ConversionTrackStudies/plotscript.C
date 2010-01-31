#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std; 

void plotscript(TString indata, TString inmc="", TString ind2="", TString inmc2="") {

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
  TString f8a = "pt_tcmetx";
  TString f8b = "pt_tcmety";
  TString f8c = "mcpt_tcmetx";
  TString f8d = "mcpt_tcmety";
  
  TFile* dafile_ = new TFile(indata.Data(), "READ");
  TFile* mcfile_ = new TFile(inmc.Data(), "READ");
  TFile* d2file_ = new TFile(ind2.Data(), "READ");
  //TFile* m2file_ = new TFile(inmc2.Data(), "READ");

  TCanvas *c1 = new TCanvas();
  //cout << "getting " << f1a.Data() << endl;
  TH1F* h1a = (TH1F*)(dafile_->Get(f1a.Data())->Clone());
  TH1F* h1b = (TH1F*)(mcfile_->Get(f1b.Data())->Clone());
  h1a->GetXaxis()->SetTitleSize(0.06);
  h1a->SetMarkerStyle(20);
  h1a->SetMarkerSize(0.6);
  h1b->SetLineColor(2);
  h1b->Scale(185939./1988370.);
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
  h2b->GetYaxis()->SetTitle("caloMET #phi");
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
  //nmc evts=1988370, ndata evts=185939
  gPad->SetLogy();
  h6a->GetXaxis()->SetTitleSize(0.06);
  //h6b->Scale(185939./1988370.);
  h6b->SetLineColor(2); //red
  h6b->SetMarkerStyle(22);
  h6b->SetMarkerSize(1.2);
  h6b->SetMarkerColor(2);
  h6c->SetLineColor(4); //blue
  for(int i=0; i<40; i++) h6b->SetBinError(i,0.001);
  h6a->Draw("hist");
  h6b->Draw("same");
  h6c->Draw("same");
  //c6a->SaveAs("pas_note/plots/calomet_datamc.eps");
  c6a->SaveAs("pas_note/plots/calomet_all.eps");

  TH1F* h7a = (TH1F*)(dafile_->Get(f7a.Data())->Clone());
  TH1F* h7b = (TH1F*)(mcfile_->Get(f7b.Data())->Clone());
  TCanvas *c7a = new TCanvas(); //10
  //nmc evts=1988370, ndata evts=185939
  gPad->SetLogy();
  h7a->GetXaxis()->SetTitleSize(0.06);
  h7a->Draw();
  h7b->Scale(185939./1988370.);
  h7b->SetLineColor(2);
  h7b->Draw("same,hist");
  c7a->SaveAs("pas_note/plots/tcmet_datamc.eps");

  TH1F* h8a = (TH1F*)(dafile_->Get(f8a.Data())->Clone()); //data
  TH1F* h8b = (TH1F*)(dafile_->Get(f8b.Data())->Clone());
  h8a->Add( h8b );
  h8a->GetXaxis()->SetTitle("tcMET components");
  h8a->GetXaxis()->SetTitleSize(0.06);
  h8a->SetMarkerStyle(20);
  h8a->SetMarkerSize(0.6);
  TH1F* h8c = (TH1F*)(mcfile_->Get(f8c.Data())->Clone()); //mc
  TH1F* h8d = (TH1F*)(mcfile_->Get(f8d.Data())->Clone());
  h8c->Add( h8d );
  h8c->Scale(185939./1988370.);
  TCanvas *c8a = new TCanvas(); //11
  h8a->Draw("");
  h8c->SetLineColor(2);
  h8c->Draw("same,hist");
  c8a->SaveAs("pas_note/plots/tcmet_comps.eps");

  //7 and 8 are repeated as 9 and 10 at 2360gev

  //fig 10 (fig 9 is actually from the other script)
  TH1F* h9a = (TH1F*)(d2file_->Get(f7a.Data())->Clone());
  //TH1F* h9b = (TH1F*)(m2file_->Get(f7b.Data())->Clone());
  TCanvas *c9a = new TCanvas(); //12
  h9a->GetXaxis()->SetTitleSize(0.06);
  gPad->SetLogy();
  h9a->Draw();
  //h9b->SetLineColor(2);
  //h9b->Draw("same,hist");
  c9a->SaveAs("pas_note/plots/tcmet_datamc_2tev.eps");

  //fig 11
  TH1F* h10a = (TH1F*)(d2file_->Get(f8a.Data())->Clone());
  TH1F* h10b = (TH1F*)(d2file_->Get(f8b.Data())->Clone());
  h10a->Add( h10b );
  h10a->GetXaxis()->SetTitle("tcMET components");
  h10a->GetXaxis()->SetTitleSize(0.06);
  h10a->SetMarkerStyle(20);
  h10a->SetMarkerSize(0.6);
  //h10a->GetXaxis()->SetTitleOffset(1.0);
  //TH1F* h10c = (TH1F*)(mcfile_->Get(f8c.Data())->Clone()); //mc
  //TH1F* h10d = (TH1F*)(mcfile_->Get(f8d.Data())->Clone());
  //h10c->Add( h10d );
  //h10c->Scale(h10a->Integral()/h10c->Integral()); //ok to normalize by integral?
  TCanvas *c10a = new TCanvas(); //13
  //h10c->SetLineColor(2);
  //h10c->Draw("hist");
  //h10a->Draw("same");
  h10a->Draw();
  c10a->SaveAs("pas_note/plots/tcmet_comps_2tev.eps");
}
