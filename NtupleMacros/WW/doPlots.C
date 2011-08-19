// Make plots from histograms
#include "TFile.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "THStack.h"
#include <iostream>

TProfile* plotIsoPerformance( TFile* ftt, 
			      const char* signal,     // histogram name
			      const char* background, // histogram name
			      const char* name,       // unique name
			      bool reverse = false,   // normally signal near zero bin, reverse means signal is around max bin
			      double bkg_eff_min = 0,
			      double bkg_eff_max = 1,
			      double sig_eff_min = 0,
			      double sig_eff_max = 0
			      )
{ 
  TH1F* S = dynamic_cast<TH1F*>(ftt->Get(signal));
  if ( ! S ) {
    std::cout << "Error: histogram not found " << signal << std::endl;
    return 0;
  }
  TH1F* B = dynamic_cast<TH1F*>(ftt->Get(background));
  if ( ! B ) {
    std::cout << "Error: histogram not found " << background << std::endl;
    return 0;
  }
  char buf[1024];
  sprintf(buf,"c_%s",name);
  // TCanvas* c = new TCanvas(buf,buf,500,500);
  sprintf(buf,"p_%s",name);
  TProfile* p = new TProfile(buf,buf,50,bkg_eff_min,bkg_eff_max,sig_eff_min,sig_eff_max);
  p->SetLineColor(kBlue);
  p->SetLineWidth(2);
  p->SetMarkerStyle(20);
  p->SetMarkerSize(1);
  p->GetXaxis()->SetTitle("Background Efficiency");
  p->GetYaxis()->SetTitle("Signal Efficiency");
  p->SetStats(kFALSE);
  for( int i=0;i<=S->GetNbinsX()+1; ++i )
    if ( reverse ) 
      p->Fill(B->Integral(i,B->GetNbinsX()+1)/B->Integral(0,B->GetNbinsX()+1),
	      S->Integral(i,S->GetNbinsX()+1)/S->Integral(0,S->GetNbinsX()+1));
    else
      p->Fill(B->Integral(0,i)/B->Integral(0,B->GetNbinsX()+1),
	      S->Integral(0,i)/S->Integral(0,S->GetNbinsX()+1));
  // p->Draw();
  return p;
}

void doPlots()
{
  TFile *ftt = TFile::Open("processed_data_tag.root");
  assert(ftt);
  /*
  // trk isolation
  TProfile* p01 = plotIsoPerformance(ftt,"ww_helTrkIsoPassId","wjets_helTrkIsoFailId","TrkIso");
  p01->SetTitle("Electron tracker isolation: pat(blue) vs cms2(red)");
  p01->SetLineColor(kRed);
  p01->SetMarkerColor(kRed);
  TProfile* p02 = plotIsoPerformance(ftt,"ww_helTrkPatIsoPassId","wjets_helTrkPatIsoFailId","TrkPatIso");
  p02->SetTitle("Electron tracker isolation: pat(blue) vs cms2(red)");
  p02->SetLineColor(kBlue);
  p02->SetMarkerColor(kBlue);
  TCanvas* c01 = new TCanvas("c01","c01",500,500);
  p01->Draw();
  p02->Draw("same");

  // ecal isolation
  TProfile* p1 = plotIsoPerformance(ftt,"ww_helEcalJuraIsoPassId","wjets_helEcalJuraIsoFailId","EcalJuraIso");
  p1->SetTitle("Electron ECAL isolation: rechits(blue) vs basic clusters(red)");
  p1->SetLineColor(kRed);
  p1->SetMarkerColor(kRed);
  TProfile* p2 = plotIsoPerformance(ftt,"ww_helEcalPatIsoPassId","wjets_helEcalPatIsoFailId","EcalPatIso");
  p2->SetTitle("Electron ECAL isolation: rechits(blue) vs basic clusters(red)");
  p2->SetLineColor(kBlue);
  p2->SetMarkerColor(kBlue);
  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  p1->Draw();
  p2->Draw("same");
  
  // hcal isolation
  TProfile* p3 = plotIsoPerformance(ftt,"ww_helHcalConeIsoPassId","wjets_helHcalConeIsoFailId","HcalJuraIso");
  p3->SetTitle("Electron HCAL isolation: pat(blue) vs calotower(red)");
  p3->SetLineColor(kRed);
  p3->SetMarkerColor(kRed);
  TProfile* p4 = plotIsoPerformance(ftt,"ww_helHcalPatIsoPassId","wjets_helHcalPatIsoFailId","HcalPatIso");
  p4->SetTitle("Electron HCAL isolation: pat(blue) vs calotower(red)");
  p4->SetLineColor(kBlue);
  p4->SetMarkerColor(kBlue);
  TCanvas* c2 = new TCanvas("c2","c2",500,500);
  p3->Draw();
  p4->Draw("same");


  TCanvas* c3 = new TCanvas("c3","c3",500,500);
  TH1F* h3_1 = dynamic_cast<TH1F*>(ftt->Get("ww_helRelPatIsoPassId"));
  assert(h3_1);
  TH1F* h3_2 = dynamic_cast<TH1F*>(ftt->Get("wjets_helRelPatIsoFailId"));
  assert(h3_2);
  THStack* stack = new THStack("stack","Relative electron isolation scaled to unit area: signal(red), background(grey)");
  h3_1->Scale(1/(1e-8+h3_1->Integral()));
  h3_2->Scale(1/(1e-8+h3_2->Integral()));
  stack->Add(h3_2);
  stack->Add(h3_1);
  stack->Draw("hist");
  
  // relative isolation
  TProfile* p41 = plotIsoPerformance(ftt,"ww_helRelIsoPassId","wjets_helRelIsoFailId","RelIso",true);
  p41->SetTitle("Electron relative isolation: pat(blue) vs cms2(red)");
  p41->SetLineColor(kRed);
  p41->SetMarkerColor(kRed);
  TProfile* p42 = plotIsoPerformance(ftt,"ww_helRelPatIsoPassId","wjets_helRelPatIsoFailId","RelPatIso",true);
  p42->SetTitle("Electron relative isolation: pat(blue) vs cms2(red)");
  p42->SetLineColor(kBlue);
  p42->SetMarkerColor(kBlue);
  TCanvas* c4 = new TCanvas("c4","c4",500,500);
  p41->Draw();
  p42->Draw("same");
  */
  // Jet veto
  TProfile* p51 = plotIsoPerformance(ftt,"ww_hmaxCaloJetEt","ttbar_hmaxCaloJetEt","MaxCaloJetEt",false,0,0.02);
  p51->SetTitle("Most energetic jet Et: CaloJet(pink) vs GenJet(red) vs JPT(blue)");
  p51->SetLineColor(kMagenta);
  p51->SetMarkerColor(kMagenta);
  TProfile* p52 = plotIsoPerformance(ftt,"ww_hmaxGenJetEt","ttbar_hmaxGenJetEt","MaxGenJetEt",false,0,0.02);
  p52->SetTitle("Most energetic jet Et: CaloJet(pink) vs GenJet(red) vs JPT(blue)");
  p52->SetLineColor(kRed);
  p52->SetMarkerColor(kRed);
  TProfile* p53 = plotIsoPerformance(ftt,"ww_hmaxJPTEt","ttbar_hmaxJPTEt","MaxJPTEt",false,0,0.02);
  p53->SetTitle("Most energetic jet Et: CaloJet(pink) vs GenJet(red) vs JPT(blue)");
  p53->SetLineColor(kBlue);
  p53->SetMarkerColor(kBlue);
  TProfile* p54 = plotIsoPerformance(ftt,"ww_hmaxCaloTrkJetEt","ttbar_hmaxCaloTrkJetEt","MaxCaloTrkJetEt",false,0,0.02);
  p54->SetTitle("Most energetic jet Et: CaloJet(pink) vs GenJet(red) vs JPT(blue)");
  p54->SetLineColor(kYellow);
  p54->SetMarkerColor(kYellow);
  TCanvas* c5 = new TCanvas("c5","c5",500,500);
  p51->Draw();
  p54->Draw("same");
  p52->Draw("same");
  p53->Draw("same");

}

