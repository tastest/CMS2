#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "smurfAnalysis.h"
#include <fstream>
#include "TRandom3.h"
#include "TCanvas.h"

void rOutIn(const double met=20)
{
  const unsigned int cut = 
    SmurfTree::Lep1FullSelection|
    SmurfTree::Lep2FullSelection|
    SmurfTree::BaseLine|
    SmurfTree::ChargeMatch|
    SmurfTree::TopVeto|
    SmurfTree::ExtraLeptonVeto;
  
  SmurfTree tree;
  // tree.LoadTree("/smurf/data/Run2011_Spring11_SmurfV6_42X/tas-TightLooseFullMET-alljets/data-met20-715ipb.root");
  tree.LoadTree("/smurf/data/Run2011_Spring11_SmurfV6_42X/tas-TightLooseFullMET-alljets/data-met20.root");
  // tree.LoadTree("/smurf/data/Run2011_Spring11_SmurfV6/mitf-alljets/dymm.root");
  tree.InitTree();
  Long64_t nDataEntries = tree.tree_->GetEntries();
  Double_t bins[6]={17,20,23,26,29,35};
  TH1F* hin_mm = new TH1F("hin_mm","hin_mm",5,bins);
  TH1F* hin_ee = new TH1F("hin_ee","hin_ee",5,bins);
  TH1F* hin_em = new TH1F("hin_em","hin_em",5,bins);
  hin_mm->SetDirectory(0);
  hin_mm->Sumw2();
  hin_ee->SetDirectory(0);
  hin_ee->Sumw2();
  hin_em->SetDirectory(0);
  hin_em->Sumw2();
  TH1F* hout_mm = new TH1F("hout_mm","hout_mm",5,bins);
  TH1F* hout_ee = new TH1F("hout_ee","hout_ee",5,bins);
  TH1F* hout_em = new TH1F("hout_em","hout_em",5,bins);
  hout_mm->SetDirectory(0);
  hout_mm->Sumw2();
  hout_ee->SetDirectory(0);
  hout_ee->Sumw2();
  hout_em->SetDirectory(0);
  hout_em->Sumw2();
  unsigned nin_mm(0);
  unsigned nin_em(0);
  unsigned nin_ee(0);
  unsigned nout_mm(0);
  unsigned nout_em(0);
  unsigned nout_ee(0);

  for (Long64_t i = 0; i < nDataEntries; i++){
    tree.tree_->GetEntry(i);
    const double mll = 75.9;
    // if ( ! (tree.lep1_.pt()>40 && tree.lep2_.pt()>25 && fabs(tree.dPhi_)<M_PI*100/180) ) continue;
    // if ( ! (tree.lep1_.pt()>25 && tree.lep2_.pt()>15 && fabs(tree.dPhi_)<M_PI*90/180) ) continue;
    // if ( ! (tree.lep1_.pt()>40 && tree.lep2_.pt()>25 && fabs(tree.dPhi_)<M_PI*100/180) ) continue;
    
    // if ( ! (tree.mt_>80 && tree.mt_<130) ) continue;
    
    // if ( !(tree.jet1_.pt()<15||tree.dPhiDiLepJet1_<M_PI/180*165) ) continue;
    // if ( tree.jet1_.pt()>25 ) continue;
    double weight = tree.scale1fb_;

    if ( !(tree.njets_==0 && (tree.cuts_&cut)==cut)  ) continue;
    // if ( !((tree.cuts_&cut)==cut)  ) continue;

    if (tree.type_==0 && fabs(tree.dilep_.mass()-91)<15 ) 
      hin_mm->Fill(std::min(tree.pmet_,tree.pTrackMet_), weight);
    if (tree.type_==3 && fabs(tree.dilep_.mass()-91)<15 ) 
      hin_ee->Fill(std::min(tree.pmet_,tree.pTrackMet_), weight);
    if ((tree.type_==1 || tree.type_==2) && fabs(tree.dilep_.mass()-91)<15 ) 
      hin_em->Fill(std::min(tree.pmet_,tree.pTrackMet_), weight);

    if (mll<76) {
      if (tree.type_==0 && tree.dilep_.mass()<mll ) 
	hout_mm->Fill(std::min(tree.pmet_,tree.pTrackMet_), weight);
      if (tree.type_==3 && tree.dilep_.mass()<mll ) 
	hout_ee->Fill(std::min(tree.pmet_,tree.pTrackMet_), weight);
      if ((tree.type_==1 || tree.type_==2) && tree.dilep_.mass()<mll ) 
	hout_em->Fill(std::min(tree.pmet_,tree.pTrackMet_), weight);
    } else {
      if (tree.type_==0 && fabs(tree.dilep_.mass()-91)>15 ) 
	hout_mm->Fill(std::min(tree.pmet_,tree.pTrackMet_), weight);
      if (tree.type_==3 && fabs(tree.dilep_.mass()-91)>15 ) 
	hout_ee->Fill(std::min(tree.pmet_,tree.pTrackMet_), weight);
      if ((tree.type_==1 || tree.type_==2) && fabs(tree.dilep_.mass()-91)>15 ) 
	hout_em->Fill(std::min(tree.pmet_,tree.pTrackMet_), weight);
    }      
    
    if ( std::min(tree.pmet_,tree.pTrackMet_)>met ){
      // if ( tree.dilep_.pt()>met ){
      if (tree.type_==0 && fabs(tree.dilep_.mass()-91)<15 ) nin_mm++;
      if (tree.type_==3 && fabs(tree.dilep_.mass()-91)<15 ) nin_ee++;
      if ((tree.type_==1 || tree.type_==2) && fabs(tree.dilep_.mass()-91)<15 ) nin_em++;
      
      if (mll<76) {
	if (tree.type_==0 && tree.dilep_.mass()<mll ) nout_mm++;
	if (tree.type_==3 && tree.dilep_.mass()<mll ) nout_ee++;
	if ((tree.type_==1 || tree.type_==2) && tree.dilep_.mass()<mll ) nout_em++;
      } else {
	if (tree.type_==0 && fabs(tree.dilep_.mass()-91)>15 ) nout_mm++;
	if (tree.type_==3 && fabs(tree.dilep_.mass()-91)>15 ) nout_ee++;
	if ((tree.type_==1 || tree.type_==2) && fabs(tree.dilep_.mass()-91)>15 ) nout_em++;
      }      
    }

  }
  TCanvas* c1 = new TCanvas("c1","c1",600,900);
  c1->Divide(2,3);
  c1->cd(1);
  hin_mm->Draw();
  c1->cd(2);
  hout_mm->Draw();
  c1->cd(3);
  hin_ee->Draw();
  c1->cd(4);
  hout_ee->Draw();
  c1->cd(5);
  hin_em->Draw();
  c1->cd(6);
  hout_em->Draw();
  new TCanvas("c2","c2",500,500);
  TH1F* hout = (TH1F*)hout_mm->Clone("hout");
  hout->Add(hout_ee,pow(1/0.77,2));
  hout->Add(hout_em,-1/0.77);
  TH1F* hin = (TH1F*)hin_mm->Clone("hin");
  hin->Add(hin_ee,pow(1/0.8,2));
  hin->Add(hin_em,-1/0.8);
  TH1F* hr = (TH1F*)hout->Clone("hr");
  hr->SetTitle("1/fb of data using opposite flavor subtraction");
  hr->GetYaxis()->SetTitle("rOutIn");
  hr->GetXaxis()->SetTitle("minMet");
  hr->Divide(hout,hin);
  hr->Draw();

  double Nout = (nout_mm+nout_ee/pow(0.77,2)-nout_em/0.77);
  double Nin  = (nin_mm+nin_ee/pow(0.77,2)-nin_em/0.77);
  std::cout << "R for special study: " << Nout/Nin << "\tNout/Nin: " << Nout << " / " << Nin << std::endl;

}
