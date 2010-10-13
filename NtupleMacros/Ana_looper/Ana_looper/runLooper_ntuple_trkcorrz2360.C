
//TCut baseline_minbias_trk   = "TMath::Abs(ran_trksp4.eta())<2.4 && TMath::Abs(ran_trksp4.eta())>1.0" ;
//TCut baseline_minbias_trk   = "TMath::Abs(ran_trksp4.eta())<1.0" ;
TCut baseline_minbias_EB   = "TMath::Abs(ran_trksp4.eta())<1.0" ;
TCut baseline_minbias_EE   = "TMath::Abs(ran_trksp4.eta())<2.4  && TMath::Abs(ran_trksp4.eta())>1.0" ;


//TCut baseline_minbias_trk   = "TMath::Abs(ran_trksp4.eta())<0.6 && TMath::Abs(ran_trksp4.eta())>0.4 && TMath::Abs(ran_trksp4.phi())<1.1 && TMath::Abs(ran_trksp4.phi())>0.9" ;
void setStats(TH1* s,TH1* r, double startingY, double startingX = .1,bool fit){
  if (startingY<0){
    s->SetStats(0);
    r->SetStats(0);
  } else {
    //gStyle->SetOptStat(1001);                                                                                                                                                          
    s->SetStats(1);
    r->SetStats(1);

    if (fit){
      s->Fit("gaus");

      TF1* f1 = (TF1*) s->GetListOfFunctions()->FindObject("gaus");
      f1->SetLineColor(2);
      f1->SetLineWidth(1);
    }
    s->Draw();
    gPad->Update();
    TPaveStats* st1 = (TPaveStats*) s->GetListOfFunctions()->FindObject("stats");
    if (fit) {st1->SetOptFit(0010);    st1->SetOptStat(1001);}
    st1->SetX1NDC(startingX);
    st1->SetX2NDC(startingX+0.30);
    st1->SetY1NDC(startingY+0.20);
    st1->SetY2NDC(startingY+0.35);
    st1->SetTextColor(s->GetLineColor());
    if (fit) {
      r->Fit("gaus");
      TF1* f2 = (TF1*) r->GetListOfFunctions()->FindObject("gaus");
      f2->SetLineColor(4);
      f2->SetLineWidth(1);
    }
    r->Draw();
    gPad->Update();
    TPaveStats* st2 = (TPaveStats*) r->GetListOfFunctions()->FindObject("stats");
    if (fit) {st2->SetOptFit(0010);    st2->SetOptStat(1001);}
    st2->SetX1NDC(startingX);
    st2->SetX2NDC(startingX+0.30);
    st2->SetY1NDC(startingY);
    st2->SetY2NDC(startingY+0.15);
    st2->SetTextColor(r->GetLineColor());
  }
}


void setBin_int_above(TH1* hist){
  //hist->Sumw2();

  int   nBins = hist->GetNbinsX();
  TH1* hist_copy = hist->Clone();
  // float nIntegral = hist->Integral(0, nBins+1);
  // hist->Scale(1/nIntegral);
  double eff = 0;
  double err = 0;
  for (int i = 1; i < = nBins; i++){
    //  hist->SetBinContent(i, hist->Integral(i,nBins+1 ));
    eff =  hist_copy->Integral(0,i );
    err = sqrt(eff*(1-eff)/hist_copy->GetEntries());
    hist->SetBinContent(i, eff);
    hist->SetBinError(i, err);
  }

}


runLooper_Iso_MinBias_skimdata_mc_EE( TString iso_opt = "ran_ecalIso03_mu",  bool jura = false){
 

  gStyle->SetHistLineWidth(2);
  gROOT->ForceStyle();
  gStyle->SetOptStat(1101);

  //gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("can", "can");
 
  TChain *MinBias_skimdata = new TChain("Events");
 
  TChain *MinBias_mc = new TChain("Events");
  MinBias_mc->Add("MCCorrZ2360_skimmednTuple.root"); 
  //MinBias_mc->Add("/data/tmp/yanjuntu/cms2/cms2-V03-00-23/MinBias_Summer09-STARTUP3X_V8D_2360GeV-v2/filtered*.root"); 
  MinBias_skimdata->Add("DataCorrZ2360_skimmednTuple.root"); 
  
  
  TString draw_opt_MinBias_skimdata = Form("%s%s", iso_opt.Data(), ">> _MinBias_skimdata(50,-0.5,2)");
  MinBias_skimdata->Draw(draw_opt_MinBias_skimdata.Data(), baseline_minbias_EE);
  TH1F* Iso03_MinBias_skimdata = (TH1F*)gDirectory->Get("_MinBias_skimdata");
 
  Iso03_MinBias_skimdata->Sumw2();
  Iso03_MinBias_skimdata->SetTitle(iso_opt.Data());
  Iso03_MinBias_skimdata->SetLineColor(kRed);
  //Iso03_MinBias_skimdata->SetAxisRange(0,1,"Y");
  Iso03_MinBias_skimdata->GetYaxis()->SetRangeUser(-1000, 1);
//   Iso03_MinBias_skimdata->SetLineWidth(3);
 
  int   nbins_MinBias_skimdata = Iso03_MinBias_skimdata->GetNbinsX();
  float n_MinBias_skimdata = Iso03_MinBias_skimdata->Integral(0, nbins_MinBias_skimdata+1);
  Iso03_MinBias_skimdata->Scale(1/n_MinBias_skimdata);
  Iso03_MinBias_skimdata->SetBinContent(nbins_MinBias_skimdata, Iso03_MinBias_skimdata->Integral(nbins_MinBias_skimdata, nbins_MinBias_skimdata+1));
  Iso03_MinBias_skimdata->SetBinContent(nbins_MinBias_skimdata+1, 0);

  TString draw_opt_MinBias_mc = Form("%s%s", iso_opt.Data(), ">> _MinBias_mc(50,-0.5,2)");
  MinBias_mc->Draw(draw_opt_MinBias_mc.Data(), baseline_minbias_EE, "sames");
  TH1F* Iso03_MinBias_mc = (TH1F*)gDirectory->Get("_MinBias_mc");
  Iso03_MinBias_mc->SetTitle(iso_opt.Data());
  Iso03_MinBias_mc->GetYaxis()->SetTitle("Events");
  Iso03_MinBias_mc->GetXaxis()->SetTitle("GeV");
  
  //Iso03_MinBias_mc->SetAxisRange(0,1,"Y");
  // Iso03_MinBias_mc->GetYaxis()->SetRangeUser(-0.5, 1.2);
  //  Iso03_MinBias_mc->Sumw2();
  Iso03_MinBias_mc->SetLineColor(kBlue);
//   Iso03_MinBias_mc->SetLineWidth(3);

  int   nbins_MinBias_mc = Iso03_MinBias_mc->GetNbinsX();
  float n_MinBias_mc = Iso03_MinBias_mc->Integral(0, nbins_MinBias_mc+1);
  Iso03_MinBias_mc->Scale(1/n_MinBias_mc);
  Iso03_MinBias_mc->SetBinContent(nbins_MinBias_mc, Iso03_MinBias_mc->Integral(nbins_MinBias_mc, nbins_MinBias_mc+1));
  Iso03_MinBias_mc->SetBinContent(nbins_MinBias_mc+1,0);
  c1->SetLogy(1);
  Double_t xl1=.2, yl1=0.75, xl2=xl1+.2, yl2=yl1+.125;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  //leg->SetLineStyle(0);
  // leg->SetLineColor(0);
 
  leg->AddEntry(Iso03_MinBias_skimdata,"MinBias_skimdata","L");   // h1 and h2 are histogram pointers
  leg->AddEntry(Iso03_MinBias_mc,"MinBias_mc","L");   // h1 and h2 are histogram pointers

  

  setStats(Iso03_MinBias_skimdata,Iso03_MinBias_mc, 0.6, 0.6, false);
 
 
  Iso03_MinBias_mc->Draw("");
  Iso03_MinBias_skimdata->Draw("esames");

 
 
 //  // leg->Draw();
     if(!jura){
     TString save_opt1= Form("%s%s%s","../../../plots/2360GeV/EE_skimdata_mc_CorrZ",iso_opt.Data(),"_V03-00-23.gif" );
     TString save_opt2= Form("%s%s%s","../../../plots/2360GeV/EE_skimdata_mc_",iso_opt.Data(),"_V03-00-23.eps" );
   }
   else{
     TString save_opt1= Form("%s%s%s","../../../plots/2360GeV/EE_skimdata_mc_CorrZ",iso_opt.Data(),"_V03-00-23_jura.gif" );
     TString save_opt2= Form("%s%s%s","../../../plots/2360GeV/EE_skimdata_mc_CorrZ",iso_opt.Data(),"__V03-00-23_jura.eps" ); 

   }
   c1->SaveAs(save_opt1.Data());
   c1->SaveAs(save_opt2.Data());
   //   can->Range(-0.5142857,-3.438716,6.0,0.4504303);
 
 
}



runLooper_Iso_MinBias_skimdata_mc_EE_int( TString iso_opt = "ran_ecalIso03_mu",  bool jura = false){
 

  gStyle->SetHistLineWidth(2);
  gROOT->ForceStyle();
  //gStyle->SetOptStat(1101);
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("can", "can");
 
  TChain *MinBias_skimdata = new TChain("Events");
 
  TChain *MinBias_mc = new TChain("Events");
 
  MinBias_mc->Add("MCCorrZ2360_skimmednTuple.root"); 
  //MinBias_mc->Add("/data/tmp/yanjuntu/cms2/cms2-V03-00-23/MinBias_Summer09-STARTUP3X_V8D_2360GeV-v2/filtered*.root"); 
  MinBias_skimdata->Add("DataCorrZ2360_skimmednTuple.root"); 

  TString draw_opt_MinBias_skimdata = Form("%s%s", iso_opt.Data(), ">> _MinBias_skimdata(40,0,2)");
  MinBias_skimdata->Draw(draw_opt_MinBias_skimdata.Data(), baseline_minbias_EE);
  TH1F* Iso03_MinBias_skimdata = (TH1F*)gDirectory->Get("_MinBias_skimdata");
 
  Iso03_MinBias_skimdata->Sumw2();
  Iso03_MinBias_skimdata->SetTitle(iso_opt.Data());
  Iso03_MinBias_skimdata->SetLineColor(kRed);
  //Iso03_MinBias_skimdata->SetAxisRange(0,1,"Y");
  Iso03_MinBias_skimdata->GetYaxis()->SetRangeUser(-1000, 1);
  Iso03_MinBias_skimdata->GetYaxis()->SetTitle("Efficiency");
  Iso03_MinBias_skimdata->GetXaxis()->SetTitle("GeV");
 
//   Iso03_MinBias_skimdata->SetLineWidth(3);
 
  int   nbins_MinBias_skimdata = Iso03_MinBias_skimdata->GetNbinsX();
  float n_MinBias_skimdata = Iso03_MinBias_skimdata->Integral(0, nbins_MinBias_skimdata+1);
  Iso03_MinBias_skimdata->Scale(1/n_MinBias_skimdata);
  // Iso03_MinBias_skimdata->SetBinContent(nbins_MinBias_skimdata, Iso03_MinBias_skimdata->Integral(nbins_MinBias_skimdata, nbins_MinBias_skimdata+1));
  //Iso03_MinBias_skimdata->SetBinContent(nbins_MinBias_skimdata+1, 0);
  setBin_int_above(Iso03_MinBias_skimdata);

  TString draw_opt_MinBias_mc = Form("%s%s", iso_opt.Data(), ">> _MinBias_mc(40,0,2)");
  MinBias_mc->Draw(draw_opt_MinBias_mc.Data(), baseline_minbias_EE, "sames");
  TH1F* Iso03_MinBias_mc = (TH1F*)gDirectory->Get("_MinBias_mc");
  Iso03_MinBias_mc->Sumw2();
  Iso03_MinBias_mc->SetTitle(iso_opt.Data());
  //Iso03_MinBias_mc->SetAxisRange(0,1,"Y");
  // Iso03_MinBias_mc->GetYaxis()->SetRangeUser(-0.5, 1.2);
  
  Iso03_MinBias_mc->SetLineColor(kBlue);

  Iso03_MinBias_mc->GetYaxis()->SetTitle("Efficiency");
  Iso03_MinBias_mc->GetXaxis()->SetTitle("GeV");
//   Iso03_MinBias_mc->SetLineWidth(3);

  int   nbins_MinBias_mc = Iso03_MinBias_mc->GetNbinsX();
  float n_MinBias_mc = Iso03_MinBias_mc->Integral(0, nbins_MinBias_mc+1);
  Iso03_MinBias_mc->Scale(1/n_MinBias_mc);
  //Iso03_MinBias_mc->SetBinContent(nbins_MinBias_mc, Iso03_MinBias_mc->Integral(nbins_MinBias_mc, nbins_MinBias_mc+1));
  //Iso03_MinBias_mc->SetBinContent(nbins_MinBias_mc+1,0);

  setBin_int_above(Iso03_MinBias_mc);
  // c1->SetLogy(1);
  Double_t xl1=.5, yl1=0.35, xl2=xl1+.3, yl2=yl1+.225;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  //leg->SetLineStyle(0);
  // leg->SetLineColor(0);
 
  leg->AddEntry(Iso03_MinBias_skimdata,"MinBias_Data","L");   // h1 and h2 are histogram pointers
  leg->AddEntry(Iso03_MinBias_mc,"MinBias_MC","L");   // h1 and h2 are histogram pointers

  

  //setStats(Iso03_MinBias_skimdata,Iso03_MinBias_mc, 0.6, 0.6, false);
 
 
  Iso03_MinBias_mc->Draw("e");
  Iso03_MinBias_skimdata->Draw("esames");

 
 
  leg->Draw();
  if(!jura){
    TString save_opt1= Form("%s%s%s","../../../plots/2360GeV/EE_skimdata_mc_intCorrZ",iso_opt.Data(),"_V03-00-23.gif" );
    TString save_opt2= Form("%s%s%s","../../../plots/2360GeV/EE_skimdata_mc_intCorrZ",iso_opt.Data(),"_V03-00-23.eps" );
  }
  else{
    TString save_opt1= Form("%s%s%s","../../../plots/2360GeV/EE_skimdata_mc_intCorrZ",iso_opt.Data(),"_V03-00-23_jura.gif" );
    TString save_opt2= Form("%s%s%s","../../../plots/2360GeV/EE_skimdata_mc_intCorrZ",iso_opt.Data(),"_V03-00-23_jura.eps" ); 

  }
  c1->SaveAs(save_opt1.Data());
  c1->SaveAs(save_opt2.Data());
  //can->Range(-0.5142857,-3.438716,6.0,0.4504303);
 
 
}


draw(TString iso_opt_1 = "ran_trkIso03_mu"){
  runLooper_Iso_MinBias_skimdata_mc_EE(iso_opt_1);
  runLooper_Iso_MinBias_skimdata_mc_EE_int(iso_opt_1);

 
}




/*
runLooper_Iso_MinBias_skimdata_mc( TString iso_opt = "ran_trkIso03_mu", bool jura = false){
 

  gStyle->SetHistLineWidth(2);
  gROOT->ForceStyle();
  gStyle->SetOptStat(1101);

  TCanvas *c1 = new TCanvas("can", "can");
  TFile *file = TFile::Open("myHist.root");
   TChain *MinBias_skimdata = new TChain("Events");
  TChain *MinBias_mc = new TChain("Events");
  MinBias_mc->Add("/data/tmp/yanjuntu/cms2/cms2-V03-00-23/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/filtered*.root"); 
  //MinBias_mc->Add("MC_skimmednTuple.root"); 
  MinBias_skimdata->Add("Data_skimmednTuple.root"); 
  // TH1F* Iso03_MinBias_skimdata = (TH1F*)file->FindObjectAny("Data_ran_trkIso03_mu_00");
  
  TString draw_opt_MinBias_skimdata = Form("%s%s", iso_opt.Data(), ">> _MinBias_skimdata(100,-1,3)");
  MinBias_skimdata->Draw(draw_opt_MinBias_skimdata.Data(), baseline_minbias_trk);
  TH1F* Iso03_MinBias_skimdata = (TH1F*)gDirectory->Get("_MinBias_skimdata");
 
  Iso03_MinBias_skimdata->Sumw2();
  Iso03_MinBias_skimdata->SetTitle(iso_opt.Data());
  Iso03_MinBias_skimdata->SetLineColor(kRed);
  //Iso03_MinBias_skimdata->SetAxisRange(0,2,"X");
//   Iso03_MinBias_skimdata->SetLineWidth(3);
 
  int   nbins_MinBias_skimdata = Iso03_MinBias_skimdata->GetNbinsX();
  float n_MinBias_skimdata = Iso03_MinBias_skimdata->Integral(0, nbins_MinBias_skimdata+1);
  Iso03_MinBias_skimdata->Scale(1/n_MinBias_skimdata);
  Iso03_MinBias_skimdata->SetBinContent(nbins_MinBias_skimdata, Iso03_MinBias_skimdata->Integral(nbins_MinBias_skimdata, nbins_MinBias_skimdata+1));
  Iso03_MinBias_skimdata->SetBinContent(nbins_MinBias_skimdata+1, 0);

  
 
  TString draw_opt_MinBias_mc = Form("%s%s", iso_opt.Data(), ">> _MinBias_mc(100,-1,3)");
  MinBias_mc->Draw(draw_opt_MinBias_mc.Data(), baseline_minbias_trk, "sames");
  TH1F* Iso03_MinBias_mc = (TH1F*)gDirectory->Get("_MinBias_mc");

  // TH1F* Iso03_MinBias_mc = (TH1F*)file->FindObjectAny("MC_ran_trkIso03_mu_00");
  Iso03_MinBias_mc->SetAxisRange(0,2,"Y");
  //  Iso03_MinBias_mc->Sumw2();
  Iso03_MinBias_mc->SetLineColor(kBlue);
//   Iso03_MinBias_mc->SetLineWidth(3);

  int   nbins_MinBias_mc = Iso03_MinBias_mc->GetNbinsX();
  float n_MinBias_mc = Iso03_MinBias_mc->Integral(0, nbins_MinBias_mc+1);
  Iso03_MinBias_mc->Scale(1/n_MinBias_mc);
  Iso03_MinBias_mc->SetBinContent(nbins_MinBias_mc, Iso03_MinBias_mc->Integral(nbins_MinBias_mc, nbins_MinBias_mc+1));
  Iso03_MinBias_mc->SetBinContent(nbins_MinBias_mc+1,0);
  Iso03_MinBias_mc->SetTitle(iso_opt.Data());
  c1->SetLogy(1);
  Double_t xl1=.2, yl1=0.75, xl2=xl1+.2, yl2=yl1+.125;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  //leg->SetLineStyle(0);
  // leg->SetLineColor(0);
 
  leg->AddEntry(Iso03_MinBias_skimdata,"MinBias_skimdata","L");   // h1 and h2 are histogram pointers
  leg->AddEntry(Iso03_MinBias_mc,"MinBias_mc","L");   // h1 and h2 are histogram pointers
  setStats(Iso03_MinBias_skimdata,Iso03_MinBias_mc, 0.6, 0.6, false);
 

  Iso03_MinBias_mc->Draw("");
  Iso03_MinBias_skimdata->Draw("esames");
 
//   // leg->Draw();
//     if(!jura){
//     TString save_opt1= Form("%s%s%s","plots/900GeV/skimdata_mc_",iso_opt.Data(),"_V03-00-20.gif" );
//     TString save_opt2= Form("%s%s%s","plots/900GeV/skimdata_mc_",iso_opt.Data(),"_V03-00-20.eps" );
//   }
//   else{
//     TString save_opt1= Form("%s%s%s","plots/900GeV/skimdata_mc_",iso_opt.Data(),"_V03-00-20_jura.gif" );
//     TString save_opt2= Form("%s%s%s","plots/900GeV/skimdata_mc_",iso_opt.Data(),"_V03-00-20_jura.eps" ); 

//   }
//   c1->SaveAs(save_opt1.Data());
//   c1->SaveAs(save_opt2.Data());
  //can->Range(-0.5142857,-3.438716,6.0,0.4504303);
 
 
}


command{

 TChain *MinBias_skimdata = new TChain("Events");
  TChain *MinBias_mc = new TChain("Events");
  //MinBias_mc->Add("/data/tmp/yanjuntu/cms2/cms2-V03-00-23/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/filtered*.root"); 
  MinBias_mc->Add("MC_skimmednTuple.root"); 
  MinBias_skimdata->Add("Data_skimmednTuple.root"); 

  MinBias_skimdata->Draw("ran_trksp4");
  

}


runLooper_Iso_MinBias_skimdata_mc_EE_int( TString iso_opt = "ran_trkIso03_mu",  bool jura = false){
 

  gStyle->SetHistLineWidth(2);
  gROOT->ForceStyle();
  gStyle->SetOptStat(1101);
  //gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("can", "can");
 
  TChain *MinBias_skimdata = new TChain("Events");
 
  TChain *MinBias_mc = new TChain("Events");
 
  MinBias_mc->Add("MC_skimmednTuple.root"); 
 
  //MinBias_mc->Add("/data/tmp/yanjuntu/cms2/cms2-V03-00-23/MinBias_Summer09-STARTUP3X_V8I_900GeV-v2/filtered*.root"); 
  //MinBias_mc->Add("MC_skimmednTuple.root"); 
  MinBias_skimdata->Add("Data_skimmednTuple.root"); 

  
  
  TString draw_opt_MinBias_skimdata = Form("%s%s", iso_opt.Data(), ">> _MinBias_skimdata(40,0,2)");
  MinBias_skimdata->Draw(draw_opt_MinBias_skimdata.Data(), baseline_minbias_trk);
  TH1F* Iso03_MinBias_skimdata = (TH1F*)gDirectory->Get("_MinBias_skimdata");
 
  Iso03_MinBias_skimdata->Sumw2();
  Iso03_MinBias_skimdata->SetTitle(iso_opt.Data());
  Iso03_MinBias_skimdata->SetLineColor(kRed);
  //Iso03_MinBias_skimdata->SetAxisRange(0,1,"Y");
  Iso03_MinBias_skimdata->GetYaxis()->SetRangeUser(-1000, 1);
//   Iso03_MinBias_skimdata->SetLineWidth(3);
 
  int   nbins_MinBias_skimdata = Iso03_MinBias_skimdata->GetNbinsX();
  float n_MinBias_skimdata = Iso03_MinBias_skimdata->Integral(0, nbins_MinBias_skimdata+1);
  Iso03_MinBias_skimdata->Scale(1/n_MinBias_skimdata);
  // Iso03_MinBias_skimdata->SetBinContent(nbins_MinBias_skimdata, Iso03_MinBias_skimdata->Integral(nbins_MinBias_skimdata, nbins_MinBias_skimdata+1));
  //Iso03_MinBias_skimdata->SetBinContent(nbins_MinBias_skimdata+1, 0);
  setBin_int_above(Iso03_MinBias_skimdata);

  TString draw_opt_MinBias_mc = Form("%s%s", iso_opt.Data(), ">> _MinBias_mc(40,0,2)");
  MinBias_mc->Draw(draw_opt_MinBias_mc.Data(), baseline_minbias_trk, "sames");
  TH1F* Iso03_MinBias_mc = (TH1F*)gDirectory->Get("_MinBias_mc");
  Iso03_MinBias_mc->Sumw2();
  Iso03_MinBias_mc->SetTitle(iso_opt.Data());
  //Iso03_MinBias_mc->SetAxisRange(0,1,"Y");
  // Iso03_MinBias_mc->GetYaxis()->SetRangeUser(-0.5, 1.2);
  //  Iso03_MinBias_mc->Sumw2();
  Iso03_MinBias_mc->SetLineColor(kBlue);
//   Iso03_MinBias_mc->SetLineWidth(3);

  int   nbins_MinBias_mc = Iso03_MinBias_mc->GetNbinsX();
  float n_MinBias_mc = Iso03_MinBias_mc->Integral(0, nbins_MinBias_mc+1);
  Iso03_MinBias_mc->Scale(1/n_MinBias_mc);
  //Iso03_MinBias_mc->SetBinContent(nbins_MinBias_mc, Iso03_MinBias_mc->Integral(nbins_MinBias_mc, nbins_MinBias_mc+1));
  //Iso03_MinBias_mc->SetBinContent(nbins_MinBias_mc+1,0);

  setBin_int_above(Iso03_MinBias_mc);
  // c1->SetLogy(1);
  Double_t xl1=.2, yl1=0.75, xl2=xl1+.2, yl2=yl1+.125;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  //leg->SetLineStyle(0);
  // leg->SetLineColor(0);
 
  leg->AddEntry(Iso03_MinBias_skimdata,"MinBias_skimdata","L");   // h1 and h2 are histogram pointers
  leg->AddEntry(Iso03_MinBias_mc,"MinBias_mc","L");   // h1 and h2 are histogram pointers

  

  setStats(Iso03_MinBias_skimdata,Iso03_MinBias_mc, 0.6, 0.6, false);
 
 
  Iso03_MinBias_mc->Draw("e");
  Iso03_MinBias_skimdata->Draw("esames");

 
 
 //  // leg->Draw();
//   if(!jura){
//     TString save_opt1= Form("%s%s%s","plots/900GeV/EE_skimdata_mc_",iso_opt.Data(),"_V03-00-23.gif" );
//     TString save_opt2= Form("%s%s%s","plots/900GeV/EE_skimdata_mc_",iso_opt.Data(),"_V03-00-23.eps" );
//   }
//   else{
//     TString save_opt1= Form("%s%s%s","plots/900GeV/EE_skimdata_mc_",iso_opt.Data(),"_V03-00-23_jura.gif" );
//     TString save_opt2= Form("%s%s%s","plots/900GeV/EE_skimdata_mc_",iso_opt.Data(),"__V03-00-23_jura.eps" ); 

//   }
//   c1->SaveAs(save_opt1.Data());
//   c1->SaveAs(save_opt2.Data());
  //can->Range(-0.5142857,-3.438716,6.0,0.4504303);
 
 
}
*/
