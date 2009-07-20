#ifndef IOStruct
#define IOStruct 1
struct IOStringStruct {
  std::string selection;

  std::string mcMatchSig;
  std::string mcMatchBkg;
  std::string isoSumSt;
  std::string isoCalSt;
  std::string isoTrkSt;
  std::string evtScaleSt;

};

#endif

void isoPerformanceElectronsHyp(TChain* sigCh, TChain* bkgCh,
		    double minX, double maxX, double minY, double maxY){

  IOStringStruct config; config.evtScaleSt = std::string("evt_scale1fb");
  config.selection = std::string("els_tightId[hyp_lt_index]&&abs(els_d0[hyp_lt_index])<0.04");
  config.mcMatchSig = std::string("hyp_type==2&&abs(hyp_lt_mc_motherid)==0&&abs(hyp_lt_mc_id)==11");
  config.mcMatchBkg = std::string("hyp_type==2&&abs(hyp_ll_mc_motherid)==0&&abs(hyp_ll_mc_id)==13");
  
  //  AutoLibraryLoader::enable();
  config.isoSumSt = std::string("min(els_p4[hyp_lt_index].pt()/(els_p4[hyp_lt_index].pt()+els_tkIso[hyp_lt_index]+max(0,els_tq_caloIso[hyp_lt_index])),0.995)");
  config.isoCalSt = std::string("min(els_p4[hyp_lt_index].pt()/(els_p4[hyp_lt_index].pt()+max(0,els_tq_caloIso[hyp_lt_index])),0.995)");
  config.isoTrkSt = std::string("min(els_p4[hyp_lt_index].pt()/(els_p4[hyp_lt_index].pt()+els_tkIso[hyp_lt_index]),0.995)");
  isoPerformance(sigCh, bkgCh, config, minX, maxX, minY, maxY);  
}

void isoPerformanceElectronsNoHypWj(TChain* sigCh, TChain* bkgCh,
		    double minX, double maxX, double minY, double maxY){

  IOStringStruct config; config.evtScaleSt = std::string("evt_scale1fb");
  config.selection = std::string("els_tightId&&abs(els_d0)<0.04&&els_p4.pt()>20");
  config.mcMatchSig = std::string("abs(els_mc_id)==11&&(abs(els_mc_motherid)==24||abs(els_mc_motherid)==0||abs(els_mc_motherid)==23)");
  config.mcMatchBkg = std::string("Sum$(abs(mus_mc_motherid)==0||abs(mus_mc_motherid)==23||abs(mus_mc_motherid)==24)>0");
  
  //  AutoLibraryLoader::enable();
  config.isoSumSt = std::string("min(els_p4.pt()/(els_p4.pt()+els_tkIso+max(0,els_tq_caloIso)),0.995)");
  config.isoCalSt = std::string("min(els_p4.pt()/(els_p4.pt()+max(0,els_tq_caloIso)),0.995)");
  config.isoTrkSt = std::string("min(els_p4.pt()/(els_p4.pt()+els_tkIso),0.995)");
  isoPerformance(sigCh, bkgCh, config, minX, maxX, minY, maxY);  
}
void isoPerformanceElectronsNoHypWjCMS2(TChain* sigCh, TChain* bkgCh,
		    double minX, double maxX, double minY, double maxY){

  IOStringStruct config; config.evtScaleSt = std::string("evt_scale1fb");
  config.selection = std::string("els_tightId&&abs(els_d0)<0.04&&els_p4.pt()>20");
  config.mcMatchSig = std::string("abs(els_mc_id)==11&&(abs(els_mc_motherid)==24||abs(els_mc_motherid)==0||abs(els_mc_motherid)==23)");
  config.mcMatchBkg = std::string("Sum$(abs(mus_mc_motherid)==0||abs(mus_mc_motherid)==23||abs(mus_mc_motherid)==24)>0");
  
  //  AutoLibraryLoader::enable();
  config.isoSumSt = std::string("min(els_p4.pt()/(els_p4.pt()+els_tkIso+max(0,els_ecalJuraIso+els_hcalConeIso)),0.995)");
  config.isoCalSt = std::string("min(els_p4.pt()/(els_p4.pt()+max(0,els_ecalJuraIso+els_hcalConeIso)),0.995)");
  config.isoTrkSt = std::string("min(els_p4.pt()/(els_p4.pt()+els_tkIso),0.995)");
  isoPerformance(sigCh, bkgCh, config, minX, maxX, minY, maxY);  
}

void isoPerformanceElectronsNoHypWjCMS2SW2pat(TChain* sigCh, TChain* bkgCh,
		    double minX, double maxX, double minY, double maxY){

  IOStringStruct config; config.evtScaleSt = std::string("evt_scale1fb");
  config.selection = std::string("els_looseId&&abs(els_d0corr)<0.04&&els_p4.pt()>20");
  config.mcMatchSig = std::string("abs(els_mc_id)==11&&(abs(els_mc_motherid)==24||abs(els_mc_motherid)==0||abs(els_mc_motherid)==23)");
  //pick electrons from QCD here
  config.mcMatchBkg = std::string("1<2");
  //protect just in case
  if (sigCh == bkgCh){
    std::cout<<"These better be different input chains"<<std::endl;
    return;
  }
  //  config.mcMatchBkg = std::string("Sum$(abs(mus_mc_motherid)==0||abs(mus_mc_motherid)==23||abs(mus_mc_motherid)==24)>0");
  
  //  AutoLibraryLoader::enable();
  config.isoSumSt = std::string("min(els_p4.pt()/(els_p4.pt()+els_pat_trackIso+max(0,els_pat_ecalIso+els_pat_hcalIso)),0.995)");
  config.isoCalSt = std::string("min(els_p4.pt()/(els_p4.pt()+max(0,els_pat_ecalIso+els_pat_hcalIso)),0.995)");
  config.isoTrkSt = std::string("min(els_p4.pt()/(els_p4.pt()+els_pat_trackIso),0.995)");
  config.evtScaleSt = std::string("evt_scale1fb");
  isoPerformance(sigCh, bkgCh, config, minX, maxX, minY, maxY);  
}

void isoPerformanceMuonsNoHypWjCMS2SW2pat(TChain* sigCh, TChain* bkgCh,
		    double minX, double maxX, double minY, double maxY){

  IOStringStruct config; config.evtScaleSt = std::string("evt_scale1fb");
  config.selection = std::string("mus_gfit_chi2/mus_gfit_ndof<10&&mus_validHits>10&&(mus_type&2)==2&&mus_p4.pt()>20");
  config.mcMatchSig = std::string("abs(mus_mc_id)==13&&(abs(mus_mc_motherid)==24||abs(mus_mc_motherid)==0||abs(mus_mc_motherid)==23)");
  //pick muons from QCD here
  config.mcMatchBkg = std::string("1<2");
  //protect just in case
  if (sigCh == bkgCh){
    std::cout<<"These better be different input chains"<<std::endl;
    return;
  }
  //  config.mcMatchBkg = std::string("Sum$(abs(mus_mc_motherid)==0||abs(mus_mc_motherid)==23||abs(mus_mc_motherid)==24)>0");
  
  //  AutoLibraryLoader::enable();
  config.isoSumSt = std::string("min(mus_p4.pt()/(mus_p4.pt()+mus_pat_trackIso+max(0,mus_pat_ecalIso+mus_pat_hcalIso)),0.995)");
  config.isoCalSt = std::string("min(mus_p4.pt()/(mus_p4.pt()+max(0,mus_pat_ecalIso+mus_pat_hcalIso)),0.995)");
  config.isoTrkSt = std::string("min(mus_p4.pt()/(mus_p4.pt()+mus_pat_trackIso),0.995)");
  config.evtScaleSt = std::string("evt_scale1fb");
  isoPerformance(sigCh, bkgCh, config, minX, maxX, minY, maxY);  
}


void isoPerformanceMuonsHyp(TChain* sigCh, TChain* bkgCh,
		    double minX, double maxX, double minY, double maxY){

  IOStringStruct config; config.evtScaleSt = std::string("evt_scale1fb");
  config.selection = std::string("mus_gfit_chi2[hyp_ll_index]/mus_gfit_ndof[hyp_ll_index]<5&&abs(mus_d0[hyp_ll_index])<0.25&&mus_validHits[hyp_ll_index]>7");
  config.mcMatchSig = std::string("hyp_type==2&&abs(hyp_ll_mc_motherid)==0&&abs(hyp_ll_mc_id)==13");
  config.mcMatchBkg = std::string("hyp_type==2&&abs(hyp_lt_mc_motherid)==0&&abs(hyp_lt_mc_id)==11");
  
  //  AutoLibraryLoader::enable();
  config.isoSumSt = std::string("min(mus_p4[hyp_ll_index].pt()/(mus_p4[hyp_ll_index].pt()+mus_iso03_sumPt[hyp_ll_index]+max(0,mus_iso03_hadEt[hyp_ll_index]+mus_iso03_emEt[hyp_ll_index])),0.995)");
  config.isoCalSt = std::string("min(mus_p4[hyp_ll_index].pt()/(mus_p4[hyp_ll_index].pt()+max(0,mus_iso03_hadEt[hyp_ll_index]+mus_iso03_emEt[hyp_ll_index])),0.995)");
  config.isoTrkSt = std::string("min(mus_p4[hyp_ll_index].pt()/(mus_p4[hyp_ll_index].pt()+mus_iso03_sumPt[hyp_ll_index]),0.995)");
  isoPerformance(sigCh, bkgCh, config, minX, maxX, minY, maxY);  
}

void isoPerformanceMuonsNoHyp(TChain* sigCh, TChain* bkgCh,
			      double minX, double maxX, double minY, double maxY){

  IOStringStruct config; config.evtScaleSt = std::string("evt_scale1fb");
  config.selection = std::string("mus_gfit_chi2/mus_gfit_ndof<5&&abs(mus_d0)<0.25&&mus_validHits>7&&mus_trk_p4.pt()>20&&abs(mus_p4.eta())<2.5");
  config.mcMatchSig = std::string("abs(mus_mc_id)==13&&(abs(mus_mc_motherid)==24||abs(mus_mc_motherid)==0||abs(mus_mc_motherid)==23)");
  config.mcMatchBkg = std::string("1<2");
  
  //  AutoLibraryLoader::enable();
  config.isoSumSt = std::string("min(mus_p4.pt()/(mus_p4.pt()+mus_iso03_sumPt+max(0,mus_iso03_hadEt+mus_iso03_emEt)),0.995)");
  config.isoCalSt = std::string("min(mus_p4.pt()/(mus_p4.pt()+max(0,mus_iso03_hadEt+mus_iso03_emEt)),0.995)");
  config.isoTrkSt = std::string("min(mus_p4.pt()/(mus_p4.pt()+mus_iso03_sumPt),0.995)");
  isoPerformance(sigCh, bkgCh, config, minX, maxX, minY, maxY);  
}

void isoPerformance(TChain* sigCh, TChain* bkgCh, IOStringStruct& cfg,
		    double minX, double maxX, double minY, double maxY){
  gROOT->SetStyle("Plain");
  //  gSystem->Load("libFWCoreFWLite");

  bkgCh->Draw(Form("%s:%s>>bkg2d(100,0,1,100,0,1)", cfg.isoTrkSt.c_str(), cfg.isoCalSt.c_str()), 
	      Form("%s*(%s&&%s)", cfg.evtScaleSt.c_str(), cfg.mcMatchBkg.c_str(), cfg.selection.c_str()), 
	      "box");
  sigCh->Draw(Form("%s:%s>>sig2d(100,0,1,100,0,1)", cfg.isoTrkSt.c_str(), cfg.isoCalSt.c_str()), 
	      Form("%s*(%s&&%s)", cfg.evtScaleSt.c_str(), cfg.mcMatchSig.c_str(), cfg.selection.c_str()), 
	      "box");
  bkg2d = (TH2F*)gDirectory->Get("bkg2d"); bkg2d->Scale(1./bkg2d->Integral());
  sig2d = (TH2F*)gDirectory->Get("sig2d"); sig2d->Scale(1./sig2d->Integral());


  bkgCh->Draw(Form("%s>>bkgTrk(100,0,1)", cfg.isoTrkSt.c_str()), 
	      Form("%s*(%s&&%s)", cfg.evtScaleSt.c_str(), cfg.mcMatchBkg.c_str(), cfg.selection.c_str()), 
	      "box");
  sigCh->Draw(Form("%s>>sigTrk(100,0,1)", cfg.isoTrkSt.c_str()), 
	      Form("%s*(%s&&%s)", cfg.evtScaleSt.c_str(), cfg.mcMatchSig.c_str(), cfg.selection.c_str()), 
	      "box");  
  bkgTrk = (TH1F*)gDirectory->Get("bkgTrk"); bkgTrk->Scale(1./bkgTrk->Integral());
  sigTrk = (TH1F*)gDirectory->Get("sigTrk"); sigTrk->Scale(1./sigTrk->Integral());

  bkgCh->Draw(Form("%s>>bkgCal(100,0,1)", cfg.isoCalSt.c_str()), 
	      Form("%s*(%s&&%s)", cfg.evtScaleSt.c_str(), cfg.mcMatchBkg.c_str(), cfg.selection.c_str()), 
	      "box");
  sigCh->Draw(Form("%s>>sigCal(100,0,1)", cfg.isoCalSt.c_str()), 
	      Form("%s*(%s&&%s)", cfg.evtScaleSt.c_str(), cfg.mcMatchSig.c_str(), cfg.selection.c_str()), 
	      "box");  
  bkgCal = (TH1F*)gDirectory->Get("bkgCal"); bkgCal->Scale(1./bkgCal->Integral());
  sigCal = (TH1F*)gDirectory->Get("sigCal"); sigCal->Scale(1./sigCal->Integral());

  bkgCh->Draw(Form("%s:0.995>>bkgsum(100,0,1,100,0,1)", cfg.isoSumSt.c_str()),
	      Form("%s*(%s&&%s)", cfg.evtScaleSt.c_str(), cfg.mcMatchBkg.c_str(), cfg.selection.c_str()), 
	      "box");
  sigCh->Draw(Form("%s:0.995>>sigsum(100,0,1,100,0,1)", cfg.isoSumSt.c_str()),
	      Form("%s*(%s&&%s)", cfg.evtScaleSt.c_str(), cfg.mcMatchSig.c_str(), cfg.selection.c_str()), 
	      "box");
  bkgsum = (TH1F*)gDirectory->Get("bkgsum"); bkgsum->Scale(1./bkgsum->Integral());
  sigsum = (TH1F*)gDirectory->Get("sigsum"); sigsum->Scale(1./sigsum->Integral());

  isoPerformance(bkg2d, sig2d, bkgTrk, sigTrk, bkgCal, sigCal, bkgsum, sigsum,
		 minX, maxX, minY, maxY);
}

void isoPerformance(TH2F* bkg2d, TH2F* sig2d, 
		    TH1F* bkgTrk, TH1F* sigTrk, TH1F* bkgCal, TH1F* sigCal, 
		    TH1F* bkgsum, TH1F* sigsum,
		    double minX, double maxX, double minY, double maxY){

  int nBinsXTrk = sigTrk->GetNbinsX();
  if (nBinsXTrk != bkgTrk->GetNbinsX()){ std::cout<<"Incompatible hists in Trk "<<std::endl; return;}
  int nBinsXCal = sigCal->GetNbinsX();
  if (nBinsXCal != bkgCal->GetNbinsX()){ std::cout<<"Incompatible hists in Cal "<<std::endl; return;}


  TH1F* lhSigTrk = new TH1F("lhSigTrk", "lhSigTrk", 1000,0, 1);
  TH1F* lhBkgTrk = new TH1F("lhBkgTrk", "lhBkgTrk", 1000,0, 1);
  for(int i=1; i<= nBinsXTrk; ++i){
    double fS = sigTrk->GetBinContent(i); double fB = bkgTrk->GetBinContent(i);
    double lhVal = 0.5;
    if ((fS + fB) != 0){ lhVal = fS/(fS+fB);}
    lhVal = lhVal > 0.9999 ? 0.9999 : lhVal;
    lhSigTrk->Fill(lhVal, fS);
    lhBkgTrk->Fill(lhVal, fB);
  }
  std::cout<<lhSigTrk->Integral()<<std::endl;
  lhSigTrk->Scale(1./lhSigTrk->Integral());
  lhBkgTrk->Scale(1./lhBkgTrk->Integral());

  TH1F* lhSigProd = new TH1F("lhSigProd", "lhSigProd", 1000,0, 1);
  TH1F* lhBkgProd = new TH1F("lhBkgProd", "lhBkgProd", 1000,0, 1);
  for(int i=1; i<= nBinsXTrk; ++i) for(int j=1; j<= nBinsXCal; ++j){
    double fS = sigTrk->GetBinContent(i)*sigCal->GetBinContent(j); 
    double fB = bkgTrk->GetBinContent(i)*bkgCal->GetBinContent(j);
    double lhVal = 0.5;
    if ((fS + fB) != 0) lhVal = fS/(fS+fB);
    lhVal = lhVal > 0.9999 ? 0.9999 : lhVal;
    lhSigProd->Fill(lhVal, sig2d->GetBinContent(i,j));
    lhBkgProd->Fill(lhVal, bkg2d->GetBinContent(i,j));
  }
  lhSigProd->Scale(1./lhSigProd->Integral());
  lhBkgProd->Scale(1./lhBkgProd->Integral());

  int nBins2DX = sig2d->GetNbinsX(); int nBins2DY = sig2d->GetNbinsY();
  if (nBins2DX != bkg2d->GetNbinsX() || nBins2DY != bkg2d->GetNbinsY()){
    std::cout<<"Incompatible bins in 2D "<<std::endl;
    return;
  }
  TH1F* lhSig2D = new TH1F("lhSig2D", "lhSig2D", 1000,0, 1);
  TH1F* lhBkg2D = new TH1F("lhBkg2D", "lhBkg2D", 1000,0, 1);
  for(int i=1; i<= nBins2DX; ++i) for(int j=1; j<= nBins2DY; ++j){
    double fS = sig2d->GetBinContent(i,j); 
    double fB = bkg2d->GetBinContent(i,j);
    double lhVal = 0.5;
    if ((fS + fB) != 0) lhVal = fS/(fS+fB);
    lhVal = lhVal > 0.9999 ? 0.9999 : lhVal;
    lhSig2D->Fill(lhVal, sig2d->GetBinContent(i,j));
    lhBkg2D->Fill(lhVal, bkg2d->GetBinContent(i,j));
  }
  lhSig2D->Scale(1./lhSig2D->Integral());
  lhBkg2D->Scale(1./lhBkg2D->Integral());


  TH2F* rocAND = new TH2F("rocAND", "rocAND", 16000, 0, 1+1e-6, 1000, 0, 1+1e-6);
  for(int i=1;i<=100;++i)for(int j=1;j<=100;++j){
    rocAND->Fill(bkg2d->Integral(i,100,j,100), sig2d->Integral(i,100,j,100));
  }

  TH2F* rocAND0909 = new TH2F("rocAND0909", "rocAND0909", 16000, 0, 1+1e-6, 1000, 0, 1+1e-6);
  //FIXME: assumes 100 bins here now, need to figure out for a general case
  rocAND0909->Fill(bkg2d->Integral(90,nBins2DX,90,nBins2DY), 
		   sig2d->Integral(90,nBins2DX,90,nBins2DY));
  rocAND0909->Fill(bkg2d->Integral(80,nBins2DX,90,nBins2DY), 
		   sig2d->Integral(80,nBins2DX,90,nBins2DY));
  rocAND0909->Fill(bkg2d->Integral(85,nBins2DX,90,nBins2DY), 
		   sig2d->Integral(85,nBins2DX,90,nBins2DY));
  rocAND0909->Fill(bkg2d->Integral(92,nBins2DX,92,nBins2DY), 
		   sig2d->Integral(92,nBins2DX,92,nBins2DY));

  TH2F* rocNoCal = new TH2F("rocNoCal", "rocNoCal", 16000, 0, 1+1e-6, 1000, 0, 1+1e-6);
  for(int i=1;i<=100;++i){
    rocNoCal->Fill(bkg2d->Integral(1,nBins2DX,i,nBins2DY), 
		   sig2d->Integral(1,nBins2DX,i,nBins2DY));
  }

  TH2F* rocSum = new TH2F("rocSum", "rocSum", 16000, 0, 1+1e-6, 1000, 0, 1+1e-6);
  for(int i=1;i<=100;++i){
    rocSum->Fill(bkgsum->Integral(1,100,i,100), 
		 sigsum->Integral(1,100,i,100));
  }
  rocSum092 = new TH2F("rocSum092", "rocSum092", 16000, 0, 1+1e-6, 1000, 0, 1+1e-6);
  rocSum092->Fill(bkgsum->Integral(1,100,85,100), sigsum->Integral(1,100,85,100));
  rocSum092->Fill(bkgsum->Integral(1,100,90,100), sigsum->Integral(1,100,90,100));
  rocSum092->Fill(bkgsum->Integral(1,100,92,100), sigsum->Integral(1,100,92,100));

  TH2F* rocDiag = new TH2F("rocDiag", "rocDiag", 16000, 0, 1+1e-6, 1000, 0, 1+1e-6);
  for(int i=1;i<=100;++i){
    rocDiag->Fill(bkg2d->Integral(i,100,i,100), sig2d->Integral(i,100,i,100));
  }

  TH2F* rocLHTrk = new TH2F("rocLHTrk", "rocLHTrk", 16000, 0, 1+1e-6, 1000, 0, 1+1e-6);
  for(int i=1;i<=1000;++i){
    rocLHTrk->Fill(lhBkgTrk->Integral(i,1000), lhSigTrk->Integral(i,1000));
  }

  TH2F* rocLHProd = new TH2F("rocLHProd", "rocLHProd", 16000, 0, 1+1e-6, 1000, 0, 1+1e-6);
  for(int i=1;i<=1000;++i){
    rocLHProd->Fill(lhBkgProd->Integral(i,1000), lhSigProd->Integral(i,1000));
  }
  TH2F* rocLH2D = new TH2F("rocLH2D", "rocLH2D", 16000, 0, 1+1e-6, 1000, 0, 1+1e-6);
  for(int i=1;i<=1000;++i){
    rocLH2D->Fill(lhBkg2D->Integral(i,1000), lhSig2D->Integral(i,1000));
  }

  rocAND->SetMarkerStyle(20);
  rocAND->SetMarkerSize(0.5);
  xax = rocAND->GetXaxis();
  xax->SetRangeUser(minX,maxX);
  yax = rocAND->GetYaxis();
  yax->SetRangeUser(minY, maxY);
  rocAND->SetStats(0);
  rocAND->Draw("*");

  rocSum->SetMarkerColor(2);
  rocSum->SetMarkerStyle(21);
  rocSum->SetMarkerSize(0.75);
  rocSum->Draw("* same");

  rocNoCal->SetMarkerSize(1.5);
  rocNoCal->SetMarkerStyle(23);
  rocNoCal->SetMarkerColor(6);
  //  rocNoCal->Draw("* same");

  rocLHTrk->SetMarkerSize(1.);
  rocLHTrk->SetMarkerStyle(29);
  rocLHTrk->SetMarkerColor(5);
  //  rocLHTrk->Draw("* same");

  rocLHProd->SetMarkerSize(1.);
  rocLHProd->SetMarkerStyle(30);
  rocLHProd->SetMarkerColor(8);
  //  rocLHProd->Draw("* same");

  rocLH2D->SetMarkerSize(1.);
  rocLH2D->SetMarkerStyle(31);
  rocLH2D->SetMarkerColor(4);
  //  rocLH2D->Draw("* same");

  rocDiag->SetMarkerSize(1.5);
  rocDiag->SetMarkerColor(3);
  rocDiag->SetMarkerStyle(20);
  rocDiag->Draw("* same");

  rocAND0909->SetMarkerStyle(22);
  rocAND0909->SetMarkerSize(1.75);
  rocAND0909->SetMarkerColor(4);
  rocAND0909->Draw("* same");

  rocSum092->SetMarkerColor(4);
  rocSum092->SetMarkerSize(2.25);
  rocSum092->SetMarkerStyle(29);
  rocSum092->Draw("* same");

  gPad->SaveAs("perfIso.png");
}
