void plotHcalVetoFromTrueMu(){
  gSystem->Load("libMiniFWLite.so");
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  
  e = new TChain("Events");
  e->Add("/data/tmp/cms2/Zmumu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");

  TCut musIsoCut = "(mus_iso03_sumPt+mus_iso03_emEt+mus_iso03_hadEt)/max(mus_p4.pt(),20)<0.1";
  TCut mcMatchCutWorZ = "(abs(mus_mc_id)==13||(abs(mus_mc_id)==22&&abs(mus_mc_motherid)==13))&&(abs(mus_mc_motherid)==13||abs(mus_mc_motherid)==23||abs(mus_mc_motherid)==24)";
  TCut muQualCut = "(mus_type&2)==2"; // keep this simple
  
  const unsigned int nBinsX = 20;
  double ptBinsD[nBinsX+1] = {10, 12, 14, 16, 18, 20, 23, 26, 29, 33, 37, 42, 48, 54, 60, 70, 80, 100,
			120, 150, 200};
  unsigned int nBinsY = 50;
  double maxY = 15;

  hhcalVmcptZ = new TH2D("hhcalVmcptZ", "hhcalVmcptZ", nBinsX, ptBinsD, nBinsY, 0,maxY);
  cv = new TCanvas("hhcalVmcptZ");
  e->Draw("min(mus_iso_hcalvetoDep,14.99):min(mus_mc_p4.pt(),199)>>hhcalVmcptZ",
	  muQualCut && musIsoCut && mcMatchCutWorZ, "colz");
  std::cout<<"Got "<<hhcalVmcptZ->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hhcalVmcptZ->SetStats(0);
  hhcalVmcptZ->SetTitle("E_{T} in veto cone vs p_{T}^{MC} in Z#rightarrow#mu#mu sample");
  hhcalVmcptZ->Draw("colz");
  gPad->SaveAs("hhcalVmcptZmm.png");
  gPad->SaveAs("hhcalVmcptZmm.eps");

  hhcalVmcptCumulZ = new TH2D("hhcalVmcptCumulZ","#epsilon(E_{T}^{HCAL}<y, p_{T}^{MC})  in Z#rightarrow#mu#mu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hhcalVmcptZ->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hhcalVmcptZ->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hhcalVmcptCumulZ->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hhcalVmcptCumulZ");
  hhcalVmcptCumulZ->SetMaximum(1.0);
  hhcalVmcptCumulZ->SetMinimum(0.9);
  hhcalVmcptCumulZ->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hhcalVmcptCumulZ->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  pls->Draw();
  gPad->SaveAs("hhcalVmcptCumulZmm.png");
  gPad->SaveAs("hhcalVmcptCumulZmm.eps");

  cv = new TCanvas("hhcalVmcptZON");
  hhcalVmcptZON = new TH2D("hhcalVmcptZON", "hhcalVmcptZON", nBinsX, ptBinsD, nBinsY, 0,maxY);
  TCut genpZONCut = "Sum$(genps_id==23&&abs(genps_p4.mass()-91)<15)==1";
  e->Draw("min(mus_iso_hcalvetoDep,14.99):min(mus_mc_p4.pt(),199)>>hhcalVmcptZON",
	  muQualCut && musIsoCut && mcMatchCutWorZ && genpZONCut, "colz");
  std::cout<<"Got "<<hhcalVmcptZON->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hhcalVmcptZON->SetStats(0);
  hhcalVmcptZON->SetTitle("E_{T} in veto cone vs p_{T}^{MC} in ZON#rightarrow#mu#mu sample");
  hhcalVmcptZON->Draw("colz");
  gPad->SaveAs("hhcalVmcptZONmm.png");
  gPad->SaveAs("hhcalVmcptZONmm.eps");

  hhcalVmcptCumulZON = new TH2D("hhcalVmcptCumulZON","#epsilon(E_{T}^{HCAL}<y, p_{T}^{MC})  in ZON#rightarrow#mu#mu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hhcalVmcptZON->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hhcalVmcptZON->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hhcalVmcptCumulZON->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hhcalVmcptCumulZON");
  hhcalVmcptCumulZON->SetMaximum(1.0);
  hhcalVmcptCumulZON->SetMinimum(0.9);
  hhcalVmcptCumulZON->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hhcalVmcptCumulZON->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  pls->Draw();
  gPad->SaveAs("hhcalVmcptCumulZONmm.png");
  gPad->SaveAs("hhcalVmcptCumulZONmm.eps");

  cv = new TCanvas("hhcalVmcptZOFF");
  hhcalVmcptZOFF = new TH2D("hhcalVmcptZOFF", "hhcalVmcptZOFF", nBinsX, ptBinsD, nBinsY, 0,maxY);
  TCut genpZOFFCut = "Sum$(genps_id==23&&abs(genps_p4.mass()-91)>15)==1";
  e->Draw("min(mus_iso_hcalvetoDep,14.99):min(mus_mc_p4.pt(),199)>>hhcalVmcptZOFF",
	  muQualCut && musIsoCut && mcMatchCutWorZ && genpZOFFCut, "colz");
  std::cout<<"Got "<<hhcalVmcptZOFF->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hhcalVmcptZOFF->SetStats(0);
  hhcalVmcptZOFF->SetTitle("E_{T} in veto cone vs p_{T}^{MC} in ZOFF#rightarrow#mu#mu sample");
  hhcalVmcptZOFF->Draw("colz");
  gPad->SaveAs("hhcalVmcptZOFFmm.png");
  gPad->SaveAs("hhcalVmcptZOFFmm.eps");

  hhcalVmcptCumulZOFF = new TH2D("hhcalVmcptCumulZOFF","#epsilon(E_{T}^{HCAL}<y, p_{T}^{MC})  in ZOFF#rightarrow#mu#mu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hhcalVmcptZOFF->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hhcalVmcptZOFF->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hhcalVmcptCumulZOFF->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hhcalVmcptCumulZOFF");
  hhcalVmcptCumulZOFF->SetMaximum(1.0);
  hhcalVmcptCumulZOFF->SetMinimum(0.9);
  hhcalVmcptCumulZOFF->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hhcalVmcptCumulZOFF->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  pls->Draw();
  gPad->SaveAs("hhcalVmcptCumulZOFFmm.png");
  gPad->SaveAs("hhcalVmcptCumulZOFFmm.eps");


  eW = new TChain("Events");
  eW->Add("/data/tmp/cms2/Wmunu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
  hhcalVmcptW = new TH2D("hhcalVmcptW", "hhcalVmcptW", nBinsX, ptBinsD, nBinsY, 0,maxY);
  cv = new TCanvas("hhcalVmcptW");
  eW->Draw("min(mus_iso_hcalvetoDep,14.99):min(mus_mc_p4.pt(),199)>>hhcalVmcptW",
	   muQualCut && musIsoCut && mcMatchCutWorZ, "colz");
  std::cout<<"Got "<<hhcalVmcptW->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hhcalVmcptW->SetStats(0);
  hhcalVmcptW->SetTitle("E_{T} in veto cone vs p_{T}^{MC} in W#rightarrow#mu#nu sample");
  hhcalVmcptW->Draw("colz");
  gPad->SaveAs("hhcalVmcptWm.png");
  gPad->SaveAs("hhcalVmcptWm.eps");

  hhcalVmcptCumulW = new TH2D("hhcalVmcptCumulW","#epsilon(E_{T}^{HCAL}<y, p_{T}^{MC})  in W#rightarrow#mu#nu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hhcalVmcptW->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hhcalVmcptW->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hhcalVmcptCumulW->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hhcalVmcptCumulW");
  hhcalVmcptCumulW->SetMaximum(1.0);
  hhcalVmcptCumulW->SetMinimum(0.9);
  hhcalVmcptCumulW->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hhcalVmcptCumulW->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  pls->Draw();
  gPad->SaveAs("hhcalVmcptCumulWm.png");
  gPad->SaveAs("hhcalVmcptCumulWm.eps");


  eTT = new TChain("Events");
  eTT->Add("/data/tmp/cms2/TTbar_Summer09-MC_31X_V3_7TeV-v1/V03-00-34/merged_ntuple*.root");

  hhcalVmcptTT = new TH2D("hhcalVmcptTT", "hhcalVmcptTT", nBinsX, ptBinsD, nBinsY, 0,maxY);
  cv = new TCanvas("hhcalVmcptTT");
  eTT->Draw("min(mus_iso_hcalvetoDep,14.99):min(mus_mc_p4.pt(),199)>>hhcalVmcptTT",
	   muQualCut && musIsoCut && mcMatchCutWorZ, "colz");
  std::cout<<"Got "<<hhcalVmcptTT->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hhcalVmcptTT->SetStats(0);
  hhcalVmcptTT->SetTitle("E_{T} in veto cone vs p_{T}^{MC} in TT#rightarrow#mu#nu sample");
  hhcalVmcptTT->Draw("colz");
  gPad->SaveAs("hhcalVmcptTTm.png");
  gPad->SaveAs("hhcalVmcptTTm.eps");

  hhcalVmcptCumulTT = new TH2D("hhcalVmcptCumulTT","#epsilon(E_{T}^{HCAL}<y, p_{T}^{MC})  in TT#rightarrow#mu#nu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hhcalVmcptTT->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hhcalVmcptTT->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hhcalVmcptCumulTT->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hhcalVmcptCumulTT");
  hhcalVmcptCumulTT->SetMaximum(1.0);
  hhcalVmcptCumulTT->SetMinimum(0.9);
  hhcalVmcptCumulTT->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hhcalVmcptCumulTT->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  pls->Draw();
  gPad->SaveAs("hhcalVmcptCumulTTm.png");
  gPad->SaveAs("hhcalVmcptCumulTTm.eps");
  //in dileptons now
  hhcalVmcptTTdil = new TH2D("hhcalVmcptTTdil", "hhcalVmcptTTdil", nBinsX, ptBinsD, nBinsY, 0,maxY);
  cv = new TCanvas("hhcalVmcptTTdil");
  TCut genpDilCut = "Sum$(abs(genps_id)==11||abs(genps_id)==13||abs(genps_id)==15)==2";
  eTT->Draw("min(mus_iso_hcalvetoDep,14.99):min(mus_mc_p4.pt(),199)>>hhcalVmcptTTdil",
	    genpDilCut && muQualCut && musIsoCut && mcMatchCutWorZ, "colz");
  std::cout<<"Got "<<hhcalVmcptTTdil->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hhcalVmcptTTdil->SetStats(0);
  hhcalVmcptTTdil->SetTitle("E_{T} in veto cone vs p_{T}^{MC} in TTdil#rightarrow#mu#nu sample");
  hhcalVmcptTTdil->Draw("colz");
  gPad->SaveAs("hhcalVmcptTTdilm.png");
  gPad->SaveAs("hhcalVmcptTTdilm.eps");

  hhcalVmcptCumulTTdil = new TH2D("hhcalVmcptCumulTTdil","#epsilon(E_{T}^{HCAL}<y, p_{T}^{MC})  in TTdil#rightarrow#mu#nu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hhcalVmcptTTdil->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hhcalVmcptTTdil->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hhcalVmcptCumulTTdil->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hhcalVmcptCumulTTdil");
  hhcalVmcptCumulTTdil->SetMaximum(1.0);
  hhcalVmcptCumulTTdil->SetMinimum(0.9);
  hhcalVmcptCumulTTdil->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hhcalVmcptCumulTTdil->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  pls->Draw();
  gPad->SaveAs("hhcalVmcptCumulTTdilm.png");
  gPad->SaveAs("hhcalVmcptCumulTTdilm.eps");

}

void plotEcalVetoFromTrueMu(){
  gSystem->Load("libMiniFWLite.so");
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  
  e = new TChain("Events");
  e->Add("/data/tmp/cms2/Zmumu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");

  TCut musIsoCut = "(mus_iso03_sumPt+mus_iso03_emEt+mus_iso03_hadEt)/max(mus_p4.pt(),20)<0.1";
  TCut mcMatchCutWorZ = "(abs(mus_mc_id)==13||(abs(mus_mc_id)==22&&abs(mus_mc_motherid)==13))&&(abs(mus_mc_motherid)==13||abs(mus_mc_motherid)==23||abs(mus_mc_motherid)==24)";
  TCut muQualCut = "(mus_type&2)==2"; // keep this simple
  
  const unsigned int nBinsX = 20;
  double ptBinsD[nBinsX+1] = {10, 12, 14, 16, 18, 20, 23, 26, 29, 33, 37, 42, 48, 54, 60, 70, 80, 100,
			120, 150, 200};
  unsigned int nBinsY = 50;
  double maxY = 15;

  hecalVmcptZ = new TH2D("hecalVmcptZ", "hecalVmcptZ", nBinsX, ptBinsD, nBinsY, 0,maxY);
  e->Draw("min(mus_iso_ecalvetoDep,14.99):min(mus_mc_p4.pt(),199)>>hecalVmcptZ",
	  muQualCut && musIsoCut && mcMatchCutWorZ, "colz");
  std::cout<<"Got "<<hecalVmcptZ->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hecalVmcptZ->SetStats(0);
  hecalVmcptZ->SetTitle("E_{T} in veto cone vs p_{T}^{MC} in Z#rightarrow#mu#mu sample");
  hecalVmcptZ->Draw("colz");
  gPad->SaveAs("hecalVmcptZmm.png");
  gPad->SaveAs("hecalVmcptZmm.eps");

  hecalVmcptCumulZ = new TH2D("hecalVmcptCumulZ","#epsilon(E_{T}^{ECAL}<y, p_{T}^{MC})  in Z#rightarrow#mu#mu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hecalVmcptZ->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hecalVmcptZ->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hecalVmcptCumulZ->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas();
  hecalVmcptCumulZ->SetMaximum(1.0);
  hecalVmcptCumulZ->SetMinimum(0.9);
  hecalVmcptCumulZ->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hecalVmcptCumulZ->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  pls->Draw();
  gPad->SaveAs("hecalVmcptCumulZmm.png");
  gPad->SaveAs("hecalVmcptCumulZmm.eps");

  cv = new TCanvas("hecalVmcptZON");
  hecalVmcptZON = new TH2D("hecalVmcptZON", "hecalVmcptZON", nBinsX, ptBinsD, nBinsY, 0,maxY);
  TCut genpZONCut = "Sum$(genps_id==23&&abs(genps_p4.mass()-91)<15)==1";
  e->Draw("min(mus_iso_ecalvetoDep,14.99):min(mus_mc_p4.pt(),199)>>hecalVmcptZON",
	  muQualCut && musIsoCut && mcMatchCutWorZ && genpZONCut, "colz");
  std::cout<<"Got "<<hecalVmcptZON->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hecalVmcptZON->SetStats(0);
  hecalVmcptZON->SetTitle("E_{T} in veto cone vs p_{T}^{MC} in ZON#rightarrow#mu#mu sample");
  hecalVmcptZON->Draw("colz");
  gPad->SaveAs("hecalVmcptZONmm.png");
  gPad->SaveAs("hecalVmcptZONmm.eps");

  hecalVmcptCumulZON = new TH2D("hecalVmcptCumulZON","#epsilon(E_{T}^{ECAL}<y, p_{T}^{MC})  in ZON#rightarrow#mu#mu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hecalVmcptZON->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hecalVmcptZON->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hecalVmcptCumulZON->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hecalVmcptCumulZON");
  hecalVmcptCumulZON->SetMaximum(1.0);
  hecalVmcptCumulZON->SetMinimum(0.9);
  hecalVmcptCumulZON->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hecalVmcptCumulZON->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  pls->Draw();
  gPad->SaveAs("hecalVmcptCumulZONmm.png");
  gPad->SaveAs("hecalVmcptCumulZONmm.eps");

  cv = new TCanvas("hecalVmcptZOFF");
  hecalVmcptZOFF = new TH2D("hecalVmcptZOFF", "hecalVmcptZOFF", nBinsX, ptBinsD, nBinsY, 0,maxY);
  TCut genpZOFFCut = "Sum$(genps_id==23&&abs(genps_p4.mass()-91)>15)==1";
  e->Draw("min(mus_iso_ecalvetoDep,14.99):min(mus_mc_p4.pt(),199)>>hecalVmcptZOFF",
	  muQualCut && musIsoCut && mcMatchCutWorZ && genpZOFFCut, "colz");
  std::cout<<"Got "<<hecalVmcptZOFF->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hecalVmcptZOFF->SetStats(0);
  hecalVmcptZOFF->SetTitle("E_{T} in veto cone vs p_{T}^{MC} in ZOFF#rightarrow#mu#mu sample");
  hecalVmcptZOFF->Draw("colz");
  gPad->SaveAs("hecalVmcptZOFFmm.png");
  gPad->SaveAs("hecalVmcptZOFFmm.eps");

  hecalVmcptCumulZOFF = new TH2D("hecalVmcptCumulZOFF","#epsilon(E_{T}^{ECAL}<y, p_{T}^{MC})  in ZOFF#rightarrow#mu#mu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hecalVmcptZOFF->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hecalVmcptZOFF->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hecalVmcptCumulZOFF->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hecalVmcptCumulZOFF");
  hecalVmcptCumulZOFF->SetMaximum(1.0);
  hecalVmcptCumulZOFF->SetMinimum(0.9);
  hecalVmcptCumulZOFF->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hecalVmcptCumulZOFF->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  pls->Draw();
  gPad->SaveAs("hecalVmcptCumulZOFFmm.png");
  gPad->SaveAs("hecalVmcptCumulZOFFmm.eps");


  eW = new TChain("Events");
  eW->Add("/data/tmp/cms2/Wmunu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
  hecalVmcptW = new TH2D("hecalVmcptW", "hecalVmcptW", nBinsX, ptBinsD, nBinsY, 0,maxY);
  cv = new TCanvas("hecalVmcptW");
  eW->Draw("min(mus_iso_ecalvetoDep,14.99):min(mus_mc_p4.pt(),199)>>hecalVmcptW",
	   muQualCut && musIsoCut && mcMatchCutWorZ, "colz");
  std::cout<<"Got "<<hecalVmcptW->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hecalVmcptW->SetStats(0);
  hecalVmcptW->SetTitle("E_{T} in veto cone vs p_{T}^{MC} in W#rightarrow#mu#nu sample");
  hecalVmcptW->Draw("colz");
  gPad->SaveAs("hecalVmcptWm.png");
  gPad->SaveAs("hecalVmcptWm.eps");

  hecalVmcptCumulW = new TH2D("hecalVmcptCumulW","#epsilon(E_{T}^{ECAL}<y, p_{T}^{MC})  in W#rightarrow#mu#nu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hecalVmcptW->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hecalVmcptW->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hecalVmcptCumulW->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hecalVmcptCumulW");
  hecalVmcptCumulW->SetMaximum(1.0);
  hecalVmcptCumulW->SetMinimum(0.9);
  hecalVmcptCumulW->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hecalVmcptCumulW->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  pls->Draw();
  gPad->SaveAs("hecalVmcptCumulWm.png");
  gPad->SaveAs("hecalVmcptCumulWm.eps");


  eTT = new TChain("Events");
  eTT->Add("/data/tmp/cms2/TTbar_Summer09-MC_31X_V3_7TeV-v1/V03-00-34/merged_ntuple*.root");

  hecalVmcptTT = new TH2D("hecalVmcptTT", "hecalVmcptTT", nBinsX, ptBinsD, nBinsY, 0,maxY);
  cv = new TCanvas("hecalVmcptTT");
  eTT->Draw("min(mus_iso_ecalvetoDep,14.99):min(mus_mc_p4.pt(),199)>>hecalVmcptTT",
	   muQualCut && musIsoCut && mcMatchCutWorZ, "colz");
  std::cout<<"Got "<<hecalVmcptTT->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hecalVmcptTT->SetStats(0);
  hecalVmcptTT->SetTitle("E_{T} in veto cone vs p_{T}^{MC} in TT#rightarrow#mu#nu sample");
  hecalVmcptTT->Draw("colz");
  gPad->SaveAs("hecalVmcptTTm.png");
  gPad->SaveAs("hecalVmcptTTm.eps");

  hecalVmcptCumulTT = new TH2D("hecalVmcptCumulTT","#epsilon(E_{T}^{ECAL}<y, p_{T}^{MC})  in TT#rightarrow#mu#nu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hecalVmcptTT->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hecalVmcptTT->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hecalVmcptCumulTT->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hecalVmcptCumulTT");
  hecalVmcptCumulTT->SetMaximum(1.0);
  hecalVmcptCumulTT->SetMinimum(0.9);
  hecalVmcptCumulTT->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hecalVmcptCumulTT->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  pls->Draw();
  gPad->SaveAs("hecalVmcptCumulTTm.png");
  gPad->SaveAs("hecalVmcptCumulTTm.eps");
  //in dileptons now
  hecalVmcptTTdil = new TH2D("hecalVmcptTTdil", "hecalVmcptTTdil", nBinsX, ptBinsD, nBinsY, 0,maxY);
  cv = new TCanvas("hecalVmcptTTdil");
  TCut genpDilCut = "Sum$(abs(genps_id)==11||abs(genps_id)==13||abs(genps_id)==15)==2";
  eTT->Draw("min(mus_iso_ecalvetoDep,14.99):min(mus_mc_p4.pt(),199)>>hecalVmcptTTdil",
	    genpDilCut && muQualCut && musIsoCut && mcMatchCutWorZ, "colz");
  std::cout<<"Got "<<hecalVmcptTTdil->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hecalVmcptTTdil->SetStats(0);
  hecalVmcptTTdil->SetTitle("E_{T} in veto cone vs p_{T}^{MC} in TTdil#rightarrow#mu#nu sample");
  hecalVmcptTTdil->Draw("colz");
  gPad->SaveAs("hecalVmcptTTdilm.png");
  gPad->SaveAs("hecalVmcptTTdilm.eps");

  hecalVmcptCumulTTdil = new TH2D("hecalVmcptCumulTTdil","#epsilon(E_{T}^{ECAL}<y, p_{T}^{MC})  in TTdil#rightarrow#mu#nu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hecalVmcptTTdil->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hecalVmcptTTdil->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hecalVmcptCumulTTdil->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hecalVmcptCumulTTdil");
  hecalVmcptCumulTTdil->SetMaximum(1.0);
  hecalVmcptCumulTTdil->SetMinimum(0.9);
  hecalVmcptCumulTTdil->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hecalVmcptCumulTTdil->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  pls->Draw();
  gPad->SaveAs("hecalVmcptCumulTTdilm.png");
  gPad->SaveAs("hecalVmcptCumulTTdilm.eps");
}

void plotCombIsoFromTrueMu(){
  gSystem->Load("libMiniFWLite.so");
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  
  e = new TChain("Events");
  e->Add("/data/tmp/cms2/Zmumu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");

  TCut musIsoCut = "1<2";
  TCut mcMatchCutWorZ = "(abs(mus_mc_id)==13||(abs(mus_mc_id)==22&&abs(mus_mc_motherid)==13))&&(abs(mus_mc_motherid)==13||abs(mus_mc_motherid)==23||abs(mus_mc_motherid)==24)";
  TCut muQualCut = "(mus_type&2)==2&&(mus_type&4)==4&&abs(mus_d0corr)<0.02"; // keep this simple
  
  const unsigned int nBinsX = 20;
  double ptBinsD[nBinsX+1] = {10, 12, 14, 16, 18, 20, 23, 26, 29, 33, 37, 42, 48, 54, 60, 70, 80, 100,
			120, 150, 200};
  unsigned int nBinsY = 50;
  double maxY = 0.5;

  hisoVmcptZ = new TH2D("hisoVmcptZ", "hisoVmcptZ", nBinsX, ptBinsD, nBinsY, 0,maxY);
  e->Draw("min((mus_iso03_sumPt+mus_iso03_emEt+mus_iso03_hadEt)/max(mus_p4.pt(),20),0.499):min(mus_mc_p4.pt(),199)>>hisoVmcptZ",
	  muQualCut && musIsoCut && mcMatchCutWorZ, "colz");
  std::cout<<"Got "<<hisoVmcptZ->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hisoVmcptZ->SetStats(0);
  hisoVmcptZ->SetTitle("sum(Iso)/pt_{reco} vs p_{T}^{MC} in Z#rightarrow#mu#mu sample");
  hisoVmcptZ->Draw("colz");
  gPad->SaveAs("hisoVmcptZmm.png");
  gPad->SaveAs("hisoVmcptZmm.eps");

  hisoVmcptCumulZ = new TH2D("hisoVmcptCumulZ","#epsilon(sum(Iso)/pt_{reco}<y, p_{T}^{MC})  in Z#rightarrow#mu#mu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hisoVmcptZ->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hisoVmcptZ->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hisoVmcptCumulZ->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas();
  hisoVmcptCumulZ->SetMaximum(1.0);
  hisoVmcptCumulZ->SetMinimum(0.8);
  hisoVmcptCumulZ->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hisoVmcptCumulZ->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  // pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  // pls->Draw();
  gPad->SaveAs("hisoVmcptCumulZmm.png");
  gPad->SaveAs("hisoVmcptCumulZmm.eps");

  cv = new TCanvas("hisoVmcptZON");
  hisoVmcptZON = new TH2D("hisoVmcptZON", "hisoVmcptZON", nBinsX, ptBinsD, nBinsY, 0,maxY);
  TCut genpZONCut = "Sum$(genps_id==23&&abs(genps_p4.mass()-91)<15)==1";
  e->Draw("min((mus_iso03_sumPt+mus_iso03_emEt+mus_iso03_hadEt)/max(mus_p4.pt(),20),0.499):min(mus_mc_p4.pt(),199)>>hisoVmcptZON",
	  muQualCut && musIsoCut && mcMatchCutWorZ && genpZONCut, "colz");
  std::cout<<"Got "<<hisoVmcptZON->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hisoVmcptZON->SetStats(0);
  hisoVmcptZON->SetTitle("sum(Iso)/pt_{reco} vs p_{T}^{MC} in ZON#rightarrow#mu#mu sample");
  hisoVmcptZON->Draw("colz");
  gPad->SaveAs("hisoVmcptZONmm.png");
  gPad->SaveAs("hisoVmcptZONmm.eps");

  hisoVmcptCumulZON = new TH2D("hisoVmcptCumulZON","#epsilon(sum(Iso)/pt_{reco}<y, p_{T}^{MC})  in ZON#rightarrow#mu#mu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hisoVmcptZON->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hisoVmcptZON->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hisoVmcptCumulZON->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hisoVmcptCumulZON");
  hisoVmcptCumulZON->SetMaximum(1.0);
  hisoVmcptCumulZON->SetMinimum(0.8);
  hisoVmcptCumulZON->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hisoVmcptCumulZON->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  // pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  // pls->Draw();
  gPad->SaveAs("hisoVmcptCumulZONmm.png");
  gPad->SaveAs("hisoVmcptCumulZONmm.eps");

  cv = new TCanvas("hisoVmcptZOFF");
  hisoVmcptZOFF = new TH2D("hisoVmcptZOFF", "hisoVmcptZOFF", nBinsX, ptBinsD, nBinsY, 0,maxY);
  TCut genpZOFFCut = "Sum$(genps_id==23&&abs(genps_p4.mass()-91)>15)==1";
  e->Draw("min((mus_iso03_sumPt+mus_iso03_emEt+mus_iso03_hadEt)/max(mus_p4.pt(),20),0.499):min(mus_mc_p4.pt(),199)>>hisoVmcptZOFF",
	  muQualCut && musIsoCut && mcMatchCutWorZ && genpZOFFCut, "colz");
  std::cout<<"Got "<<hisoVmcptZOFF->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hisoVmcptZOFF->SetStats(0);
  hisoVmcptZOFF->SetTitle("sum(Iso)/pt_{reco} vs p_{T}^{MC} in ZOFF#rightarrow#mu#mu sample");
  hisoVmcptZOFF->Draw("colz");
  gPad->SaveAs("hisoVmcptZOFFmm.png");
  gPad->SaveAs("hisoVmcptZOFFmm.eps");

  hisoVmcptCumulZOFF = new TH2D("hisoVmcptCumulZOFF","#epsilon(sum(Iso)/pt_{reco}<y, p_{T}^{MC})  in ZOFF#rightarrow#mu#mu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hisoVmcptZOFF->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hisoVmcptZOFF->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hisoVmcptCumulZOFF->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hisoVmcptCumulZOFF");
  hisoVmcptCumulZOFF->SetMaximum(1.0);
  hisoVmcptCumulZOFF->SetMinimum(0.8);
  hisoVmcptCumulZOFF->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hisoVmcptCumulZOFF->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  // pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  // pls->Draw();
  gPad->SaveAs("hisoVmcptCumulZOFFmm.png");
  gPad->SaveAs("hisoVmcptCumulZOFFmm.eps");


  eW = new TChain("Events");
  eW->Add("/data/tmp/cms2/Wmunu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
  hisoVmcptW = new TH2D("hisoVmcptW", "hisoVmcptW", nBinsX, ptBinsD, nBinsY, 0,maxY);
  cv = new TCanvas("hisoVmcptW");
  eW->Draw("min((mus_iso03_sumPt+mus_iso03_emEt+mus_iso03_hadEt)/max(mus_p4.pt(),20),0.499):min(mus_mc_p4.pt(),199)>>hisoVmcptW",
	   muQualCut && musIsoCut && mcMatchCutWorZ, "colz");
  std::cout<<"Got "<<hisoVmcptW->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hisoVmcptW->SetStats(0);
  hisoVmcptW->SetTitle("sum(Iso)/pt_{reco} vs p_{T}^{MC} in W#rightarrow#mu#nu sample");
  hisoVmcptW->Draw("colz");
  gPad->SaveAs("hisoVmcptWm.png");
  gPad->SaveAs("hisoVmcptWm.eps");

  hisoVmcptCumulW = new TH2D("hisoVmcptCumulW","#epsilon(sum(Iso)/pt_{reco}<y, p_{T}^{MC})  in W#rightarrow#mu#nu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hisoVmcptW->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hisoVmcptW->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hisoVmcptCumulW->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hisoVmcptCumulW");
  hisoVmcptCumulW->SetMaximum(1.0);
  hisoVmcptCumulW->SetMinimum(0.8);
  hisoVmcptCumulW->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hisoVmcptCumulW->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  // pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  // pls->Draw();
  gPad->SaveAs("hisoVmcptCumulWm.png");
  gPad->SaveAs("hisoVmcptCumulWm.eps");


  eTT = new TChain("Events");
  eTT->Add("/data/tmp/cms2/TTbar_Summer09-MC_31X_V3_7TeV-v1/V03-00-34/merged_ntuple*.root");

  hisoVmcptTT = new TH2D("hisoVmcptTT", "hisoVmcptTT", nBinsX, ptBinsD, nBinsY, 0,maxY);
  cv = new TCanvas("hisoVmcptTT");
  eTT->Draw("min((mus_iso03_sumPt+mus_iso03_emEt+mus_iso03_hadEt)/max(mus_p4.pt(),20),0.499):min(mus_mc_p4.pt(),199)>>hisoVmcptTT",
	   muQualCut && musIsoCut && mcMatchCutWorZ, "colz");
  std::cout<<"Got "<<hisoVmcptTT->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hisoVmcptTT->SetStats(0);
  hisoVmcptTT->SetTitle("sum(Iso)/pt_{reco} vs p_{T}^{MC} in TT#rightarrow#mu#nu sample");
  hisoVmcptTT->Draw("colz");
  gPad->SaveAs("hisoVmcptTTm.png");
  gPad->SaveAs("hisoVmcptTTm.eps");

  hisoVmcptCumulTT = new TH2D("hisoVmcptCumulTT","#epsilon(sum(Iso)/pt_{reco}<y, p_{T}^{MC})  in TT#rightarrow#mu#nu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hisoVmcptTT->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hisoVmcptTT->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hisoVmcptCumulTT->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hisoVmcptCumulTT");
  hisoVmcptCumulTT->SetMaximum(1.0);
  hisoVmcptCumulTT->SetMinimum(0.8);
  hisoVmcptCumulTT->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hisoVmcptCumulTT->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  // pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  // pls->Draw();
  gPad->SaveAs("hisoVmcptCumulTTm.png");
  gPad->SaveAs("hisoVmcptCumulTTm.eps");
  //in dileptons now
  hisoVmcptTTdil = new TH2D("hisoVmcptTTdil", "hisoVmcptTTdil", nBinsX, ptBinsD, nBinsY, 0,maxY);
  cv = new TCanvas("hisoVmcptTTdil");
  TCut genpDilCut = "Sum$(abs(genps_id)==11||abs(genps_id)==13||abs(genps_id)==15)==2";
  eTT->Draw("min((mus_iso03_sumPt+mus_iso03_emEt+mus_iso03_hadEt)/max(mus_p4.pt(),20),0.499):min(mus_mc_p4.pt(),199)>>hisoVmcptTTdil",
	    genpDilCut && muQualCut && musIsoCut && mcMatchCutWorZ, "colz");
  std::cout<<"Got "<<hisoVmcptTTdil->GetEntries()<<" entries passing the cuts"<<std::endl;

  gPad->SetLogz();
  hisoVmcptTTdil->SetStats(0);
  hisoVmcptTTdil->SetTitle("sum(Iso)/pt_{reco} vs p_{T}^{MC} in TTdil#rightarrow#mu#nu sample");
  hisoVmcptTTdil->Draw("colz");
  gPad->SaveAs("hisoVmcptTTdilm.png");
  gPad->SaveAs("hisoVmcptTTdilm.eps");

  hisoVmcptCumulTTdil = new TH2D("hisoVmcptCumulTTdil","#epsilon(sum(Iso)/pt_{reco}<y, p_{T}^{MC})  in TTdil#rightarrow#mu#nu sample",
			     nBinsX, ptBinsD,nBinsY,0,maxY);
  for (unsigned int iX=0;iX<20+2;++iX){
    for (unsigned int iY=0;iY<nBinsY+2;++iY){
      double sumX = hisoVmcptTTdil->Integral(iX,iX, 0,nBinsY+2);      
      double sumXiY = hisoVmcptTTdil->Integral(iX,iX, 0,iY);
      double rat = sumX > 0 ? sumXiY/sumX : 0;
      hisoVmcptCumulTTdil->SetBinContent(iX,iY, rat);
    }    
  }
  cv = new TCanvas("hisoVmcptCumulTTdil");
  hisoVmcptCumulTTdil->SetMaximum(1.0);
  hisoVmcptCumulTTdil->SetMinimum(0.8);
  hisoVmcptCumulTTdil->SetStats(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hisoVmcptCumulTTdil->Draw("colz");
  TLine* pl = new TLine(10,6,30,6);
  pl->SetLineWidth(2);
  // pl->Draw();
  TLine* pls = new TLine(30,6,120,15);
  pls->SetLineWidth(2);
  // pls->Draw();
  gPad->SaveAs("hisoVmcptCumulTTdilm.png");
  gPad->SaveAs("hisoVmcptCumulTTdilm.eps");
}
