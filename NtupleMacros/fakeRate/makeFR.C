{
  cout << "Be careful about abs vs fabs" << endl;
  gROOT->LoadMacro("eff2.C");
  gROOT->LoadMacro("histio.cc");
  //gStyle->SetOptStat(0);

  TChain *ch1 = new TChain("tree");
  TChain *ch2 = new TChain("tree");
  ch1->Add("JMTMonitor.root");
  ch2->Add("./august9/JMT_600nb.root");

  // The bins for the electron FR
  double ptEleBin[4] ={10.,15.,20.,40.};
  double etaEleBin[5]={0.,1., 1.479, 2., 2.5}

  // Book the electron FR
  TH2F* yyy = new TH2F("yyyy","yyyy",4,etaEleBin,3,ptEleBin);
  TH2F* enum6u = new TH2F("enum6u","enum6u",4,etaEleBin,3,ptEleBin);
  TH2F* v16u  = new TH2F("v16u", "v16u", 4,etaEleBin,3,ptEleBin);
  TH2F* v26u  = new TH2F("v26u", "v26u", 4,etaEleBin,3,ptEleBin);
  TH2F* v36u  = new TH2F("v36u", "v36u", 4,etaEleBin,3,ptEleBin);

  TH2F* enum10u = new TH2F("enum10u","enum10u",4,etaEleBin,3,ptEleBin);
  TH2F* v110u  = new TH2F("v110u", "v110u", 4,etaEleBin,3,ptEleBin);
  TH2F* v210u  = new TH2F("v210u", "v210u", 4,etaEleBin,3,ptEleBin);
  TH2F* v310u  = new TH2F("v310u", "v310u", 4,etaEleBin,3,ptEleBin);

  TH2F* enum15u = new TH2F("enum15u","enum15u",4,etaEleBin,3,ptEleBin);
  TH2F* v115u  = new TH2F("v115u", "v115u", 4,etaEleBin,3,ptEleBin);
  TH2F* v215u  = new TH2F("v215u", "v215u", 4,etaEleBin,3,ptEleBin);
  TH2F* v315u  = new TH2F("v315u", "v315u", 4,etaEleBin,3,ptEleBin);

  // The bins for the muon FR
  double ptMuBin[3] ={10.,15.,40.};
  double etaMuBin[5]={0., 1., 1.5, 2., 2.5};
  
  // Book the muon FR
  TH2F* mnum6u = new TH2F("mnum6u","num6u",4,etaMuBin,2,ptMuBin);
  TH2F* fo6u  = new TH2F("fo6u", "fo6u", 4,etaMuBin,2,ptMuBin);

  TH2F* mnum10u = new TH2F("mnum10u","num10u",4,etaMuBin,2,ptMuBin);
  TH2F* fo10u  = new TH2F("fo10u", "fo10u", 4,etaMuBin,2,ptMuBin);

  TH2F* mnum15u = new TH2F("mnum15u","num15u",4,etaMuBin,2,ptMuBin);
  TH2F* fo15u  = new TH2F("fo15u", "fo15u", 4,etaMuBin,2,ptMuBin);


  // Fill the muon histograms
  // 6u
  ch1->Draw("min(pt,39.99):abs(eta)>>mnum6u","numAug9&&l16u>1&&abs(id)==13&&tcmet<20");
  ch1->Draw("min(pt,39.99):abs(eta)>>fo6u",      "l16u>1&&abs(id)==13&&tcmet<20");
  // 10u
  ch1->Draw("min(pt,39.99):abs(eta)>>mnum10u","numAug9&&l110u>1&&abs(id)==13&&tcmet<20");
  ch1->Draw("min(pt,39.99):abs(eta)>>fo10u",      "l110u>1&&abs(id)==13&&tcmet<20");
  // 15u
  ch2->Draw("min(pt,39.99):abs(eta)>>mnum15u","numAug9&&hlt15u>1&&abs(id)==13&&tcmet<20");
  ch2->Draw("min(pt,39.99):abs(eta)>>fo15u",      "hlt15u>1&&abs(id)==13&&tcmet<20");


  // Get the muon FR
  TH2F* muFR6u  = eff2(fo6u, mnum6u, "muFR6u");
  TH2F* muFR10u = eff2(fo10u,mnum10u,"muFR10u");
  TH2F* muFR15u = eff2(fo15u,mnum15u,"muFR15u");


  // Fill the electron histograms
  // 6u
  ch1->Draw("min(pt,39.99):abs(eta)>>enum6u","numAug9&&l16u>1&&abs(id)==11&&tcmet<20");
  ch1->Draw("min(pt,39.99):abs(eta)>>v16u","v1Aug9&&l16u>1&&abs(id)==11&&tcmet<20");
  ch1->Draw("min(pt,39.99):abs(eta)>>v26u","v2Aug9&&l16u>1&&abs(id)==11&&tcmet<20");
  ch1->Draw("min(pt,39.99):abs(eta)>>v36u","v3Aug9&&l16u>1&&abs(id)==11&&tcmet<20");
  // 10u
  ch1->Draw("min(pt,39.99):abs(eta)>>enum10u","numAug9&&l110u>1&&abs(id)==11&&tcmet<20");
  ch1->Draw("min(pt,39.99):abs(eta)>>v110u","v1Aug9&&l110u>1&&abs(id)==11&&tcmet<20");
  ch1->Draw("min(pt,39.99):abs(eta)>>v210u","v2Aug9&&l110u>1&&abs(id)==11&&tcmet<20");
  ch1->Draw("min(pt,39.99):abs(eta)>>v310u","v3Aug9&&l110u>1&&abs(id)==11&&tcmet<20");
  // 15u
  ch2->Draw("min(pt,39.99):abs(eta)>>enum15u","numAug9&&hlt15u>1&&abs(id)==11&&tcmet<20");
  ch2->Draw("min(pt,39.99):abs(eta)>>v115u","v1Aug9&&hlt15u>1&&abs(id)==11&&tcmet<20");
  ch2->Draw("min(pt,39.99):abs(eta)>>v215u","v2Aug9&&hlt15u>1&&abs(id)==11&&tcmet<20");
  ch2->Draw("min(pt,39.99):abs(eta)>>v315u","v3Aug9&&hlt15u>1&&abs(id)==11&&tcmet<20");


  // Get the electron FR
  // 6u
  TH2F* eFRv16u = eff2(v16u,enum6u,"eFRv16u");
  TH2F* eFRv26u = eff2(v26u,enum6u,"eFRv26u");
  TH2F* eFRv36u = eff2(v36u,enum6u,"eFRv36u");
  //10u
  TH2F* eFRv110u = eff2(v110u,enum10u,"eFRv110u");
  TH2F* eFRv210u = eff2(v210u,enum10u,"eFRv210u");
  TH2F* eFRv310u = eff2(v310u,enum10u,"eFRv310u");
  //15u
  TH2F* eFRv115u = eff2(v115u,enum15u,"eFRv115u");
  TH2F* eFRv215u = eff2(v215u,enum15u,"eFRv215u");
  TH2F* eFRv315u = eff2(v315u,enum15u,"eFRv315u");
}
