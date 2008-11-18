{
// Load and compile something to allow proper treatment of vectors
// Not clear that it is needed
//gROOT->LoadMacro("loader.C+");

// Load various tools
gROOT->ProcessLine(".x ../Tools/setup.C");

//hist file:
TFile *ftt = TFile::Open("processed_data_tag.root");

 float DYee[4], DYmm[4], DYtt[4], tt[4], wjets[4], wz[4], zz[4], ww[4], h140[4], tw[4];

for ( int njet=1; njet<4; ++njet) {// actual njet is one less than this counter!
  for (int i=0; i<4; i++){
    DYee[i] = 0.0;
    DYmm[i] = 0.0;
    DYtt[i] = 0.0;
    tt[i] = 0.0;
    wjets[i] = 0.0;
    wz[i] = 0.0;
    zz[i] = 0.0;
    ww[i] = 0.0;
    h140[i] = 0.0;
    tw[i] = 0.0;
  }

  //do the emu row:
  DYee[0] = dyee_hnJetLepVeto_em->GetBinContent(njet);
  DYmm[0] = dymm_hnJetLepVeto_em->GetBinContent(njet);
  DYtt[0] = dytt_hnJetLepVeto_em->GetBinContent(njet);
  tt[0] = ttbar_hnJetLepVeto_em->GetBinContent(njet);
  wjets[0] = wjets_hnJetLepVeto_em->GetBinContent(njet);
  wz[0] = wz_hnJetLepVeto_em->GetBinContent(njet);
  zz[0] = zz_hnJetLepVeto_em->GetBinContent(njet);
  ww[0] = ww_hnJetLepVeto_em->GetBinContent(njet);
  // h140[0] = H140_hnJetLepVeto_em->GetBinContent(njet);
  tw[0] = tw_hnJetLepVeto_em->GetBinContent(njet);
  
  //do the mumu row:
  DYee[1] = dyee_hnJetLepVeto_mm->GetBinContent(njet);
  DYmm[1] = dymm_hnJetLepVeto_mm->GetBinContent(njet);
  DYtt[1] = dytt_hnJetLepVeto_mm->GetBinContent(njet);
  tt[1] = ttbar_hnJetLepVeto_mm->GetBinContent(njet);
  wjets[1] = wjets_hnJetLepVeto_mm->GetBinContent(njet);
  wz[1] = wz_hnJetLepVeto_mm->GetBinContent(njet);
  zz[1] = zz_hnJetLepVeto_mm->GetBinContent(njet);
  ww[1] = ww_hnJetLepVeto_mm->GetBinContent(njet);
  // h140[1] = H140_hnJetLepVeto_mm->GetBinContent(njet);
  tw[1] = tw_hnJetLepVeto_mm->GetBinContent(njet);
  
  //do the ee row:
  DYee[2] = dyee_hnJetLepVeto_ee->GetBinContent(njet);
  DYmm[2] = dymm_hnJetLepVeto_ee->GetBinContent(njet);
  DYtt[2] = dytt_hnJetLepVeto_ee->GetBinContent(njet);
  tt[2] = ttbar_hnJetLepVeto_ee->GetBinContent(njet);
  wjets[2] = wjets_hnJetLepVeto_ee->GetBinContent(njet);
  wz[2] = wz_hnJetLepVeto_ee->GetBinContent(njet);
  zz[2] = zz_hnJetLepVeto_ee->GetBinContent(njet);
  ww[2] = ww_hnJetLepVeto_ee->GetBinContent(njet);
  // h140[2] = H140_hnJetLepVeto_ee->GetBinContent(njet);
  tw[2] = tw_hnJetLepVeto_ee->GetBinContent(njet);
  
  char* finalState[4];
  finalState[0] = "emu";
  finalState[1] = "mumu";
  finalState[2] = "ee";
  finalState[3] = "total";
  
  cout << "  " << endl;
  cout << " Njet = " << njet-1 << endl;
  cout << "| | DY ee| DY mumu | DY tautau | ttbar | Wjets | WZ | ZZ | WW | TW | " << endl;

  for (int i=0; i<4; i++){
    
    cout << "| " << finalState[i] << " | " << DYee[i] << " | " << DYmm[i] << " | " << DYtt[i] << " | " << tt[i] << " | " << wjets[i] << " | " << wz[i] << " | " << zz[i] << " | " << ww[i] << " | " << tw[i] << " | " << endl;

    DYee[3] = DYee[3]+DYee[i];
    DYmm[3] = DYmm[3]+DYmm[i];
    DYtt[3] = DYtt[3]+DYtt[i];
    tt[3] = tt[3]+tt[i];
    wjets[3] = wjets[3]+wjets[i];
    wz[3] = wz[3]+wz[i];
    zz[3] = zz[3]+zz[i];
    ww[3] = ww[3]+ww[i];
    //   h140[3] = h140[3]+h140[i];
    tw[3] = tw[3]+tw[i];
    
  }
}//end of for loop over njet
}
