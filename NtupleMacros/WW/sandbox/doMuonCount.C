{
// Load and compile something to allow proper treatment of vectors
// Not clear that it is needed
//gROOT->LoadMacro("loader.C+");

// Load various tools
gROOT->ProcessLine(".x ../Tools/setup.C");

//hist file:
TFile *ftt = TFile::Open("myHist_WW_fast_bjetstudy_allmcs_alljets.root");

 float DYee[4], DYmm[4], DYtt[4], tt[4], wjets[4], wz[4], zz[4], ww[4], tw[4], tot[4];

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
    tw[i] = 0.0;
    tot[i] = 0.0;
  }

  //do the total:
  DYee[0] = dyee_extramuonsvsnjet_all->Integral(0,10,njet,njet);
  DYmm[0] = dymm_extramuonsvsnjet_all->Integral(0,10,njet,njet);
  DYtt[0] = dytt_extramuonsvsnjet_all->Integral(0,10,njet,njet);
  tt[0] = ttbar_extramuonsvsnjet_all->Integral(0,10,njet,njet);
  wjets[0] = wjets_extramuonsvsnjet_all->Integral(0,10,njet,njet);
  wz[0] = wz_extramuonsvsnjet_all->Integral(0,10,njet,njet);
  zz[0] = zz_extramuonsvsnjet_all->Integral(0,10,njet,njet);
  ww[0] = ww_extramuonsvsnjet_all->Integral(0,10,njet,njet);
  tw[0] = tw_extramuonsvsnjet_all->Integral(0,10,njet,njet);
  tot[0] = DYee[0]+DYee[0]+DYee[0]+tt[0]+wjets[0]+wz[0]+zz[0]+ww[0]+tw[0];
  
  //do the muon tags:
  DYee[1] = dyee_extramuonsvsnjet_all->Integral(2,10,njet,njet);
  DYmm[1] = dymm_extramuonsvsnjet_all->Integral(2,10,njet,njet);
  DYtt[1] = dytt_extramuonsvsnjet_all->Integral(2,10,njet,njet);
  tt[1] = ttbar_extramuonsvsnjet_all->Integral(2,10,njet,njet);
  wjets[1] = wjets_extramuonsvsnjet_all->Integral(2,10,njet,njet);
  wz[1] = wz_extramuonsvsnjet_all->Integral(2,10,njet,njet);
  zz[1] = zz_extramuonsvsnjet_all->Integral(2,10,njet,njet);
  ww[1] = ww_extramuonsvsnjet_all->Integral(2,10,njet,njet);
  tw[1] = tw_extramuonsvsnjet_all->Integral(2,10,njet,njet);
  tot[1] = DYee[1]+DYee[1]+DYee[1]+tt[1]+wjets[1]+wz[1]+zz[1]+ww[1]+tw[1];
  
  char* finalState[4];
  finalState[0] = "Njet = total";
  finalState[1] = "Njet = muon tag";
  finalState[3] = "Njet = tag prob";
  
  cout << "  " << endl;
  cout << " Njet = " << njet-1 << endl;
  cout << "| | DY ee| DY mumu | DY tautau | ttbar | Wjets | WZ | ZZ | WW | TW | total | " << endl;

  for (int i=0; i<2; i++){
    
    cout << "| " << finalState[i] << " | " << DYee[i] << " | " << DYmm[i] << " | " << DYtt[i] << " | " << tt[i] << " | " << wjets[i] << " | " << wz[i] << " | " << zz[i] << " | " << ww[i] << " | " << tw[i] << " | " << tot[i] << " | " << endl;

  }

    DYee[3] = DYee[1]/DYee[0];
    DYmm[3] = DYmm[1]/DYmm[0];
    DYtt[3] = DYtt[1]/DYtt[0];
    tt[3] = tt[1]/tt[0];
    wjets[3] = wjets[1]/wjets[0];
    wz[3] = wz[1]/wz[0];
    zz[3] = zz[1]/zz[0];
    ww[3] = ww[1]/ww[0];
    tw[3] = tw[1]/tw[0];
    tot[3] = tot[1]/tot[0];

    cout << "| " << finalState[3] << " | " << DYee[3] << " | " << DYmm[3] << " | " << DYtt[3] << " | " << tt[3] << " | " << wjets[3] << " | " << wz[3] << " | " << zz[3] << " | " << ww[3] << " | " << tw[3] << " | " << tot[3] << " | " << endl;


}//end of for loop over njet
}
