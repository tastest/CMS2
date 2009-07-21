void doTTconstituents() {
  TFile *_file0 = TFile::Open("processed_data_tag.root");

  float tt[4];
  float ttww[4];
  float ttwo[4];
  float ttwosemi[4];
  float ttwoother[4];
  float other[4];

  for (int i = 0; i < 4; i++){
    tt[i] = 0.0;
    ttww[i] = 0.0;
    ttwosemi[i] = 0.0;
    ttwoother[i] = 0.0;
    other[i] = 0.0;
  }
  
  tt[0]=ttbar_hnJet_ee->Integral();
  ttww[0] = ttbar_hnJetSemiTop_ee->Integral();
  ttwo[0]=ttbar_hnJetSemiWTop_ee->Integral();
  ttwosemi[0]=ttbar_hnJetSemitrueTop_ee->Integral();
  ttwoother[0]= ttbar_hnJetSemiOtherTop_ee->Integral();
  other[0]=ttbar_hnJetfakeLep_ee->Integral();

  tt[1]=ttbar_hnJet_mm->Integral();
  ttww[1] = ttbar_hnJetSemiTop_mm->Integral();
  ttwo[1]=ttbar_hnJetSemiWTop_mm->Integral();
  ttwosemi[1]=ttbar_hnJetSemitrueTop_mm->Integral();
  ttwoother[1]= ttbar_hnJetSemiOtherTop_mm->Integral();
  other[1]=ttbar_hnJetfakeLep_mm->Integral();

  tt[2]=ttbar_hnJet_em->Integral();
  ttww[2] = ttbar_hnJetSemiTop_em->Integral();
  ttwo[2]=ttbar_hnJetSemiWTop_em->Integral();
  ttwosemi[2]=ttbar_hnJetSemitrueTop_em->Integral();
  ttwoother[2]= ttbar_hnJetSemiOtherTop_em->Integral();
  other[2]=ttbar_hnJetfakeLep_em->Integral();

  tt[3]=ttbar_hnJet_all->Integral();
  ttww[3] = ttbar_hnJetSemiTop_all->Integral();
  ttwo[3]=ttbar_hnJetSemiWTop_all->Integral();
  ttwosemi[3]=ttbar_hnJetSemitrueTop_all->Integral();
  ttwoother[3]= ttbar_hnJetSemiOtherTop_all->Integral();
  other[3]=ttbar_hnJetfakeLep_all->Integral();

  cout << "| Same Sign | Type I+II+III | Type-I | Type-II | Type-III | Type-II-SemiLep | Type-II-Fakes | " <<  endl;
  
  for (int i = 0; i < 4 ; i++) {
    cout << "| | " << ttww[i]+ttwo[i]+other[i] << " | " << ttww[i] << " | " << ttwo[i] << " | " << other[i] << " | " << ttwosemi[i] << " | " << ttwoother[i] << " | " << endl;
  }


}
