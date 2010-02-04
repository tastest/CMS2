void mcDilEffs(char* fNam){
  gSystem->Load("libMiniFWLite.so");

  TChain* e = new TChain("Events");
  e->Add(fNam);

  e->Draw("evt_nEvts:evt_scale1fb","","",1);
  double nEvts = e->GetV1()[0];
  double scale1pb = e->GetV2()[];
  scale1pb*=1e-3;

  std::cout<<"Processed "<<nEvts<<" expect "<<nEvts*scale1pb<<" pb"<<std::endl;
  double nAll = nEvts*scale1pb;

  //Now, obviously this sample should have dilepton gen filter
  double nEEAll = e->GetEntries("Sum$(abs(genps_id)==11)==2"); nEEAll*=scale1pb;
  double nMMAll = e->GetEntries("Sum$(abs(genps_id)==13)==2"); nMMAll*=scale1pb;
  double nEMAll = e->GetEntries("Sum$(abs(genps_id)==11||abs(genps_id)==13)==2");
  nEMAll*=scale1pb; nEMAll -= (nEEAll+nMMAll);
  std::cout<<"EE:MM:EM fEE:fMM:fEM "
	   << nEEAll <<" "<< nMMAll <<" "<< nEMAll <<" "
	   <<nEEAll/nAll<<" "<<nMMAll/nAll<<" "<<nEMAll/nAll
	   <<std::endl;

  // now count fiducial dilepton events
  double nEE2020 = e->GetEntries("Sum$(abs(genps_id)==11&&abs(genps_p4.eta())<2.4&&genps_p4.pt()>20.)==2"); nEE2020*=scale1pb;
  double nMM2020 = e->GetEntries("Sum$(abs(genps_id)==13&&abs(genps_p4.eta())<2.4&&genps_p4.pt()>20.)==2"); nMM2020*=scale1pb;
  double nEM2020 = e->GetEntries("Sum$((abs(genps_id)==11||abs(genps_id)==13)&&abs(genps_p4.eta())<2.4&&genps_p4.pt()>20.)==2");
  nEM2020*=scale1pb; nEM2020 -= (nEE2020+nMM2020);
  std::cout<<"   * 2020 eta<2.4 EE:MM:EM fEE:fMM:fEM frEE:frMM:frEM "<<std::endl;
  std::cout<<"      * "<< nEE2020 <<" "<< nMM2020 <<" "<< nEM2020 <<std::endl;
  std::cout<<"      * "<< nEE2020/nAll<<" "<<nMM2020/nAll<<" "<<nEM2020/nAll<<std::endl;
  std::cout<<"      * "<< nEE2020/nEEAll<<" "<<nMM2020/nMMAll<<" "<<nEM2020/nEMAll<<std::endl;

  //look at events with two reco leptons with pt>20, eta<2.4
  double nEE2020W2Reco = e->GetEntries("Sum$(abs(genps_id)==11&&abs(genps_p4.eta())<2.4&&genps_p4.pt()>20.)==2&&Sum$(els_p4.pt()>20&&abs(els_p4.eta())<2.4)>1");
  nEE2020W2Reco*=scale1pb;
  double nMM2020W2Reco = e->GetEntries("Sum$(abs(genps_id)==13&&abs(genps_p4.eta())<2.4&&genps_p4.pt()>20.)==2&&Sum$(mus_p4.pt()>20&&abs(mus_p4.eta())<2.4)>1");
  nMM2020W2Reco*=scale1pb;
  double nEM2020W2Reco = e->GetEntries("Sum$((abs(genps_id)==11||abs(genps_id)==13)&&abs(genps_p4.eta())<2.4&&genps_p4.pt()>20.)==2&&(Sum$(mus_p4.pt()>20&&abs(mus_p4.eta())<2.4)+Sum$(els_p4.pt()>20&&abs(els_p4.eta())<2.4&&els_closestMuon==-1))>1");
  nEM2020W2Reco*=scale1pb;  nEM2020W2Reco-= (nEE2020W2Reco+nMM2020W2Reco);
  std::cout<<"   * 2020 eta<2.4 w2Reco EE:MM:EM fEE:fMM:fEM (no bratio)  frEE:frMM:frEM (wrt previous sels)"<<std::endl;
  std::cout<<"      * "<< nEE2020W2Reco <<" "<< nMM2020W2Reco <<" "<< nEM2020W2Reco <<std::endl;
  std::cout<<"      * "<<nEE2020W2Reco/nEEAll<<" "<<nMM2020W2Reco/nMMAll<<" "<<nEM2020W2Reco/nEMAll<<std::endl;
  std::cout<<"      * "<<nEE2020W2Reco/nEE2020<<" "<<nMM2020W2Reco/nMM2020<<" "<<nEM2020W2Reco/nEM2020<<std::endl;

  //look at events with two reco leptons with pt>20, eta<2.4 and a match to MC with parent<50
  double nEE2020W2RecoMatch = e->GetEntries("Sum$(abs(genps_id)==11&&abs(genps_p4.eta())<2.4&&genps_p4.pt()>20.)==2&&Sum$(els_p4.pt()>20&&abs(els_p4.eta())<2.4&&abs(els_mc_motherid)<50&&(abs(els_mc_id)==11||(abs(els_mc_id)==22&&abs(els_mc_motherid)==11)))>1");
  nEE2020W2RecoMatch*=scale1pb;
  double nMM2020W2RecoMatch = e->GetEntries("Sum$(abs(genps_id)==13&&abs(genps_p4.eta())<2.4&&genps_p4.pt()>20.)==2&&Sum$(mus_p4.pt()>20&&abs(mus_p4.eta())<2.4&&abs(mus_mc_motherid)<50&&(abs(mus_mc_id)==13||(abs(mus_mc_id)==22&&abs(mus_mc_motherid)==13)))>1");
  nMM2020W2RecoMatch*=scale1pb;
  double nEM2020W2RecoMatch = e->GetEntries("Sum$((abs(genps_id)==11||abs(genps_id)==13)&&abs(genps_p4.eta())<2.4&&genps_p4.pt()>20.)==2&&(Sum$(mus_p4.pt()>20&&abs(mus_p4.eta())<2.4&&abs(mus_mc_motherid)<50&&(abs(mus_mc_id)==13||(abs(mus_mc_id)==22&&abs(mus_mc_motherid)==13)))+Sum$(els_p4.pt()>20&&abs(els_p4.eta())<2.4&&els_closestMuon==-1&&abs(els_mc_motherid)<50&&(abs(els_mc_id)==11||(abs(els_mc_id)==22&&abs(els_mc_motherid)==11))))>1");
  nEM2020W2RecoMatch*=scale1pb;  nEM2020W2RecoMatch-= (nEE2020W2RecoMatch+nMM2020W2RecoMatch);
  std::cout<<"   * 2020 eta<2.4 w2Reco matching MC w mother<50 EE:MM:EM fEE:fMM:fEM (no bratio)  frEE:frMM:frEM (wrt previous sels)"<<std::endl;
  std::cout<<"      * "<< nEE2020W2RecoMatch <<" "<< nMM2020W2RecoMatch <<" "<< nEM2020W2RecoMatch <<std::endl;
  std::cout<<"      * "<<nEE2020W2RecoMatch/nEEAll<<" "<<nMM2020W2RecoMatch/nMMAll<<" "<<nEM2020W2RecoMatch/nEMAll<<std::endl;
  std::cout<<"      * "<<nEE2020W2RecoMatch/nEE2020W2Reco<<" "<<nMM2020W2RecoMatch/nMM2020W2Reco<<" "<<nEM2020W2RecoMatch/nEM2020W2Reco<<std::endl;


  //look at events with two reco leptons with pt>20, eta<2.4 and a match to MC with parent<50
  double nEEW2Reco2020Match = e->GetEntries("Sum$(abs(genps_id)==11)==2&&Sum$(els_p4.pt()>20&&abs(els_p4.eta())<2.4&&abs(els_mc_motherid)<50&&(abs(els_mc_id)==11||(abs(els_mc_id)==22&&abs(els_mc_motherid)==11)))>1");
  nEEW2Reco2020Match*=scale1pb;
  double nMMW2Reco2020Match = e->GetEntries("Sum$(abs(genps_id)==13)==2&&Sum$(mus_p4.pt()>20&&abs(mus_p4.eta())<2.4&&abs(mus_mc_motherid)<50&&(abs(mus_mc_id)==13||(abs(mus_mc_id)==22&&abs(mus_mc_motherid)==13)))>1");
  nMMW2Reco2020Match*=scale1pb;
  double nEMW2Reco2020Match = e->GetEntries("Sum$((abs(genps_id)==11||abs(genps_id)==13))==2&&(Sum$(mus_p4.pt()>20&&abs(mus_p4.eta())<2.4&&abs(mus_mc_motherid)<50&&(abs(mus_mc_id)==13||(abs(mus_mc_id)==22&&abs(mus_mc_motherid)==13)))+Sum$(els_p4.pt()>20&&abs(els_p4.eta())<2.4&&els_closestMuon==-1&&abs(els_mc_motherid)<50&&(abs(els_mc_id)==11||(abs(els_mc_id)==22&&abs(els_mc_motherid)==11))))>1");
  nEMW2Reco2020Match*=scale1pb;  nEMW2Reco2020Match-= (nEEW2Reco2020Match+nMMW2Reco2020Match);
  std::cout<<"   * Reco 2020 eta<2.4 (ee/mumu genp) w matching MC w mother<50 EE:MM:EM fEE:fMM:fEM (this/genp2020)"<<std::endl;
  std::cout<<"      * "<< nEEW2Reco2020Match <<" "<< nMMW2Reco2020Match <<" "<< nEMW2Reco2020Match <<std::endl;
  std::cout<<"      * "<< nEEW2Reco2020Match/nEE2020 <<" "<< nMMW2Reco2020Match/nMM2020 <<" "<< nEMW2Reco2020Match/nEM2020 <<std::endl;

  //look at events with two reco leptons with pt>20, eta<2.4 and a match to MC with parent<50
  double nEE1515LW2Reco2020Match = e->GetEntries("Sum$(abs(genps_id)==11&&abs(genps_p4.eta())<2.7&&genps_p4.pt()>15.)==2&&Sum$(els_p4.pt()>20&&abs(els_p4.eta())<2.4&&abs(els_mc_motherid)<50&&(abs(els_mc_id)==11||(abs(els_mc_id)==22&&abs(els_mc_motherid)==11)))>1");
  nEE1515LW2Reco2020Match*=scale1pb;
  double nMM1515LW2Reco2020Match = e->GetEntries("Sum$(abs(genps_id)==13&&abs(genps_p4.eta())<2.7&&genps_p4.pt()>15.)==2&&Sum$(mus_p4.pt()>20&&abs(mus_p4.eta())<2.4&&abs(mus_mc_motherid)<50&&(abs(mus_mc_id)==13||(abs(mus_mc_id)==22&&abs(mus_mc_motherid)==13)))>1");
  nMM1515LW2Reco2020Match*=scale1pb;
  double nEM1515LW2Reco2020Match = e->GetEntries("Sum$((abs(genps_id)==11||abs(genps_id)==13)&&abs(genps_p4.eta())<2.7&&genps_p4.pt()>15.)==2&&(Sum$(mus_p4.pt()>20&&abs(mus_p4.eta())<2.4&&abs(mus_mc_motherid)<50&&(abs(mus_mc_id)==13||(abs(mus_mc_id)==22&&abs(mus_mc_motherid)==13)))+Sum$(els_p4.pt()>20&&abs(els_p4.eta())<2.4&&els_closestMuon==-1&&abs(els_mc_motherid)<50&&(abs(els_mc_id)==11||(abs(els_mc_id)==22&&abs(els_mc_motherid)==11))))>1");
  nEM1515LW2Reco2020Match*=scale1pb;  nEM1515LW2Reco2020Match-= (nEE1515LW2Reco2020Match+nMM1515LW2Reco2020Match);
  std::cout<<"   * 2020 eta<2.4 w2Reco matching MC w mother<50 EE:MM:EM fEE:fMM:fEM (this/genp2020)"<<std::endl;
  std::cout<<"      * "<< nEE1515LW2Reco2020Match <<" "<< nMM1515LW2Reco2020Match <<" "<< nEM1515LW2Reco2020Match <<std::endl;
  std::cout<<"      * "<< nEE1515LW2Reco2020Match/nEE2020 <<" "<< nMM1515LW2Reco2020Match/nMM2020 <<" "<< nEM1515LW2Reco2020Match/nEM2020 <<std::endl;


}
