void skimMusMET(TChain* eq){
eq->SetBranchStatus("*", 0);
eq->SetBranchStatus("*_muonMaker_musp4_*", 1);
eq->SetBranchStatus("*_muonMaker_muse*_*", 1);
eq->SetBranchStatus("*_muonMaker_mustype*_*", 1);
eq->SetBranchStatus("*_muonMaker_musvalidHits*_*", 1);
eq->SetBranchStatus("*_muonMaker_musgfit*_*", 1);
eq->SetBranchStatus("*_muonMaker_musd0corr*_*", 1);
eq->SetBranchStatus("*_muonMaker_musz0corr*_*", 1);
eq->SetBranchStatus("*_candToGenAssMaker_musmc*_*", 1);
eq->SetBranchStatus("*_patMuonMaker_*_*", 1);
eq->SetBranchStatus("*_eventMaker_*_*", 1);
eq->SetBranchStatus("*_eventMaker_*HLT*_*", 0);
eq->SetBranchStatus("*_eventMaker_*L1*_*", 0);
eq->SetBranchStatus("*_patMETMaker_*metpatmet*Cor_*", 1);
}
