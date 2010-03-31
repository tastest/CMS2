{
    gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");
    TChain *chain_v0 = new TChain("Events");
    chain_v0->Add("/tas03/disk02/slava77/reltestdata/CMSSW_3_5_6-cms2-data/*.root.v0");                          
    TString cut = "evt_nels == 1 && els_egamma_looseId[0] == 1 && els_p4[0].Pt() > 5.0";
//    chain_v0->Scan("evt_event:evt_lumiBlock:evt_run:els_p4[0].Eta():els_p4[0].Pt():els_type[0]:(els_eMax[0]/els_e5x5[0]):els_fbrem:els_eOverPIn:els_egamma_looseId[0]", cut);  
//    chain_v0->Scan("evt_event:evt_lumiBlock:evt_run:els_p4[0].Pt():els_hOverE:els_dEtaIn:els_dPhiIn:els_eSeedOverPIn:els_sigmaIEtaIEta:els_e2x5Max/els_e5x5:els_d0corr:els_tkIso:els_hcalIso:els_ecalIso:els_egamma_looseId[0]", cut);  
    chain_v0->Scan("evt_event:evt_lumiBlock:evt_run:els_p4[0].Pt():els_conv_dist:els_conv_dcot:els_exp_innerlayers", cut);



}
