{
    gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");
    TChain *chain_v0 = new TChain("Events");
    //chain_v0->Add("/tas03/disk02/slava77/reltestdata/CMSSW_3_5_6-cms2-data/*.root.v0");                          
    // slavas 15 gev electron
    chain_v0->Add("/tas03/disk02/slava77/reltestdata/CMSSW_3_5_6-cms2-data/_cms__store_express_Commissioning10_ExpressPhysics_FEVT_v7_000_132_442_FACFAE07-033C-DF11-8746-0019B9F70607.root.v0");
    TString cut = "evt_run == 132442 && evt_lumiBlock == 72 && evt_event == 1739397 && evt_nels == 2 && (els_p4[0].Pt() > 5.0 || els_p4[0].Pt() > 5.0)";
//    chain_v0->Scan("evt_event:evt_lumiBlock:evt_run:els_p4[0].Eta():els_p4[0].Pt():els_type[0]:(els_eMax[0]/els_e5x5[0]):els_fbrem:els_eOverPIn:els_egamma_looseId[0]", cut);  
    chain_v0->Scan("evt_event:evt_lumiBlock:evt_run:els_p4[0].Pt():els_hOverE[0]:els_dEtaIn[0]:els_dPhiIn[0]:els_eSeedOverPIn[0]:els_sigmaIEtaIEta[0]:els_e2x5Max[0]/els_e5x5[0]:els_d0corr[0]:els_tkIso[0]:els_hcalIso[0]:els_ecalIso[0]:els_egamma_looseId[0]", cut);  
    //chain_v0->Scan("evt_event:evt_lumiBlock:evt_run:els_p4[0].Pt():els_p4[1].Pt()", cut);



}
