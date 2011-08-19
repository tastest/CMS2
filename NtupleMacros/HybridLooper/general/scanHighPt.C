{
    gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");
    TChain *chain_v0 = new TChain("Events");
    chain_v0->Add("/tas03/disk02/slava77/reltestdata/CMSSW_3_5_6-cms2-data/*.root.v0");                          
    TString cut = "evt_nels == 1 && els_p4[0].Pt() > 30";
    chain_v0->Scan("evt_event:evt_lumiBlock:evt_run:evt_pfmet:evt_tcmet:els_p4[0].Pt():acos(cos(evt_tcmetPhi - els_p4[0].Phi())):els_type[0]:(els_eMax[0]/els_e5x5[0])", cut);  

}
