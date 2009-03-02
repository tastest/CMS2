{
chain = new TChain("Events");
chain->Add("/data/tmp/cms2-V01-02-06/WJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*root");
chain->Draw(">>mufakes","abs(genps_id)==11");
chain->SetEventList(mufakes);
cout << "Looking only at events with electorn in gen block (W->enu)" << endl;

TCut base_cut("mus_p4.pt()>20&&(mus_type&0x2)&&abs(mus_d0corr)<0.2&&mus_gfit_chi2/mus_gfit_ndof<10&&mus_validHits>=11");
TCut bad_global_fit("abs(mus_p4.pt()/mus_trk_p4.pt()-1)>0.1||mus_type!=14");
TCut mcmuon("abs(mus_mc_id)==13");

h_mcmuons = new TH1F("h_mcmuons","h_mcmuons",20,0,1.);
chain->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_mcmuons",base_cut&&mcmuon,"goff");

h_badgfit = new TH1F("h_badgfit","h_badgfit",20,0,1.);
chain->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_badgfit",base_cut&&(!mcmuon)&&bad_global_fit,"goff");

h_fakes = new TH1F("h_fakes","h_fakes",20,0,1.);
chain->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_fakes",base_cut&&(!mcmuon)&&(!bad_global_fit),"goff");

stack = new THStack("mufakes","Muon Fakes");
h_mcmuons->SetFillColor(kBlue);
h_badgfit->SetFillColor(kMagenta);
h_fakes->SetFillColor(kYellow);
stack->Add(h_fakes);
stack->Add(h_badgfit);
stack->Add(h_mcmuons);
stack->Draw();

}
