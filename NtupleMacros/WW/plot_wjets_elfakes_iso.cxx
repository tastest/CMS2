{
chain = new TChain("Events");
chain->Add("/data/tmp/cms2-V01-02-06/WJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*root");
chain->Draw(">>elfakes","abs(genps_id)==13");
chain->SetEventList(elfakes);
cout << "Looking only at events with muon in gen block (W->munu)" << endl;

TCut base_cut("els_p4.pt()>20&&els_tightId22XMinMatteo&&abs(els_d0corr)<0.025");
TCut wmuon(     "abs(els_mc_id)==13&&abs(els_mc_motherid)==24");
TCut mugamma(   "abs(els_mc_id)==22&&abs(els_mc_motherid)==13");
TCut conversion("abs(els_mc_id)==22&&abs(els_mc_motherid)!=13");
TCut heavy(     "abs(els_mc_id)==11&&els_mc_motherid!=111");
TCut piz3body(  "abs(els_mc_id)==11&&els_mc_motherid==111");



chain->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_wmuons(20,0,1)",
	    (base_cut && wmuon)*"evt_scale1fb*0.1","goff");
chain->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_mugamma(20,0,1)",
	    (base_cut && mugamma)*"evt_scale1fb*0.1","goff");
chain->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_conver(20,0,1)",
	    (base_cut && conversion)*"evt_scale1fb*0.1","goff");
chain->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_heavy(20,0,1)",
	    (base_cut && heavy)*"evt_scale1fb*0.1","goff");
chain->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_piz(20,0,1)",
	    (base_cut && piz3body)*"evt_scale1fb*0.1","goff");
chain->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_fakes(20,0,1)",
	    (base_cut && !wmuon && !mugamma && !heavy && !piz3body && !conversion)*"evt_scale1fb*0.1","goff");

stack = new THStack("elfakes","Electron Fakes (no muon veto)");
h_wmuons->SetFillColor(kBlue);
h_mugamma->SetFillColor(kMagenta);
h_fakes->SetFillColor(kYellow);
h_heavy->SetFillColor(kRed);
h_conver->SetFillColor(kGray);
h_piz->SetFillColor(kBlack);

stack->Add(h_piz);
stack->Add(h_heavy);
stack->Add(h_conver);
stack->Add(h_fakes);
stack->Add(h_mugamma);
stack->Add(h_wmuons);
stack->Draw();

}
