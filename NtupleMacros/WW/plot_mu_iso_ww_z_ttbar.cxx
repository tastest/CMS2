{
z = new TChain("Events");
z->Add("/data/tmp/cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*root");

ttbar = new TChain("Events");
ttbar->Add("/data/tmp/cms2-V01-02-06/TTJets-madgraph_Fall08_IDEAL_V9_v2/merged_ntuple*root");

ww = new TChain("Events");
ww->Add("/data/tmp/cms2-V01-02-06/WW_2l_Summer08_IDEAL_V9_v2/merged_ntuple*root");

TCut base_cut("mus_p4.pt()>20&&(mus_type&0x2)&&abs(mus_d0corr)<0.2&&mus_gfit_chi2/mus_gfit_ndof<10&&mus_validHits>=11");
TCut zee(     "abs(mus_mc_id)==13&&mus_mc_motherid==23");
TCut ttbar_e( "abs(mus_mc_id)==13&&abs(mus_mc_motherid)==24");
TCut ww_e(    "abs(mus_mc_id)==13&&abs(mus_mc_motherid)==24");
c1 = new TCanvas("plot_mu_iso_ww_z_ttbar","plot_mu_iso_ww_z_ttbar",1200,400);

z->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_z(100,0,1)",
	base_cut && zee,"goff");
ttbar->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_ttbar(100,0,1)",
	    base_cut && ttbar_e,"goff");
ww->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_ww(100,0,1)",
	 base_cut && ww_e,"goff");

c1->Divide(3,1);
c1->cd(1);
h_z->SetTitle("Isolated muons from Z events (truth matched)");
h_z->GetXaxis()->SetTitle("relative isolation");
h_z->SetFillColor(kMagenta);
h_z->Draw();
c1->cd(2);
h_ttbar->SetTitle("Isolated muons from TTbar events (truth matched)");
h_ttbar->GetXaxis()->SetTitle("relative isolation");
h_ttbar->SetFillColor(kYellow);
h_ttbar->Draw();
c1->cd(3);
h_ww->SetTitle("Isolated muons from WW events (truth matched)");
h_ww->GetXaxis()->SetTitle("relative isolation");
h_ww->SetFillColor(kBlue);
h_ww->Draw();

}
