{
qcd30 = new TChain("Events");
qcd30->Add("/data/tmp/cms2-V01-02-06/QCDpt30/merged_ntuple*root");

qcd80 = new TChain("Events");
qcd80->Add("/data/tmp/cms2-V01-02-06/QCDpt80/merged_ntuple*root");

qcd170 = new TChain("Events");
qcd170->Add("/data/tmp/cms2-V01-02-06/QCDpt170/merged_ntuple*root");

qcd300 = new TChain("Events");
qcd300->Add("/data/tmp/cms2-V01-02-06/QCDpt300/merged_ntuple*root");

TCut base_cut("mus_p4.pt()>20&&(mus_type&0x2)&&abs(mus_d0corr)<0.2&&mus_gfit_chi2/mus_gfit_ndof<10&&mus_validHits>=11");
TCut bad_global_fit("abs(mus_p4.pt()/mus_trk_p4.pt()-1)>0.1||mus_type!=14");
TCut mcmuon("abs(mus_mc_id)==13");

qcd30->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_qcd30_mc(20,0.5,1)",
	    base_cut+mcmuon,"goff");
qcd30->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_qcd30_bad(20,0.5,1)",
	    base_cut+bad_global_fit,"goff");
qcd30->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_qcd30_fake(20,0.5,1)",
	    base_cut+!bad_global_fit+!mcmuon,"goff");

qcd80->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_qcd80_mc(20,0.5,1)",
	    base_cut+mcmuon,"goff");
qcd80->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_qcd80_bad(20,0.5,1)",
	    base_cut+bad_global_fit,"goff");
qcd80->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_qcd80_fake(20,0.5,1)",
	    base_cut+!bad_global_fit+!mcmuon,"goff");

qcd170->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_qcd170_mc(20,0.5,1)",
	    base_cut+mcmuon,"goff");
qcd170->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_qcd170_bad(20,0.5,1)",
	    base_cut+bad_global_fit,"goff");
qcd170->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_qcd170_fake(20,0.5,1)",
	    base_cut+!bad_global_fit+!mcmuon,"goff");

qcd300->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_qcd300_mc(20,0.5,1)",
	    base_cut+mcmuon,"goff");
qcd300->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_qcd300_bad(20,0.5,1)",
	    base_cut+bad_global_fit,"goff");
qcd300->Draw("mus_p4.pt()/(mus_p4.pt()+mus_pat_ecalIso+mus_pat_hcalIso+mus_pat_trackIso+1e-5)>>h_qcd300_fake(20,0.5,1)",
	    base_cut+!bad_global_fit+!mcmuon,"goff");

c1 = new TCanvas("plot_mu_iso_qcd","plot_mu_iso_qcd",800,800);
c1->Divide(2,2);
c1->cd(1);
stack30 = new THStack("s_qcd30","Muons from QCD events (#hat{pt}>30GeV)");
h_qcd30_fake->SetFillColor(kYellow);
h_qcd30_mc->SetFillColor(kBlue);
h_qcd30_bad->SetFillColor(kMagenta);
stack30->Add(h_qcd30_fake);
stack30->Add(h_qcd30_bad);
stack30->Add(h_qcd30_mc);
stack30->Draw();

c1->cd(2);
stack80 = new THStack("s_qcd80","Muons from QCD events (#hat{pt}>80GeV)");
h_qcd80_fake->SetFillColor(kYellow);
h_qcd80_mc->SetFillColor(kBlue);
h_qcd80_bad->SetFillColor(kMagenta);
stack80->Add(h_qcd80_fake);
stack80->Add(h_qcd80_bad);
stack80->Add(h_qcd80_mc);
stack80->Draw();

c1->cd(3);
stack170 = new THStack("s_qcd170","Muons from QCD events (#hat{pt}>170GeV)");
h_qcd170_fake->SetFillColor(kYellow);
h_qcd170_mc->SetFillColor(kBlue);
h_qcd170_bad->SetFillColor(kMagenta);
stack170->Add(h_qcd170_fake);
stack170->Add(h_qcd170_bad);
stack170->Add(h_qcd170_mc);
stack170->Draw();

c1->cd(4);
stack300 = new THStack("s_qcd300","Muons from QCD events (#hat{pt}>300GeV)");
h_qcd300_fake->SetFillColor(kYellow);
h_qcd300_mc->SetFillColor(kBlue);
h_qcd300_bad->SetFillColor(kMagenta);
stack300->Add(h_qcd300_fake);
stack300->Add(h_qcd300_bad);
stack300->Add(h_qcd300_mc);
stack300->Draw();

}
