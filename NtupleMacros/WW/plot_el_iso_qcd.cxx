{
qcd30 = new TChain("Events");
qcd30->Add("/data/tmp/cms2-V01-02-06/QCDpt30/merged_ntuple*root");

qcd80 = new TChain("Events");
qcd80->Add("/data/tmp/cms2-V01-02-06/QCDpt80/merged_ntuple*root");

qcd170 = new TChain("Events");
qcd170->Add("/data/tmp/cms2-V01-02-06/QCDpt170/merged_ntuple*root");

qcd300 = new TChain("Events");
qcd300->Add("/data/tmp/cms2-V01-02-06/QCDpt300/merged_ntuple*root");

TCut base_cut("els_p4.pt()>20&&els_tightId22XMinMatteo&&abs(els_d0corr)<0.025");
TCut conversion("abs(els_mc_id)==22");
TCut heavy(     "abs(els_mc_id)==11&&els_mc_motherid!=111");
TCut piz3body(  "abs(els_mc_id)==11&&els_mc_motherid==111");

qcd30->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd30_piz(20,0.5,1)",
	    base_cut+piz3body,"goff");
qcd30->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd30_conv(20,0.5,1)",
	    base_cut+conversion,"goff");
qcd30->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd30_heavy(20,0.5,1)",
	    base_cut+heavy,"goff");
qcd30->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd30_fake(20,0.5,1)",
	    base_cut+!heavy+!piz3body+!conversion,"goff");

qcd80->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd80_piz(20,0.5,1)",
	    base_cut+piz3body,"goff");
qcd80->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd80_conv(20,0.5,1)",
	    base_cut+conversion,"goff");
qcd80->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd80_heavy(20,0.5,1)",
	    base_cut+heavy,"goff");
qcd80->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd80_fake(20,0.5,1)",
	    base_cut+!heavy+!piz3body+!conversion,"goff");

qcd170->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd170_piz(20,0.5,1)",
	    base_cut+piz3body,"goff");
qcd170->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd170_conv(20,0.5,1)",
	    base_cut+conversion,"goff");
qcd170->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd170_heavy(20,0.5,1)",
	    base_cut+heavy,"goff");
qcd170->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd170_fake(20,0.5,1)",
	    base_cut+!heavy+!piz3body+!conversion,"goff");

qcd300->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd300_piz(20,0.5,1)",
	    base_cut+piz3body,"goff");
qcd300->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd300_conv(20,0.5,1)",
	    base_cut+conversion,"goff");
qcd300->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd300_heavy(20,0.5,1)",
	    base_cut+heavy,"goff");
qcd300->Draw("els_p4.pt()/(els_p4.pt()+els_pat_ecalIso+els_pat_hcalIso+els_pat_trackIso+1e-5)>>h_qcd300_fake(20,0.5,1)",
	    base_cut+!heavy+!piz3body+!conversion,"goff");


c1 = new TCanvas("plot_el_iso_qcd","plot_el_iso_qcd",800,800);
c1->Divide(2,2);
c1->cd(1);
stack30 = new THStack("s_qcd30","Electrons from QCD events (#hat{pt}>30GeV)");
h_qcd30_fake->SetFillColor(kYellow);
h_qcd30_heavy->SetFillColor(kRed);
h_qcd30_conv->SetFillColor(kGray);
h_qcd30_piz->SetFillColor(kBlack);
stack30->Add(h_qcd30_piz);
stack30->Add(h_qcd30_heavy);
stack30->Add(h_qcd30_conv);
stack30->Add(h_qcd30_fake);
stack30->Draw();

c1->cd(2);
stack80 = new THStack("s_qcd80","Electrons from QCD events (#hat{pt}>80GeV)");
h_qcd80_fake->SetFillColor(kYellow);
h_qcd80_heavy->SetFillColor(kRed);
h_qcd80_conv->SetFillColor(kGray);
h_qcd80_piz->SetFillColor(kBlack);
stack80->Add(h_qcd80_piz);
stack80->Add(h_qcd80_heavy);
stack80->Add(h_qcd80_conv);
stack80->Add(h_qcd80_fake);
stack80->Draw();

c1->cd(3);
stack170 = new THStack("s_qcd170","Electrons from QCD events (#hat{pt}>170GeV)");
h_qcd170_fake->SetFillColor(kYellow);
h_qcd170_heavy->SetFillColor(kRed);
h_qcd170_conv->SetFillColor(kGray);
h_qcd170_piz->SetFillColor(kBlack);
stack170->Add(h_qcd170_piz);
stack170->Add(h_qcd170_heavy);
stack170->Add(h_qcd170_conv);
stack170->Add(h_qcd170_fake);
stack170->Draw();

c1->cd(4);
stack300 = new THStack("s_qcd300","Electrons from QCD events (#hat{pt}>300GeV)");
h_qcd300_fake->SetFillColor(kYellow);
h_qcd300_heavy->SetFillColor(kRed);
h_qcd300_conv->SetFillColor(kGray);
h_qcd300_piz->SetFillColor(kBlack);
stack300->Add(h_qcd300_piz);
stack300->Add(h_qcd300_heavy);
stack300->Add(h_qcd300_conv);
stack300->Add(h_qcd300_fake);
stack300->Draw();

}
