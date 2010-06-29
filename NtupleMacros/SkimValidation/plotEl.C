{
gSystem->Load("../validation/CMSSW_3_7_0_patch2_V03-05-01/src/CMS2/NtupleMacros/Tools/MiniFWLite/libMiniFWLite.so");
gStyle->SetOptStat("nemou");

int   N_nbins = 2;
float N_min = 0.5;
float N_max = 2.5;

int pt_nbins = 40;
float pt_min = 10;
float pt_max = 50;

int eta_nbins = 25;
float eta_min = -2.5;
float eta_max = 2.5;

int phi_nbins = 50;
float phi_min = -3.15;
float phi_max = 3.15;

int d0_nbins = 50;
float d0_min = -.1;
float d0_max = .1;

int Iso_nbins = 50;
float Iso_min = 0;
float Iso_max = 1;


TH1F* el1_N        = new TH1F("el1_N", "el1_N", N_nbins, N_min, N_max );
TH1F* el1_pt       = new TH1F("el1_pt", "el1_pt", pt_nbins, pt_min, pt_max );
TH1F* el1_eta      = new TH1F("el1_eta", "el1_eta", eta_nbins, eta_min, eta_max );
TH1F* el1_phi      = new TH1F("el1_phi", "el1_phi", phi_nbins, phi_min, phi_max );
TH1F* el1_d0corr   = new TH1F("el1_d0corr", "el1_d0corr", d0_nbins, d0_min, d0_max );
TH1F* el1_Iso      = new TH1F("el1_Iso", "el1_Iso", Iso_nbins, Iso_min, Iso_max );
TH2F* el1_phi_d0   = new TH2F("el1_d0phi", "el1_d0phi", phi_nbins, phi_min, phi_max, d0_nbins, d0_min, d0_max );

TH1F* el2_N        = new TH1F("el2_N", "el2_N", N_nbins, N_min, N_max );
TH1F* el2_pt       = new TH1F("el2_pt", "el2_pt", pt_nbins, pt_min, pt_max );
TH1F* el2_eta      = new TH1F("el2_eta", "el2_eta", eta_nbins, eta_min, eta_max );
TH1F* el2_phi      = new TH1F("el2_phi", "el2_phi", phi_nbins, phi_min, phi_max );
TH1F* el2_d0corr   = new TH1F("el2_d0corr", "el2_d0corr", d0_nbins, d0_min, d0_max );
TH1F* el2_Iso      = new TH1F("el2_Iso", "el2_Iso", Iso_nbins, Iso_min, Iso_max );
TH2F* el2_phi_d0   = new TH2F("el2_d0phi", "el2_d0phi", phi_nbins, phi_min, phi_max, d0_nbins, d0_min, d0_max );

TChain *chain1 = new TChain("tree");
chain1->Add("validate_els_before.root");

TChain *chain2 = new TChain("tree");
chain2->Add("validate_els_after.root");


chain1->Draw("nels >> el1_N", "elsId");
chain1->Draw("elp4.pt() >> el1_pt", "elsId");
chain1->Draw("elp4.eta() >> el1_eta", "elsId");
chain1->Draw("elp4.phi() >> el1_phi", "elsId");
//chain1->Draw("eld0corr >> el1_d0corr", "elsId");
chain1->Draw("elreliso >> el1_Iso");
//add the three individual isos
//chain1->Draw("els.phi():elsd0corr >> el_phi_d0", "elsId");

chain2->Draw("nEl >> el2_N", "elsId");
chain2->Draw("elp4.pt() >> el2_pt", "elsId");
chain2->Draw("elp4.eta() >> el2_eta", "elsId");
chain2->Draw("elp4.phi() >> el2_phi", "elsId");
//chain2->Draw("eld0corr >> el2_d0corr", "elsId");
chain2->Draw("elreliso >> el2_Iso");
//chain2->Draw("els.phi():elsd0corr >> el_phi_d0", "elsId");

TCanvas *c1 = new TCanvas();
c1->SetWindowSize(1100,850);
c1->Divide(3,2);
c1->cd(1);
el1_pt->Draw();
c1->cd(2);
el1_eta->Draw();
c1->cd(3);
el1_phi->Draw();
c1->cd(4);
el2_pt->Draw();
c1->cd(5);
el2_eta->Draw();
c1->cd(6);
el2_phi->Draw();

TCanvas *c2 = new TCanvas();
c2->SetWindowSize(1100,850);
c2->Divide(3,2);
c2->cd(1);
el1_N->Draw();
c2->cd(2);
el1_d0corr->Draw();
c2->cd(3);
el1_Iso->Draw();
c2->cd(4);
el2_N->Draw();
c2->cd(5);
el2_d0corr->Draw();
c2->cd(6);
el2_Iso->Draw();

}
