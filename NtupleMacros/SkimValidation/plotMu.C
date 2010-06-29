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


TH1F* mu1_N        = new TH1F("mu1_N", "mu1_N", N_nbins, N_min, N_max );
TH1F* mu1_pt       = new TH1F("mu1_pt", "mu1_pt", pt_nbins, pt_min, pt_max );
TH1F* mu1_eta      = new TH1F("mu1_eta", "mu1_eta", eta_nbins, eta_min, eta_max );
TH1F* mu1_phi      = new TH1F("mu1_phi", "mu1_phi", phi_nbins, phi_min, phi_max );
TH1F* mu1_d0corr   = new TH1F("mu1_d0corr", "mu1_d0corr", d0_nbins, d0_min, d0_max );
TH1F* mu1_Iso      = new TH1F("mu1_Iso", "mu1_Iso", Iso_nbins, Iso_min, Iso_max );
TH2F* mu1_phi_d0   = new TH2F("mu1_d0phi", "mu1_d0phi", phi_nbins, phi_min, phi_max, d0_nbins, d0_min, d0_max );

TH1F* mu2_N        = new TH1F("mu2_N", "mu2_N", N_nbins, N_min, N_max );
TH1F* mu2_pt       = new TH1F("mu2_pt", "mu2_pt", pt_nbins, pt_min, pt_max );
TH1F* mu2_eta      = new TH1F("mu2_eta", "mu2_eta", eta_nbins, eta_min, eta_max );
TH1F* mu2_phi      = new TH1F("mu2_phi", "mu2_phi", phi_nbins, phi_min, phi_max );
TH1F* mu2_d0corr   = new TH1F("mu2_d0corr", "mu2_d0corr", d0_nbins, d0_min, d0_max );
TH1F* mu2_Iso      = new TH1F("mu2_Iso", "mu2_Iso", Iso_nbins, Iso_min, Iso_max );
TH2F* mu2_phi_d0   = new TH2F("mu2_d0phi", "mu2_d0phi", phi_nbins, phi_min, phi_max, d0_nbins, d0_min, d0_max );

TChain *chain1 = new TChain("tree");
chain1->Add("validate_mus_before.root");

TChain *chain2 = new TChain("tree");
chain2->Add("validate_mus_after.root");

chain1->Draw("nMu >> mu1_N", "musId");
chain1->Draw("musp4.pt() >> mu1_pt", "musId");
chain1->Draw("musp4.eta() >> mu1_eta", "musId");
chain1->Draw("musp4.phi() >> mu1_phi", "musId");
chain1->Draw("musd0corr >> mu1_d0corr", "musId");
chain1->Draw("musIso >> mu1_Iso");
//chain1->Draw("mus.phi():musd0corr >> mu_phi_d0", "musId");

chain2->Draw("nMu >> mu2_N", "musId");
chain2->Draw("musp4.pt() >> mu2_pt", "musId");
chain2->Draw("musp4.eta() >> mu2_eta", "musId");
chain2->Draw("musp4.phi() >> mu2_phi", "musId");
chain2->Draw("musd0corr >> mu2_d0corr", "musId");
chain2->Draw("musIso >> mu2_Iso");
//chain2->Draw("mus.phi():musd0corr >> mu_phi_d0", "musId");

TCanvas *c1 = new TCanvas();
c1->SetWindowSize(1100,850);
c1->Divide(3,2);
c1->cd(1);
mu1_pt->Draw();
c1->cd(2);
mu1_eta->Draw();
c1->cd(3);
mu1_phi->Draw();
c1->cd(4);
mu2_pt->Draw();
c1->cd(5);
mu2_eta->Draw();
c1->cd(6);
mu2_phi->Draw();

TCanvas *c2 = new TCanvas();
c2->SetWindowSize(1100,850);
c2->Divide(3,2);
c2->cd(1);
mu1_N->Draw();
c2->cd(2);
mu1_d0corr->Draw();
c2->cd(3);
mu1_Iso->Draw();
c2->cd(4);
mu2_N->Draw();
c2->cd(5);
mu2_d0corr->Draw();
c2->cd(6);
mu2_Iso->Draw();

}
