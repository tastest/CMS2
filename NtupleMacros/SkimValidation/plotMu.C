{
gSystem->Load("./CMSSW_3_7_0_patch2_V03-05-01/src/CMS2/NtupleMacros/Tools/MiniFWLite/libMiniFWLite.so");
gStyle->SetOptStat("nemou");
TChain *chain = new TChain("tree");
chain->Add("validate.root");

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

TH1F* mu_pt       = new TH1F("mu_pt", "mu_pt", pt_nbins, pt_min, pt_max );
TH1F* mu_eta      = new TH1F("mu_eta", "mu_eta", eta_nbins, eta_min, eta_max );
TH1F* mu_phi      = new TH1F("mu_phi", "mu_phi", phi_nbins, phi_min, phi_max );
TH1F* mu_d0corr   = new TH1F("mu_d0corr", "mu_d0corr", d0_nbins, d0_min, d0_max );
TH1F* mu_Iso      = new TH1F("mu_Iso", "mu_Iso", Iso_nbins, Iso_min, Iso_max );

TH2F* mu_phi_d0   = new TH2F("mu_d0phi", "mu_d0phi", phi_nbins, phi_min, phi_max, d0_nbins, d0_min, d0_max );

chain->Draw("mus.pt() >> mu_pt", "musId");
chain->Draw("mus.eta() >> mu_eta", "musId");
chain->Draw("mus.phi() >> mu_phi", "musId");
chain->Draw("musd0corr >> mu_d0corr", "musId");
chain->Draw("musIso >> mu_Iso", "musId");
chain->Draw("mus.phi():musd0corr >> mu_phi_d0", "musId");

mu_pt->Draw();
new TCanvas();
mu_eta->Draw();
new TCanvas();
mu_phi->Draw();
new TCanvas();
mu_d0corr->Draw();
new TCanvas();
mu_Iso->Draw();
new TCanvas();
mu_phi_d0->Draw();

}
