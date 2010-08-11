{
gSystem->Load("../validation/CMSSW_3_7_0_patch2_V03-05-01/src/CMS2/NtupleMacros/Tools/MiniFWLite/libMiniFWLite.so");
gROOT->SetStyle("Plain");
gStyle->SetHistFillColor(kRed);
gStyle->SetHistMinielmZero();
gStyle->SetOptStat("nemoui");

const string suffix = ".png";
const string dir = "plots/";

int   N_nbins = 2;
float N_min = 0.5;
float N_max = 2.5;

int pt_nbins = 40;
float pt_min = 10;
float pt_max = 50;

int eta_nbins = 25;
float eta_min = -2.5;
float eta_max = 2.5;

int phi_nbins = 30;
float phi_min = -3.15;
float phi_max = 3.15;

int d0_nbins = 50;
float d0_min = -.04;
float d0_max = .04;

int Iso_nbins = 100;
float Iso_min = 0;
float Iso_max = 1;

int   met_nbins = 75;
float met_min   = 0;
float met_max   = 75;

int jet_pt_nbins  = 70;
float jet_pt_min  = 30;
float jet_pt_max  = 100;

int jet_eta_nbins = 25;
float jet_eta_min = -2.4;
float jet_eta_max = 2.4;

int Njet_nbins  = 6;
int Njet_min    = 0;
int Njet_max    = 6; 

// elon quantities - before
TH1F* el1_N        = new TH1F("el1_Nels", "el1_Nels", N_nbins, N_min, N_max );
TH1F* el1_Niso     = new TH1F("el1_Nelsiso", "el1_Nels_isolated", N_nbins, N_min, N_max );
TH1F* el1_pt       = new TH1F("el1_pt", "el1_pt", pt_nbins, pt_min, pt_max );
TH1F* el1_eta      = new TH1F("el1_eta", "el1_eta", eta_nbins, eta_min, eta_max );
TH1F* el1_phi      = new TH1F("el1_phi", "el1_phi", phi_nbins, phi_min, phi_max );
TH1F* el1_d0corr   = new TH1F("el1_d0corr", "el1_d0corr", d0_nbins, d0_min, d0_max );
TH1F* el1_Iso      = new TH1F("el1_Iso", "el1_Iso", Iso_nbins, Iso_min, Iso_max );
TH2F* el1_phi_d0   = new TH2F("el1_phi_d0", "el1_phi_d0", phi_nbins, phi_min, phi_max, d0_nbins, d0_min, d0_max );

// met & met phi - before
TH1F* el1_clmet        = new TH1F("el1_clmet", "el1_clmet", met_nbins, met_min, met_max );
TH1F* el1_pfmet        = new TH1F("el1_pfmet", "el1_pfmet", met_nbins, met_min, met_max );
TH1F* el1_tcmet        = new TH1F("el1_tcmet", "el1_tcmet", met_nbins, met_min, met_max );
TH1F* el1_clmetphi     = new TH1F("el1_clmetphi", "el1_clmetphi", phi_nbins, phi_min, phi_max );
TH1F* el1_pfmetphi     = new TH1F("el1_pfmetphi", "el1_pfmetphi", phi_nbins, phi_min, phi_max );
TH1F* el1_tcmetphi     = new TH1F("el1_tcmetphi", "el1_tcmetphi", phi_nbins, phi_min, phi_max );

// jets - before

TH1F* el1_Njets     = new TH1F("el1_Njets", "el1_Ncalojets", Njet_nbins, Njet_min, Njet_max );
TH1F* el1_Npfjets   = new TH1F("el1_Npfjets", "el1_Npfjets", Njet_nbins, Njet_min, Njet_max );
TH1F* el1_Ntrkjets  = new TH1F("el1_Ntrkjets", "el1_Ntrkjets", Njet_nbins, Njet_min, Njet_max );

TH1F* el1_jets_pt   	= new TH1F("el1_jets_pt", "el1_calojets_pt", jet_pt_nbins, jet_pt_min, jet_pt_max );
TH1F* el1_pfjets_pt   	= new TH1F("el1_pfjets_pt", "el1_pfjets_pt", jet_pt_nbins, jet_pt_min, jet_pt_max );
TH1F* el1_trkjets_pt   	= new TH1F("el1_trkjets_pt", "el1_trkjets_pt", jet_pt_nbins, jet_pt_min, jet_pt_max );

TH1F* el1_jets_eta   	= new TH1F("el1_jets_eta", "el1_calojets_eta", jet_eta_nbins, jet_eta_min, jet_eta_max );
TH1F* el1_pfjets_eta   	= new TH1F("el1_pfjets_eta", "el1_pfjets_eta", jet_eta_nbins, jet_eta_min, jet_eta_max );
TH1F* el1_trkjets_eta   = new TH1F("el1_trkjets_eta", "el1_trkjets_eta", jet_eta_nbins, jet_eta_min, jet_eta_max );

TH1F* el1_jets_phi   	= new TH1F("el1_jets_phi", "el1_calojets_phi", phi_nbins, phi_min, phi_max );
TH1F* el1_pfjets_phi   	= new TH1F("el1_pfjets_phi", "el1_pfjets_phi", phi_nbins, phi_min, phi_max );
TH1F* el1_trkjets_phi   = new TH1F("el1_trkjets_phi", "el1_trkjets_phi", phi_nbins, phi_min, phi_max );

// elons - after
TH1F* el2_N        = new TH1F("el2_Nels", "el2_Nels", N_nbins, N_min, N_max );
TH1F* el2_Niso     = new TH1F("el2_Nelsiso", "el2_Nels_isolated", N_nbins, N_min, N_max );
TH1F* el2_pt       = new TH1F("el2_pt", "el2_pt", pt_nbins, pt_min, pt_max );
TH1F* el2_eta      = new TH1F("el2_eta", "el2_eta", eta_nbins, eta_min, eta_max );
TH1F* el2_phi      = new TH1F("el2_phi", "el2_phi", phi_nbins, phi_min, phi_max );
TH1F* el2_d0corr   = new TH1F("el2_d0corr", "el2_d0corr", d0_nbins, d0_min, d0_max );
TH1F* el2_Iso      = new TH1F("el2_Iso", "el2_Iso", Iso_nbins, Iso_min, Iso_max );
TH2F* el2_phi_d0   = new TH2F("el2_phi_d0", "el2_phi_d0", phi_nbins, phi_min, phi_max, d0_nbins, d0_min, d0_max );

// met & met phi - after
TH1F* el2_clmet        = new TH1F("el2_clmet", "el2_clmet", met_nbins, met_min, met_max );
TH1F* el2_pfmet        = new TH1F("el2_pfmet", "el2_pfmet", met_nbins, met_min, met_max );
TH1F* el2_tcmet        = new TH1F("el2_tcmet", "el2_tcmet", met_nbins, met_min, met_max );
TH1F* el2_clmetphi     = new TH1F("el2_clmetphi", "el2_clmetphi", phi_nbins, phi_min, phi_max );
TH1F* el2_pfmetphi     = new TH1F("el2_pfmetphi", "el2_pfmetphi", phi_nbins, phi_min, phi_max );
TH1F* el2_tcmetphi     = new TH1F("el2_tcmetphi", "el2_tcmetphi", phi_nbins, phi_min, phi_max );

// jets - after

TH1F* el2_Njets     = new TH1F("el2_Njets", "el2_Ncalojets", Njet_nbins, Njet_min, Njet_max );
TH1F* el2_Npfjets   = new TH1F("el2_Npfjets", "el2_Npfjets", Njet_nbins, Njet_min, Njet_max );
TH1F* el2_Ntrkjets  = new TH1F("el2_Ntrkjets", "el2_Ntrkjets", Njet_nbins, Njet_min, Njet_max );

TH1F* el2_jets_pt   	= new TH1F("el2_jets_pt", "el2_calojets_pt", jet_pt_nbins, jet_pt_min, jet_pt_max );
TH1F* el2_pfjets_pt   	= new TH1F("el2_pfjets_pt", "el2_pfjets_pt", jet_pt_nbins, jet_pt_min, jet_pt_max );
TH1F* el2_trkjets_pt   	= new TH1F("el2_trkjets_pt", "el2_trkjets_pt", jet_pt_nbins, jet_pt_min, jet_pt_max );
TH1F* el2_jets_eta   	= new TH1F("el2_jets_eta", "el2_calojets_eta", jet_eta_nbins, jet_eta_min, jet_eta_max );
TH1F* el2_pfjets_eta   	= new TH1F("el2_pfjets_eta", "el2_pfjets_eta", jet_eta_nbins, jet_eta_min, jet_eta_max );
TH1F* el2_trkjets_eta   = new TH1F("el2_trkjets_eta", "el2_trkjets_eta", jet_eta_nbins, jet_eta_min, jet_eta_max );
TH1F* el2_jets_phi   	= new TH1F("el2_jets_phi", "el2_calojets_phi", phi_nbins, phi_min, phi_max );
TH1F* el2_pfjets_phi   	= new TH1F("el2_pfjets_phi", "el2_pfjets_phi", phi_nbins, phi_min, phi_max );
TH1F* el2_trkjets_phi   = new TH1F("el2_trkjets_phi", "el2_trkjets_phi", phi_nbins, phi_min, phi_max );

// Chains
TChain *chain1 = new TChain("tree");
TChain *chain2 = new TChain("tree");
chain1->Add("/tas03/home/warren/CMSSW_3_7_0_patch2/src/CMS2/NtupleMacros/SkimValidation/validate_els_before.root");
chain2->Add("/tas03/home/warren/CMSSW_3_7_0_patch2/src/CMS2/NtupleMacros/SkimValidation/validate_els_after.root");

// Fill Before
TCanvas *ctemp = new TCanvas();
chain1->Draw("nel >> el1_Nels");
chain1->Draw("Sum$(elsid) >> el1_Nelsiso"); //isolated bc elsid is isolated
chain1->Draw("elsp4.pt() >> el1_pt", "elsid");
chain1->Draw("elsp4.eta() >> el1_eta", "elsid");
chain1->Draw("elsp4.phi() >> el1_phi", "elsid");
chain1->Draw("elsd0corr >> el1_d0corr", "elsid");
chain1->Draw("elsiso >> el1_Iso");
chain1->Draw("elsd0corr:elsp4.phi() >> el1_phi_d0", "elsid", "BOX");

chain1->Draw("clmet >> el1_clmet");
chain1->Draw("pfmet >> el1_pfmet");
chain1->Draw("tcmet >> el1_tcmet");

chain1->Draw("clmetphi >> el1_clmetphi", "clmet>10");
chain1->Draw("pfmetphi >> el1_pfmetphi", "pfmet>10");
chain1->Draw("tcmetphi >> el1_tcmetphi", "tcmet>10");

chain1->Draw("jets@.size() >> el1_Njets");
chain1->Draw("pfjets@.size() >> el1_Npfjets");
chain1->Draw("trkjets@.size() >> el1_Ntrkjets");

chain1->Draw("jets.pt() >> el1_jets_pt");
chain1->Draw("pfjets.pt() >> el1_pfjets_pt");
chain1->Draw("trkjets.pt() >> el1_trkjets_pt");
chain1->Draw("jets.eta() >> el1_jets_eta");
chain1->Draw("pfjets.eta() >> el1_pfjets_eta");
chain1->Draw("trkjets.eta() >> el1_trkjets_eta");
chain1->Draw("jets.phi() >> el1_jets_phi");
chain1->Draw("pfjets.phi() >> el1_pfjets_phi");
chain1->Draw("trkjets.phi() >> el1_trkjets_phi");

// Fill After
chain2->Draw("nel >> el2_Nels");
chain2->Draw("Sum$(elsid) >> el2_Nelsiso");
chain2->Draw("elsp4.pt() >> el2_pt", "elsid");
chain2->Draw("elsp4.eta() >> el2_eta", "elsid");
chain2->Draw("elsp4.phi() >> el2_phi", "elsid");
chain2->Draw("elsd0corr >> el2_d0corr", "elsid");
chain2->Draw("elsiso >> el2_Iso");
chain2->Draw("elsd0corr:elsp4.phi() >> el2_phi_d0", "elsid", "BOX PROJ");

chain2->Draw("clmet >> el2_clmet");
chain2->Draw("pfmet >> el2_pfmet");
chain2->Draw("tcmet >> el2_tcmet");

chain2->Draw("clmetphi >> el2_clmetphi", "clmet>10");
chain2->Draw("pfmetphi >> el2_pfmetphi", "pfmet>10");
chain2->Draw("tcmetphi >> el2_tcmetphi", "tcmet>10");

chain2->Draw("jets@.size() >> el2_Njets");
chain2->Draw("pfjets@.size() >> el2_Npfjets");
chain2->Draw("trkjets@.size() >> el2_Ntrkjets");

chain2->Draw("jets.pt() >> el2_jets_pt");
chain2->Draw("pfjets.pt() >> el2_pfjets_pt");
chain2->Draw("trkjets.pt() >> el2_trkjets_pt");
chain2->Draw("jets.eta() >> el2_jets_eta");
chain2->Draw("pfjets.eta() >> el2_pfjets_eta");
chain2->Draw("trkjets.eta() >> el2_trkjets_eta");
chain2->Draw("jets.phi() >> el2_jets_phi");
chain2->Draw("pfjets.phi() >> el2_pfjets_phi");
chain2->Draw("trkjets.phi() >> el2_trkjets_phi");

delete ctemp;

TCanvas *c1 = new TCanvas();
c1->SetWindowSize(1100,850);
c1->Divide(3,2);
c1->cd(1)->SetLogy();
el1_pt->Draw();
c1->cd(2);
el1_eta->Draw();
c1->cd(3);
el1_phi->Draw();
c1->cd(4)->SetLogy();
el2_pt->Draw();
c1->cd(5);
el2_eta->Draw();
c1->cd(6);
el2_phi->Draw();
c1->SaveAs((dir+"compare_el_kinematics"+suffix).c_str());

TCanvas *c2 = new TCanvas();
 c2->SetWindowSize(1450,850);//wider bc 4 plots
c2->Divide(4,2);
c2->cd(1)->SetLogy();
el1_N->Draw();
c2->cd(2)->SetLogy();
el1_Niso->Draw();
c2->cd(3);
el1_d0corr->Draw();
c2->cd(4);
el1_Iso->Draw();
c2->cd(5)->SetLogy();
el2_N->Draw();
c2->cd(6)->SetLogy();
el2_Niso->Draw();
c2->cd(7);
el2_d0corr->Draw();
c2->cd(8);
el2_Iso->Draw();
c1->SaveAs((dir+"compare_el_N_do_iso"+suffix).c_str());

TCanvas *c3 = new TCanvas();
c3->SetWindowSize(1100,900);
c3->Divide(2,2);
c3->cd(1);
el1_phi_d0->Draw("BOX");
c3->cd(2);
el2_phi_d0->Draw("BOX");
c3->cd(3);
el1_phi_d0->ProfileX()->Draw();
c3->cd(4);
el2_phi_d0->ProfileX()->Draw();
c1->SaveAs((dir+"compare_el_2dkin"+suffix).c_str());

TCanvas *c4 = new TCanvas();
c4->SetWindowSize(1100,850);
c4->Divide(3,2);
c4->cd(1)->SetLogy();
el1_clmet->Draw();
c4->cd(2)->SetLogy();
el1_pfmet->Draw();
c4->cd(3)->SetLogy();
el1_tcmet->Draw();
c4->cd(4)->SetLogy();
el2_clmet->Draw();
c4->cd(5)->SetLogy();
el2_pfmet->Draw();
c4->cd(6)->SetLogy();
el2_tcmet->Draw();
c1->SaveAs((dir+"compare_el_met"+suffix).c_str());

TCanvas *c5 = new TCanvas();
c5->SetWindowSize(1100,850);
c5->Divide(3,2);
c5->cd(1);
el1_clmetphi->Draw();
c5->cd(2);
el1_pfmetphi->Draw();
c5->cd(3);
el1_tcmetphi->Draw();
c5->cd(4);
el2_clmetphi->Draw();
c5->cd(5);
el2_pfmetphi->Draw();
c5->cd(6);
el2_tcmetphi->Draw();
c1->SaveAs((dir+"compare_el_metphi"+suffix).c_str());

TCanvas *c = new TCanvas();
c->SetWindowSize(1100,850);
c->Divide(3,2);
c->cd(1);
el1_Njets->Draw();
c->cd(2);
el1_Npfjets->Draw();
c->cd(3);
el1_Ntrkjets->Draw();
c->cd(4);
el2_Njets->Draw();
c->cd(5);
el2_Npfjets->Draw();
c->cd(6);
el2_Ntrkjets->Draw();
c1->SaveAs((dir+"compare_el_jet"+suffix).c_str());

TCanvas *c6 = new TCanvas();
c6->SetWindowSize(1100,850);
c6->Divide(3,2);
c6->cd(1);
el1_jets_pt->Draw();
c6->cd(2);
el1_jets_eta->Draw();
c6->cd(3);
el1_jets_phi->Draw();
c6->cd(4);
el2_jets_pt->Draw();
c6->cd(5);
el2_jets_eta->Draw();
c6->cd(6);
el2_jets_phi->Draw();
c1->SaveAs((dir+"compare_el_cljetkin"+suffix).c_str());

TCanvas *c7 = new TCanvas();
c7->SetWindowSize(1100,850);
c7->Divide(3,2);
c7->cd(1);
el1_pfjets_pt->Draw();
c7->cd(2);
el1_pfjets_eta->Draw();
c7->cd(3);
el1_pfjets_phi->Draw();
c7->cd(4);
el2_pfjets_pt->Draw();
c7->cd(5);
el2_pfjets_eta->Draw();
c7->cd(6);
el2_pfjets_phi->Draw();
c1->SaveAs((dir+"compare_el_pfjetkin"+suffix).c_str());

TCanvas *c8 = new TCanvas();
c8->SetWindowSize(1100,850);
c8->Divide(3,2);
c8->cd(1);
el1_trkjets_pt->Draw();
c8->cd(2);
el1_trkjets_eta->Draw();
c8->cd(3);
el1_trkjets_phi->Draw();
c8->cd(4);
el2_trkjets_pt->Draw();
c8->cd(5);
el2_trkjets_eta->Draw();
c8->cd(6);
el2_trkjets_phi->Draw();
c1->SaveAs((dir+"compare_el_tkjetkin"+suffix).c_str());

}
