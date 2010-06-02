
{

    TFile f_mc("../histos_mc_spike.root", "READ");
    TFile f_data("../histos_data_spike.root", "READ");

    TH2F *h2_scatteret_mc = (TH2F*)f_mc.Get("h2_spike_scatteret");
    h2_scatteret_mc->SetMarkerStyle(5);
    TH2F *h2_scatteret_data = (TH2F*)f_data.Get("h2_spike_scatteret");
    h2_scatteret_data->SetMarkerColor(kRed);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    h2_scatteret_mc->Draw();
    h2_scatteret_mc->GetXaxis()->SetRangeUser(0, 0.15);
    h2_scatteret_data->Draw("SAMES");

    TH2F *h2_scattermet_mc = (TH2F*)f_mc.Get("h2_spike_scattermet");
    h2_scattermet_mc->SetMarkerStyle(5);
    TH2F *h2_scattermet_data = (TH2F*)f_data.Get("h2_spike_scattermet");
    h2_scattermet_data->SetMarkerColor(kRed);

    TCanvas *c2 = new TCanvas();
    c2->cd();
    h2_scattermet_mc->Draw();
    h2_scattermet_mc->GetXaxis()->SetRangeUser(0, 0.15);
    h2_scattermet_data->Draw("SAMES");

    TH2F *h2_scattermet_fixed_mc = (TH2F*)f_mc.Get("h2_spike_scattermet_fixed");
    h2_scattermet_fixed_mc->SetMarkerStyle(5);
    TH2F *h2_scattermet_fixed_data = (TH2F*)f_data.Get("h2_spike_scattermet_fixed");
    h2_scattermet_fixed_data->SetMarkerColor(kRed);

    TCanvas *c3 = new TCanvas();
    c3->cd();
    h2_scattermet_fixed_mc->Draw();
    h2_scattermet_fixed_mc->GetXaxis()->SetRangeUser(0, 0.15);
    h2_scattermet_fixed_data->Draw("SAMES");

    TH2F *h2_scatteret_eid_mc = (TH2F*)f_mc.Get("h2_spike_scatteret_eid");
    h2_scatteret_eid_mc->SetMarkerStyle(5);
    TH2F *h2_scatteret_eid_data = (TH2F*)f_data.Get("h2_spike_scatteret_eid");
    h2_scatteret_eid_data->SetMarkerColor(kRed);

    TCanvas *c4 = new TCanvas();
    c4->cd();
    h2_scatteret_eid_mc->Draw();
    h2_scatteret_eid_mc->GetXaxis()->SetRangeUser(0, 0.15);
    h2_scatteret_eid_data->Draw("SAMES");



/*
    gSystem->Load("../../../Tools/MiniFWLite/libMiniFWLite.so");

    TChain *chain_may6eg = new TChain("Events");
    chain_may6eg->Add("/tas07/disk00/cms2/MinimumBias_Commissioning10-May6thPDSkim2_SD_EG-v1/V03-04-09/merged*.root");
    TChain *chain_qcd30 = new TChain("Events");
    chain_qcd30->Add("/tas07/disk00/cms2/QCD_Pt30_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*.root");

    TString cut = "(els_type & (1<<2)) > 0 && abs(els_etaSC) < 1.5 && (els_eSC/cosh(els_etaSC)) > 10";
    gStyle->SetOptStat(111111);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    TH2F *h2_scatter10_data = new TH2F("h2_scatter10_data", "h2_scatter10_data;eMax/e5x5;tcMET", 150, 0, 1.5, 100, 0, 100);
    h2_scatter10_data->SetMarkerColor(kRed);
    h2_scatter10_data->SetLineColor(kRed);
    TH2F *h2_scatter10_mc = new TH2F("h2_scatter10_mc", "h2_scatter10_mc;eMax/e5x5;tcMET", 150, 0, 1.5, 100, 0, 100);
    chain_may6eg->Draw("evt_tcmet:(els_eMax/els_e5x5)>>h2_scatter10_data", cut + "&& (els_eSC/cosh(els_etaSC)) > 10");
    chain_qcd30->Draw("evt_tcmet:(els_eMax/els_e5x5)>>h2_scatter10_mc", cut + "&& (els_eSC/cosh(els_etaSC)) > 10");
    h2_scatter10_mc->Draw("HIST");
    h2_scatter10_data->Draw("SAMES");

    TCanvas *c2 = new TCanvas();
    c2->cd();
    TH2F *h2_scatter20_data = new TH2F("h2_scatter20_data", "h2_scatter20_data;eMax/e5x5;tcMET", 150, 0, 1.5, 100, 0, 100);
    h2_scatter20_data->SetMarkerColor(kRed);
    h2_scatter20_data->SetLineColor(kRed);
    TH2F *h2_scatter20_mc = new TH2F("h2_scatter20_mc", "h2_scatter20_mc;eMax/e5x5;tcMET", 150, 0, 1.5, 100, 0, 100); 
    chain_may6eg->Draw("evt_tcmet:(els_eMax/els_e5x5)>>h2_scatter20_data", cut + "&& (els_eSC/cosh(els_etaSC)) > 20");
    chain_qcd30->Draw("evt_tcmet:(els_eMax/els_e5x5)>>h2_scatter20_mc", cut + "&& (els_eSC/cosh(els_etaSC)) > 20");
    h2_scatter20_mc->Draw("HIST");
    h2_scatter20_data->Draw("SAMES");

    TCanvas *c3 = new TCanvas();
    c3->cd();
    TH2F *h2_scatter30_data = new TH2F("h2_scatter30_data", "h2_scatter30_data;eMax/e5x5;tcMET", 150, 0, 1.5, 100, 0, 100);
    h2_scatter30_data->SetMarkerColor(kRed);
    h2_scatter30_data->SetLineColor(kRed);
    TH2F *h2_scatter30_mc = new TH2F("h2_scatter30_mc", "h2_scatter30_mc;eMax/e5x5;tcMET", 150, 0, 1.5, 100, 0, 100); 
    chain_may6eg->Draw("evt_tcmet:(els_eMax/els_e5x5) >> h2_scatter30_data", cut + "&& (els_eSC/cosh(els_etaSC)) > 30");
    chain_qcd30->Draw("evt_tcmet:(els_eMax/els_e5x5) >> h2_scatter30_mc", cut + "&& (els_eSC/cosh(els_etaSC)) > 30");
    h2_scatter30_mc->Draw("HIST");
    h2_scatter30_data->Draw("SAMES");

    TCanvas *c4 = new TCanvas();
    c4->cd();
    TH2F *h2_scatteret_data = new TH2F("h2_scatteret_data", "h2_scatteret_data;eMax/e5x5;ET_{SC}", 150, 0, 1.5, 200, 0, 200);
    h2_scatteret_data->SetMarkerColor(kRed);
    h2_scatteret_data->SetLineColor(kRed);
    TH2F *h2_scatteret_mc = new TH2F("h2_scatteret_mc", "h2_scatteret_mc;eMax/e5x5;ET_{SC}", 150, 0, 1.5, 200, 0, 200);
    chain_may6eg->Draw("(els_eSC/cosh(els_etaSC)):(els_eMax/els_e5x5) >> h2_scatteret_data", cut);
    chain_qcd30->Draw("(els_eSC/cosh(els_etaSC)):(els_eMax/els_e5x5) >> h2_scatteret_mc", cut);
    h2_scatteret_mc->Draw("HIST");
    h2_scatteret_data->Draw("SAMES");

    //  
    // look at ele id criteria
    //

    TCanvas *c5 = new TCanvas();
    c5->cd();
    TH2F *h2_scattermatch_data = new TH2F("h2_scattermatch_data", "h2_scattermatch_data;dEtaIn;dPhiIn", 200, -0.1, 0.1, 200, -0.1, 0.1);
    h2_scattermatch_data->SetMarkerColor(kRed);
    h2_scattermatch_data->SetLineColor(kRed);
    TH2F *h2_scattermatch_mc = new TH2F("h2_scattermatch_mc", "h2_scattermatch_mc;dEtaIn;dPhiIn", 200, -0.1, 0.1, 200, -0.1, 0.1);
    chain_may6eg->Draw("els_dEtaIn:els_dPhiIn >> h2_scattermatch_data", cut + "&& (els_eMax/els_e5x5) > 0.95");
    chain_qcd30->Draw("els_dEtaIn:els_dPhiIn >> h2_scattermatch_mc", cut);
    h2_scattermatch_mc->Draw("BOX");
    h2_scattermatch_data->Draw("SAMES");

    TCanvas *c6 = new TCanvas();
    c6->cd();
    TH1F *h1_hoe_data = new TH1F("h1_hoe_data", "h1_hoe_data;H/E", 100, 0, 0.5);
    h1_hoe_data->SetMarkerColor(kRed);
    h1_hoe_data->SetLineColor(kRed);
    h1_hoe_data->SetFillColor(kRed);
    chain_may6eg->Draw("els_hOverE >> h1_hoe_data", cut + "&& (els_eMax/els_e5x5) > 0.95");
    h1_hoe_data->Draw("HIST");

    TCanvas *c7 = new TCanvas();
    c7->cd();
    TH1F *h1_sigmaIEtaIEta_data = new TH1F("h1_sigmaIEtaIEta_data", "h1_sigmaIEtaIEta_data;sigmaIEtaIEta", 100, 0, 0.05);
    h1_sigmaIEtaIEta_data->SetMarkerColor(kRed);
    h1_sigmaIEtaIEta_data->SetLineColor(kRed);
    h1_sigmaIEtaIEta_data->SetFillColor(kRed);
    chain_may6eg->Draw("els_sigmaIEtaIEta >> h1_sigmaIEtaIEta_data", cut + "&& (els_eMax/els_e5x5) > 0.95");
    h1_sigmaIEtaIEta_data->Draw("HIST");

    TCanvas *c8 = new TCanvas();
    c8->cd();
    TH1F *h1_eopin_data = new TH1F("h1_eopin_data", "h1_eopin_data;E/p_{IN}", 100, 0, 10.0);
    h1_eopin_data->SetMarkerColor(kRed);
    h1_eopin_data->SetLineColor(kRed);
    h1_eopin_data->SetFillColor(kRed);
    chain_may6eg->Draw("els_eOverPIn >> h1_eopin_data", cut + "&& (els_eMax/els_e5x5) > 0.95");
    h1_eopin_data->Draw("HIST");
*/
}


