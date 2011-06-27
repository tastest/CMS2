
void format(TH1F *h1, bool data, unsigned int rebin, float xmin, float xmax)
{

    h1->Rebin(rebin);
    h1->GetXaxis()->SetRangeUser(xmin, xmax);

    if (data) {
        h1->SetMarkerStyle(20);
    }
    else {
        h1->SetFillStyle(3002);
        h1->SetFillColor(kBlue);
        h1->SetLineColor(kBlue);
    }

}

void overlay(TString data, TString title, TString titleX, TString name, TString type)
{

    TFile f("histos_baby.root");
    gROOT->cd();
    gStyle->SetOptTitle(1);
 
    // get norm for data-mc
    // DoubleElectronRun2011AApr22ReReco_h1_hyp_pt_etaincl_allj_ee
    TH1F *h1_ptnorm_data = (TH1F*)f.Get(TString(data+"_h1_hyp_pt_etaincl_allj_"+type))->Clone();
    TH1F *h1_ptnorm_mc = (TH1F*)f.Get(TString("dyeemm_h1_hyp_pt_etaincl_allj_"+type))->Clone();
    float n_data = h1_ptnorm_data->Integral(0, h1_ptnorm_data->FindBin(40.0));
    float n_mc = h1_ptnorm_mc->Integral(0, h1_ptnorm_mc->FindBin(40.0));
    float scale = n_data/n_mc;

    // get histogram
    TH1F *h1_data = (TH1F*)f.Get(TString(data+"_h1_" + name + "_" + type))->Clone();
    TH1F *h1_mc = (TH1F*)f.Get(TString("dyeemm_h1_" + name + "_" + type))->Clone();
    format(h1_data, true, 1, 40.0, 200.0);
    format(h1_mc, false, 1, 40.0, 200.0);

    // scale mc to data
    h1_mc->Scale(scale);

    // draw
    TLegend *l1 = new TLegend(0.6, 0.8, 0.9, 0.9);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->SetShadowColor(0);
    l1->AddEntry(h1_mc, "Madgraph DY", "f");
    l1->AddEntry(h1_data, data, "lp");

    TCanvas *c1 = new TCanvas();
    c1->cd();
    h1_data->Draw("E1");
    h1_mc->Draw("SAME HIST");
    h1_data->SetTitle(title);
    h1_data->GetXaxis()->SetTitle(titleX);
    l1->Draw();
    c1->SaveAs(data+"_"+type+"_"+name+".png");


    // tidy up
    f.Close();
    delete h1_ptnorm_data;
    delete h1_ptnorm_mc;
    delete h1_data;
    delete h1_mc;

}

void plot(TString data, TString hyp)
{

    // pt
    overlay(data, "Inclusive", "dilepton p_{T}", "hyp_pt_etaincl_allj", hyp);
    overlay(data, "Inclusive, |#eta_{Z}| < 1", "dilepton p_{T}", "hyp_pt_eta1_allj", hyp);
    overlay(data, "1-Jet", "dilepton p_{T}", "hyp_pt_etaincl_1j", hyp);
    overlay(data, "1-Jet, |#eta_{Z}| < 1", "dilepton p_{T}", "hyp_pt_eta1_1j", hyp);
    // njets in "signal" region
    overlay(data, "80 < dilepton p_{T} < 100", "#Jets", "hyp_njets_sig_etaincl", hyp);
    overlay(data, "80 < dilepton p_{T} < 100, |#eta_{Z}| < 1", "#Jets", "hyp_njets_sig_eta1", hyp);
    // delta-phi between leptons in "signal" region
    overlay(data, "Inclusive, 80 < dilepton p_{T} < 100", "#Delta#phi", "hyp_dphi_sig_etaincl_allj", hyp);
    overlay(data, "Inclusive, 80 < dilepton p_{T} < 100, |#eta_{Z}| < 1", "#Delta#phi", "hyp_dphi_sig_eta1_allj", hyp);
    overlay(data, "1-Jet, 80 < dilepton p_{T} < 100", "#Delta#phi", "hyp_dphi_sig_etaincl_1j", hyp);
    overlay(data, "1-Jet 80 < dilepton p_{T} < 100, |#eta_{Z}| < 1", "#Delta#phi", "hyp_dphi_sig_eta1_1j", hyp);

}

