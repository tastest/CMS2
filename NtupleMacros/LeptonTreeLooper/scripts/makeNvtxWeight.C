{

    TChain *ch_data = new TChain("leptons");
    ch_data->Add("/smurf/dlevans/LeptonTree/V00-01-07/DoubleMuRun2012APromptV1/merged_Cert_190456-191859_8TeV_PromptReco_Collisions12_JSON.root");

    TChain *ch_dyee = new TChain("leptons");
    ch_dyee->Add("/smurf/dlevans/LeptonTree/V00-01-07/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/merged.root");

    gROOT->cd();

    TFile outfile("nvtxweight.root", "RECREATE");
    TH1F *h1_nvtx_data = new TH1F("h1_nvtx_data", "h1_nvtx_data", 20, -0.5, 39.5);
    TH1F *h1_nvtx_mc = new TH1F("h1_nvtx_mc", "h1_nvtx_mc", 20, -0.5, 39.5);

    ch_dyee->Draw("TMath::Min(39, nvtx) >> h1_nvtx_mc",    "eventSelection & (1<<0)");
    ch_data->Draw("TMath::Min(39, nvtx) >> h1_nvtx_data",  "eventSelection & (1<<0)");

    h1_nvtx_data->Scale(1.0/h1_nvtx_data->Integral(0, 20));
    h1_nvtx_mc->Scale(1.0/h1_nvtx_mc->Integral(0, 20));

    TH1F *h1_weight = (TH1F*)h1_nvtx_data->Clone("weight");
    h1_weight->Divide(h1_nvtx_mc);

    outfile.cd();
    h1_weight->Write();
    outfile.Close();


}


