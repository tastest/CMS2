
{

    TFile f_signal("../baby_lm0.root", "READ");
    TFile f_bg("../baby_wjets.root", "READ");

    TTree *t_signal = (TTree*)f_signal.Get("tree");
    TTree *t_bg = (TTree*)f_bg.Get("tree");

    TH1F *h1_signal_pfmva = new TH1F("h1_signal_pfmva", "signal_pfmva", 100, -1.0, 1.0);
    TH1F *h1_bg_pfmva = new TH1F("h1_bg_pfmva", "bg_pfmva", 100, -1.0, 1.0);

    TCut denominator("abs(reco_typelt) == 11 && reco_ptlt > 15.0");
    TCut denominator_signal = denominator + TCut("abs(reco_mctypelt) == 11");
    TCut denominator_bg = denominator + TCut("abs(reco_mctypelt) != 11");

    t_signal->Draw("reco_pfmvalt >> h1_signal_pfmva", "weight * " + denominator_signal);
    t_bg->Draw("reco_pfmvalt >> h1_bg_pfmva", "weight * " + denominator_bg);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    h1_signal_pfmva->Draw();
    h1_bg_pfmva->Draw("SAME");

}

