
int bin_njets = 2;

void printIntegral(std::string text, const TH1F *histSignal, const TH1F *histBackground)
{

    float s = histSignal->Integral(bin_njets, histSignal->GetNbinsX() + 1);
    float b = histBackground->Integral(bin_njets, histSignal->GetNbinsX() + 1);

    std::cout << text << " \t& ";
    printf("%4.2f", s);
    std::cout << " \t& ";
    printf("%4.2f", b);
    std::cout << " \t& ";
    printf("%4.2f", s/sqrt(b));
    std::cout << " \\\\ \\hline" << std::endl;

}

void formatHist(TH1F *hist, bool signal) {

    hist->SetDirectory(gDirectory);
    hist->SetLineWidth(2);

        if (signal) {
            hist->SetFillColor(kYellow);
        }
        else {
            hist->SetLineColor(kRed);
            hist->SetMarkerStyle(20);
            hist->SetMarkerColor(kRed);
        }
}

void setError(TH1F *h1_numer, TH1F *h1_denom, TH1F *h1_eff) {

    Float_t scale = h1_numer->GetEntries() / h1_numer->Integral(0, h1_numer->GetNbinsX() + 1);
    for (Int_t i = 0; i < h1_numer->GetNbinsX() + 1; ++i) {
        Float_t nNumer = h1_numer->GetBinContent(i) * scale;
        Float_t nDenom = h1_denom->GetBinContent(i) * scale;
        Float_t eff = h1_eff->GetBinContent(i);
        Float_t err = 0.0;
        if (nDenom != 0) {
            err = sqrt( eff * (1 - eff) / nDenom);
        }
        h1_eff->SetBinError(i, err);
    }

}

void makeQuickTables() {


    gROOT->ProcessLine(".L ../tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    TFile f("../histos_mc.root", "READ");
    gROOT->cd();    

    //
    // just sani loose
    //
    TH1F *ttbar_hyp_njets_ee = (TH1F*)f.Get("ttbar_hyp_njets_ee");
    formatHist(ttbar_hyp_njets_ee, true);
    TH1F *wjets_hyp_njets_ee = (TH1F*)f.Get("wjets_hyp_njets_ee");
    formatHist(wjets_hyp_njets_ee, true);

    TH1F *ttbar_hyp_njets_em = (TH1F*)f.Get("ttbar_hyp_njets_em");
    formatHist(ttbar_hyp_njets_em, true);
    TH1F *wjets_hyp_njets_em = (TH1F*)f.Get("wjets_hyp_njets_em");
    formatHist(wjets_hyp_njets_em, true);

    TH1F *ttbar_hyp_njets_mm = (TH1F*)f.Get("ttbar_hyp_njets_mm");
    formatHist(ttbar_hyp_njets_mm, true);
    TH1F *wjets_hyp_njets_mm = (TH1F*)f.Get("wjets_hyp_njets_mm");
    formatHist(wjets_hyp_njets_mm, true);

    //
    // pt distributions before iso loose
    //

    TCanvas *c1 = new TCanvas();
    c1->cd();

    TLegend *l1 = new TLegend(0.5, 0.5, 0.9, 0.9);
    l1->SetFillColor(kWhite);
    l1->SetLineColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(ttbar_hyp_njets_ee, "TTBar", "lf");
    l1->AddEntry(wjets_hyp_njets_ee, "W+Jets", "lf");


    ttbar_hyp_njets_ee->Draw("HIST");
    wjets_hyp_njets_ee->Draw("HIST SAME E1");
    l1->Draw();

    c1->SaveAs("../results/study_njets.png");

    // print table
    std::cout << "\\begin{table}[ht]" << std::endl;
    std::cout << "\\caption{Yields for two or more jets}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|l*{4}{|c}|r|}\\hline" << std::endl;
    std::cout << "Final state   & S & B & S/$\\sqrt{B}$  \\\\ \\hline" << std::endl;
    printIntegral("EE", ttbar_hyp_njets_ee, wjets_hyp_njets_ee);
    printIntegral("EM", ttbar_hyp_njets_em, wjets_hyp_njets_em);
    printIntegral("EM", ttbar_hyp_njets_mm, wjets_hyp_njets_mm);
    std::cout <<"\\end{tabular}" << std::endl;
    std::cout <<"\\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;

}

