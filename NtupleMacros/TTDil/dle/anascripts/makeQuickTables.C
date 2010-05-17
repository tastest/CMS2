
void printIntegral(std::string text, const TH1F *histSignalee, const TH1F *histBackgroundee, 
                                        const TH1F *histSignalem, const TH1F *histBackgroundem, int njets)
{

    Double_t errsee = 0.0;
    Double_t errbee = 0.0;

    Double_t see = histSignalee->IntegralAndError(njets + 1, histSignalee->GetNbinsX() + 1, errsee);
    Double_t bee = histBackgroundee->IntegralAndError(njets + 1, histSignalee->GetNbinsX() + 1, errbee);

    Double_t errsem = 0.0;
    Double_t errbem = 0.0;

    Double_t sem = histSignalem->IntegralAndError(njets + 1, histSignalem->GetNbinsX() + 1, errsem);
    Double_t bem = histBackgroundem->IntegralAndError(njets + 1, histSignalem->GetNbinsX() + 1, errbem);

    Double_t s = see + sem;
    Double_t b = bee + bem;

    Double_t errs = sqrt ( pow(errsee, 2) + pow(errsem, 2) );
    Double_t errb = 0.0;
    if (bee != 0 && bem != 0) errb = sqrt ( pow(errbee, 2) + pow(errbem, 2) );
    if (bee == 0 && bem != 0) errb = errbem;
    if (bee != 0 && bem == 0) errb = errbee;


    //std::cout << text << " \t& ";
    //printf("%4.2f", s);
    //std::cout << " \t& ";
    //printf("%4.2f", b);
    //std::cout << " \t& ";
    //if (b > 0)
    //    printf("%4.2f", s/sqrt(b));
    //else
    //    printf("N/A");
    //std::cout << " \\\\ \\hline" << std::endl;

    //if (text == "ALL") {
        std::cout << "gr->SetPoint(P, " << b << ", " << s << ");" << std::endl;
        std::cout << "gr->SetPointError(P, " << errb << ", " << errs << ");" << std::endl;
    //}

}

void printIntegral(std::string text, const TH1F *histSignal, const TH1F *histBackground, int njets)
{


    //Double_t IntegralAndError(Int_t binx1, Int_t binx2, Double_t& err, Option_t* option = "") const
    Double_t errs = 0.0;
    Double_t errb = 0.0;

    Double_t s = histSignal->IntegralAndError(njets + 1, histSignal->GetNbinsX() + 1, errs);
    Double_t b = histBackground->IntegralAndError(njets + 1, histSignal->GetNbinsX() + 1, errb);

    std::cout << text << " \t& ";
    printf("%4.2f", s);
    std::cout << " \t& ";
    printf("%4.2f", b);
    std::cout << " \t& ";
    if (b > 0)
        printf("%4.2f", s/sqrt(b));
    else 
        printf("N/A");
    std::cout << " \\\\ \\hline" << std::endl;

    if (text == "ALL") {
        std::cout << "gr->SetPoint(P, " << b << ", " << s << ");" << std::endl;
        std::cout << "gr->SetPointError(P, " << errb << ", " << errs << ");" << std::endl;
    }

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

void makeQuickTables(TString id, TString str_njets) {

    int njets = 0;
    if (str_njets == "all") njets = 0; // 0 jets and up
    if (str_njets == "2j") njets = 2; // 2 jets and up

    gROOT->ProcessLine(".L ../tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    TFile f("../histos_mc_" + id + ".root", "READ");
    gROOT->cd();    

    //
    // just sani loose
    //
    TH1F *ttbar_hyp_njets_ee = (TH1F*)f.Get("ttbar_hyp_njets_ee");
    formatHist(ttbar_hyp_njets_ee, true);
    ttbar_hyp_njets_ee->GetYaxis()->SetRangeUser(0, 5);
    TH1F *wjets_hyp_njets_ee = (TH1F*)f.Get("wjets_hyp_njets_ee");
    formatHist(wjets_hyp_njets_ee, false);

    TH1F *ttbar_hyp_njets_em = (TH1F*)f.Get("ttbar_hyp_njets_em");
    formatHist(ttbar_hyp_njets_em, true);
    ttbar_hyp_njets_em->GetYaxis()->SetRangeUser(0, 15);
    TH1F *wjets_hyp_njets_em = (TH1F*)f.Get("wjets_hyp_njets_em");
    formatHist(wjets_hyp_njets_em, false);

    TH1F *ttbar_hyp_njets_mm = (TH1F*)f.Get("ttbar_hyp_njets_mm");
    formatHist(ttbar_hyp_njets_mm, true);
    ttbar_hyp_njets_mm->GetYaxis()->SetRangeUser(0, 10);
    TH1F *wjets_hyp_njets_mm = (TH1F*)f.Get("wjets_hyp_njets_mm");
    formatHist(wjets_hyp_njets_mm, false);

    TH1F *ttbar_hyp_njets_all = (TH1F*)f.Get("ttbar_hyp_njets_all");
    formatHist(ttbar_hyp_njets_all, true);
    ttbar_hyp_njets_all->GetYaxis()->SetRangeUser(0, 25);
    TH1F *wjets_hyp_njets_all = (TH1F*)f.Get("wjets_hyp_njets_all");
    formatHist(wjets_hyp_njets_all, false);

    //
    // make a simple plot and table
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
    c1->SaveAs("../results/study_njets_ee_" + id + ".png");

    ttbar_hyp_njets_em->Draw("HIST");
    wjets_hyp_njets_em->Draw("HIST SAME E1");
    l1->Draw();
    c1->SaveAs("../results/study_njets_em_" + id + ".png");

    ttbar_hyp_njets_all->Draw("HIST");
    wjets_hyp_njets_all->Draw("HIST SAME E1");
    l1->Draw();
    c1->SaveAs("../results/study_njets_all_" + id + ".png");

    // print table
//    std::cout << "\\begin{table}[ht]" << std::endl;
//    std::cout << "\\caption{Yields for two or more jets (" + id + ")}" << std::endl;
//    std::cout << "\\begin{center}" << std::endl;
//    std::cout << "\\begin{tabular}{|l*{4}{|c}|r|}\\hline" << std::endl;
//    std::cout << "Final state   & S & B & S/$\\sqrt{B}$  \\\\ \\hline" << std::endl;
//    printIntegral("EE", ttbar_hyp_njets_ee, wjets_hyp_njets_ee, njets);
//    printIntegral("EM", ttbar_hyp_njets_em, wjets_hyp_njets_em, njets);
//    printIntegral("MM", ttbar_hyp_njets_mm, wjets_hyp_njets_mm, njets);
//    printIntegral("ALL", ttbar_hyp_njets_all, wjets_hyp_njets_all, njets);
//    std::cout <<"\\end{tabular}" << std::endl;
//    std::cout <<"\\end{center}" << std::endl;
//    std::cout << "\\end{table}" << std::endl;

    printIntegral("EE+EM", ttbar_hyp_njets_ee, wjets_hyp_njets_ee, ttbar_hyp_njets_em, wjets_hyp_njets_em, njets);


}

