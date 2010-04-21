

void printIntegral(std::string text, const TH1F *histSignal, const TH1F *histBackground)
{

    float s = histSignal->Integral(0, histSignal->GetNbinsX() + 1);
    float b = histBackground->Integral(0, histSignal->GetNbinsX() + 1);

    std::cout << text << " \t& ";
    printf("%4.2f", s);
    std::cout << " \t& ";
    printf("%4.2f", b);
    std::cout << " \t& ";
    printf("%4.2f", s/sqrt(b));
    std::cout << " \\\\ \\hline" << std::endl;

}

void formatHist(TH1F *hist, bool signal, unsigned int idType) {
    hist->SetDirectory(gDirectory);
    hist->Rebin(10);
    hist->SetLineWidth(2);

    if (idType == 0) {
        if (signal) {
            hist->SetLineColor(kBlack);
            hist->SetMarkerStyle(20);
            hist->SetMarkerColor(kBlack);
        }
        else {
            hist->SetLineColor(kRed);
            hist->SetMarkerStyle(20);
            hist->SetMarkerColor(kRed);
        }
    }
    if (idType == 1) {
        if (signal) {
            hist->SetLineColor(kBlack);
            hist->SetMarkerStyle(4);
            hist->SetMarkerColor(kBlack);
        }
        else {
            hist->SetLineColor(kRed);
            hist->SetMarkerStyle(4);
            hist->SetMarkerColor(kRed);
        }
    }
    if (idType == 2) {
        if (signal) {
            hist->SetLineColor(kGray +2);
            hist->SetFillColor(kGray +2);
            hist->SetFillStyle(3004);
        }
        else {
            hist->SetLineColor(kRed - 9);
            hist->SetFillColor(kRed - 9);
            hist->SetFillStyle(3005);
        }
    }
    if (idType == 3) {
        if (signal) {
            hist->SetLineColor(kBlue);
            hist->SetMarkerStyle(20);
            hist->SetMarkerColor(kBlue);
        }
        else {
            hist->SetLineColor(kMagenta);
            hist->SetMarkerStyle(20);
            hist->SetMarkerColor(kMagenta);
        }
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

void printSanityStudy(TString det) {


    gROOT->ProcessLine(".L tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    TFile f("histos_eleid_pt20up.root", "READ");
    gROOT->cd();    


    //
    // just sani loose
    //
    TH1F *s_h1_hyp_idstudy_after_classExpLoose_pt = (TH1F*)f.Get("ttbar_h1_hyp_idstudy_after_classExpLoose_pt_" + det + "_ee");
    formatHist(s_h1_hyp_idstudy_after_classExpLoose_pt, true, 0);
    TH1F *b_h1_hyp_idstudy_after_classExpLoose_pt = (TH1F*)f.Get("wm_h1_hyp_idstudy_after_classExpLoose_pt_" + det + "_ee");
    formatHist(b_h1_hyp_idstudy_after_classExpLoose_pt, false, 0);

    TH1F *s_h1_hyp_idstudy_after_classExpLoose_reliso = (TH1F*)f.Get("ttbar_h1_hyp_idstudy_after_classExpLoose_reliso_" + det + "_ee");
    formatHist(s_h1_hyp_idstudy_after_classExpLoose_reliso, true, 0);
    TH1F *b_h1_hyp_idstudy_after_classExpLoose_reliso = (TH1F*)f.Get("wm_h1_hyp_idstudy_after_classExpLoose_reliso_" + det + "_ee");
    formatHist(b_h1_hyp_idstudy_after_classExpLoose_reliso, false, 0);


    //
    // sani loose with sani iso
    //
    TH1F *s_h1_hyp_idstudy_after_classExpLooseFull_pt = (TH1F*)f.Get("ttbar_h1_hyp_idstudy_after_classExpLooseFull_pt_" + det + "_ee");
    formatHist(s_h1_hyp_idstudy_after_classExpLooseFull_pt, true, 3);
    TH1F *b_h1_hyp_idstudy_after_classExpLooseFull_pt = (TH1F*)f.Get("wm_h1_hyp_idstudy_after_classExpLooseFull_pt_" + det + "_ee");
    formatHist(b_h1_hyp_idstudy_after_classExpLooseFull_pt, false, 3);

    //
    // sani loose with reliso01
    //
    TH1F *s_h1_hyp_idstudy_after_classExpLooseRel01_pt = (TH1F*)f.Get("ttbar_h1_hyp_idstudy_after_classExpLooseRel01_pt_" + det + "_ee");
    formatHist(s_h1_hyp_idstudy_after_classExpLooseRel01_pt, true, 0);
    TH1F *b_h1_hyp_idstudy_after_classExpLooseRel01_pt = (TH1F*)f.Get("wm_h1_hyp_idstudy_after_classExpLooseRel01_pt_" + det + "_ee");
    formatHist(b_h1_hyp_idstudy_after_classExpLooseRel01_pt, false, 0);

    //
    // just sani tight
    //
    TH1F *s_h1_hyp_idstudy_after_classExpTight_pt = (TH1F*)f.Get("ttbar_h1_hyp_idstudy_after_classExpTight_pt_" + det + "_ee");
    formatHist(s_h1_hyp_idstudy_after_classExpTight_pt, true, 1);
    TH1F *b_h1_hyp_idstudy_after_classExpTight_pt = (TH1F*)f.Get("wm_h1_hyp_idstudy_after_classExpTight_pt_" + det + "_ee");
    formatHist(b_h1_hyp_idstudy_after_classExpTight_pt, false, 1);

    TH1F *s_h1_hyp_idstudy_after_classExpTight_reliso = (TH1F*)f.Get("ttbar_h1_hyp_idstudy_after_classExpTight_reliso_" + det + "_ee");
    formatHist(s_h1_hyp_idstudy_after_classExpTight_reliso, true, 1);
    TH1F *b_h1_hyp_idstudy_after_classExpTight_reliso = (TH1F*)f.Get("wm_h1_hyp_idstudy_after_classExpTight_reliso_" + det + "_ee");
    formatHist(b_h1_hyp_idstudy_after_classExpTight_reliso, false, 1);


    //
    // sani tight with sani iso
    //
    TH1F *s_h1_hyp_idstudy_after_classExpTightFull_pt = (TH1F*)f.Get("ttbar_h1_hyp_idstudy_after_classExpTightFull_pt_" + det + "_ee");
    formatHist(s_h1_hyp_idstudy_after_classExpTightFull_pt, true, 3);
    TH1F *b_h1_hyp_idstudy_after_classExpTightFull_pt = (TH1F*)f.Get("wm_h1_hyp_idstudy_after_classExpTightFull_pt_" + det + "_ee");
    formatHist(b_h1_hyp_idstudy_after_classExpTightFull_pt, false, 3);

    //
    // sani tight with reliso01
    //
    TH1F *s_h1_hyp_idstudy_after_classExpTightRel01_pt = (TH1F*)f.Get("ttbar_h1_hyp_idstudy_after_classExpTightRel01_pt_" + det + "_ee");
    formatHist(s_h1_hyp_idstudy_after_classExpTightRel01_pt, true, 1);
    TH1F *b_h1_hyp_idstudy_after_classExpTightRel01_pt = (TH1F*)f.Get("wm_h1_hyp_idstudy_after_classExpTightRel01_pt_" + det + "_ee");
    formatHist(b_h1_hyp_idstudy_after_classExpTightRel01_pt, false, 1);


    //
    // just cand01
    //
    TH1F *s_h1_hyp_idstudy_after_cand01_pt = (TH1F*)f.Get("ttbar_h1_hyp_idstudy_after_cand01_pt_" + det + "_ee");
    formatHist(s_h1_hyp_idstudy_after_cand01_pt, true, 2);
    TH1F *b_h1_hyp_idstudy_after_cand01_pt = (TH1F*)f.Get("wm_h1_hyp_idstudy_after_cand01_pt_" + det + "_ee");
    formatHist(b_h1_hyp_idstudy_after_cand01_pt, false, 2);

    TH1F *s_h1_hyp_idstudy_after_cand01_reliso = (TH1F*)f.Get("ttbar_h1_hyp_idstudy_after_cand01_reliso_" + det + "_ee");
    formatHist(s_h1_hyp_idstudy_after_cand01_reliso, true, 2);
    TH1F *b_h1_hyp_idstudy_after_cand01_reliso = (TH1F*)f.Get("wm_h1_hyp_idstudy_after_cand01_reliso_" + det + "_ee");
    formatHist(b_h1_hyp_idstudy_after_cand01_reliso, false, 2);


    //
    // cand01 with reliso01
    //
    TH1F *s_h1_hyp_idstudy_after_cand01Rel01_pt = (TH1F*)f.Get("ttbar_h1_hyp_idstudy_after_cand01Rel01_pt_" + det + "_ee"); 
    formatHist(s_h1_hyp_idstudy_after_cand01Rel01_pt, true, 2); 
    TH1F *b_h1_hyp_idstudy_after_cand01Rel01_pt = (TH1F*)f.Get("wm_h1_hyp_idstudy_after_cand01Rel01_pt_" + det + "_ee"); 
    formatHist(b_h1_hyp_idstudy_after_cand01Rel01_pt, false, 2);

    //
    // Efficiencies
    //

    // efficiency of sani loose reliso01
    TH1F *s_eff_classExpLooseRel01_pt = (TH1F*)s_h1_hyp_idstudy_after_classExpLooseRel01_pt->Clone("s_eff_classExpLooseRel01_pt");
    s_eff_classExpLooseRel01_pt->Divide(s_h1_hyp_idstudy_after_classExpLoose_pt);
    setError(s_h1_hyp_idstudy_after_classExpLooseRel01_pt, s_h1_hyp_idstudy_after_classExpLoose_pt, s_eff_classExpLooseRel01_pt);
    TH1F *b_eff_classExpLooseRel01_pt = (TH1F*)b_h1_hyp_idstudy_after_classExpLooseRel01_pt->Clone("b_eff_classExpLooseRel01_pt");
    b_eff_classExpLooseRel01_pt->Divide(b_h1_hyp_idstudy_after_classExpLoose_pt);
    setError(b_h1_hyp_idstudy_after_classExpLooseRel01_pt, b_h1_hyp_idstudy_after_classExpLoose_pt, b_eff_classExpLooseRel01_pt);

    // efficiency of sani loose full
    TH1F *s_eff_classExpLooseFull_pt = (TH1F*)s_h1_hyp_idstudy_after_classExpLooseFull_pt->Clone("s_eff_classExpLooseFull_pt");
    s_eff_classExpLooseFull_pt->Divide(s_h1_hyp_idstudy_after_classExpLoose_pt);
    setError(s_h1_hyp_idstudy_after_classExpLooseFull_pt, s_h1_hyp_idstudy_after_classExpLoose_pt, s_eff_classExpLooseFull_pt);
    TH1F *b_eff_classExpLooseFull_pt = (TH1F*)b_h1_hyp_idstudy_after_classExpLooseFull_pt->Clone("b_eff_classExpLooseFull_pt");
    b_eff_classExpLooseFull_pt->Divide(b_h1_hyp_idstudy_after_classExpLoose_pt);
    setError(b_h1_hyp_idstudy_after_classExpLooseFull_pt, s_h1_hyp_idstudy_after_classExpLoose_pt, b_eff_classExpLooseFull_pt);

    // efficiency of cand0 full
    TH1F *s_eff_cand01Rel01_pt = (TH1F*)s_h1_hyp_idstudy_after_cand01Rel01_pt->Clone("s_eff_cand01Rel01_pt");
    s_eff_cand01Rel01_pt->Divide(s_h1_hyp_idstudy_after_cand01_pt);
    setError(s_h1_hyp_idstudy_after_cand01Rel01_pt, s_h1_hyp_idstudy_after_cand01_pt, s_eff_cand01Rel01_pt);
    TH1F *b_eff_cand01Rel01_pt = (TH1F*)b_h1_hyp_idstudy_after_cand01Rel01_pt->Clone("b_eff_cand01Rel01_pt");
    b_eff_cand01Rel01_pt->Divide(b_h1_hyp_idstudy_after_cand01_pt);
    setError(b_h1_hyp_idstudy_after_cand01Rel01_pt, b_h1_hyp_idstudy_after_cand01_pt, b_eff_cand01Rel01_pt);

    //
    //

    // efficiency of sani tight reliso01
    TH1F *s_eff_classExpTightRel01_pt = (TH1F*)s_h1_hyp_idstudy_after_classExpTightRel01_pt->Clone("s_eff_classExpTightRel01_pt");
    s_eff_classExpTightRel01_pt->Divide(s_h1_hyp_idstudy_after_classExpTight_pt);
    setError(s_h1_hyp_idstudy_after_classExpTightRel01_pt, s_h1_hyp_idstudy_after_classExpTight_pt, s_eff_classExpTightRel01_pt);
    TH1F *b_eff_classExpTightRel01_pt = (TH1F*)b_h1_hyp_idstudy_after_classExpTightRel01_pt->Clone("b_eff_classExpTightRel01_pt");
    b_eff_classExpTightRel01_pt->Divide(b_h1_hyp_idstudy_after_classExpTight_pt);
    setError(b_h1_hyp_idstudy_after_classExpTightRel01_pt, b_h1_hyp_idstudy_after_classExpTight_pt, b_eff_classExpTightRel01_pt);
    
    // efficiency of sani tight full
    TH1F *s_eff_classExpTightFull_pt = (TH1F*)s_h1_hyp_idstudy_after_classExpTightFull_pt->Clone("s_eff_classExpTightFull_pt");
    s_eff_classExpTightFull_pt->Divide(s_h1_hyp_idstudy_after_classExpTight_pt);
    setError(s_h1_hyp_idstudy_after_classExpTightFull_pt, s_h1_hyp_idstudy_after_classExpTight_pt, s_eff_classExpTightFull_pt);
    TH1F *b_eff_classExpTightFull_pt = (TH1F*)b_h1_hyp_idstudy_after_classExpTightFull_pt->Clone("b_eff_classExpTightFull_pt");
    b_eff_classExpTightFull_pt->Divide(b_h1_hyp_idstudy_after_classExpTight_pt);
    setError(b_h1_hyp_idstudy_after_classExpTightFull_pt, s_h1_hyp_idstudy_after_classExpTight_pt, b_eff_classExpTightFull_pt);

    // efficiency of cand0 full
    TH1F *s_eff_cand01Rel01_pt = (TH1F*)s_h1_hyp_idstudy_after_cand01Rel01_pt->Clone("s_eff_cand01Rel01_pt");
    s_eff_cand01Rel01_pt->Divide(s_h1_hyp_idstudy_after_cand01_pt);
    setError(s_h1_hyp_idstudy_after_cand01Rel01_pt, s_h1_hyp_idstudy_after_cand01_pt, s_eff_cand01Rel01_pt);
    TH1F *b_eff_cand01Rel01_pt = (TH1F*)b_h1_hyp_idstudy_after_cand01Rel01_pt->Clone("b_eff_cand01Rel01_pt");
    b_eff_cand01Rel01_pt->Divide(b_h1_hyp_idstudy_after_cand01_pt);
    setError(b_h1_hyp_idstudy_after_cand01Rel01_pt, b_h1_hyp_idstudy_after_cand01_pt, b_eff_cand01Rel01_pt);



    //
    // pt distributions before iso loose
    //

    TCanvas *c1 = new TCanvas();
    c1->cd();

    TLegend *l1 = new TLegend(0.5, 0.5, 0.9, 0.9);
    l1->SetFillColor(kWhite);
    l1->SetLineColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(s_h1_hyp_idstudy_after_cand01_pt, "S (cand01)", "lf");
    l1->AddEntry(b_h1_hyp_idstudy_after_cand01_pt, "BG (cand01)", "lf");
    l1->AddEntry(s_h1_hyp_idstudy_after_classExpLoose_pt, "S (class loose)", "lp");
    l1->AddEntry(b_h1_hyp_idstudy_after_classExpLoose_pt, "BG (class loose)", "lp");


    s_h1_hyp_idstudy_after_cand01_pt->Draw("HIST");
    b_h1_hyp_idstudy_after_cand01_pt->Draw("HIST SAME");
    s_h1_hyp_idstudy_after_classExpLoose_pt->Draw("HIST SAME E1");
    b_h1_hyp_idstudy_after_classExpLoose_pt->Draw("HIST SAME E1");
    s_h1_hyp_idstudy_after_cand01_pt->GetYaxis()->SetRangeUser(0, 10);
    s_h1_hyp_idstudy_after_cand01_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
    l1->Draw();

    c1->SaveAs("results/study_" + det + "_classExpLoose.png");


    //
    // pt distributions before iso tight
    //
    TCanvas *c2 = new TCanvas();
    c2->cd();

    TLegend *l2 = new TLegend(0.5, 0.5, 0.9, 0.9);
    l2->SetFillColor(kWhite);
    l2->SetLineColor(kWhite);
    l2->SetShadowColor(kWhite);
    l2->AddEntry(s_h1_hyp_idstudy_after_cand01_pt, "S (cand01)", "lf");
    l2->AddEntry(b_h1_hyp_idstudy_after_cand01_pt, "BG (cand01)", "lf");
    l2->AddEntry(s_h1_hyp_idstudy_after_classExpTight_pt, "S (class tight)", "lp");
    l2->AddEntry(b_h1_hyp_idstudy_after_classExpTight_pt, "BG (class tight)", "lp");

    
    s_h1_hyp_idstudy_after_cand01_pt->Draw("HIST");
    b_h1_hyp_idstudy_after_cand01_pt->Draw("HIST SAME");
    s_h1_hyp_idstudy_after_classExpTight_pt->Draw("HIST SAME E1");
    b_h1_hyp_idstudy_after_classExpTight_pt->Draw("HIST SAME E1");
    s_h1_hyp_idstudy_after_cand01_pt->GetYaxis()->SetRangeUser(0, 6);
    s_h1_hyp_idstudy_after_cand01_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
    l2->Draw();
    
    c2->SaveAs("results/study_" + det + "_classExpTight.png");

    // print table of IDs before Iso
    std::cout << "\\begin{table}[ht]" << std::endl;
    std::cout << "\\caption{Electron ID before isolation " + det + "}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|l*{4}{|c}|r|}\\hline" << std::endl;
    std::cout << "ID type   & S & B & S/$\\sqrt{B}$  \\\\ \\hline" << std::endl;
    printIntegral("cand01", s_h1_hyp_idstudy_after_cand01_pt, b_h1_hyp_idstudy_after_cand01_pt);
    printIntegral("class loose", s_h1_hyp_idstudy_after_classExpLoose_pt, b_h1_hyp_idstudy_after_classExpLoose_pt);
    printIntegral("class tight", s_h1_hyp_idstudy_after_classExpTight_pt, b_h1_hyp_idstudy_after_classExpTight_pt);
    std::cout <<"\\end{tabular}" << std::endl;
    std::cout <<"\\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;

    //
    // pt distributions after iso loose
    //
    TCanvas *c3 = new TCanvas();
    c3->cd();

    TLegend *l3 = new TLegend(0.5, 0.5, 0.9, 0.9);
    l3->SetFillColor(kWhite);
    l3->SetLineColor(kWhite);
    l3->SetShadowColor(kWhite);
    l3->AddEntry(s_h1_hyp_idstudy_after_cand01Rel01_pt, "S (cand01 rel01)", "lf");
    l3->AddEntry(b_h1_hyp_idstudy_after_cand01Rel01_pt, "BG (cand01 rel01)", "lf");
    l3->AddEntry(s_h1_hyp_idstudy_after_classExpLooseRel01_pt, "S (class loose rel01)", "lp");
    l3->AddEntry(b_h1_hyp_idstudy_after_classExpLooseRel01_pt, "BG (class loose rel01)", "lp");
    l3->AddEntry(s_h1_hyp_idstudy_after_classExpLooseFull_pt, "S (class loose full)", "lp");
    l3->AddEntry(b_h1_hyp_idstudy_after_classExpLooseFull_pt, "BG (class loose full)", "lp");
    

    s_h1_hyp_idstudy_after_cand01Rel01_pt->Draw("HIST");
    b_h1_hyp_idstudy_after_cand01Rel01_pt->Draw("HIST SAME");
    s_h1_hyp_idstudy_after_classExpLooseRel01_pt->Draw("HIST SAME E1");
    b_h1_hyp_idstudy_after_classExpLooseRel01_pt->Draw("HIST SAME E1");
    s_h1_hyp_idstudy_after_classExpLooseFull_pt->Draw("HIST SAME E1");
    b_h1_hyp_idstudy_after_classExpLooseFull_pt->Draw("HIST SAME E1");
    s_h1_hyp_idstudy_after_cand01Rel01_pt->GetYaxis()->SetRangeUser(0, 6.5);
    s_h1_hyp_idstudy_after_cand01Rel01_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
    l3->Draw();

    c3->SaveAs("results/study_" + det + "_classExpLooseRel01.png");
    
    
    //
    // pt distributions after iso tight
    //
    TCanvas *c4 = new TCanvas();
    c4->cd();
    
    TLegend *l4 = new TLegend(0.5, 0.5, 0.9, 0.9);
    l4->SetFillColor(kWhite);
    l4->SetLineColor(kWhite);
    l4->SetShadowColor(kWhite);
    l4->AddEntry(s_h1_hyp_idstudy_after_cand01Rel01_pt, "S (cand01 rel01)", "lf");
    l4->AddEntry(b_h1_hyp_idstudy_after_cand01Rel01_pt, "BG (cand01 rel01)", "lf");
    l4->AddEntry(s_h1_hyp_idstudy_after_classExpTightRel01_pt, "S (class tight rel01)", "lp");
    l4->AddEntry(b_h1_hyp_idstudy_after_classExpTightRel01_pt, "BG (class tight rel01)", "lp");
    l4->AddEntry(s_h1_hyp_idstudy_after_classExpTightFull_pt, "S (class tight full)", "lp");
    l4->AddEntry(b_h1_hyp_idstudy_after_classExpTightFull_pt, "BG (class tight full)", "lp");
    
    s_h1_hyp_idstudy_after_cand01Rel01_pt->Draw("HIST");
    b_h1_hyp_idstudy_after_cand01Rel01_pt->Draw("HIST SAME");
    s_h1_hyp_idstudy_after_classExpTightRel01_pt->Draw("HIST SAME E1");
    b_h1_hyp_idstudy_after_classExpTightRel01_pt->Draw("HIST SAME E1");
    s_h1_hyp_idstudy_after_classExpTightFull_pt->Draw("HIST SAME E1");
    b_h1_hyp_idstudy_after_classExpTightFull_pt->Draw("HIST SAME E1");
    s_h1_hyp_idstudy_after_cand01Rel01_pt->GetYaxis()->SetRangeUser(0, 10);
    s_h1_hyp_idstudy_after_cand01Rel01_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
    l4->Draw();
    
    c4->SaveAs("results/study_" + det + "_classExpTightRel01.png");

    // print table of IDs after Iso
    std::cout << "\\begin{table}[ht]" << std::endl;
    std::cout << "\\caption{Electron ID after isolation " + det + "}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|l*{4}{|c}|r|}\\hline" << std::endl;
    std::cout << "ID type   & S & B & S/$\\sqrt{B}$  \\\\ \\hline" << std::endl;
    printIntegral("cand01 + reliso", s_h1_hyp_idstudy_after_cand01Rel01_pt, b_h1_hyp_idstudy_after_cand01Rel01_pt);
    printIntegral("class loose + reliso", s_h1_hyp_idstudy_after_classExpLooseRel01_pt, b_h1_hyp_idstudy_after_classExpLooseRel01_pt);
    printIntegral("class loose + full", s_h1_hyp_idstudy_after_classExpLooseFull_pt, b_h1_hyp_idstudy_after_classExpLooseFull_pt);
    printIntegral("class tight + reliso", s_h1_hyp_idstudy_after_classExpTightRel01_pt, b_h1_hyp_idstudy_after_classExpTightRel01_pt);
    printIntegral("class tight + full", s_h1_hyp_idstudy_after_classExpTightFull_pt, b_h1_hyp_idstudy_after_classExpTightFull_pt);
    std::cout <<"\\end{tabular}" << std::endl;
    std::cout <<"\\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;


    //
    // Efficiencies for loose signal
    //
    
    TCanvas *c5 = new TCanvas();
    c5->cd();
    
    TLegend *l5 = new TLegend(0.5, 0.2, 0.9, 0.6);
    l5->SetFillColor(kWhite);
    l5->SetLineColor(kWhite);
    l5->SetShadowColor(kWhite);
    l5->AddEntry(s_eff_cand01Rel01_pt, "S (cand01 rel01)", "lf");
    l5->AddEntry(s_eff_classExpLooseRel01_pt, "S (class loose rel01)", "lp");
    l5->AddEntry(s_eff_classExpLooseFull_pt, "S (class loose full)", "lp");
    
    s_eff_cand01Rel01_pt->Draw("HIST");
    s_eff_classExpLooseRel01_pt->Draw("HIST SAME E1");
    s_eff_classExpLooseFull_pt->Draw("HIST SAME E1");

    s_eff_cand01Rel01_pt->GetYaxis()->SetRangeUser(0, 1.1);
    s_eff_cand01Rel01_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
    l5->Draw();

    c5->SaveAs("results/study_" + det + "_classExpLoose_effS.png");

    //
    // Efficiencies for loose signal
    //

    TCanvas *c6 = new TCanvas();
    c6->cd();

    TLegend *l6 = new TLegend(0.5, 0.2, 0.9, 0.6);
    l6->SetFillColor(kWhite);
    l6->SetLineColor(kWhite);
    l6->SetShadowColor(kWhite);
    l6->AddEntry(s_eff_cand01Rel01_pt, "S (cand01 rel01)", "lf");
    l6->AddEntry(s_eff_classExpTightRel01_pt, "S (class tight rel01)", "lp");
    l6->AddEntry(s_eff_classExpTightFull_pt, "S (class tight full)", "lp");
    
    s_eff_cand01Rel01_pt->Draw("HIST");
    s_eff_classExpTightRel01_pt->Draw("HIST SAME E1");
    s_eff_classExpTightFull_pt->Draw("HIST SAME E1");
    
    s_eff_cand01Rel01_pt->GetYaxis()->SetRangeUser(0, 1.1);
    s_eff_cand01Rel01_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
    l6->Draw();
    
    c6->SaveAs("results/study_" + det + "_classExpTight_effS.png");

    //
    // Efficiencies for loose bg
    //

    TCanvas *c7 = new TCanvas();
    c7->cd();

    TLegend *l7 = new TLegend(0.65, 0.2, 0.95, 0.6);
    l7->SetFillColor(kWhite);
    l7->SetLineColor(kWhite);
    l7->SetShadowColor(kWhite);
    l7->AddEntry(b_eff_cand01Rel01_pt, "BG (cand01 rel01)", "lf");
    l7->AddEntry(b_eff_classExpLooseRel01_pt, "BG (class loose rel01)", "lp");
    l7->AddEntry(b_eff_classExpLooseFull_pt, "BG (class loose full)", "lp");

    b_eff_cand01Rel01_pt->Draw("HIST");
    b_eff_classExpLooseRel01_pt->Draw("HIST SAME E1");
    b_eff_classExpLooseFull_pt->Draw("HIST SAME E1");

    s_eff_cand01Rel01_pt->GetYaxis()->SetRangeUser(0, 1.1);
    s_eff_cand01Rel01_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
    l7->Draw();

    c7->SaveAs("results/study_" + det + "_classExpLoose_effB.png");

    //
    // Efficiencies for loose bg
    //

    TCanvas *c8 = new TCanvas();
    c8->cd();

    TLegend *l8 = new TLegend(0.65, 0.2, 0.95, 0.6);
    l8->SetFillColor(kWhite);
    l8->SetLineColor(kWhite);
    l8->SetShadowColor(kWhite);
    l8->AddEntry(b_eff_cand01Rel01_pt, "BG (cand01 rel01)", "lf");
    l8->AddEntry(b_eff_classExpTightRel01_pt, "BG (class tight rel01)", "lp");
    l8->AddEntry(b_eff_classExpTightFull_pt, "BG (class tight full)", "lp");

    b_eff_cand01Rel01_pt->Draw("HIST");
    b_eff_classExpTightRel01_pt->Draw("HIST SAME E1");
    b_eff_classExpTightFull_pt->Draw("HIST SAME E1");

    s_eff_cand01Rel01_pt->GetYaxis()->SetRangeUser(0, 1.1);
    s_eff_cand01Rel01_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
    l8->Draw();

    c8->SaveAs("results/study_" + det + "_classExpTight_effB.png");


    //
    // reliso distributions after iso tight
    //

    TCanvas *c9 = new TCanvas();
    c9->cd();
    c9->SetLogy();

    TLegend *l9 = new TLegend(0.6, 0.6, 1.0, 1.0);
    l9->SetFillColor(kWhite);
    l9->SetLineColor(kWhite);
    l9->SetShadowColor(kWhite);
    l9->AddEntry(s_h1_hyp_idstudy_after_cand01_reliso, "S (cand01 rel01)", "lf");
    l9->AddEntry(b_h1_hyp_idstudy_after_cand01_reliso, "BG (cand01 rel01)", "lf");
    l9->AddEntry(s_h1_hyp_idstudy_after_classExpTight_reliso, "S (class tight rel01)", "lp");
    l9->AddEntry(b_h1_hyp_idstudy_after_classExpTight_reliso, "BG (class tight rel01)", "lp");

    s_h1_hyp_idstudy_after_cand01_reliso->Draw("HIST");
    b_h1_hyp_idstudy_after_cand01_reliso->Draw("HIST SAME");
    s_h1_hyp_idstudy_after_classExpTight_reliso->Draw("HIST SAME E1");
    b_h1_hyp_idstudy_after_classExpTight_reliso->Draw("HIST SAME E1");
    s_h1_hyp_idstudy_after_cand01_reliso->GetYaxis()->SetRangeUser(0.01, 20);
    s_h1_hyp_idstudy_after_cand01_reliso->GetXaxis()->SetTitle("RelIso");
    l9->Draw();

    c9->SaveAs("results/study_" + det + "_classExpTight_reliso.png");

    //
    // reliso distributions after iso loose
    //

    TCanvas *c10 = new TCanvas();
    c10->cd();
    c10->SetLogy();

    TLegend *l10 = new TLegend(0.6, 0.6, 1.0, 1.0);
    l10->SetFillColor(kWhite);
    l10->SetLineColor(kWhite);
    l10->SetShadowColor(kWhite);
    l10->AddEntry(s_h1_hyp_idstudy_after_cand01_reliso, "S (cand01 rel01)", "lf");
    l10->AddEntry(b_h1_hyp_idstudy_after_cand01_reliso, "BG (cand01 rel01)", "lf");
    l10->AddEntry(s_h1_hyp_idstudy_after_classExpLoose_reliso, "S (class loose rel01)", "lp");
    l10->AddEntry(b_h1_hyp_idstudy_after_classExpLoose_reliso, "BG (class loose rel01)", "lp");

    s_h1_hyp_idstudy_after_cand01_reliso->Draw("HIST");
    b_h1_hyp_idstudy_after_cand01_reliso->Draw("HIST SAME");
    s_h1_hyp_idstudy_after_classExpLoose_reliso->Draw("HIST SAME E1");
    b_h1_hyp_idstudy_after_classExpLoose_reliso->Draw("HIST SAME E1");
    s_h1_hyp_idstudy_after_cand01_reliso->GetYaxis()->SetRangeUser(0.01, 20);
    s_h1_hyp_idstudy_after_cand01_reliso->GetXaxis()->SetTitle("RelIso");
    l10->Draw();

    c10->SaveAs("results/study_" + det + "_classExpLoose_reliso.png");


}

