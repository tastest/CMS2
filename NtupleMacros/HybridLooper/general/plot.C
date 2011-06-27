void plot(TString lm, TString t1, TString t2, TString var) 
{

    TFile f("histos_mc.root", "read");
    gROOT->cd();

    Int_t rebin = 20;

    TH1F *h1_lm0_t1 = (TH1F*)f.Get(lm + "_" + var + "_" + t1)->Clone();
    TH1F *h1_lm0_t2 = (TH1F*)f.Get(lm + "_" + var + "_" + t2)->Clone();
    h1_lm0_t1->SetDirectory(gDirectory);
    h1_lm0_t2->SetDirectory(gDirectory);
    h1_lm0_t1->Sumw2();
    h1_lm0_t2->Sumw2();
    h1_lm0_t1->Rebin(rebin);
    h1_lm0_t2->Rebin(rebin);

    TH1F *h1_dyee_t1 = (TH1F*)f.Get("dyee_" + var + "_" + t1)->Clone();
    TH1F *h1_dyee_t2 = (TH1F*)f.Get("dyee_" + var + "_" + t2)->Clone();
    h1_dyee_t1->SetDirectory(gDirectory);
    h1_dyee_t2->SetDirectory(gDirectory);
    h1_dyee_t1->Sumw2();
    h1_dyee_t2->Sumw2();
    h1_dyee_t1->Rebin(rebin);
    h1_dyee_t2->Rebin(rebin);
    h1_dyee_t1->SetFillColor(kYellow);
    h1_dyee_t2->SetFillColor(kYellow);

    TH1F *h1_dyee = (TH1F*) h1_dyee_t1->Clone();
    h1_dyee->Divide(h1_dyee_t2);
    h1_dyee->SetFillColor(kYellow);
    h1_dyee->SetLineWidth(2);

    TH1F *h1_lm0 = (TH1F*) h1_lm0_t1->Clone();
    h1_lm0->Divide(h1_lm0_t2);
    h1_lm0->SetMarkerStyle(20);
    h1_lm0->SetLineStyle(kDashed);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    h1_dyee->Draw("HIST E1");
    h1_lm0->Draw("SAME E1");

    h1_dyee->GetYaxis()->SetRangeUser(0, 10);
    h1_dyee->SetTitle("ratio of T1 events to T2 events as fn of " + var);
    h1_dyee->GetXaxis()->SetTitle(var);

    TLegend *l1 = new TLegend(0.5, 0.7, 0.9, 0.9);
    l1->AddEntry(h1_dyee, "DYEE", "fp");
    l1->AddEntry(h1_dyee, "LM0", "lp");

    l1->Draw();


    TCanvas *c2 = new TCanvas();
    c2->Divide(2, 1);
    c2->cd(1);
    h1_lm0_t1->Draw("HIST");
    h1_lm0_t2->Draw("SAME E1");
    c2->cd(2);
    h1_dyee_t1->Draw("HIST");
    h1_dyee_t2->Draw("SAME E1");




}


