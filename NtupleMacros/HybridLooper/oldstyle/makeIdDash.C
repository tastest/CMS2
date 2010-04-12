
void format(TH1F *hist, bool isSignal)
{

        hist->SetDirectory(gDirectory);

        if (isSignal) {
                hist->SetMarkerStyle(22);
                hist->SetMarkerColor(kRed);
                hist->SetLineWidth(2);
                hist->SetLineColor(kRed);
        } else {
                hist->SetFillColor(kYellow);
        }

}

void format_mc(TH1F *hist, bool isSignal, float luminorm)
{

        hist->SetDirectory(gDirectory);
        hist->Scale(luminorm);
        hist->SetLineWidth(2.0);

    if (isSignal) {
        hist->SetLineColor(kBlue);
    } else {
        hist->SetLineColor(kGreen);
        
    }   
}

void makeIdDash(TString det)
{

        gROOT->ProcessLine(".L tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");
        TFile f("histos_data.root", "READ");
        TFile f_ref("histos_reference.root", "READ");
        float lumi = 0.201; // lumi integrated (nb)
        float luminorm = 25*lumi/1e+06; // reference being scaled up 10 times

        gROOT->cd();

        gStyle->SetOptTitle(1);
        TCanvas *c1 = new TCanvas();
        c1->Divide(2,3);

        c1->cd(1);
        TH1F *antiselected_fbrem = (TH1F*)f.Get("whunt_ele_antiselected_fbrem_" + det);
        TH1F *selected_fbrem = (TH1F*)f.Get("whunt_ele_selected_fbrem_" + det);
        TH1F *selected_ref_fbrem = (TH1F*)f_ref.Get("wjets_ele_selected_fbrem_" + det);
        TH1F *antiselected_ref_fbrem = (TH1F*)f_ref.Get("QCDpt30_ele_antiselected_fbrem_" + det);
        format_mc(selected_ref_fbrem, true, luminorm);
        format_mc(antiselected_ref_fbrem, false, luminorm);
        format(antiselected_fbrem, false);
        format(selected_fbrem, true);
        antiselected_fbrem->Draw();
        selected_fbrem->Draw("SAME E1");
        selected_ref_fbrem->Draw("SAME HIST");
        antiselected_ref_fbrem->Draw("SAME HIST");

        TLegend *l1 = new TLegend(0.7, 0.7, 0.9, 0.9);
        l1->SetFillColor(kWhite);
        l1->SetShadowColor(kWhite);
        l1->SetLineColor(kWhite);
        l1->AddEntry(antiselected_fbrem, "MET < 15", "f");
        l1->AddEntry(selected_fbrem, "MET > 20", "lp");
        l1->Draw();

        c1->cd(2);
        TH1F *antiselected_eOverPIn = (TH1F*)f.Get("whunt_ele_antiselected_eOverPIn_" + det);
        TH1F *selected_eOverPIn = (TH1F*)f.Get("whunt_ele_selected_eOverPIn_" + det);
        TH1F *selected_ref_eOverPIn = (TH1F*)f_ref.Get("wjets_ele_selected_eOverPIn_" + det);
        TH1F *antiselected_ref_eOverPIn = (TH1F*)f_ref.Get("QCDpt30_ele_antiselected_eOverPIn_" + det);
        format_mc(antiselected_ref_eOverPIn, false, luminorm);
        format_mc(selected_ref_eOverPIn, true, luminorm);
        format(antiselected_eOverPIn, false);
        format(selected_eOverPIn, true);
        antiselected_eOverPIn->Draw();
        selected_eOverPIn->Draw("SAME E1");
        selected_ref_eOverPIn->Draw("SAME");
        antiselected_ref_eOverPIn->Draw("SAME HIST");
        l1->Draw();

        c1->cd(3);
        TH1F *antiselected_sigmaIEtaIEta = (TH1F*)f.Get("whunt_ele_antiselected_sigmaIEtaIEta_" + det);
        TH1F *selected_sigmaIEtaIEta = (TH1F*)f.Get("whunt_ele_selected_sigmaIEtaIEta_" + det);
        TH1F *selected_ref_sigmaIEtaIEta = (TH1F*)f_ref.Get("wjets_ele_selected_sigmaIEtaIEta_" + det);
        TH1F *antiselected_ref_sigmaIEtaIEta = (TH1F*)f_ref.Get("QCDpt30_ele_antiselected_sigmaIEtaIEta_" + det);
        format_mc(antiselected_ref_sigmaIEtaIEta, false, luminorm);
        format_mc(selected_ref_sigmaIEtaIEta, true, luminorm);
        format(antiselected_sigmaIEtaIEta, false);
        format(selected_sigmaIEtaIEta, true);
        antiselected_sigmaIEtaIEta->Draw();
        selected_sigmaIEtaIEta->Draw("SAME E1");
        selected_ref_sigmaIEtaIEta->Draw("SAME");
        antiselected_ref_sigmaIEtaIEta->Draw("SAME HIST");
        l1->Draw();

        c1->cd(4);
        TH1F *antiselected_dEtaIn = (TH1F*)f.Get("whunt_ele_antiselected_dEtaIn_" + det);
        TH1F *selected_dEtaIn = (TH1F*)f.Get("whunt_ele_selected_dEtaIn_" + det);
        TH1F *selected_ref_dEtaIn = (TH1F*)f_ref.Get("wjets_ele_selected_dEtaIn_" + det);
        TH1F *antiselected_ref_dEtaIn = (TH1F*)f_ref.Get("QCDpt30_ele_antiselected_dEtaIn_" + det);
        format_mc(antiselected_ref_dEtaIn, false, luminorm);
        format_mc(selected_ref_dEtaIn, true, luminorm);
        format(antiselected_dEtaIn, false);
        format(selected_dEtaIn, true);
        antiselected_dEtaIn->Draw();
        selected_dEtaIn->Draw("SAME E1");
        selected_ref_dEtaIn->Draw("SAME");
        antiselected_ref_dEtaIn->Draw("SAME HIST");
        l1->Draw();

        c1->cd(5);
        TH1F *antiselected_dPhiIn = (TH1F*)f.Get("whunt_ele_antiselected_dPhiIn_" + det);
        TH1F *selected_dPhiIn = (TH1F*)f.Get("whunt_ele_selected_dPhiIn_" + det);
        TH1F *selected_ref_dPhiIn = (TH1F*)f_ref.Get("wjets_ele_selected_dPhiIn_" + det);
        TH1F *antiselected_ref_dPhiIn = (TH1F*)f_ref.Get("QCDpt30_ele_antiselected_dPhiIn_" + det);
        format_mc(antiselected_ref_dPhiIn, false, luminorm);
        format_mc(selected_ref_dPhiIn, true, luminorm);
        format(antiselected_dPhiIn, false);
        format(selected_dPhiIn, true);
        antiselected_dPhiIn->Draw();
        selected_dPhiIn->Draw("SAME E1");
        selected_ref_dPhiIn->Draw("SAME");
        antiselected_ref_dPhiIn->Draw("SAME HIST");
        l1->Draw();

        c1->cd(6);
        TH1F *antiselected_hOverE = (TH1F*)f.Get("whunt_ele_antiselected_hOverE_" + det);
        TH1F *selected_hOverE = (TH1F*)f.Get("whunt_ele_selected_hOverE_" + det);
        TH1F *selected_ref_hOverE = (TH1F*)f_ref.Get("wjets_ele_selected_hOverE_" + det);
        TH1F *antiselected_ref_hOverE = (TH1F*)f_ref.Get("QCDpt30_ele_antiselected_hOverE_" + det);
        format_mc(antiselected_ref_hOverE, false, luminorm);
        format_mc(selected_ref_hOverE, true, luminorm);
        format(antiselected_hOverE, false);
        format(selected_hOverE, true);
        antiselected_hOverE->Draw();
        selected_hOverE->Draw("SAME E1");
        selected_ref_hOverE->Draw("SAME");
        antiselected_ref_hOverE->Draw("SAME HIST");
        l1->Draw();

}

