#include "Utilities.h"
#include "processDYEstResults.C"

void printMetStacks(TString region)
{

        gROOT->ProcessLine(".L ~/tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");
        gStyle->SetOptTitle(1);

//        gROOT->ProcessLine(".L processDYEstResults.C");
//        gROOT->ProcessLine(".L Utilities.h");

        TFile f("DYEstResults_ForWW_MET35.root");

	bool setLog = 1;
	gROOT->cd();

        AnaHist anaHist;
        TLegend *lg_tmp = anaHist.getLegend(f);

        THStack *st_tmp = anaHist.getMetStack(f, "0j", region, "mm");
        TCanvas *c0j_mm = new TCanvas();
        c0j_mm->cd();
	c0j_mm->SetLogy(setLog);
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("MET_{#mu#mu} (GeV)");
        st_tmp->SetTitle("nJets = 0, (" + region + ") hyp_type = mm");
        lg_tmp->Draw();
        saveCanvas(c0j_mm, "st_met_0j_" + region + "_mm");

        THStack *st_tmp = anaHist.getMetStack(f, "1j", region, "mm");
        TCanvas *c1j_mm = new TCanvas();
        c1j_mm->cd();
        c1j_mm->SetLogy(setLog);
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("MET_{#mu#mu} (GeV)");
        st_tmp->SetTitle("nJets = 1, (" + region + ") hyp_type = mm");
        lg_tmp->Draw();
        saveCanvas(c1j_mm, "st_met_1j_" + region + "_mm");

        THStack *st_tmp = anaHist.getMetStack(f, "2j", region, "mm");
        TCanvas *c2j_mm = new TCanvas();
        c2j_mm->cd();
        c2j_mm->SetLogy(setLog);
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("MET_{#mu#mu} (GeV)");
        st_tmp->SetTitle("nJets #geq 2, (" + region + ") hyp_type = mm");
        lg_tmp->Draw();
        saveCanvas(c2j_mm, "st_met_2j_" + region + "_mm");

        THStack *st_tmp = anaHist.getMetStack(f, "0j", region, "ee");
        TCanvas *c0j_ee = new TCanvas();
        c0j_ee->cd();
        c0j_ee->SetLogy(setLog);
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("MET_{ee} (GeV)");
        st_tmp->SetTitle("nJets = 0, (" + region + ") hyp_type = ee");
        lg_tmp->Draw();
        saveCanvas(c0j_ee, "st_met_0j_" + region + "_ee");

        THStack *st_tmp = anaHist.getMetStack(f, "1j", region, "ee");
        TCanvas *c1j_ee = new TCanvas();
        c1j_ee->cd();
        c1j_ee->SetLogy(setLog);
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("MET_{ee} (GeV)");
        st_tmp->SetTitle("nJets = 1, (" + region + ") hyp_type = ee");
        lg_tmp->Draw();
        saveCanvas(c1j_ee, "st_met_1j_" + region + "_ee");

        THStack *st_tmp = anaHist.getMetStack(f, "2j", region, "ee");
        TCanvas *c2j_ee = new TCanvas(); 
        c2j_ee->cd();
        c2j_ee->SetLogy(setLog);
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("MET_{ee} (GeV)");
        st_tmp->SetTitle("nJets #geq 2, (" + region + ") hyp_type = ee");
        lg_tmp->Draw();
        saveCanvas(c2j_ee, "st_met_2j_" + region + "_ee");

        // now for em

        THStack *st_tmp = anaHist.getMetStack(f, "0j", region, "em");
        TCanvas *c0j_em = new TCanvas();
        c0j_em->cd();
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("MET_{e#mu} (GeV)");
        st_tmp->SetTitle("nJets = 0, (" + region + ") hyp_type = em");
        lg_tmp->Draw();
        saveCanvas(c0j_em, "st_met_0j_" + region + "_em");

        THStack *st_tmp = anaHist.getMetStack(f, "1j", region, "em");
        TCanvas *c1j_em = new TCanvas();
        c1j_em->cd();
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("MET_{e#mu} (GeV)");
        st_tmp->SetTitle("nJets = 1, (" + region + ") hyp_type = em");
        lg_tmp->Draw();
        saveCanvas(c1j_em, "st_met_1j_" + region + "_em");

        THStack *st_tmp = anaHist.getMetStack(f, "2j", region, "em");
        TCanvas *c2j_em = new TCanvas();
        c2j_em->cd();
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("MET_{e#mu} (GeV)");
        st_tmp->SetTitle("nJets #geq 2, (" + region + ") hyp_type = em");
        lg_tmp->Draw();
        saveCanvas(c2j_em, "st_met_2j_" + region + "_em");

}

