#include "Utilities.h"
#include "HistogramUtilities.h"
#include "DataSource.h"

void printMassStacks(TString met)
{

        gROOT->ProcessLine(".L ~/tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");
        gStyle->SetOptTitle(1);

	TString metName = "_MET35";
	if (met == "45") metName = "_MET45";
        if (met == "20") metName = "_MET20";
        if (met == "25") metName = "_MET25";
        if (met == "30") metName = "_MET30";
        if (met == "45_INCL") metName = "_MET45_INCL";

	TFile f("DYEstResults_ForWW" + metName + ".root");
	HistogramUtilities *hUtil = new HistogramUtilities("DYEstResults_ForWW" + metName + ".root", 20.0, true, 0.1);

	gROOT->cd();

        // now make the mass stacks

        AnaHist anaHist;
        TLegend *lg_tmp = anaHist.getLegend(f);

        THStack *st_tmp = hUtil->Stack(f, "0j", "mm");
        TCanvas *c0j_mm = new TCanvas();
        c0j_mm->cd();
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("M_{#mu#mu} (GeV/c^{2})");
        st_tmp->SetTitle("nJets = 0, hyp_type = mm");
        lg_tmp->Draw();
        saveCanvas(c0j_mm, "st_mass_v2_0j_" + met + "_mm");

        THStack *st_tmp = anaHist.getMassStack(f, "1j", "mm");
        TCanvas *c1j_mm = new TCanvas();
        c1j_mm->cd();
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("M_{#mu#mu} (GeV/c^{2})");
        st_tmp->SetTitle("nJets = 1, hyp_type = mm");
        lg_tmp->Draw();
        saveCanvas(c1j_mm, "st_mass_v2_1j_" + met + "_mm");

        THStack *st_tmp = anaHist.getMassStack(f, "2j", "mm");
        TCanvas *c2j_mm = new TCanvas();
        c2j_mm->cd();
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("M_{#mu#mu} (GeV/c^{2})");
        st_tmp->SetTitle("nJets #geq 2, hyp_type = mm");
        lg_tmp->Draw();
        saveCanvas(c2j_mm, "st_mass_v2_2j_" + met + "_mm");

        THStack *st_tmp = anaHist.getMassStack(f, "0j", "ee");
        TCanvas *c0j_ee = new TCanvas();
        c0j_ee->cd();
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
        st_tmp->SetTitle("nJets = 0, hyp_type = ee");
        lg_tmp->Draw();
        saveCanvas(c0j_ee, "st_mass_v2_0j_" + met + "_ee");

        THStack *st_tmp = anaHist.getMassStack(f, "1j", "ee");
        TCanvas *c1j_ee = new TCanvas();
        c1j_ee->cd();
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
        st_tmp->SetTitle("nJets = 1, hyp_type = ee");
        lg_tmp->Draw();
        saveCanvas(c1j_ee, "st_mass_v2_1j_" + met + "_ee");

        THStack *st_tmp = anaHist.getMassStack(f, "2j", "ee");
        TCanvas *c2j_ee = new TCanvas(); 
        c2j_ee->cd();
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
        st_tmp->SetTitle("nJets #geq 2, hyp_type = ee");
        lg_tmp->Draw();
        saveCanvas(c2j_ee, "st_mass_v2_2j_" + met + "_ee");

        // now for em

        THStack *st_tmp = anaHist.getMassStack(f, "0j", "em");
        TCanvas *c0j_em = new TCanvas();
        c0j_em->cd();
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("M_{e#mu} (GeV/c^{2})");
        st_tmp->SetTitle("nJets = 0, hyp_type = em");
        lg_tmp->Draw();
        saveCanvas(c0j_em, "st_mass_v2_0j_" + met + "_em");

        THStack *st_tmp = anaHist.getMassStack(f, "1j", "em");
        TCanvas *c1j_em = new TCanvas();
        c1j_em->cd();
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("M_{e#mu} (GeV/c^{2})");
        st_tmp->SetTitle("nJets = 1, hyp_type = em");
        lg_tmp->Draw();
        saveCanvas(c1j_em, "st_mass_v2_1j_" + met + "_em");

        THStack *st_tmp = anaHist.getMassStack(f, "2j", "em");
        TCanvas *c2j_em = new TCanvas();
        c2j_em->cd();
        st_tmp->Draw("HIST");
        st_tmp->GetXaxis()->SetTitle("M_{e#mu} (GeV/c^{2})");
        st_tmp->SetTitle("nJets #geq 2, hyp_type = em");
        lg_tmp->Draw();
        saveCanvas(c2j_em, "st_mass_v2_2j_" + met + "_em");

	delete hUtil;

}


