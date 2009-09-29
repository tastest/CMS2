
#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
#include "Tools/Utilities.h"
#include "Tools/histtools.cc"

#include "TROOT.h"

#include "TGraphAsymmErrors.h"

        const static sources_t theSources_31X =
//		(1ll << H_ZZ);
                (1ll << H_WENU_7TeV);        
//                (1ll << H_EM30_80)       |
//                (1ll << H_BC30_80);

        const static sources_t theSignal_31X =
//		(1ll << H_ZZ);
                (1ll << H_WENU_7TeV);

        const static sources_t theBackground_31X =
		(1ll << H_WENU_7TeV);
//		(1ll << H_ZZ);
//                (1ll << H_EM30_80)       |
//                (1ll << H_BC30_80);

        const static sources_t theSources_22X =
                (1ll << H_QCD30) 	|
                (1ll << H_QCD80)        |
		(1ll << H_WJET_ALP);

	const static sources_t theSignal_22X = 
                (1ll << H_WJET_ALP);

	const static sources_t theBackground_22X = 
                (1ll << H_QCD30) 	|
                (1ll << H_QCD80);

	// for 2_2_X
	const static sources_t &theSignal = theSignal_31X;
	const static sources_t &theBackground = theBackground_31X;
	const static sources_t &theSources = theSources_31X;

void plotEff(HistogramUtilities &h1, TString name, TString det, bool ascending)
{

	TH1F *h1_signal = h1.getHistogram(theSignal, name, "", det);
        TH1F *h1_background = h1.getHistogram(theBackground, name, "", det);

	// latter bool is "ascending"
        TGraph *gr = (TGraph*)(eff_rej(*h1_signal, *h1_background, true, ascending).Clone());
        gr->SetTitle(name + ";Signal;Background");
	gr->SetMarkerStyle(23);

        TCanvas *c = new TCanvas();
        c->cd();
        gr->Draw("AP");
	Utilities::saveCanvas(c, "results/eff_" + name + "_" + det);

	delete gr;
	delete h1_signal;
	delete h1_background;
	delete c;

}

void plotEffVar(HistogramUtilities &h1, TString name, Int_t rebin = 1)
{

        TCanvas *c = new TCanvas();
        TH1F *h1_eff = 0;
        TH1F *h1_total = 0;
        TH1F *h1_pass = 0;
	//TGraphAsymmErrors *gr_eff = new TGraphAsymmErrors();

        h1_total = h1.getHistogram(theSignal, name, "", "denom", 1, "");
        h1_pass = h1.getHistogram(theSignal, name, "", "numer", 1, "");
        h1_total->Rebin(rebin);
        h1_pass->Rebin(rebin);
	//gr_eff->BayesDivide(h1_pass, h1_total);
        h1_eff = (TH1F*)h1_pass->Clone();
        h1_eff->Reset();
        h1_eff->Divide(h1_pass, h1_total, 1.0, 1.0, "B");
        h1_eff->SetName(h1_pass->GetName());
        c->SetName(TString("c_s_") + h1_pass->GetName());
        c->cd();
        h1_eff->SetLineWidth(2);
        h1_eff->SetMarkerStyle(20);
	h1_eff->GetYaxis()->SetRangeUser(0, 1.1);
        h1_eff->Draw("E1");
	//gr_eff->Draw("AP");
        //gr_eff->GetYaxis()->SetRangeUser(0, 1.1);
        Utilities::saveCanvas(c, "results/effVar_s_" + name);


        h1_total = h1.getHistogram(theBackground, name, "", "denom", 1, "");
        h1_pass = h1.getHistogram(theBackground, name, "", "numer", 1, "");
	h1_total->Rebin(rebin);
	h1_pass->Rebin(rebin);
        //gr_eff->BayesDivide(h1_pass, h1_total);
        h1_eff->Reset();
        h1_eff->Divide(h1_pass, h1_total, 1.0, 1.0, "B");
        h1_eff->SetName(h1_pass->GetName());
        c->SetName(TString("c_bg_") + h1_pass->GetName());
        c->cd();
	h1_eff->SetLineWidth(2);
	h1_eff->SetMarkerStyle(20);
        h1_eff->GetYaxis()->SetRangeUser(0, 1.1);
        h1_eff->Draw("E1");
	//gr_eff->Draw("AP");
        //gr_eff->GetYaxis()->SetRangeUser(0, 1.1);
        Utilities::saveCanvas(c, "results/effVar_bg_" + name);

        delete c;
        delete h1_eff;
        delete h1_pass;
        delete h1_total;
	//delete gr_eff;
}

void test()
{
        gROOT->ProcessLine(".L ~/tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");
        HistogramUtilities h1("Results.root", 0.001);

        plotEffVar(h1, "dEtaIn_pt_ee", 4);
        plotEffVar(h1, "dEtaIn_eta_ee");
        plotEffVar(h1, "dEtaIn_phi_ee");
        plotEffVar(h1, "dEtaIn_pt_eb", 4);
        plotEffVar(h1, "dEtaIn_eta_eb");
        plotEffVar(h1, "dEtaIn_phi_eb");

        plotEffVar(h1, "dPhiIn_pt_ee", 4);
        plotEffVar(h1, "dPhiIn_eta_ee");
        plotEffVar(h1, "dPhiIn_phi_ee");
        plotEffVar(h1, "dPhiIn_pt_eb", 4);
        plotEffVar(h1, "dPhiIn_eta_eb");
        plotEffVar(h1, "dPhiIn_phi_eb");

        plotEffVar(h1, "hoe_pt_ee", 4);
        plotEffVar(h1, "hoe_eta_ee");
        plotEffVar(h1, "hoe_phi_ee");
        plotEffVar(h1, "hoe_pt_eb", 4);
        plotEffVar(h1, "hoe_eta_eb");
        plotEffVar(h1, "hoe_phi_eb");

        plotEffVar(h1, "sieie_pt_ee", 4);
        plotEffVar(h1, "sieie_eta_ee");
        plotEffVar(h1, "sieie_phi_ee");
        plotEffVar(h1, "sieie_pt_eb", 4);
        plotEffVar(h1, "sieie_eta_eb");
        plotEffVar(h1, "sieie_phi_eb");
}

void plotResults(TString det)
{

	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");

	// luminorm for 1pb-1
	HistogramUtilities h1("Results.root", 0.001);
	THStack *st_pt = h1.getStack(theSources, "h1_pt", "", det, 2);

        THStack *st_eta = h1.getStack(theSources, "h1_eta", "", det);
	TLegend *lg_all = h1.getLegend(theSources, "h1_pt", "", det);

        THStack *st_ecalIso03 = h1.getStack(theSources, "h1_ecalIso03", "", det);
        THStack *st_hcalIso03 = h1.getStack(theSources, "h1_hcalIso03", "", det);
        THStack *st_tkIso03 = h1.getStack(theSources, "h1_tkIso03", "", det);
        //THStack *st_esJuraIso03 = h1.getStack(theSources, "h1_esJuraIso03", "", det);
        THStack *st_wwIso = h1.getStack(theSources, "h1_wwIso", "", det, 4);

	//plotEff(h1, "h1_esJuraIso03", det, true);
        //plotEff(h1, "h1_ecalIso03", det, true);
        //plotEff(h1, "h1_hcalIso03", det, true);
        //plotEff(h1, "h1_tkIso03", det, true);
	//plotEff(h1, "h1_wwIso", det, false);

	TCanvas *c1 = new TCanvas();
	c1->cd();
	st_pt->Draw();
	lg_all->Draw();
	Utilities::saveCanvas(c1, "results/pt_" + det);	

        TCanvas *c2 = new TCanvas();
        c2->cd();
        st_eta->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c2, "results/eta_" + det);

        TCanvas *c3 = new TCanvas();
        c3->cd();
        st_ecalIso03->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c3, "results/ecalIso03_" + det);

        TCanvas *c4 = new TCanvas();
        c4->cd();
        st_hcalIso03->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c4, "results/hcalIso03_" + det);

        TCanvas *c5 = new TCanvas();
        c5->cd();
        st_tkIso03->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c5, "results/tkIso03_" + det);

        TCanvas *c6 = new TCanvas();
        c6->cd();
        st_wwIso->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c6, "results/wwIso_" + det);

        //TCanvas *c7 = new TCanvas();
        //c7->cd();
        //st_esJuraIso03->Draw();
        //lg_all->Draw();
        //Utilities::saveCanvas(c7, "results/esJuraIso_" + det);

	// electron id related
	//
        THStack *st_dEtaIn = h1.getStack(theSources, "h1_dEtaIn", "", det, 2);
        THStack *st_dPhiIn = h1.getStack(theSources, "h1_dPhiIn", "", det, 2);
        THStack *st_hoe = h1.getStack(theSources, "h1_hoe", "", det, 4);
        THStack *st_eopIn = h1.getStack(theSources, "h1_eopIn", "", det, 4);
        THStack *st_sigmaIEtaIEta = h1.getStack(theSources, "h1_sigmaIEtaIEta", "", det);
        THStack *st_sigmaIPhiIPhi = h1.getStack(theSources, "h1_sigmaIPhiIPhi", "", det);
        THStack *st_E2x5Norm5x5 = h1.getStack(theSources, "h1_E2x5Norm5x5", "", det, 4);
        THStack *st_E1x5Norm5x5 = h1.getStack(theSources, "h1_E1x5Norm5x5", "", det, 4);

        //plotEff(h1, "h1_hoe", det, true);
	
        TCanvas *c_id1 = new TCanvas();
        c_id1->cd();
        st_dEtaIn->Draw();
        lg_all->Draw();
	c_id1->Update();
        Utilities::saveCanvas(c_id1, "results/dEtaIn_" + det);
        
	TCanvas *c_id2 = new TCanvas();
        c_id2->cd();
        st_dPhiIn->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c_id2, "results/dPhiIn_" + det);

        TCanvas *c_id3 = new TCanvas();
        c_id3->cd();
        st_hoe->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c_id3, "results/hoe_" + det);

        TCanvas *c_id4 = new TCanvas();
        c_id4->cd();
        st_sigmaIEtaIEta->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c_id4, "results/sigmaIEtaIEta_" + det);

        TCanvas *c_id5 = new TCanvas();
        c_id5->cd();
        st_sigmaIPhiIPhi->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c_id5, "results/sigmaIPhiIPhi_" + det);

        TCanvas *c_id6 = new TCanvas();
        c_id6->cd();
        st_eopIn->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c_id6, "results/eopIn_" + det);

        TCanvas *c_id7 = new TCanvas();
        c_id7->cd();
        st_E2x5Norm5x5->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c_id7, "results/E2x5Norm5x5_" + det);

        TCanvas *c_id8 = new TCanvas();
        c_id8->cd();
        st_E1x5Norm5x5->Draw();
        lg_all->Draw();
        Utilities::saveCanvas(c_id8, "results/E1x5Norm5x5_" + det);

}

