
{

	TFile f_wjets("Looper_wjets.root");
	TTree *t_wjets = (TTree*)f_wjets.Get("T1");

	TFile f_em30to80("Looper_QCDEMenrichedPt30to80.root");
        TTree *t_em30to80 = (TTree*)f_em30to80.Get("T1");

	Int_t nBins = 100;
	Float_t binMin = 0.0;
	Float_t binMax = 100.0;
	TH1F *h1_wjets_tcmet = new TH1F("h1_wjets_tcmet", "h1_wjets_tcmet", nBins, binMin, binMax);
		h1_wjets_tcmet->SetLineColor(kRed);
        TH1F *h1_em80to80_tcmet = new TH1F("h1_em30to80_tcmet", "h1_em30to80_tcmet", nBins, binMin, binMax);
		h1_em80to80_tcmet->SetLineColor(kBlack);


	TCanvas *c1 = new TCanvas();
	c1->cd();

	// weight divided by 100 should be 10pb-1
	t_wjets->Draw("evt_tcmet >> h1_wjets_tcmet", "(ele_count > 0 && ele_pt[0] > 20 && ele_hOverE[0] < 0.05)*(evt_weight/100.)");
	t_em30to80->Draw("evt_tcmet >> h1_em30to80_tcmet", "(ele_count > 0 && ele_pt[0] > 20 && ele_hOverE[0] < 0.05)*(evt_weight/100.)");

	h1_em30to80_tcmet->Draw();
	h1_wjets_tcmet->Draw("SAME");



}

