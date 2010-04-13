{

     gROOT->ProcessLine("gStyle->SetOptStat(1100111)");

     TChain* chain = new TChain("tree");
     chain->Add("/tas03/disk01/whunt/baby/*.root");

     TCut wselection("iso < 0.2 && pfmet > 25 && pfmet < 80 && abs(d0corr) < 0.4");
     TCut welSelection( wselection + "eormu == 11");
     TCut wmuSelection( wselection + "eormu == 13");

     TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 400);
     canvas->cd();

     // make some general plots without additional selections that help guide further investigation
     unsigned int baselineRunNumber = 133000;
     chain->Draw("run-133000>>candsPerRun(10001, -500.5, 500.5)");
     canvas->Print("plots/candsPerRun.png");

     chain->Draw("pfmet:iso>>isomet(100, 0., 1., 50, 0., 100.)", "","box");
     canvas->Print("plots/pfmetVsIso.png");

     chain->Draw("pfmet:pt>>ptmet(50, 0., 100., 50, 0., 100.)", "","box");
     canvas->Print("plots/pfmetVsPt.png");

     chain->Draw("eormu>>numLepCands(21, -0.5, 20.5)");
     canvas->Print("plots/numLepCands.png");

     // now make some basic plots for events after selection
     chain->Draw("pt>>pt(50, 0., 100.)", wselection);
     canvas->Print("plots/pt.png");

     chain->Draw("d0corr>>d0corr(20, -0.1, 0.1)", wselection);
     canvas->Print("plots/d0corr.png");

     chain->Draw("njets-Sum$(drjet<0.5)>>njets(4, -0.5, 3.5)", wselection + "drjet > 0.5");
     canvas->Print("plots/njets.png");

     // now make some basic electron plots for events passing selection
     chain->Draw("e_cand01>>cand01(2, -0.5, 1.5)", welSelection);
     canvas->Print("plots/e_cand01.png");

     chain->Draw("e_eopin>>eopin(50, 0., 5.)", welSelection);
     canvas->Print("plots/e_eopin.png");

     chain->Draw("e_hoe>>hoe(50, 0., 1.)", welSelection);
     canvas->Print("plots/e_hoe.png");

     chain->Draw("e_dphiin>>dphiin(100, -0.5., 0.5.)", welSelection);
     canvas->Print("plots/e_dphiin.png");

     chain->Draw("e_detain>>detain(100, -0.1., 0.1.)", welSelection);
     canvas->Print("plots/e_detain.png");
   
     chain->Draw("e_eMe55>>eMe55(50, 0., 1.5)", welSelection);
     canvas->Print("plots/e_eMe55.png");

     chain->Draw("e_nmHits>>nmHits(11, -0.5, 10.5)", welSelection);
     canvas->Print("plots/e_nmHits.png");

     chain->Draw("type>>eltype(21, -0.5, 20.5)", welSelection);
     canvas->Print("plots/el_type.png");

     // now make some basic muon plots for events passing selection
     chain->Draw("mu_muonid>>muonid(2, -0.5, 1.5)", wmuSelection);
     canvas->Print("plots/mu_muonid.png");

     chain->Draw("min(max(mu_gfitchi2, -1.), 199.999)>>gfitchi2(100, 0., 200.)", wmuSelection);
     canvas->Print("plots/mu_gfitchi2.png");

     chain->Draw("type>>mutype(21, -0.5, 20.5)", wmuSelection);
     canvas->Print("plots/mu_type.png");

     TTreePlayer *tp = (TTreePlayer*)chain->GetPlayer();
     tp->SetScanRedirect(kTRUE);
     tp->SetScanFileName("wcands.txt");
     chain->Scan("run:ls:evt:pfmet:njets:jet1pt:dphimetjet:eormu:type:pt:iso:d0corr:dphimet:drjet:mt:mu_muonid:mu_goodmask:mu_gfitchi2:e_cand01:e_eopin:e_hoe:e_dphiin:e_detain:e_eMe55:e_nmHits", wselection);
}
