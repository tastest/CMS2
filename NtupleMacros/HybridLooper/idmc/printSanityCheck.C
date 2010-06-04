
void printSanityCheck()
{

    TFile f("histos_eleid_pt20up.root");

    TCanvas *c1 = new TCanvas();
    c1->cd();

    // loose   
    eleidval_h1_hyp_idstudy_classExpSaniLooseId_EB_ee->Draw("HIST");
    eleidval_h1_hyp_idstudy_classExpLooseRecompId_EB_ee->Draw("SAME E1");
    eleidval_h1_hyp_idstudy_classExpLooseRecompId_EB_ee->SetMarkerStyle(24);
    c1->SaveAs("results/val_loose_id_eb.png");

    eleidval_h1_hyp_idstudy_classExpSaniLooseId_EE_ee->Draw("HIST");
    eleidval_h1_hyp_idstudy_classExpLooseRecompId_EE_ee->Draw("SAME E1");
    eleidval_h1_hyp_idstudy_classExpLooseRecompId_EE_ee->SetMarkerStyle(24);
    c1->SaveAs("results/val_loose_id_ee.png");

    eleidval_h1_hyp_idstudy_classExpSaniLooseIso_EB_ee->Draw("HIST");
    eleidval_h1_hyp_idstudy_classExpLooseRecompIso_EB_ee->Draw("SAME E1");
    eleidval_h1_hyp_idstudy_classExpLooseRecompIso_EB_ee->SetMarkerStyle(24);
    c1->SaveAs("results/val_loose_iso_eb.png");

    eleidval_h1_hyp_idstudy_classExpSaniLooseIso_EE_ee->Draw("HIST");
    eleidval_h1_hyp_idstudy_classExpLooseRecompIso_EE_ee->Draw("SAME E1");
    eleidval_h1_hyp_idstudy_classExpLooseRecompIso_EE_ee->SetMarkerStyle(24);
    c1->SaveAs("results/val_loose_iso_ee.png");


    // tight
    eleidval_h1_hyp_idstudy_classExpSaniTightId_EB_ee->Draw("HIST");
    eleidval_h1_hyp_idstudy_classExpTightRecompId_EB_ee->Draw("SAME E1");
    eleidval_h1_hyp_idstudy_classExpTightRecompId_EB_ee->SetMarkerStyle(24);
    c1->SaveAs("results/val_tight_id_eb.png");

    eleidval_h1_hyp_idstudy_classExpSaniTightId_EE_ee->Draw("HIST");
    eleidval_h1_hyp_idstudy_classExpTightRecompId_EE_ee->Draw("SAME E1");
    eleidval_h1_hyp_idstudy_classExpTightRecompId_EE_ee->SetMarkerStyle(24);
    c1->SaveAs("results/val_tight_id_ee.png");

    eleidval_h1_hyp_idstudy_classExpSaniTightIso_EB_ee->Draw("HIST");
    eleidval_h1_hyp_idstudy_classExpTightRecompIso_EB_ee->Draw("SAME E1");
    eleidval_h1_hyp_idstudy_classExpTightRecompIso_EB_ee->SetMarkerStyle(24);
    c1->SaveAs("results/val_tight_iso_eb.png");

    eleidval_h1_hyp_idstudy_classExpSaniTightIso_EE_ee->Draw("HIST");
    eleidval_h1_hyp_idstudy_classExpTightRecompIso_EE_ee->Draw("SAME E1");
    eleidval_h1_hyp_idstudy_classExpTightRecompIso_EE_ee->SetMarkerStyle(24);
    c1->SaveAs("results/val_tight_iso_ee.png");




}

