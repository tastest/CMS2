
void printSanityCheck()
{

    TFile f("histos_eleid_saniv02_pt20up.root");

    TCanvas *c1 = new TCanvas();
    c1->cd();   
    eleidval_h1_hyp_idstudy_classExpSaniLooseId_EB_ee->Draw("HIST");
    eleidval_h1_hyp_idstudy_classExpRecompId_EB_ee->Draw("SAME E1");
    eleidval_h1_hyp_idstudy_classExpRecompId_EB_ee->SetMarkerStyle(24);
    c1->SaveAs("results/val_loose_id_eb.png");

    eleidval_h1_hyp_idstudy_classExpSaniLooseId_EE_ee->Draw("HIST");
    eleidval_h1_hyp_idstudy_classExpRecompId_EE_ee->Draw("SAME E1");
    eleidval_h1_hyp_idstudy_classExpRecompId_EE_ee->SetMarkerStyle(24);
    c1->SaveAs("results/val_loose_id_ee.png");

    eleidval_h1_hyp_idstudy_classExpSaniLooseIso_EB_ee->Draw("HIST");
    eleidval_h1_hyp_idstudy_classExpRecompIso_EB_ee->Draw("SAME E1");
    eleidval_h1_hyp_idstudy_classExpRecompIso_EB_ee->SetMarkerStyle(24);
    c1->SaveAs("results/val_loose_iso_eb.png");

    eleidval_h1_hyp_idstudy_classExpSaniLooseIso_EE_ee->Draw("HIST");
    eleidval_h1_hyp_idstudy_classExpRecompIso_EE_ee->Draw("SAME E1");
    eleidval_h1_hyp_idstudy_classExpRecompIso_EE_ee->SetMarkerStyle(24);
    c1->SaveAs("results/val_loose_iso_ee.png");

}

