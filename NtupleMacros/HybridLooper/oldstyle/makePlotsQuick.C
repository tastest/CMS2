
{

        TFile f("histos_data.root", "READ");
        
        TCanvas *c1 = new TCanvas();
        c1->cd();

        // electrons
        //

        whunt_ele_nm1_r19_all->Draw();
        c1->SaveAs("r19_nm1.png");

        whunt_ele_nm1nor19_pfmet_all->Draw();
        whunt_ele_nm1nor19_tcmet_all->Draw("SAME");
        whunt_ele_nm1nor19_tcmet_all->SetLineColor(kBlue);
        c1->SaveAs("met_nm1nor19.png");

        whunt_ele_nm1_pfmet_all->Draw();
        whunt_ele_nm1_tcmet_all->Draw("SAME");
        whunt_ele_nm1_tcmet_all->SetLineColor(kBlue);
        c1->SaveAs("met_nm1.png");

        whunt_ele_nm1nor19_pfmetratio_all->Draw();
        whunt_ele_nm1nor19_tcmetratio_all->Draw("SAME");
        whunt_ele_nm1nor19_tcmetratio_all->SetLineColor(kBlue);
        c1->SaveAs("metratio_nm1nor19.png");

        whunt_ele_nm1_pfmetratio_all->Draw();
        whunt_ele_nm1_tcmetratio_all->Draw("SAME");
        whunt_ele_nm1_tcmetratio_all->SetLineColor(kBlue);
        c1->SaveAs("metratio_nm1.png");

        whunt_ele_selected_pfmetratio_all->Draw();
        whunt_ele_selected_tcmetratio_all->Draw("SAME");
        whunt_ele_selected_tcmetratio_all->SetLineColor(kBlue);
        c1->SaveAs("metratio_selected.png");

        whunt_ele_selected_pfmet_all->Draw();
        whunt_ele_selected_tcmet_all->Draw("SAME");
        whunt_ele_selected_tcmet_all->SetLineColor(kBlue);
        c1->SaveAs("met_selected.png");

        whunt_ele_selected_pt_all->Draw();
        c1->SaveAs("pt_selected.png");

        whunt_ele_selected_eta_all->Draw();
        c1->SaveAs("eta_selected.png");

        whunt_ele_selected_phi_all->Draw();
        c1->SaveAs("phi_selected.png");

        // 
        // muons

        whunt_mu_selected_pfmetratio_all->Draw();
        whunt_mu_selected_tcmetratio_all->Draw("SAME");
        whunt_mu_selected_tcmetratio_all->SetLineColor(kBlue);
        c1->SaveAs("mu_metratio_selected.png");

        whunt_mu_selected_pfmet_all->Draw();
        whunt_mu_selected_tcmet_all->Draw("SAME");
        whunt_mu_selected_tcmet_all->SetLineColor(kBlue);
        c1->SaveAs("mu_met_selected.png");

        whunt_mu_selected_pt_all->Draw();
        c1->SaveAs("mu_pt_selected.png");

        whunt_mu_selected_eta_all->Draw();
        c1->SaveAs("mu_eta_selected.png");

        whunt_mu_selected_phi_all->Draw();
        c1->SaveAs("mu_phi_selected.png");

        

}


