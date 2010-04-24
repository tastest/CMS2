{
    gROOT->ProcessLine("gStyle->SetOptStat(1100111)");

    TChain* chain = new TChain("tree");
    chain->Add("/tas03/disk01/dilephunt/baby/*.root");

    TCut zselection("pt1 > 20. && pt2 > 20.");                  // two pt > 20 leptons
    zselection   += "iso1 < 0.2 || iso2 < 0.2";                 // at least one is "well" isolated
    zselection   += "iso1 < 0.4 && iso2 < 0.4";                 // and both are somewhat isolated
    zselection   += "abs(eormu1) == abs(eormu2)";               // same flavor
    zselection   += "abs(d0corr1) < 0.4 && abs(d0corr2) < 0.4"; // humble impact parameters
    //zselection   += "(abs(eormu1) == 13 && (type1&6) == 6) || abs(eormu1) == 11"; // global && tracker muon
    //zselection   += "(abs(eormu2) == 13 && (type2&6) == 6) || abs(eormu2) == 11"; // global && tracker muon

    TCut topselection("(pt1 > 20. || pt2 > 20.) && pt1 > 10. && pt2 > 10."); // 20/10
    topselection   += "iso1 < 0.2 || iso2 < 0.2";                            // at least one is "well" isolated
    topselection   += "iso1 < 0.4 && iso2 < 0.4           ";                 // and both are somewhat isolated
    topselection   += "abs(d0corr1) < 0.4 && abs(d0corr2) < 0.4";            // humble impact parameters
    topselection   += "((hyp_type == 0 || hyp_type == 3) && (pfmet > 10. || tcmet > 10.)) || (hyp_type == 1 || hyp_type == 2)"; // met cuts for ee,mumu only
    topselection   += "(abs(eormu1) == 13 && (type1&6) == 6) || abs(eormu1) == 11"; // global && tracker muon
    topselection   += "(abs(eormu2) == 13 && (type2&6) == 6) || abs(eormu2) == 11"; // global && tracker muon

    TCut zeeSelection  = zselection   + "hyp_type == 3";
    TCut zmmSelection  = zselection   + "hyp_type == 0";

    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 400);
    canvas->cd();

    // make some general plots without additional selections that help guide further investigation
    unsigned int baselineRunNumber = 133000;
    chain->Draw("run-133000>>candsPerRun(10001, -500.5, 500.5)");
    canvas->Print("plots/candsPerRun.png");

    chain->Draw("pfmet:iso1>>isomet1(100, 0., 1., 50, 0., 100.)", "","box");
    chain->Draw("pfmet:iso2>>isomet2(100, 0., 1., 50, 0., 100.)", "","box");
    TH2F* h2 = (TH2F*)isomet1->Clone();
    h2->Add(isomet2, 1);
    h2->Draw("box");
    canvas->Print("plots/pfmetVsIso.png");

    chain->Draw("pfmet:pt1>>pt1(100, 0., 100., 50, 0., 100.)", "","box");
    chain->Draw("pfmet:pt2>>pt2(100, 0., 100., 50, 0., 100.)", "","box");
    h2 = (TH2F*)pt1->Clone();
    h2->Add(pt2, 1);
    h2->Draw("box");
    canvas->Print("plots/pfmetVsPt.png");

    chain->Draw("hyp_type>>ztype(4, -0.5, 3.5)", zselection);
    canvas->Print("plots/zhyptype.png");

    chain->Draw("hyp_type>>ttype(4, -0.5, 3.5)", topselection);
    canvas->Print("plots/tophyptype.png");

    chain->Draw("pfmet>>zpfmet(50, 0., 100.)", zselection);
    canvas->Print("plots/zpfmet.png");

    chain->Draw("pfmet>>tpfmet(50, 0., 100.)", topselection);
    canvas->Print("plots/toppfmet.png");

    chain->Draw("tcmet>>ztcmet(50, 0., 100.)", zselection);
    canvas->Print("plots/ztcmet.png");

    chain->Draw("tcmet>>ttcmet(50, 0., 100.)", topselection);
    canvas->Print("plots/toptcmet.png");

    chain->Draw("pt1>>zpt1(50, 0., 100.)", zselection);
    chain->Draw("pt2>>zpt2(50, 0., 100.)", zselection);
    TH1F* h1 = (TH1F*)zpt1->Clone();
    h1->Add(zpt2, 1);
    h1->Draw();
    canvas->Print("plots/zpt.png");

    chain->Draw("pt1>>tpt1(50, 0., 100.)", topselection);
    chain->Draw("pt2>>tpt2(50, 0., 100.)", topselection);
    h1 = (TH1F*)tpt1->Clone();
    h1->Add(tpt2, 1);
    h1->Draw();
    canvas->Print("plots/toppt.png");

    chain->Draw("eta1>>zeta1(60, -3., 3.)", zselection);
    chain->Draw("eta2>>zeta2(60, -3., 3.)", zselection);
    h1 = (TH1F*)zeta1->Clone();
    h1->Add(zeta2, 1);
    h1->Draw();
    canvas->Print("plots/zeta.png");
    
    chain->Draw("eta1>>teta1(60, -3., 3.)", topselection);
    chain->Draw("eta2>>teta2(60, -3., 3.)", topselection);
    h1 = (TH1F*)teta1->Clone();
    h1->Add(teta2, 1);
    h1->Draw();
    canvas->Print("plots/topeta.png");
        
    chain->Draw("phi1>>zphi1(64, -3.2, 3.2)", zselection);
    chain->Draw("phi2>>zphi2(64, -3.2, 3.2)", zselection);
    h1 = (TH1F*)zphi1->Clone();
    h1->Add(zphi2, 1);
    h1->Draw();
    canvas->Print("plots/zphi.png");

    chain->Draw("phi1>>tphi1(64, -3.2, 3.2)", topselection);
    chain->Draw("phi2>>tphi2(64, -3.2, 3.2)", topselection);
    h1 = (TH1F*)tphi1->Clone();
    h1->Add(tphi2, 1);
    h1->Draw();
    canvas->Print("plots/topphi.png");
    
    chain->Draw("iso1>>ziso1(50, 0., 1.)", zselection);
    chain->Draw("iso2>>ziso2(50, 0., 1.)", zselection);
    h1 = (TH1F*)ziso1->Clone();
    h1->Add(ziso2, 1);
    h1->Draw();
    canvas->Print("plots/ziso.png");
    
    chain->Draw("iso1>>tiso1(50, 0., 1.)", topselection);
    chain->Draw("iso2>>tiso2(50, 0., 1.)", topselection);
    h1 = (TH1F*)tiso1->Clone();
    h1->Add(tiso2, 1);
    h1->Draw();
    canvas->Print("plots/topiso.png");
    
    chain->Draw("d0vtx1>>zd0vtx1(50, -0.5, 0.5)", zselection);
    chain->Draw("d0vtx2>>zd0vtx2(50, -0.5, 0.5)", zselection);
    h1 = (TH1F*)zd0vtx1->Clone();
    h1->Add(zd0vtx2, 1);
    h1->Draw();
    canvas->Print("plots/zd0vtx.png");
    
    chain->Draw("d0vtx1>>td0vtx1(50, -0.5, 0.5)", topselection);
    chain->Draw("d0vtx2>>td0vtx2(50, -0.5, 0.5)", topselection);
    h1 = (TH1F*)td0vtx1->Clone();
    h1->Add(td0vtx2, 1);
    h1->Draw();
    canvas->Print("plots/topd0vtx.png");
    
    chain->Draw("njets-Sum$(drjet1<0.5)-Sum$(drjet2<0.5)>>znjets(4, -0.5, 3.5)", zselection);
    canvas->Print("plots/znjets.png");

    chain->Draw("njets-Sum$(drjet1<0.5)-Sum$(drjet2<0.5)>>topnjets(4, -0.5, 3.5)", topselection);
    canvas->Print("plots/topnjets.png");

    chain->Draw("mass>>zmass(40, 36., 116.)", zselection);
    canvas->Print("plots/zmass.png");

    chain->Draw("mass>>topmass(80, 0., 160.)", topselection);
    canvas->Print("plots/topmass.png");

    chain->Draw("mu1_muonid>>zmuonid1(2, -0.5, 1.5)", zmmSelection);
    chain->Draw("mu2_muonid>>zmuonid2(2, -0.5, 1.5)", zmmSelection);
    h1 = (TH1F*)zmuonid1->Clone();
    h1->Add(zmuonid2, 1);
    h1->Draw();
    canvas->Print("plots/zmuonid.png");
    
    chain->Draw("mu1_muonid>>tmuonid1(2, -0.5, 1.5)", topselection + "hyp_type == 0 || hyp_type == 1");
    chain->Draw("mu2_muonid>>tmuonid2(2, -0.5, 1.5)", topselection + "hyp_type == 0 || hyp_type == 2");
    h1 = (TH1F*)tmuonid1->Clone();
    h1->Add(tmuonid2, 1);
    h1->Draw();
    canvas->Print("plots/topmuonid.png");
    
    chain->Draw("type1>>ztype1(21, -0.5, 20.5)", zmmSelection);
    chain->Draw("type2>>ztype2(21, -0.5, 20.5)", zmmSelection);
    h1 = (TH1F*)ztype1->Clone();
    h1->Add(ztype2, 1);
    h1->Draw();
    canvas->Print("plots/zmutype.png");
    
    chain->Draw("type1>>ttype1(21, -0.5, 20.5)", topselection + "hyp_type == 0 || hyp_type == 1");
    chain->Draw("type2>>ttype2(21, -0.5, 20.5)", topselection + "hyp_type == 0 || hyp_type == 2");
    h1 = (TH1F*)ttype1->Clone();
    h1->Add(ttype2, 1);
    h1->Draw();
    canvas->Print("plots/topmutype.png");

    chain->Draw("e1_eopin>>zeopin1(50, 0., 5.)", zeeSelection);
    chain->Draw("e2_eopin>>zeopin2(50, 0., 5.)", zeeSelection);
    h1 = (TH1F*)zeopin1->Clone();
    h1->Add(zeopin2, 1);
    h1->Draw();
    canvas->Print("plots/ze_eopin.png");

    chain->Draw("e1_eopin>>teopin1(50, 0., 5.)", topselection + "hyp_type == 2 || hyp_type == 3");
    chain->Draw("e2_eopin>>teopin2(50, 0., 5.)", topselection + "hyp_type == 1 || hyp_type == 3");
    h1 = (TH1F*)teopin1->Clone();
    h1->Add(teopin2, 1);
    h1->Draw();
    canvas->Print("plots/tope_eopin.png");

    chain->Draw("e1_hoe>>zhoe1(50, 0., 1.)", zeeSelection);
    chain->Draw("e2_hoe>>zhoe2(50, 0., 1.)", zeeSelection);
    h1 = (TH1F*)zhoe1->Clone();
    h1->Add(zhoe2, 1);
    h1->Draw();
    canvas->Print("plots/ze_hoe.png");

    chain->Draw("e1_hoe>>thoe1(50, 0., 1.)", topselection + "hyp_type == 2 || hyp_type == 3");
    chain->Draw("e2_hoe>>thoe2(50, 0., 1.)", topselection + "hyp_type == 1 || hyp_type == 3");
    h1 = (TH1F*)thoe1->Clone();
    h1->Add(thoe2, 1);
    h1->Draw();
    canvas->Print("plots/tope_hoe.png");

    chain->Draw("e1_dphiin>>zdphiin1(100, -0.5, 0.5)", zeeSelection);
    chain->Draw("e2_dphiin>>zdphiin2(100, -0.5, 0.5)", zeeSelection);
    h1 = (TH1F*)zdphiin1->Clone();
    h1->Add(zdphiin2, 1);
    h1->Draw();
    canvas->Print("plots/ze_dphiin.png");

    chain->Draw("e1_dphiin>>tdphiin1(100, -0.5, 0.5)", topselection + "hyp_type == 2 || hyp_type == 3");
    chain->Draw("e2_dphiin>>tdphiin2(100, -0.5, 0.5)", topselection + "hyp_type == 1 || hyp_type == 3");
    h1 = (TH1F*)tdphiin1->Clone();
    h1->Add(tdphiin2, 1);
    h1->Draw();
    canvas->Print("plots/tope_dphiin.png");

    chain->Draw("e1_detain>>zdetain1(100, -0.1, 0.1)", zeeSelection);
    chain->Draw("e2_detain>>zdetain2(100, -0.1, 0.1)", zeeSelection);
    h1 = (TH1F*)zdetain1->Clone();
    h1->Add(zdetain2, 1);
    h1->Draw();
    canvas->Print("plots/ze_detain.png");

    chain->Draw("e1_detain>>tdetain1(100, -0.1, 0.1)", topselection + "hyp_type == 2 || hyp_type == 3");
    chain->Draw("e2_detain>>tdetain2(100, -0.1, 0.1)", topselection + "hyp_type == 1 || hyp_type == 3");
    h1 = (TH1F*)tdetain1->Clone();
    h1->Add(tdetain2, 1);
    h1->Draw();
    canvas->Print("plots/tope_detain.png");

    chain->Draw("e1_eMe55>>zeMe551(50, 0., 1.5)", zeeSelection);
    chain->Draw("e2_eMe55>>zeMe552(50, 0., 1.5)", zeeSelection);
    h1 = (TH1F*)zeMe551->Clone();
    h1->Add(zeMe552, 1);
    h1->Draw();
    canvas->Print("plots/ze_eMe55.png");

    chain->Draw("e1_eMe55>>teMe551(50, 0., 1.5)", topselection + "hyp_type == 2 || hyp_type == 3");
    chain->Draw("e2_eMe55>>teMe552(50, 0., 1.5)", topselection + "hyp_type == 1 || hyp_type == 3");
    h1 = (TH1F*)teMe551->Clone();
    h1->Add(teMe552, 1);
    h1->Draw();
    canvas->Print("plots/tope_eMe55.png");

    chain->Draw("e1_nmHits>>znmHits1(11, -0.5, 10.5)", zeeSelection);
    chain->Draw("e2_nmHits>>znmHits2(11, -0.5, 10.5)", zeeSelection);
    h1 = (TH1F*)znmHits1->Clone();
    h1->Add(znmHits2, 1);
    h1->Draw();
    canvas->Print("plots/ze_nmHits.png");

    chain->Draw("e1_nmHits>>tnmHits1(11, -0.5, 10.5)", topselection + "hyp_type == 2 || hyp_type == 3");
    chain->Draw("e2_nmHits>>tnmHits2(11, -0.5, 10.5)", topselection + "hyp_type == 1 || hyp_type == 3");
    h1 = (TH1F*)tnmHits1->Clone();
    h1->Add(tnmHits2, 1);
    h1->Draw();
    canvas->Print("plots/tope_nmHits.png");

    chain->Draw("e1_cand01>>z1cand01(2, -0.5, 1.5)", zeeSelection);
    chain->Draw("e2_cand01>>z2cand01(2, -0.5, 1.5)", zeeSelection);
    h1 = (TH1F*)z1cand01->Clone();
    h1->Add(z2cand01, 1);
    h1->Draw();
    canvas->Print("plots/ze_cand01.png");

    chain->Draw("e1_cand01>>t1cand01(2, -0.5, 1.5)", topselection + "hyp_type == 2 || hyp_type == 3");
    chain->Draw("e2_cand01>>t2cand01(2, -0.5, 1.5)", topselection + "hyp_type == 1 || hyp_type == 3");
    h1 = (TH1F*)t1cand01->Clone();
    h1->Add(t2cand01, 1);
    h1->Draw();
    canvas->Print("plots/tope_cand01.png");

    chain->Draw("type1>>zeltype1(21, -0.5, 20.5)", zeeSelection);
    chain->Draw("type2>>zeltype2(21, -0.5, 20.5)", zeeSelection);
    h1 = (TH1F*)zeltype1->Clone();
    h1->Add(zeltype2, 1);
    h1->Draw();
    canvas->Print("plots/ze_type.png");

    chain->Draw("type1>>teltype1(21, -0.5, 20.5)", topselection + "hyp_type == 2 || hyp_type == 3");
    chain->Draw("type2>>teltype2(21, -0.5, 20.5)", topselection + "hyp_type == 1 || hyp_type == 3");
    h1 = (TH1F*)teltype1->Clone();
    h1->Add(teltype2, 1);
    h1->Draw();
    canvas->Print("plots/tope_type.png");

    TTreePlayer *tp = (TTreePlayer*)chain->GetPlayer();
    tp->SetScanRedirect(kTRUE);
    tp->SetScanFileName("zcands.txt");
    chain->Scan("*", zselection, "colsize=50");
    tp->SetScanFileName("topcands.txt");
    chain->Scan("*", topselection, "colsize=50");
}
