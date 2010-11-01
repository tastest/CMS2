void CompareWWFRwithPrediction(bool doels = true, bool domus = true) {
    gROOT->LoadMacro("eff2.C");
    gStyle->SetOptStat(0);

    //-------------------------------------------------------------------------
    // Load the baby ntuple
    //-------------------------------------------------------------------------

    TChain *ch_el = new TChain("tree");
    if (doels) {
        ch_el->Add("EG.root");
        //ch_el->Add("qcd30.root");
    }

    TChain *ch_mu = new TChain("tree");
    if (domus) {
        ch_mu->Add("Mu.root");
        //ch_mu->Add("inclMu.root");
    }

    //-------------------------------------------------------------------------
    // Load the fake rate histograms
    // Need to know the name of the 
    // fake rate files and histograms
    //-------------------------------------------------------------------------

    loadHist("ww_el_fr_EG.root",0,0,"*fr*");
    //loadHist("ww_el_fr_qcd30.root",0,0,"*fr*");
    loadHist("ww_mu_fr_Mu.root",0,0,"*fr*");
    //loadHist("ww_mu_fr_inclMu.root",0,0,"*fr*");

    TH2F* el_fr_v1_wwV0b   = (TH2F*)gDirectory->Get("el_fr_v1_wwV0b");
    TH2F* el_fr_v2_wwV0b   = (TH2F*)gDirectory->Get("el_fr_v2_wwV0b");
    TH2F* el_fr_v3_wwV0b   = (TH2F*)gDirectory->Get("el_fr_v3_wwV0b");
    TH2F* el_fr_v4_wwV0b   = (TH2F*)gDirectory->Get("el_fr_v4_wwV0b");
    TH2F* mu_fr_fo_wwV0_10 = (TH2F*)gDirectory->Get("mu_fr_fo_wwV0_10");

    //-------------------------------------------------------------------------
    // These are the cuts that define the new selection
    //-------------------------------------------------------------------------

    // A cut against Ws
    TCut notWCut  = "tcmet<20 && mt<25";

    // A pt cut...
    // Remember, we use 35 (at least for muons) to
    // minimize the impact of Ws
    TCut ptCut    = "pt>20";

    // Only consider events with >0 pfjet with pt
    // above a threshold separated by at least dR
    TCut jetCut   = "ptpfj1>5";

    // The trigger selections
    TCut trgCutEl = "el10_lw>1 || el10_sw>1 || el15_lw>1 || el15_sw>1 || el17_sw>1 || el20_sw>1";
    TCut trgCutMu = "mu9>1 || mu11>1 || mu15>1";

    // Numerator selections
    TCut is_el_num_wwV0  = "num_wwV0&&abs(id)==11"+trgCutEl+jetCut+ptCut+notWCut;
    TCut is_el_num_wwV0b = "num_wwV0b&&abs(id)==11"+trgCutEl+jetCut+ptCut+notWCut;

    TCut is_mu_num_wwV0  = "num_wwV0&&abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;

    // Denominator selections
    TCut is_el_v1_wwV0    = "v1_wwV0&&abs(id)==11"+trgCutEl+jetCut+ptCut+notWCut;
    TCut is_el_v2_wwV0    = "v2_wwV0&&abs(id)==11"+trgCutEl+jetCut+ptCut+notWCut;
    TCut is_el_v3_wwV0    = "v3_wwV0&&abs(id)==11"+trgCutEl+jetCut+ptCut+notWCut;
    TCut is_el_v4_wwV0    = "v4_wwV0&&abs(id)==11"+trgCutEl+jetCut+ptCut+notWCut;

    TCut is_el_v1_wwV0b   = "v1_wwV0b&&abs(id)==11"+trgCutEl+jetCut+ptCut+notWCut;
    TCut is_el_v2_wwV0b   = "v2_wwV0b&&abs(id)==11"+trgCutEl+jetCut+ptCut+notWCut;
    TCut is_el_v3_wwV0b   = "v3_wwV0b&&abs(id)==11"+trgCutEl+jetCut+ptCut+notWCut;
    TCut is_el_v4_wwV0b   = "v4_wwV0b&&abs(id)==11"+trgCutEl+jetCut+ptCut+notWCut;

    TCut is_mu_fo_wwV0_04 = "fo_wwV0_04&&abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;
    TCut is_mu_fo_wwV0_10 = "fo_wwV0_10&&abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;

    //-------------------------------------------------------------------------
    // Make three empty histograms with the same
    // binning as the fake rate histogram
    //-------------------------------------------------------------------------

    TH2F* el_num_wwV0b = el_fr_v1_wwV0b->Clone("el_num_wwV0b");
    el_num_wwV0b->Reset();

    TH2F* el_v1_wwV0b_NotNum = el_fr_v1_wwV0b->Clone("el_v1_wwV0b_NotNum");
    el_v1_wwV0b_NotNum->Reset();
    TH2F* el_v2_wwV0b_NotNum = el_fr_v2_wwV0b->Clone("el_v2_wwV0b_NotNum");
    el_v2_wwV0b_NotNum->Reset();
    TH2F* el_v3_wwV0b_NotNum = el_fr_v3_wwV0b->Clone("el_v3_wwV0b_NotNum");
    el_v3_wwV0b_NotNum->Reset();
    TH2F* el_v4_wwV0b_NotNum = el_fr_v4_wwV0b->Clone("el_v4_wwV0b_NotNum");
    el_v4_wwV0b_NotNum->Reset();

    TH2F* pred_el_v1_wwV0b = el_fr_v1_wwV0b->Clone("pred_el_v1_wwV0b");
    pred_el_v1_wwV0b->Reset();
    TH2F* pred_el_v2_wwV0b = el_fr_v2_wwV0b->Clone("pred_el_v2_wwV0b");
    pred_el_v2_wwV0b->Reset();
    TH2F* pred_el_v3_wwV0b = el_fr_v3_wwV0b->Clone("pred_el_v3_wwV0b");
    pred_el_v3_wwV0b->Reset();
    TH2F* pred_el_v4_wwV0b = el_fr_v4_wwV0b->Clone("pred_el_v4_wwV0b");
    pred_el_v4_wwV0b->Reset();

    TH2F* mu_num_wwV0 = mu_fr_fo_wwV0_10->Clone("mu_num_wwV0");
    mu_num_wwV0->Reset();

    TH2F* mu_fo_wwV0_10_NotNum = mu_fr_fo_wwV0_10->Clone("mu_fo_wwV0_10_NotNum");
    mu_fo_wwV0_10_NotNum->Reset();

    TH2F* pred_mu_fo_wwV0_10 = mu_fr_fo_wwV0_10->Clone("pred_mu_fo_wwV0_10");
    pred_mu_fo_wwV0_10->Reset();

    //-------------------------------------------------------------------------
    // Fill two of the three histograms.
    // The 1st one is filled with the numerator
    // The 2nd one is filled with the denominator but not numerator
    // The 3rd one will be filled with the prediction
    //-------------------------------------------------------------------------

    float ptmax = el_num_wwV0b->GetYaxis()->GetXmax();
    cout << "Leptons with pt > " << ptmax   << " will be put in the last bin" << endl; 

    ch_el->Draw(Form("min(pt,%f-0.1):abs(eta)>>el_num_wwV0b",ptmax),is_el_num_wwV0b);
    ch_el->Draw(Form("min(pt,%f-0.1):abs(eta)>>el_v1_wwV0b_NotNum",ptmax),is_el_v1_wwV0b&&!is_el_num_wwV0b);
    ch_el->Draw(Form("min(pt,%f-0.1):abs(eta)>>el_v2_wwV0b_NotNum",ptmax),is_el_v2_wwV0b&&!is_el_num_wwV0b);
    ch_el->Draw(Form("min(pt,%f-0.1):abs(eta)>>el_v3_wwV0b_NotNum",ptmax),is_el_v3_wwV0b&&!is_el_num_wwV0b);
    ch_el->Draw(Form("min(pt,%f-0.1):abs(eta)>>el_v4_wwV0b_NotNum",ptmax),is_el_v4_wwV0b&&!is_el_num_wwV0b);

    ch_mu->Draw(Form("min(pt,%f-0.1):abs(eta)>>mu_num_wwV0",ptmax),is_mu_num_wwV0);
    ch_mu->Draw(Form("min(pt,%f-0.1):abs(eta)>>mu_fo_wwV0_10_NotNum",ptmax),is_mu_fo_wwV0_10&&!is_mu_num_wwV0);

    //-------------------------------------------------------------------------
    // Now we need to get a histogram full of ones
    //-------------------------------------------------------------------------

    TH2F* ones = el_fr_v1_wwV0b->Clone();
    ones->Reset();
    ones->SetTitle("ones");
    ones->SetName("ones");
    for (int ix=1; ix<ones->GetNbinsX()+1; ix++) {
        for (int iy=1; iy<ones->GetNbinsY()+1; iy++) {
            ones->SetBinContent(ix,iy,1.);
        }
    }

    //-------------------------------------------------------------------------
    // Now do the actual calculation and output the results
    //-------------------------------------------------------------------------

    TH2F *junk1 = 0, *junk2 = 0;

    junk1 = (TH2F*)el_fr_v1_wwV0b->Clone("junk1");
    junk2 = (TH2F*)el_fr_v1_wwV0b->Clone("junk2");
    junk1->Reset();
    junk2->Reset();

    junk1->Add(ones,el_fr_v1_wwV0b,1,-1);  // = 1-FR
    junk2->Divide(el_fr_v1_wwV0b,junk1);   // = FR/(1-FR)
    pred_el_v1_wwV0b->Multiply(el_v1_wwV0b_NotNum,junk2);

    cout << "el_fr_v1_wwV0b:\n";
    cout << "observed  = " << el_num_wwV0b->Integral() << endl;
    cout << "predicted = " << pred_el_v1_wwV0b->Integral() << endl << endl;

    junk1 = (TH2F*)el_fr_v2_wwV0b->Clone("junk1");
    junk2 = (TH2F*)el_fr_v2_wwV0b->Clone("junk2");
    junk1->Reset();
    junk2->Reset();

    junk1->Add(ones,el_fr_v2_wwV0b,1,-1);  // = 1-FR
    junk2->Divide(el_fr_v2_wwV0b,junk1);   // = FR/(1-FR)
    pred_el_v2_wwV0b->Multiply(el_v2_wwV0b_NotNum,junk2);

    cout << "el_fr_v2_wwV0b:\n";
    cout << "observed  = " << el_num_wwV0b->Integral() << endl;
    cout << "predicted = " << pred_el_v2_wwV0b->Integral() << endl << endl;

    junk1 = (TH2F*)el_fr_v3_wwV0b->Clone("junk1");
    junk2 = (TH2F*)el_fr_v3_wwV0b->Clone("junk2");
    junk1->Reset();
    junk2->Reset();

    junk1->Add(ones,el_fr_v3_wwV0b,1,-1);  // = 1-FR
    junk2->Divide(el_fr_v3_wwV0b,junk1);   // = FR/(1-FR)
    pred_el_v3_wwV0b->Multiply(el_v3_wwV0b_NotNum,junk2);

    cout << "el_fr_v3_wwV0b:\n";
    cout << "observed  = " << el_num_wwV0b->Integral() << endl;
    cout << "predicted = " << pred_el_v3_wwV0b->Integral() << endl << endl;

    junk1 = (TH2F*)el_fr_v4_wwV0b->Clone("junk1");
    junk2 = (TH2F*)el_fr_v4_wwV0b->Clone("junk2");
    junk1->Reset();
    junk2->Reset();

    junk1->Add(ones,el_fr_v4_wwV0b,1,-1);  // = 1-FR
    junk2->Divide(el_fr_v4_wwV0b,junk1);   // = FR/(1-FR)
    pred_el_v4_wwV0b->Multiply(el_v4_wwV0b_NotNum,junk2);

    cout << "el_fr_v4_wwV0b:\n";
    cout << "observed  = " << el_num_wwV0b->Integral() << endl;
    cout << "predicted = " << pred_el_v4_wwV0b->Integral() << endl << endl;

    junk1 = (TH2F*)mu_fr_fo_wwV0_10->Clone("junk1");
    junk2 = (TH2F*)mu_fr_fo_wwV0_10->Clone("junk2");
    junk1->Reset();
    junk2->Reset();

    junk1->Add(ones,mu_fr_fo_wwV0_10,1,-1);  // = 1-FR
    junk2->Divide(mu_fr_fo_wwV0_10,junk1);   // = FR/(1-FR)
    pred_mu_fo_wwV0_10->Multiply(mu_fo_wwV0_10_NotNum,junk2);

    cout << "mu_fr_fo_wwV0_10:\n";
    cout << "observed  = " << mu_num_wwV0->Integral() << endl;
    cout << "predicted = " << pred_mu_fo_wwV0_10->Integral() << endl << endl;
}
