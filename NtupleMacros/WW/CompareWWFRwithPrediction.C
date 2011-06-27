void CompareWWFRwithPrediction(bool doels = true, bool domus = true);

// Method by pointer
TH2F* eff2(TH2F* h1, TH2F* h2, const char* name="eff"){

  // first, verify that all histograms have same binning
  // nx is the number of visible bins
  // nxtot = nx+2 includes underflow and overflow
  Int_t nx = h1->GetNbinsX();
  if (h2->GetNbinsX() != nx) {
    cout << "Histograms must have same number of bins" << endl;
    return 0;
  }

  // get the new histogram
  TH2F* temp = (TH2F*) h1->Clone(name);
  temp->SetTitle(name);
  temp->Reset();
  temp->Sumw2();

  // Do the calculation
  temp->Divide(h2,h1,1.,1.,"B");

  // Done
  return temp;
}


// Method by name
TH2F* eff2(const char* name1, const char* name2, const char* name="eff"){

  // Get a list of object and their iterator
  TList* list = gDirectory->GetList() ;
  TIterator* iter = list->MakeIterator();

  // Loop over objects, set the pointers
  TObject* obj;
  TH2F* h1=0;
  TH2F* h2=0;
  TString str1 = Form("%s",name1);
  TString str2 = Form("%s",name2);
  while((obj=iter->Next())) {
    TString objName = obj->GetName();
    if (objName == str1) h1 = (TH2F*) obj;
    if (objName == str2) h2 = (TH2F*) obj;
  }

  // quit if not found
  if (h1 == 0) {
    cout << "Histogram " << name1 << " not found" << endl;
    return 0;
  }
  if (h2 == 0) {
    cout << "Histogram " << name2 << " not found" << endl;
    return 0;
  }

  // Call the method by pointer
  TH2F* temp = eff2(h1, h2, name);
  return temp;
}

void CompareWWFRwithPrediction(bool doels, bool domus) {
    gStyle->SetOptStat(0);

    //-------------------------------------------------------------------------
    // Load the baby ntuple
    //-------------------------------------------------------------------------

    TChain *ch_el = new TChain("tree");
    TChain *ch_egmon = new TChain("tree");
    if (doels) {
        ch_el   ->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/EG.root");
        ch_egmon->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/EGMon.root");
        //ch_el   ->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/qcd_pt_30to50_fall10.root");
    }

    TChain *ch_mu = new TChain("tree");
    if (domus) {
        ch_mu->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/Mu.root");
        //ch_mu->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/mu15.root");
        //ch_mu->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/mu15_1.root");
        //ch_mu->Add("/tas/cms2/FRBabies/FakeRates10November2010-v2/mu10.root");
    }

    //-------------------------------------------------------------------------
    // Load the fake rate histograms
    // Need to know the name of the 
    // fake rate files and histograms
    //-------------------------------------------------------------------------

    loadHist("ww_el_fr_EGandEGMon.root",0,0,"*fr*");
    //loadHist("ww_el_fr_qcd_pt_30to50_fall10.root",0,0,"*fr*");
    loadHist("ww_mu_fr_Mu.root",0,0,"*fr*");
    //loadHist("ww_mu_fr_mu15.root",0,0,"*fr*");
    //loadHist("ww_mu_fr_mu10.root",0,0,"*fr*");

    TH2F* el_fr_v1_wwV1   = (TH2F*)gDirectory->Get("el_fr_v1_wwV1");
    TH2F* el_fr_v2_wwV1   = (TH2F*)gDirectory->Get("el_fr_v2_wwV1");
    TH2F* el_fr_v3_wwV1   = (TH2F*)gDirectory->Get("el_fr_v3_wwV1");
    TH2F* el_fr_v4_wwV1   = (TH2F*)gDirectory->Get("el_fr_v4_wwV1");
    TH2F* mu_fr_fo_wwV1_10 = (TH2F*)gDirectory->Get("mu_fr_fo_wwV1_10");
    TH2F* mu_fr_fo_wwV1_10_d0 = (TH2F*)gDirectory->Get("mu_fr_fo_wwV1_10_d0");

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

    //
    // The trigger selections
    //

    // For EG
    TCut trgCutEl_eg = "((el10_lw>1 || el10_sw>1 || el15_lw>1) && run < 141956) ||\
                        ((el15_sw>1 || el20_sw>1) && run <= 145761)";
    // For EGMon
    TCut trgCutEl_egmon = "((el10_sw>1 || el17_sw>1) && run < 149181) ||\
                           (el10_sw_v2>1 || el17_sw_v2>1)";
    // For Mu
    TCut trgCutMu = "mu9>1 || mu11>1 || mu15>1";

    // Numerator selections
    TCut is_el_num_wwV1  = "num_wwV1&&abs(id)==11"+jetCut+ptCut+notWCut;
    TCut is_el_num_wwV1 = "num_wwV1&&abs(id)==11"+jetCut+ptCut+notWCut;

    TCut is_mu_num_wwV1  = "num_wwV1&&abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;

    // Denominator selections
    TCut is_el_v1_wwV1    = "v1_wwV1&&abs(id)==11"+jetCut+ptCut+notWCut;
    TCut is_el_v2_wwV1    = "v2_wwV1&&abs(id)==11"+jetCut+ptCut+notWCut;
    TCut is_el_v3_wwV1    = "v3_wwV1&&abs(id)==11"+jetCut+ptCut+notWCut;
    TCut is_el_v4_wwV1    = "v4_wwV1&&abs(id)==11"+jetCut+ptCut+notWCut;

    TCut is_mu_fo_wwV1_04 = "fo_wwV1_04&&abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;
    TCut is_mu_fo_wwV1_10 = "fo_wwV1_10&&abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;
    TCut is_mu_fo_wwV1_10_d0 = "fo_wwV1_10_d0&&abs(id)==13"+trgCutMu+jetCut+ptCut+notWCut;

    //-------------------------------------------------------------------------
    // Make three empty histograms with the same
    // binning as the fake rate histogram
    //-------------------------------------------------------------------------

    TH2F* el_num_wwV1 = el_fr_v1_wwV1->Clone("el_num_wwV1");
    el_num_wwV1->Reset();

    TH2F* el_v1_wwV1_NotNum = el_fr_v1_wwV1->Clone("el_v1_wwV1_NotNum");
    el_v1_wwV1_NotNum->Reset();
    TH2F* el_v2_wwV1_NotNum = el_fr_v2_wwV1->Clone("el_v2_wwV1_NotNum");
    el_v2_wwV1_NotNum->Reset();
    TH2F* el_v3_wwV1_NotNum = el_fr_v3_wwV1->Clone("el_v3_wwV1_NotNum");
    el_v3_wwV1_NotNum->Reset();
    TH2F* el_v4_wwV1_NotNum = el_fr_v4_wwV1->Clone("el_v4_wwV1_NotNum");
    el_v4_wwV1_NotNum->Reset();

    TH2F* pred_el_v1_wwV1 = el_fr_v1_wwV1->Clone("pred_el_v1_wwV1");
    pred_el_v1_wwV1->Reset();
    TH2F* pred_el_v2_wwV1 = el_fr_v2_wwV1->Clone("pred_el_v2_wwV1");
    pred_el_v2_wwV1->Reset();
    TH2F* pred_el_v3_wwV1 = el_fr_v3_wwV1->Clone("pred_el_v3_wwV1");
    pred_el_v3_wwV1->Reset();
    TH2F* pred_el_v4_wwV1 = el_fr_v4_wwV1->Clone("pred_el_v4_wwV1");
    pred_el_v4_wwV1->Reset();

    TH2F* mu_num_wwV1 = mu_fr_fo_wwV1_10->Clone("mu_num_wwV1");
    mu_num_wwV1->Reset();

    TH2F* mu_fo_wwV1_10_NotNum = mu_fr_fo_wwV1_10->Clone("mu_fo_wwV1_10_NotNum");
    mu_fo_wwV1_10_NotNum->Reset();
    TH2F* mu_fo_wwV1_10_d0_NotNum = mu_fr_fo_wwV1_10_d0->Clone("mu_fo_wwV1_10_d0_NotNum");
    mu_fo_wwV1_10_d0_NotNum->Reset();

    TH2F* pred_mu_fo_wwV1_10 = mu_fr_fo_wwV1_10->Clone("pred_mu_fo_wwV1_10");
    pred_mu_fo_wwV1_10->Reset();
    TH2F* pred_mu_fo_wwV1_10_d0 = mu_fr_fo_wwV1_10_d0->Clone("pred_mu_fo_wwV1_10_d0");
    pred_mu_fo_wwV1_10_d0->Reset();

    //-------------------------------------------------------------------------
    // Fill two of the three histograms.
    // The 1st one is filled with the numerator
    // The 2nd one is filled with the denominator but not numerator
    // The 3rd one will be filled with the prediction
    //-------------------------------------------------------------------------

    float ptmax = el_num_wwV1->GetYaxis()->GetXmax();
    cout << "Leptons with pt > " << ptmax   << " will be put in the last bin" << endl; 

    ch_el->Draw(Form("min(pt,%f-0.1):abs(eta)>>el_num_wwV1",ptmax),trgCutEl_eg&&is_el_num_wwV1);
    ch_el->Draw(Form("min(pt,%f-0.1):abs(eta)>>el_v1_wwV1_NotNum",ptmax),trgCutEl_eg&&is_el_v1_wwV1&&!is_el_num_wwV1);
    ch_el->Draw(Form("min(pt,%f-0.1):abs(eta)>>el_v2_wwV1_NotNum",ptmax),trgCutEl_eg&&is_el_v2_wwV1&&!is_el_num_wwV1);
    ch_el->Draw(Form("min(pt,%f-0.1):abs(eta)>>el_v3_wwV1_NotNum",ptmax),trgCutEl_eg&&is_el_v3_wwV1&&!is_el_num_wwV1);
    ch_el->Draw(Form("min(pt,%f-0.1):abs(eta)>>el_v4_wwV1_NotNum",ptmax),trgCutEl_eg&&is_el_v4_wwV1&&!is_el_num_wwV1);
    ch_egmon->Draw(Form("min(pt,%f-0.1):abs(eta)>>+el_num_wwV1",ptmax),(trgCutEl_egmon&&is_el_num_wwV1));
    ch_egmon->Draw(Form("min(pt,%f-0.1):abs(eta)>>+el_v1_wwV1_NotNum",ptmax),trgCutEl_egmon&&is_el_v1_wwV1&&!is_el_num_wwV1);
    ch_egmon->Draw(Form("min(pt,%f-0.1):abs(eta)>>+el_v2_wwV1_NotNum",ptmax),trgCutEl_egmon&&is_el_v2_wwV1&&!is_el_num_wwV1);
    ch_egmon->Draw(Form("min(pt,%f-0.1):abs(eta)>>+el_v3_wwV1_NotNum",ptmax),trgCutEl_egmon&&is_el_v3_wwV1&&!is_el_num_wwV1);
    ch_egmon->Draw(Form("min(pt,%f-0.1):abs(eta)>>+el_v4_wwV1_NotNum",ptmax),trgCutEl_egmon&&is_el_v4_wwV1&&!is_el_num_wwV1);

    ch_mu->Draw(Form("min(pt,%f-0.1):abs(eta)>>mu_num_wwV1",ptmax),is_mu_num_wwV1);
    ch_mu->Draw(Form("min(pt,%f-0.1):abs(eta)>>mu_fo_wwV1_10_NotNum",ptmax),is_mu_fo_wwV1_10&&!is_mu_num_wwV1);
    ch_mu->Draw(Form("min(pt,%f-0.1):abs(eta)>>mu_fo_wwV1_10_d0_NotNum",ptmax),is_mu_fo_wwV1_10_d0&&!is_mu_num_wwV1);

    //-------------------------------------------------------------------------
    // Now we need to get a histogram full of ones
    //-------------------------------------------------------------------------

    TH2F* ones = el_fr_v1_wwV1->Clone();
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

    junk1 = (TH2F*)el_fr_v1_wwV1->Clone("junk1");
    junk2 = (TH2F*)el_fr_v1_wwV1->Clone("junk2");
    junk1->Reset();
    junk2->Reset();

    junk1->Add(ones,el_fr_v1_wwV1,1,-1);  // = 1-FR
    junk2->Divide(el_fr_v1_wwV1,junk1);   // = FR/(1-FR)
    pred_el_v1_wwV1->Multiply(el_v1_wwV1_NotNum,junk2);

    cout << "el_fr_v1_wwV1:\n";
    cout << "observed  = " << el_num_wwV1->Integral() << endl;
    cout << "predicted = " << pred_el_v1_wwV1->Integral() << endl << endl;

    junk1 = (TH2F*)el_fr_v2_wwV1->Clone("junk1");
    junk2 = (TH2F*)el_fr_v2_wwV1->Clone("junk2");
    junk1->Reset();
    junk2->Reset();

    junk1->Add(ones,el_fr_v2_wwV1,1,-1);  // = 1-FR
    junk2->Divide(el_fr_v2_wwV1,junk1);   // = FR/(1-FR)
    pred_el_v2_wwV1->Multiply(el_v2_wwV1_NotNum,junk2);

    cout << "el_fr_v2_wwV1:\n";
    cout << "observed  = " << el_num_wwV1->Integral() << endl;
    cout << "predicted = " << pred_el_v2_wwV1->Integral() << endl << endl;

    junk1 = (TH2F*)el_fr_v3_wwV1->Clone("junk1");
    junk2 = (TH2F*)el_fr_v3_wwV1->Clone("junk2");
    junk1->Reset();
    junk2->Reset();

    junk1->Add(ones,el_fr_v3_wwV1,1,-1);  // = 1-FR
    junk2->Divide(el_fr_v3_wwV1,junk1);   // = FR/(1-FR)
    pred_el_v3_wwV1->Multiply(el_v3_wwV1_NotNum,junk2);

    cout << "el_fr_v3_wwV1:\n";
    cout << "observed  = " << el_num_wwV1->Integral() << endl;
    cout << "predicted = " << pred_el_v3_wwV1->Integral() << endl << endl;

    junk1 = (TH2F*)el_fr_v4_wwV1->Clone("junk1");
    junk2 = (TH2F*)el_fr_v4_wwV1->Clone("junk2");
    junk1->Reset();
    junk2->Reset();

    junk1->Add(ones,el_fr_v4_wwV1,1,-1);  // = 1-FR
    junk2->Divide(el_fr_v4_wwV1,junk1);   // = FR/(1-FR)
    pred_el_v4_wwV1->Multiply(el_v4_wwV1_NotNum,junk2);

    cout << "el_fr_v4_wwV1:\n";
    cout << "observed  = " << el_num_wwV1->Integral() << endl;
    cout << "predicted = " << pred_el_v4_wwV1->Integral() << endl << endl;

    junk1 = (TH2F*)mu_fr_fo_wwV1_10->Clone("junk1");
    junk2 = (TH2F*)mu_fr_fo_wwV1_10->Clone("junk2");
    junk1->Reset();
    junk2->Reset();

    junk1->Add(ones,mu_fr_fo_wwV1_10,1,-1);  // = 1-FR
    junk2->Divide(mu_fr_fo_wwV1_10,junk1);   // = FR/(1-FR)
    pred_mu_fo_wwV1_10->Multiply(mu_fo_wwV1_10_NotNum,junk2);

    cout << "mu_fr_fo_wwV1_10:\n";
    cout << "observed  = " << mu_num_wwV1->Integral() << endl;
    cout << "predicted = " << pred_mu_fo_wwV1_10->Integral() << endl << endl;

    junk1 = (TH2F*)mu_fr_fo_wwV1_10_d0->Clone("junk1");
    junk2 = (TH2F*)mu_fr_fo_wwV1_10_d0->Clone("junk2");
    junk1->Reset();
    junk2->Reset();

    junk1->Add(ones,mu_fr_fo_wwV1_10_d0,1,-1);  // = 1-FR
    junk2->Divide(mu_fr_fo_wwV1_10_d0,junk1);   // = FR/(1-FR)
    pred_mu_fo_wwV1_10_d0->Multiply(mu_fo_wwV1_10_d0_NotNum,junk2);

    cout << "mu_fr_fo_wwV1_10_d0:\n";
    cout << "observed  = " << mu_num_wwV1->Integral() << endl;
    cout << "predicted = " << pred_mu_fo_wwV1_10_d0->Integral() << endl << endl;
}
