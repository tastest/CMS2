void makeWWFakeRates(bool doels = true, bool domus = true) {
    gROOT->LoadMacro("eff2.C");
    gStyle->SetOptStat(0);

    //--------------------------------------------------------------------------
    // Load lepton data that you want to use to
    // make your fake rate. It should be a baby
    // ntuple.  
    //--------------------------------------------------------------------------

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

    //--------------------------------------------------------------------------
    // Define numerator and denominator cuts
    //--------------------------------------------------------------------------

    // A cut against Ws
    TCut notWCut  = "tcmet<20 && mt<25"; 

    // A pt cut...
    // Remember, we use 35 (at least for muons) to
    // minimize the impact of Ws
    TCut ptCut    = "pt>10 && pt<35";

    // Only consider events with >0 pfjet with pt
    // above a threshold separated by at least dR
    TCut jetCut   = "ptpfj1>15";

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

    //--------------------------------------------------------------------------
    // Define pt and eta bins of fake rate histograms
    //--------------------------------------------------------------------------

    double ybin[6]={10.,15.,20.,25.,30.,35.};
    int nbinsy = 5;
    //double ybin[4]={20.,25.,30.,35.};
    //int nbinsy = 3;
    double xbin[5]={0.0,1.0,1.479,2.0,2.5};
    int nbinsx = 4;

    //--------------------------------------------------------------------------
    // Book numerator and denominator histograms
    //--------------------------------------------------------------------------

    TH2F* el_num_wwV0   = new TH2F("el_num_wwV0","el_num_wwV0",nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v1_wwV0    = new TH2F("el_v1_wwV0","el_v1_wwV0",nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v2_wwV0    = new TH2F("el_v2_wwV0","el_v2_wwV0",nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v3_wwV0    = new TH2F("el_v3_wwV0","el_v3_wwV0",nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v4_wwV0    = new TH2F("el_v4_wwV0","el_v4_wwV0",nbinsx,xbin,nbinsy,ybin);

    TH2F* el_num_wwV0b  = new TH2F("el_num_wwV0b","el_num_wwV0b",nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v1_wwV0b   = new TH2F("el_v1_wwV0b","el_v1_wwV0b",nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v2_wwV0b   = new TH2F("el_v2_wwV0b","el_v2_wwV0b",nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v3_wwV0b   = new TH2F("el_v3_wwV0b","el_v3_wwV0b",nbinsx,xbin,nbinsy,ybin);
    TH2F* el_v4_wwV0b   = new TH2F("el_v4_wwV0b","el_v4_wwV0b",nbinsx,xbin,nbinsy,ybin);

    TH2F* mu_num_wwV0   = new TH2F("mu_num_wwV0","mu_num_wwV0",nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_fo_wwV0_04 = new TH2F("mu_fo_wwV0_04","mu_fo_wwV0_04",nbinsx,xbin,nbinsy,ybin);
    TH2F* mu_fo_wwV0_10 = new TH2F("mu_fo_wwV0_10","mu_fo_wwV0_10",nbinsx,xbin,nbinsy,ybin);

    //--------------------------------------------------------------------------
    // Fill Histograms
    //--------------------------------------------------------------------------

    if (doels) {
        ch_el->Draw("pt:abs(eta)>>el_num_wwV0",is_el_num_wwV0);
        ch_el->Draw("pt:abs(eta)>>el_v1_wwV0",is_el_v1_wwV0);
        ch_el->Draw("pt:abs(eta)>>el_v2_wwV0",is_el_v2_wwV0);
        ch_el->Draw("pt:abs(eta)>>el_v3_wwV0",is_el_v3_wwV0);
        ch_el->Draw("pt:abs(eta)>>el_v4_wwV0",is_el_v4_wwV0);

        ch_el->Draw("pt:abs(eta)>>el_num_wwV0b",is_el_num_wwV0b);
        ch_el->Draw("pt:abs(eta)>>el_v1_wwV0b",is_el_v1_wwV0b);
        ch_el->Draw("pt:abs(eta)>>el_v2_wwV0b",is_el_v2_wwV0b);
        ch_el->Draw("pt:abs(eta)>>el_v3_wwV0b",is_el_v3_wwV0b);
        ch_el->Draw("pt:abs(eta)>>el_v4_wwV0b",is_el_v4_wwV0b);
    }

    if (domus) {
        ch_mu->Draw("pt:abs(eta)>>mu_num_wwV0",is_mu_num_wwV0);
        ch_mu->Draw("pt:abs(eta)>>mu_fo_wwV0_04",is_mu_fo_wwV0_04);
        ch_mu->Draw("pt:abs(eta)>>mu_fo_wwV0_10",is_mu_fo_wwV0_10);
    }

    //--------------------------------------------------------------------------
    // Make fake rate histograms. Last argument is
    // the histogram name
    //--------------------------------------------------------------------------

    TH2F* el_fr_v1_wwV0    = doels ? eff2(el_v1_wwV0,el_num_wwV0,"el_fr_v1_wwV0") : 0;
    TH2F* el_fr_v2_wwV0    = doels ? eff2(el_v2_wwV0,el_num_wwV0,"el_fr_v2_wwV0") : 0;
    TH2F* el_fr_v3_wwV0    = doels ? eff2(el_v3_wwV0,el_num_wwV0,"el_fr_v3_wwV0") : 0;
    TH2F* el_fr_v4_wwV0    = doels ? eff2(el_v4_wwV0,el_num_wwV0,"el_fr_v4_wwV0") : 0;

    TH2F* el_fr_v1_wwV0b   = doels ? eff2(el_v1_wwV0b,el_num_wwV0b,"el_fr_v1_wwV0b") : 0;
    TH2F* el_fr_v2_wwV0b   = doels ? eff2(el_v2_wwV0b,el_num_wwV0b,"el_fr_v2_wwV0b") : 0;
    TH2F* el_fr_v3_wwV0b   = doels ? eff2(el_v3_wwV0b,el_num_wwV0b,"el_fr_v3_wwV0b") : 0;
    TH2F* el_fr_v4_wwV0b   = doels ? eff2(el_v4_wwV0b,el_num_wwV0b,"el_fr_v4_wwV0b") : 0;

    TH2F* mu_fr_fo_wwV0_04 = domus ? eff2(mu_fo_wwV0_04,mu_num_wwV0,"mu_fr_fo_wwV0_04") : 0;
    TH2F* mu_fr_fo_wwV0_10 = domus ? eff2(mu_fo_wwV0_10,mu_num_wwV0,"mu_fr_fo_wwV0_10") : 0;

    if (doels) {
        saveHist("ww_el_fr_EG.root", "el_*");
        //saveHist("ww_el_fr_qcd30.root", "el_*");
    }

    if (domus) {
        saveHist("ww_mu_fr_Mu.root", "mu_*");
        //saveHist("ww_mu_fr_inclMu.root", "mu_*");
    }
}
