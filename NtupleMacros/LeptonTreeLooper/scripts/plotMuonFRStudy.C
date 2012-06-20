
void plotMuonFRStudy(TString det, TString pt)
{

    //
    // load lins
    //

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gROOT->ProcessLine(".L ../libLeptonTreeLooper.so");

    gROOT->ProcessLine(".L ~/tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    gStyle->SetOptTitle(0);
    gROOT->ForceStyle();

    //
    // open file
    //

    TString inputFile = "/smurf/dlevans/LeptonTree/V00-02-00_frtests2/DoubleMuRun2012APromptV1/merged_Cert_190456-193336_8TeV_PromptReco_Collisions12_JSON.root";
    TFile *f = new TFile(inputFile, "READ");
    TTree *leptons = (TTree*)f->Get("leptons");
    gROOT->cd();

    TCut kinematics("");
    TCut numerator("(leptonSelection&((1<<18)|(1<<19)))==((1<<18)|(1<<19))");
    TCut denominator("(eventSelection&(1<<5))==(1<<5) && HLT_Mu8_probe > 0");
    if (det == "EE" && pt == "PT1020") {
        kinematics = TCut("abs(probe.Eta()) > 1.479 && probe.Pt() < 20.0");
    }
    if (det == "EE" && pt == "PT20UP") {
        kinematics = TCut("abs(probe.Eta()) > 1.479 && probe.Pt() > 20.0");
    }
    if (det == "EB" && pt == "PT1020") {
        kinematics = TCut("abs(probe.Eta()) < 1.479 && probe.Pt() < 20.0");
    }
    if (det == "EB" && pt == "PT20UP") {
        kinematics = TCut("abs(probe.Eta()) < 1.479 && probe.Pt() > 20.0");
    }

    TH1F *h1_pfiso = new TH1F("h1_pfiso", "h1_pfiso; PFIso/p_{T}", 100, 0.0, 2.0);
    TH1F *h1_mva = new TH1F("h1_mva", "h1_pfiso; MVA", 100, -1.0, 1.0);

    leptons->Draw("TMath::Min(iso2011, 1.99) >> h1_pfiso", kinematics + denominator);
    leptons->Draw("TMath::Min(muonHZZ2012IsoRingsMVA, 0.99) >> h1_mva", kinematics + denominator);

    TCanvas *c1_pfiso = new TCanvas();
    c1_pfiso->SetLogy();
    c1_pfiso->cd();
    h1_pfiso->Draw();
    c1_pfiso->SaveAs("muonFRStudy_pfiso_"+det+"_"+pt+".png");

    TCanvas *c1_mva = new TCanvas();
    c1_mva->SetLogy();
    c1_mva->cd();
    h1_mva->Draw();
    c1_mva->SaveAs("muonFRStudy_mva_"+det+"_"+pt+".png");

    float n_denominator = leptons->GetEntries(kinematics + denominator + "iso2011<0.4");
    float n_numerator = leptons->GetEntries(numerator + kinematics + denominator + "iso2011<0.4");
    printf("%s %s\t Fake Rate = %4.3f\n", det.Data(), pt.Data(), n_numerator/n_denominator);

    float n_denominator = leptons->GetEntries(kinematics + denominator);
    float n_numerator = leptons->GetEntries(numerator + kinematics + denominator);
    printf("%s %s\t Fake Rate (no iso) = %4.3f\n", det.Data(), pt.Data(), n_numerator/n_denominator);

    n_denominator = leptons->GetEntries(kinematics + denominator + "muonHZZ2012IsoRingsMVA > -0.6");
    n_numerator = leptons->GetEntries(numerator + kinematics + denominator + "muonHZZ2012IsoRingsMVA > -0.6");
    printf("%s %s\t Fake Rate (new) = %4.3f\n", det.Data(), pt.Data(), n_numerator/n_denominator);

    f->Close();

}

