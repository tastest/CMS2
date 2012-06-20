
void plotElectronFRStudy(TString det, TString pt)
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

    TString inputFile = "/smurf/dlevans/LeptonTree/V00-02-00_frtests2/DoubleElectronRun2012APromptV1/merged_Cert_190456-193336_8TeV_PromptReco_Collisions12_JSON.root";
    TFile *f = new TFile(inputFile, "READ");
    TTree *leptons = (TTree*)f->Get("leptons");
    gROOT->cd();

    TCut kinematics("");
    TCut detIso  = ("hcaliso/probe.Pt() < 0.2 && ecaliso/probe.Pt() < 0.2 && trkiso/probe.Pt() < 0.2");
    TCut pfIso04Denom("(pfchiso04 + TMath::Max(0.0, pfnhiso04 + pfemiso04 - ea04 * TMath::Max(0.0, rhoIsoAll)))/probe.Pt() < 0.40");
    TCut pfIso04Numer("(pfchiso04 + TMath::Max(0.0, pfnhiso04 + pfemiso04 - ea04 * TMath::Max(0.0, rhoIsoAll)))/probe.Pt() < 0.15");
    TCut denominator("(eventSelection&(1<<4))==(1<<4) && HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe > 0");
    TCut numerator("(leptonSelection&(1<<6))==(1<<6)"+pfIso04Numer);
    if (det == "EE" && pt == "PT1020") {
        kinematics = TCut("abs(sceta) > 1.479 && probe.Pt() < 20.0");
    }
    if (det == "EE" && pt == "PT20UP") {
        kinematics = TCut("abs(sceta) > 1.479 && probe.Pt() > 20.0");
    }
    if (det == "EB" && pt == "PT1020") {
        kinematics = TCut("abs(sceta) < 1.479 && probe.Pt() < 20.0");
    }
    if (det == "EB" && pt == "PT20UP") {
        kinematics = TCut("abs(sceta) < 1.479 && probe.Pt() > 20.0");
    }

    float n_denominator = leptons->GetEntries(kinematics + denominator + detIso);
    float n_numerator = leptons->GetEntries(numerator + kinematics + denominator + detIso);
    printf("%4.3f, %4.3f\n", n_denominator, n_numerator);
    printf("%s %s\t Fake Rate (DetIso Denominator) = %4.3f +/- %4.3f\n", det.Data(), pt.Data(), n_numerator/n_denominator, (n_numerator/n_denominator) * sqrt(1/n_numerator + 1/n_denominator));

    float n_denominator = leptons->GetEntries(kinematics + denominator + pfIso04Denom);
    float n_numerator = leptons->GetEntries(numerator + kinematics + denominator + pfIso04Denom);
    printf("%4.3f, %4.3f\n", n_denominator, n_numerator);
    printf("%s %s\t Fake Rate (PFIso Denominator) = %4.3f +/- %4.3f\n", det.Data(), pt.Data(), n_numerator/n_denominator, (n_numerator/n_denominator) * sqrt(1/n_numerator + 1/n_denominator));

    f->Close();

}

