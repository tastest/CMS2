
void doMC() {

    //
    // the looper
    //
    gSystem->SetIncludePath(Form("%s -I../../Tools", gSystem->GetIncludePath()));

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gSystem->Load("libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");
    gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");

    // for likelihood id
    gSystem->Load("../../Tools/EgammaAnalysisTools/lib/libElectronLikelihoodId.so");

    gROOT->ProcessLine(".L histtools.C+");

    //
    // output file for histograms
    //

    MyScanChain *looper = new MyScanChain();

    //
    // chains for input files
    //TString ntuple_location = "/nfs-3/userdata/cms2/";
    TString ntuple_location = "/tas/cms2/";

    TChain *chain_hww130 = new TChain("Events");
    chain_hww130->Add(ntuple_location + "GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall10-START38_V12-v1/V03-06-18/merged*.root");

    //
    // run it
    //
    // k factor for pythia samples
    looper->ScanChain(false, "hww130", chain_hww130);

    //
    // tidy up
    //
    delete looper;
    delete chain_hww130;

}

