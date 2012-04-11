{

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gROOT->ProcessLine(".L ../../../Smurf/Core/LeptonTree.h+");
    gROOT->ProcessLine(".L libLeptonTreeLooper.so");

    LeptonTreeLooper *looper = new LeptonTreeLooper();

    //
    // data
    //

    TChain *ch_data_2011B = new TChain("leptons");
    ch_data_2011B->Add("/smurf/dlevans/LeptonTree/Wednesday110412_00/DoubleElectronRun2011BPromptV1/merged.root");
    looper->setGoodRunList("runlists/HWW_2011.jmu");
    looper->loop(ch_data_2011B, "data_2011B");

    TChain *ch_data_2012A = new TChain("leptons");
    ch_data_2012A->Add("/smurf/dlevans/LeptonTree/Wednesday110412_00/DoubleElectronRun2012APromptV1/merged.root");
    looper->setGoodRunList("runlists/DCSOnly_2012A.jmu");
    looper->loop(ch_data_2012A, "data_2012A");

    //
    // mc
    //

    //TChain *ch_dyee = new TChain("leptons");
    //ch_dyee->Add("/smurf/dlevans/LeptonTree/V00-00-08/dyee.root");
    //looper->loop(ch_dyee, "dyee");


    const std::string outFile = Form("histos.root");
    saveHist(outFile.c_str());  
    deleteHistos();

    delete looper;
    delete ch_data_2011B;
    delete ch_data_2012A;
    //delete ch_dyee;

}


