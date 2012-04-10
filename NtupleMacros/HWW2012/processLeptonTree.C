
#include "SmurfDataTypes.h"

void test()
{
//    processLeptonTree("dyee_merged_ntuple.root", 45, "/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple.root", 0, "nogoodrunlist");

    processLeptonTree("test.root", SmurfDataType::data, 
        "/smurf/dlevans/CMSSW_5_2_3_patch3_V05-02-07/DoubleElectron_Run2012A-PromptReco-v1_AOD/merged/merged_ntuple_190459_0.root", 
        true, "");

}

void processLeptonTree(TString outfileid, enum SmurfDataType::DataType sample, TString file, bool realData, TString goodrunlist)
{

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gSystem->Load("libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");
    gSystem->Load("/tas/dlevans/HWW2012/CMSSW_5_2_3/src/CMS2/NtupleMacros/HWW2012/libMiniFWLite.so");

    //
    // const config parameters
    //

    const bool lockToCoreSelectors = false;
    const bool useLHeleId = false;
    const bool useMVAeleId = true;
    const bool doDYNNLOw = true;
    const unsigned int prescale = 1;
    const double integratedLumi = 1000.0; // pb^1  

    //
    // create looper
    //

    std::cout << "going to loop" << std::endl;
    LeptonTreeMaker *looper = new LeptonTreeMaker(lockToCoreSelectors, useLHeleId, useMVAeleId, doDYNNLOw, prescale, realData);
    looper->SetBatchMode(true);

    //
    // set up chain
    //

    TChain *chain = new TChain("Events");
    chain->Add(file);

    //
    // loop
    //

    looper->ScanChain(outfileid, chain, sample, integratedLumi, -1, -1, false, realData, goodrunlist);

    //
    // tidy up
    //

    delete looper;
    delete chain;

}



