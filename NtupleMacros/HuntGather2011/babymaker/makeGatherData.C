
void makeGatherBaby(const char *inputFileName, const char *outputFileName)
{

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gSystem->Load("libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");
    gSystem->Load("libCMS2NtupleMacrosMT2.so");
    gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

    TChain *chain_data = new TChain("Events");
    chain_data->Add("/tas/cms2/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/skimmed_ntuple*.root");

    dilepbabymaker *dilepbaby = new dilepbabymaker();
    dilepbaby->ScanChain(inputFileName, outputFileName);

    delete chain_data;
    delete dilepbabymaker;

}


