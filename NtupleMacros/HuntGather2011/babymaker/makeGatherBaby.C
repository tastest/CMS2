
void makeGatherBaby(const char *inputFileName, const char *outputFileName)
{

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gSystem->Load("../libs/libCMS2NtupleMacrosCORE.so");
    gSystem->Load("../libs/libHuntGather2011Babymaker.so");
    gSystem->Load("../libs/libCMS2NtupleMacrosCOREMT2.so");
    gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");

    dilepbabymaker *dilepbaby = new dilepbabymaker();
    std::cout << "about to do scanchain" << std::endl;
    dilepbaby->ScanChain(inputFileName, outputFileName);
}


