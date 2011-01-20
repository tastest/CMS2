
void makeGatherBaby(const char *inputFileName, const char *outputFileName)
{

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gSystem->Load("libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");
    gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

    emubabymaker *emubaby = new emubabymaker();
    std::cout << "about to do scanchain" << std::endl;
    emubaby->ScanChain(inputFileName, outputFileName);
    delete emubaby;

}


