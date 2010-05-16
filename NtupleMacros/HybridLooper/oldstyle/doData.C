
void doData() {

    //
    // the looper
    //
    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gSystem->Load("libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");
    gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");
    gROOT->ProcessLine(".L histtools.C+");

    //
    // output file for histograms
    //

    MyScanChain *looper = new MyScanChain();

    // DATA
    TChain *chain_whunt_skim = new TChain("Events");
    //chain_whunt_skim->Add("/tas03/disk01/whunt/skim/emuskim_*.root");
    // on dle laptop/other
    //chain_whunt_skim->Add("/Users/dlevans/tas03/disk01/whunt/skim/emuskim_merged.root");
    chain_whunt_skim->Add("/tmp/emuskim_merged_1400_160510.root");
    looper->ScanChain(true, "whunt", chain_whunt_skim);

    //
    // write histograms
    // 
    const char* outFile = "histos_data.root";
    hist::saveHist(outFile); 
    hist::deleteHistos();

    //
    // tidy up
    //
    delete looper;
    delete chain_whunt_skim;

}

