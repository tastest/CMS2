
void doMC() {

    //
    // the looper
    //
    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gSystem->Load("/tas07/disk00/jribnik/lhapdf-5.8.3/lib/.libs/libLHAPDF.so");

    gSystem->Load("libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");

    gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
    gROOT->ProcessLine(".L histtools.C+");

    //
    // output file for histograms
    //

    MyScanChain *looper = new MyScanChain();

    //
    // specify the PDF set to use
    // this should be a valid LHgrid file
    //
    //looper->specifyPDF("/tas07/disk00/jribnik/lhapdf/share/lhapdf/PDFsets/cteq61");
    //looper->specifyPDF("pdfs/MSTW2008nlo90cl");
    looper->specifyPDF("pdfs/cteq61", "pdfs/cteq66");
    //looper->specifyPDF("pdfs/cteq61", "pdfs/MSTW2008nlo68cl");
    //looper->specifyPDF("pdfs/cteq61", "pdfs/NNPDF20_100");

    //
    // chains for input files
    //TString ntuple_location = "/nfs-3/userdata/cms2/";
    TString ntuple_location = "/tas/cms2/";

    // ttbar
    TChain *chain_ttbar = new TChain("Events");
    chain_ttbar->Add(ntuple_location + "TTbar_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*.root");

    TChain *chain_ttbarmg = new TChain("Events");
    chain_ttbarmg->Add(ntuple_location + "TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*.root");

    // zee
    TChain *chain_zee = new TChain("Events");
    chain_zee->Add(ntuple_location + "EarlyDataSamples/Zee_Spring10-START3X_V26_S09-v1/merged_ntuple*.root");

    // WW
    TChain *chain_ww = new TChain("Events");
    chain_ww->Add(ntuple_location + "WW_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple.root");

    // zee
    TChain *chain_zee = new TChain("Events");
    //chain_zee->Add(ntuple_location +"Zee_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*.root");
    chain_zee->Add("/tas/dlevans/data_skim/zee_skim_2020.root");
    chain_zee->Add(ntuple_location +"DYee_M10to20_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");
    chain_zee->Add(ntuple_location +"ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");

    // zmumu
    TChain *chain_zmumu = new TChain("Events");
    //chain_zmumu->Add(ntuple_location +"Zmumu_Spring10-START3X_V26_S09-v1/V03-04-13-07/*.root");
    chain_zmumu->Add("/tas/dlevans/data_skim/zmumu_skim_2020.root");
    chain_zmumu->Add(ntuple_location +"DYmumu_M10to20_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");
    chain_zmumu->Add(ntuple_location +"ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");


    // dy combined
    TChain *chain_dy = new TChain("Events");
    chain_dy->Add("/tas/dlevans/data_skim/zee_skim_2020.root");
    chain_dy->Add(ntuple_location +"DYee_M10to20_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");
    chain_dy->Add(ntuple_location +"ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");
    chain_dy->Add("/tas/dlevans/data_skim/zmumu_skim_2020.root");
    chain_dy->Add(ntuple_location +"DYmumu_M10to20_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");
    chain_dy->Add(ntuple_location +"ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");


    //
    // run it
    //

    looper->ScanChain("ttbar", chain_ttbarmg, (157.5/165.0));
    looper->ScanChain("dyee", chain_zee, (1666.0/1300.0));
    looper->ScanChain("dymm", chain_zmumu, (1666.0/1300.0));

    //
    // write histograms
    // 

    const char* outFile = "histos_mc.root";
    hist::saveHist(outFile); 
    hist::deleteHistos();

    //
    // tidy up
    //

    delete looper;
    delete chain_ttbar;
    delete chain_ttbarmg;
    delete chain_zee;
    delete chain_ww;
    delete chain_zmumu;
    delete chain_zee;
    delete chain_dy;

}

