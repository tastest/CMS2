
void doMC() {

    //
    // the looper
    //
    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gROOT->ProcessLine(".L SmurfTree.h+");

    gSystem->Load("/tas03/home/dlevans/LHAPDF/lib/libLHAPDF.so");
    //gSystem->Load("/afs/cern.ch/cms/slc5_amd64_gcc434/external/lhapdf/5.8.5/full/lib/libLHAPDF.so");

    gROOT->ProcessLine(".L histtools.C+");

    gSystem->Load("libCMS2NtupleMacrosLooper.so");


    //
    // output_V5 file for histograms
    //

    MyScanChain *looper = new MyScanChain();

    //
    // run it
    //

    std::vector<std::string> pdfSets;

    // CTEQ6
    pdfSets.push_back("cteq6mE");

    // CT10
    pdfSets.push_back("CT10");
    pdfSets.push_back("CT10as");

    // MSTW
    pdfSets.push_back("MSTW2008nlo68cl");
    pdfSets.push_back("MSTW2008nlo68cl_asmz+68cl");
    pdfSets.push_back("MSTW2008nlo68cl_asmz+68clhalf");
    pdfSets.push_back("MSTW2008nlo68cl_asmz-68cl");
    pdfSets.push_back("MSTW2008nlo68cl_asmz-68clhalf");

    // NNPDF
    pdfSets.push_back("NNPDF20_as_0116_100");
    pdfSets.push_back("NNPDF20_as_0117_100");
    pdfSets.push_back("NNPDF20_as_0118_100");
    pdfSets.push_back("NNPDF20_100");
    pdfSets.push_back("NNPDF20_as_0120_100");    
    pdfSets.push_back("NNPDF20_as_0121_100");
    pdfSets.push_back("NNPDF20_as_0122_100");

    for (unsigned int i = 0; i < pdfSets.size(); ++i) {
        looper->ScanChain("nn_hww115_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_115train_0jets_hww115.root");
        //looper->ScanChain("nn_hww120_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_120train_0jets_hww120.root");
        //looper->ScanChain("nn_hww130_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_130train_0jets_hww130.root");
        //looper->ScanChain("nn_hww140_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_140train_0jets_hww140.root");
        //looper->ScanChain("bdtg_hww150_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_150train_0jets_hww150.root");
        //looper->ScanChain("bdtg_hww160_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_160train_0jets_hww160.root");
        //looper->ScanChain("nn_hww170_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_170train_0jets_hww170.root");
        //looper->ScanChain("nn_hww180_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_180train_0jets_hww180.root");
        //looper->ScanChain("bdtg_hww190_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_190train_0jets_hww190.root");
        //looper->ScanChain("nn_hww200_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_200train_0jets_hww200.root");
        //looper->ScanChain("nn_hww250_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_250train_0jets_hww250.root");
        //looper->ScanChain("nn_hww300_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_300train_0jets_hww300.root");
        //looper->ScanChain("nn_hww350_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_350train_0jets_hww350.root");
        //looper->ScanChain("nn_hww400_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_400train_0jets_hww400.root");
        //looper->ScanChain("nn_hww450_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_450train_0jets_hww450.root");
        //looper->ScanChain("nn_hww500_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_500train_0jets_hww500.root");
        //looper->ScanChain("nn_hww550_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_550train_0jets_hww550.root");
        //looper->ScanChain("nn_hww600_ww", pdfSets[i], "/smurf/ceballos/tmva/output_V5/ntuples_600train_0jets_hww600.root");

    }

    /*
       looper->ScanChain("bdtg_hww115_ww", "/smurf/ceballos/tmva/output_V5/ntuples_115train_1jets_hww115.root");
       looper->ScanChain("bdtg_hww120_ww", "/smurf/ceballos/tmva/output_V5/ntuples_120train_1jets_hww120.root");
       looper->ScanChain("nn_hww130_ww", "/smurf/ceballos/tmva/output_V5/ntuples_130train_1jets_hww130.root");
       looper->ScanChain("bdtg_hww140_ww", "/smurf/ceballos/tmva/output_V5/ntuples_140train_1jets_hww140.root");
       looper->ScanChain("bdtg_hww150_ww", "/smurf/ceballos/tmva/output_V5/ntuples_150train_1jets_hww150.root");
       looper->ScanChain("bdtg_hww160_ww", "/smurf/ceballos/tmva/output_V5/ntuples_160train_1jets_hww160.root");
       looper->ScanChain("bdtg_hww170_ww", "/smurf/ceballos/tmva/output_V5/ntuples_170train_1jets_hww170.root");
       looper->ScanChain("bdtg_hww180_ww", "/smurf/ceballos/tmva/output_V5/ntuples_180train_1jets_hww180.root");
       looper->ScanChain("bdtg_hww190_ww", "/smurf/ceballos/tmva/output_V5/ntuples_190train_1jets_hww190.root");
       looper->ScanChain("bdtg_hww200_ww", "/smurf/ceballos/tmva/output_V5/ntuples_200train_1jets_hww200.root");
       looper->ScanChain("bdtg_hww250_ww", "/smurf/ceballos/tmva/output_V5/ntuples_250train_1jets_hww250.root");
       looper->ScanChain("bdtg_hww300_ww", "/smurf/ceballos/tmva/output_V5/ntuples_300train_1jets_hww300.root");
       looper->ScanChain("bdtg_hww350_ww", "/smurf/ceballos/tmva/output_V5/ntuples_350train_1jets_hww350.root");
       looper->ScanChain("bdtg_hww400_ww", "/smurf/ceballos/tmva/output_V5/ntuples_400train_1jets_hww400.root");
       looper->ScanChain("bdtg_hww450_ww", "/smurf/ceballos/tmva/output_V5/ntuples_450train_1jets_hww450.root");
       looper->ScanChain("bdtg_hww500_ww", "/smurf/ceballos/tmva/output_V5/ntuples_500train_1jets_hww500.root");
       looper->ScanChain("bdtg_hww550_ww", "/smurf/ceballos/tmva/output_V5/ntuples_550train_1jets_hww550.root");
       looper->ScanChain("bdtg_hww600_ww", "/smurf/ceballos/tmva/output_V5/ntuples_600train_1jets_hww600.root");
     */

    //looper->ScanChain("LR[10]", "/tas03/home/yygao/WWAnalysis/CMSSW_3_11_3_ME_SmurfV5/src/Smurf/ME/scripts/output_V5/hww160_LR_hww160.root");


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

}

