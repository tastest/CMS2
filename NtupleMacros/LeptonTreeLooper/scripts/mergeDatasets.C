void mergeDatasets() {

    TString base = "/smurf/dlevans/LeptonTree/V00-02-07/";
    TString file = "merged_HCP.root";

    std::vector<TString> pds;
    //pds.push_back("DoubleElectron");
    pds.push_back("DoubleMu");
    //pds.push_back("SingleElectron");
    //pds.push_back("SingleMu");

    std::vector<TString> eras;
    eras.push_back("Run2012A-13Jul2012-v1_AOD_190456_193621");
    eras.push_back("Run2012A-recover-06Aug2012-v1_AOD_190782_190949");
    eras.push_back("Run2012B-13Jul2012-v1_AOD_193834_196531");
    eras.push_back("Run2012C-24Aug2012-v1_AOD_198022_198523");
    eras.push_back("Run2012C-PromptReco-v2_AOD_198934_202950");

    for (unsigned int i = 0; i < pds.size(); ++i) {

        std::cout << "Doing " << pds[i] << std::endl;
        gSystem->mkdir(base+pds[i]);
        TChain *chain = new TChain("leptons");
        chain->SetMaxTreeSize(5000000000LL);

        for (unsigned int j = 0; j < eras.size(); ++j) {

            TString era = eras[j];
            if (pds[i] == "DoubleMu" && era == "Run2012B-13Jul2012-v1_AOD_193834_196531")
                era = "Run2012B-13Jul2012-v4_AOD_193834_196531";
            std::cout << "\t" << eras[j] << std::endl;
            //std::cout << TString(base+pds[i]+"_"+eras[j]+"/"+file) << std::endl;
            chain->Add(base+pds[i]+"_"+era+"/"+file);
        }

        std::cout << "\tMerging" << std::endl;
        //std::cout << TString(base+pds[i]+"/"+file) << std::endl;
        chain->Merge(base+pds[i]+"/"+file);
        delete chain;
    }

}

