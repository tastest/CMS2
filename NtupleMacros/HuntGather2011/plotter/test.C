
{

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gROOT->ProcessLine(".L ../libs/libHuntGather2011Plotter.so");
    gROOT->ProcessLine(".L ../libs/libCMS2NtupleMacrosCORE.so");
    gROOT->ProcessLine(".L ../libs/libCMS2NtupleMacrosTools.so");
    gROOT->ProcessLine(".L ../libs/libHuntGather2011Babymaker.so");

    //gROOT->ProcessLine(".L ../../Tools/goodrun.cc++");

    const char *goodrunlist = "../runlists/Cert_TopNov5_Merged_135821-149442_allPVT.txt";
    set_goodrun_file(goodrunlist);

    TFile *_file0 = TFile::Open("/nfs-3/userdata/cms2/gather/data/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather_skimmed_ntuple_146430_0.root");

    tree->Draw("pt1", "(pt1>20.0&&pt2>20.0&&goodrun(run,ls))");

    //tree->Draw("pt1", "(pt1>20.0&&pt2>20.0)");

}

