#ifndef __CINT__
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "histtools.h"
#include "ossusy_looper.h"

#include <iostream>
#endif

void pickSkimIfExists( TChain *ch, const std::string& base, const std::string& skimPrefix ) {
    TChain *dummy = new TChain("Events");

    /*
       if (skimPrefix != "") {
       std::string skimName = Form("data/skim_v2/%s_skimmednTuple.root", skimPrefix.c_str());
       if (dummy->Add(skimName.c_str())) {
       int nFiles = ch->Add(skimName.c_str());
       std::cout << "Skim " << skimName.c_str() << " exists: use it. Loaded "<< nFiles << " files" << std::endl;
       return;
       } else
       std::cout << "Skim " << skimName.c_str() << " does not exist: use " << base.c_str() << std::endl;
       }
     */
    if (dummy->Add(base.c_str())) {
        int nFiles = ch->Add(base.c_str());
        std::cout << "Main " <<base.c_str() << " exists: use it. Loaded " 
            << nFiles << " files" << std::endl;
    } else
        std::cout << "FUCK SHIT DAMN!" << std::endl;

    // be paranoid
    if (nFiles == 0) {
        std::cout << "ERROR: expected to read files " 
            << base.c_str() << "  but found none" << std::endl;
        assert(0);
    }

    return;
}

void doAll_ossusy_looper(bool skipFWLite = true) {
    // Load various tools  
    gROOT->ProcessLine(Form(".x setup.C(%d)", skipFWLite));

    // Load and compile the looping code
    gSystem->CompileMacro("ossusy_looper.C", 
            "++k", "libossusy_looper");

    ossusy_looper* looper = new ossusy_looper();
    //use OS/SS baseline selection as documented in:
    //http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/SusyStudies3x
    looper->set_susybaseline(0);
    //make baby ntuple
    looper->set_createTree(0);
    //use bitmask selection
    looper->set_useBitMask(0);

    // K-factors
    // these have been k-factors NLO/LO before
    // now using them as sample normalizations to NLO

    // these two are taken from Ceballos's pdf. 
    // It looks like the top x-section is for mtop = 175 GeV
    float kttdil    = 1.;  // 375pb, 127000 events processed
    float kttotr    = 1.;  // 375pb, 127000 events processed
    float kWW       = 1.;
    float kWZ       = 1.;
    float kZZ       = 1.;
    float kWjets    = 1.;  // 11850 pb, 980000 events processed
    float kWcharm   = 1.1;
    float kZjets    = 1.;  
    float kDYee     = 1.;  // 1230 pb,  970360 events processed
    float kDYmm     = 1.;  // 1230 pb,  970360 events processed
    float kDYtautau = 1.;  // 1230 pb,  970360 events processed
    float kppMuX    = 1.;  // xsec/nevents
    float kEM       = 1.;
    float ktW       = 1.;  // the evtScale is all negative for some reason
    float kVQQ      = 1.;
    float kLM0      = 1.;
    float kLM1      = 1.;
    float kLM2      = 1.;
    float kLM3      = 1.;
    float kLM4      = 1.;
    float kLM5      = 1.;
    float kLM6      = 1.;
    float kLM7      = 1.;
    float kLM8      = 1.;
    float kLM9      = 1.;
    float kLM10     = 1.;
    float kLM11     = 1.;
    float kLM12     = 1.;

    // Prescales
    int prettdil    = 1;
    int prettotr    = 1;
    int preWW       = 1;
    int preWZ       = 1;
    int preZZ       = 1;
    int preWjets    = 1;
    int preWcharm   = 1;
    int preZjets    = 1;
    int preDYee     = 1;
    int preDYmm     = 1;
    int preDYtautau = 1;
    int preppMuX    = 1;
    int preEM       = 1;
    int pretW       = 1;
    int preVQQ      = 1;
    int preLM0      = 1;
    int preLM1      = 1;
    int preLM2      = 1;
    int preLM3      = 1;
    int preLM4      = 1;
    int preLM5      = 1;
    int preLM6      = 1;
    int preLM7      = 1;
    int preLM8      = 1;
    int preLM9      = 1;
    int preLM10     = 1;
    int preLM11     = 1;
    int preLM12     = 1;

    //Flags for files to run over
    bool runttdil    = 1;
    bool runttotr    = 1;
    bool runWW       = 1;
    bool runWZ       = 1;
    bool runZZ       = 1;
    bool runWjets    = 1;
    bool runWcharm   = 0;
    bool runZjets    = 1;
    bool runDYee     = 0;
    bool runDYmm     = 0;
    bool runDYtautau = 0;
    bool runppMuX    = 0;
    bool runEM       = 0;
    bool runtW       = 1;
    bool runVQQ      = 0;
    bool runLM0      = 1;
    bool runLM1      = 1;
    bool runLM2      = 1;
    bool runLM3      = 1;
    bool runLM4      = 1;
    bool runLM5      = 1;
    bool runLM6      = 1;
    bool runLM7      = 1;
    bool runLM8      = 1;
    bool runLM9      = 1;
    bool runLM10     = 1;
    bool runLM11     = 1;
    bool runLM12     = 1;

    TChain* chZjets = new  TChain("Events");
    if(runZjets){
        pickSkimIfExists(chZjets, 
                //"data3x/ZJets-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged*.root", 
                "data3x/ZJets-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/skimDilepton-20-10/skimmed*.root", 
                "Zjets");
    }

    TChain* chtopdil = new TChain("Events");
    if (runttdil) {
        pickSkimIfExists(chtopdil, 
                "data3x/TTbarJets-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged*.root",
                "TTJets");
    }

    TChain* chtopotr = new TChain("Events");
    if(runttotr){
        pickSkimIfExists(chtopotr, 
                "data3x/TTbarJets-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged*.root",
                "TTJets");
    }

    TChain* chww = new TChain("Events");
    if(runWW){
        pickSkimIfExists(chww, 
                "data3x/WW_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root",
                "WW");
    }

    TChain* chWZ = new TChain("Events");
    if(runWZ){
        pickSkimIfExists(chWZ, 
                "data3x/WZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "WZ_incl"); // can try WZ_3l-Pythia
    }

    TChain* chZZ = new TChain("Events");
    if(runZZ){
        pickSkimIfExists(chZZ, 
                "data3x/ZZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "ZZ");
    }

    TChain* chWjets = new  TChain("Events");
    if(runWjets){
        pickSkimIfExists(chWjets, 
                //"data3x/WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "data3x/WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/skimDilepton-20-10/skimmed*.root",
                "WJets");
    }

    //SAMPLE NOT YET AVAILABLE
    //   TChain* chWcharm = new TChain("Events");
    //   if(runWcharm){
    //     pickSkimIfExists(chWcharm, 
    // 		     "data/Wc-madgraph_Fall08_IDEAL_V11_redigi_v1/merged*.root", 
    // 		     "Wc");
    //   }

    TChain* chDYtautau = new  TChain("Events");
    if(runDYtautau){
        pickSkimIfExists(chDYtautau, 
                "data3x/Ztautau_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "DYtautau");
    }

    // the low-mass splice has no choice other than the skim
    //    pickSkimIfExists(chDYtautau, "data/Ztautau_M20_Summer08_IDEAL_V11_redigi_v1/merged*.root", "Ztautau_M20");

    TChain* chDYee = new  TChain("Events");
    if(runDYee){
        pickSkimIfExists(chDYee, 
                "data3x/Zee_Summer09-MC_31X_V3_7TeV_TrackingParticles-v1/V03-00-35/merged*.root", 
                "DYee");
    }
    // the low-mass splice has no choice other than the skim
    //    pickSkimIfExists(chDYee, "data/Zee_M20_Summer08_IDEAL_V11_redigi_v1/merged*.root", "Zee_M20");

    if(runDYmm){
        TChain* chDYmm = new  TChain("Events");
        pickSkimIfExists(chDYmm, 
                "data3x/Zmumu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "DYmm");
    }
    // the low-mass splice has no choice other than the skim
    //    pickSkimIfExists(chDYmm, "data/Zmumu_M20_Summer08_IDEAL_V11_redigi_v1/merged*.root", "Zmumu_M20");

    // ppMuX
    TChain* chppMuX = new  TChain("Events");
    if (runppMuX) {
        pickSkimIfExists(chppMuX, 
                "data3x/InclusiveMu15_Summer09-MC_31X_V3_7TeV-v1_dilepfilt/merged*.root", 
                "InclusiveMuPt15"); 
        // can try InclusiveMu5Pt50 .. figure out how to merge later
    }

    // ppEM
    TChain* chEM =  new  TChain("Events");
    if (runEM) {
        pickSkimIfExists(chEM, 
                "data3x/QCD_EMEnriched_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "_skimSimple2020");
        pickSkimIfExists(chEM, 
                "data3x/QCD_EMEnriched_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "_skimSimple2020");
        pickSkimIfExists(chEM, 
                "data3x/QCD_EMEnriched_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "_skimSimple2020");
        pickSkimIfExists(chEM, 
                "data3x/QCD_BCtoE_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "_skimSimple2020");
        pickSkimIfExists(chEM, 
                "data3x/QCD_BCtoE_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "_skimSimple2020");
        pickSkimIfExists(chEM, 
                "data3x/QCD_BCtoE_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "_skimSimple2020");
    }

    // tW
    TChain* chtW = new  TChain("Events");
    if (runtW) {
        pickSkimIfExists(chtW, 
                "data3x/SingleTop_sChannel-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SingleTop_sChannel"); 
        pickSkimIfExists(chtW, 
                "data3x/SingleTop_tChannel-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged*.root", 
                "SingleTop_tChannel"); 
        pickSkimIfExists(chtW, 
                "data3x/SingleTop_tWChannel-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged*.root",
                "SingleTop_tWChannel"); 
    }

    //SAMPLE NOT YET AVAILABLE
    //   // VQQ
    //   TChain* chVQQ = new TChain("Events");
    //   if (runVQQ) {
    //     pickSkimIfExists(chVQQ, 
    // 		     "data/VQQ-madgraph_Summer08_IDEAL_V11_redigi_v2/merged*.root", 
    // 		     "VQQ");
    //   }

    // LM0
    TChain *chLM0 = new TChain("Events");
    if (runLM0) {
        pickSkimIfExists(chLM0, 
                "data3x/LM0_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM0");
    }

    // LM1
    TChain *chLM1 = new TChain("Events");
    if (runLM1) {
        pickSkimIfExists(chLM1, 
                "data3x/LM1_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM1");
    }

    // LM2
    TChain *chLM2 = new TChain("Events");
    if (runLM2) {
        pickSkimIfExists(chLM2, 
                "data3x/LM2_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM2");
    }

    // LM3
    TChain *chLM3 = new TChain("Events");
    if (runLM3) {
        pickSkimIfExists(chLM3, 
                "data3x/LM3_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM3");
    }

    // LM4
    TChain *chLM4 = new TChain("Events");
    if (runLM4) {
        pickSkimIfExists(chLM4, 
                "data3x/LM4_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM4");
    }

    // LM5
    TChain *chLM5 = new TChain("Events");
    if (runLM5) {
        pickSkimIfExists(chLM5, 
                "data3x/LM5_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM5");
    }

    // LM6
    TChain *chLM6 = new TChain("Events");
    if (runLM6) {
        pickSkimIfExists(chLM6, 
                "data3x/LM6_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM6");
    }

    // LM7
    TChain *chLM7 = new TChain("Events");
    if (runLM7) {
        pickSkimIfExists(chLM7, 
                "data3x/LM7_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM7");
    }

    // LM8
    TChain *chLM8 = new TChain("Events");
    if (runLM8) {
        pickSkimIfExists(chLM8, 
                "data3x/LM8_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM8");
    }

    // LM9
    TChain *chLM9 = new TChain("Events");
    if (runLM9) {
        pickSkimIfExists(chLM9, 
                "data3x/LM9_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM9");
    }

    // LM10
    TChain *chLM10 = new TChain("Events");
    if (runLM10) {
        pickSkimIfExists(chLM10, 
                "data3x/LM10_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM10");
    }

    // LM11 
    TChain *chLM11 = new TChain("Events");
    if (runLM11) {
        pickSkimIfExists(chLM11, 
                "data3x/LM11_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM11");
    }

    // LM12
    TChain *chLM12 = new TChain("Events");
    if (runLM12) {
        pickSkimIfExists(chLM12, 
                "data3x/LM12_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                "SUSY_LM12");
    }

    char jetTypeStrings[2][5] = {"JPT", "calo"};
    char metTypeStrings[3][8] = {"tcmet", "muon", "muonjes"};
    bool doFakeApp            = false;

    // Process files one at a time, and color them as needed
    for (int jetTypeIdx = 1; jetTypeIdx < 2; ++jetTypeIdx)
    {
        for (int metTypeIdx = 0; metTypeIdx < 1; ++metTypeIdx)
        {
            ossusy_looper::JetTypeEnum jetType(jetTypeIdx);
            ossusy_looper::MetTypeEnum metType(metTypeIdx);
            if (runZjets) {
                cout << "Processing Zjets" << endl;
                looper->ScanChain(chZjets,"Zjets", kZjets, preZjets, jetType, metType, doFakeApp);
                cout << "Done processing Zjets" << endl;
                hist::color("Zjets", kBlack);
            }
            if (runttdil) {
                cout << "Processing ttbar dileptonic.. " << endl;
                looper->ScanChain(chtopdil,"ttdil", kttdil, prettdil, jetType, metType, doFakeApp);
                cout << "Done processing ttbar dileptonic.. " << endl;
                hist::color("ttdil", kYellow);
            }
            if (runttotr) {
                cout << "Processing ttbar no-dileptons.. " << endl;
                looper->ScanChain(chtopotr,"ttotr", kttotr, prettotr, jetType, metType, doFakeApp);
                cout << "Done processing ttbar no-dileptons.. " << endl;
                hist::color("ttotr", 30);
            }
            if (runWW) {
                cout << "Processing WW.." << endl;
                looper->ScanChain(chww,"ww", kWW, preWW, jetType, metType, doFakeApp);
                cout << "Done processing WW.." << endl;
                hist::color("ww", kRed);
            }
            if (runWZ) {
                cout << "Processing WZ.." << endl;
                looper->ScanChain(chWZ,"wz", kWZ, preWZ, jetType, metType, doFakeApp);
                cout << "Done processing WZ.." << endl;
                hist::color("wz", kBlue);
            }
            if (runZZ) {
                cout << "Processing ZZ.." << endl;
                looper->ScanChain(chZZ,"zz", kZZ, preZZ, jetType, metType, doFakeApp);
                cout << "Done processing ZZ.." << endl;
                hist::color("zz", kGreen);
            }
            if (runWjets) {
                cout << "Processing Wjets.." << endl;
                looper->ScanChain(chWjets,"wjets", kWjets, preWjets, jetType, metType, doFakeApp);
                cout << "Done processing Wjets.." << endl;
                hist::color("wjets", 40);
            }
            // 	  if (runWcharm) { //SAMPLE NOT YET AVAILABLE
            // 	    cout << "Processing Wcharm.." << endl;
            // 	    looper->ScanChain(chWcharm, "wcharm", kWcharm, preWcharm, jetType, metType, doFakeApp);
            // 	    cout << "Done processing Wcharm.." << endl;
            // 	    hist::color("wcharm", 50);
            // 	  }
            if (runDYtautau) {
                cout << "Processing DY->tautau" << endl;
                looper->ScanChain(chDYtautau,"DYtautau", kDYtautau, preDYtautau, jetType, metType, doFakeApp);
                cout << "Done processing DY->tautau" << endl;
                hist::color("DYtautau", kBlack);
            }
            if (runDYee) {
                cout << "Processing DY->ee" << endl;
                looper->ScanChain(chDYee,"DYee", kDYee, preDYee, jetType, metType, doFakeApp);
                cout << "Done rocessing DY->ee" << endl;
                hist::color("DYee", kMagenta);
            }
            if (runDYmm) {
                cout << "Processing DY->mm" << endl;
                looper->ScanChain(chDYmm,"DYmm", kDYmm, preDYmm, jetType, metType, doFakeApp);
                cout << "Done processing DY->mm" << endl;
                hist::color("DYmm", kCyan);
            }
            if (runppMuX) {
                cout << "Processing ppMuX" << endl;
                looper->ScanChain(chppMuX,"ppMuX", kppMuX, preppMuX, jetType, metType, doFakeApp);
                cout << "Done processing ppMuX" << endl;
                hist::color("ppMuX", 51);
            }
            if (runEM) {
                cout << "Processing EM" << endl;
                looper->ScanChain(chEM,"EM", kEM, preEM, jetType, metType, doFakeApp);
                cout << "Done processing EM" << endl;
                hist::color("EM", 49);
            }
            if (runtW) {
                cout << "Processing tW" << endl;
                looper->ScanChain(chtW,"tW", ktW, pretW, jetType, metType, doFakeApp);
                cout << "Done processing tW" << endl;
                hist::color("tW", 63);
            }
            // 	  if (runVQQ) { //SAMPLE NOT YET AVAILABLE
            // 	    cout << "Processing VQQ" << endl;
            // 	    looper->ScanChain(chVQQ,"VQQ", kVQQ, preVQQ, jetType, metType, doFakeApp);
            // 	    cout << "Done processing VQQ" << endl;
            // 	    hist::color("VQQ", 45);
            // 	  }
            if (runLM0) {
                cout << "Processing LM0" << endl;
                looper->ScanChain(chLM0, "LM0", kLM0, preLM0, jetType, metType, doFakeApp);
                cout << "Done processing LM0" << endl;
                hist::color("LM0", kOrange);
            }
            if (runLM1) {
                cout << "Processing LM1" << endl;
                looper->ScanChain(chLM1, "LM1", kLM1, preLM1, jetType, metType, doFakeApp);
                cout << "Done processing LM1" << endl;
                hist::color("LM1", kOrange+1);
            }
            if (runLM2) {
                cout << "Processing LM2" << endl;
                looper->ScanChain(chLM2, "LM2", kLM2, preLM2, jetType, metType, doFakeApp);
                cout << "Done processing LM2" << endl;
                hist::color("LM2", kOrange+2);
            }
            if (runLM3) {
                cout << "Processing LM3" << endl;
                looper->ScanChain(chLM3, "LM3", kLM3, preLM3, jetType, metType, doFakeApp);
                cout << "Done processing LM3" << endl;
                hist::color("LM3", kOrange+3);
            }
            if (runLM4) {
                cout << "Processing LM4" << endl;
                looper->ScanChain(chLM4, "LM4", kLM4, preLM4, jetType, metType, doFakeApp);
                cout << "Done processing LM4" << endl;
                hist::color("LM4", kOrange+4);
            }
            if (runLM5) {
                cout << "Processing LM5" << endl;
                looper->ScanChain(chLM5, "LM5", kLM5, preLM5, jetType, metType, doFakeApp);
                cout << "Done processing LM5" << endl;
                hist::color("LM5", kOrange+5);
            }
            if (runLM6) {
                cout << "Processing LM6" << endl;
                looper->ScanChain(chLM6, "LM6", kLM6, preLM6, jetType, metType, doFakeApp);
                cout << "Done processing LM6" << endl;
                hist::color("LM6", kOrange+6);
            }
            if (runLM7) {
                cout << "Processing LM7" << endl;
                looper->ScanChain(chLM7, "LM7", kLM7, preLM7, jetType, metType, doFakeApp);
                cout << "Done processing LM7" << endl;
                hist::color("LM7", kOrange+7);
            }
            if (runLM8) {
                cout << "Processing LM8" << endl;
                looper->ScanChain(chLM8, "LM8", kLM8, preLM8, jetType, metType, doFakeApp);
                cout << "Done processing LM8" << endl;
                hist::color("LM8", kOrange+8);
            }
            if (runLM9) {
                cout << "Processing LM9" << endl;
                looper->ScanChain(chLM9, "LM9", kLM9, preLM9, jetType, metType, doFakeApp);
                cout << "Done processing LM9" << endl;
                hist::color("LM9", kOrange+9);
            }
            if (runLM10) {
                cout << "Processing LM10" << endl;
                looper->ScanChain(chLM10, "LM10", kLM10, preLM10, jetType, metType, doFakeApp);
                cout << "Done processing LM10" << endl;
                hist::color("LM10", kOrange+10);
            }
            if (runLM11) { 
                cout << "Processing LM11" << endl;
                looper->ScanChain(chLM11, "LM11", kLM11, preLM11, jetType, metType, doFakeApp);
                cout << "Done processing LM11" << endl;
                hist::color("LM11", kOrange-7);
            }
            if (runLM12) {
                cout << "Processing LM12" << endl;
                looper->ScanChain(chLM12, "LM12", kLM12, preLM12, jetType, metType, doFakeApp);
                cout << "Done processing LM12" << endl;
                hist::color("LM12", kOrange-7);
            }

            // save all the histograms
            if(doFakeApp) {
                const char* outFile = Form("victory_baseline_metgt50_sumjetptgt200_cand01_%s_%s_FakeApp_3x.root", 
                        jetTypeStrings[jetTypeIdx], metTypeStrings[metTypeIdx]);
            }
            else {
                const char* outFile = Form("victory_baseline_metgt50_sumjetptgt200_cand01_%s_%s_3x.root", 
                        jetTypeStrings[jetTypeIdx], metTypeStrings[metTypeIdx]);
            }
            //const char* outFile = Form("root3x/victory_susybaseline_%s_%s_3x_newemselectors.root", 
            //			     jetTypeStrings[jetTypeIdx], metTypeStrings[metTypeIdx]);

            //const char* outFile = Form("victory_baseline_genmetgt50_nosumjetptcut_%s_%s_pleasework_varbins.root", 
            //jetTypeStrings[jetTypeIdx], metTypeStrings[metTypeIdx]);
            TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
            rootdir->cd();
            saveHist(outFile);
            deleteHistos();
        } // metTypeIdx
    } // jetTypeIdx

    gSystem->Exit(0);
}
