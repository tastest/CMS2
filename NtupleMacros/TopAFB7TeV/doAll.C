void doAll(TString outputDir="results", bool rundata=true, bool runsig=false, bool runmc = true, bool requireBTag=true, bool require2BTag=false, bool usePF = true, 
	   bool doFRestimation = false, bool sendOutputToLogFile = true, bool BTagAlgTCHE = false, bool createBabyNtuples = true, bool doBFR = false)
{
  //gSystem->Load("/home/users/yanjuntu/MiniFWlib/libMiniFWLite_v5.28.00.so");
  // gSystem->Load("/home/users/yanjuntu/MiniFWlib/libMiniFWLite_5.27.06b-cms10.so");
   gSystem->Load("/home/users/yanjuntu/MiniFWlib/libMiniFWLite_CMSSW_5_3_2_patch4_V05-03-13.so");
  //gSystem->Load("/nfs-3/userdata/yanjuntu/lhapdf/lib/libLHAPDF.so");
  gSystem->Load("/nfs-6/userdata/yanjuntu/LHAPDF/lib/libLHAPDF.so");
  
  gSystem->AddIncludePath(" -w -I../CORE/topmass -I/nfs-6/userdata/yanjuntu/LHAPDF/include");
  gROOT->ProcessLine(".L ../CORE/topmass/ttdilepsolve.cpp+"); 
  gROOT->ProcessLine(".L ../CORE/CMS2.cc+");
  gROOT->ProcessLine(".L ../CORE/utilities.cc+");
  gROOT->ProcessLine(".L ../CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/MITConversionUtilities.cc+");
  gROOT->ProcessLine(".L ../CORE/muonSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L ../CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/metSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/SimpleFakeRate.cc+");
  gROOT->ProcessLine(".L ../CORE/mcSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/MT2/MT2.cc+");
  gROOT->ProcessLine(".L ../CORE/triggerUtils.cc+");  
  gROOT->ProcessLine(".L ../CORE/susySelections.cc+");
  //gROOT->ProcessLine(".L ../CORE/mcSUSYkfactor.cc+");
  gROOT->ProcessLine(".L ../CORE/triggerSuperModel.cc+");
  gROOT->ProcessLine(".L ../CORE/jetSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/ttbarSelections.cc+");
  gROOT->ProcessLine(".L ../CORE/triggerUtils.cc+");
  

  gSystem->CompileMacro("histtools.C", "++k", "libhisttools");
  gSystem->CompileMacro("topAFB_looper.C","++k", "libtopAFB_looper");

  //float lumiToNormalizeTo   = 1144.5e-3; //1090.0e-3; //349.0e-3; //204.0e-3;// 191.0e-3 ; //36.1e-3; 
  float lumiToNormalizeTo   = 4980.0e-3; //4684.0e-3; //3230.0e-3; //1090.0e-3; //349.0e-3; //204.0e-3;// 191.0e-3 ; //36.1e-3; 

  topAFB_looper * baby = new topAFB_looper();
  // Flags for files to run over
  //string cms2_skim_location="/nfs-6/userdata/yanjuntu/TPrimeSkim";
  string cms2_skim_location="/nfs-6/userdata/yanjuntu/AfbSkimv2";
  string cms2_location="/nfs-6/userdata/cms2";

  bool runskim          = true;
  bool runSMS           = false;
  bool runwprime        = runsig;
  bool runAxigluon      = runsig;
  bool runttdil         = runmc;
  bool runttdilsys      = false;

  bool runttotr         = runmc ; //false;
  bool runWjets         = runmc ; //false;
  bool runDYee          = runmc ; //false;
  bool runDYmm          = runmc ; //false;
  bool runDYtautau      = runmc ; //false;
  bool  runVV           = runmc ; //false;
  bool  runtW           = runmc ; //false;


  bool  runQCDPt15      = false;
  bool  runQCDPt30      = false;
  bool  runVqq          = false;
  bool  runphotonJet15  = false;
  bool  runphotonVJets  = false;
  //NLO cross-sections
  float kWprime400    = 1.;
  float kWprime600    = 1.;
  float kWprime1000    = 1.;
  float kWprime600w300    = 1.;
  float kAxigluonR    = 1.;
  float ksms       = 1.;
  float kttdil    = 1.;
  float kttotr    = 1.;
  float kWjets    = 1.;
  float kDYee     = 1.;
  float kDYmm     = 1.;
  float kDYtautau = 1.;
  float kVV       = 1.;
  float ktW       = 1.;
  float kqcd15    = 1.;
  float kqcd30    = 1.;
  float kVqq      = 1.;
  float kphoton15 = 1.;

    vector<TString> v_baseCuts;
  //v_baseCuts.push_back("applyNoCuts");         // no cuts, set runskim=false too
  v_baseCuts.push_back("usePtGt2020");         // use leptons with pt > 20
  //v_baseCuts.push_back("applyTriggers");       // apply triggers
  v_baseCuts.push_back("hypDisamb");           // do hyp. disambiguation
  //v_baseCuts.push_back("usePtGt2010");         // one lepton > 20, other > 10
  //v_baseCuts.push_back("excludePtGt2020");     // one lepton > 10, < 20 other >20

  if(!usePF) {
    v_baseCuts.push_back("usejptJets");          // use jpt jets for jet counting      
    v_baseCuts.push_back("usetcMET");            // use tcMET
  } else {
    v_baseCuts.push_back("usepfMET");   //use PFMET   
    v_baseCuts.push_back("usepfJets");  // use pf jets for jet counting
  }


  //v_baseCuts.push_back("requireEcalEls");
  v_baseCuts.push_back("useOS");
  //if(rundata)
  //v_baseCuts.push_back("applyAlignmentCorrection");
  // v_baseCuts.push_back("vetoHypMassLt10");
  v_baseCuts.push_back("vetoHypMassLt12");

  //possible cuts
  //v_baseCuts.push_back("applyFOv1Cuts");
  //v_baseCuts.push_back("applyFOv2Cuts");  
  //v_baseCuts.push_back("applyFOv3Cuts");
  //v_baseCuts.push_back("applyLooseIDCuts");      // apply loose ID cuts
  //v_baseCuts.push_back("applylepLooseIsoCuts");// loose iso cuts
  //v_baseCuts.push_back("requireZmass");        // leptons only in zmass
  //v_baseCuts.push_back("useCorMET");           // use corrected calo MET ---> NOT SUPPORTED RIGHT NOW
  //v_baseCuts.push_back("vetoProjectedMET");    // cut on projected MET
  //v_baseCuts.push_back("usecaloJets");         // use caloJETs for jet counting
  //v_baseCuts.push_back("useSS");
  
  // -----------
  // top pt reweighting
  // activate to reweight events according to top pt
  // or alterniative scenarios
  //
  // v_baseCuts.push_back("applyTopPtWeighting");
  // v_baseCuts.push_back("applyLeptonPtWeighting");
  // v_baseCuts.push_back("applyJetPtWeighting");
  //
  // default is using sqrt(weight(object_1)*weight(object_2)) for 100% correlated objects
  // other possibilities are
  // v_baseCuts.push_back("useReweightingUncorrelated"); // for completely uncorrelated objects
  // v_baseCuts.push_back("useReweightingLeadingObject"); // to only reweight according to the leading object
  
  if(requireBTag && !doBFR)
    v_baseCuts.push_back("requireBTag");  
  if(require2BTag && !doBFR)
    v_baseCuts.push_back("require2BTag");
  if(doBFR)
    v_baseCuts.push_back("doBFR");
  v_baseCuts.push_back("sortJetCandidatesbyPt");
  //v_baseCuts.push_back("sortJetCandidatesbyDR");
  // v_baseCuts.push_back("matchLeptonJetbyMaxDR");
  //v_baseCuts.push_back("applyLeptonJetInvMassCut450");
  //v_baseCuts.push_back("applyTopSystEta");
  // v_baseCuts.push_back("applyMinMassLBCut");
  //v_baseCuts.push_back("requireExact2BTag");
  
  v_baseCuts.push_back("generalLeptonVeto");
  //v_baseCuts.push_back("applyHTCut");
  if(BTagAlgTCHE)
    v_baseCuts.push_back("BTagAlgTCHE");
  if(createBabyNtuples)
    v_baseCuts.push_back("createBabyNtuples");

  // systematics
  // activate if needed
  // v_baseCuts.push_back("scaleJER");
  // v_baseCuts.push_back("scaleJESMETUp");
  // v_baseCuts.push_back("scaleJESMETDown");
  // v_baseCuts.push_back("scaleLeptonEnergyUp");
  // v_baseCuts.push_back("scaleLeptonEnergyDown");
  // v_baseCuts.push_back("scaleBTAGSFup");
  // v_baseCuts.push_back("scaleBTAGSFdown");
  // v_baseCuts.push_back("scaleTrigSFup");
  // v_baseCuts.push_back("scaleTrigSFdown");
  // v_baseCuts.push_back("noVertexReweighting");
  // v_baseCuts.push_back("weighttaudecay");
  
  vector<TString> v_Cuts = v_baseCuts;
  vector<TString> v_otherCuts;
  if(!doFRestimation) {

    //base set of cuts that will be used for the yields
    //v_otherCuts.push_back("applylepIDCuts;applylepIsoCuts");
    //v_otherCuts.push_back("applylepIDCuts;applylepIsoCuts;vetoMET;veto2Jets");
    
    //reduced cuts for DYEst
    //v_otherCuts.push_back("applylepIDCuts;applylepIsoCuts;veto2Jets");
    
    //full cuts
    v_otherCuts.push_back("applylepIDCuts;applylepIsoCuts;vetoZmass;veto2Jets;vetoMET");

	// MET > 50 
    //v_otherCuts.push_back("applylepIDCuts;applylepIsoCuts;vetoZmass;veto2Jets;vetoMET50");

    // //this is to make NJet plots since we don't want to make plots with the njet cut in there
    //v_otherCuts.push_back("applylepIDCuts;applylepIsoCuts;vetoZmass");
    //v_otherCuts.push_back("applylepIDCuts;applylepIsoCuts;vetoZmass;vetoMET");

  }
  else {
    
    v_otherCuts.push_back("vetoZmass;veto2Jets;vetoMET");
    
  }
  for(unsigned int index = 0; index < v_otherCuts.size(); index++) {

   
    vector<TString> v_Cuts = v_baseCuts;
    
    if(v_otherCuts.at(index).Contains(";")) {
      TString cut = v_otherCuts.at(index);
      TObjArray *ob = cut.Tokenize(";");
      for(int j = 0; j < ob->GetEntries(); j++) {       
	v_Cuts.push_back(TString(ob->At(j)->GetName()));
      }
    } else {
      v_Cuts.push_back(v_otherCuts.at(index));
    }

  
  if(rundata) {
    TChain  *ch_data= new TChain("Events");  
    cout << "Doing data" << endl;
    if(runskim){
      
      ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_4_V04-02-20/DoubleElectron_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
      ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_4_V04-02-20/DoubleMu_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
      // ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_4_V04-02-20/MuEG_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/MuEG_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
     
      ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_4_V04-02-20/DoubleElectron_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
      ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_4_V04-02-20/DoubleMu_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
      //ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_4_V04-02-20/MuEG_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/MuEG_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
      
      ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_7_patch1_V04-02-30/DoubleElectron_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_7_patch1_V04-02-30/DoubleMu_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      // ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_7_patch1_V04-02-30//MuEG_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30//MuEG_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");

      ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_7_patch1_V04-02-30/DoubleElectron_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_7_patch1_V04-02-30/DoubleMu_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      //ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_7_patch1_V04-02-30/MuEG_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/MuEG_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");

      ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_7_patch1_V04-02-30/DoubleElectron_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_7_patch1_V04-02-30/DoubleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      //ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_7_patch1_V04-02-30/MuEG_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/MuEG_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      
      ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_7_patch1_V04-02-34/DoubleElectron_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged*root");
      ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_7_patch1_V04-02-34/DoubleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged*root");
      // ch_data->Add("/nfs-6/userdata/yanjuntu/AfbSkimv2/CMSSW_4_2_7_patch1_V04-02-34/MuEG_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-34/MuEG_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged*root");
    }
    else{
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/DoubleElectron_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/DoubleMu_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/MuEG_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
     
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/DoubleElectron_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/DoubleMu_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/MuEG_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/merged*root");
      
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/DoubleElectron_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/DoubleMu_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30//MuEG_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
     
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/DoubleElectron_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/DoubleMu_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/MuEG_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/DoubleElectron_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/DoubleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/MuEG_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/merged*root");
      
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-34/DoubleElectron_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-34/DoubleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged*root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-34/MuEG_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/merged*root");
    
    }
    baby->ScanChain(ch_data, v_Cuts, "data", doFRestimation, lumiToNormalizeTo, 1, false);//, 0.01, 1, true);
    hist::color("data", kRed);
    delete ch_data;
  }
  if(runwprime){
    TChain  *ch_wprime400= new TChain("Events");
    TChain  *ch_wprime600= new TChain("Events");
    TChain  *ch_wprime1000= new TChain("Events");
    TChain  *ch_wprime600w300= new TChain("Events");

    //always run unskimmed

	cout << "Doing the MadGraph wprime 400 " << endl;
	ch_wprime400->Add("/nfs-6/userdata/yanjuntu/Wprime_SM_400_Madgraph_v2_yanjuntu-Wprime_SM_400_Madgraph_v2-f3d3f52ad6235ba5a3ccb05162c152b9_USER/VB04-02-29_Fastsim/merged_ntuple*.root");
	cout << "Doing the MadGraph wprime 600 " << endl;
        //ch_wprime600->Add("/nfs-6/userdata/yanjuntu/Wprime_ttbar_600_Madgraph-f3d3f52ad6235ba5a3ccb05162c152b9_USER/VB04-02-29_Fastsim/merged_ntuple*.root");
	ch_wprime600->Add("/nfs-6/userdata/yanjuntu/Wprime_ttbar_600_inclusive_Madgraph-6bfefe7ba7e693b7a5695fa468de1212_USER/VB04-02-29_Fastsim/merged_ntuple*.root");

        cout << "Doing the MadGraph wprime 1000 " << endl;
        ch_wprime1000->Add("/nfs-6/userdata/linacre/Wprime_ttbar_1000_Madgraph_linacre-Wprime_ttbar_1000_Madgraph-f0fcb122c9c1934e58d5f27499313ac3/VB04-02-29_Fastsim/merged_ntuple*.root");
        cout << "Doing the MadGraph wprime 600w300 " << endl;
        ch_wprime600w300->Add("/nfs-6/userdata/linacre/Wprime_ttbar_600_w300_2_Madgraph_linacre-Wprime_ttbar_600_w300_2_Madgraph-f3d3f52ad6235ba5a3ccb05162c152b9/VB04-02-29_Fastsim/merged_ntuple*.root");

  
    baby->ScanChain(ch_wprime400, v_Cuts,"wprime400", doFRestimation, lumiToNormalizeTo*6.2539/1.0, kWprime400, false);
    hist::color("wprime400", kRed); 
    baby->ScanChain(ch_wprime600, v_Cuts,"wprime600", doFRestimation, lumiToNormalizeTo*79997/99997, kWprime600, false);
    hist::color("wprime600", kRed);
    baby->ScanChain(ch_wprime1000, v_Cuts,"wprime1000", doFRestimation, lumiToNormalizeTo, kWprime1000, false);
    hist::color("wprime1000", kRed);
    baby->ScanChain(ch_wprime600w300, v_Cuts,"wprime600w300", doFRestimation, lumiToNormalizeTo, kWprime600w300, false);
    hist::color("wprime600w300", kRed);
  }


 if(runAxigluon){
    TChain  *ch_axigluonR= new TChain("Events");

    //always run unskimmed

	cout << "Doing the MadGraph Axigluon rght-handed " << endl;
	ch_axigluonR->Add("/nfs-3/userdata/cms2/AxigluonR_2TeV_ttbar_MadGraph_sergo-AxigluonR_2TeV_ttbar_MadGraph/VB04-02-29_Fastsim/merged_ntuple*.root");
      
    baby->ScanChain(ch_axigluonR, v_Cuts,"axigluonR", doFRestimation, lumiToNormalizeTo, kAxigluonR, false);
    hist::color("axigluonR", kRed); 
 }
  
if(runttdil) {
    TChain  *ch_ttbar= new TChain("Events");
    if(runskim){
      // cout << "Doing the MadGraph ttbar sample" << endl; ch_ttbar->Add(Form("%s/%s", cms2_skim_location.c_str(),"TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*root"));

      //cout << "Doing the powheg ttbar sample" << endl; ch_ttbar->Add(Form("%s/%s", cms2_skim_location.c_str(),"TTTo2L2Nu2B_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*root"));

      //cout << "Doing the powheg tauola ttbar sample" << endl; ch_ttbar->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/TT_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
      
      cout << "Doing the Fall11 MC@NLO ttbar sample" << endl; ch_ttbar->Add(Form("%s/%s", cms2_skim_location.c_str(),"TT_TuneZ2_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1/V04-02-29_fix_dilepton/skimmed*root"));
    } 
    else{
      //cout << "Doing the MadGraph ttbar sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    
      //cout << "Doing the powheg ttbar sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/TTTo2L2Nu2B_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
      
      //cout << "Doing the powheg tauola ttbar sample" << endl; ch_ttbar->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/TT_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
      
      cout << "Doing the Fall11 MC@NLO ttbar sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/TT_TuneZ2_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1/V04-02-29_fix_dilepton/merged*root");
    }

    //for samples with no negative weights
    //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil",doFRestimation, lumiToNormalizeTo*154./157.5, kttdil, false);
    //for mc@NLO
    baby->ScanChain(ch_ttbar, v_Cuts, "ttdil",doFRestimation, lumiToNormalizeTo*(154./157.5)*(190.41256/147.4), kttdil, false); //must mutliply by ratio of per-event xsec and PREP xsec to account for negative weights

    hist::color("ttdil", kGreen);
    //delete ch_ttbar;
  }
 if(runttdilsys) {
   TChain  *ch_ttbar= new TChain("Events");
   //  cout << "Doing the MadGraph ttbar sys sample" << endl; ch_ttbar->Add("/nfs-4/userdata/cms2/TTJets_TuneZ2_mass178_5_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v3/V04-02-29/merged*root");
   //    baby->ScanChain(ch_ttbar, v_Cuts, "ttdil_m178",doFRestimation, lumiToNormalizeTo*154./157.5, kttdil, false);
   //    hist::color("ttdil_m178", kGreen);
   //  cout << "Doing the MadGraph ttbar sys sample" << endl; ch_ttbar->Add("/nfs-4/userdata/cms2/TTJets_TuneZ2_mass166_5_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v3/V04-02-29/merged*root");
   //    baby->ScanChain(ch_ttbar, v_Cuts, "ttdil_m166",doFRestimation, lumiToNormalizeTo*154./157.5, kttdil, false);
   //    hist::color("ttdil_m166", kGreen);
   //   cout << "Doing the MadGraph ttbar sys sample" << endl; ch_ttbar->Add("/nfs-4/userdata/cms2/TTjets_TuneZ2_matchingup_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
   // cout << "Doing the MadGraph ttbar sys sample" << endl; ch_ttbar->Add("/nfs-4/userdata/cms2/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
   //cout << "Doing the MadGraph ttbar sys sample" << endl; ch_ttbar->Add("/hadoop/cms/store/group/snt/papers2011/Fall11MC/TTjets_TuneZ2_scaleup_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-29/merged*root");
   //cout << "Doing the MadGraph ttbar sys sample" << endl; ch_ttbar->Add("/nfs-6/userdata/cms2/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v2/V04-02-29/merged*root");

   //cout << "Doing the scale down mc@nlo ttbar sys FastSim sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/FastSim_mcatnlo-private-SC-fac-0.5-ren-0.5_LHE2EDM_v1_RWTH_0614_fhohle-TT_mcatnlo_private_fac-0.5-ren-0.5_withSC_FastSim_coherent_Summer11FullSim_FastSimJetFixTest_AOD_0616-98893edf89c918d9b6d63453d739e0c0/VB04-02-29_Fastsim/merged*root");
   //cout << "Doing the scale up mc@nlo ttbar sys FastSim sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/FastSim_mcatnlo-private-SC-fac-2-ren-2_LHE2EDM_v1_RWTH_0614_fhohle-TT_mcatnlo_private_fac-2-ren-2_withSC_FastSim_coherent_Summer11FullSim_FastSimJetFixTest_AOD_0616-98893edf89c918d9b6d63453d739e0c0/VB04-02-29_Fastsim/merged*root");

   cout << "Doing the mass down mc@nlo ttbar sys FastSim sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/TT-7TeV_MCatNLO_private_Mtop167.5_withSC_LHE_0827_fhohle-TTbar_MCatNLO_private_Mtop167.5_withSC_herwig_FastSim_coherent_Summer11FullSim_AOD_DESY_0827-98893edf89c918d9b6d63453d739e0c0/VB04-02-29_Fastsim/merged_ntuple*root");
   //cout << "Doing the mass up mc@nlo ttbar sys FastSim sample"   << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/TT-7TeV_MCatNLO_private_Mtop177.5_withSC_LHE_0827_fhohle-TTbar_MCatNLO_private_Mtop177.5_withSC_herwig_FastSim_coherent_Summer11FullSim_AOD_DESY_0827-98893edf89c918d9b6d63453d739e0c0/VB04-02-29_Fastsim/merged_ntuple*root");
   //cout << "Doing the nominal mc@nlo ttbar sys FastSim sample"   << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/TT-7TeV_MCatNLO-5217-SC_LHE_0822_fhohle-TTbar_mcatnlo-5217_herwig_FastSim_coherent_Summer11FullSim_AOD_DESY_0822-98893edf89c918d9b6d63453d739e0c0/VB04-02-29_Fastsim/merged_ntuple*root");


   //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil",doFRestimation, lumiToNormalizeTo*154./157.5, kttdil, false);
   
   //for mc@NLO
   //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil",doFRestimation, lumiToNormalizeTo*(154./157.5)*(190.41256/147.4), kttdil, false); //must mutliply by ratio of per-event xsec and PREP xsec to account for negative weights
   baby->ScanChain(ch_ttbar, v_Cuts, "ttdil",doFRestimation, lumiToNormalizeTo*(154./157.5)*(190.41256/147.4), kttdil, false, 1, 167.5); //must mutliply by ratio of per-event xsec and PREP xsec to account for negative weights
   //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil",doFRestimation, lumiToNormalizeTo*(154./157.5)*(190.41256/147.4), kttdil, false, 1, 177.5); //must mutliply by ratio of per-event xsec and PREP xsec to account for negative weights
   hist::color("ttdil", kGreen);
    
 }
 if(runSMS) {
   TChain  *ch_sms= new TChain("Events");
   if(runskim){
     ch_sms->Add(Form("%s/%s", cms2_skim_location.c_str(),"SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/skimmed*.root"));
   }
   else{
     cout << "Doing the fast sim SMS.. " << endl;
     // ch_sms->Add("/nfs-7/userdata/cms2/SMS-T2blnu_x-0p25to0p75_mStop-50to850_mLSP-50to800_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v2/VB04-02-29_Fastsim/merged*.root");
     ch_sms->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged*.root");
   }
   
   baby->ScanChain(ch_sms, v_Cuts,"sms", doFRestimation, lumiToNormalizeTo, ksms, false);
   hist::color("sms", kYellow);
   //delete ch_sms;                                                                                                                                                                   
 }
  
  if(runttotr) {
    TChain  *ch_ttor= new TChain("Events");
    if(runskim){
      cout << "Doing the Fall11 MC@NLO ttbar no-dileptons" << endl; ch_ttor->Add(Form("%s/%s", cms2_skim_location.c_str(),"TT_TuneZ2_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1/V04-02-29_fix_dilepton/skimmed*root"));
            
    }
    else{
      //cout << "Doing the MadGraph ttbar no-dileptons " << endl; ch_ttor->Add("/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
      cout << "Doing the Fall11 MC@NLO ttbar no-dileptons" << endl; ch_ttor->Add("/nfs-6/userdata/cms2/TT_TuneZ2_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1/V04-02-29_fix_dilepton/merged*root");
    }
    //for samples with no negative weights 
    //baby->ScanChain(ch_ttor, v_Cuts,"ttotr", doFRestimation, lumiToNormalizeTo*154./157.5, kttotr, false);
    //for MC@NLO
    baby->ScanChain(ch_ttor, v_Cuts,"ttotr", doFRestimation, lumiToNormalizeTo*(154./157.5)*(190.41256/147.4), kttotr, false); //must mutliply by ratio of per-event xsec and PREP xsec to account for negative weights
    hist::color("ttotr", kYellow);
    //delete ch_ttor;
  }
     
  if (runWjets) {
    TChain  *ch_wjets= new TChain("Events");
    cout << "Processing Wjets.."<<endl;
  
    if(runskim){
      ch_wjets->Add(Form("%s/%s", cms2_skim_location.c_str(),"WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/DileptonHyp/skimmed*root")); 
    }
    else{
      ch_wjets->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton/merged*root"); 
    }
    baby->ScanChain(ch_wjets,v_Cuts, "wjets", doFRestimation, lumiToNormalizeTo, kWjets, false);
    hist::color("wjets", kViolet);
  }
  if (runDYee) {
    TChain  *ch_dyee= new TChain("Events");
    cout << "Processing DY->ee" << endl;
    if(runskim){
      ch_dyee->Add(Form("%s/%s", cms2_skim_location.c_str(),"DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_dyee->Add(Form("%s/%s", cms2_skim_location.c_str(),"DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_dyee->Add(Form("%s/%s", cms2_skim_location.c_str(),"DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root")); 
    }
    else{
      ch_dyee->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
      ch_dyee->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
      ch_dyee->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
    }
    baby->ScanChain(ch_dyee, v_Cuts, "DYee",doFRestimation, lumiToNormalizeTo, kDYee, false);
    hist::color("DYee", kMagenta);
    //delete ch_dyee;
  }

  if (runDYmm) {
    TChain  *ch_dymm= new TChain("Events");
    cout << "Processing DY->mm" << endl;
    if(runskim){
      ch_dymm->Add(Form("%s/%s", cms2_skim_location.c_str(),"DYToMuMu_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_dymm->Add(Form("%s/%s", cms2_skim_location.c_str(),"DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_dymm->Add(Form("%s/%s", cms2_skim_location.c_str(),"DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root")); 
    }
    else{
      ch_dymm->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYToMuMu_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
      ch_dymm->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
      ch_dymm->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root"); 
    }
 
    baby->ScanChain(ch_dymm, v_Cuts, "DYmm",doFRestimation, lumiToNormalizeTo, kDYmm, false);
    hist::color("DYmm", kCyan);
    //delete ch_dymm;
  }
    
  if (runDYtautau) {
    TChain  *ch_dytt= new TChain("Events");
    cout << "Processing DY->tautau" << endl;
    if(runskim){
      ch_dytt->Add(Form("%s/%s", cms2_skim_location.c_str(),"DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2/V04-02-29/skimmed*.root"));
      ch_dytt->Add(Form("%s/%s", cms2_skim_location.c_str(),"DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_dytt->Add(Form("%s/%s", cms2_skim_location.c_str(),"DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
    }
    else{
      ch_dytt->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2/V04-02-29/merged*.root");

      ch_dytt->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
      ch_dytt->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
    }
    baby->ScanChain(ch_dytt,v_Cuts, "DYtautau",doFRestimation, lumiToNormalizeTo, kDYtautau, false);
    hist::color("DYtautau", kBlack);
    //delete ch_dytt;
  }
  if(runVV) {
    TChain  *ch_vv= new TChain("Events");
    cout << "Processing VV" << endl;
    if(runskim){
      ch_vv->Add(Form("%s/%s", cms2_skim_location.c_str(),"WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_vv->Add(Form("%s/%s", cms2_skim_location.c_str(),"WZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_vv->Add(Form("%s/%s", cms2_skim_location.c_str(),"WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_vv->Add(Form("%s/%s", cms2_skim_location.c_str(),"ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_vv->Add(Form("%s/%s", cms2_skim_location.c_str(),"ZZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/skimmed*.root"));
      ch_vv->Add(Form("%s/%s", cms2_skim_location.c_str(),"ZZJetsTo4L_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      /*
      ch_vv->Add(Form("%s/%s", cms2_skim_location.c_str(),"WWTo2L2Nu_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_vv->Add(Form("%s/%s", cms2_skim_location.c_str(),"ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_vv->Add(Form("%s/%s", cms2_skim_location.c_str(),"WZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      */
    }
    else{
      ch_vv->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");                                                     
      ch_vv->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/WZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");                                                     
      ch_vv->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");                                                     
      ch_vv->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root"); 
      ch_vv->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/ZZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root"); 
      ch_vv->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/ZZJetsTo4L_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root"); 
     /*
      ch_vv->Add("/nfs-7/userdata/cms2/WWTo2L2Nu_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
      ch_vv->Add("/nfs-7/userdata/cms2/ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
      ch_vv->Add("/nfs-7/userdata/cms2/WZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
      */
      }
    baby->ScanChain(ch_vv, v_Cuts, "VV" ,doFRestimation, lumiToNormalizeTo, kVV, false);
    hist::color("vv", 31);
    //delete ch_vv;
  }

  if(runtW) {
    TChain  *ch_tw= new TChain("Events");
    cout << "Processing tW" << endl;
    if(runskim){
      
      ch_tw->Add(Form("%s/%s", cms2_skim_location.c_str(),"T_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_tw->Add(Form("%s/%s", cms2_skim_location.c_str(),"Tbar_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_tw->Add(Form("%s/%s", cms2_skim_location.c_str(),"T_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_tw->Add(Form("%s/%s", cms2_skim_location.c_str(),"Tbar_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_tw->Add(Form("%s/%s", cms2_skim_location.c_str(),"T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      ch_tw->Add(Form("%s/%s", cms2_skim_location.c_str(),"Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/skimmed*.root"));
      /*
      ch_tw->Add(Form("%s/%s", "/nfs-4/userdata/yanjuntu/TPrimeSkim","TToBLNu_TuneZ2_t-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/skimmed*.root"));
      ch_tw->Add(Form("%s/%s", "/nfs-4/userdata/yanjuntu/TPrimeSkim","TToBLNu_TuneZ2_s-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/skimmed*.root"));
      ch_tw->Add(Form("%s/%s", "/nfs-4/userdata/yanjuntu/TPrimeSkim","TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/skimmed*.root")); 
      */
    }
    else{
      ch_tw->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/T_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1//V04-02-29/merged*.root");
      ch_tw->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
      ch_tw->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/T_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
      ch_tw->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
      ch_tw->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
      ch_tw->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*.root");
    
    }
    baby->ScanChain(ch_tw, v_Cuts, "tw", doFRestimation, lumiToNormalizeTo, ktW, false);
    hist::color("tw", kRed-3);
    //delete ch_tw;
  }
     

  if(runQCDPt15) {
    TChain  *ch_qcdpt15= new TChain("Events");
    cout << "Processing QCDPt15" << endl;
    ch_qcdpt15->Add("/nfs-3/userdata/cms2/QCD_Pt15_Spring10-START3X_V26_S09-v1/V03-04-13-07/diLepPt2010Skim/skimmed_ntuple.root");
    baby->ScanChain(ch_qcdpt15, v_Cuts, "qcd15", doFRestimation, lumiToNormalizeTo, kQCDPt15, false);
    hist::color("qcd15", 28);
    // delete ch_qcdpt15;
  }
     

  if(runQCDPt30) {
    TChain  *ch_qcdpt30= new TChain("Events");
    cout << "Processing QCDPt30" << endl;
    ch_qcdpt30->Add("/nfs-3/userdata/cms2/QCD_Pt30_Spring10-START3X_V26_S09-v1/V03-04-13-07/diLepPt2010Skim/skimmed_ntuple.root");
    baby->ScanChain(ch_qcdpt30, v_Cuts, "qcd30", doFRestimation, lumiToNormalizeTo, kQCDPt30, false);
    hist::color("qcd30", 29);
    //delete ch_qcdpt30;
  }


  if(runVqq) {
    TChain  *ch_vqq             = new TChain("Events");
    cout << "Processing Vqq" << endl;
    ch_vqq->Add("/nfs-3/userdata/cms2/VqqJets-madgraph_Spring10-START3X_V26_S09-v1/merged*.root" );
    baby->ScanChain(ch_vqq, v_Cuts, "vqq", doFRestimation, lumiToNormalizeTo, kVqq, false);
    //delete ch_vqq;
  }
      
     
  if(runphotonJet15) {

    cout << "USING OLD PHOTON JET SAMPLE" << endl;
    TChain  *ch_photonjet15= new TChain("Events"); 
    cout << "Processing PhotonJet15" << endl;
    ch_photonjet15->Add("/nfs-3/userdata/cms2/PhotonJet_Pt15_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*.root");
    baby->ScanChain(ch_photonjet15, v_Cuts, "photonjet15", doFRestimation, lumiToNormalizeTo, 1, false);
    hist::color("photonjet15", 30);
    //delete ch_photonjet15;
  }


  if(runphotonVJets) {
      
    cout << "Using Photon VJets sample" << endl;
    TChain *ch_photonVJets = new TChain("Events");
    cout << "Processing PhotonVJets" << endl;
    ch_photonVJets->Add("/nfs-3/userdata/cms2/PhotonVJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*.root");
    baby->ScanChain(ch_photonVJets, v_Cuts, "photonVJets", doFRestimation, lumiToNormalizeTo, 1, false);
    //delete ch_photonVJets;

  }
  TString cutstring = "";
  for(unsigned int i = 0; i < v_Cuts.size(); i++) 
    cutstring   = cutstring + "_" + v_Cuts.at(i);

  if(doFRestimation) 
    cutstring = outputDir + "/FRhist" + cutstring + ".root";
  else  
    cutstring   = outputDir + "/hist" + cutstring + ".root";    
    if( cutstring.Contains("applyNoCuts") ) cutstring  = outputDir + "/hist_noCuts.root";
    
  cout << "Saving histograms to: " << cutstring << endl;
  hist::saveHist(cutstring.Data());
  cout << "done saving" << endl;

  hist::deleteHistos();
  


  }


}
