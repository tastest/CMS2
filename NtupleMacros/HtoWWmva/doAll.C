{
  gROOT->ProcessLine(".L CORE/CMS2.cc+");
  gROOT->ProcessLine(".L CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L CORE/metSelections.cc+");
  gROOT->ProcessLine(".L CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L CORE/muonSelections.cc+");
  gROOT->ProcessLine(".L CORE/jetSelections.cc+");
  gROOT->ProcessLine(".L CORE/mcSelections.cc+");
  gROOT->ProcessLine(".L histtools.C+");
  gROOT->ProcessLine(".L runLooper.C+");
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");

  //choose samples to run over----------------------
/*************** BACKGROUND SAMPLES ***************/
  runLooper("TTJets_PU");
  runLooper("WJetsToLNu_PU");
  runLooper("WWTo2L2Nu_PU");
  runLooper("tW_PU");
  runLooper("ZZ_PU");
  runLooper("WZ_PU");
  runLooper("GluGluToWWTo4L_PU");
  runLooper("VVJetsTo4L_PU");
  runLooper("DYToEEM10To20_PU");
  runLooper("DYToEEM20_PU");
  runLooper("DYToMuMuM10To20_PU");
  runLooper("DYToMuMuM20_PU");
  runLooper("DYToTauTauM10To20_PU");
  runLooper("DYToTauTauM20_PU");
  /*************** SIGNAL SAMPLES (E+MU) ***************/
  runLooper("HToWWTo2L2NuM130_PU");
  runLooper("HToWWTo2L2NuM160_PU");
  runLooper("HToWWTo2L2NuM200_PU");
/*************** SIGNAL SAMPLES (TAUS) ***************/
  runLooper("HToWWTo2Tau2NuM130_PU");
  runLooper("HToWWToLNuTauNuM130_PU");
  runLooper("HToWWTo2Tau2NuM160_PU");
  runLooper("HToWWToLNuTauNuM160_PU");
  runLooper("HToWWTo2Tau2NuM200_PU");
  runLooper("HToWWToLNuTauNuM200_PU");


/*************** SAMPLES NO PU ***************/
//   runLooper("HToWWTo2L2NuM130");
//   runLooper("HToWWTo2L2NuM160");
//   runLooper("HToWWTo2L2NuM200");
//   runLooper("TTJets");
//   runLooper("WJetsToLNu");
//   runLooper("WWTo2L2Nu");

}
