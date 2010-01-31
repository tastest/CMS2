 {
   //copied from setup.C
	gSystem->Load("libTree.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libEG.so");
	gSystem->Load("libMathCore.so");
	//gSystem->Load("libCMS2NtupleMacrosCORE.so");
	//gSystem->Load("libCMS2NtupleMacrosTools.so");
	//gSystem->Load("libCMS2NtupleMacrosLooper.so");



   gROOT->ProcessLine(".L plotscript.C++");

   gStyle->SetOptStat(0);
   gStyle->SetPadRightMargin ( 0.1);
   gStyle->SetLabelSize  ( 0.040,"X");
   gStyle->SetLabelSize  ( 0.040,"Y");
   gStyle->SetTitleSize  ( 0.040,"X");//was 0.050
   gStyle->SetTitleOffset( 0.600,"X");//was 1.2

   plotscript("ExclSpike_Runs_allRuns_withTrackingCuts_withTriggers_jan23_6.root",
			  "ExclSpike_Runs_allRuns_withTrackingCuts_withTriggers_mc_2.root",
			  "ExclSpike_Runs_allRuns_withTrackingCuts_withTriggers_jan23_2tev_3.root");

 }
