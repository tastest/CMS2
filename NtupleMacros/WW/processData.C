//==============================================================
//
// This runs over the skimmed files (see AAREADME.txt)
//
// To run on unskimmed files, change the file names, but be
// careful about Drell Yan.
//
// The DY skims have separated out the three generated final states,
// while the unskimmed file has them all together, so it is a bit 
// more complicated.  You should:
// (1) Uncomment the block following the "Full Drell Yan file" comment
// (2) Optionally comment out blocks where the skimmed DY files are opened
// (3) Replace these three statements
//        ScanTree(tDYtautau,"DYtautau", -1, 1.2);
//        ScanTree(tDYmm,"DYmm", -1, 1.2);
//        ScanTree(tDYee,"DYee", -1, 1.2);
//     by
//        ScanTree(tDY,"DYtautau", 2, 1.2);
//        ScanTree(tDY,"DYmm",     1, 1.2);
//        ScanTree(tDY,"DYee",     0, 1.2);
//     Note the change in the 3rd calling parameter!
//
//==============================================================
{
  using namespace std;
// Output file
  const char* outFile = "processed_data.root";
  const bool identifyVVEvents = false; // careful with this option. You don't need it for not mixed samples
  const bool identifyDYEvents = false;

// Flags for files to run over
  bool runWW    = true;
  bool runWZ    = true;
  bool runZZ    = false;
  bool runWjets = false;
  bool runDYee  = false;
  bool runDYmm  = false;
  bool runDYtt  = false;
  bool runttbar = false;
  bool runtW    = true;
  bool runQCD   = false; 

// Load various tools
 gROOT->ProcessLine(".x init.C");
 gSystem->Load("libCMS2NtupleMacrosCORE");
 gSystem->Load("libRooFit.so");
 gSystem->Load("libCMS2NtupleMacrosLooper");
 gROOT->LoadMacro("../Tools/getMyHistosNames.C");
 gROOT->LoadMacro("../Tools/histtools.C+");
 gROOT->LoadMacro("../Tools/browseStacks.C");
 
 // Define colors numbers:
 gStyle->SetPalette(1);
 enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };

 RooDataSet *fullDataSet(0);

// read dataset prefix
 string dataset = "data";
 
 if (runWW)
   ProcessSample(dataset+"/WW_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple.root", WW, 1.65, fullDataSet, kRed);

 if (runWZ)
   ProcessSample(dataset+"/WZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple.root", WZ, 1.84, fullDataSet, kBlue);

 if (runZZ)
   ProcessSample(dataset+"/ZZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple.root", ZZ, 1.47, fullDataSet, kGreen);
 
 if (runWjets)
   ProcessSample(dataset+"/WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root", Wjets, 1.0, fullDataSet, 40);

 if (runDYee)
   ProcessSample(dataset+"/Zee_Summer09-MC_31X_V3_7TeV_TrackingParticles-v1/V03-00-35/merged_ntuple*.root", DYee, 1.14, fullDataSet, kMagenta, identifyDYEvents);
 
 if (runDYmm)
   ProcessSample(dataset+"/Zmumu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root", DYmm, 1.14, fullDataSet, kCyan, identifyDYEvents);
 
 if (runDYtt)
   ProcessSample(dataset+"/Ztautau_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root", DYtt, 1.14, fullDataSet, kBlack, identifyDYEvents);

 if (runttbar)
   ProcessSample(dataset+"/TTbarJets-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged_ntuple*.root", ttbar, 1.0, fullDataSet, kYellow);
 
 if (runtW)
   ProcessSample(dataset+"/SingleTop_tWChannel-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged_ntuple*.root", tW, 1.0, fullDataSet, 63);

 //QCD file
TChain *fqcd = new TChain("Events");
if (runQCD) {
  fqcd->Add((dataset+"/cms2-V01-02-06/InclusiveMuPt15/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/InclusiveMu5Pt50/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/QCD_EMenriched_Pt20to30/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/QCD_EMenriched_Pt30to80/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/QCD_EMenriched_Pt80to170/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/QCD_BCtoE_Pt20to30/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/QCD_BCtoE_Pt30to80/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/QCD_BCtoE_Pt80to170/merged_ntuple*.root").c_str());
}


if (runQCD) {
  cout << "Processing QCD.."<<endl;
  RooDataSet* data = ScanChain(fWjets, qcd, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("qcd", 40);
}

//save all the histograms
//saveHist(outFile);
 cout << "got up to here" << endl;
 TList* list = gDirectory->GetList() ;
 TIterator* iter = list->MakeIterator();
 
 TRegexp re("*",kTRUE) ;
 
 TFile outf(outFile,"RECREATE") ;
 while(obj=iter->Next()) {
   if (TString(obj->GetName()).Index(re)>=0) {
     obj->Write() ;
     cout << "." ;
     cout.flush() ;
   }
 }
 cout << endl ;
 if (fullDataSet) {
   std::string description("Full N-1 dataset:");
   if (runWW)     description+=" WW";
   if (runWZ)     description+=" WZ";
   if (runZZ)     description+=" ZZ";
   if (runWjets)  description+=" Wjets";
   if (runDYee)   description+=" DYee";
   if (runDYmm)   description+=" DYmm";
   if (runDYtt)   description+=" DYtt";
   if (runttbar)  description+=" ttbar";
   if (runtW)     description+=" tW";
   if (runQCD)    description+=" QCD";
   fullDataSet->SetName("fulldataset");
   fullDataSet->SetTitle(description.c_str());
   fullDataSet->Write();
 }
 outf.Close() ;
 
 delete iter ;

}
