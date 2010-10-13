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
// Output file
char* outFile = "myHist_caloIso_new.root";

// Flags for files to run over
bool runWW    = true;
bool runWZ    = false;
bool runZZ    = false;
bool runWjets = true;
bool runDYee  = false;
bool runDYmm  = false;
bool runDYtt  = false;
bool runttbar = false;
bool runtW    = false;

// Load various tools
gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

gROOT->ProcessLine(".x setup.C");

// Load and compile the looping code
gROOT->ProcessLine(".L caloIsoLooper.C+");

//WW file
if (runWW) {
     TChain *fWW = new TChain("Events");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_1.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_10.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_11.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_12.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_13.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_14.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_15.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_16.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_17.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_18.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_19.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_2.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_20.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_21.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_22.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_23.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_24.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_25.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_26.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_27.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_28.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_29.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_3.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_30.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_31.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_32.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_33.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_34.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_35.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_36.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_37.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_38.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_39.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_4.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_40.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_41.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_42.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_43.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_44.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_45.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_46.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_47.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_48.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_49.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_5.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_50.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_51.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_52.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_53.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_54.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_55.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_56.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_57.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_58.root");
//      fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_59.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_6.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_60.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_61.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_62.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_63.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_64.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_65.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_66.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_67.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_68.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_69.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_7.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_70.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_71.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_72.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_73.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_74.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_75.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_76.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_8.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_CMSSW_1_6_7-CSA07-1198096665_080901_d2c33c2823e290f7aa24569a0c68265e/ntuple_signal_9.root");
}

//WZ file
if (runWZ) {
     TChain *fWZ = new TChain("Events");
     fWZ->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WZ_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_1.root");
}

//ZZ file
if (runZZ) {
     TChain *fZZ = new TChain("Events");
     fZZ->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_ZZ_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_1.root");
}

//Wjets file
if (runWjets) {
     TChain *fWjets = new TChain("Events");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_1.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_10.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_11.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_12.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_13.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_14.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_15.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_16.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_17.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_18.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_19.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_2.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_20.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_21.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_22.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_23.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_24.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_25.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_26.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_27.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_28.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_29.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_3.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_30.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_31.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_32.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_33.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_34.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_35.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_36.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_37.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_38.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_39.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_4.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_40.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_41.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_42.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_43.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_44.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_45.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_46.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_47.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_48.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_49.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_5.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_50.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_51.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_52.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_53.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_54.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_55.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_56.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_57.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_58.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_59.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_6.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_60.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_61.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_62.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_63.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_64.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_65.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_66.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_67.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_68.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_69.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_7.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_70.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_71.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_72.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_73.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_74.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_75.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_76.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_8.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_CMSSW_1_6_7-CSA07-Chowder-P1-PDMuon-Skims4-topDiLeptonMuonX_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_9.root");

     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_1.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_10.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_11.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_12.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_13.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_14.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_15.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_16.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_17.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_18.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_19.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_2.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_20.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_21.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_22.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_23.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_24.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_25.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_26.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_27.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_28.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_29.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_3.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_30.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_31.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_32.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_33.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_34.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_35.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_36.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_37.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_38.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_39.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_4.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_40.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_41.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_42.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_43.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_44.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_45.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_46.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_47.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_48.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_49.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_5.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_50.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_51.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_52.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_53.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_54.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_55.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_56.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_57.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_58.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_59.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_6.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_60.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_61.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_62.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_63.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_64.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_65.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_66.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_7.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_8.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_CMSSW_1_6_7-CSA07-Chowder-B1-PDElectron-Skims6-topDiLepton2Electron_080901_5a13a81119dfebc5ee781d4849bcdb7f/ntuple_soup_9.root");
}

//DYee file
if (runDYee) {
     TChain *fDYee = new TChain("Events");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_1.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_2.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_3.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_4.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_5.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_6.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_1.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_2.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_3.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_4.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_5.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_6.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_7.root");
}

//DYmm file
if (runDYmm) {
     TChain *fDYmm = new TChain("Events");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_1.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_2.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_3.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_4.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_5.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_6.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_1.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_2.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_3.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_4.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_5.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_6.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_7.root");
}

//DYtt file
if (runDYtt) {
     TChain *fDYtt = new TChain("Events");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_1.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_2.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_3.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_4.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_5.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_6.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_1.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_2.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_3.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_4.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_5.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_6.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_7.root");
}

//ttbar file
if (runttbar) {
     TChain *fttbar = new TChain("Events");
//   fttbar->Add("data/electron_soup/ntuple*.root");
     fttbar->Add("/data/tmp/jmuelmen/ntuple_ttbar_1.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_ttbar_a8faccabfc12f0755896310bdda19928/ntuple_ttbar_1.root");
}

if (runtW) {
     TChain *ftW = new TChain("Events");
//   ftW->Add("data/electron_soup/ntuple*.root");
     ftW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_tW_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_1.root");
}

// Define colors numbers:
gStyle->SetPalette(1);
enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };

// Process files one at a time, and color them as needed
if (runWW) {
  cout << "Processing WW.."<< endl;
  ScanChain(fWW, WW);
  hist::color("ww", kRed);
}

if (runWZ) {
  cout << "Processing WZ.."<< endl;
  ScanChain(fWZ, WZ);
  hist::color("wz", kBlue);
}

if (runZZ) {
  cout << "Processing ZZ.."<< endl;
  ScanChain(fZZ, ZZ);
  hist::color("zz", kGreen);
}

if (runWjets) {
  cout << "Processing Wjets.."<<endl;
  ScanChain(fWjets, Wjets);
  hist::color("wjets", 40);
}

if (runDYee) {
  cout << "Processing DYee.."<<endl;
  ScanChain(fDYee, DYee);
  hist::color("dyee", kMagenta);
}

if (runDYmm) {
  cout << "Processing DYmm.."<<endl;
  ScanChain(fDYmm, DYmm);
  hist::color("dymm", kCyan);
}

if (runDYtt) {
  cout << "Processing DYtt.."<<endl;
  ScanChain(fDYtt, DYtt);
  hist::color("dytt", kBlack);
}

if (runttbar) {
  cout << "Processing ttbar.."<<endl;
  ScanChain(fttbar, ttbar);
  hist::color("ttbar", kYellow);
}

if (runtW) {
  cout << "Processing tW.."<<endl;
  ScanChain(ftW, tW);
  hist::color("ttbar", 63);
}

//save all the histograms
saveHist(outFile);
}
