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
char* outFile = "myHist_WW_fast.root";

// Flags for files to run over
bool runWW    = true;
bool runWZ    = true;
bool runZZ    = true;
bool runWjets = true;
bool runDYee  = true;
bool runDYmm  = true;
bool runDYtt  = true;
bool runttbar = true;

// Load various tools
gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

gROOT->ProcessLine(".x setup.C");

// Load and compile the looping code
gROOT->ProcessLine(".L ClaudioLoopingFunctionFast.C+");

//WW file
if (runWW) {
     TChain *fWW = new TChain("Events");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_1.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_2.root");
     fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_3.root");
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
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_1.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_2.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_3.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_4.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_5.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_6.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_7.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_1.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_2.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_3.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_4.root");
     fWjets->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_5.root");
}

//DYee file
if (runDYee) {
     TChain *fDYee = new TChain("Events");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_1.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_10.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_11.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_12.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_13.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_14.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_15.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_16.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_17.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_18.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_19.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_2.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_20.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_21.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_22.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_23.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_24.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_25.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_26.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_27.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_28.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_29.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_3.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_30.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_31.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_32.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_33.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_34.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_35.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_36.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_37.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_38.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_39.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_4.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_40.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_41.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_42.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_43.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_44.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_45.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_46.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_47.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_48.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_49.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_5.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_50.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_51.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_52.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_53.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_54.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_55.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_56.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_57.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_58.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_59.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_6.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_60.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_61.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_62.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_63.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_64.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_65.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_66.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_67.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_68.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_69.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_7.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_8.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_9.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_1.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_10.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_11.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_12.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_13.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_14.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_15.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_16.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_17.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_18.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_19.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_2.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_20.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_21.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_22.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_23.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_24.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_25.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_26.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_27.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_28.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_29.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_3.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_30.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_31.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_32.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_33.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_34.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_35.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_36.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_37.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_38.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_39.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_4.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_40.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_41.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_42.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_43.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_44.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_45.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_46.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_47.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_48.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_49.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_5.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_50.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_51.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_52.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_53.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_54.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_55.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_56.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_57.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_6.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_7.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_8.root");
     fDYee->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_9.root");
}

//DYmm file
if (runDYmm) {
  TChain *fDYmm = new TChain("Events");
//   fDYmm->Add("data/electron_soup/ntuple*.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_1.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_10.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_11.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_12.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_13.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_14.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_15.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_16.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_17.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_18.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_19.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_2.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_20.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_21.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_22.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_23.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_24.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_25.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_26.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_27.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_28.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_29.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_3.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_30.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_31.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_32.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_33.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_34.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_35.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_36.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_37.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_38.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_39.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_4.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_40.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_41.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_42.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_43.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_44.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_45.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_46.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_47.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_48.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_49.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_5.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_50.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_51.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_52.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_53.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_54.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_55.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_56.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_57.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_58.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_59.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_6.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_60.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_61.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_62.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_63.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_64.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_65.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_66.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_67.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_68.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_69.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_7.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_8.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_9.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_1.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_10.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_11.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_12.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_13.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_14.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_15.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_16.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_17.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_18.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_19.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_2.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_20.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_21.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_22.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_23.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_24.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_25.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_26.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_27.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_28.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_29.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_3.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_30.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_31.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_32.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_33.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_34.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_35.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_36.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_37.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_38.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_39.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_4.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_40.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_41.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_42.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_43.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_44.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_45.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_46.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_47.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_48.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_49.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_5.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_50.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_51.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_52.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_53.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_54.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_55.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_56.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_57.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_6.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_7.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_8.root");
     fDYmm->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_9.root");
}

//DYtt file
if (runDYtt) {
  TChain *fDYtt = new TChain("Events");
//   fDYtt->Add("data/electron_soup/ntuple*.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_1.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_10.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_11.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_12.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_13.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_14.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_15.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_16.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_17.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_18.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_19.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_2.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_20.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_21.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_22.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_23.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_24.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_25.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_26.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_27.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_28.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_29.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_3.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_30.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_31.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_32.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_33.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_34.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_35.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_36.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_37.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_38.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_39.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_4.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_40.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_41.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_42.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_43.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_44.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_45.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_46.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_47.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_48.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_49.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_5.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_50.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_51.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_52.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_53.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_54.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_55.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_56.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_57.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_58.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_59.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_6.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_60.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_61.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_62.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_63.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_64.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_65.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_66.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_67.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_68.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_69.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_7.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_8.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_9.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_1.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_10.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_11.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_12.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_13.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_14.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_15.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_16.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_17.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_18.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_19.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_2.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_20.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_21.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_22.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_23.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_24.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_25.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_26.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_27.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_28.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_29.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_3.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_30.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_31.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_32.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_33.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_34.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_35.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_36.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_37.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_38.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_39.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_4.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_40.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_41.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_42.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_43.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_44.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_45.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_46.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_47.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_48.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_49.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_5.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_50.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_51.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_52.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_53.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_54.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_55.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_56.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_57.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_6.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_7.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_8.root");
     fDYtt->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_electron_soup_postprocessed_b93843bed05e4e5b456b0e788b5de036/ntuple_soup_9.root");
}

//ttbar file
if (runttbar) {
  TChain *fttbar = new TChain("Events");
//   fttbar->Add("data/electron_soup/ntuple*.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_1.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_10.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_11.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_12.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_13.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_14.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_15.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_16.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_17.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_18.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_19.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_2.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_20.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_21.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_22.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_23.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_24.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_25.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_26.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_27.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_28.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_29.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_3.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_30.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_31.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_32.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_33.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_34.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_35.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_36.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_37.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_38.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_39.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_4.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_40.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_41.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_42.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_43.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_44.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_45.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_46.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_47.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_48.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_49.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_5.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_50.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_51.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_52.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_53.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_54.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_55.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_56.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_57.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_58.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_59.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_6.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_60.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_61.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_62.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_63.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_64.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_65.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_66.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_67.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_68.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_69.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_7.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_8.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_muon_soup_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_soup_9.root");
     fttbar->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/jmuelmen/cms2_electron_soup_postprocessed_split_ttbar_a8faccabfc12f0755896310bdda19928/ntuple_ttbar_1.root");
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

//save all the histograms
saveHist(outFile);
}
