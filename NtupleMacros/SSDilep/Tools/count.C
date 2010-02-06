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
  TChain *fE = new TChain("Events");
  fE->Add("/uscms_data/d1/fgolf/sntuples/soups/topDilepton2Electron/ntuple*.root");
  int nE = fE->GetEntries();
  TChain *fM = new TChain("Events");
  fM->Add("/uscms_data/d1/fgolf/sntuples/soups/topDileptonMuonX/ntuple*.root");
  int nM = fM->GetEntries();
  std::cout << nE << " topDilepton2Electron events were ntuplized." << endl;
  std::cout << nM << " topDileptonMuonX events were ntuplized." << endl;
}
