#include <fstream>

void doTable() {
  // Load and compile something to allow proper treatment of vectors
  // Not clear that it is needed
  // gROOT->LoadMacro("loader.C+");

  // Load various tools
  // gROOT->ProcessLine(".x setup.C");

  // Load and compile the looping code
  gROOT->ProcessLine(".L ../Tools/tableUtilities.C");

  // histogram name
  const char* histname = "hNjetsBothLeptonsVeto";
  unsigned int nBins = 5;

  cout << "oli" << endl;

  ofstream file(Form("tables_%s.txt",histname));

  file << " NjetsBothLeptonsVeto " << endl;
  file << printTrilepTable("WZ",Form("wz_%s",histname),nBins);
  file << printTrilepTable("DYee",Form("DYee_%s",histname),nBins);
  file << printTrilepTable("DYmm",Form("DYmm_%s",histname),nBins);
  file << printTrilepTable("DYtautau",Form("DYtautau_%s",histname),nBins);
  file << printTrilepTable("TT",Form("tt_%s",histname),nBins);
  file << printTrilepTable("WJets",Form("wjets_%s",histname),nBins);
  file << printTrilepTable("WW",Form("ww_%s",histname),nBins);
  file << printTrilepTable("ZZ",Form("zz_%s",histname),nBins);
  file << printTrilepTable("tW",Form("tW_%s",histname),nBins);
   
  file.close();
}
