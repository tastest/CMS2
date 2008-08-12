#include <iomanip>
#include <sstream>
#include <fstream>

const char* printTable(const char* pattern) {

  const unsigned int allBuckets = 20;
  char *suffix[allBuckets+1];
  suffix[0]  = "mpmpmp";
  suffix[1]  = "mpmpmm";
  suffix[2]  = "mpmmmm";
  suffix[3]  = "mmmmmm";
  suffix[4]  = "mpmpep";
  suffix[5]  = "mpmpem";
  suffix[6]  = "mpmmep";
  suffix[7]  = "mpmmem";
  suffix[8]  = "mmmmep";
  suffix[9]  = "mmmmem";
  suffix[10] = "mpepep";
  suffix[11] = "mpepem";
  suffix[12] = "mpemem";
  suffix[13] = "mmepep";
  suffix[14] = "mmepem";
  suffix[15] = "mmemem";
  suffix[16] = "epepep";
  suffix[17] = "epepem";
  suffix[18] = "epemem";
  suffix[19] = "ememem";
  suffix[20] = "all";

  float wz_5[num_suffix];
  float wz_10[num_suffix];
  float wz_20[num_suffix];
   
  TH1* histos[num_suffix];
   
  for ( int i = 0; i < num_suffix; ++i ) {
    TH1* hist = (TH1*)gFile->Get(Form("%s_%s",pattern,suffix[i]));
    // cout << "Reading histogram: " << Form("%s_%s",pattern,suffix[i]) << endl;
    if ( hist != 0 ) {
      histos[i] = hist;
    } else {
      cout << "Histogram: " << Form("%s_%s",pattern,suffix[i]) << "could not be found in input file" << endl;
      return;
    }
  }

  for ( int i = 0; i < num_suffix; ++i ) {
    wz_5[i]  = histos[i]->GetBinContent(2) + histos[i]->GetBinContent(3) + histos[i]->GetBinContent(4);
    wz_10[i] = histos[i]->GetBinContent(3) + histos[i]->GetBinContent(4);
    wz_20[i] = histos[i]->GetBinContent(4);
  }
   
  ostringstream stream;

  stream << "| " << pattern << " | > 5 GeV | > 10 GeV | > 20 GeV |" << endl;
  for (int i=0; i<num_suffix; i++){
    stream << "| " << suffix[i] << " | " << wz_5[i] << " | " << wz_10[i] << " | " << wz_20[i] << " | " << endl;
  }
  stream << endl;
   
  cout << stream.str();
   
  return stream.str().c_str();

}

void doTable() {
  // Load and compile something to allow proper treatment of vectors
  // Not clear that it is needed
  gROOT->LoadMacro("loader.C+");

  // Load various tools
  gROOT->ProcessLine(".x setup.C");

  ofstream file("tables.txt");

  file << printTable("wz_hLowestPt");
  file << printTable("wz_hLowestPtGood");
  file << printTable("wz_hLowestPtTrue");
  file << printTable("wz_hMET");
  file << printTable("wz_hMETTrue");
  file << printTable("wz_hMETGood");
  file << printTable("wz_hMETAll");
  file << printTable("wz_hMETAllGood");
  file << printTable("wz_hMETAllTrue");

  file << printTable("DY_hLowestPt");
  file << printTable("DY_hLowestPtGood");
  file << printTable("DY_hLowestPtTrue");
  file << printTable("DY_hMET");
  file << printTable("DY_hMETTrue");
  file << printTable("DY_hMETGood");
  file << printTable("DY_hMETAll");
  file << printTable("DY_hMETAllGood");
  file << printTable("DY_hMETAllTrue");
   
  file.close();
}
