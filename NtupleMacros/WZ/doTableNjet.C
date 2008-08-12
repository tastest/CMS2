#include <iomanip>
#include <sstream>
#include <fstream>

const char* printTable(const char* mode, const char* pattern, const int numBins) {

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

  TH1* histos[allBuckets+1];
   
  for ( int i = 0; i < allBuckets+1; ++i ) {
    TH1* hist = (TH1*)gFile->Get(Form("%s_%s",pattern,suffix[i]));
    // cout << "Reading histogram: " << Form("%s_%s",pattern,suffix[i]) << endl;
    if ( hist != 0 ) {
      histos[i] = hist;
    } else {
      cout << "Histogram: " << Form("%s_%s",pattern,suffix[i]) << "could not be found in input file" << endl;
      return;
    }
  }

  ostringstream stream;
  stream.setf(ios::fixed);

  stream << "| *Sample* | *Njets<0* | *Njets=0* | *Njets=1* | *Njets=2* | *Njets=3* | *Njets=4* | *Njets>4* | " << endl; 

  int bins = histos[1]->GetNbinsX();
  for (int i=0; i<allBuckets+1; i++){
    stream << "| " << mode << " " << suffix[i] ;
    for (int j=0; j<=numBins; j++){
      if ( histos[i]->GetBinContent(j) > 0 ) {
	stream << " |  %RED%" << setprecision(3) << histos[i]->GetBinContent(j) << "%ENDCOLOR%";
      } else {
	stream << " |  " << setprecision(3) << histos[i]->GetBinContent(j) ;
      }
    }
    if ( histos[i]->GetBinContent(numBins+1,bins+1) > 0 ) {
      stream << " |  %RED%" << setprecision(3) << histos[i]->GetBinContent(numBins+1,bins+1) << "%ENDCOLOR% | " << endl;
    } else {
      stream << " |  " << setprecision(3) << histos[i]->GetBinContent(numBins+1,bins+1) << " | " << endl;
    }
  }
  stream << endl;
   
  cout << stream.str();
   
  return stream.str().c_str();

}

void doTableNjet() {
  // Load and compile something to allow proper treatment of vectors
  // Not clear that it is needed
  //gROOT->LoadMacro("loader.C+");

  // Load various tools
  //gROOT->ProcessLine(".x setup.C");

  ofstream file("tables_final.txt");

  file << " NjetsFinal " << endl;
  file << printTable("WZ","wz_hNjetsFinal",5);
//   file << printTable("DYee","DYee_hNjetsFinal",5);
//   file << printTable("DYmm","DYmm_hNjetsFinal",5);
//   file << printTable("DYtautau","DYtautau_hNjetsFinal",5);
//   file << printTable("tt","tt_hNjetsFinal",5);
//   file << printTable("wjets","wjets_hNjetsFinal",5);
//   file << printTable("ww","ww_hNjetsFinal",5);
//   file << printTable("zz","zz_hNjetsFinal",5);
   
  file.close();
}

