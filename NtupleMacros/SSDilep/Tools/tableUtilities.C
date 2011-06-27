#include <sstream>
#include <iomanip>

const char* printTrilepTable(const char* mode, const char* pattern, const int numBins) {

  const int allBuckets = 20;
  char *suffix[allBuckets+1];
  suffix[0]  = "epepep";
  suffix[1]  = "epepem";
  suffix[2]  = "epemem";
  suffix[3]  = "ememem";
  suffix[4]  = "mpmpmp";
  suffix[5]  = "mpmpmm";
  suffix[6]  = "mpmmmm";
  suffix[7]  = "mmmmmm";
  suffix[8]  = "mpepep";
  suffix[9]  = "mmepep";
  suffix[10] = "mmepem";
  suffix[11] = "mpepem";
  suffix[12] = "mmemem";
  suffix[13] = "mpemem";
  suffix[14] = "mpmpep";
  suffix[15] = "mpmpem";
  suffix[16] = "mpmmem";
  suffix[17] = "mpmmep";
  suffix[18] = "mmmmem";
  suffix[19] = "mmmmep";
  suffix[20] = "all";

   TH1* histos[allBuckets+1];
   
   for ( int i = 0; i < allBuckets+1; ++i ) {
      TH1* hist = (TH1*)gFile->Get(Form("%s_%s",pattern,suffix[i]));
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
