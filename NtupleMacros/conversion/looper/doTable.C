#include <iostream>
#include <iomanip>


void doTable(){
gROOT->ProcessLine(".x setup.C");
TFile *file = TFile::Open("myHist.root");
 float electron[4], gamma[4];
 
 for(int i=0; i<4; i++){
   electron[i] = 0.0;
   gamma [i] = 0.0;
 }
 char* finalState[4];
 finalState[0] = "analysis cut";
 finalState[1] = "d0 cut";
 finalState[2] = "cutting on number of missing layers";
 finalState[3] = "cutting on dist and dcot";

 //for channel 1
 electron[0]     = Electron_elsEta_Ch0H1->Integral(0,  Electron_elsEta_Ch0H1->GetNbinsX()+1);
 gamma[0] = Gamma_elsEta_Ch0H1->Integral(0,  Gamma_elsEta_Ch0H1->GetNbinsX()+1);

 //for channel 2
 
 electron[1]     = Electron_elsEta_Ch0H3->Integral(0,  Electron_elsEta_Ch0H1->GetNbinsX()+1);
 gamma[1] = Gamma_elsEta_Ch0H3->Integral(0,  Gamma_elsEta_Ch0H1->GetNbinsX()+1);
 
 electron[2]     = Electron_elsEta_Ch0H4->Integral(0,  Electron_elsEta_Ch0H1->GetNbinsX()+1);
 gamma[2] = Gamma_elsEta_Ch0H4->Integral(0,  Gamma_elsEta_Ch0H1->GetNbinsX()+1);
 
 electron[3]     = Electron_elsEta_Ch0H5->Integral(0,  Electron_elsEta_Ch0H1->GetNbinsX()+1);
 gamma[3] = Gamma_elsEta_Ch0H5->Integral(0,  Gamma_elsEta_Ch0H1->GetNbinsX()+1);
 
 cout << "| | electron|gamma" <<"|"<< endl;

for (int i=0; i<4; i++){
  
  cout << setiosflags(ios::fixed)<< "| " << finalState[i] << " | " <<electron[i] <<"|"<<gamma[i] <<"|"<< endl;
  cout << setiosflags(ios::fixed)<< "| " << "efficiency" << " | " <<electron[i]/electron[0]  <<"|"<<gamma[i]/gamma[0]  <<"|"<< endl;
}

}
