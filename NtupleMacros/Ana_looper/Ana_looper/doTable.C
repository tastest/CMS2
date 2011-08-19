void doTable(){
gROOT->ProcessLine(".x setup.C");
TFile *file = TFile::Open("myHist.root");
 float DYee[2], DYtautau[2];
 
 for(int i=0; i<2; i++){
   DYee[i] = 0.0;
   DYtautau [i] = 0.0;
 }
 char* finalState[2];
 finalState[0] = "channel 1";
 finalState[1] = "channel 2";
 

 //for channel 1
 DYee[0]     = DYee_nJets_Ch0H0->Integral();
 DYtautau[0] = DYtautau_nJets_Ch0H0->Integral();

 //for channel 2
 DYee[1]     = DYee_nJets_Ch1H0->Integral();
 DYtautau[1] = DYtautau_nJets_Ch1H0->Integral();

 cout << "| | DYee|DYtautau" <<"|"<< endl;

for (int i=0; i<2; i++){
  
  cout << "| " << finalState[i] << " | " <<DYee[i] <<"|"<<DYtautau[i] <<"|"<< endl;
}

}
