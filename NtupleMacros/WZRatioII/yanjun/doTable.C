void doTable(){
gROOT->ProcessLine(".x setup.C");
TFile *file = TFile::Open("myHist.root");
 float DYee[4], DYtautau[4], WJets[4], ZJets[4], TTBar[4], QCD[4];
 
 for(int i=0; i<4; i++){
   DYee[i] = 0.0;
   DYtautau [i] = 0.0;
   WJets[i]=0.0;
   ZJets[i]=0.0;
   TTBar[i]=0.0;
   QCD[i]=0.0;
 }
 char* finalState[4];
 finalState[0] = "channel e";
 finalState[1] = "channel m";
 finalState[2] = "channel ee";
 finalState[3] = "channel mm";

 //for channel 1
 // DYee[0]     = DYee_nJets_e0->Integral();
//  DYtautau[0] = DYtautau_nJets_e0->Integral();
 WJets[0]     = WJets_nJets_e0->Integral();
 ZJets[0]     = ZJets_nJets_e0->Integral();
 TTBar[0]     = TTBar_nJets_e0->Integral();
 QCD[0]       = QCD_nJets_e0->Integral();

 //for channel 2
//  DYee[1]     = DYee_nJets_m0->Integral();
//  DYtautau[1] = DYtautau_nJets_m0->Integral();

 WJets[1]     = WJets_nJets_m0->Integral();
 ZJets[1]     = ZJets_nJets_m0->Integral();
 TTBar[1]     = TTBar_nJets_m0->Integral();
 QCD[1]       = QCD_nJets_m0->Integral();

 WJets[2]     = WJets_nJets_ee0->Integral();
 ZJets[2]     = ZJets_nJets_ee0->Integral();
 TTBar[2]     = TTBar_nJets_ee0->Integral();
 QCD[2]       = QCD_nJets_ee0->Integral();

 
 WJets[3]     = WJets_nJets_mm0->Integral();
 ZJets[3]     = ZJets_nJets_mm0->Integral();
 TTBar[3]     = TTBar_nJets_mm0->Integral();
 QCD[3]       = QCD_nJets_mm0->Integral();

 cout << "| | W+0jets|Z+0jets|TTBar|QCD" <<"|"<< endl;

for (int i=0; i<4; i++){
  
  cout << "| " << finalState[i] << " | " <<WJets[i] <<"|"<<ZJets[i] <<"|"<< TTBar[i] <<"|"<<QCD[i] <<"|"<<endl;
}

}
