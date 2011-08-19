{
// Load and compile something to allow proper treatment of vectors
// Not clear that it is needed
//gROOT->LoadMacro("loader.C+");

// Load various tools
//gROOT->ProcessLine(".x setup.C");

//hist file:
 TFile *ftt = TFile::Open("./processed_data_tagtcmet.root");

 float DY[4], tt[4], wjets[4], wz[4], zz[4], ww[4], tw[4], lm0x[4], lm1x[4], lm2x[4], lm3x[4], lm4x[4], lm5x[4], lm6x[4], lm7x[4], lm8x[4], lm9x[4];
 float totSM[4];

 for (int i=0; i<4; i++){
   DY[i] = 0.0;
   tt[i] = 0.0;
   wjets[i] = 0.0;
   wz[i] = 0.0;
   zz[i] = 0.0;
   ww[i] = 0.0;
   tw[i] = 0.0;
   lm0x[i] = 0.0;
   lm1x[i] = 0.0;
   lm2x[i] = 0.0;
   lm3x[i] = 0.0;
   lm4x[i] = 0.0;
   lm5x[i] = 0.0;
   lm6x[i] = 0.0;
   lm7x[i] = 0.0;
   lm8x[i] = 0.0;
   lm9x[i] = 0.0;
   totSM[i] = 0.0;
 }

//do the emu row:
 DY[2] = dy_hnJet_em->Integral();
 tt[2] = ttbar_hnJet_em->Integral();
 wjets[2] = wjets_hnJet_em->Integral();
 wz[2] = wz_hnJet_em->Integral();
 zz[2] = zz_hnJet_em->Integral();
 ww[2] = ww_hnJet_em->Integral();
 tw[2] = tw_hnJet_em->Integral();
 totSM[2] = DY[2] + tt[2] + wjets[2] + wz[2] + zz[2] + ww[2] + tw[2]; 
 lm0x[2] = lm0x_hnJet_em->Integral();
 lm1x[2] = lm1x_hnJet_em->Integral();
 lm2x[2] = lm2x_hnJet_em->Integral();
 lm3x[2] = lm3x_hnJet_em->Integral();
 lm4x[2] = lm4x_hnJet_em->Integral();
 lm5x[2] = lm5x_hnJet_em->Integral();
 lm6x[2] = lm6x_hnJet_em->Integral();
 lm7x[2] = lm7x_hnJet_em->Integral();
 lm8x[2] = lm8x_hnJet_em->Integral();
 lm9x[2] = lm9x_hnJet_em->Integral();

//do the mumu row:
 DY[1] = dy_hnJet_mm->Integral();
 tt[1] = ttbar_hnJet_mm->Integral();
 wjets[1] = wjets_hnJet_mm->Integral();
 wz[1] = wz_hnJet_mm->Integral();
 zz[1] = zz_hnJet_mm->Integral();
 ww[1] = ww_hnJet_mm->Integral();
 tw[1] = tw_hnJet_mm->Integral();
 totSM[1] = DY[1] + tt[1] + wjets[1] + wz[1] + zz[1] + ww[1] + tw[1]; 
 lm0x[1] = lm0x_hnJet_mm->Integral();
 lm1x[1] = lm1x_hnJet_mm->Integral();
 lm2x[1] = lm2x_hnJet_mm->Integral();
 lm3x[1] = lm3x_hnJet_mm->Integral();
 lm4x[1] = lm4x_hnJet_mm->Integral();
 lm5x[1] = lm5x_hnJet_mm->Integral();
 lm6x[1] = lm6x_hnJet_mm->Integral();
 lm7x[1] = lm7x_hnJet_mm->Integral();
 lm8x[1] = lm8x_hnJet_mm->Integral();
 lm9x[1] = lm9x_hnJet_mm->Integral();


//do the ee row:
 DY[0] = dy_hnJet_ee->Integral();
 tt[0] = ttbar_hnJet_ee->Integral();
 wjets[0] = wjets_hnJet_ee->Integral();
 wz[0] = wz_hnJet_ee->Integral();
 zz[0] = zz_hnJet_ee->Integral();
 ww[0] = ww_hnJet_ee->Integral();
 tw[0] = tw_hnJet_ee->Integral();
 totSM[0] = DY[0] + tt[0] + wjets[0] + wz[0] + zz[0] + ww[0] + tw[0]; 
 lm0x[0] = lm0x_hnJet_ee->Integral();
 lm1x[0] = lm1x_hnJet_ee->Integral();
 lm2x[0] = lm2x_hnJet_ee->Integral();
 lm3x[0] = lm3x_hnJet_ee->Integral();
 lm4x[0] = lm4x_hnJet_ee->Integral();
 lm5x[0] = lm5x_hnJet_ee->Integral();
 lm6x[0] = lm6x_hnJet_ee->Integral();
 lm7x[0] = lm7x_hnJet_ee->Integral();
 lm8x[0] = lm8x_hnJet_ee->Integral();
 lm9x[0] = lm9x_hnJet_ee->Integral();



 char* finalState[4];
 finalState[2] = "e&mu;";
 finalState[1] = "&mu;&mu;";
 finalState[0] = "ee";
 finalState[3] = "total";

 cout << "|Same Sign | Total SM | ttbar | TW | WZ | ZZ | WW | DY | Wjets | LM0 | LM1 | LM2 | LM3 | LM4 | LM5 | LM6 | LM7 | LM8 | LM9 |" << endl;
// cout.setf(ios::fixed, ios::floatfield);
// cout.precision(2);

 for (int i=0; i<4; i++){

//   cout << "| " << finalState[i] << " | " << totSM[i] <<  " | "  << DY[i] << " | " << tt[i] << " | " << wjets[i] << " | " << ww[i] << " | " << wz[i] << " | " << zz[i] << " | " << tw[i] << " | " << lm0x[i] << "|" << lm1x[i] << "|" << lm2x[i] << "|" << lm3x[i] << "|" << lm4x[i] << "|" << lm5x[i] << "|" << lm6x[i] << "|" << lm7x[i] << "|" << lm8x[i] << "|" << lm9x[i] << "|" << endl;
   cout << "| " << finalState[i] << " | " << totSM[i] <<  " | "  << tt[i] << " | " << tw[i] << " | " << wz[i] << " | " << zz[i] << " | " << ww[i] << " | " << DY[i] << " | " << wjets[i] << " | " << lm0x[i] << "|" << lm1x[i] << "|" << lm2x[i] << "|" << lm3x[i] << "|" << lm4x[i] << "|" << lm5x[i] << "|" << lm6x[i] << "|" << lm7x[i] << "|" << lm8x[i] << "|" << lm9x[i] << "|" << endl;

  totSM[3] = totSM[3] + totSM[i]; 
   DY[3] = DY[3]+DY[i];
   tt[3] = tt[3]+tt[i];
   wjets[3] = wjets[3]+wjets[i];
   wz[3] = wz[3]+wz[i];
   zz[3] = zz[3]+zz[i];
   ww[3] = ww[3]+ww[i];
   tw[3] = tw[3]+tw[i];
   lm0x[3] = lm0x[3]+lm0x[i];
   lm1x[3] = lm1x[3]+lm1x[i];
   lm2x[3] = lm2x[3]+lm2x[i];
   lm3x[3] = lm3x[3]+lm3x[i];
   lm4x[3] = lm4x[3]+lm4x[i];
   lm5x[3] = lm5x[3]+lm5x[i];
   lm6x[3] = lm6x[3]+lm6x[i];
   lm7x[3] = lm7x[3]+lm7x[i];
   lm8x[3] = lm8x[3]+lm8x[i];
   lm9x[3] = lm9x[3]+lm9x[i];
 }

}
