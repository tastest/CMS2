#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <vector>
#include <string>
#include "TMath.h"
#include "TDirectory.h"
#include "topAFB_looper.h"
void topAFB_looper::bookHistos(const char *prefix, int nchannel, int nhists) {
  //  Book histograms...
  //  Naming Convention:
  //  Prefix comes from the sample and it is passed to the scanning function
  //  Suffix is "ee" "em" "em" "all" which depends on the final state
  //  For example: histogram named tt_hnJet_ee would be the Njet distribution
  //  for the ee final state in the ttbar sample.
  
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  
  cout << "Begin book histos..." << endl;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
    
  char jetbins[5][7]    = {"0", "1", "2", "3", "#geq 4"};
  char suffixall[4][4]  = {"ee", "mm", "em", "all"};
  char njetCh[4][5]     = {"0j", "1j", "2j", "allj"};
  

  
  //the plots not differ by number of jets
  for (int i=0; i<4; i++) {

    for (int j = 3; j < 4; j++) {
      char suffix[7];
      sprintf(suffix, "%s_%s", njetCh[j], suffixall[i]);

      hnJet[i][j] = new TH1D(Form("%s_hnJet_%s",prefix,suffix),Form("%s_nJet_%s",prefix,suffix),9,0,9);	
      hnJet[i][j]->GetXaxis()->SetTitle("Number of jets");
      hnJet[i][j]->Sumw2();

      hnBtagJet[i][j] = new TH1D(Form("%s_hnBtagJet_%s",prefix,suffix),Form("%s_nBtagJet_%s",prefix,suffix),6,0,6);	
      hnBtagJet[i][j]->GetXaxis()->SetTitle("Number of b tagged jets");
      hnBtagJet[i][j]->Sumw2();

      hnVtx[i][j] = new TH1D(Form("%s_hnVtx_%s",prefix,suffix),Form("%s_nVtx_%s",prefix,suffix),21,0.,21.);	
      hnVtx[i][j]->GetXaxis()->SetTitle("Number of vertices");
      hnVtx[i][j]->Sumw2();

      hNsolns[i][j] = new TH1D(Form("%s_hNsolns_%s",prefix,suffix),Form("%s_Nsolns_%s",prefix,suffix),120,-19,101);
      hNsolns[i][j]->GetXaxis()->SetTitle("Jet smearing solution multiplicity");
      hNsolns[i][j]->Sumw2();

      hmaxAMWTweight[i][j] = new TH1D(Form("%s_hmaxAMWTweight_%s",prefix,suffix),Form("%s_maxAMWTweight_%s",prefix,suffix),120,-0.2,1.0);
      hmaxAMWTweight[i][j]->GetXaxis()->SetTitle("Maximum AMWT weight");
      hmaxAMWTweight[i][j]->Sumw2();

      haveAMWTweight[i][j] = new TH1D(Form("%s_haveAMWTweight_%s",prefix,suffix),Form("%s_aveAMWTweight_%s",prefix,suffix),120,-0.2,1.0);
      haveAMWTweight[i][j]->GetXaxis()->SetTitle("Average AMWT weight");
      haveAMWTweight[i][j]->Sumw2();

      hsumAMWTweight[i][j] = new TH1D(Form("%s_hsumAMWTweight_%s",prefix,suffix),Form("%s_sumAMWTweight_%s",prefix,suffix),120,-20.,100.);
      hsumAMWTweight[i][j]->GetXaxis()->SetTitle("Sum of AMWT weights");
      hsumAMWTweight[i][j]->Sumw2();

      hAMWTweightnojetsmear[i][j] = new TH1D(Form("%s_hAMWTweightnojetsmear_%s",prefix,suffix),Form("%s_AMWTweightnojetsmear_%s",prefix,suffix),120,-0.2,1.0);
      hAMWTweightnojetsmear[i][j]->GetXaxis()->SetTitle("AMWT weight (no jet smearing)");
      hAMWTweightnojetsmear[i][j]->Sumw2();


      
      hlepChargeAsym_2d[i][j] = new TH2D(Form("%s_hlepChargeAsym2d_%s",prefix,suffix),Form("%s_lepChargeAsym2d_%s",prefix,suffix),80,-4,4, 80,-4,4);
      hlepChargeAsym_2d[i][j]->GetXaxis()->SetTitle("Charge_Asymmetry_lep_gen");
      hlepChargeAsym_2d[i][j]->GetYaxis()->SetTitle("Charge_Asymmetry_lep_rec");
      hlepChargeAsym_2d[i][j]->Sumw2();

      hlepAzimAsym_2d[i][j] = new TH2D(Form("%s_hlepAzimAsym2d_%s",prefix,suffix),Form("%s_lepAzimAsym2d_%s",prefix,suffix),80,-1,1, 80,-1,1);
      hlepAzimAsym_2d[i][j]->GetXaxis()->SetTitle("Azimuthal_Asymmetry_lep_gen");
      hlepAzimAsym_2d[i][j]->GetYaxis()->SetTitle("Azimuthal_Asymmetry_lep_rec");
      hlepAzimAsym_2d[i][j]->Sumw2();

      htopSpinCorr_2d[i][j] = new TH2D(Form("%s_htopSpinCorr2d_%s",prefix,suffix),Form("%s_topSpinCorr2d_%s",prefix,suffix),80,-1,1,80,-1,1);
      htopSpinCorr_2d[i][j]->GetXaxis()->SetTitle("Spin_Correlation_top_gen");
      htopSpinCorr_2d[i][j]->GetYaxis()->SetTitle("Spin_Correlation_top_rec");
      htopSpinCorr_2d[i][j]->Sumw2();
      
      htopCosTheta_2d[i][j] = new TH2D(Form("%s_htopCosTheta2d_%s",prefix,suffix),Form("%s_topCosTheta2d_%s",prefix,suffix),80,-1,1, 80,-1,1);
      htopCosTheta_2d[i][j]->GetXaxis()->SetTitle("Cos(theta)_top_gen");
      htopCosTheta_2d[i][j]->GetYaxis()->SetTitle("Cos(theta)_top_rec");
      htopCosTheta_2d[i][j]->Sumw2();

      hlepCosTheta_2d[i][j] = new TH2D(Form("%s_hlepCosTheta2d_%s",prefix,suffix),Form("%s_lepCosTheta2d_%s",prefix,suffix),80,-1,1, 80,-1,1);
      hlepCosTheta_2d[i][j]->GetXaxis()->SetTitle("Cos(theta)_lep_gen");
      hlepCosTheta_2d[i][j]->GetYaxis()->SetTitle("Cos(theta)_lep_rec");
      hlepCosTheta_2d[i][j]->Sumw2();
      
      
      
      hlepChargeAsym_gen[i][j] = new TH1D(Form("%s_hlepChargeAsymGen_%s",prefix,suffix),Form("%s_lepChargeAsymGen_%s",prefix,suffix),80,-4,4);
      hlepChargeAsym_gen[i][j]->GetXaxis()->SetTitle("|#eta_{l^{+}}| - |#eta_{l^{-}}|");
      hlepChargeAsym_gen[i][j]->Sumw2();

      hlepAzimAsym_gen[i][j] = new TH1D(Form("%s_hlepAzimAsymGen_%s",prefix,suffix),Form("%s_lepAzimAsymGen_%s",prefix,suffix),80,-1,1);
      hlepAzimAsym_gen[i][j]->GetXaxis()->SetTitle("cos(#Delta #phi_{l^{+}l^{-}})");
      hlepAzimAsym_gen[i][j]->Sumw2();

      hlepAzimAsym2_gen[i][j] = new TH1D(Form("%s_hlepAzimAsym2Gen_%s",prefix,suffix),Form("%s_lepAzimAsym2Gen_%s",prefix,suffix),80,0,3.141592653589793);
      hlepAzimAsym2_gen[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym2_gen[i][j]->Sumw2();

      htopSpinCorr_gen[i][j] = new TH1D(Form("%s_htopSpinCorrGen_%s",prefix,suffix),Form("%s_topSpinCorrGen_%s",prefix,suffix),80,-1,1);
      htopSpinCorr_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) #times cos(#theta_{l^{-}}^{#bar{t}} )");
      htopSpinCorr_gen[i][j]->Sumw2();
      
      htopCosTheta_gen[i][j] = new TH1D(Form("%s_htopCosThetaGen_%s",prefix,suffix),Form("%s_topCosThetaGen_%s",prefix,suffix),80,-1,1);
      htopCosTheta_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{t}^{t#bar{t}})");
      htopCosTheta_gen[i][j]->Sumw2();

      hlepCosTheta_gen[i][j] = new TH1D(Form("%s_hlepCosThetaGen_%s",prefix,suffix),Form("%s_lepCosThetaGen_%s",prefix,suffix),80,-1,1);
      hlepCosTheta_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{l}^{t})");
      hlepCosTheta_gen[i][j]->Sumw2();

      hlepPlusCosTheta_gen[i][j] = new TH1D(Form("%s_hlepPlusCosThetaGen_%s",prefix,suffix),Form("%s_lepPlusCosThetaGen_%s",prefix,suffix),80,-1,1);
      hlepPlusCosTheta_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t})");
      hlepPlusCosTheta_gen[i][j]->Sumw2();
      
      hlepMinusCosTheta_gen[i][j] = new TH1D(Form("%s_hlepMinusCosThetaGen_%s",prefix,suffix),Form("%s_lepMinusCosThetaGen_%s",prefix,suffix),80,-1,1);
      hlepMinusCosTheta_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#bar{t}})");
      hlepMinusCosTheta_gen[i][j]->Sumw2();
            
      hpseudorapiditydiff_gen[i][j] = new TH1D(Form("%s_hpseudorapiditydiffGen_%s",prefix,suffix),Form("%s_pseudorapiditydiffGen_%s",prefix,suffix),80,-4,4);
      hpseudorapiditydiff_gen[i][j]->GetXaxis()->SetTitle("|#eta_{t}| - |#eta_{#bar{t}}|");
      hpseudorapiditydiff_gen[i][j]->Sumw2();
      
      hrapiditydiff_gen[i][j] = new TH1D(Form("%s_hrapiditydiffGen_%s",prefix,suffix),Form("%s_rapiditydiffGen_%s",prefix,suffix),80,-4,4);
      hrapiditydiff_gen[i][j]->GetXaxis()->SetTitle("(y_{t}-y_{#bar{t}}) #times (y_{t}+y_{#bar{t}})");
      hrapiditydiff_gen[i][j]->Sumw2();
      
      hrapiditydiffMarco_gen[i][j] = new TH1D(Form("%s_hrapiditydiffMarcoGen_%s",prefix,suffix),Form("%s_rapiditydiffMarcoGen_%s",prefix,suffix),80,-4,4);
      hrapiditydiffMarco_gen[i][j]->GetXaxis()->SetTitle("|y_{t}|-|y_{#bar{t}}|");
      hrapiditydiffMarco_gen[i][j]->Sumw2();      
      
      
      
      
      //2d histos for asyms vs mass at gen level
      hlepChargeAsym_gen2d[i][j] = new TH2D(Form("%s_hlepChargeAsymGen2d_%s",prefix,suffix),Form("%s_lepChargeAsymGen2d_%s",prefix,suffix),80,-4,4,120,0.,1200.);
      hlepChargeAsym_gen2d[i][j]->GetXaxis()->SetTitle("|#eta_{l^{+}}| - |#eta_{l^{-}}|");
      hlepChargeAsym_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepChargeAsym_gen2d[i][j]->Sumw2();

      hlepAzimAsym_gen2d[i][j] = new TH2D(Form("%s_hlepAzimAsymGen2d_%s",prefix,suffix),Form("%s_lepAzimAsymGen2d_%s",prefix,suffix),80,-1,1,120,0.,1200.);
      hlepAzimAsym_gen2d[i][j]->GetXaxis()->SetTitle("cos(#Delta #phi_{l^{+}l^{-}})");
      hlepAzimAsym_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepAzimAsym_gen2d[i][j]->Sumw2();

      hlepAzimAsym2_gen2d[i][j] = new TH2D(Form("%s_hlepAzimAsym2Gen2d_%s",prefix,suffix),Form("%s_lepAzimAsym2Gen2d_%s",prefix,suffix),80,0,3.141592653589793,120,0.,1200.);
      hlepAzimAsym2_gen2d[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym2_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepAzimAsym2_gen2d[i][j]->Sumw2();

      htopSpinCorr_gen2d[i][j] = new TH2D(Form("%s_htopSpinCorrGen2d_%s",prefix,suffix),Form("%s_topSpinCorrGen2d_%s",prefix,suffix),80,-1,1,120,0.,1200.);
      htopSpinCorr_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) #times cos(#theta_{l^{-}}^{#bar{t}} )");
      htopSpinCorr_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      htopSpinCorr_gen2d[i][j]->Sumw2();
      
      htopCosTheta_gen2d[i][j] = new TH2D(Form("%s_htopCosThetaGen2d_%s",prefix,suffix),Form("%s_topCosThetaGen2d_%s",prefix,suffix),80,-1,1,120,0.,1200.);
      htopCosTheta_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{t}^{t#bar{t}})");
      htopCosTheta_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      htopCosTheta_gen2d[i][j]->Sumw2();

      hlepCosTheta_gen2d[i][j] = new TH2D(Form("%s_hlepCosThetaGen2d_%s",prefix,suffix),Form("%s_lepCosThetaGen2d_%s",prefix,suffix),80,-1,1,120,0.,1200.);
      hlepCosTheta_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l}^{t})");
      hlepCosTheta_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepCosTheta_gen2d[i][j]->Sumw2();

      hlepPlusCosTheta_gen2d[i][j] = new TH2D(Form("%s_hlepPlusCosThetaGen2d_%s",prefix,suffix),Form("%s_lepPlusCosThetaGen2d_%s",prefix,suffix),80,-1,1,120,0.,1200.);
      hlepPlusCosTheta_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t})");
      hlepPlusCosTheta_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepPlusCosTheta_gen2d[i][j]->Sumw2();
      
      hlepMinusCosTheta_gen2d[i][j] = new TH2D(Form("%s_hlepMinusCosThetaGen2d_%s",prefix,suffix),Form("%s_lepMinusCosThetaGen2d_%s",prefix,suffix),80,-1,1,120,0.,1200.);
      hlepMinusCosTheta_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#bar{t}})");
      hlepMinusCosTheta_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepMinusCosTheta_gen2d[i][j]->Sumw2();
            
      hpseudorapiditydiff_gen2d[i][j] = new TH2D(Form("%s_hpseudorapiditydiffGen2d_%s",prefix,suffix),Form("%s_pseudorapiditydiffGen2d_%s",prefix,suffix),80,-4,4,120,0.,1200.);
      hpseudorapiditydiff_gen2d[i][j]->GetXaxis()->SetTitle("|#eta_{t}| - |#eta_{#bar{t}}|");
      hpseudorapiditydiff_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hpseudorapiditydiff_gen2d[i][j]->Sumw2();
      
      hrapiditydiff_gen2d[i][j] = new TH2D(Form("%s_hrapiditydiffGen2d_%s",prefix,suffix),Form("%s_rapiditydiffGen2d_%s",prefix,suffix),80,-4,4,120,0.,1200.);
      hrapiditydiff_gen2d[i][j]->GetXaxis()->SetTitle("(y_{t}-y_{#bar{t}}) #times (y_{t}+y_{#bar{t}})");
      hrapiditydiff_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hrapiditydiff_gen2d[i][j]->Sumw2();
      
      hrapiditydiffMarco_gen2d[i][j] = new TH2D(Form("%s_hrapiditydiffMarcoGen2d_%s",prefix,suffix),Form("%s_rapiditydiffMarcoGen2d_%s",prefix,suffix),80,-4,4,120,0.,1200.);
      hrapiditydiffMarco_gen2d[i][j]->GetXaxis()->SetTitle("|y_{t}|-|y_{#bar{t}}|");
      hrapiditydiffMarco_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hrapiditydiffMarco_gen2d[i][j]->Sumw2();



      hlepChargeAsym_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepChargeAsymttpTGen2d_%s",prefix,suffix),Form("%s_lepChargeAsymttpTGen2d_%s",prefix,suffix),80,-4,4,300.,0.,300.);
      hlepChargeAsym_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("|#eta_{l^{+}}| - |#eta_{l^{-}}|");
      hlepChargeAsym_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepChargeAsym_ttpT_gen2d[i][j]->Sumw2();

      hlepAzimAsym_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepAzimAsymttpTGen2d_%s",prefix,suffix),Form("%s_lepAzimAsymttpTGen2d_%s",prefix,suffix),80,-1,1,300.,0.,300.);
      hlepAzimAsym_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("cos(#Delta #phi_{l^{+}l^{-}})");
      hlepAzimAsym_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepAzimAsym_ttpT_gen2d[i][j]->Sumw2();

      hlepAzimAsym2_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepAzimAsym2ttpTGen2d_%s",prefix,suffix),Form("%s_lepAzimAsym2ttpTGen2d_%s",prefix,suffix),80,0,3.141592653589793,300.,0.,300.);
      hlepAzimAsym2_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym2_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepAzimAsym2_ttpT_gen2d[i][j]->Sumw2();

      htopSpinCorr_ttpT_gen2d[i][j] = new TH2D(Form("%s_htopSpinCorrttpTGen2d_%s",prefix,suffix),Form("%s_topSpinCorrttpTGen2d_%s",prefix,suffix),80,-1,1,300.,0.,300.);
      htopSpinCorr_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) #times cos(#theta_{l^{-}}^{#bar{t}} )");
      htopSpinCorr_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      htopSpinCorr_ttpT_gen2d[i][j]->Sumw2();
      
      htopCosTheta_ttpT_gen2d[i][j] = new TH2D(Form("%s_htopCosThetattpTGen2d_%s",prefix,suffix),Form("%s_topCosThetattpTGen2d_%s",prefix,suffix),80,-1,1,300.,0.,300.);
      htopCosTheta_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{t}^{t#bar{t}})");
      htopCosTheta_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      htopCosTheta_ttpT_gen2d[i][j]->Sumw2();

      hlepCosTheta_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepCosThetattpTGen2d_%s",prefix,suffix),Form("%s_lepCosThetattpTGen2d_%s",prefix,suffix),80,-1,1,300.,0.,300.);
      hlepCosTheta_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l}^{t})");
      hlepCosTheta_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepCosTheta_ttpT_gen2d[i][j]->Sumw2();

      hlepPlusCosTheta_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepPlusCosThetattpTGen2d_%s",prefix,suffix),Form("%s_lepPlusCosThetattpTGen2d_%s",prefix,suffix),80,-1,1,300.,0.,300.);
      hlepPlusCosTheta_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t})");
      hlepPlusCosTheta_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepPlusCosTheta_ttpT_gen2d[i][j]->Sumw2();
      
      hlepMinusCosTheta_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepMinusCosThetattpTGen2d_%s",prefix,suffix),Form("%s_lepMinusCosThetattpTGen2d_%s",prefix,suffix),80,-1,1,300.,0.,300.);
      hlepMinusCosTheta_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#bar{t}})");
      hlepMinusCosTheta_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepMinusCosTheta_ttpT_gen2d[i][j]->Sumw2();
            
      hpseudorapiditydiff_ttpT_gen2d[i][j] = new TH2D(Form("%s_hpseudorapiditydiffttpTGen2d_%s",prefix,suffix),Form("%s_pseudorapiditydiffttpTGen2d_%s",prefix,suffix),80,-4,4,300.,0.,300.);
      hpseudorapiditydiff_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("|#eta_{t}| - |#eta_{#bar{t}}|");
      hpseudorapiditydiff_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hpseudorapiditydiff_ttpT_gen2d[i][j]->Sumw2();
      
      hrapiditydiff_ttpT_gen2d[i][j] = new TH2D(Form("%s_hrapiditydiffttpTGen2d_%s",prefix,suffix),Form("%s_rapiditydiffttpTGen2d_%s",prefix,suffix),80,-4,4,300.,0.,300.);
      hrapiditydiff_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("(y_{t}-y_{#bar{t}}) #times (y_{t}+y_{#bar{t}})");
      hrapiditydiff_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hrapiditydiff_ttpT_gen2d[i][j]->Sumw2();
      
      hrapiditydiffMarco_ttpT_gen2d[i][j] = new TH2D(Form("%s_hrapiditydiffMarcottpTGen2d_%s",prefix,suffix),Form("%s_rapiditydiffMarcottpTGen2d_%s",prefix,suffix),80,-4,4,300.,0.,300.);
      hrapiditydiffMarco_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("|y_{t}|-|y_{#bar{t}}|");
      hrapiditydiffMarco_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hrapiditydiffMarco_ttpT_gen2d[i][j]->Sumw2();


    
      hlepChargeAsym_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepChargeAsymttRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepChargeAsymttRapidity2Gen2d_%s",prefix,suffix),80,-4,4,90,0.,3.0);
      hlepChargeAsym_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("|#eta_{l^{+}}| - |#eta_{l^{-}}|");
      hlepChargeAsym_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepChargeAsym_ttRapidity2_gen2d[i][j]->Sumw2();

      hlepAzimAsym_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepAzimAsymttRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepAzimAsymttRapidity2Gen2d_%s",prefix,suffix),80,-1,1,90,0.,3.0);
      hlepAzimAsym_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("cos(#Delta #phi_{l^{+}l^{-}})");
      hlepAzimAsym_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepAzimAsym_ttRapidity2_gen2d[i][j]->Sumw2();

      hlepAzimAsym2_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepAzimAsym2ttRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepAzimAsym2ttRapidity2Gen2d_%s",prefix,suffix),80,0,3.141592653589793,90,0.,3.0);
      hlepAzimAsym2_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym2_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepAzimAsym2_ttRapidity2_gen2d[i][j]->Sumw2();

      htopSpinCorr_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_htopSpinCorrttRapidity2Gen2d_%s",prefix,suffix),Form("%s_topSpinCorrttRapidity2Gen2d_%s",prefix,suffix),80,-1,1,90,0.,3.0);
      htopSpinCorr_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) #times cos(#theta_{l^{-}}^{#bar{t}} )");
      htopSpinCorr_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      htopSpinCorr_ttRapidity2_gen2d[i][j]->Sumw2();
      
      htopCosTheta_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_htopCosThetattRapidity2Gen2d_%s",prefix,suffix),Form("%s_topCosThetattRapidity2Gen2d_%s",prefix,suffix),80,-1,1,90,0.,3.0);
      htopCosTheta_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{t}^{t#bar{t}})");
      htopCosTheta_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      htopCosTheta_ttRapidity2_gen2d[i][j]->Sumw2();

      hlepCosTheta_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepCosThetattRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepCosThetattRapidity2Gen2d_%s",prefix,suffix),80,-1,1,90,0.,3.0);
      hlepCosTheta_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l}^{t})");
      hlepCosTheta_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepCosTheta_ttRapidity2_gen2d[i][j]->Sumw2();

      hlepPlusCosTheta_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepPlusCosThetattRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepPlusCosThetattRapidity2Gen2d_%s",prefix,suffix),80,-1,1,90,0.,3.0);
      hlepPlusCosTheta_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t})");
      hlepPlusCosTheta_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepPlusCosTheta_ttRapidity2_gen2d[i][j]->Sumw2();
      
      hlepMinusCosTheta_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepMinusCosThetattRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepMinusCosThetattRapidity2Gen2d_%s",prefix,suffix),80,-1,1,90,0.,3.0);
      hlepMinusCosTheta_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#bar{t}})");
      hlepMinusCosTheta_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepMinusCosTheta_ttRapidity2_gen2d[i][j]->Sumw2();
            
      hpseudorapiditydiff_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hpseudorapiditydiffttRapidity2Gen2d_%s",prefix,suffix),Form("%s_pseudorapiditydiffttRapidity2Gen2d_%s",prefix,suffix),80,-4,4,90,0.,3.0);
      hpseudorapiditydiff_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("|#eta_{t}| - |#eta_{#bar{t}}|");
      hpseudorapiditydiff_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hpseudorapiditydiff_ttRapidity2_gen2d[i][j]->Sumw2();
      
      hrapiditydiff_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hrapiditydiffttRapidity2Gen2d_%s",prefix,suffix),Form("%s_rapiditydiffttRapidity2Gen2d_%s",prefix,suffix),80,-4,4,90,0.,3.0);
      hrapiditydiff_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("(y_{t}-y_{#bar{t}}) #times (y_{t}+y_{#bar{t}})");
      hrapiditydiff_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hrapiditydiff_ttRapidity2_gen2d[i][j]->Sumw2();
      
      hrapiditydiffMarco_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hrapiditydiffMarcottRapidity2Gen2d_%s",prefix,suffix),Form("%s_rapiditydiffMarcottRapidity2Gen2d_%s",prefix,suffix),80,-4,4,90,0.,3.0);
      hrapiditydiffMarco_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("|y_{t}|-|y_{#bar{t}}|");
      hrapiditydiffMarco_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hrapiditydiffMarco_ttRapidity2_gen2d[i][j]->Sumw2();



            
      
      //Reco-Gen hists for asyms
      hlepChargeAsymGenDiff[i][j] = new TH1D(Form("%s_hlepChargeAsymGenDiff_%s",prefix,suffix),Form("%s_lepChargeAsymGenDiff_%s",prefix,suffix),80,-4,4);
      hlepChargeAsymGenDiff[i][j]->GetXaxis()->SetTitle("|#eta_{l^{+}}| - |#eta_{l^{-}}| (reco-gen)");
      hlepChargeAsymGenDiff[i][j]->Sumw2();

      hlepAzimAsymGenDiff[i][j] = new TH1D(Form("%s_hlepAzimAsymGenDiff_%s",prefix,suffix),Form("%s_lepAzimAsymGenDiff_%s",prefix,suffix),80,-2,2);
      hlepAzimAsymGenDiff[i][j]->GetXaxis()->SetTitle("cos(#Delta #phi_{l^{+}l^{-}}) (reco-gen)");
      hlepAzimAsymGenDiff[i][j]->Sumw2();

      hlepAzimAsym2GenDiff[i][j] = new TH1D(Form("%s_hlepAzimAsym2GenDiff_%s",prefix,suffix),Form("%s_lepAzimAsym2GenDiff_%s",prefix,suffix),80,-3.141592653589793,3.141592653589793);
      hlepAzimAsym2GenDiff[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}} (reco-gen)");
      hlepAzimAsym2GenDiff[i][j]->Sumw2();

      htopSpinCorrGenDiff[i][j] = new TH1D(Form("%s_htopSpinCorrGenDiff_%s",prefix,suffix),Form("%s_topSpinCorrGenDiff_%s",prefix,suffix),80,-2,2);
      htopSpinCorrGenDiff[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) #times cos(#theta_{l^{-}}^{#bar{t}} ) (reco-gen)");
      htopSpinCorrGenDiff[i][j]->Sumw2();
      
      htopCosThetaGenDiff[i][j] = new TH1D(Form("%s_htopCosThetaGenDiff_%s",prefix,suffix),Form("%s_topCosThetaGenDiff_%s",prefix,suffix),80,-2,2);
      htopCosThetaGenDiff[i][j]->GetXaxis()->SetTitle("cos(#theta_{t}^{t#bar{t}}) (reco-gen)");
      htopCosThetaGenDiff[i][j]->Sumw2();

      hlepCosThetaGenDiff[i][j] = new TH1D(Form("%s_hlepCosThetaGenDiff_%s",prefix,suffix),Form("%s_lepCosThetaGenDiff_%s",prefix,suffix),80,-2,2);
      hlepCosThetaGenDiff[i][j]->GetXaxis()->SetTitle("cos(#theta_{l}^{t}) (reco-gen)");
      hlepCosThetaGenDiff[i][j]->Sumw2();

      hlepPlusCosThetaGenDiff[i][j] = new TH1D(Form("%s_hlepPlusCosThetaGenDiff_%s",prefix,suffix),Form("%s_lepPlusCosThetaGenDiff_%s",prefix,suffix),80,-2,2);
      hlepPlusCosThetaGenDiff[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) (reco-gen)");
      hlepPlusCosThetaGenDiff[i][j]->Sumw2();
      
      hlepMinusCosThetaGenDiff[i][j] = new TH1D(Form("%s_hlepMinusCosThetaGenDiff_%s",prefix,suffix),Form("%s_lepMinusCosThetaGenDiff_%s",prefix,suffix),80,-2,2);
      hlepMinusCosThetaGenDiff[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#bar{t}}) (reco-gen)");
      hlepMinusCosThetaGenDiff[i][j]->Sumw2();
            
      hpseudorapiditydiffGenDiff[i][j] = new TH1D(Form("%s_hpseudorapiditydiffGenDiff_%s",prefix,suffix),Form("%s_pseudorapiditydiffGenDiff_%s",prefix,suffix),80,-4,4);
      hpseudorapiditydiffGenDiff[i][j]->GetXaxis()->SetTitle("|#eta_{t}| - |#eta_{#bar{t}}| (reco-gen)");
      hpseudorapiditydiffGenDiff[i][j]->Sumw2();
      
      hrapiditydiffGenDiff[i][j] = new TH1D(Form("%s_hrapiditydiffGenDiff_%s",prefix,suffix),Form("%s_rapiditydiffGenDiff_%s",prefix,suffix),80,-4,4);
      hrapiditydiffGenDiff[i][j]->GetXaxis()->SetTitle("(y_{t}-y_{#bar{t}}) #times (y_{t}+y_{#bar{t}}) (reco-gen)");
      hrapiditydiffGenDiff[i][j]->Sumw2();
      
      hrapiditydiffMarcoGenDiff[i][j] = new TH1D(Form("%s_hrapiditydiffMarcoGenDiff_%s",prefix,suffix),Form("%s_rapiditydiffMarcoGenDiff_%s",prefix,suffix),80,-4,4);
      hrapiditydiffMarcoGenDiff[i][j]->GetXaxis()->SetTitle("|y_{t}|-|y_{#bar{t}}| (reco-gen)");
      hrapiditydiffMarcoGenDiff[i][j]->Sumw2();
            
      


      
 
      hlepChargeAsym[i][j] = new TH1D(Form("%s_hlepChargeAsym_%s",prefix,suffix),Form("%s_lepChargeAsym_%s",prefix,suffix),80,-4,4);
      hlepChargeAsym[i][j]->GetXaxis()->SetTitle(" |#eta_{l^{+}}| - |#eta_{l^{-}}| ");
      hlepChargeAsym[i][j]->Sumw2();
 
      hlepRapDiff[i][j] = new TH1D(Form("%s_hlepRapDiff_%s",prefix,suffix),Form("%s_lepRapDiff_%s",prefix,suffix),100,-4,4);
      hlepRapDiff[i][j]->GetXaxis()->SetTitle(" #eta_{l^{+}} - #eta_{l^{-}} ");
      hlepRapDiff[i][j]->Sumw2();

      hlepAngleBetween[i][j] = new TH1D(Form("%s_hlepAngleBetween_%s",prefix,suffix),Form("%s_lepAngleBetween_%s",prefix,suffix),80,-1,1);
      hlepAngleBetween[i][j]->GetXaxis()->SetTitle("cos( #alpha_{l^{+}l^{-}}^{lab}) ");
      hlepAngleBetween[i][j]->Sumw2();
      
      hlepAngleBetweenCMS[i][j] = new TH1D(Form("%s_hlepAngleBetweenCMS_%s",prefix,suffix),Form("%s_lepAngleBetweenCMS_%s",prefix,suffix),80,-1,1);
      hlepAngleBetweenCMS[i][j]->GetXaxis()->SetTitle("cos( #alpha_{l^{+}l^{-}}^{t#bar{t}}) ");
      hlepAngleBetweenCMS[i][j]->Sumw2();


      hpseudorapiditydiff2[i][j] = new TH1D(Form("%s_hpseudorapiditydiff2_%s",prefix,suffix),Form("%s_pseudorapiditydiff2_%s",prefix,suffix),100,-6,6);
      hpseudorapiditydiff2[i][j]->GetXaxis()->SetTitle("#eta_{t} - #eta_{#bar{t}}");
      hpseudorapiditydiff2[i][j]->Sumw2();
      
      hrapiditydiff2[i][j] = new TH1D(Form("%s_hrapiditydiff2_%s",prefix,suffix),Form("%s_rapiditydiff2_%s",prefix,suffix),100,-4,4);
      hrapiditydiff2[i][j]->GetXaxis()->SetTitle("y_{t}-y_{#bar{t}} ");
      hrapiditydiff2[i][j]->Sumw2();

      hlepPlusCosTheta[i][j] = new TH1D(Form("%s_hlepPlusCosTheta_%s",prefix,suffix),Form("%s_lepPlusCosTheta_%s",prefix,suffix),80,-1,1);
      hlepPlusCosTheta[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t})");
      hlepPlusCosTheta[i][j]->Sumw2();
      
      hlepMinusCosTheta[i][j] = new TH1D(Form("%s_hlepMinusCosTheta_%s",prefix,suffix),Form("%s_lepMinusCosTheta_%s",prefix,suffix),80,-1,1);
      hlepMinusCosTheta[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#bar{t}})");
      hlepMinusCosTheta[i][j]->Sumw2();


      hjetAzimAsym[i][j] = new TH1D(Form("%s_hjetAzimAsym_%s",prefix,suffix),Form("%s_jetAzimAsym_%s",prefix,suffix),80,-1,1);
      hjetAzimAsym[i][j]->GetXaxis()->SetTitle("cos(#Delta #phi_{j1j2})");
      hjetAzimAsym[i][j]->Sumw2();
      
      hjetRapDiff[i][j] = new TH1D(Form("%s_hjetRapDiff_%s",prefix,suffix),Form("%s_jetRapDiff_%s",prefix,suffix),100,-4,4);
      hjetRapDiff[i][j]->GetXaxis()->SetTitle(" #eta_{j1} - #eta_{j2} ");
      hjetRapDiff[i][j]->Sumw2();
      
      hjetAngleBetween[i][j] = new TH1D(Form("%s_hjetAngleBetween_%s",prefix,suffix),Form("%s_jetAngleBetween_%s",prefix,suffix),80,-1,1);
      hjetAngleBetween[i][j]->GetXaxis()->SetTitle("cos( #alpha_{j1j2}^{lab}) ");
      hjetAngleBetween[i][j]->Sumw2();
      
      hjetAngleBetweenCMS[i][j] = new TH1D(Form("%s_hjetAngleBetweenCMS_%s",prefix,suffix),Form("%s_jetAngleBetweenCMS_%s",prefix,suffix),80,-1,1);
      hjetAngleBetweenCMS[i][j]->GetXaxis()->SetTitle("cos( #alpha_{j1j2}^{t#bar{t}}) ");
      hjetAngleBetweenCMS[i][j]->Sumw2();

      hlepPhi[i][j] = new TH1D(Form("%s_hlepPhi_%s",prefix,suffix),Form("%s_lepPhi_%s",prefix,suffix),80,-3.141592653589793,3.141592653589793);
      hlepPhi[i][j]->GetXaxis()->SetTitle("Lepton #phi");
      hlepPhi[i][j]->Sumw2();
            
      hlepPlusPhi[i][j] = new TH1D(Form("%s_hlepPlusPhi_%s",prefix,suffix),Form("%s_lepPlusPhi_%s",prefix,suffix),80,-3.141592653589793,3.141592653589793);
      hlepPlusPhi[i][j]->GetXaxis()->SetTitle("#phi_{l^{+}}");
      hlepPlusPhi[i][j]->Sumw2();
            
      hlepMinusPhi[i][j] = new TH1D(Form("%s_hlepMinusPhi_%s",prefix,suffix),Form("%s_lepMinusPhi_%s",prefix,suffix),80,-3.141592653589793,3.141592653589793);
      hlepMinusPhi[i][j]->GetXaxis()->SetTitle("#phi_{l^{-}}");
      hlepMinusPhi[i][j]->Sumw2();

      hjetPhi[i][j] = new TH1D(Form("%s_hjetPhi_%s",prefix,suffix),Form("%s_jetPhi_%s",prefix,suffix),80,-3.141592653589793,3.141592653589793);
      hjetPhi[i][j]->GetXaxis()->SetTitle("Jet #phi");
      hjetPhi[i][j]->Sumw2();

      hlepPlusEta[i][j] = new TH1D(Form("%s_hlepPlusEta_%s",prefix,suffix),Form("%s_lepPlusEta_%s",prefix,suffix),60,-3.,3.);
      hlepPlusEta[i][j]->GetXaxis()->SetTitle("#eta_{l^{+}}");
      hlepPlusEta[i][j]->Sumw2();

      hlepMinusEta[i][j] = new TH1D(Form("%s_hlepMinusEta_%s",prefix,suffix),Form("%s_lepMinusEta_%s",prefix,suffix),60,-3.,3.);
      hlepMinusEta[i][j]->GetXaxis()->SetTitle("#eta_{l^{-}}");
      hlepMinusEta[i][j]->Sumw2();

      hlepPlusPt[i][j] = new TH1D(Form("%s_hlepPlusPt_%s",prefix,suffix),Form("%s_lepPlusPt_%s",prefix,suffix),60,0.,240.);
      hlepPlusPt[i][j]->GetXaxis()->SetTitle("l^{+} p_{T} (GeV/c)");
      hlepPlusPt[i][j]->Sumw2();
      
      hlepMinusPt[i][j] = new TH1D(Form("%s_hlepMinusPt_%s",prefix,suffix),Form("%s_lepMinusPt_%s",prefix,suffix),60,0.,240.);
      hlepMinusPt[i][j]->GetXaxis()->SetTitle("l^{-} p_{T} (GeV/c)");
      hlepMinusPt[i][j]->Sumw2();








      hlepAzimAsym[i][j] = new TH1D(Form("%s_hlepAzimAsym_%s",prefix,suffix),Form("%s_lepAzimAsym_%s",prefix,suffix),80,-1,1);
      hlepAzimAsym[i][j]->GetXaxis()->SetTitle("cos(#Delta #phi_{l^{+}l^{-}})");
      hlepAzimAsym[i][j]->Sumw2();

      hlepAzimAsym2[i][j] = new TH1D(Form("%s_hlepAzimAsym2_%s",prefix,suffix),Form("%s_lepAzimAsym2_%s",prefix,suffix),80,0,3.141592653589793);
      hlepAzimAsym2[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym2[i][j]->Sumw2();

      htopSpinCorr[i][j] = new TH1D(Form("%s_htopSpinCorr_%s",prefix,suffix),Form("%s_topSpinCorr_%s",prefix,suffix),80,-1,1);
      htopSpinCorr[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) #times cos(#theta_{l^{-}}^{#bar{t}} )");
      htopSpinCorr[i][j]->Sumw2();
      
      htopCosTheta[i][j] = new TH1D(Form("%s_htopCosTheta_%s",prefix,suffix),Form("%s_topCosTheta_%s",prefix,suffix),80,-1,1);
      htopCosTheta[i][j]->GetXaxis()->SetTitle("cos(#theta_{t}^{t#bar{t}})");
      htopCosTheta[i][j]->Sumw2();

      hpseudorapiditydiff[i][j] = new TH1D(Form("%s_hpseudorapiditydiff_%s",prefix,suffix),Form("%s_pseudorapiditydiff_%s",prefix,suffix),80,-4,4);
      hpseudorapiditydiff[i][j]->GetXaxis()->SetTitle("|#eta_{t}| - |#eta_{#bar{t}}|");
      hpseudorapiditydiff[i][j]->Sumw2();
      
      hrapiditydiff[i][j] = new TH1D(Form("%s_hrapiditydiff_%s",prefix,suffix),Form("%s_rapiditydiff_%s",prefix,suffix),80,-4,4);
      hrapiditydiff[i][j]->GetXaxis()->SetTitle("(y_{t}-y_{#bar{t}}) #times (y_{t}+y_{#bar{t}})");
      hrapiditydiff[i][j]->Sumw2();
      
      hrapiditydiffMarco[i][j] = new TH1D(Form("%s_hrapiditydiffMarco_%s",prefix,suffix),Form("%s_rapiditydiffMarco_%s",prefix,suffix),80,-4,4);
      hrapiditydiffMarco[i][j]->GetXaxis()->SetTitle("|y_{t}|-|y_{#bar{t}}|");
      hrapiditydiffMarco[i][j]->Sumw2();

      hlepCosTheta[i][j] = new TH1D(Form("%s_hlepCosTheta_%s",prefix,suffix),Form("%s_lepCosTheta_%s",prefix,suffix),80,-1,1);
      hlepCosTheta[i][j]->GetXaxis()->SetTitle("cos(#theta_{l}^{t})");
      hlepCosTheta[i][j]->Sumw2();

      httpT[i][j] = new TH1D(Form("%s_httpT_%s",prefix,suffix),Form("%s_ttpT_%s",prefix,suffix),101,-4.,400.);
      httpT[i][j]->GetXaxis()->SetTitle("t#bar{t} p_{T} estimate (GeV/c)");
      httpT[i][j]->Sumw2();

      httpT_2d[i][j] = new TH2D(Form("%s_httpT2d_%s",prefix,suffix),Form("%s_ttpT2d_%s",prefix,suffix),101,-4.,400., 101,-4.,400.);
      httpT_2d[i][j]->GetXaxis()->SetTitle("t#bar{t} p_{T} gen (GeV/c)");
      httpT_2d[i][j]->GetYaxis()->SetTitle("t#bar{t} p_{T} estimate (GeV/c)");
      httpT_2d[i][j]->Sumw2();

      htopP_2d[i][j] = new TH2D(Form("%s_htopP2d_%s",prefix,suffix),Form("%s_topP2d_%s",prefix,suffix),100,0.,800., 100,0.,800.);
      htopP_2d[i][j]->GetXaxis()->SetTitle("top CM momentum gen (GeV/c)");
      htopP_2d[i][j]->GetYaxis()->SetTitle("top CM momentum estimate (GeV/c)");
      htopP_2d[i][j]->Sumw2();

      httMass[i][j] = new TH1D(Form("%s_httMass_%s",prefix,suffix),Form("%s_ttMass_%s",prefix,suffix),100,100.,1100.);
      httMass[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} estimate (GeV/c^{2})");
      httMass[i][j]->Sumw2();
      
      httMass_pull[i][j] = new TH1D(Form("%s_httMasspull_%s",prefix,suffix),Form("%s_ttMasspull_%s",prefix,suffix),200,-3,3);
      httMass_pull[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} (reco-gen)/gen");
      httMass_pull[i][j]->Sumw2();

      httMass_diff[i][j] = new TH1D(Form("%s_httMassdiff_%s",prefix,suffix),Form("%s_ttMassdiff_%s",prefix,suffix),100,-300,300);
      httMass_diff[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} (reco-gen)");
      httMass_diff[i][j]->Sumw2();

      htopMass_diff_plus[i][j] = new TH1D(Form("%s_htopMassdiff_plus_%s",prefix,suffix),Form("%s_topMassdiff_plus_%s",prefix,suffix),150,-150,150);
      htopMass_diff_plus[i][j]->GetXaxis()->SetTitle("M_{t} (reco-gen)");
      htopMass_diff_plus[i][j]->Sumw2();

      htopMass_diff_minus[i][j] = new TH1D(Form("%s_htopMassdiff_minus_%s",prefix,suffix),Form("%s_topMassdiff_minus_%s",prefix,suffix),150,-150,150);
      htopMass_diff_minus[i][j]->GetXaxis()->SetTitle("M_{#bar{t}} (reco-gen)");
      htopMass_diff_minus[i][j]->Sumw2();

      htopPCM_diff_plus[i][j] = new TH1D(Form("%s_htopPCMdiff_plus_%s",prefix,suffix),Form("%s_topPCMdiff_plus_%s",prefix,suffix),80,-300,300);
      htopPCM_diff_plus[i][j]->GetXaxis()->SetTitle("top momentum in CM (reco-gen)");
      htopPCM_diff_plus[i][j]->Sumw2();

      htopPCM_diff_minus[i][j] = new TH1D(Form("%s_htopPCMdiff_minus_%s",prefix,suffix),Form("%s_topPCMdiff_minus_%s",prefix,suffix),80,-300,300);
      htopPCM_diff_minus[i][j]->GetXaxis()->SetTitle("#bar{t} momentum in CM (reco-gen)");
      htopPCM_diff_minus[i][j]->Sumw2();

      httMass_2d[i][j] = new TH2D(Form("%s_httMass2d_%s",prefix,suffix),Form("%s_ttMass2d_%s",prefix,suffix),100,100.,1100., 100,100.,1100.);
      httMass_2d[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} gen (GeV/c^{2})");
      httMass_2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} estimate (GeV/c^{2})");
      httMass_2d[i][j]->Sumw2();

      htopMass[i][j] = new TH1D(Form("%s_htopMass_%s",prefix,suffix),Form("%s_topMass_%s",prefix,suffix),102,98,302);
      htopMass[i][j]->GetXaxis()->SetTitle("M_{t} estimate (GeV/c^{2})");
      htopMass[i][j]->Sumw2();

      httpT_nojetsmear[i][j] = new TH1D(Form("%s_httpT_nojetsmear_%s",prefix,suffix),Form("%s_ttpT_nojetsmear_%s",prefix,suffix),101,-4.,400.);
      httpT_nojetsmear[i][j]->GetXaxis()->SetTitle("t#bar{t} p_{T} estimate (GeV/c) (no jet smearing)");
      httpT_nojetsmear[i][j]->Sumw2();

      httpT_nojetsmear_2d[i][j] = new TH2D(Form("%s_httpT_nojetsmear2d_%s",prefix,suffix),Form("%s_ttpT_nojetsmear2d_%s",prefix,suffix),101,-4.,400., 101,-4.,400.);
      httpT_nojetsmear_2d[i][j]->GetXaxis()->SetTitle("t#bar{t} p_{T} gen (GeV/c)");
      httpT_nojetsmear_2d[i][j]->GetYaxis()->SetTitle("t#bar{t} p_{T} estimate (GeV/c) (no jet smearing)");
      httpT_nojetsmear_2d[i][j]->Sumw2();

      htopP_nojetsmear_2d[i][j] = new TH2D(Form("%s_htopP_nojetsmear2d_%s",prefix,suffix),Form("%s_topP_nojetsmear2d_%s",prefix,suffix),100,0.,800., 100,0.,800.);
      htopP_nojetsmear_2d[i][j]->GetXaxis()->SetTitle("top CM momentum gen (GeV/c)");
      htopP_nojetsmear_2d[i][j]->GetYaxis()->SetTitle("top CM momentum estimate (GeV/c) (no jet smearing)");
      htopP_nojetsmear_2d[i][j]->Sumw2();

      httMass_nojetsmear[i][j] = new TH1D(Form("%s_httMass_nojetsmear_%s",prefix,suffix),Form("%s_ttMass_nojetsmear_%s",prefix,suffix),100,100.,1100.);
      httMass_nojetsmear[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} estimate (GeV/c^{2}) (no jet smearing)");
      httMass_nojetsmear[i][j]->Sumw2();
      
      httMass_nojetsmear_pull[i][j] = new TH1D(Form("%s_httMass_nojetsmearpull_%s",prefix,suffix),Form("%s_ttMass_nojetsmearpull_%s",prefix,suffix),200,-3,3);
      httMass_nojetsmear_pull[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} (reco-gen)/gen (no jet smearing)");
      httMass_nojetsmear_pull[i][j]->Sumw2();

      httMass_nojetsmear_diff[i][j] = new TH1D(Form("%s_httMass_nojetsmeardiff_%s",prefix,suffix),Form("%s_ttMass_nojetsmeardiff_%s",prefix,suffix),100,-300,300);
      httMass_nojetsmear_diff[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} (reco-gen) (no jet smearing)");
      httMass_nojetsmear_diff[i][j]->Sumw2();

      httMass_nojetsmear_2d[i][j] = new TH2D(Form("%s_httMass_nojetsmear2d_%s",prefix,suffix),Form("%s_ttMass_nojetsmear2d_%s",prefix,suffix),100,100.,1100., 100,100.,1100.);
      httMass_nojetsmear_2d[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} gen (GeV/c^{2})");
      httMass_nojetsmear_2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} estimate (GeV/c^{2}) (no jet smearing)");
      httMass_nojetsmear_2d[i][j]->Sumw2();

      htopMass_nojetsmear_diff_plus[i][j] = new TH1D(Form("%s_htopMass_nojetsmeardiff_plus_%s",prefix,suffix),Form("%s_topMass_nojetsmeardiff_plus_%s",prefix,suffix),150,-150,150);
      htopMass_nojetsmear_diff_plus[i][j]->GetXaxis()->SetTitle("M_{t} (reco-gen) (no jet smearing)");
      htopMass_nojetsmear_diff_plus[i][j]->Sumw2();

      htopMass_nojetsmear_diff_minus[i][j] = new TH1D(Form("%s_htopMass_nojetsmeardiff_minus_%s",prefix,suffix),Form("%s_topMass_nojetsmeardiff_minus_%s",prefix,suffix),150,-150,150);
      htopMass_nojetsmear_diff_minus[i][j]->GetXaxis()->SetTitle("M_{#bar{t}} (reco-gen) (no jet smearing)");
      htopMass_nojetsmear_diff_minus[i][j]->Sumw2();

      htopMass_nojetsmear[i][j] = new TH1D(Form("%s_htopMass_nojetsmear_%s",prefix,suffix),Form("%s_topMass_nojetsmear_%s",prefix,suffix),100,0.,500.);
      htopMass_nojetsmear[i][j]->GetXaxis()->SetTitle("M_{t} estimate (GeV/c^{2}) (no jet smearing)");
      htopMass_nojetsmear[i][j]->Sumw2();



      httpT_gen[i][j] = new TH1D(Form("%s_httpTGen_%s",prefix,suffix),Form("%s_ttpTGen_%s",prefix,suffix),101,-4.,400.);
      httpT_gen[i][j]->GetXaxis()->SetTitle("t#bar{t} p_{T} gen (GeV/c)");
      httpT_gen[i][j]->Sumw2();

      httMass_gen[i][j] = new TH1D(Form("%s_httMassGen_%s",prefix,suffix),Form("%s_ttMassGen_%s",prefix,suffix),100,100.,1100.);
      httMass_gen[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} gen (GeV/c^{2})");
      httMass_gen[i][j]->Sumw2();

      htopMass_plus_gen[i][j] = new TH1D(Form("%s_htopMass_plusGen_%s",prefix,suffix),Form("%s_topMass_plusGen_%s",prefix,suffix),100,100.,250.);
      htopMass_plus_gen[i][j]->GetXaxis()->SetTitle("M_{t} gen (GeV/c^{2})");
      htopMass_plus_gen[i][j]->Sumw2();

      htopMass_minus_gen[i][j] = new TH1D(Form("%s_htopMass_minusGen_%s",prefix,suffix),Form("%s_topMass_minusGen_%s",prefix,suffix),100,100.,250.);
      htopMass_minus_gen[i][j]->GetXaxis()->SetTitle("M_{#bar{t}} gen (GeV/c^{2})");
      htopMass_minus_gen[i][j]->Sumw2();


      hllbbMass[i][j] = new TH1D(Form("%s_hllbbMass_%s",prefix,suffix),Form("%s_llbbMass_%s",prefix,suffix),100,100.,1100.);
      hllbbMass[i][j]->GetXaxis()->SetTitle("Mass(llbb) Estimate (GeV/c^{2})");
      hllbbMass[i][j]->Sumw2();
      

      httRapidity[i][j] = new TH1D(Form("%s_httRapidity_%s",prefix,suffix),Form("%s_ttRapidity_%s",prefix,suffix),120,-6,6);
      httRapidity[i][j]->GetXaxis()->SetTitle("y_{t} + y_{#bar{t}} estimate");
      httRapidity[i][j]->Sumw2();
      
      httRapidity2[i][j] = new TH1D(Form("%s_httRapidity2_%s",prefix,suffix),Form("%s_ttRapidity2_%s",prefix,suffix),120,-6,6);
      httRapidity2[i][j]->GetXaxis()->SetTitle("y_{t#bar{t}} estimate");
      httRapidity2[i][j]->Sumw2();
	
      hmassltb[i][j] =  new TH1D(Form("%s_hmassltb_%s",prefix,suffix),Form("%s_massltb_%s",prefix,suffix),75,0.,510.);
      hmassltb[i][j]->GetXaxis()->SetTitle("M_{l1b1} (GeV/c^{2})");
      hmassltb[i][j]->Sumw2();

      hmassllb[i][j] =  new TH1D(Form("%s_hmassllb_%s",prefix,suffix),Form("%s_massllb_%s",prefix,suffix),75,0.,510.);
      hmassllb[i][j]->GetXaxis()->SetTitle("M_{l2b2} (GeV/c^{2})");
      hmassllb[i][j]->Sumw2();
      
      hmassltb1Dmasscut[i][j] =  new TH1D(Form("%s_hmassltb1Dmasscut_%s",prefix,suffix),Form("%s_massltb1Dmasscut_%s",prefix,suffix),75,0.,510.);
      hmassltb1Dmasscut[i][j]->GetXaxis()->SetTitle("M_{l1b1} (GeV/c^{2}) for M_{l2b2} > 170 GeV/c^{2}");
      hmassltb1Dmasscut[i][j]->Sumw2();

      hmassllb1Dmasscut[i][j] =  new TH1D(Form("%s_hmassllb1Dmasscut_%s",prefix,suffix),Form("%s_massllb1Dmasscut_%s",prefix,suffix),75,0.,510.);
      hmassllb1Dmasscut[i][j]->GetXaxis()->SetTitle("M_{l2b2} (GeV/c^{2}) for M_{l1b1} > 170 GeV/c^{2}");
      hmassllb1Dmasscut[i][j]->Sumw2();

      htheSumJetPt[i][j] =  new TH1D(Form("%s_theSumJetPt_%s",prefix,suffix),Form("%s_theSumJetPt_%s",prefix,suffix),100,0.,600.);
      htheSumJetPt[i][j]->GetXaxis()->SetTitle("Sum of jet p_{T} (GeV/c)");
      htheSumJetPt[i][j]->Sumw2();

      htheSumBtagJetPt[i][j] =  new TH1D(Form("%s_theSumBtagJetPt_%s",prefix,suffix),Form("%s_theSumBtagJetPt_%s",prefix,suffix),100,0.,600.);
      htheSumBtagJetPt[i][j]->GetXaxis()->SetTitle("Sum of b tagged jet p_{T} (GeV/c)");
      htheSumBtagJetPt[i][j]->Sumw2();

      hthefirstJetPt[i][j] =  new TH1D(Form("%s_thefirstJetPt_%s",prefix,suffix),Form("%s_thefirstJetPt_%s",prefix,suffix),60,0.,360.);
      hthefirstJetPt[i][j]->GetXaxis()->SetTitle("First jet p_{T} (GeV/c)");
      hthefirstJetPt[i][j]->Sumw2();

      hthesecondJetPt[i][j] =  new TH1D(Form("%s_thesecondJetPt_%s",prefix,suffix),Form("%s_thesecondJetPt_%s",prefix,suffix),60,0.,360.);
      hthesecondJetPt[i][j]->GetXaxis()->SetTitle("Second jet p_{T} (GeV/c)");
      hthesecondJetPt[i][j]->Sumw2();
      
      htheleadinglepPt[i][j] = new TH1D(Form("%s_htheleadinglepPt_%s",prefix,suffix),Form("%s_theleadinglepPt_%s",prefix,suffix),60,0.,240.);
	  htheleadinglepPt[i][j]->GetXaxis()->SetTitle("Leading lepton p_{T} (GeV/c)");
      htheleadinglepPt[i][j]->Sumw2();

      hthesecondlepPt[i][j] = new TH1D(Form("%s_hthesecondlepPt_%s",prefix,suffix),Form("%s_thesecondlepPt_%s",prefix,suffix),60,0.,240.);
	  hthesecondlepPt[i][j]->GetXaxis()->SetTitle("Second lepton p_{T} (GeV/c)");
      hthesecondlepPt[i][j]->Sumw2();

      htheSumLepPt[i][j] =  new TH1D(Form("%s_theSumLepPt_%s",prefix,suffix),Form("%s_theSumLepPt_%s",prefix,suffix),120,0.,480.);
      htheSumLepPt[i][j]->GetXaxis()->SetTitle("Sum of lepton p_{T} (GeV/c)");
      htheSumLepPt[i][j]->Sumw2();

      hlepEta[i][j] = new TH1D(Form("%s_hlepEta_%s",prefix,suffix),Form("%s_lepEta_%s",prefix,suffix),60,-3.,3.);
      hlepEta[i][j]->GetXaxis()->SetTitle("Lepton #eta");
      hlepEta[i][j]->Sumw2();

      hlepPt[i][j] = new TH1D(Form("%s_hlepPt_%s",prefix,suffix),Form("%s_lepPt_%s",prefix,suffix),60,0.,240.);
      hlepPt[i][j]->GetXaxis()->SetTitle("Lepton p_{T} (GeV/c)");
      hlepPt[i][j]->Sumw2();
      
      hjetPt[i][j] = new TH1D(Form("%s_hjetPt_%s",prefix,suffix),Form("%s_jetPt_%s",prefix,suffix),60,0.,360.);
      hjetPt[i][j]->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      hjetPt[i][j]->Sumw2();

      hjetEta[i][j] = new TH1D(Form("%s_hjetEta_%s",prefix,suffix),Form("%s_jetEta_%s",prefix,suffix),60,-3.,3.);
      hjetEta[i][j]->GetXaxis()->SetTitle("Jet #eta");
      hjetEta[i][j]->Sumw2();
      
     
      hMET[i][j] = new TH1D(Form("%s_hMET_%s",prefix,suffix),Form("%s_MET_%s",prefix,suffix),100,0.,400.);
      hMET[i][j]->GetXaxis()->SetTitle("E_{T}^{miss} (GeV)");
      hMET[i][j]->Sumw2();

       
      hmasslb_2d[i][j] = new TH2D(Form("%s_hmasslb2d_%s",prefix,suffix),  Form("%s_masslb2d_%s" ,prefix,suffix),100,0,500,100,0,500);
      hmasslb_2d[i][j]->GetXaxis()->SetTitle("M_{l2b2}(GeV/c^{2})");
      hmasslb_2d[i][j]->GetYaxis()->SetTitle("M_{l1b1}(GeV/c^{2})");
      hmasslb_2d[i][j]->Sumw2();

      
      habcd_2d[i][j] = new TH2D(Form("%s_habcd2d_%s",prefix,suffix), Form("%s_habcd2d_%s" ,prefix,suffix),100,0,500,100,0,500);
      habcd_2d[i][j]->GetXaxis()->SetTitle("M_{l2b2}(GeV/c^{2})");
      habcd_2d[i][j]->GetYaxis()->SetTitle("M_{l1b1}(GeV/c^{2})");
      habcd_2d[i][j]->Sumw2();


      httmasssm_2d[i][j] = new TH2D(Form("%s_httmasssm2d_%s",prefix,suffix), Form("%s_httmasssm2d_%s" ,prefix,suffix),100,100,1100,100,100,1100);
      httmasssm_2d[i][j]->GetXaxis()->SetTitle("M_{t#bar{t},1} (GeV/c^{2})");
      httmasssm_2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t},100} (GeV/c^{2})");
      httmasssm_2d[i][j]->Sumw2();
      htopmasssm_2d[i][j] = new TH2D(Form("%s_htopmasssm2d_%s",prefix,suffix), Form("%s_htopmasssm2d_%s" ,prefix,suffix),102,98,302,102,98,302);
      htopmasssm_2d[i][j]->GetXaxis()->SetTitle("M_{t,1} (GeV/c^{2})");
      htopmasssm_2d[i][j]->GetYaxis()->SetTitle("M_{t,100} (GeV/c^{2})");
      htopmasssm_2d[i][j]->Sumw2();

      htop1pCMsm_2d[i][j] = new TH2D(Form("%s_htop1pCMsm2d_%s",prefix,suffix), Form("%s_htop1pCMsm2d_%s" ,prefix,suffix),101,-8,800,101,-8,800);
      htop1pCMsm_2d[i][j]->GetXaxis()->SetTitle("p t CM,1 (GeV/c)");
      htop1pCMsm_2d[i][j]->GetYaxis()->SetTitle("p t CM,100 (GeV/c)");
      htop1pCMsm_2d[i][j]->Sumw2();
      htop2pCMsm_2d[i][j] = new TH2D(Form("%s_htop2pCMsm2d_%s",prefix,suffix), Form("%s_htop2pCMsm2d_%s" ,prefix,suffix),101,-8,800,101,-8,800);
      htop2pCMsm_2d[i][j]->GetXaxis()->SetTitle("p #bar{t} CM,1 (GeV/c)");
      htop2pCMsm_2d[i][j]->GetYaxis()->SetTitle("p #bar{t} CM,100 (GeV/c)");
      htop2pCMsm_2d[i][j]->Sumw2();

      htop1pTsm_2d[i][j] = new TH2D(Form("%s_htop1pTsm2d_%s",prefix,suffix), Form("%s_htop1pTsm2d_%s" ,prefix,suffix),101,-8,800,101,-8,800);
      htop1pTsm_2d[i][j]->GetXaxis()->SetTitle("p_{T} t lab,1 (GeV/c)");
      htop1pTsm_2d[i][j]->GetYaxis()->SetTitle("p_{T} t lab,100 (GeV/c)");
      htop1pTsm_2d[i][j]->Sumw2();
      //htop1etasm_2d[i][j] = new TH2D(Form("%s_htop1etasm2d_%s",prefix,suffix), Form("%s_htop1etasm2d_%s" ,prefix,suffix),100,-4,4,100,-4,4);
      //htop1etasm_2d[i][j]->GetXaxis()->SetTitle("#eta t lab,1");
      //htop1etasm_2d[i][j]->GetYaxis()->SetTitle("#eta t lab,100");
      //htop1etasm_2d[i][j]->Sumw2();
      //htop1phism_2d[i][j] = new TH2D(Form("%s_htop1phism2d_%s",prefix,suffix), Form("%s_htop1phism2d_%s" ,prefix,suffix),100,-3.141592653589793,3.141592653589793,100,-3.141592653589793,3.141592653589793);
      //htop1phism_2d[i][j]->GetXaxis()->SetTitle("#phi t lab,1 (GeV/c)");
      //htop1phism_2d[i][j]->GetYaxis()->SetTitle("#phi t lab,100 (GeV/c)");
      //htop1phism_2d[i][j]->Sumw2();
      htop2pTsm_2d[i][j] = new TH2D(Form("%s_htop2pTsm2d_%s",prefix,suffix), Form("%s_htop2pTsm2d_%s" ,prefix,suffix),101,-8,800,101,-8,800);
      htop2pTsm_2d[i][j]->GetXaxis()->SetTitle("p_{T} #bar{t} lab,1 (GeV/c)");
      htop2pTsm_2d[i][j]->GetYaxis()->SetTitle("p_{T} #bar{t} lab,100 (GeV/c)");
      htop2pTsm_2d[i][j]->Sumw2();
      //htop2etasm_2d[i][j] = new TH2D(Form("%s_htop2etasm2d_%s",prefix,suffix), Form("%s_htop2etasm2d_%s" ,prefix,suffix),100,-4,4,100,-4,4);
      //htop2etasm_2d[i][j]->GetXaxis()->SetTitle("#eta #bar{t} lab,1");
      //htop2etasm_2d[i][j]->GetYaxis()->SetTitle("#eta #bar{t} lab,100");
      //htop2etasm_2d[i][j]->Sumw2();
      //htop2phism_2d[i][j] = new TH2D(Form("%s_htop2phism2d_%s",prefix,suffix), Form("%s_htop2phism2d_%s" ,prefix,suffix),100,-3.141592653589793,3.141592653589793,100,-3.141592653589793,3.141592653589793);
      //htop2phism_2d[i][j]->GetXaxis()->SetTitle("#phi #bar{t} lab,1");
      //htop2phism_2d[i][j]->GetYaxis()->SetTitle("#phi #bar{t} lab,100");
      //htop2phism_2d[i][j]->Sumw2();

      // generator level distributions

      httMassGluongenp[i][j] = new TH1D(Form("%s_httMassGluongenp_%s",prefix,suffix),Form("%s_ttMassGluongenp_%s",prefix,suffix),200, 0 , 1000);
      httMassGluongenp[i][j]->GetXaxis()->SetTitle("ttMassGluongenp ");
      httMassGluongenp[i][j]->Sumw2();
      
      httMassQuarkgenp[i][j] = new TH1D(Form("%s_httMassQuarkgenp_%s",prefix,suffix),Form("%s_ttMassQuarkgenp_%s",prefix,suffix),200, 0 , 1000);
      httMassQuarkgenp[i][j]->GetXaxis()->SetTitle("ttMassQuarkgenp ");
      httMassQuarkgenp[i][j]->Sumw2();

      httRapidityGluongenp[i][j] = new TH1D(Form("%s_httRapidityGluongenp_%s",prefix,suffix),Form("%s_ttRapidityGluongenp_%s",prefix,suffix),200, -10, 10);
      httRapidityGluongenp[i][j]->GetXaxis()->SetTitle("ttRapidityGluongenp ");
      httRapidityGluongenp[i][j]->Sumw2();
      
      httRapidityQuarkgenp[i][j] = new TH1D(Form("%s_httRapidityQuarkgenp_%s",prefix,suffix),Form("%s_ttRapidityQuarkgenp_%s",prefix,suffix),200,-10, 10);
      httRapidityQuarkgenp[i][j]->GetXaxis()->SetTitle("ttRapidityQuarkgenp ");
      httRapidityQuarkgenp[i][j]->Sumw2();

      hllbbRapidityGluongenp[i][j] = new TH1D(Form("%s_hllbbRapidityGluongenp_%s",prefix,suffix),Form("%s_llbbRapidityGluongenp_%s",prefix,suffix),200, -10, 10);
      hllbbRapidityGluongenp[i][j]->GetXaxis()->SetTitle("llbbRapidityGluongenp ");
      hllbbRapidityGluongenp[i][j]->Sumw2();
      
      hllbbRapidityQuarkgenp[i][j] = new TH1D(Form("%s_hllbbRapidityQuarkgenp_%s",prefix,suffix),Form("%s_llbbRapidityQuarkgenp_%s",prefix,suffix),200,-10, 10);
      hllbbRapidityQuarkgenp[i][j]->GetXaxis()->SetTitle("llbbRapidityQuarkgenp ");
      hllbbRapidityQuarkgenp[i][j]->Sumw2();


      htheSumBtagJetPtgenp[i][j] =  new TH1D(Form("%s_theSumBtagJetPtgenp_%s",prefix,suffix),Form("%s_theSumBtagJetPtgenp_%s",prefix,suffix),100,0.,500.);
      htheSumBtagJetPtgenp[i][j]->GetXaxis()->SetTitle("SumBtagJetPtgenp (GeV/c)");
      htheSumBtagJetPtgenp[i][j]->Sumw2();

      hthefirstBtagJetPtgenp[i][j] =  new TH1D(Form("%s_thefirstBtagJetPtgenp_%s",prefix,suffix),Form("%s_thefirstBtagJetPtgenp_%s",prefix,suffix),100,0.,500.);
      hthefirstBtagJetPtgenp[i][j]->GetXaxis()->SetTitle("firstBtagJetPtgenp (GeV/c)");
      hthefirstBtagJetPtgenp[i][j]->Sumw2();

      hthesecondBtagJetPtgenp[i][j] =  new TH1D(Form("%s_thesecondBtagJetPtgenp_%s",prefix,suffix),Form("%s_thesecondBtagJetPtgenp_%s",prefix,suffix),100,0.,500.);
      hthesecondBtagJetPtgenp[i][j]->GetXaxis()->SetTitle("secondBtagJetPtgenp (GeV/c)");
      hthesecondBtagJetPtgenp[i][j]->Sumw2();
      
      htheleadinglepPtgenp[i][j] = new TH1D(Form("%s_htheleadinglepPtgenp_%s",prefix,suffix),Form("%s_theleadinglepPtgenp_%s",prefix,suffix),50,0.,250.);
      htheleadinglepPtgenp[i][j]->GetXaxis()->SetTitle("leadingLeptongenp p_{T} (GeV/c)");
      htheleadinglepPtgenp[i][j]->Sumw2();

      hthesecondlepPtgenp[i][j] = new TH1D(Form("%s_hthesecondlepPtgenp_%s",prefix,suffix),Form("%s_thesecondlepPtgenp_%s",prefix,suffix),50,0.,250.);
      hthesecondlepPtgenp[i][j]->GetXaxis()->SetTitle("secondLeptongenp p_{T} (GeV/c)");
      hthesecondlepPtgenp[i][j]->Sumw2();

      htheSumLepPtgenp[i][j] =  new TH1D(Form("%s_theSumLepPtgenp_%s",prefix,suffix),Form("%s_theSumLepPtgenp_%s",prefix,suffix),100,0.,500.);
      htheSumLepPtgenp[i][j]->GetXaxis()->SetTitle("SumLepPtgenp (GeV/c)");
      htheSumLepPtgenp[i][j]->Sumw2();

      htheleadingNuPtgenp[i][j] = new TH1D(Form("%s_htheleadingNuPtgenp_%s",prefix,suffix),Form("%s_theleadingNuPtgenp_%s",prefix,suffix),50,0.,250.);
      htheleadingNuPtgenp[i][j]->GetXaxis()->SetTitle("leadingNugenp p_{T} (GeV/c)");
      htheleadingNuPtgenp[i][j]->Sumw2();

      hthesecondNuPtgenp[i][j] = new TH1D(Form("%s_hthesecondNuPtgenp_%s",prefix,suffix),Form("%s_thesecondNuPtgenp_%s",prefix,suffix),50,0.,250.);
      hthesecondNuPtgenp[i][j]->GetXaxis()->SetTitle("secondNugenp p_{T} (GeV/c)");
      hthesecondNuPtgenp[i][j]->Sumw2();

      hMETgenp[i][j] = new TH1D(Form("%s_hMETgenp_%s",prefix,suffix),Form("%s_METgenp_%s",prefix,suffix),50,0.,250.);
      hMETgenp[i][j]->GetXaxis()->SetTitle("METgenp (GeV)");
      hMETgenp[i][j]->Sumw2();

      htheSumLBPtgenp[i][j] =  new TH1D(Form("%s_theSumLBPtgenp_%s",prefix,suffix),Form("%s_theSumLBPtgenp_%s",prefix,suffix),100,0.,500.);
      htheSumLBPtgenp[i][j]->GetXaxis()->SetTitle("SumLBPtgenp (GeV/c)");
      htheSumLBPtgenp[i][j]->Sumw2();

      htheSumLBNPtgenp[i][j] =  new TH1D(Form("%s_theSumLBNPtgenp_%s",prefix,suffix),Form("%s_theSumLBNPtgenp_%s",prefix,suffix),100,0.,500.);
      htheSumLBNPtgenp[i][j]->GetXaxis()->SetTitle("SumLBNPtgenp (GeV/c)");
      htheSumLBNPtgenp[i][j]->Sumw2();


      hdRlbtruegenp[i][j] = new TH1D(Form("%s_hdRlbtruegenp_%s",prefix,suffix),Form("%s_dRlbtruegenp_%s",prefix,suffix),35,0.,3.5);
      hdRlbtruegenp[i][j]->GetXaxis()->SetTitle("dRlbtruegenp ");
      hdRlbtruegenp[i][j]->Sumw2();

      
      hdRlbfalsegenp[i][j] = new TH1D(Form("%s_hdRlbfalsegenp_%s",prefix,suffix),Form("%s_dRlbfalsegenp_%s",prefix,suffix),35,0.,3.5);
      hdRlbfalsegenp[i][j]->GetXaxis()->SetTitle("dRlbfalsegenp ");
      hdRlbfalsegenp[i][j]->Sumw2();

      hdRlbratiogenp[i][j] = new TH1D(Form("%s_hdRlbratiogenp_%s",prefix,suffix),Form("%s_dRlbratiogenp_%s",prefix,suffix),40,-2.,2);
      hdRlbratiogenp[i][j]->GetXaxis()->SetTitle("dRlbratiogenp ");
      hdRlbratiogenp[i][j]->Sumw2();
      
      htopptgenp[i][j] = new TH1D(Form("%s_htopptgenp_%s",prefix,suffix),Form("%s_topptgenp_%s",prefix,suffix),100, 0 , 500);
      htopptgenp[i][j]->GetXaxis()->SetTitle("topptgenp ");
      htopptgenp[i][j]->Sumw2();

      htopMassgenp[i][j] = new TH1D(Form("%s_htopMassgenp_%s",prefix,suffix),Form("%s_topMassgenp_%s",prefix,suffix),500, 0 , 500);
      htopMassgenp[i][j]->GetXaxis()->SetTitle("topMassgenp ");
      htopMassgenp[i][j]->Sumw2();
      
      htopptdrgenp_2d[i][j] = new TH2D(Form("%s_htopptdrgenp2d_%s",prefix,suffix), Form("%s_htopptdrgenp2d_%s" ,prefix,suffix),100,0,500,35,0,3.5);
      htopptdrgenp_2d[i][j]->GetXaxis()->SetTitle("top p_{T} (GeV/c)");
      htopptdrgenp_2d[i][j]->GetYaxis()->SetTitle("DR_lb");
      htopptdrgenp_2d[i][j]->Sumw2();
    
        
      hmasslbgenp_2d[i][j] = new TH2D(Form("%s_hmasslbgenp2d_%s",prefix,suffix), Form("%s_masslbgenp2d_%s" ,prefix,suffix),100,0,500,100,0,500);
      hmasslbgenp_2d[i][j]->GetXaxis()->SetTitle("Gen M_{l2b2}(GeV/c^{2})");
      hmasslbgenp_2d[i][j]->GetYaxis()->SetTitle("Gen M_{l1b1}(GeV/c^{2})");
      hmasslbgenp_2d[i][j]->Sumw2();

      hmasslbgenmatch1_2d[i][j] = new TH2D(Form("%s_hmasslbgenmatch12d_%s",prefix,suffix), Form("%s_masslbgenmatch12d_%s" ,prefix,suffix),100,0,500,100,0,500);
      hmasslbgenmatch1_2d[i][j]->GetXaxis()->SetTitle("Gen M_{l2b2}(GeV/c^{2})");
      hmasslbgenmatch1_2d[i][j]->GetYaxis()->SetTitle("Gen M_{l1b1}(GeV/c^{2})");
      hmasslbgenmatch1_2d[i][j]->Sumw2();
     
      hmasslbgenmatch_2d[i][j] = new TH2D(Form("%s_hmasslbgenmatch2d_%s",prefix,suffix), Form("%s_masslbgenmatch2d_%s" ,prefix,suffix),100,0,500,100,0,500);
      hmasslbgenmatch_2d[i][j]->GetXaxis()->SetTitle("Gen M_{l2b2}(GeV/c^{2})");
      hmasslbgenmatch_2d[i][j]->GetYaxis()->SetTitle("Gen M_{l1b1}(GeV/c^{2})");
      hmasslbgenmatch_2d[i][j]->Sumw2();
      
      
      //daughter lepton angle in tau rest frame to check if MC is correctly using the tau polarisation
      hlepPlusCosThetaTau_gen[i][j] = new TH1D(Form("%s_hlepPlusCosThetaTauGen_%s",prefix,suffix),Form("%s_lepPlusCosThetaTauGen_%s",prefix,suffix),80,-1,1);
      hlepPlusCosThetaTau_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{#tau})");
      hlepPlusCosThetaTau_gen[i][j]->Sumw2();
      
      hlepMinusCosThetaTau_gen[i][j] = new TH1D(Form("%s_hlepMinusCosThetaTauGen_%s",prefix,suffix),Form("%s_lepMinusCosThetaTauGen_%s",prefix,suffix),80,-1,1);
      hlepMinusCosThetaTau_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#tau})");
      hlepMinusCosThetaTau_gen[i][j]->Sumw2();

      //daughter lepton E/Emax in tau rest frame to check if MC is correctly using the tau polarisation
      hlepPlusxTau_gen[i][j] = new TH1D(Form("%s_hlepPlusxTauGen_%s",prefix,suffix),Form("%s_lepPlusxTauGen_%s",prefix,suffix),80,0,1);
      hlepPlusxTau_gen[i][j]->GetXaxis()->SetTitle("x");
      hlepPlusxTau_gen[i][j]->Sumw2();
      
      hlepMinusxTau_gen[i][j] = new TH1D(Form("%s_hlepMinusxTauGen_%s",prefix,suffix),Form("%s_lepMinusxTauGen_%s",prefix,suffix),80,0,1);
      hlepMinusxTau_gen[i][j]->GetXaxis()->SetTitle("x");
      hlepMinusxTau_gen[i][j]->Sumw2();

            
      //DYEst histos
      hdilMassWithMetDYEst[i][j] = new TH1D(Form("%s_hdilMassWithMetDYEst_%s",  prefix,suffix), "Di-lepton mass with MET for DY Estimation", 40, 0., 200.);
      hdilMassWithMetDYEst[i][j]->GetXaxis()->SetTitle("M_{ll}(GeV/c^{2}), with MET > 30 GeV cut");
      hdilMassWithMetDYEst[i][j]->Sumw2();

      hdilMassNoMetDYEst[i][j] = new TH1D(Form("%s_hdilMassNoMetDYEst_%s",  prefix,suffix), "Di-lepton mass without MET for DY Estimation", 40, 0., 200.);
      hdilMassNoMetDYEst[i][j]->GetXaxis()->SetTitle("M_{ll}(GeV/c^{2}), no MET cut");
      hdilMassNoMetDYEst[i][j]->Sumw2();

      hmetInDYEst[i][j] = new TH1D(Form("%s_hmetInDYEst_%s",  prefix,suffix), "MET in Z mass for DY Estimation", 40, 0., 200.);
      hmetInDYEst[i][j]->GetXaxis()->SetTitle("MET (GeV), inside Z mass window");
      hmetInDYEst[i][j]->Sumw2();

      hmetOutDYEst[i][j] = new TH1D(Form("%s_hmetOutDYEst_%s",  prefix,suffix), "MET outside Z mass for DY Estimation", 40, 0., 200.);
      hmetOutDYEst[i][j]->GetXaxis()->SetTitle("MET (GeV), outside Z mass window");
      hmetOutDYEst[i][j]->Sumw2();
    }
    
  }

 

}
