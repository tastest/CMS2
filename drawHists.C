#include <iostream>
#include <vector>
#include "TStyle.h"
//#include "histtools.C"
#include "browseStacks.C"
#include "getMyHistosNames.C"
#include "TSystem.h"
#include "tdrStyle.C"
#include "DYest.C"
#include "CommonFunctions.C"
#include "EstimateSignalSpillage.C"

/*                                                                                                                                                                                                                                                            
  replacewFRest = replace the yields and errors for WJets, QCD with the FR estimate results                                                                                                                                                                   
  Does the Spillage subtraction as well                                                                                                                                                                                                                       
  replacewDYest = replace the yields and errors for DY with the data driven estimate                                                                                                                                                                          
  drawFullErrors = draw the errors in all their glory. Includes the systematic errors:                                                                                                                                                                        
  50% on the WJets estimate                                                                                                                                                                                                                                   
  100% on the QCD estimate                                                                                                                                                                                                                                    
  50% on DY estimate                                                                                                                                                                                                                                          
  50% on the backgrounds from MC                                                                                                                                                                                                                              
*/

//void makePSFile(const TString dataFName="results_data/hist_usePtGt2020_applyTriggers_hypDisamb_usejptJets_usetcMET_requireEcalEls_useOS_vetoHypMassLt12_require2BTag_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root", const TString mcFName="results/hist_usePtGt2020_applyTriggers_hypDisamb_usejptJets_usetcMET_requireEcalEls_useOS_vetoHypMassLt12_require2BTag_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root",         float scaleMC=36.1, bool drawLogY = true, bool drawDiffs = false,

void makePSFile(const TString FName="results/hist_usePtGt2020_applyTriggers_hypDisamb_usepfMET_usepfJets_requireEcalEls_useOS_vetoHypMassLt12_require2BTag_sortJetCandidatesbyDR_applylepIDCuts_applylepIsoCuts_vetoZmass_veto2Jets_vetoMET.root",     
		int rebin=1, bool drawtprime = false,
		const TString pat="default", 
		float scaleMC=4.98, bool drawLogY = false, bool allhypOnly = false, bool drawDiffs = false,
                bool replacewFRest = false, bool replacewDYest = false,
                bool drawFullErrors = false, bool scalebTagMCbyDataDrivenEsts = false) {

  setTDRStyle();

  // hist::loadHist(FName.Data(),0,"*hnJet*");
  
  if(pat!="default") hist::loadHist(FName.Data(),0,pat.Data());
  else {

  hist::loadHist(FName.Data(),0,"*hnJet_allj*");
  hist::loadHist(FName.Data(),0,"*hnVtx_allj*");
  hist::loadHist(FName.Data(),0,"*hnBtagJet_allj*");

  hist::loadHist(FName.Data(),0,"*hlepPt_allj*");
  //hist::loadHist(FName.Data(),0,"*hlepPlusPt_allj*");
  //hist::loadHist(FName.Data(),0,"*hlepMinusPt_allj*");
  hist::loadHist(FName.Data(),0,"*hlepEta_allj*");
  //hist::loadHist(FName.Data(),0,"*hlepPlusEta_allj*");
  //hist::loadHist(FName.Data(),0,"*hlepMinusEta_allj*");
  hist::loadHist(FName.Data(),0,"*theSumLepPt_allj*");
  //hist::loadHist(FName.Data(),0,"*hlepPhi_allj*");
  //hist::loadHist(FName.Data(),0,"*hlepPlusPhi_allj*");
  //hist::loadHist(FName.Data(),0,"*hlepMinusPhi_allj*");
  //hist::loadHist(FName.Data(),0,"*hlepAzimAsym_allj*");
  
  hist::loadHist(FName.Data(),0,"*hjetPt_allj*");
  hist::loadHist(FName.Data(),0,"*hjetEta_allj*");
  //hist::loadHist(FName.Data(),0,"*hjetPhi_allj*");
  //hist::loadHist(FName.Data(),0,"*thefirstJetPt_allj*");
  //hist::loadHist(FName.Data(),0,"*thesecondJetPt_allj*");
  hist::loadHist(FName.Data(),0,"*theSumJetPt_allj*");
  //hist::loadHist(FName.Data(),0,"*theSumBtagJetPt_allj*");
 
  //hist::loadHist(FName.Data(),0,"*htheleadinglepPt_allj*");
  //hist::loadHist(FName.Data(),0,"*hthesecondlepPt_allj*");

  hist::loadHist(FName.Data(),0,"*hMET_allj*");
  //hist::loadHist(FName.Data(),0,"*hmassltb_allj*");
  //hist::loadHist(FName.Data(),0,"*hmassllb_allj*");

  hist::loadHist(FName.Data(),0,"*hNsolns_allj*");
  hist::loadHist(FName.Data(),0,"*hmaxAMWTweight_allj*");
  hist::loadHist(FName.Data(),0,"*hsumAMWTweight_allj*");
  hist::loadHist(FName.Data(),0,"*haveAMWTweight_allj*");
  //hist::loadHist(FName.Data(),0,"*hAMWTweight_nojetsmear_allj*");

  hist::loadHist(FName.Data(),0,"*htopMass_allj*");
  hist::loadHist(FName.Data(),0,"*httMass_allj*"); 
  hist::loadHist(FName.Data(),0,"*httpT_allj*"); 
  hist::loadHist(FName.Data(),0,"*httRapidity_allj*");   
  hist::loadHist(FName.Data(),0,"*httRapidity2_allj*");   

  hist::loadHist(FName.Data(),0,"*hlepChargeAsym_allj*");
  hist::loadHist(FName.Data(),0,"*hlepAzimAsym_allj*");
  hist::loadHist(FName.Data(),0,"*hlepAzimAsym2_allj*");
  hist::loadHist(FName.Data(),0,"*htopCosTheta_allj*");
  hist::loadHist(FName.Data(),0,"*hpseudorapiditydiff_allj*");
  hist::loadHist(FName.Data(),0,"*hrapiditydiff_allj*");
  hist::loadHist(FName.Data(),0,"*hrapiditydiffMarco_allj*");
  hist::loadHist(FName.Data(),0,"*hlepPlusCosTheta_allj*");
  hist::loadHist(FName.Data(),0,"*hlepMinusCosTheta_allj*");
  hist::loadHist(FName.Data(),0,"*hlepCosTheta_allj*");
  hist::loadHist(FName.Data(),0,"*htopSpinCorr_allj*");

  //hist::loadHist(FName.Data(),0,"*hpseudorapiditydiff2_allj*");
  //hist::loadHist(FName.Data(),0,"*hrapiditydiff2_allj*");

  //hist::loadHist(FName.Data(),0,"*hlepRapDiff_allj*");
  //hist::loadHist(FName.Data(),0,"*hlepAngleBetween_allj*");
  //hist::loadHist(FName.Data(),0,"*hlepAngleBetweenCMS_allj*");  
  
  //hist::loadHist(FName.Data(),0,"*hjetAzimAsym_allj*");
  //hist::loadHist(FName.Data(),0,"*hjetRapDiff_allj*");
  //hist::loadHist(FName.Data(),0,"*hjetAngleBetween_allj*");
  //hist::loadHist(FName.Data(),0,"*hjetAngleBetweenCMS_allj*");
    
  //  hist::loadHist(FName.Data(),0,"*htopCosThetaGen_allj*");
  //  hist::loadHist(FName.Data(),0,"*hlepChargeAsymGen_allj*");
  //  hist::loadHist(FName.Data(),0,"*hlepCosThetaGen_allj*");
  //  hist::loadHist(FName.Data(),0,"*htopSpinCorrGen_allj*");
  }
  //if(drawtprime && ! drawLogY) hist::scale("ttprime*", 10);

  //hist::scale("ttdil_*", 1.+(9824. - 10086.77)/9344.25); //scale mc@nlo to data
  hist::scale("ttdil_*", 1.+(9824. - 10063.47)/9323.84); //scale mc@nlo to data, after pT reweighting


  //if(scalebTagMCbyDataDrivenEsts) {                                                                                                                                                                                                
  vector<TString> v_samples;
  vector<Color_t> v_colors;
  vector<TString> v_legEntries;
  vector<Style_t> v_styles;

  //the order you put the sample into the vector is the order                                                                                                                                                                                                 
  //that the histogram gets drawn. So if you put ttdil in first,                                                                                                                                                                                              
  //it will be drawn first, i.e, it will be drawn on the bottom                                                                                                                                                                                               
  //of the stack!                                                                                                                                                                                                                                             

  
  
  if(drawLogY) {
    v_samples.push_back("ttdil");
    v_colors.push_back(kRed+1);
    v_legEntries.push_back("t#bar{t}");
    v_styles.push_back(1001);
    
  }
  
   
  if(!drawLogY) {
  v_samples.push_back("wjets");
  v_colors.push_back(kGreen-3);
  if(!scalebTagMCbyDataDrivenEsts)
    v_legEntries.push_back("W#rightarrowl#nu");
  else
    v_legEntries.push_back("Non-W/Z prediction");
  v_styles.push_back(1001);    
  }


  //  v_samples.push_back("FR");
  // v_colors.push_back(kGreen-3);
  //v_legEntries.push_back("Non-W/Z prediction");
  //v_styles.push_back(1001);


  v_samples.push_back("VV");
  v_colors.push_back(kYellow-10);
  v_legEntries.push_back("VV");
  v_styles.push_back(1001);

  v_samples.push_back("tw");
  v_colors.push_back(kMagenta);
  v_legEntries.push_back("Single top");
  v_styles.push_back(1001);


  
  v_samples.push_back("DYtautau");
  v_colors.push_back(kAzure+8);
  v_legEntries.push_back("Z/#gamma*#rightarrow#tau^{+}#tau^{-}");
  v_styles.push_back(1001);


  std::vector<TString> v_prfxsToCombine;
  v_prfxsToCombine.clear();
  v_prfxsToCombine.push_back("DYee");
  v_prfxsToCombine.push_back("DYmm");
  hist::combineHists(v_prfxsToCombine, "DYeemm");


    
  v_samples.push_back("DYeemm");
  v_colors.push_back(kAzure - 2);
  //if(replacewDYest || scalebTagMCbyDataDrivenEsts)                                                                                                                                                                                                          
  //  v_legEntries.push_back("Z/#gamma*#rightarrowl^{+}l^{-} prediction");                                                                                                                                                                                    
  //else                                                                                                                                                                                                                                                      
  v_legEntries.push_back("Z/#gamma*#rightarrowl^{+}l^{-}");
  v_styles.push_back(1001);

  v_samples.push_back("ttotr");
  v_colors.push_back(kRed-1);
  v_legEntries.push_back("t#bar{t} (other)");
  v_styles.push_back(1001);
  
  
  if(!drawLogY) {
    v_samples.push_back("ttdil");
    v_colors.push_back(kRed+1);
    v_legEntries.push_back("t#bar{t} (dilepton)");
    v_styles.push_back(1001);
    
  }
    

  if(drawtprime) {
    v_samples.push_back("wprime400");
    v_colors.push_back(kBlue-2);
    if(drawLogY) v_legEntries.push_back("W' 400 GeV");
    else v_legEntries.push_back("W' 400 GeV");
    v_styles.push_back(0);
  }

  
 
 
  //data should be the last thing put into the vectors                                                                                                                                                                                                        
  v_samples.push_back("data");
  v_colors.push_back(kBlack);
  v_legEntries.push_back("Data");
  v_styles.push_back(20);
 

  
  browseStacks(v_samples, v_colors, FName, v_legEntries, drawLogY, v_styles, drawFullErrors, drawDiffs, scaleMC, rebin, allhypOnly);

  hist::deleteHistos();

}
  
