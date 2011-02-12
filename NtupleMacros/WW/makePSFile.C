#include <iostream>
#include <vector>
#include "TStyle.h"
#include "histtools.C"
#include "TSystem.h"
#include "tdrstyle.C"

/*
  drawFullErrors = draw the errors in all their glory. Includes the systematic errors:
                   50% on the WJets estimate
		   100% on the QCD estimate
		   50% on DY estimate
		   50% on the backgrounds from MC
*/

bool sampleIsPresent(const char* pattern){
  TRegexp re(pattern, kTRUE);
  TIter nextobj(gDirectory->GetList());
  while ( TObject *obj = (TObject *) nextobj() ) {
    TString s = obj->GetName();
    if (s.Index(re) == kNPOS) continue;
    return true;
  }
  cout << "No histogram for sample " << pattern << " is found" << endl;
  return false;
}

void makePSFile(const TString dataFName,
		float scaleMC=0.355 /*with respect to 100/pb*/,  
		bool drawLogY = true, bool drawFullErrors = false) {
  
  setTDRStyle();

  hist::loadHist(dataFName.Data(), 0, "*_hmetVal_*");
  hist::loadHist(dataFName.Data(), 0, "*_hmetProjVal_*");
  hist::loadHist(dataFName.Data(), 0, "*_hmaxPFJetPtVal_*");
  hist::loadHist(dataFName.Data(), 0, "*_hdilMassVal_*");

  hist::loadHist(dataFName.Data(), 0, "*_hmetNM1_*");
  hist::loadHist(dataFName.Data(), 0, "*_hmetProjNM1_*");
  hist::loadHist(dataFName.Data(), 0, "*_hmaxPFJetPtNM1_*");
  hist::loadHist(dataFName.Data(), 0, "*_hdilMassNM1_*");

  hist::loadHist(dataFName.Data(), 0, "*_hmet_*");
  hist::loadHist(dataFName.Data(), 0, "*_hmetProj_*");
  hist::loadHist(dataFName.Data(), 0, "*_hdilMass_*");

  /*
  hist::scale("ww_*", scaleMC);     
  hist::scale("dyee_*", scaleMC);     
  hist::scale("dymm_*", scaleMC);     
  hist::scale("dytt_*", scaleMC);     
  hist::scale("ttbar_*", scaleMC);     
  hist::scale("tw_*", scaleMC);     
  hist::scale("wjets_*", scaleMC);     
  hist::scale("wz_*", scaleMC);     
  hist::scale("zz_*", scaleMC);     
  */
  vector<TString> v_prfxsToCombine;
  if (sampleIsPresent("dyee_*")&&sampleIsPresent("dymm_*")){
    v_prfxsToCombine.clear();
    v_prfxsToCombine.push_back("dyee");
    v_prfxsToCombine.push_back("dymm");
    hist::combineHists(v_prfxsToCombine, "DYeemm");
  }

  if (sampleIsPresent("wz_*")&&sampleIsPresent("zz_*")){
    v_prfxsToCombine.clear();
    v_prfxsToCombine.push_back("wz");
    v_prfxsToCombine.push_back("zz");
    hist::combineHists(v_prfxsToCombine, "VZ");
  }

  vector<TString> v_samples;    
  vector<Color_t> v_colors;
  vector<TString> v_legEntries;
  vector<Style_t> v_styles;
  
  //the order you put the sample into the vector is the order
  //that the histogram gets drawn. So if you put ttdil in first,
  //it will be drawn first, i.e, it will be drawn on the bottom
  //of the stack!
  
  if (sampleIsPresent("ww_*")){
    v_samples.push_back("ww");
    // v_colors.push_back(kRed);
    v_colors.push_back(kYellow+2);
    v_legEntries.push_back("#font[12]{WW} signal");
    v_styles.push_back(1001);
  }

  if (sampleIsPresent("wjets_*")){
    v_samples.push_back("wjets");
    // v_colors.push_back(kGreen);
    v_colors.push_back(kCyan);
    v_legEntries.push_back("#font[12]{W}#rightarrowl#nu");
    v_styles.push_back(1001);    
  }

  if (sampleIsPresent("ttbar_*")){
    v_samples.push_back("ttbar");
    v_colors.push_back(kMagenta);
    v_legEntries.push_back("#font[12]{t#bar{t}}");
    v_styles.push_back(1001);
  }
    
  if (sampleIsPresent("tw_*")){
    v_samples.push_back("tw");
    v_colors.push_back(kMagenta+1);
    v_legEntries.push_back("single top");
    v_styles.push_back(1001);
  }

  if (sampleIsPresent("dytt_*")){
    v_samples.push_back("dytt");
    v_colors.push_back(kAzure+8);
    v_legEntries.push_back("Z/#gamma*#rightarrow#tau^{+}#tau^{-}");
    v_styles.push_back(1001);
  }

  if (sampleIsPresent("DYeemm_*")){
    v_samples.push_back("DYeemm");
    // v_colors.push_back(kAzure - 2);
    v_colors.push_back(kBlue);
    v_legEntries.push_back("Z/#gamma*#rightarrowl^{+}l^{-}");
    v_styles.push_back(1001);
  }

  if (sampleIsPresent("VZ_*")){
    v_samples.push_back("VZ");
    // v_colors.push_back(kWhite);
    v_colors.push_back(kGreen);
    v_legEntries.push_back("(W/Z)Z#rightarrowl^{+}l^{-}");
    v_styles.push_back(1001);
  }

  //data should be the last thing put into the vectors
  if (sampleIsPresent("data_*")){
    v_samples.push_back("data");
    v_colors.push_back(kBlack);
    v_legEntries.push_back("Data");
    v_styles.push_back(20);
  }
  /*
    void browseStacks(vector<TString> v_samples, vector<Color_t> v_colors, 
    TString outfile, vector<TString> v_legEntries, bool drawLogY = false, 
    vector<Style_t> v_style = vector<Style_t>(), bool drawFullErrors) {
  */
  browseStacks(v_samples, v_colors, dataFName, v_legEntries, drawLogY, v_styles, drawFullErrors, scaleMC);
  
}
