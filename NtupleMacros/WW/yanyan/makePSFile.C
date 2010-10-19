#include <iostream>
#include <vector>
#include "TStyle.h"
#include "histtools.C"
#include "browseStacks.C"
#include "getMyHistosNames.C"
#include "TSystem.h"
#include "tdrStyle.C"
#include "CommonFunctions.C"

/*
  drawFullErrors = draw the errors in all their glory. Includes the systematic errors:
                   50% on the WJets estimate
		   100% on the QCD estimate
		   50% on DY estimate
		   50% on the backgrounds from MC
*/

void makePSFile(const TString dataFName, const TString mcFName, 
		float scaleMC,  bool drawLogY = true, bool drawFullErrors = false) {
  
  setTDRStyle();

  hist::loadHist(mcFName.Data(), 0, "*_hmetVal_*");
  hist::loadHist(mcFName.Data(), 0, "*_hmetProjVal_*");
  hist::loadHist(mcFName.Data(), 0, "*_hmaxPFJetPtVal_*");
  hist::loadHist(mcFName.Data(), 0, "*_hdilMassVal_*");
  hist::scale("*_*", scaleMC);     

  vector<TString> v_prfxsToCombine;
  v_prfxsToCombine.clear();
  v_prfxsToCombine.push_back("dyee");
  v_prfxsToCombine.push_back("dymm");
  hist::combineHists(v_prfxsToCombine, "DYeemm");


  v_prfxsToCombine.clear();
  v_prfxsToCombine.push_back("wz");
  v_prfxsToCombine.push_back("zz");
  hist::combineHists(v_prfxsToCombine, "VZ");

  hist::loadHist(dataFName.Data(), 0, "data_hmetVal_*");
  hist::loadHist(dataFName.Data(), 0, "*_hmetProjVal_*");
  hist::loadHist(dataFName.Data(), 0, "data_hmaxPFJetPtVal_*");
  hist::loadHist(dataFName.Data(), 0, "*_hdilMassVal_*");

  vector<TString> v_samples;    
  vector<Color_t> v_colors;
  vector<TString> v_legEntries;
  vector<Style_t> v_styles;
  
  //the order you put the sample into the vector is the order
  //that the histogram gets drawn. So if you put ttdil in first,
  //it will be drawn first, i.e, it will be drawn on the bottom
  //of the stack!
  
  v_samples.push_back("ww");
  v_colors.push_back(kRed);
  v_legEntries.push_back("#font[12]{WW} signal");
  v_styles.push_back(1001);

  v_samples.push_back("wjets");
  v_colors.push_back(kGreen);
  v_legEntries.push_back("#font[12]{W}#rightarrowl#nu");
  v_styles.push_back(1001);    
  
  v_samples.push_back("ttbar");
  v_colors.push_back(kMagenta);
  v_legEntries.push_back("#font[12]{t#bar{t}}");
  v_styles.push_back(1001);
  
  v_samples.push_back("tw");
  v_colors.push_back(kMagenta+1);
  v_legEntries.push_back("single top");
  v_styles.push_back(1001);
  
  v_samples.push_back("dytt");
  v_colors.push_back(kAzure+8);
  v_legEntries.push_back("Z/#gamma*#rightarrow#tau^{+}#tau^{-}");
  v_styles.push_back(1001);
  
  v_samples.push_back("DYeemm");
  v_colors.push_back(kAzure - 2);
  v_legEntries.push_back("Z/#gamma*#rightarrowl^{+}l^{-}");
  v_styles.push_back(1001);
    
  v_samples.push_back("VZ");
  v_colors.push_back(kWhite);
  v_legEntries.push_back("(W/Z)Z#rightarrowl^{+}l^{-}");
  v_styles.push_back(1001);

  //data should be the last thing put into the vectors
  v_samples.push_back("data");
  v_colors.push_back(kBlack);
  v_legEntries.push_back("Data");
  v_styles.push_back(20);

  /*
    void browseStacks(vector<TString> v_samples, vector<Color_t> v_colors, 
    TString outfile, vector<TString> v_legEntries, bool drawLogY = false, 
    vector<Style_t> v_style = vector<Style_t>(), bool drawFullErrors) {
  */
  browseStacks(v_samples, v_colors, dataFName, v_legEntries, drawLogY, v_styles, drawFullErrors, scaleMC);
  
}
