void makeAllPlots(char* fname, bool logScale=false,char* hPatt = 0, char* refName, int titleStyle = 1, bool noLegend = false, char* bsmName=0){
  gROOT->SetStyle("Plain");
  gROOT->ProcessLine(".x setup.C(true)"); //don't need FWLite to make plots
  hist::loadHist(fname, 0, hPatt);
  
  if (refName!=0){
    hist::loadHist(refName, "ref", hPatt);
  }  

  browseStacks(true,                         //makePictures
	       false,                        //wait
	       titleStyle,                   //addHistName
	       1.1,                          // maxYrescale
	       logScale,                     // logScale
	       logScale ? false : true,      // set min to 0
	       0,//mod                            // color scheme (0 original tas, 1 as in pas 09-002, 2 almost like 1 from the wheel
	       noLegend,                     // no legend on the plot
	       0,//mod                             // order scheme
	       bsmName);
  gSystem->Exit(0);
}
