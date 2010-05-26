void makeAllPlots(char* fname, bool logScale=false,char* hPatt = 0, char* refName=0, 
		  int titleStyle = 1, int colorStyle=0, int orderStyle=0,
		  bool noLegend = false, char* bsmName=0){
  gROOT->SetStyle("Plain");
  ROOT->gErrorIgnoreLevel = 0;
  gROOT->ProcessLine(".x setup.C(1)"); //don't need FWLite to make plots
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
	       colorStyle,//mod              // color scheme (0 original tas, 1 as in pas 09-002, 2 almost like 1 from the wheel
	       noLegend,                     // no legend on the plot
	       orderStyle,//mod              // order scheme
	       bsmName);
  gSystem->Exit(0);
}
