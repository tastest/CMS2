void makeAllPlots(char* fname, bool logScale=false){
  gROOT->SetStyle("Plain");
  gROOT->LoadMacro("loader.C+");
  gROOT->ProcessLine(".x setup.C");
  hist::loadHist(fname);
  
  browseStacks(true, false, true, 1.1, logScale, logScale ? false : true);
  gSystem->Exit(0);
}
