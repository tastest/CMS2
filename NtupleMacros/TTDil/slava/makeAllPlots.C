void makeAllPlots(char* fname, bool logScale=false){
  gROOT->SetStyle("Plain");
  gSystem->CompileMacro("loader.C", "++k", "libloader");
  gROOT->ProcessLine(".x setup.C(true)"); //don't need FWLite to make plots
  hist::loadHist(fname);
  
  browseStacks(true, false, true, 1.1, logScale, logScale ? false : true);
  gSystem->Exit(0);
}
