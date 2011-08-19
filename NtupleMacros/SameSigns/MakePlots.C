void MakePlots(TString histofilename){
  gROOT->ProcessLine(".x ../Tools/setup.C");
  loadHist(histofilename);
  histofilename+="_out";
  browseStacks(true, false, histofilename);
}
