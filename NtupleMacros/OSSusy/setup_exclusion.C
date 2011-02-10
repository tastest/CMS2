{
  gROOT->ProcessLine(".L tdrstyle_SUSY.C");
  //gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();
  gROOT->ProcessLine(".L ExclusionPlot.C+");
}

