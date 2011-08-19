{
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L xsecLoop.C+");
  gROOT->ProcessLine(".x plotMu.C+");
}
