{
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L xsecLoop.C+");
  gROOT->ProcessLine(".x plotEl.C+");
}
