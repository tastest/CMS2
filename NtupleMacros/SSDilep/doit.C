{
gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
gSystem->SetAclicMode(TSystem::kDebug);
gROOT->ProcessLine(".x processData.C");
}
