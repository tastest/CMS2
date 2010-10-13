{
TH1::AddDirectory(true);
gSystem->Load("libPhysics.so");  
gSystem->Load("libEG.so");

gROOT->LoadMacro("doAnalysis.C+");
gROOT->ProcessLine(".x processData.C");

}
