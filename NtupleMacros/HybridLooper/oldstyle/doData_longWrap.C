{
//
// the looper
//
gSystem->Load("libTree.so");
gSystem->Load("libPhysics.so");
gSystem->Load("libEG.so");
gSystem->Load("libMathCore.so");
gSystem->Load("libCMS2NtupleMacrosCORE.so");
gSystem->Load("libCMS2NtupleMacrosLooper.so");

gROOT->ProcessLine(".x doData_long.C+");



}
