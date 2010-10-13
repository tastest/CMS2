// Loads a bunch of tools in interpretative mode
{
  gSystem->Load("libGui.so");
  gSystem->Load("libPhysics.so");
  
  gROOT->LoadMacro("getMyHistosNames.C+");
  gROOT->LoadMacro("histtools.C+");
  gROOT->LoadMacro("browseStacks.C+");
  
}
