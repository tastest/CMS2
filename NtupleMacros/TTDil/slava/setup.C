// Loads a bunch of tools in interpretative mode
void setup(bool skipFWLite = false){
  if (!skipFWLite){
    gSystem->Load("libFWCoreFWLite");
    AutoLibraryLoader::enable();
  }

  gSystem->Load("libGui.so");
  gSystem->Load("libPhysics.so");
  
  gROOT->LoadMacro("getMyHistosNames.C+");
  gROOT->LoadMacro("histtools.C+");
  gROOT->LoadMacro("browseStacks.C+");
  
}
