// Loads a bunch of tools in interpretative mode
void setup(bool skipFWLite = false){
  if (!skipFWLite){
    gSystem->Load("libFWCoreFWLite");
    AutoLibraryLoader::enable();
  }

  gSystem->Load("libGui.so");
  gSystem->Load("libPhysics.so");
  
  gSystem->CompileMacro("getMyHistosNames.C", "++k", "libgetMyHistosNames");
  gSystem->CompileMacro("histtools.C", "++k", "libhisttools");
  gSystem->CompileMacro("browseStacks.C", "++k", "libbrowseStacks");
  
}
