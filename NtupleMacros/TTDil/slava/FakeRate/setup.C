// Loads a bunch of tools in interpretative mode
void setup(){
  
  gSystem->Load("libGui.so");
  gSystem->Load("libPhysics.so");
  
 
  // Load and compile something to allow proper treatment of vectors
  // Not clear that it is needed
  gSystem->CompileMacro("../loader.C", "++k", "libloader");
  
  gSystem->CompileMacro("../histscripts/getMyHistosNames.C", "++k", "libgetMyHistosNames");
  gSystem->CompileMacro("../histscripts/histtools.C", "++k", "libhisttools");
  gSystem->CompileMacro("../histscripts/browseStacks.C", "++k", "libbrowseStacks");
  gSystem->CompileMacro("CORE/CMS2.cc", "++k", "libCMS2");

}
