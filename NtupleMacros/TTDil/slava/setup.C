// Loads a bunch of tools in interpretative mode
void setup(){
  gSystem->Load("libMiniFWLite.so");

  gSystem->Load("libGui.so");
  gSystem->Load("libPhysics.so");
  
  // Load and compile something to allow proper treatment of vectors
  // Not clear that it is needed
  gSystem->CompileMacro("loader.C", "++k", "libloader");

  gSystem->CompileMacro("histscripts/getMyHistosNames.C", "++k", "libgetMyHistosNames");
  gSystem->CompileMacro("histscripts/histtools.C", "++k", "libhisttools");
  gSystem->CompileMacro("histscripts/browseStacks.C", "++k", "libbrowseStacks");
  gSystem->CompileMacro("CORE/CMS2.cc", "++k", "libCMS2");
  gSystem->CompileMacro("ttDilCounts_looper.C", "++k", "libttDilCounts_looper");
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();

  //punch in more guesses where NtupleMacros are -- don't keep multiple copies ;)
  gSystem->AddIncludePath(Form(" -w -I./ " ));
}
