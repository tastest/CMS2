// Loads a bunch of tools in interpretative mode
void setup(bool noCompile=false){
  gSystem->Load("libMiniFWLite.so");

  gSystem->Load("libGui.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMathMore.so");
  
  // Load and compile something to allow proper treatment of vectors
  // Not clear that it is needed
  if (noCompile){
    gSystem->Load("libloader.so");
    
    gSystem->Load("libgetMyHistosNames.so");
    gSystem->Load("libhisttools.so");
    gSystem->Load("libbrowseStacks.so");
    gSystem->Load("libCMS2.so");
    gSystem->Load("libMT2.so");
    gSystem->Load("libttDilCounts_looper.so");
    gSystem->Load("libttDilRefS10_looper.so");
  } else {
    gSystem->CompileMacro("loader.C", "++k", "libloader");
    
    gSystem->CompileMacro("histscripts/getMyHistosNames.C", "++k", "libgetMyHistosNames");
    gSystem->CompileMacro("histscripts/histtools.C", "++k", "libhisttools");
    gSystem->CompileMacro("histscripts/browseStacks.C", "++k", "libbrowseStacks");
    gSystem->CompileMacro("CORE/CMS2.cc", "++k", "libCMS2");
    gSystem->CompileMacro("CORE/MT2.cc", "++k", "libMT2");
    gSystem->CompileMacro("ttDilCounts_looper.C", "++k", "libttDilCounts_looper");
    gSystem->CompileMacro("ttDilRefS10_looper.C", "++k", "libttDilRefS10_looper");
  } 
    gROOT->ProcessLine(".L tdrstyle.C");
    setTDRStyle();

  //punch in more guesses where NtupleMacros are -- don't keep multiple copies ;)
  gSystem->AddIncludePath(Form(" -w -I./ " ));
}
