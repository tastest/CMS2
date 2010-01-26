// Loads a bunch of tools in interpretative mode
void setup(bool skipFWLite = false){
  if (!skipFWLite){
    gSystem->Load("libFWCoreFWLite");
    AutoLibraryLoader::enable();
  }

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

  /*
  std::string cms2Location = "";
  if (getenv("CMS2_LOCATION") ) cms2Location = Form("%s/NtupleMacros", getenv("CMS2_LOCATION"));
  else if (1>2){
    std::cout<<"CMS2_LOCATION is not set: it should point to what you where a nominal /pathOrRelativePath/CMS2 is"<<std::endl;
    if (getenv("CMSSW_BASE")){
      cms2Location = Form("%s/src/CMS2/NtupleMacros", getenv("CMSSW_BASE"));
      std::cout<<"You have CMSSW_BASE -- will try "<<cms2Location.c_str()<<std::endl;
    }
  }
  */

  //punch in more guesses where NtupleMacros are -- don't keep multiple copies ;)
  gSystem->AddIncludePath(Form(" -w -I./ " ));
}
