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

  std::string cms2Location = "";
  if (getenv("CMS2_LOCATION") ) cms2Location = Form("%s/NtupleMacros", getenv("CMS2_LOCATION"));
  else {
    std::cout<<"CMS2_LOCATION is not set: it should point to what you where a nominal /pathOrRelativePath/CMS2 is"<<std::endl;
    if (getenv("CMSSW_BASE")){
      cms2Location = Form("%s/src/CMS2/NtupleMacros", getenv("CMSSW_BASE"));
      std::cout<<"You have CMSSW_BASE -- will try "<<cms2Location.c_str()<<std::endl;
    }
  }

  //punch in more quesses where NtupleMacros are -- don't keep multiple copies ;)
  gSystem->AddIncludePath(Form(" -w -I./ -I%s/CORE -I../CMS2/NtupleMacros/CORE -I../../CMS2/NtupleMacros/CORE\
  -I%s -I../CMS2/NtupleMacros -I../../CMS2/NtupleMacros ", cms2Location.c_str(),cms2Location.c_str() ));
}
