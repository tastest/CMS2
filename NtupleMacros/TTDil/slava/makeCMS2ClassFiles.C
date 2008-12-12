// Original author: Puneeth Kalavase (UCSB)
// 
/*Root macro to make a header and .C file for basic analysis
  of CMS2 ntuples. Usage:

  [kalavase@stau ~/rootmacros]$ root
  root [0] .L makeCMS2ClassFiles.C
  //second string is optional. The classname is CMS2 by default
  root [1] makeCMS2Header("tablemaker_Zmumu_ntuple.root","classname")

*/

makeCMS2Files(std::string fname, std::string className = "") {

  using namespace std;
  
  TFile *f = TFile::Open( fname.c_str() );
  if(f->IsZombie()) { 
    cout << "File is not a valid root file, or root file is corruped" << endl;
    cout << "Exiting..." << endl;
  }
  
  //class is CMS2 by default
  std::string Classname = className=="" ? "CMS2" : className;
  
  ofstream headerf;
  ofstream headerBLOBf;
  ofstream headcodef;
  ofstream codef;

  headerf.open((Classname+".h").c_str());
  headerBLOBf.open((Classname+"_BLOB.h").c_str());
  headcodef.open((Classname+"_looper.h").c_str());
  codef.open((Classname+"_looper.C").c_str());
  headerf << "// -*- C++ -*-" << endl;
  headerf << "#ifndef " << Classname << "_H" << endl;
  headerf << "#define " << Classname << "_H" << endl;
  headerf << "#include \"Math/LorentzVector.h\"" << endl;
  headerf << "#include \"Math/Point3D.h\"" << endl;
  headerf << "#include \"TMath.h\"" << endl;
  headerf << "#include \"TBranch.h\"" << endl;
  headerf << "#include \"TTree.h\"" << endl << endl;
  headerf << "#include <vector> " << endl;
  headerf << "using namespace std; " << endl;
  //  headerf << "class " << Classname << " { " << endl;
  //  headerf << "private: " << endl;
  //  headerf << "protected: " << endl;
  headerf << "\tunsigned int index_;" << endl;
  TTree *ev = (TTree*)f->Get("Events");
  TList *fullarray =  ev->GetListOfAliases();
  TList *aliasarray = new TList();
    
  //for(Int_t i = 0; i < aliasarray->GetSize(); ++i) {
  for(Int_t i = 0; i < fullarray->GetSize(); ++i) {
    TString aliasname(fullarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    TString branchname(branch->GetName());
    if(!branchname.BeginsWith("int") && 
       !branchname.BeginsWith("uint") && 
       !branchname.BeginsWith("float") &&
       !branchname.BeginsWith("double") ) continue;
    aliasarray->Add(fullarray->At(i));
  }
  
  headerf << "#include \""<< (Classname+"_BLOB.h\"")<< endl;
  headerf << "// Included in the BLOB void Init(TTree *tree) {" << endl;
  headerf << "// Included in the BLOB void GetEntry(unsigned int idx) {" << endl;
  
  for(Int_t i = 0; i< aliasarray->GetSize(); ++i) {
    
    //Class name is blank for a int of float
    TString aliasname(aliasarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    
    TString classname = branch->GetClassName();
    if ( classname.Contains("vector") ) {
      classname = classname(0,classname.Length()-2);
      classname.ReplaceAll("edm::Wrapper<","");
      headerBLOBf << "\t" << classname << "\t" << aliasname << "_;" << endl;
    } else {
      classname = classname(0,classname.Length()-1);
      classname.ReplaceAll("edm::Wrapper<","");
      headerBLOBf << "\t" << classname << "\t" << aliasname << "_;" << endl;
    }
    headerBLOBf << "\tTBranch *" << Form("%s_branch",aliasname.Data()) << ";" << endl;
    headerBLOBf << "\tbool " << Form("%s_isLoaded",aliasname.Data()) << ";" << endl;
  }
  
  
  headerBLOBf << "void Init(TTree *tree) {" << endl;

  // SetBranchAddresses for LorentzVectors
  for(Int_t i = 0; i< aliasarray->GetSize(); i++) {
    TString aliasname(aliasarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    TString classname = branch->GetClassName();
    if ( !classname.Contains("vector<vector") ) {
      if ( classname.Contains("Lorentz") ) {
	   headerBLOBf << "\t" << Form("%s_branch",aliasname.Data()) << " = 0;" << endl;
	   headerBLOBf << "\t" << "if (tree->GetAlias(\"" << aliasname << "\") != 0) {" << endl;
	   headerBLOBf << "\t\t" << Form("%s_branch",aliasname.Data()) << " = tree->GetBranch(tree->GetAlias(\"" << aliasname << "\"));" << endl;
	   headerBLOBf << "\t\t" << Form("%s_branch",aliasname.Data()) << "->SetAddress(&" << aliasname << "_);" << endl << "\t}" << endl;
	   headerBLOBf << "\t" << "if(" << Form("%s_branch",aliasname.Data()) << " == 0 ) {" << endl;
	   headerBLOBf << "\t" << "cout << \"Branch " << aliasname.Data() << " does not exist.\" << endl;" << endl;
	   headerBLOBf << "\t" << "}" << endl;
      }
    }
  }


  // SetBranchAddresses for everything else
  headerBLOBf << "  tree->SetMakeClass(1);" << endl;
  for(Int_t i = 0; i< aliasarray->GetSize(); i++) {
    TString aliasname(aliasarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    TString classname = branch->GetClassName();
    if ( !classname.Contains("Lorentz") || classname.Contains("vector<vector") ) {
	 headerBLOBf << "\t" << Form("%s_branch",aliasname.Data()) << " = 0;" << endl;
	 headerBLOBf << "\t" << "if (tree->GetAlias(\"" << aliasname << "\") != 0) {" << endl;
	 headerBLOBf << "\t\t" << Form("%s_branch",aliasname.Data()) << " = tree->GetBranch(tree->GetAlias(\"" << aliasname << "\"));" << endl;
	 headerBLOBf << "\t\t" << Form("%s_branch",aliasname.Data()) << "->SetAddress(&" << aliasname << "_);" << endl << "\t}" << endl;
	 headerBLOBf << "\t" << "if(" << Form("%s_branch",aliasname.Data()) << " == 0 ) {" << endl;
	 headerBLOBf << "\t" << "cout << \"Branch " << aliasname.Data() << " does not exist.\" << endl;" << endl;
	 headerBLOBf << "\t" << "}" << endl;
    }
  }

  headerBLOBf << "  tree->SetMakeClass(0);" << endl;
  headerBLOBf << "}" << endl;

  // GetEntry
  headerBLOBf << "void GetEntry(unsigned int idx) " << endl;
  headerBLOBf << "\t// this only marks branches as not loaded, saving a lot of time" << endl << "\t{" << endl;
  headerBLOBf << "\t\tindex_ = idx;" << endl;
  for(Int_t i = 0; i< aliasarray->GetSize(); i++) {
       TString aliasname(aliasarray->At(i)->GetName());
       headerBLOBf << "\t\t" << Form("%s_isLoaded",aliasname.Data()) << " = false;" << endl;
  }
  headerBLOBf << "\t}" << endl << endl;

  // accessor functions
  for (Int_t i = 0; i< aliasarray->GetSize(); i++) {
       TString aliasname(aliasarray->At(i)->GetName());
       TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
       TString classname = branch->GetClassName();
       if ( classname.Contains("vector") ) {
	    classname = classname(0,classname.Length()-2);
	    classname.ReplaceAll("edm::Wrapper<","");
	    headerBLOBf << "\t" << classname << " &" << aliasname << "()" << endl;
       } else {
	    classname = classname(0,classname.Length()-1);
	    classname.ReplaceAll("edm::Wrapper<","");
	    headerBLOBf << "\t" << classname << " &" << aliasname << "()" << endl;
       }
       TString aliasname(aliasarray->At(i)->GetName());
       headerBLOBf << "\t{" << endl;
       headerBLOBf << "\t\t" << "if (not " << Form("%s_isLoaded) {",aliasname.Data()) << endl;
       headerBLOBf << "\t\t\t" << "if (" << Form("%s_branch",aliasname.Data()) << " != 0) " << endl;
       headerBLOBf << "\t\t\t\t" << Form("%s_branch",aliasname.Data()) << "->GetEntry(index_);" << endl;
       headerBLOBf << "\t\t\t" << Form("%s_isLoaded",aliasname.Data()) << " = true;" << endl;
       headerBLOBf << "\t\t" << "}" << endl;
       headerBLOBf << "\t\t" << "return " << aliasname << "_;" << endl << "\t}" << endl;
  }
  //  headerf << "};" << endl << endl;

  headerf << "#endif" << endl;

  //  headerf << "public: " << endl;
  headcodef << "// -*- C++ -*-" << endl;
  headcodef << "#ifndef " << Classname << "_looper_H" << endl;
  headcodef << "#define " << Classname << "_looper_H" << endl;
  headcodef << "#include \""<< (Classname+".h\"")<< endl;
  headcodef << "#include \"TH1F.h\""  << endl << endl;
  headcodef << "class "<<Classname<<"_looper {"  << endl << endl;
  headcodef << "\t public: " << endl;
  headcodef << "int ScanChain( TChain* chain, int nEvents=-1);" << endl;
  headcodef << "\t TH1F *samplehisto;" << endl;
  headcodef << "};" << endl;
  headcodef << "#endif" << endl;

  codef << "/* Usage:" << endl;
  codef << "   root[0] .L " << Classname << "_looper.C++" << endl;
  codef << "   root [1] TFile *_file0 = TFile::Open(\"ntuple_file.root\")" << endl;
  codef << "   root [2] TChain *chain = new TChain(\"Events\")" << endl;
  codef << "   root [3] chain->Add(\"ntuple_file.root\")" << endl;
  codef << "   root [4] " << Classname << "_looper a " << endl; 
  codef << "   root [5] a.ScanChain(chain) // will give the same results" << endl;
  codef << "*/" << endl;
  codef << "#include <iostream>" << endl;
  codef << "#include <vector>" << endl;
  codef << "" << endl;
  codef << "#include \"TChain.h\"" << endl;
  codef << "#include \"TFile.h\"" << endl;
  codef << "#include \"TDirectory.h\"" << endl;
  codef << "#include \"TROOT.h\"" << endl;
  codef << "" << endl;
  codef << "#include \"" + Classname+"_looper.h\"" << endl;
  codef << "" << endl;
  codef << "" << endl;
  
  codef << "int "<<Classname<<"_looper::ScanChain( TChain* chain, int nEvents) {" << endl;
  codef << "" << endl;
  codef << "  TObjArray *listOfFiles = chain->GetListOfFiles();" << endl;
  codef << "" << endl;
  codef << "  unsigned int nEventsChain=0;" << endl;
  codef << "  if(nEvents==-1) " << endl << "     nEvents = chain->GetEntries();" << endl;
  codef << "  else nEventsChain = nEvents;" << endl;
  
  codef << "  unsigned int nEventsTotal = 0;" << endl;
  codef << "  TDirectory *rootdir = gDirectory->GetDirectory(\"Rint:\");" << endl << endl;
  codef << "  TH1F *samplehisto = new TH1F(\"samplehisto\", \"Example histogram\", 200,0,200);" << endl;
  codef << "  samplehisto->SetDirectory(rootdir);" << endl;

  codef << "  // file loop" << endl;
  codef << "  TIter fileIter(listOfFiles);" << endl;
  codef << "  TFile *currentFile = 0;" << endl;
  codef << "  while ( currentFile = (TFile*)fileIter.Next() ) {" << endl;
  codef << "    TFile f(currentFile->GetTitle());" << endl;
  codef << "    TTree *tree = (TTree*)f.Get(\"Events\");" << endl;
  codef << "    Init(tree);" << endl;
  codef << "    " << endl;
  codef << "    //Event Loop" << endl;
  codef << "    unsigned int nEvents = tree->GetEntries();" << endl;
  codef << "    for( unsigned int event = 0; event < nEvents; ++event) {" << endl;
  codef << "      GetEntry(event);" << endl;
  codef << "      ++nEventsTotal;" << endl;
  codef << "      std::cout << \"els size: \" << els_p4().size() << \" \";" << endl;
  codef << "      std::cout << \"mus size: \" << mus_p4().size() << std::endl;" << endl;
  codef << "      for (unsigned int mus = 0; " << endl;
  codef << "           mus < mus_p4().size(); mus++) " << endl << endl;
  codef << "         samplehisto->Fill(mus_p4().at(mus).Pt());" << endl << endl;
  codef << "      for (unsigned int hyp = 0;" << endl;
  codef << "           hyp < hyp_jets_p4().size();" << endl;
  codef << "           ++hyp) {" << endl;
  codef << "        std::cout << \"hyp: \" << hyp << \"jet corrections:\";" << endl;
  codef << "        for ( unsigned int jet = 0;" << endl;
  codef << "              jet < hyp_jets_p4()[hyp].size();" << endl;
  codef << "              ++jet ) {" << endl;
  codef << "          std::cout << \" \" << hyp_jets_p4()[hyp][jet].pt();" << endl;
  codef << "        }" << endl;
  codef << "        std::cout << endl;" << endl;
  codef << "      }" << endl;
  codef << "      if ( hyp_jets_p4().size() == 0 ) {" << endl;
  codef << "        std::cout << \"no hypothesis!\" << std::endl;" << endl;
  codef << "      }" << endl;
  codef << "    }" << endl;
  codef << "  }" << endl;
  codef << "" << endl;
  codef << "  if ( nEventsChain != nEventsTotal ) {" << endl;
  codef << "    std::cout << \"ERROR: number of events from files is not equal to total number of events\" << std::endl;" << endl;
  codef << "  }" << endl;
  codef << "" << endl;
  codef << "  samplehisto->Draw();" << endl;
  codef << "  return 0;" << endl;
  codef << "}" << endl;
   
  headerf.close();
  headerBLOBf.close();
  headcodef.close();
  codef.close();
  f->Close();
}
