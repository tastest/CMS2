// Original author: Puneeth Kalavase (UCSB)
// 
/*Root macro to make a header and .C file for basic analysis
  of CMS2 ntuples. Usage:

  [kalavase@stau ~/rootmacros]$ root
  root [0] .L makeCMS2Header.C
  root [1] makeCMS2Header("tablemaker_Zmumu_ntuple.root")

*/

makeCMS2Files(std::string fname) {

  using namespace std;
  
  TFile *f = TFile::Open( fname.c_str() );
  if(f->IsZombie()) { 
    cout << "File is not a valid root file, or root file is corruped" << endl;
    cout << "Exiting..." << endl;
  }
  

  ofstream headerf;
  ofstream codef;
  headerf.open("CMS2.h");
  codef.open("CMS2.C");
  headerf << "#include \"Math/LorentzVector.h\"" << endl;
  headerf << "#include \"TMath.h\"" << endl;
  headerf << "#include \"TBranch.h\"" << endl;
  headerf << "#include \"TTree.h\"" << endl;


  
  TTree *ev = (TTree*)f->Get("Events");
  TObjArray *aliasarray = ev->GetListOfAliases();

  for(Int_t i = 0; i< aliasarray->GetSize(); ++i) {

    //Class name is blank for a int of float
    TString aliasname(aliasarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    TString classname = branch->GetClassName();
    if ( classname.Contains("vector") ) {
      classname = classname(0,classname.Length()-2);
      classname.ReplaceAll("edm::Wrapper<","");
      headerf << "\t" << classname << "\t" << aliasname << ";" << endl;
    } else {
      classname = classname(0,classname.Length()-1);
      classname.ReplaceAll("edm::Wrapper<","");
      headerf << "\t" << classname << "\t" << aliasname << ";" << endl;
    }
    headerf << "\tTBranch *" << Form("%s_branch",aliasname.Data()) << ";" << endl;
  }

  headerf << "void Init(TTree *tree) {" << endl;
  

  // SetBranchAddresses for LorentzVectors
  for(Int_t i = 0; i< aliasarray->GetSize(); i++) {
    TString aliasname(aliasarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    TString classname = branch->GetClassName();
    if ( classname.Contains("Lorentz") ) {
      headerf << "\t" << Form("%s_branch",aliasname.Data()) << " = tree->GetBranch(tree->GetAlias(\"" << aliasname << "\"));" << endl;
      headerf << "\t" << Form("%s_branch",aliasname.Data()) << "->SetAddress(&" << aliasname << ");" << endl;
    }
  }
  
  // SetBranchAddresses for everything else
  headerf << "  tree->SetMakeClass(1);" << endl;
  for(Int_t i = 0; i< aliasarray->GetSize(); i++) {
    TString aliasname(aliasarray->At(i)->GetName());
    TBranch *branch = ev->GetBranch(ev->GetAlias(aliasname.Data()));
    TString classname = branch->GetClassName();
    if ( !classname.Contains("Lorentz") ) {
      headerf << "\t" << Form("%s_branch",aliasname.Data()) << " = tree->GetBranch(tree->GetAlias(\"" << aliasname << "\"));" << endl;
      headerf << "\t" << Form("%s_branch",aliasname.Data()) << "->SetAddress(&" << aliasname << ");" << endl;
    }
  }
  headerf << "  tree->SetMakeClass(0);" << endl;

  headerf << "}" << endl;
  headerf << "void GetEntry(unsigned int index) {" << endl;
  

  //SetBranchAddresses
  for(Int_t i = 0; i< aliasarray->GetSize(); i++) {
    TString aliasname(aliasarray->At(i)->GetName());
    headerf << "\t" << Form("%s_branch",aliasname.Data()) << "->GetEntry(index);" << endl;
  }

  headerf << "}" << endl;



  codef << "//now make the source file" << endl;
  codef << "#include <iostream>" << endl;
  codef << "#include <vector>" << endl;
  codef << "" << endl;
  codef << "#include \"TChain.h\"" << endl;
  codef << "#include \"TFile.h\"" << endl;
  codef << "" << endl;
  codef << "#include \"TH1F.h\"" << endl;
  codef << "#include \"TH2F.h\"" << endl;
  codef << "" << endl;
  codef << "#include \"CMS2.h\"" << endl;
  codef << "" << endl;
  codef << "int ScanChain( TChain* chain) {" << endl;
  codef << "" << endl;
  codef << "  TObjArray *listOfFiles = chain->GetListOfFiles();" << endl;
  codef << "" << endl;
  codef << "  unsigned int nEventsChain = chain->GetEntries();" << endl;
  codef << "  unsigned int nEventsTotal = 0;" << endl;
  codef << "" << endl;
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
  codef << "	  std::cout << \"els size: \" << els_p4.size() << std::endl;" << endl;
  codef << "	  std::cout << \"mus size: \" << mus_p4.size() << std::endl;" << endl;
  codef << "    }" << endl;
  codef << "  }" << endl;
  codef << "" << endl;
  codef << "  if ( nEventsChain != nEventsTotal ) {" << endl;
  codef << "    std::cout << \"ERROR: number of events from files is not equal to total number of events\" << std::endl;" << endl;
  codef << "  }" << endl;
  codef << "" << endl;
  codef << "  return 0;" << endl;
  codef << "}" << endl;
  
   
  headerf.close();
  codef.close();
  f->Close();
}
