{
   gROOT->SetStyle("Plain");
   gStyle->SetPalette(1);
   gStyle->SetOptStat("nemruoi");
   if ( gSystem->Getenv("CMSSW_VERSION") ){
     TString line;
     FILE *fp = gSystem->OpenPipe("scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/'", "r");
     line.Gets(fp);
     gSystem->ClosePipe(fp);
     gSystem->AddIncludePath(line.Data());
   }
   // cout << "loading..." <<endl;
   // gSystem->Load("libCintex");
   // Cintex::Enable();
   // gSystem->Load("libFWCoreFWLite");
   // AutoLibraryLoader::enable();
}
