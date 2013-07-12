void macro_drawHists(const char* rootfile){
gROOT->ProcessLine(".L drawHists.C+");
makePSFile(Form("results/%s",rootfile),5,0,"default",4.98,0,1,1);
}
