void macro_getYields(const char* rootfile){
gROOT->ProcessLine(".L getYields.C+");
getYields(Form("results/%s",rootfile),1.0,false,"%6.2f",false);
getYields(Form("results/%s",rootfile),1.0,false,"%6.2f",true);
}
