{
gROOT->ProcessLine(".x init.C");
gROOT->ProcessLine(".L doAnalysis.C+");
gROOT->ProcessLine(".x processData.C");
// gROOT->ProcessLine(".x showResults.C+");
}
