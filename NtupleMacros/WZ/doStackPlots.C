{

  gROOT->SetBatch(kTRUE);

  gROOT->LoadMacro("../Tools/histtools.C");

  bool goOn = true;

  goOn = makeTrilepStackPlots("hNjetsBothLeptonsVeto");
  if ( !goOn ) return;

  gROOT->SetBatch(kFALSE);

}
