{

  gROOT->SetBatch(kTRUE);

  gROOT->LoadMacro("../Tools/histtools.C");

  bool goOn = true;

  goOn = makeTrilepStackPlots("hNjets");
  if ( !goOn ) return;
  goOn = makeTrilepStackPlots("hPtFirst");
  if ( !goOn ) return;
  goOn = makeTrilepStackPlots("hPtSecond");
  if ( !goOn ) return;
  goOn = makeTrilepStackPlots("hPtThird");
  if ( !goOn ) return;
  goOn = makeTrilepStackPlots("hMET");

  gROOT->SetBatch(kFALSE);

}
