void doAll_oneFile(const char* fileName){


  gROOT->ProcessLine(".x setup.C(0)");
  gSystem->CompileMacro("ttDilCounts_looper.C", "++k", "libttDilCounts_looper");
  
  // 1957888 -- baseline
  // 35512320 -- baseline using tcmet
  // 1695744 -- baseline without MET
  // 1433600 -- baseline without zveto
  // 1926144 -- baseline without tight iso (only loose iso)
  // 1171456 -- baseline without MET, without zveto
  // 269606912 -- baseline without MET, without zveto, with MC truth match
  // 268435456 -- truthmatch only
  // 122880  -- baseline without MET, without zveto, without trigg reqm
  //  268566528  truthmatch && pt20Eta2.4
  // 269615104  truthmatch && pt20Eta2.4 && trigger
  // 269647872  truthmatch && pt20Eta2.4 && trigger && isolation

  unsigned int bitms[14] = {   1957888  ,  35512320 ,    1695744,       1433600,     1926144, 
			       1171456  ,          0,  269606912,     268435456,      122880,
			       268566528,  269615104,  269647872,     269647872};
  
  for (unsigned int iM=0; iM< 14; ++iM){
    TChain* ch = new TChain("Events");
    ch->Add(fileName);
    ttDilCounts_looper* looper = new ttDilCounts_looper();
    unsigned int bitmask =  bitms[iM];
    looper->ScanChain(ch,"ttdil", 1, 1, false, bitmask);
    hist::saveHist(Form("testTTHists_%d_%s.root", bitmask, looper->compactConfig.c_str()));
    hist::deleteHistos();
    cout << "Finished with bitmask "<<bitmask << endl;
    delete looper;
    delete ch;
  }
  gSystem->Exit(0);
}
