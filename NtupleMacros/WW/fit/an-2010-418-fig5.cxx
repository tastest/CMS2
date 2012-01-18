{
  gSystem->CompileMacro("RooATGCPdf.C","k");
  gSystem->CompileMacro("wwfit.C","k");
  setDefaults();
  setSigPdf_LZ_GZ();
  // setSigPdf_LZ_KG();
  // setBkgPdf();
  // sensitivity();
  // ww1DFits();
  wwATGC1DFits();
  // wwATGC1DLzKgFits();
  // fitData();
  // fitTop();
}
