//
{
  
  //cout << "Starting" << endl;
  //cout << //(return of lines below are 0 if success)
  gROOT->ProcessLine(".x setup.C");
  gSystem->Load("histtools_C.so");

  gROOT->ProcessLine(".L doScanChain.C++");

  //first bool is data, second is MC
  doScanChain("output", true, false);
  //doScanChain("output_MC", false, true);  

  //both
  //doScanChain("output", true, true);

}
