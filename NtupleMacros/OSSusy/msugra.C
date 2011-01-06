{

 //---------------
 // the limit
 //---------------
 float nev = 4.7;

 //-------------------
 // a kfactor fudged
 //-------------------
 float kfact = 1.4;

 //----------------------------
 // Benny and/or Sanjay messed up the normalization by a factor of 1000
 //----------------------------
 float fudge=1000.;

 //---------------------------
 // Load some tools
 //---------------------------
 //gROOT->LoadMacro("histio.cc");
 //gStyle->SetOptStat(0);

 //-----------------------------------
 // Here we load the yield histogram
 //-----------------------------------
 //loadHist("ossusy_JPT_tcmet_bitmask_LMscan.root", 0, "lmgridyield");
 TFile *f = TFile::Open("ossusy_JPT_tcmet_bitmask_LMscan.root");
 TH2F* yield = (TH2F*) f->Get("lmgridyield");

 yield->RebinX(2);
 yield->RebinY(2);
 yield->Scale(0.25);

 //--------------------------------------------
 // Get the parameters of the yield histogram
 //--------------------------------------------
 float xmin = yield->GetXaxis()->GetXmin();
 float xmax = yield->GetXaxis()->GetXmax();
 float ymin = yield->GetYaxis()->GetXmin();
 float ymax = yield->GetYaxis()->GetXmax();
 int nx     = yield->GetXaxis()->GetNbins();
 int ny     = yield->GetYaxis()->GetNbins();

 cout << "xmin = " << xmin <<endl;
 cout << "xmax = " << xmax <<endl;
 cout << "ymin = " << ymin <<endl;
 cout << "ymax = " << ymax <<endl;
 cout << "nx  = " << nx << endl;
 cout << "ny  = " << ny << endl;

 //---------------------------------------------------------------
 // We book a 1D histogram to keep the results in
 //---------------------------------------------------------------
 TH1F* limit = new TH1F("limit","limit",nx,xmin,xmax);
 //---------------------------------------------------------------
 // Use Sanjay's method, scan from the top and hope for the best
 //---------------------------------------------------------------
 float ybinsize = (ymax-ymin)/ny;
 float xbinsize = (xmax-xmin)/nx;
 for (int ix=1; ix<=nx; ix++) {
   float x = xmin + ix*xbinsize - 0.5*xbinsize;
   bool foundOne = false;

   for (int iy=ny; iy>0; iy--) {
     float this_ = yield->GetBinContent(ix,iy);
     this_ = kfact*fudge*this_;
     if (this_ > nev) {
       float yupperedge = ymin + iy*ybinsize;
       limit->Fill(x,yupperedge);
       // cout << yupperedge << " " << iy << endl;
       foundOne=true;
       break;
     }
   } //close iy loop
   //    if (!foundOne) limit->Fill(x,ymin);

 }   //close ix loop


 //----------------------------
 // Plot it
 //----------------------------
 limit->SetMinimum(90.);
 limit->SetTitle("Exclusion Curve in mSUGRA Space");
 //limit->SetTitle("Exclusion Curve on the Camel's Ass");
 limit->GetXaxis()->SetTitle("m_{0} (GeV)");
 limit->GetYaxis()->SetTitle("m_{1/2} (GeV)");
 limit->SetLineWidth(2);
 limit->SetLineColor(2);
 limit->GetXaxis()->SetRangeUser(0,500);
 limit->GetYaxis()->SetRangeUser(100,400);
 limit->Draw("c");
 //limit->Draw();

}
