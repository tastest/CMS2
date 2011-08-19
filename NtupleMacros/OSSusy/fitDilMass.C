

void fitDilMass(){
  
  using namespace RooFit;
  gROOT->SetStyle("Plain");

  //declare variables
  RooRealVar mll("mll",      "dilepton mass",              0,100 ,"GeV");
  RooRealVar edge("edge",    "dilepton mass edge",         0,100 ,"GeV");
  RooRealVar min("min",      "min dilepton mass",          0,100 ,"GeV");
  RooRealVar k("k",          "constant",                   0,100 ,"GeV");
  RooRealVar mean("mean",    "mean of mll resolution",     0,100 ,"GeV");
  RooRealVar sigma("sigma",  "sigma of mll resolution",    0,10  ,"GeV");

  //set variables
  mll.setBins(100);
  edge.setVal(50.);
  mean.setVal(0.);
  sigma.setVal(2.);
  min.setVal(10.);
  min.setConstant();
  k.setVal(1.);
  mean.setConstant();
  sigma.setConstant();
  mll.setMin(0);
  mll.setMax(100);

  //get data
  TFile* file = TFile::Open("LM0_dilmass.root");
  TH1D*  hist = static_cast<TH1D *>((file)->Get("dilmass"));
  RooDataHist data("data","data histogram",RooArgSet(mll),hist);

  //declare and make PDFs
  RooGenericPdf tri("tri","( (mll < min && mll > 0) * k + (mll > min && mll < edge ) * (mll - min) + (mll > edge && mll < 200) * k)", RooArgList(mll , edge, min, k) );
  RooGaussian gaus("gaus","Detector Gaussian",mll,mean,sigma); 
  RooNumConvPdf trigaus("trigaus","triangle + gaussian function",mll,tri,gaus); 
  trigaus.setConvolutionWindow(mean,sigma,25);

  //do fit
  RooFitResult *result = trigaus.fitTo( data , Save() , Minos(kFALSE) , SumW2Error(kFALSE));

  //make plot
  TCanvas *c1=new TCanvas("c1","",800,600);
  c1->cd();

  RooPlot* frame = mll.frame();
  frame->SetXTitle("dilepton mass (GeV)");
  frame->SetYTitle("");
  frame->SetTitle("");
  data.plotOn(frame);
  tri.plotOn(frame , LineStyle(kDashed));
  trigaus.plotOn( frame );
  trigaus.paramOn( frame );
  frame->Draw();
  //frame->SetMinimum(-10);
  //frame->SetMaximum(100);

 
  int nfloatpars = result->floatParsFinal().getSize();
  int ndf        = 20 - nfloatpars;
  float chi2     = frame->chiSquare(nfloatpars)*ndf;
  float prob     = TMath::Prob(chi2,ndf);
 
  cout << endl << endl;
  cout << "Fit Results-------------------------------------------"     << endl;
  cout << "chi2/ndf = " << chi2 << " / " << ndf                        << endl;
  cout << "prob     = " << prob                                        << endl;  
  cout << "edge     = " << edge.getVal() << " +/- " << edge.getError() << endl;
  cout << "k        = " << k.getVal()    << " +/- " << k.getError()    << endl;
  cout << "------------------------------------------------------"     << endl;
  cout << endl << endl;
}
