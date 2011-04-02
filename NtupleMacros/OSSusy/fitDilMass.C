

void fitDilMass(){
  
  using namespace RooFit;
  gROOT->SetStyle("Plain");

  //declare variables
  RooRealVar mll("mll",    "dilepton mass",              0,100 ,"GeV");
  RooRealVar edge("edge",  "dilepton mass edge",         45,55 ,"GeV");
  RooRealVar mean("mean",  "mean of mll resolution",     0,100 ,"GeV");
  RooRealVar sigma("sigma",  "sigma of mll resolution",  0,10  ,"GeV");

  //set variables
  mll.setBins(100);
  edge.setVal(50.);
  mean.setVal(0.);
  sigma.setVal(2.);
  mean.setConstant();
  sigma.setConstant();

  //get data
  TFile* file = TFile::Open("data.root");
  TH1D*  hist = static_cast<TH1D *>((file)->Get("dilmass"));
  RooDataHist data("data","data histogram",RooArgSet(mll),hist);

  //declare and make PDFs
  RooGenericPdf tri("tri","( (mll < edge ) * mll )", RooArgList(mll , edge) );
  RooGaussian gaus("gaus","Detector Gaussian",mll,mean,sigma); 
  RooNumConvPdf trigaus("trigaus","triangle + gaussian function",mll,tri,gaus); 
  trigaus.setConvolutionWindow(mean,sigma,5);

  //do fit
  RooFitResult *result = trigaus.fitTo( data , Save() , Minos(kFALSE) );

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
  frame->SetMinimum(-10);
  frame->SetMaximum(100);
}
