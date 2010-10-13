#include "fitWjets.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooGlobalFunc.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "RooPolynomial.h"
#include <iostream>
#include "RooArgusBG.h"
#include "RooGenericPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooFitResult.h"
#include "TPaveText.h"
// best for asymmetrical distributions with long tails
// fake type:
//   0 - electron fakes
//   1 - muon fakes
// 
// pdf type:
//   0 - 1st degree polynomial (may lead to negative background estimates
//   1 - 0 degree polynomial (flat background)
//   2 - power law pdf
//   3 - non-parametric QCD based background pdf


TPaveText* estimate_background(int pdf_type, RooFitResult* result, float xmin, float xsig){
  RooRealVar* Nbkg   = dynamic_cast<RooRealVar*>(result->floatParsFinal().find("Nbkg"));
  RooRealVar* bkg_c  = dynamic_cast<RooRealVar*>(result->floatParsFinal().find("bkg_c"));
  RooRealVar* bkg_a1 = dynamic_cast<RooRealVar*>(result->floatParsFinal().find("bkg_a1"));
  char text[1024];
  text[0]=0;

  if ( pdf_type == 2 ){
    double c = (1-xsig)/(1-xmin);
    double central_value = pow(c,bkg_c->getVal()+1)*Nbkg->getVal();
    double a1 = Nbkg->getError()*pow(c,bkg_c->getVal()+1);
    double a2 = bkg_c->getError()*pow(c,bkg_c->getVal()+1)*log(c)*Nbkg->getVal();
    double variance = a1*a1 + a2*a2 + result->correlation("Nbkg","bkg_c")*a1*a2;
    sprintf(text, "%0.1f +/- %0.1f", central_value, sqrt(variance));
  }
  if ( pdf_type == 1 ){
    double c = (1-xsig)/(1-xmin);
    sprintf(text, "%0.1f +/- %0.1f", c*Nbkg->getVal(), c*Nbkg->getError() );
  }
  if ( pdf_type == 0 ){
    double c1 = 1-pow(xsig,2);
    double c2 = 1-pow(xmin,2);
    double cv = (1-xsig+bkg_a1->getVal()/2*c1)/(1-xmin+bkg_a1->getVal()/2*c2);
    double a1 = Nbkg->getError()*cv;
    double a2 = bkg_a1->getError()*Nbkg->getVal()*(
						   (1/2*c1/(1-xmin+bkg_a1->getVal()/2*c2)) - 
						   (1-xsig+bkg_a1->getVal()/2*c1)/pow((1-xmin+bkg_a1->getVal()/2*c2),2)*c2/2);
    double variance = a1*a1 + a2*a2 + result->correlation("Nbkg","bkg_a1")*a1*a2;
    sprintf(text, "%0.1f +/- %0.1f", cv*Nbkg->getVal(), sqrt(variance));
  }

  TPaveText* pText = new TPaveText(0.15,0.75,0.45,0.85,"NDC");
  if (text[0]!=0){
    std::cout << "\t" << text << std::endl;
    pText->AddText(text);
  }
  return pText;
}

RooFitResult* fit_isolation(RooAbsData* control_sample, 
			    RooAbsData* analysis_sample,
			    int fake_type,
			    int pdf_type,
			    const char* title,
			    TH1F* bkg_control_sample)
{
  float xmin=0.5;
  float xsig=0.92;
  using namespace RooFit;
  RooFitResult* result(0);
  
  // data
  RooRealVar x("iso","iso",0.5,1.0);
  
  // signal pdf
  RooAbsData* control_sample_reduced(0);
  RooAbsData* data = 0;
  if ( fake_type == 0 ) {
    control_sample_reduced = control_sample->reduce("fake_type==0&&iso>0.5");
    data = analysis_sample->reduce("fake_type==0&&iso>0.5");
  } else {
    control_sample_reduced = control_sample->reduce("fake_type==1&&iso>0.5");
    data = analysis_sample->reduce("fake_type==1&&iso>0.5");
  }
  
  RooDataHist control_sample_reduced_binned("ref_data","ref_data",RooArgSet(x),*control_sample_reduced);
  RooHistPdf sig_pdf("sig_pdf","sig_pdf",RooArgSet(x), control_sample_reduced_binned,2);
  
  // background pdf
  RooAbsPdf* bkg_pdf(0);
  RooRealVar bkg_a1("bkg_a1","bkg_a1",-.8);
  RooRealVar bkg_c("bkg_c","bkg_c", 1, 0.1, 10);
  RooDataHist* bkg_ref(0);
  // RooPolynomial bkg_pdf("bkg","Background",x,RooArgList(bkg_a1),0);
  switch ( pdf_type ){
  case 0:
    bkg_pdf = new RooPolynomial("bkg","Background",x,bkg_a1);
    bkg_a1.setConstant(kFALSE);
    break;
  case 1:
    bkg_pdf = new RooPolynomial("bkg","Background",x,RooArgSet());
    break;
  case 2:
    bkg_pdf = new RooGenericPdf("bkg","Background","pow((1-iso),bkg_c)",RooArgSet(bkg_c,x));
    bkg_c.setConstant(kFALSE);
    break;
  case 3:
    bkg_ref = new RooDataHist("ref_bkg","ref_bkg",RooArgSet(x),bkg_control_sample);
    bkg_pdf = new RooHistPdf("bkg","Background",RooArgSet(x), *bkg_ref,2);
    break;
  default:
    std::cout << "ERROR: unsupported backgroun pdf type" << std::endl;
  }
  
  // full pdf
  double nTotal = data->sumEntries();
  std::cout << "Total number of weighted events in the data sample: " << nTotal << std::endl;
  RooRealVar Nsig("Nsig","Nsig",0.5*nTotal,0,nTotal); 
  RooRealVar Nbkg("Nbkg","Nbkg",0.5*nTotal,0,nTotal); 
  RooAddPdf pdf("pdf","Cumulative pdf", RooArgList(sig_pdf,*bkg_pdf), RooArgList(Nsig,Nbkg));

  // fit analysis sample
  std::cout << "Fitting analysis sample" << std::endl;
  RooPlot* analysis_frame = x.frame(RooFit::Bins(20));
  Nsig.setConstant(kFALSE);
  Nbkg.setConstant(kFALSE);

  result = pdf.fitTo(*data,Extended(kTRUE),Save(kTRUE));//,RooFit::PrintLevel(2),RooFit::Verbose(kTRUE));//,RooFit::Minos(kFALSE));
  std::cout<< "Wjets background yeild estimate:" << std::endl;
  TPaveText* text = estimate_background(pdf_type, result, xmin, xsig);
  if ( pdf_type == 3 ){
    assert(bkg_control_sample);
    TH1F* h(bkg_control_sample);
    double c = h->Integral(h->FindBin(xsig),h->FindBin(1.0))/h->Integral(h->FindBin(xmin),h->FindBin(1.0));
    const char* estimate = Form("%0.1f +/- %0.1f", Nbkg.getVal()*c, Nbkg.getError()*c);
    std::cout << "\t" << estimate << std::endl;
    text->AddText(estimate);
  }

  analysis_frame->addObject(text);
  data->plotOn(analysis_frame,DataError(RooAbsData::SumW2));
  // data->statOn(analysis_frame,What("N"));
  // pdf.paramOn(analysis_frame,Format("NEA",AutoPrecision(1)));
  pdf.plotOn(analysis_frame);
  bkg_pdf->plotOn(analysis_frame,LineColor(kRed),Normalization(Nbkg.getVal(),RooAbsPdf::NumEvent));
  analysis_frame->SetTitle(title);
  analysis_frame->Draw();
  /*
    double bkg_pdf_full_integral = bkg_pdf->createIntegral(x)->getVal();
    double xmin = x.getMin();
    double xmax = x.getMax();
    x.setRange(0.92,1.0);
    double bkg_pdf_signal_integral = bkg_pdf->createIntegral(x)->getVal();
    x.setRange(xmin,xmax);
    std::cout << "\tcentral value (from integral): " << Nbkg.getVal()*bkg_pdf_signal_integral/bkg_pdf_full_integral << std::endl;
  */
  return result;
}

RooFitResult* fit_isolation(RooAbsData* full_sample, 
			    int fake_type, 
			    int pdf_type, 
			    const char* title,
			    TH1F* bkg_control_sample){
  RooAbsData* control_sample = full_sample->reduce("sample_type==1");
  if ( bkg_control_sample->GetEntries() < 1 ) {
    std::cout << "Bad background reference histogram: " << bkg_control_sample->GetName()<< std::endl;
    std::cout << "Number of entries: " << bkg_control_sample->GetEntries() << std::endl;
    exit(1);
  }
  if ( control_sample->numEntries() < 100 ){
    std::cout << "Signal control sample is too small, make sure you have Z events" << std::endl;
    exit(1);
  }
  RooAbsData* analysis_sample = full_sample->reduce("sample_type==0");
  // std::cout << "Number of entries in the signal isolation control sample: " << control_sample->numEntries() << std::endl;
  // std::cout << "Number of entries in the analysis sample: " << analysis_sample->numEntries() << std::endl;
  return fit_isolation(control_sample, analysis_sample, fake_type, pdf_type, title, bkg_control_sample);
}
