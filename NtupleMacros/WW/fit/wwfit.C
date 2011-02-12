#include "TFile.h"
#include "RooWorkspace.h"
#include "TProfile.h"
#include "TTree.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TSystem.h"
#include "RooNDKeysPdf.h"
#include "RooATGCPdf.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "TH2.h"
#include "TMarker.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"
// #include "RooLognormal.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooMsgService.h"

const unsigned int nEvents = 11;
const unsigned int Nbins = 48;
const double rangeX = 1.5;
const double rangeY = 1.5;

RooRealVar* var_pt1;
RooRealVar* x_par;
RooRealVar* y_par;
RooRealVar* var_dummy;

RooRealVar* n_ww;
RooRealVar* n_top;
RooRealVar* n_wjets;
RooRealVar* n_wz;
RooRealVar* n_zz;

RooATGCPdf* pdf;
RooAbsPdf*  pdf_bkg;
RooDataSet* data;
RooAbsPdf*  cpdf;
RooAbsPdf*  cSigPdf;
RooAbsPdf*  cBkgPdf;

RooGaussian* n_ww_con;
RooGaussian* n_top_con;
RooGaussian* n_wjets_con;
RooGaussian* n_wz_con;
RooGaussian* n_zz_con;

struct Sample{
  std::string file_name;
  std::string name;
  RooAbsData* dataset;
  double x;
  double y;
  TH1*   hist;
  TH1*   hist_keys;
  Sample( const char* file, const char* _name, double _x, double _y ){
    file_name = file;
    name = _name;
    x=_x;
    y=_y;
    TFile* f = TFile::Open(file);
    assert(f);
    RooAbsData* ds = (RooAbsData*)f->Get("ww");
    ds->SetName(_name);
    dataset = ds->reduce(*var_pt1,"selected==1&&unique==1");
    ((TTree*)dataset->tree())->Draw(Form("pt1>>h(%u,20,500)",Nbins),"","goff");
    hist = (TH1*)gDirectory->Get("h");
    hist->SetTitle(Form("%s: %0.2f, %s: %0.2f",x_par->GetTitle(),x,y_par->GetTitle(),y));
    hist->GetXaxis()->SetTitle("Pt, [GeV]");
    hist->SetDirectory(0);
    // for (Int_t i=0; i<=hist->GetNbinsX(); ++i)
    // hist->SetBinContent(i,hist->GetBinContent(i)+1e-5);
    hist->Scale(1/hist->Integral());
    hist->GetYaxis()->SetRangeUser(0.001,0.25);
    hist->SetStats(kFALSE);
    RooNDKeysPdf keys_pt1("keys_pt1","keys_pt1",*var_pt1,*((RooDataSet*)dataset),"am");
    hist_keys = (TH1F*)keys_pt1.createHistogram("hist_keys",*var_pt1);
    hist_keys->Scale(hist_keys->Integral());
    hist_keys->SetLineColor(kRed);
    hist_keys->SetLineWidth(2);
    hist_keys->SetDirectory(0);
    hist_keys->SetStats(kFALSE);
  }
  Sample():dataset(0),x(0),y(0),hist(0),hist_keys(0){}
};

Sample samples[11]; 

// http://root.cern.ch/root/html/RooLognormal.html
// http://en.wikipedia.org/wiki/Log-normal_distribution
double meanLogNormal(double mean, double sigma){
  return 1/sqrt(mean*mean+sigma*sigma);
}
double sigmaLogNormal(double mean, double sigma){
  return exp(sqrt(log(1+sigma*sigma/mean/mean)));
}

void setDefaults()
{
  gSystem->CompileMacro("RooATGCPdf.C","k");
  var_pt1   = new RooRealVar("pt1","pt1",20,500);
  var_dummy = new RooRealVar("var_dummy","var_dummy",0);
  var_pt1->setBins(Nbins);

  n_ww    = new RooRealVar("n_ww","n_ww", nEvents, 0, 1E+3);
  n_top   = new RooRealVar("n_top","n_top", 0, 0, 1E+3);
  n_wjets = new RooRealVar("n_wjets","n_wjets", 0, 0, 1E+3);
  n_wz    = new RooRealVar("n_wz","n_wz", 0, 0, 1E+3);
  n_zz    = new RooRealVar("n_zz","n_zz", 0, 0, 1E+3);


  n_ww_con    = new RooGaussian("n_ww_con",    "WW uncertainty",    *n_ww,    
				// RooFit::RooConst(10.77),RooFit::RooConst(1.4));   // 11% (lumi) + 5% (theory) + 5% (jet veto) = 13%
				RooFit::RooConst(nEvents),RooFit::RooConst(nEvents*0.13));
  n_top_con   = new RooGaussian("n_top_con",   "Top uncertainty",   *n_top,   
				RooFit::RooConst(0.85),RooFit::RooConst(0.85));
  n_wjets_con = new RooGaussian("n_wjets_con", "Wjets uncertainty", *n_wjets, 
				RooFit::RooConst(2.1),RooFit::RooConst(0.7));
  n_wz_con    = new RooGaussian("n_wz_con",    "WZ uncertainty",    *n_wz,    
				RooFit::RooConst(0.24),RooFit::RooConst(0.2));
  n_zz_con    = new RooGaussian("n_zz_con",    "ZZ uncertainty",    *n_zz,   
				RooFit::RooConst(0.10),RooFit::RooConst(0.01));
}

void setSigPdf_LZ_GZ()
{
  x_par = new RooRealVar("x_par","#lambda_{Z}",0);
  y_par = new RooRealVar("y_par","#Delta g^{Z}_{1}",0);
  double norm = nEvents/3.09163;
  
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(3,3);

  pdf = new RooATGCPdf("pdf", "pdf", *var_pt1, *n_ww, *x_par, *y_par);
  cSigPdf = new RooProdPdf("cSigPdf","model with constraint",RooArgSet(*pdf,*n_ww_con)) ;
  
  Int_t i=0;
  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz100_lz0_gz1000000.root",
		      "sm_sm",0,0);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*3.09163, norm*0.0547394), 
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz175_lz0_gz1750000.root",
		      "sm_p",0,0.75);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*6.12261, 0.112729*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz25_lz0_gz250000.root",
		      "sm_m",0,-0.75);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*6.01277, 0.112194*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz100_lz50_gz1000000.root",
		      "p_sm", 0.5, 0);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*5.53635, 0.103594*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz175_lz50_gz1750000.root",
		      "p_p",0.5,0.75);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*9.43854, 0.178215*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz25_lz50_gz250000.root",
		      "p_m",0.5,-0.75);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*7.22219, 0.13752*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg100_lgm50_gg100_kz100_lzm50_gz1000000.root",
		      "m_sm",-0.5,0);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*5.20204, 0.0981948*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg100_lgm50_gg100_kz175_lzm50_gz1750000.root",
		      "m_p",-0.5,0.75);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*6.89533, 0.131155*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg100_lgm50_gg100_kz25_lzm50_gz250000.root",
		      "m_m",-0.5,-0.75);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*9.72764, 0.183288*norm),
		samples[i].hist_keys);
  i++;			    
  pdf->build();
  
  c1->cd(1);
  x_par->setVal(0); y_par->setVal(0);
  TH1* h_sm_sm = pdf->createHistogram("h_sm_sm",*var_pt1);
  h_sm_sm->SetLineColor(kBlue);
  h_sm_sm->Scale(h_sm_sm->Integral());
  h_sm_sm->SetStats(kFALSE);
  h_sm_sm->Draw("same");

  c1->cd(2);
  x_par->setVal(0); y_par->setVal(0.75);
  TH1* h_sm_p = pdf->createHistogram("h_sm_p",*var_pt1);
  h_sm_p->SetLineColor(kBlue);
  h_sm_p->Scale(h_sm_p->Integral());
  h_sm_p->SetStats(kFALSE);
  h_sm_p->Draw("same");

  c1->cd(3);
  x_par->setVal(0); y_par->setVal(-0.75);
  TH1* h_sm_m = pdf->createHistogram("h_sm_m",*var_pt1);
  h_sm_m->SetLineColor(kBlue);
  h_sm_m->Scale(h_sm_m->Integral());
  h_sm_m->SetStats(kFALSE);
  h_sm_m->Draw("same");

  c1->cd(4);
  x_par->setVal(0.5); y_par->setVal(0);
  TH1* h_p_sm = pdf->createHistogram("h_p_sm",*var_pt1);
  h_p_sm->SetLineColor(kBlue);
  h_p_sm->Scale(h_p_sm->Integral());
  h_p_sm->SetStats(kFALSE);
  h_p_sm->Draw("same");

  c1->cd(5);
  x_par->setVal(0.5); y_par->setVal(0.75);
  TH1* h_p_p = pdf->createHistogram("h_p_p",*var_pt1);
  h_p_p->SetLineColor(kBlue);
  h_p_p->Scale(h_p_p->Integral());
  h_p_p->SetStats(kFALSE);
  h_p_p->Draw("same");

  c1->cd(6);
  x_par->setVal(0.5); y_par->setVal(-0.75);
  TH1* h_p_m = pdf->createHistogram("h_p_m",*var_pt1);
  h_p_m->SetLineColor(kBlue);
  h_p_m->Scale(h_p_m->Integral());
  h_p_m->SetStats(kFALSE);
  h_p_m->Draw("same");

  c1->cd(7);
  x_par->setVal(-0.5); y_par->setVal(0);
  TH1* h_m_sm = pdf->createHistogram("h_m_sm",*var_pt1);
  h_m_sm->SetLineColor(kBlue);
  h_m_sm->Scale(h_m_sm->Integral());
  h_m_sm->SetStats(kFALSE);
  h_m_sm->Draw("same");

  c1->cd(8);
  x_par->setVal(-0.5); y_par->setVal(0.75);
  TH1* h_m_p = pdf->createHistogram("h_m_p",*var_pt1);
  h_m_p->SetLineColor(kBlue);
  h_m_p->Scale(h_m_p->Integral());
  h_m_p->SetStats(kFALSE);
  h_m_p->Draw("same");
  
  c1->cd(9);
  x_par->setVal(-0.5); y_par->setVal(-0.75);
  TH1* h_m_m = pdf->createHistogram("h_m_m",*var_pt1);
  h_m_m->SetLineColor(kBlue);
  h_m_m->Scale(h_m_m->Integral());
  h_m_m->SetStats(kFALSE);
  h_m_m->Draw("same");

}

void setSigPdf_LZ_KG()
{
  x_par = new RooRealVar("x_par","#lambda_{Z}",0);
  y_par = new RooRealVar("y_par","#Delta#kappa_{#gamma}",0);
  double norm = nEvents/3.09163;
  
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(3,3);

  pdf = new RooATGCPdf("pdf", "pdf", *var_pt1, *n_ww, *x_par, *y_par);
  cSigPdf = new RooProdPdf("cSigPdf","model with constraint",RooArgSet(*pdf,*n_ww_con)) ;
  
  Int_t i=0;
  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg100_lg0_gg100_kz100_lz0_gz1000000.root",
		      "sm_sm",0,0);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*3.09163, norm*0.0547394), 
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg170_lg0_gg100_kz779311_lz0_gz100.root",
		      "sm_p",0,0.70);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*3.2987, 0.0592563*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg30_lg0_gg100_kz1220689_lz0_gz100.root",
		      "sm_m",0,-0.70);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*3.35685 , 0.0604608*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg100_lg50_gg100_kz100_lz50_gz1000000.root",
		      "p_sm", 0.5, 0);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*5.53635, 0.103594*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg170_lg50_gg100_kz779311_lz50_gz100.root",
		      "p_p",0.5,0.70);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*5.82764, 0.108893*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg170_lg100_gg100_kz779311_lz100_gz100.root",
		      "p_p",1.0,0.70);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  // FIXME - WRONG X-section
  // pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*12.8868, 0.247537*norm),
  // samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg30_lg50_gg100_kz1220689_lz50_gz100.root",
		      "p_m",0.5,-0.70);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*5.76679, 0.108426*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg30_lg100_gg100_kz1220689_lz100_gz100.root",
		      "p_m",1.0,-0.70);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  // FIXME - WRONG X-section
  // pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*12.9291, 0.248591*norm),
  // samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg100_lgm50_gg100_kz100_lzm50_gz1000000.root",
		      "m_sm",-0.5,0);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*5.20204, 0.0981948*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg170_lgm50_gg100_kz779311_lzm50_gz100.root",
		      "m_p",-0.5,0.70);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*5.76225, 0.107762*norm),
		samples[i].hist_keys);
  i++;			    

  c1->cd(i+1);
  gPad->SetLogy(kTRUE);
  samples[i] = Sample("jake-samples/processed_data_job_WW_1j_kg30_lgm50_gg100_kz1220689_lzm50_gz100.root",
		      "m_m",-0.5,-0.70);
  samples[i].hist->Draw();
  samples[i].hist_keys->Draw("same");
  pdf->addPoint(Measurement(samples[i].x, samples[i].y, norm*5.75791, 0.107823*norm),
		samples[i].hist_keys);
  i++;			    
  pdf->build();
  
  c1->cd(1);
  x_par->setVal(0); y_par->setVal(0);
  TH1* h_sm_sm = pdf->createHistogram("h_sm_sm",*var_pt1);
  h_sm_sm->SetLineColor(kBlue);
  h_sm_sm->Scale(h_sm_sm->Integral());
  h_sm_sm->SetStats(kFALSE);
  h_sm_sm->Draw("same");

  c1->cd(2);
  x_par->setVal(0); y_par->setVal(0.75);
  TH1* h_sm_p = pdf->createHistogram("h_sm_p",*var_pt1);
  h_sm_p->SetLineColor(kBlue);
  h_sm_p->Scale(h_sm_p->Integral());
  h_sm_p->SetStats(kFALSE);
  h_sm_p->Draw("same");

  c1->cd(3);
  x_par->setVal(0); y_par->setVal(-0.75);
  TH1* h_sm_m = pdf->createHistogram("h_sm_m",*var_pt1);
  h_sm_m->SetLineColor(kBlue);
  h_sm_m->Scale(h_sm_m->Integral());
  h_sm_m->SetStats(kFALSE);
  h_sm_m->Draw("same");

  c1->cd(4);
  x_par->setVal(0.5); y_par->setVal(0);
  TH1* h_p_sm = pdf->createHistogram("h_p_sm",*var_pt1);
  h_p_sm->SetLineColor(kBlue);
  h_p_sm->Scale(h_p_sm->Integral());
  h_p_sm->SetStats(kFALSE);
  h_p_sm->Draw("same");

  c1->cd(5);
  x_par->setVal(0.5); y_par->setVal(0.75);
  TH1* h_p_p = pdf->createHistogram("h_p_p",*var_pt1);
  h_p_p->SetLineColor(kBlue);
  h_p_p->Scale(h_p_p->Integral());
  h_p_p->SetStats(kFALSE);
  h_p_p->Draw("same");

  c1->cd(6);
  x_par->setVal(0.5); y_par->setVal(-0.75);
  TH1* h_p_m = pdf->createHistogram("h_p_m",*var_pt1);
  h_p_m->SetLineColor(kBlue);
  h_p_m->Scale(h_p_m->Integral());
  h_p_m->SetStats(kFALSE);
  h_p_m->Draw("same");

  c1->cd(7);
  x_par->setVal(-0.5); y_par->setVal(0);
  TH1* h_m_sm = pdf->createHistogram("h_m_sm",*var_pt1);
  h_m_sm->SetLineColor(kBlue);
  h_m_sm->Scale(h_m_sm->Integral());
  h_m_sm->SetStats(kFALSE);
  h_m_sm->Draw("same");

  c1->cd(8);
  x_par->setVal(-0.5); y_par->setVal(0.75);
  TH1* h_m_p = pdf->createHistogram("h_m_p",*var_pt1);
  h_m_p->SetLineColor(kBlue);
  h_m_p->Scale(h_m_p->Integral());
  h_m_p->SetStats(kFALSE);
  h_m_p->Draw("same");
  
  c1->cd(9);
  x_par->setVal(-0.5); y_par->setVal(-0.75);
  TH1* h_m_m = pdf->createHistogram("h_m_m",*var_pt1);
  h_m_m->SetLineColor(kBlue);
  h_m_m->Scale(h_m_m->Integral());
  h_m_m->SetStats(kFALSE);
  h_m_m->Draw("same");

}

void setBkgPdf()
{
  TCanvas* c9 = new TCanvas("c9","c9",600,900);
  c9->Divide(2,3);

  TFile* f = TFile::Open("samples/processed_data_final.root");
  assert(f);
  
  c9->cd(1);
  // ww
  RooAbsData* ds_ww = (RooAbsData*)f->Get("ww");
  ds_ww->SetName("ds_ww");
  RooAbsData* ds_ww_pt1 = ds_ww->reduce(*var_pt1,"selected==1&&unique==1");
  RooNDKeysPdf* pdf_ww = new RooNDKeysPdf("pdf_ww","pdf_ww",*var_pt1,*((RooDataSet*)ds_ww_pt1),"am");
  TH1F* hpdf_ww = (TH1F*)pdf_ww->createHistogram("hpdf_ww",*var_pt1);
  hpdf_ww->SetTitle("WW");
  hpdf_ww->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_ww->Scale(hpdf_ww->Integral());
  hpdf_ww->Draw();

  c9->cd(2);
  // wjets
  RooAbsData* ds_wjets = (RooAbsData*)f->Get("wjets");
  ds_wjets->SetName("ds_wjets");
  RooAbsData* ds_wjets_pt1 = ds_wjets->reduce(*var_pt1,"selected==1&&unique==1");
  RooNDKeysPdf* pdf_wjets = new RooNDKeysPdf("pdf_wjets","pdf_wjets",*var_pt1,*((RooDataSet*)ds_wjets_pt1),"am");
  TH1F* hpdf_wjets = (TH1F*)pdf_wjets->createHistogram("hpdf_wjets",*var_pt1);
  hpdf_wjets->SetTitle("WJets");
  hpdf_wjets->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_wjets->Scale(hpdf_wjets->Integral());
  hpdf_wjets->Draw();

  c9->cd(3);
  // ttbar
  RooAbsData* ds_ttbar = (RooAbsData*)f->Get("ttbar");
  ds_ttbar->SetName("ds_ttbar");
  RooAbsData* ds_ttbar_pt1 = ds_ttbar->reduce(*var_pt1,"selected==1&&unique==1");
  RooNDKeysPdf* pdf_ttbar = new RooNDKeysPdf("pdf_ttbar","pdf_ttbar",*var_pt1,*((RooDataSet*)ds_ttbar_pt1),"am");
  TH1F* hpdf_ttbar = (TH1F*)pdf_ttbar->createHistogram("hpdf_ttbar",*var_pt1);
  hpdf_ttbar->SetTitle("TTbar");
  hpdf_ttbar->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_ttbar->Scale(hpdf_ttbar->Integral());
  hpdf_ttbar->Draw();

  c9->cd(4);
  // tW
  RooAbsData* ds_tw = (RooAbsData*)f->Get("tw");
  ds_tw->SetName("ds_tw");
  RooAbsData* ds_tw_pt1 = ds_tw->reduce(*var_pt1,"selected==1&&unique==1");
  RooNDKeysPdf* pdf_tw = new RooNDKeysPdf("pdf_tw","pdf_tw",*var_pt1,*((RooDataSet*)ds_tw_pt1),"am");
  TH1F* hpdf_tw = (TH1F*)pdf_tw->createHistogram("hpdf_tw",*var_pt1);
  hpdf_tw->SetTitle("tW");
  hpdf_tw->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_tw->Scale(hpdf_tw->Integral());
  hpdf_tw->Draw();
  
  c9->cd(5);
  // wz
  RooAbsData* ds_wz = (RooAbsData*)f->Get("wz");
  ds_wz->SetName("ds_wz");
  RooAbsData* ds_wz_pt1 = ds_wz->reduce(*var_pt1,"selected==1&&unique==1");
  RooNDKeysPdf* pdf_wz = new RooNDKeysPdf("pdf_wz","pdf_wz",*var_pt1,*((RooDataSet*)ds_wz_pt1),"am");
  TH1F* hpdf_wz = (TH1F*)pdf_wz->createHistogram("hpdf_wz",*var_pt1);
  hpdf_wz->SetTitle("WZ");
  hpdf_wz->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_wz->Scale(hpdf_wz->Integral());
  hpdf_wz->Draw();
  
  c9->cd(6);
  // zz
  RooAbsData* ds_zz = (RooAbsData*)f->Get("zz");
  ds_zz->SetName("ds_zz");
  RooAbsData* ds_zz_pt1 = ds_zz->reduce(*var_pt1,"selected==1&&unique==1");
  RooNDKeysPdf* pdf_zz = new RooNDKeysPdf("pdf_zz","pdf_zz",*var_pt1,*((RooDataSet*)ds_zz_pt1),"am");
  TH1F* hpdf_zz = (TH1F*)pdf_zz->createHistogram("hpdf_zz",*var_pt1);
  hpdf_zz->SetTitle("ZZ");
  hpdf_zz->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_zz->Scale(hpdf_zz->Integral());
  hpdf_zz->Draw();

  // extended pdfs
  RooAbsPdf* epdf_wjets = new RooExtendPdf("epdf_wjets","epdf_wjets",*pdf_wjets,*n_wjets);
  RooAbsPdf* epdf_top = new RooExtendPdf("epdf_top","epdf_top",*pdf_ttbar,*n_top);
  RooAbsPdf* epdf_wz = new RooExtendPdf("epdf_wz","epdf_wz",*pdf_wz,*n_wz);
  RooAbsPdf* epdf_zz = new RooExtendPdf("epdf_zz","epdf_zz",*pdf_ttbar,*n_zz);
  pdf_bkg = new RooAddPdf("pdf_bkg","pdf_bkg",RooArgList(*epdf_wjets,*epdf_top, *epdf_wz, *epdf_zz));
  cBkgPdf = new RooProdPdf("cBkgPdf","model with constraint",RooArgSet(*pdf_bkg,*n_top_con,*n_wjets_con,*n_wz_con,*n_zz_con)) ;
}

void sensitivity()
{
  TFile *f = TFile::Open("atgc.root");
  assert(f);
  RooWorkspace *ww = (RooWorkspace*)f->Get("ww");
  assert(ww);
  
  const unsigned int nPoints = 21;
  std::vector<double> sum(nPoints,0);
  std::vector<double> sum2(nPoints,0);
    
  // loop over WW dataset and make small samples
  // perform Likelihood difference calculation for each
  // value of the parameter
  unsigned int n(0);
  TTree* iTree = const_cast<TTree*>(((RooDataSet*)ww->data("ds_ww_pt1"))->tree());
  assert(iTree);
  unsigned int size = iTree->GetEntries();
  cout << "iTree->GetEntries(): " << size << endl;
  unsigned int firstEntry = 0;
  y_par->setVal(0);
  data=0;
  while ( firstEntry <= size - nEvents ){
    n++;
    TTree* tree = iTree->CopyTree("","",nEvents,firstEntry);
    assert(tree);
    firstEntry+=nEvents;
    if (data) delete data;
    data = new RooDataSet("data","data",tree, *ww->var("pt1"));
    // continue;
    for ( unsigned int i = 0; i<nPoints; ++i ){
      x_par->setVal(0);
      RooAbsReal* nll = cSigPdf->createNLL(*data,RooFit::Extended(),RooFit::Constrain(*n_ww));
      nll->addServer(*var_dummy);
      double nll_sm(nll->getVal());
      delete nll;
      x_par->setVal(-0.5+i*0.05);
      nll = cSigPdf->createNLL(*data,RooFit::Extended(),RooFit::Constrain(*n_ww));
      nll->addServer(*var_dummy);
      double nll_atgc(nll->getVal());
      delete nll;
      double sig = sqrt(2.*fabs(nll_sm-nll_atgc));
      sum[i]+=sig;
      sum2[i]+=sig*sig;
    }
  }
  // return;

  // p->Draw();
  if (n==0) return;
  // show results
  TCanvas* c2 = new TCanvas("c2","Significance",500,500);
  //c1->SetFillColor(42);
  c2->SetGrid();
  Float_t xval[nPoints];
  Float_t yval[nPoints];
  Float_t yerr[nPoints];
  for ( unsigned int i = 0; i<nPoints; ++i ){
    xval[i] = -0.5+i*0.05;
    yval[i] = sum[i]/n;
    yerr[i] = sqrt(sum2[i]/n-pow(sum[i]/n,2));
  }
  TGraphErrors* gr = new TGraphErrors(nPoints,xval,yval,(Float_t*)0,yerr);
  // gr->SetTitle("TGraphErrors Example");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("AP");

  x_par->setRange(-rangeX,rangeX);
  y_par->setRange(-rangeY,rangeY);
  y_par->setConstant(0);
  x_par->setConstant(0);

  // fit 
//   {
//     RooAbsReal* nll = cSigPdf->createNLL(*data,RooFit::Extended(),RooFit::Constrain(*n_ww));
//     RooMinuit m(*nll);
//     m.migrad();
//     m.hesse();
//     RooPlot* p1 = m.contour(*x_par,*y_par);
//     m.migrad();
//     m.hesse();
//     TCanvas* c4 = new TCanvas("c4","c4",500,500);
//     c4->SetGrid();
//     p1->Draw();
//   }

  // fit the reference point samples
  TCanvas* c5 = new TCanvas("c5","c5",800,800);
  c5->Divide(3,3);
  TCanvas* c6 = new TCanvas("c6","c6",800,800);
  c6->Divide(3,3);
  for (unsigned int i=0; i<9; ++i){
    x_par->setVal(samples[i].x);
    y_par->setVal(samples[i].y);
    Int_t n = cSigPdf->expectedEvents(RooArgSet());
    cout << "i: " << i << " \tExpected number of events: " << n << endl;
    TTree* iTree = const_cast<TTree*>(((RooDataSet*)samples[i].dataset)->tree());
    assert(iTree);
    TTree* tree = iTree->CopyTree("","",n,0);
    assert(tree);
    RooDataSet*ds = new RooDataSet("ds","ds",tree, *ww->var("pt1"));
    c5->cd(i+1);
    gPad->SetGrid();
    // x_par->setVal(0);
    // y_par->setVal(0);
    RooAbsReal* nll = cSigPdf->createNLL(*ds,RooFit::Extended(),RooFit::Constrain(*n_ww));
    double minNLL = nll->getVal();
    RooMinuit m(*nll);
    m.migrad();
    m.hesse();
    TMarker* markerFit = new TMarker(x_par->getVal(),y_par->getVal(),3);
    RooPlot* p = m.contour(*x_par,*y_par);
    m.migrad();
    m.hesse();
    p->Draw();
    c6->cd(i+1);
    const Int_t nBins = 200;
    RooMsgService::instance().setStreamStatus(1,0);
    TH2F* nll2d = new TH2F("nll2d",Form("#sqrt{2 #Delta logL}: %s: %0.2f, %s: %0.2f",x_par->GetTitle(),samples[i].x,x_par->GetTitle(),samples[i].y),
			   nBins,-1,1,nBins,-1,1);
    nll2d->GetXaxis()->SetTitle(x_par->GetTitle());
    nll2d->GetYaxis()->SetTitle(y_par->GetTitle());
    nll2d->SetDirectory(0);
    for ( Int_t xi = 1; xi<=nBins; ++xi ){
      cout << ".";
      cout.flush();
      for ( Int_t yi = 1; yi<=nBins; ++yi )
	{
	  x_par->setVal(-0.95+2.0*(xi-1)/nBins);
	  y_par->setVal(-0.95+2.0*(yi-1)/nBins);
	  RooAbsReal* nll = cSigPdf->createNLL(*ds,RooFit::Extended(),RooFit::Constrain(*n_ww));
	  nll->addServer(*var_dummy);
	  double sig = sqrt(2*fabs(nll->getVal()-minNLL));
	  nll2d->SetBinContent(xi,yi,sig);
	  delete nll;
	}
    }
    nll2d->Draw("colz");
    RooMsgService::instance().setStreamStatus(1,1);
    TMarker* marker = new TMarker(samples[i].x,samples[i].y,20);
    // marker->SetColor(kBlue);
    marker->Draw();
    markerFit->Draw();
  }
  
}

void ww1DFits()
{
  TFile* f = TFile::Open("samples/processed_data_final.root");
  assert(f);
  RooAbsData* ds = (RooAbsData*)f->Get("ww");
  ds->SetName("ds_ww");
  RooAbsData* dataset = ds->reduce(*var_pt1,"selected==1&&unique==1");

  //TFile *f = TFile::Open("atgc.root");
  // assert(f);
  // RooWorkspace *ww = (RooWorkspace*)f->Get("ww");
  // assert(ww);
  
  // loop over WW dataset and make small samples
  // perform Likelihood difference calculation for each
  // value of the parameter
  unsigned int n(0);
  // TTree* iTree = const_cast<TTree*>(((RooDataSet*)ww->data("ds_ww_pt1"))->tree());
  TTree* iTree = const_cast<TTree*>(((RooDataSet*)dataset)->tree());
  assert(iTree);
  unsigned int size = iTree->GetEntries();
  cout << "iTree->GetEntries(): " << size << endl;
  unsigned int firstEntry = 0;
  y_par->setVal(0);
  data=0;
  // TRandom::Poisson

  TH1F* hx = new TH1F("hx","Fit of WW SM events",40,-1,1);
  hx->GetXaxis()->SetTitle("#lambda_Z");
  TH1F* hy = new TH1F("hy","Fit of WW SM events",40,-1,1);
  hy->GetXaxis()->SetTitle("#Delta#kappa_Z");
  while ( firstEntry <= size - nEvents ){
    n++;
    TTree* tree = iTree->CopyTree("","",nEvents,firstEntry);
    assert(tree);
    firstEntry+=nEvents;
    if (data) delete data;
    data = new RooDataSet("data","data",tree, *var_pt1);
    
    x_par->setVal(0);
    x_par->setConstant(0);
    y_par->setVal(0);
    y_par->setConstant(1);
    pdf->fitTo(*data);
    hx->Fill(x_par->getVal());

    x_par->setVal(0);
    x_par->setConstant(1);
    y_par->setVal(0);
    y_par->setConstant(0);
    pdf->fitTo(*data);
    hy->Fill(y_par->getVal());
  }

  // show results
  TCanvas* c7 = new TCanvas("c7","Fit results",1000,500);
  c7->Divide(2,1);
  c7->cd(1);
  hx->Draw();
  c7->cd(2);
  hy->Draw();
}

TH1F* wwATGC1DFit(const char* file, const char* name, double lz, double dkz)
{
  x_par->setRange(-rangeX,rangeX);
  y_par->setRange(-rangeY,rangeY);
  assert(lz==0||dkz==0);
  TFile* f = TFile::Open(file);
  assert(f);
  RooAbsData* ds = (RooAbsData*)f->Get("ww");
  ds->SetName(name);
  RooAbsData* dataset = ds->reduce(*var_pt1,"selected==1&&unique==1");
  
  // loop over WW dataset and make small samples
  // perform Likelihood difference calculation for each
  // value of the parameter
  unsigned int n(0);
  TTree* iTree = const_cast<TTree*>(((RooDataSet*)dataset)->tree());
  assert(iTree);
  unsigned int size = iTree->GetEntries();
  cout << "iTree->GetEntries(): " << size << endl;
  unsigned int firstEntry = 0;
  y_par->setVal(0);
  data=0;
  // TRandom::Poisson

  // TH1F* h = new TH1F(Form("h_%s",name),"Fit of on WW with anomalous couplings",40,-1,1);
  TH1F* h = new TH1F(Form("h_%s",name),"Fit of on WW with anomalous couplings",40,0,
		     std::max(rangeX,rangeY));
  h->SetDirectory(0);
  while ( firstEntry <= size - nEvents ){
    n++;
    TTree* tree = iTree->CopyTree("","",nEvents,firstEntry);
    assert(tree);
    firstEntry+=nEvents;
    if (data) delete data;
    data = new RooDataSet("data","data",tree, *var_pt1);
    
    x_par->setVal(lz);
    y_par->setVal(dkz);
    if ( lz==0 ){
      x_par->setConstant(1);
      y_par->setConstant(0);
    } else {
      x_par->setConstant(0);
      y_par->setConstant(1);
    }      
    pdf->fitTo(*data);
    if ( lz==0 )
      h->Fill(fabs(y_par->getVal()));
    else
      h->Fill(fabs(x_par->getVal()));
  }
  return h;
}

void wwATGC1DFits()
{
  TCanvas* c8 = new TCanvas("c8","Fit results",800,800);
  c8->Divide(2,2);
  
  c8->cd(1);
  TH1F* h1 = wwATGC1DFit("processed_data_WW_1j_kg100_lg0_gg100_kz175_lz0_gz175_fastsim386_v1.root",
			 "sm_p",0,0.75);
  h1->SetTitle("#Delta#kappa_{Z} 0.75");
  h1->GetXaxis()->SetTitle("|#Delta#kappa_{Z}|");
  h1->Draw();

  c8->cd(2);
  TH1F* h2 = wwATGC1DFit("processed_data_WW_1j_kg100_lg0_gg100_kz25_lz0_gz25_fastsim386_v1.root",
			 "sm_m",0,-0.75);
  h2->SetTitle("#Delta#kappa_{Z} -0.75");
  h2->GetXaxis()->SetTitle("|#Delta#kappa_{Z}|");
  h2->Draw();

  c8->cd(3);
  TH1F* h3 = wwATGC1DFit("processed_data_WW_1j_kg100_lg50_gg100_kz100_lz50_gz100_fastsim386_v1.root",
			 "p_sm",0.5,0);
  h3->SetTitle("#lambda_{Z} 0.5");
  h3->GetXaxis()->SetTitle("|#lambda_{Z}|");
  h3->Draw();
    
  c8->cd(4);
  TH1F* h4 = wwATGC1DFit("processed_data_WW_1j_kg100_lgm50_gg100_kz100_lzm50_gz100_fastsim386_v1.root",
			 "p_sm",-0.5,0);
  h4->SetTitle("#lambda_{Z} -0.5");
  h4->GetXaxis()->SetTitle("|#lambda_{Z}|");
  h4->Draw();
}

void wwATGC1DLzKgFits()
{
  TCanvas* c8 = new TCanvas("c8","Fit results",800,400);
  c8->Divide(2,1);
  
  c8->cd(1);
  TH1F* h1 = wwATGC1DFit("samples/processed_data_WW_1j_kg170_lg0_gg100_kz779311_lz0_gz100.root",
			 "sm_p",0,0.7);
  h1->SetTitle(Form("%s 0.7",y_par->GetTitle()));
  h1->GetXaxis()->SetTitle(y_par->GetTitle());
  h1->Draw();

  c8->cd(2);
  TH1F* h2 = wwATGC1DFit("samples/processed_data_WW_1j_kg30_lg0_gg100_kz1220689_lz0_gz100.root",
			 "sm_m",0,-0.7);
  h2->SetTitle(Form("%s -0.7",y_par->GetTitle()));
  h2->GetXaxis()->SetTitle(y_par->GetTitle());
  h2->Draw();
}

void fitData()
{
  TCanvas* c10 = new TCanvas("c10","c10",1000,500);
  c10->Divide(2,1);
  TFile* f = TFile::Open("samples/processed_data_final.root");
  assert(f);
  RooAbsData* ds_data = (RooAbsData*)f->Get("data");
  ds_data->SetName("ds_data");
  RooAbsData* ds_data_pt1 = ds_data->reduce(*var_pt1,"selected==1&&unique==1");
  c10->cd(1);
  ((TTree*)ds_data_pt1->tree())->Draw(Form("pt1>>h(%u,20,500)",Nbins));
  data = dynamic_cast<RooDataSet*>(ds_data_pt1);
  RooNDKeysPdf* pdf_data = new RooNDKeysPdf("pdf_data","pdf_data",*var_pt1,*((RooDataSet*)ds_data_pt1),"am");
  TH1F* hpdf_data = (TH1F*)pdf_data->createHistogram("hpdf_data",*var_pt1);
  hpdf_data->SetTitle("DATA");
  hpdf_data->GetXaxis()->SetTitle("Leading lepton pt, [GeV]");
  hpdf_data->Scale(hpdf_data->Integral());
  c10->cd(2);
  hpdf_data->Draw();
  
  cpdf = new RooAddPdf("cpdf","combined pdf",RooArgList(*cSigPdf,*cBkgPdf));
  
  RooArgSet constrainedParams(*n_ww,*n_top,*n_wjets,*n_wz,*n_zz);

  if (1) {
    n_top->setConstant(0);
    n_wjets->setConstant(0);
    n_wz->setConstant(0);
    n_zz->setConstant(0);
  }

  x_par->setRange(-rangeX,rangeX);
  y_par->setRange(-rangeY,rangeY);
  y_par->setConstant(0);
  x_par->setConstant(0);
  y_par->setVal(0);
  x_par->setVal(0);

  // signal only pdf
  {
    TCanvas* c11 = new TCanvas("c11","c11",500,500);
    c11->SetGrid();
    y_par->setVal(0);
    x_par->setVal(0);
    RooAbsReal* nll = pdf->createNLL(*ds_data_pt1,RooFit::Extended(),RooFit::Constrain(constrainedParams));
    RooMinuit m(*nll);
    m.migrad();
    m.hesse();
    RooPlot* p1 = m.contour(*x_par,*y_par,sqrt(2.3),sqrt(6.0));
    p1->SetTitle("68% and 95% C.L.");
    m.migrad();
    m.hesse();
    p1->Draw();
  }

  // signal + background pdf 
  {
    TCanvas* c12 = new TCanvas("c12","c12",500,500);
    c12->SetGrid();
    y_par->setVal(0);
    x_par->setVal(0);
    RooAbsReal* nll = cpdf->createNLL(*ds_data_pt1,RooFit::Extended(),RooFit::Constrain(constrainedParams));
    RooMinuit m(*nll);
    m.migrad();
    m.hesse();
    RooPlot* p1 = m.contour(*x_par,*y_par,sqrt(2.3),sqrt(6.0));
    p1->SetTitle("68% and 95% C.L.");
    m.migrad();
    m.hesse();
    p1->Draw();
  }

  if (0){
    RooMsgService::instance().setStreamStatus(1,0);
    // fits
    x_par->setConstant(0);
    y_par->setConstant(0);
    x_par->setVal(0);
    y_par->setVal(0);
    cpdf->fitTo(*data);
    double minNLL = cpdf->createNLL(*data,RooFit::Extended(),RooFit::Constrain(constrainedParams))->getVal();
    TMarker* marker = new TMarker(x_par->getVal(),y_par->getVal(),20);

    TCanvas* c12_2 = new TCanvas("c12_2","c12_2",500,500);
    c12_2->SetGrid();
    const Int_t nBins = 200;
    TH2F* nll2d = new TH2F("nll2d","#sqrt{2 #Delta logL} for 35.5/pb of data", nBins,-1,1,nBins,-1,1);
    nll2d->GetXaxis()->SetTitle(x_par->GetTitle());
    nll2d->GetYaxis()->SetTitle(y_par->GetTitle());
    nll2d->SetDirectory(0);
    for ( Int_t xi = 1; xi<=nBins; ++xi ){
      cout << ".";
      cout.flush();
      for ( Int_t yi = 1; yi<=nBins; ++yi )
	{
	  x_par->setVal(-0.95+2.0*(xi-1)/nBins);
	  y_par->setVal(-0.95+2.0*(yi-1)/nBins);
	  RooAbsReal* nll = cpdf->createNLL(*ds_data_pt1,RooFit::Extended(),RooFit::Constrain(constrainedParams)/*,RooFit::NumCPU(2),RooFit::Verbose(0)*/);
	  nll->addServer(*var_dummy);
	  // double sig = sqrt(2*fabs(nll->getVal()-minNLL));
	  nll2d->SetBinContent(xi,yi,2*fabs(nll->getVal()-minNLL));
	  delete nll;
	}
    }
    nll2d->SetStats(kFALSE);
    nll2d->Draw("colz");
    // marker->SetMarkerColor(kBlue);
    marker->Draw();
    RooMsgService::instance().setStreamStatus(1,1);
  }

  if (0){
    // yield fit
    x_par->setConstant(1);
    y_par->setConstant(1);
    x_par->setVal(0);
    y_par->setVal(0);
    // cpdf->fitTo(*ds_data_pt1, RooFit::Minos());
    RooAbsReal* nll = cpdf->createNLL(*ds_data_pt1,RooFit::Extended(),RooFit::Constrain(RooArgSet(*n_top,*n_wjets,*n_wz,*n_zz)));
    RooMinuit m(*nll);
    m.setErrorLevel(0.5); //68% C.L. 
    m.migrad();
    m.hesse();
    m.minos();
  }
  
  {
    // fits
    x_par->setConstant(0);
    y_par->setConstant(1);
    x_par->setVal(0);
    y_par->setVal(0);
    // cpdf->fitTo(*ds_data_pt1, RooFit::Minos());
    RooAbsReal* nll = cpdf->createNLL(*ds_data_pt1,RooFit::Extended(),RooFit::Constrain(constrainedParams));
    RooMinuit m(*nll);
    m.setErrorLevel(3.84*0.5); //95% C.L. 
    m.migrad();
    m.hesse();
    m.minos();
  }

  {
    x_par->setConstant(1);
    y_par->setConstant(0);
    x_par->setVal(0);
    y_par->setVal(0);
    // cpdf->fitTo(*ds_data_pt1, RooFit::Minos());
    RooAbsReal* nll = cpdf->createNLL(*ds_data_pt1,RooFit::Extended(),RooFit::Constrain(constrainedParams));
    RooMinuit m(*nll);
    m.setErrorLevel(3.84*0.5); //95% C.L. 
    m.migrad();
    m.hesse();
    m.minos();
  }
}

void fitTop()
{
  TFile* f = TFile::Open("samples/processed_data_final.root");
  assert(f);
  RooAbsData* ds_ttbar = (RooAbsData*)f->Get("ttbar");
  ds_ttbar->SetName("ds_ttbar");
  RooAbsData* ds_ttbar_pt1 = ds_ttbar->reduce(*var_pt1,"selected==1&&unique==1");

  TTree* iTree = const_cast<TTree*>(((RooDataSet*)ds_ttbar_pt1)->tree());
  assert(iTree);
  unsigned int size = iTree->GetEntries();
  cout << "iTree->GetEntries(): " << size << endl;
  TTree* tree = iTree->CopyTree("","",nEvents,0);
  assert(tree);
  data = new RooDataSet("data","data",tree, *var_pt1);

  RooAbsPdf* cpdf = new RooAddPdf("cpdf","combined pdf",RooArgList(*pdf,*pdf_bkg));

  x_par->setRange(-rangeX,rangeX);
  y_par->setRange(-rangeY,rangeY);
  y_par->setConstant(0);
  x_par->setConstant(0);
  y_par->setVal(0);
  x_par->setVal(0);

  // signal + background pdf 
//   {
//     TCanvas* c13 = new TCanvas("c13","c13",500,500);
//     c13->SetGrid();
//     y_par->setVal(0);
//     x_par->setVal(0);
//     RooAbsReal* nll = cpdf->createNLL(*data,RooFit::Extended());
//     RooMinuit m(*nll);
//     m.migrad();
//     m.hesse();
//     RooPlot* p1 = m.contour(*x_par,*y_par);
//     m.migrad();
//     m.hesse();
//     p1->Draw();
//   }

  {
    // fits
    x_par->setConstant(0);
    y_par->setConstant(0);
    x_par->setVal(0);
    y_par->setVal(0);
    cpdf->fitTo(*data);
    double minNLL = cpdf->createNLL(*data)->getVal();
    TMarker* marker = new TMarker(x_par->getVal(),y_par->getVal(),20);

    TCanvas* c13_2 = new TCanvas("c13_2","c13_2",500,500);
    c13_2->SetGrid();
    const Int_t nBins = 200;
    TH2F* nll2d = new TH2F("nll2d","TTbar Monte Carlo as WW aTGC", nBins,-1,1,nBins,-1,1);
    nll2d->SetTitle("2#Delta#ln L");
    nll2d->GetXaxis()->SetTitle("|#lambda_{Z}|");
    nll2d->GetYaxis()->SetTitle("|#Delta#kappa_{Z}|");
    nll2d->SetDirectory(0);
    for ( Int_t xi = 1; xi<=nBins; ++xi ){
      cout << ".";
      cout.flush();
      for ( Int_t yi = 1; yi<=nBins; ++yi )
	{
	  x_par->setVal(-0.95+2.0*(xi-1)/nBins);
	  y_par->setVal(-0.95+2.0*(yi-1)/nBins);
	  RooAbsReal* nll = cpdf->createNLL(*data,RooFit::Extended());
	  nll->addServer(*var_dummy);
	  // double sig = sqrt(2*fabs(nll->getVal()-minNLL));
	  nll2d->SetBinContent(xi,yi,2*fabs(nll->getVal()-minNLL));
	  delete nll;
	}
    }
    nll2d->Draw("colz");
    // marker->SetMarkerColor(kBlue);
    marker->Draw();
  }

  // fits
  x_par->setConstant(0);
  y_par->setConstant(1);
  x_par->setVal(0);
  y_par->setVal(0);
  cpdf->fitTo(*data, RooFit::Minos());
  
  x_par->setConstant(1);
  y_par->setConstant(0);
  x_par->setVal(0);
  y_par->setVal(0);
  cpdf->fitTo(*data, RooFit::Minos());

}
