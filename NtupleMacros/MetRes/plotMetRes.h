
//#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
#include "Tools/Utilities.h"
#include "Tools/histtools.cc"

#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "/home/users/wandrews/macros/comparison.C"

void saveStack(THStack* st, TLegend* leg, bool cpylog=false) {

  TCanvas *c = new TCanvas();
  st->Draw("hist");
  leg->Draw();
  c->SaveAs((TString)st->GetName()+".png");

  if(cpylog) {
	gPad->SetLogy();
	st->Draw("hist");
	leg->Draw();
	c->SaveAs((TString)st->GetName()+"_log.png");
  }
  
}

void makeStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suffix, bool cpylog=false) {
  THStack *st = h->getStack(theSources, title, subtitle, suffix);

  saveStack( st, leg, cpylog );
}

void makeSumStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suff1, TString suff2, bool cpylog=false) {
  THStack *st = h->getSumStack(theSources, title, subtitle, suff1, suff2);

  saveStack( st, leg, cpylog );
  /*
  TCanvas *c = new TCanvas();
  st->Draw("hist");
  leg->Draw();
  c->SaveAs((TString)st->GetName()+".png");

  if(cpylog) {
	gPad->SetLogy();
	st->Draw("hist");
	leg->Draw();
	c->SaveAs((TString)st->GetName()+"_log.png");
	}*/
}

void makeSumDifStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suff1, TString suff2, TString suff3, bool cpylog=false) {
  THStack *st = h->getSumDifStack(theSources, title, subtitle, suff1, suff2, suff3);

  saveStack( st, leg, cpylog );
}

void plotResults() {

  sources_t theSources = sources_all;
  
  //HistogramUtilities h1("Results.root");
  //HistogramUtilities h1("Susy_Results.root");
  HistogramUtilities* h1 = new HistogramUtilities("Results.root");

  TString all = "all";
  TString ee = "ee"; //these correspond to defintion in DileptonHypType.h
  TString mm = "mm";
  TString em = "em"; 

  //makeStack(h1, theSources, "tcMet", "", all);
  
  TLegend *lg_all = h1->getLegend(theSources, "sumJetPt_outz", "", all);

  makeStack(h1, lg_all, theSources, "sumJetPt", "", all, true);
  makeStack(h1, lg_all, theSources, "sumJetPt_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "sumJetPt_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet_xy", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_xy_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_xy_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet_x", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_x_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_x_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet_y", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_y_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_y_inz", "", all, true);

  for(int i=0; i<4;i++) { //nsjpbins is hardcoded as 4 here
	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_xy_sjp", i), "", all, true);
	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_xy_outz_sjp", i), "", all, true);
	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_xy_inz_sjp", i), "", all, true);

	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_xy_sjp", i), "", em, true);
	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_xy_outz_sjp", i), "", em, true);
	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_xy_inz_sjp", i), "", em, true);

	makeSumStack(h1, lg_all, theSources, Form("%s%i", "tcMet_xy_sjp", i), "", ee, mm, true);
	makeSumStack(h1, lg_all, theSources, Form("%s%i", "tcMet_xy_outz_sjp", i), "", ee, mm, true);
	makeSumStack(h1, lg_all, theSources, Form("%s%i", "tcMet_xy_inz_sjp", i), "", ee, mm, true);

	makeSumDifStack(h1, lg_all, theSources, Form("%s%i", "tcMet_xy_sjp", i), "", ee, mm, em, true);
	makeSumDifStack(h1, lg_all, theSources, Form("%s%i", "tcMet_xy_outz_sjp", i), "", ee, mm, em, true);
	makeSumDifStack(h1, lg_all, theSources, Form("%s%i", "tcMet_xy_inz_sjp", i), "", ee, mm, em, true);
  }

  //do 2d hist
  TH2F* metxy = h1->get2dHistogram(theSources, "tcMet_xvy", "", all);

  //TString filename = "Results_plot.root";
  //TFile outf(filename,"RECREATE") ;
  TCanvas *c1 = new TCanvas();
  //c1->cd();
  metxy->SetName("tcMet_xvy");
  metxy->Draw("box");
  c1->SaveAs((TString)metxy->GetName()+".png");
  //metxy->Write();

  //do sum hists
  TH1F* tcMet_xy_inz = h1->getHistogramSum(theSources, "tcMet_xy_inz", "", all);
  TH1F* tcMet_xy_ouz = h1->getHistogramSum(theSources, "tcMet_xy_outz", "", all);
  over_save(tcMet_xy_ouz, tcMet_xy_inz);
  over_save(tcMet_xy_ouz, tcMet_xy_inz, true, true); //two bools set log, and put log in name
  //same two scaled
  tcMet_xy_ouz->Scale( tcMet_xy_inz->Integral()/tcMet_xy_ouz->Integral() );
  tcMet_xy_ouz->SetName("tcMet_xy_outz_scale");
  over_save(tcMet_xy_ouz, tcMet_xy_inz);
  over_save(tcMet_xy_ouz, tcMet_xy_inz, true, true);

  //do zero sjp bin
  TH1F* tcMet_xy_inz_sjp0 = h1->getHistogramSum(theSources, "tcMet_xy_inz_sjp0", "", all);
  TH1F* tcMet_xy_ouz_sjp0 = h1->getHistogramSum(theSources, "tcMet_xy_outz_sjp0", "", all);
  over_save(tcMet_xy_ouz_sjp0, tcMet_xy_inz_sjp0);
  over_save(tcMet_xy_ouz_sjp0, tcMet_xy_inz_sjp0, true, true);
  //scale
  tcMet_xy_ouz_sjp0->Scale( tcMet_xy_inz_sjp0->Integral()/tcMet_xy_ouz_sjp0->Integral() );
  tcMet_xy_ouz_sjp0->SetName("tcMet_xy_outz_sjp0_scale");
  over_save(tcMet_xy_ouz_sjp0, tcMet_xy_inz_sjp0);
  over_save(tcMet_xy_ouz_sjp0, tcMet_xy_inz_sjp0, true, true);
  
}



/* //don't bother with this one--assume that if you want a stack, you want a legend with it...
void makeStack(HistogramUtilities* h, sources_t theSources, TString title, TString subtitle, TString suffix, bool cpylog=false) {
  THStack *st = h->getStack(theSources, title, subtitle, suffix);

  TCanvas *c = new TCanvas();
  st->Draw("hist");
  c->SaveAs((TString)st->GetName()+".png");
  //delete c;
}
*/
