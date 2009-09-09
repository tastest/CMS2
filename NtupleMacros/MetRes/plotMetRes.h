
//#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
//#include "Tools/Utilities.h"
//#include "Tools/histtools.cc"

#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "/home/users/wandrews/macros/comparison.C"

bool verbose = false;

//void saveStack(THStack* st, TLegend* leg, bool cpylog=false, double ymin=10., double ymax=-999, TString name="") {
void saveStack(THStack* st, TLegend* leg, bool cpylog, double ymin, double ymax, TString name="") {

  double oldymax = st->GetMaximum();
  TCanvas *c = new TCanvas();
  TString usename;
  if( name != "" ) {
	usename = name;//+".png";
	//st->SetName(name);
	st->SetTitle(name);
  }
  else
	usename = (TString)st->GetName();//+".png";

  //for the y axis won't give me a max...just gives 1--don't use GetYaxis for stacks, use GetMaximum
  //cout << "ymin " << ymin << "  ymax " << ymax << "  st y max " << st->GetYaxis()->GetXmax() << "  max " << st->GetMaximum() << endl;
  if( ymax != -999 ) {
	//st->SetMinimum( ymin );
	st->SetMaximum( ymax );
  }
  //else
  //st->SetMinimum( ymin );
  //st->GetYaxis()->SetRangeUser( ymin, st->GetMaximum()+0.01*st->GetMaximum() );

  st->Draw("hist"); //must draw before the axis range is defined...(uh, don't ask me)

  if( verbose )
	cout << "stack " << usename << "  max " << st->GetMaximum() << "  " << oldymax << "  min " << st->GetMinimum() << endl;
  //c->Update();
  //getc(stdin);
  leg->Draw();
  c->SaveAs(usename+".png");
  //st->SetMaximum( oldymax );
  st->SetMinimum( ymin );
  
  if(cpylog) {
	gPad->SetLogy();

	st->Draw("hist");
 	leg->Draw();
	c->Update();
	//c->SaveAs((TString)st->GetName()+"_log.png");
	c->SaveAs(usename+"_log.png");
  }
  
}

void makeStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suffix, bool cpylog=false, double ymin=10., double ymax=-999, TString name="") {
  THStack *st = h->getStack(theSources, title, subtitle, suffix);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}

void makeSumStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suff1, TString suff2, bool cpylog=false, TString title2="", double scale=1.0, TString name="", double ymin=10., double ymax=-999) {
  THStack *st = h->getSumStack(theSources, title, subtitle, suff1, suff2, 1, title2, scale);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}

void make2fileStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suff1, bool cpylog=false, double ymin=10., double ymax=-999, TString name="") {
  THStack *st = h->get2fileStack(theSources, title, subtitle, suff1);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}

void makeSumDifStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suff1, TString suff2, TString suff3, bool cpylog=false, TString name="", double ymin=10., double ymax=-999) {
  THStack *st = h->getSumDifStack(theSources, title, subtitle, suff1, suff2, suff3);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}

void makeHistOverlay( HistogramUtilities* h, sources_t theSources1,  sources_t theSources2, TString title1, TString title2, TString subtitle, TString suff1, TString suff2, bool cpylog=false, bool cpyscale=false, bool cpydiff=false, TString nameprefix="", double ymin=10., double ymax=-999) {
  TH1F* h1 = h->getHistogram(theSources1, title1, "", suff1, 1, nameprefix);
  TH1F* h2 = h->getHistogram(theSources2, title2, "", suff2, 1, nameprefix);
  over_save(h2, h1);
  if( cpylog )
	over_save(h2, h1, true, true); //two bools set log, and put log in name
  //same two scaled
  if( cpyscale ) {
	h2->Scale( h1->Integral()/h2->Integral() );
	h2->SetName((TString)h2->GetName() + "_scale");
	over_save(h2, h1);
	if( cpylog )
	  over_save(h2, h1, true, true);
  }
  //same two difference--keep the scaling if cpyscale
  if( cpydiff ) {
	h1->Add( h2, -1. );
	h2->SetName( (TString)h2->GetName() + "_diff" ); //change name on h2 for consistency with above
	TCanvas *c1 = new TCanvas();
	h1->Draw();
	c1->SaveAs( (TString)h2->GetName() + ".png" );
  }
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

  makeStack(h1, lg_all, theSources, "tcMet", "", all, true); //magnitude of met
  makeStack(h1, lg_all, theSources, "tcMet_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet_x", "", all, true); //x comp
  makeStack(h1, lg_all, theSources, "tcMet_x_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_x_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet_y", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_y_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_y_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet_xy", "", all, true); //each x and y comps seperately
  makeStack(h1, lg_all, theSources, "tcMet_xy_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_xy_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet_xy", "", em, true);
  makeStack(h1, lg_all, theSources, "tcMet_xy_outz", "", em, true);
  makeStack(h1, lg_all, theSources, "tcMet_xy_inz", "", em, true);

  makeSumStack(h1, lg_all, theSources, "tcMet_xy", "", ee, mm, true);
  makeSumStack(h1, lg_all, theSources, "tcMet_xy_outz", "", ee, mm, true);
  makeSumStack(h1, lg_all, theSources, "tcMet_xy_inz", "", ee, mm, true);

  makeSumStack(h1, lg_all, theSources, "tcMet_xy"     , "", all, all, true, "genMet_xy", -1., "tcMet-genMet_xy");
  makeSumStack(h1, lg_all, theSources, "tcMet_xy_outz", "", all, all, true, "genMet_xy", -1., "tcMet-genMet_xy_outz");
  makeSumStack(h1, lg_all, theSources, "tcMet_xy_inz" , "", all, all, true, "genMet_xy", -1., "tcMet-genMet_xy_inz"); 

  makeSumDifStack(h1, lg_all, theSources, "tcMet_xy", "", ee, mm, em, true);
  makeSumDifStack(h1, lg_all, theSources, "tcMet_xy_outz", "", ee, mm, em, true);
  makeSumDifStack(h1, lg_all, theSources, "tcMet_xy_inz", "", ee, mm, em, true);

  //verbose = true;
  //cout << "setting verbose = true" << endl;
  //h1->setVerbose( true );
  //makeSumDifStack(h1, lg_all, sources_nody, "tcMet_xy", "", ee, mm, em, true, "nody_tcMet_xy_eemm-em", -20, 370); //just for this one, set range manually
  makeSumDifStack(h1, lg_all, sources_nody, "tcMet_xy", "", ee, mm, em, true, "nody_tcMet_xy_eemm-em");
  makeSumDifStack(h1, lg_all, sources_nody, "tcMet_xy_outz", "", ee, mm, em, true, "nody_tcMet_xy_outz_eemm-em");
  makeSumDifStack(h1, lg_all, sources_nody, "tcMet_xy_inz", "", ee, mm, em, true, "nody_tcMet_xy_inz_eemm-em");
  //verbose = false;
  //cout << "setting verbose = false" << endl;
  //h1->setVerbose( false );

  for(int i=0; i<6; i++) { //nsjpbins is hardcoded as 6 here (was 4)
	TString nametcsjp = Form("%s%i", "tcMet_xy_sjp", i);
	TString nmtcousjp = Form("%s%i", "tcMet_xy_outz_sjp", i);
	TString nmtcinsjp = Form("%s%i", "tcMet_xy_inz_sjp", i);

	TString namegnsjp = Form("%s%i", "genMet_xy_sjp", i);
	TString nmgnousjp = Form("%s%i", "genMet_xy_outz_sjp", i);
	TString nmgninsjp = Form("%s%i", "genMet_xy_inz_sjp", i);

	makeStack(h1, lg_all, theSources, nametcsjp, "", all, true);
	makeStack(h1, lg_all, theSources, nmtcousjp, "", all, true);
	makeStack(h1, lg_all, theSources, nmtcinsjp, "", all, true);

	makeStack(h1, lg_all, theSources, nametcsjp, "", em, true);
	makeStack(h1, lg_all, theSources, nmtcousjp, "", em, true);
	makeStack(h1, lg_all, theSources, nmtcinsjp, "", em, true);

	makeSumStack(h1, lg_all, theSources, nametcsjp, "", ee, mm, true);
	makeSumStack(h1, lg_all, theSources, nmtcousjp, "", ee, mm, true);
	makeSumStack(h1, lg_all, theSources, nmtcinsjp, "", ee, mm, true);

	makeSumStack(h1, lg_all, theSources, nametcsjp, "", all, all, true, namegnsjp, -1., "tcMet-" + namegnsjp);
	makeSumStack(h1, lg_all, theSources, nmtcousjp, "", all, all, true, nmgnousjp, -1., "tcMet-" + nmgnousjp);
	makeSumStack(h1, lg_all, theSources, nmtcinsjp, "", all, all, true, nmgninsjp, -1., "tcMet-" + nmgninsjp);

	makeSumDifStack(h1, lg_all, theSources, nametcsjp, "", ee, mm, em, true);
	makeSumDifStack(h1, lg_all, theSources, nmtcousjp, "", ee, mm, em, true);
	makeSumDifStack(h1, lg_all, theSources, nmtcinsjp, "", ee, mm, em, true);

	makeSumDifStack(h1, lg_all, sources_nody, nametcsjp, "", ee, mm, em, true, "nody_" + nametcsjp + "_eemm-em");
	makeSumDifStack(h1, lg_all, sources_nody, nmtcousjp, "", ee, mm, em, true, "nody_" + nmtcousjp + "_eemm-em");
	makeSumDifStack(h1, lg_all, sources_nody, nmtcinsjp, "", ee, mm, em, true, "nody_" + nmtcinsjp + "_eemm-em");

	makeStack(h1, lg_all, theSources, Form("%s%i", "dphi_tcMetl_sjp", i), "", all, true);
	makeStack(h1, lg_all, theSources, Form("%s%i", "dphi_tcMetl_outz_sjp", i), "", all, true);
	makeStack(h1, lg_all, theSources, Form("%s%i", "dphi_tcMetl_inz_sjp", i), "", all, true);

	//for this one, the one without 'den' has weight as weight*tcmet, so don't look at, just use for numerator in divide
	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_mllden_sjp", i), "", all, true);
  }

  //do 2d hist
  TH2F* metxy = h1->get2dHistogram(theSources, "tcMet_xvy", "", all);

  //TString filename = "Results_plot.root";
  //TFile outf(filename,"RECREATE") ;
  TCanvas *c1 = new TCanvas();
  metxy->SetName("tcMet_xvy");
  metxy->Draw("box");
  c1->SaveAs((TString)metxy->GetName()+".png");

  //do sum hists

  //subtraction check
  THStack *tcMet_xy_nody = h1->getSumStack(sources_nody, "tcMet_xy", "", ee, mm);
  TH1F* tcMet_xy_em = h1->getHistogram(theSources, "tcMet_xy", "", em);
  tcMet_xy_nody->Draw("hist");
  tcMet_xy_em->Draw("same");
  c1->SaveAs("tcMet_xy_sf_nody_of_comp.png");
  gPad->SetLogy();
  c1->SaveAs("tcMet_xy_sf_nody_of_comp_log.png");
  gPad->SetLogy(0);

  //tcMet_x v tcMet_y
  //TH1F* tcMet_x_all = h1->getHistogram(theSources, "tcMet_x", "", all);
  //TH1F* tcMet_y_all = h1->getHistogram(theSources, "tcMet_y", "", all);
  ////move difference here
  //TH1F* tcMet_xmy = new TH1F( *tcMet_x_all );
  //tcMet_xmy->Add( tcMet_y_all, -1 );
  //tcMet_xmy->Draw();
  //c1->SaveAs((TString)tcMet_xmy->GetName()+"_diff.png");
  //over_save(tcMet_x_all, tcMet_y_all);
  //over_save(tcMet_x_all, tcMet_y_all, true, true);
  //bools are copylog, copyscale, copydiff
  makeHistOverlay( h1, theSources, theSources, "tcMet_x", "tcMet_y", "", all, all, true, false, true);
  makeHistOverlay( h1, theSources, theSources, "tcMet_x_nar", "tcMet_y_nar", "", all, all, true, false, true);

  //xy_in vs xy_out
  makeHistOverlay( h1, theSources, theSources, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true );
  makeHistOverlay( h1, theSources, theSources, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true );
  sources_t sources_dyee = (1ll << H_DYEE);
  sources_t sources_dymm = (1ll << H_DYMM);
  makeHistOverlay( h1, sources_dyee, sources_dyee, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true, true, (TString)"dyee_");
  makeHistOverlay( h1, sources_dyee, sources_dyee, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, (TString)"dyee_");
  makeHistOverlay( h1, sources_dymm, sources_dymm, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true, true, (TString)"dymm_");
  makeHistOverlay( h1, sources_dymm, sources_dymm, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, (TString)"dymm_");

  //these are wrong--ignore for now
  //makeHistOverlay( h1, sources_dyee, sources_dyee, "tcMet_xy_inz_nar", "tcMet_xy_outz_nar", "", all, all, true, true, true, (TString)"dyee_");
  //makeHistOverlay( h1, sources_dyee, sources_dyee, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, (TString)"dyee_");
  //makeHistOverlay( h1, sources_dymm, sources_dymm, "tcMet_xy_inz_nar", "tcMet_xy_outz_nar", "", all, all, true, true, true, (TString)"dymm_");
  //makeHistOverlay( h1, sources_dymm, sources_dymm, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, (TString)"dymm_");
}


/*
  //do zero sjp bin
  TH1F* tcMet_xy_inz_sjp0 = h1->getHistogram(theSources, "tcMet_xy_inz_sjp0", "", all);
  TH1F* tcMet_xy_ouz_sjp0 = h1->getHistogram(theSources, "tcMet_xy_outz_sjp0", "", all);
  over_save(tcMet_xy_ouz_sjp0, tcMet_xy_inz_sjp0);
  over_save(tcMet_xy_ouz_sjp0, tcMet_xy_inz_sjp0, true, true);
  //scale
  tcMet_xy_ouz_sjp0->Scale( tcMet_xy_inz_sjp0->Integral()/tcMet_xy_ouz_sjp0->Integral() );
  tcMet_xy_ouz_sjp0->SetName("tcMet_xy_outz_sjp0_all_scale");
  over_save(tcMet_xy_ouz_sjp0, tcMet_xy_inz_sjp0);
  over_save(tcMet_xy_ouz_sjp0, tcMet_xy_inz_sjp0, true, true);
  */

/* //don't bother with this one--assume that if you want a stack, you want a legend with it...
void makeStack(HistogramUtilities* h, sources_t theSources, TString title, TString subtitle, TString suffix, bool cpylog=false) {
  THStack *st = h->getStack(theSources, title, subtitle, suffix);

  TCanvas *c = new TCanvas();
  st->Draw("hist");
  c->SaveAs((TString)st->GetName()+".png");
  //delete c;
}
*/
