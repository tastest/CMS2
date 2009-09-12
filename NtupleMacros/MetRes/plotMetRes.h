
//#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
//#include "Tools/Utilities.h"
//#include "Tools/histtools.cc"

#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "/home/users/wandrews/macros/comparison.C"

bool verbose = false;

void saveStack(THStack* st, TLegend* leg, bool cpylog, double ymin, double ymax, TString name="", THStack* st2=0) {

  double oldymax = st->GetMaximum();
  //double oldymin = st->GetMinimum();
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
  if( st2 != 0 ) {
	st->SetMinimum( st2->GetMinimum() );
	st2->Draw("same,hist");
  }

  if( verbose )
	cout << "stack " << usename << "  max " << st->GetMaximum() << "  " << oldymax << "  min " << st->GetMinimum() << endl;
  //c->Update();
  //getc(stdin);
  leg->Draw();
  c->SaveAs(usename+".png");

  if(cpylog) {
	gPad->SetLogy();

	//st->SetMaximum( oldymax );
	//if( st2 == 0 ) { //don't reset min if st2 exists--this is wrong
	st->SetMinimum( ymin ); //need to reset regardless
	//cout << "set min to ymin " << ymin << endl;
	//}
  
	st->Draw("hist");
 	leg->Draw();
	c->Update();
	//c->SaveAs((TString)st->GetName()+"_log.png");
	c->SaveAs(usename+"_log.png");
  }
  
}

void makeStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suffix, bool cpylog=false, double ymin=1., double ymax=-999, TString name="") {
  THStack *st = h->getStack(theSources, title, subtitle, suffix);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}

void makeSumStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suff1, TString suff2, bool cpylog=false, TString title2="", double scale=1.0, TString name="", double ymin=1., double ymax=-999) {
  THStack *st = h->getSumStack(theSources, title, subtitle, suff1, suff2, 1, title2, scale);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}

void make2fileStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suff1, bool cpylog=false, double ymin=1., double ymax=-999, TString name="") {
  THStack *st = h->get2fileStack(theSources, title, subtitle, suff1);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}

void makeSumDifStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suff1, TString suff2, TString suff3, bool cpylog=false, TString name="", double ymin=1., double ymax=-999) {
  //second to last argument is scale--leave 1
  THStack *stpos = h->getSumDifStack(theSources, title, subtitle, suff1, suff2, suff3, 1, 1);
  THStack *stneg = h->getSumDifStack(theSources, title, subtitle, suff1, suff2, suff3, 1, -1); //negative
  saveStack( stpos, leg, cpylog, ymin, ymax, name, stneg );
}

void makeHistOverlay( HistogramUtilities* h, sources_t theSources1,  sources_t theSources2, TString title1, TString title2, TString subtitle, TString suff1, TString suff2, bool cpylog=false, bool cpyscale=false, bool cpydiff=false, TString nameprefix="", TString namesuffix="", double ymin=1., double ymax=-999) {
  TH1F* h1 = h->getHistogram(theSources1, title1, "", suff1, 1, nameprefix);
  TH1F* h2 = h->getHistogram(theSources2, title2, "", suff2, 1, nameprefix);
  if( namesuffix != "" ) {
	h1->SetName( (TString)h1->GetName() + namesuffix );
	h2->SetName( (TString)h2->GetName() + namesuffix );
  }
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
	h2->SetTitle( (TString)h2->GetName() + "_diff" ); //for display purposes
	TCanvas *c1 = new TCanvas();
	h1->Draw();
	c1->SaveAs( (TString)h2->GetName() + ".png" );
  }
}

//fit functions

double getHistSigma( TH1F* h ) {
  h->Fit("gaus", "Q");
  TF1 *f = h->GetFunction("gaus"); //Q for quiet
  return f->GetParameter(2); // 0 is a scale constant, 1 is mean, 2 is sigma
  //return {f->GetParameter(2), f->GetParError(2)}; 
	//f->GetChisquare(); //this is useless for now
}

double getTailRatio( TH1F* h ) {
  h->Fit("gaus");
  TF1 *f = h->GetFunction("gaus");
  double sigma = f->GetParameter(2);
  double mean  = f->GetParameter(1);
  //int binlow = h->GetBin( mean - sigma );
  //int binhgh = h->GetBin( mean + sigma );
  int binlow = h->GetXaxis()->FindFixBin( mean - sigma );
  int binhgh = h->GetXaxis()->FindFixBin( mean + sigma );
  cout << "hist " << h->GetName() << "  mean " << mean << "  sigma " << sigma << "  binlow " << binlow << "  binhigh " << binhgh << endl;
  double intsigma = h->Integral( binlow, binhgh );
  double inttaillow = h->Integral( 0, binlow-1 ); //include underflow
  double inttailhgh = h->Integral( binhgh+1, h->GetNbinsX()+1 ); //include overflow
  return intsigma/(inttaillow + inttailhgh);
}


//ratio function : return ratio of integral(low2,hgh2)/integral(low1,hgh1)
//hgh2=0 is for infinity as up range
double getIntegralRatio( TH1F* h, double low1, double hgh1, double low2, double hgh2=0 ) { 
  int binlow1 = h->GetXaxis()->FindFixBin( low1 );
  int binhgh1 = h->GetXaxis()->FindFixBin( hgh1 );
  int binlow2 = h->GetXaxis()->FindFixBin( low2 );
  int binhgh2 = 0;
  double int1 = h->Integral( binlow1, binhgh1 );
  double int2 = 0;
  if( hgh2 == 0 ) //default is infinity as upper limit
	int2 = h->Integral( binlow2, h->GetNbinsX()+1 );
  else {
	binhgh2 = h->GetXaxis()->FindFixBin( hgh2 );
	int2 = h->Integral( binlow2, binhgh2 );
  }

  if( int1 == 0 )
	return 0;
  else
	return int2/int1;
}


//results function

void plotResults() {

  sources_t theSources = sources_all;
  
  //HistogramUtilities h1("Results.root");
  //HistogramUtilities h1("Susy_Results.root");
  HistogramUtilities* h1 = new HistogramUtilities("Results.root");
  HistogramUtilities* h2 = new HistogramUtilities("Results_Nar.root");

  TString all = "all";
  TString ee = "ee"; //these correspond to defintion in DileptonHypType.h
  TString mm = "mm";
  TString em = "em"; 

  //Need this in order to not store clones, which could have same names in 2 files, thus breaking file_->Get()
  TH1::AddDirectory(false); 

  //makeStack(h1, theSources, "tcMet", "", all);
  
  TLegend *lg_all = h1->getLegend(theSources, "sumJetPt_outz", "", all);
  // 
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

  //this value actually comes from looper.h via preprocessor directive...kind of scary, kind of cool
  //int nsjpbins = 6; //nsjpbins is set as 6 here (was 4)
  //nsjpbins = 6;
  double res_ratio1[nsjpbins];
  TH1F* hres_ratio1 = new TH1F( "Ratio_200+_over_50-100", "Ratio_200+_over_50-100", nsjpbins, 0, nsjpbins );
  hres_ratio1->Sumw2();
  double res_ratio2[nsjpbins];
  TH1F* hres_ratio2 = new TH1F( "Ratio_100+_over_0-50", "Ratio_100+_over_0-50", nsjpbins, 0, nsjpbins );
  hres_ratio2->Sumw2();
  double res_ratio3[nsjpbins];
  TH1F* hres_ratio3 = new TH1F( "Ratio_150+_over_50-100", "Ratio_150+_over_50-100", nsjpbins, 0, nsjpbins );
  hres_ratio3->Sumw2();

  for(int i=0; i<nsjpbins; i++) {
	
	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_sjp", i), "", all, true);
	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_outz_sjp", i), "", all, true);
	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_inz_sjp", i), "", all, true);

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

	//integral ratio
	res_ratio1[i] = 0; //initialize here
	res_ratio1[i] = getIntegralRatio( h1->getHistogram(theSources, Form("%s%i", "tcMet_sjp", i), "", all), 50, 100, 200 );
	hres_ratio1->Fill( i, res_ratio1[i] );

	res_ratio2[i] = 0; //initialize here
	res_ratio2[i] = getIntegralRatio( h1->getHistogram(theSources, Form("%s%i", "tcMet_sjp", i), "", all), 0, 50, 100 );
	hres_ratio2->Fill( i, res_ratio2[i] );

	res_ratio3[i] = 0; //initialize here
	res_ratio3[i] = getIntegralRatio( h1->getHistogram(theSources, Form("%s%i", "tcMet_sjp", i), "", all), 50, 100, 150 );
	hres_ratio3->Fill( i, res_ratio3[i] );
  }

  //save integral ratio
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas();
  hres_ratio1->Draw();
  c1->SaveAs((TString)hres_ratio1->GetName()+".png");
  hres_ratio2->Draw();
  c1->SaveAs((TString)hres_ratio2->GetName()+".png");
  hres_ratio3->Draw();
  c1->SaveAs((TString)hres_ratio3->GetName()+".png");
  gStyle->SetOptStat(1110111);
  
  //do 2d hist
  TH2F* metxy = h1->get2dHistogram(theSources, "tcMet_xvy", "", all);
  //TString filename = "Results_plot.root";
  //TFile outf(filename,"RECREATE") ;
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

  //second subtraction check
  THStack *tcMet_xy_subpos = h1->getSumDifStack(theSources, "tcMet_xy", "", ee, mm, em, 1, 1);
  THStack *tcMet_xy_subneg = h1->getSumDifStack(theSources, "tcMet_xy", "", ee, mm, em, 1, -1);
  TH1F* tcMet_xy_dy = h1->getHistogram(sources_dy, "tcMet_xy", "", all);
  tcMet_xy_dy->SetLineColor( 1 ); //0=white, 1=black
  tcMet_xy_dy->SetMarkerColor( 1 );
  tcMet_xy_subpos->SetMinimum( tcMet_xy_subneg->GetMinimum() );
  tcMet_xy_subpos->Draw("hist");
  tcMet_xy_subneg->Draw("same,hist");
  tcMet_xy_dy->Draw("same");
  c1->SaveAs("tcMet_xy_eemm-em_dy_comp.png");
  gPad->SetLogy();
  c1->SaveAs("tcMet_xy_eemm-em_dy_comp_log.png");
  gPad->SetLogy(0);

  //tcMet_x v tcMet_y
  //bools are copylog, copyscale, copydiff
  makeHistOverlay( h1, theSources, theSources, "tcMet_x", "tcMet_y", "", all, all, true, false, true);
  //makeHistOverlay( h1, theSources, theSources, "tcMet_x_nar", "tcMet_y_nar", "", all, all, true, false, true); //replaced with below
  makeHistOverlay( h2, theSources, theSources, "tcMet_x", "tcMet_y", "", all, all, true, false, true,"", "_nar"); //add suffix for h2 
  // * / 
  //xy_in vs xy_out
  makeHistOverlay( h1, theSources, theSources, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true );
  makeHistOverlay( h1, theSources, theSources, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true );
  sources_t sources_dyee = (1ll << H_DYEE);
  sources_t sources_dymm = (1ll << H_DYMM);
  makeHistOverlay( h1, sources_dyee, sources_dyee, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true, true, "dyee_");
  makeHistOverlay( h1, sources_dyee, sources_dyee, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, "dyee_");
  makeHistOverlay( h1, sources_dymm, sources_dymm, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true, true, "dymm_");
  makeHistOverlay( h1, sources_dymm, sources_dymm, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, "dymm_");
  //same in narrow
  makeHistOverlay( h2, sources_dyee, sources_dyee, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true, true, "dyee_", "_nar");
  makeHistOverlay( h2, sources_dyee, sources_dyee, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, "dyee_", "_nar");
  makeHistOverlay( h2, sources_dymm, sources_dymm, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true, true, "dymm_", "_nar");
  makeHistOverlay( h2, sources_dymm, sources_dymm, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, "dymm_", "_nar");

  //these are wrong--ignore for now
  //makeHistOverlay( h1, sources_dyee, sources_dyee, "tcMet_xy_inz_nar", "tcMet_xy_outz_nar", "", all, all, true, true, true, (TString)"dyee_");
  //makeHistOverlay( h1, sources_dyee, sources_dyee, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, (TString)"dyee_");
  //makeHistOverlay( h1, sources_dymm, sources_dymm, "tcMet_xy_inz_nar", "tcMet_xy_outz_nar", "", all, all, true, true, true, (TString)"dymm_");
  //makeHistOverlay( h1, sources_dymm, sources_dymm, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, (TString)"dymm_");


  //fitting -- use h2 b'c finer binning

  //tcMet_xy
  TH1F* tcMet_xy_all = h2->getHistogram(theSources, "tcMet_xy", "", all);
  double sigma_xy_all = getHistSigma( tcMet_xy_all );
  double trato_xy_all = getTailRatio( tcMet_xy_all );
  cout << "name " << tcMet_xy_all->GetName() << "  sigma " << sigma_xy_all << "  ratio core/tail " << trato_xy_all << endl;

  //tcMet_x
  TH1F* tcMet_x_all = h2->getHistogram(theSources, "tcMet_x", "", all);
  double sigma_x_all = getHistSigma( tcMet_x_all );
  double trato_x_all = getTailRatio( tcMet_x_all );
  cout << "name " << tcMet_x_all->GetName() << "  sigma " << sigma_x_all << "  ratio core/tail " << trato_x_all << endl;

  //tcMet_y
  TH1F* tcMet_y_all = h2->getHistogram(theSources, "tcMet_y", "", all);
  double sigma_y_all = getHistSigma( tcMet_y_all );
  double trato_y_all = getTailRatio( tcMet_y_all );
  cout << "name " << tcMet_y_all->GetName() << "  sigma " << sigma_y_all << "  ratio core/tail " << trato_y_all << endl;

  //for dyee+mm
  TH1F* dy_tcMet_xy_all = h2->getHistogram(sources_dy, "tcMet_xy", "", all);
  double dy_sigma_xy_all = getHistSigma( dy_tcMet_xy_all );
  double dy_trato_xy_all = getTailRatio( dy_tcMet_xy_all );
  cout << "name " << dy_tcMet_xy_all->GetName() << "  sigma " << dy_sigma_xy_all << "  ratio core/tail " << dy_trato_xy_all << endl;

  cout << endl << endl;

  //twiki table
  cout << "Fit results for met xy distribution\n";
  cout << "||*DYee+DYmm*|*All Samples*|\n";
  cout << "|*sigma*|" << dy_sigma_xy_all << "|" << sigma_xy_all << "|\n";
  cout << "|*ratio*|" << dy_trato_xy_all << "|" << trato_xy_all << "|\n";
  //add note on twiki on definition of all samples, and ratio
  
}

