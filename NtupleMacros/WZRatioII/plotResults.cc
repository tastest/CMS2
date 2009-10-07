
#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
#include "Tools/Utilities.h"
#include "Tools/histtools.cc"

#include "TROOT.h"
//#include "TGraphAsymmErrors.h"

using namespace std;

//const static sources_t theSources_22X =
const static sources_t theSources =
//(1ll << H_QCD30)	|
(1ll << H_ZEEJET_ALP) 	|
(1ll << H_ZMMJET_ALP)   |
(1ll << H_ZTTJET_ALP)   |
//(1ll << H_WJET_ALP) |
(1ll << H_WEJET_ALP) |
(1ll << H_WMJET_ALP) |
(1ll << H_WTJET_ALP) |
(1ll << H_MU15_SINGLE)	|
(1ll << H_QCDEM)	|
(1ll << H_PHOTONJET)	|
(1ll << H_TTBAR)     
  ;

const static sources_t bkgSources =
(1ll << H_MU15_SINGLE)	|
(1ll << H_QCDEM)	|
//(1ll << H_QCDBCTOE)	|
(1ll << H_PHOTONJET)	|
(1ll << H_TTBAR)       ;

const static sources_t sigsSources =
(1ll << H_WEJET_ALP) |
(1ll << H_WMJET_ALP) |
(1ll << H_WTJET_ALP) ; //does tau count as signal or bkg?

const static sources_t sigdSources =
(1ll << H_ZEEJET_ALP) 	|
(1ll << H_ZMMJET_ALP)   |
(1ll << H_ZTTJET_ALP)   ;


// for 2_2_X
//const static sources_t &theSources = theSources_22X;

TString rootfilename = "wzratio.root";

void saveStack(THStack* st, TLegend* leg, bool cpylog, double ymin, double ymax, TString name="", THStack* st2=0) {

  //double oldymax = st->GetMaximum();
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

  //if( verbose )
  //cout << "stack " << usename << "  max " << st->GetMaximum() << "  " << oldymax << "  min " << st->GetMinimum() << endl;
  //c->Update();
  //getc(stdin);
  leg->Draw();
  c->SaveAs(usename+".png");
  c->SaveAs(rootfilename);

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
	c->SaveAs(rootfilename);
  }
  
}

void makeStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suffix, bool cpylog=false, double ymin=1., double ymax=-999, TString name="") {
  THStack *st = h->getStack(theSources, title, subtitle, suffix);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}



void plotResults() {

  HistogramUtilities* h1 = new HistogramUtilities("Results.root");  //normalization is in weight

  vector<DataSource> vSources;
  vSources.push_back(   fH_WEJET_ALP()   );	
  vSources.push_back(   fH_WMJET_ALP()   );	
  vSources.push_back(   fH_WTJET_ALP()   );	
  vSources.push_back(	fH_ZEEJET_ALP()	);
  vSources.push_back(   fH_ZMMJET_ALP() );
  vSources.push_back(   fH_ZTTJET_ALP() );
  //vSources.push_back(	fH_QCD30()	);
  //vSources.push_back(   fH_QCD80()      );
  vSources.push_back(	fH_QCDEM()	);
  vSources.push_back(	fH_PHOTONJET()	);
  vSources.push_back(   fH_MU15_SINGLE() );
  vSources.push_back(	fH_TTBAR()	);
  //sources_.push_back(   fH_TTBAR_SINGLE() ); //it's single, but not called that
  h1->setOrder(vSources);
  
  TLegend *lg_all = h1->getLegend(theSources, "lep_pt", "", "all");
  /*
  makeStack(h1, lg_all, theSources, "Highlep_pt", 			"", "all", true);
  makeStack(h1, lg_all, theSources, "Highlep_Met", 			"", "all", true);
  makeStack(h1, lg_all, theSources, "Highlep_RelIso", 		"", "all", true);
  makeStack(h1, lg_all, theSources, "Highlep_RelIsoPtLg20", "", "all", true);
  makeStack(h1, lg_all, theSources, "Lowlep_pt",			"", "all", true);
  makeStack(h1, lg_all, theSources, "Lowlep_Met",			"", "all", true);
  makeStack(h1, lg_all, theSources, "Lowlep_RelIso",		"", "all", true);
  makeStack(h1, lg_all, theSources, "Lowlep_RelIsoPtLg20",	"", "all", true);
  makeStack(h1, lg_all, theSources, "Lowlep_NLepGt10Lt20",	"", "all", true);
  makeStack(h1, lg_all, theSources, "Lowlep_NLepGt20",		"", "all", true);
  */
  /*
  makeStack(h1, lg_all, theSources, "lep_pt",				"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_transmass",		"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_tcmet",			"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_calomet_muon",		"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_met_dphi",			"", "all", true);
  //isos
  makeStack(h1, lg_all, theSources, "lep_trckIso",			"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_ecalIso",			"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_relIso",			"", "all", true);
  //els
  makeStack(h1, lg_all, theSources, "lep_trckIso",			"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_ecalIso",			"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_relIso",			"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_njet20",			"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_njet30",			"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_conversions",		"", "e", true);
  //mus
  makeStack(h1, lg_all, theSources, "lep_trckIso",			"", "m", true);
  makeStack(h1, lg_all, theSources, "lep_ecalIso",			"", "m", true);
  makeStack(h1, lg_all, theSources, "lep_relIso",			"", "m", true);
  makeStack(h1, lg_all, theSources, "lep_njet20",			"", "m", true);
  makeStack(h1, lg_all, theSources, "lep_njet30",			"", "m", true);

  makeStack(h1, lg_all, theSources, "lep_nlep",				"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_njet20",			"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_njet30",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_0_pt",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_1_pt",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_pt",				"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_mass",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_tcmet",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_calomet_muon",	"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_njet20",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_njet30",			"", "all", true);
*/

  //Tests for ABCD
  //get this manually (not histogramUtilities way)
  TFile* file = new TFile("Results.root" , "READ");
  //TH2Fs
  //TH2F* h = (TH2F*)(file->Get(sources_[i].getName() + histNameSuffix)->Clone());
  TH2F* h = (TH2F*)(file->Get("wejetsAlpgen_lep_tcMet_relIso_all")->Clone());
  //printTH2F( h );
  cout << endl << endl;
  compareTH2F( h, h );

  TH2F* hdytt_a = (TH2F*)(file->Get("dyttAlpgen_lep_tcMet_relIso_all")->Clone());  
  TH2F* hdytt_m = (TH2F*)(file->Get("dyttAlpgen_lep_tcMet_relIso_m")->Clone());  
  //cout << endl << endl;
  //compareTH2F( hdytt_a, hdytt_m );

  TH2F* hdytt_e = (TH2F*)(file->Get("dyttAlpgen_lep_tcMet_relIso_e")->Clone());  
  //cout << endl << endl;
  hdytt_m->Add(hdytt_e);
  compareTH2F( hdytt_a, hdytt_m );

  //cout << integrateTH2F( hdytt_a, 0, 50, 0.89, 1.0 ) << endl;
  //I expect this to be sum of:
  /*
bin 55,100  first 0.0754978   second 0
bin 60,98   first 0.00046883   second 0
bin 50,96   first 0.0008938   second 0
bin 21,95   first 0.0013502   second 0
bin 64,95   first 0.00787893   second 0
bin 66,94   first 0.00046883   second 0
bin 14,92   first 0.0000376945   3.76945e-05   second 0
bin 77,92   first 0.000284761   second 0
bin 81,92   first 0.00787893   second 0
bin 60,90   first 0.000927268   second 0

0.095687043
  */

  TCanvas* c1 = new TCanvas();

  //ABCD
  //tcmet vs reliso
  //e+m
  TH2F* tcmet_reliso_sig = h1->get2dHistogram(sigsSources, "lep_tcMet_relIso", "", "all", 1, "_sig"); //all or e/mu?
  //tcmet_reliso_sig->SetName("tcMet_relIso_signal");
  tcmet_reliso_sig->Draw("box");
  c1->SaveAs((TString)tcmet_reliso_sig->GetName()+".png");

  TH2F* tcmet_reliso_bkg = h1->get2dHistogram(bkgSources, "lep_tcMet_relIso", "", "all", 1, "_bkg"); //all or e/mu?
  //tcmet_reliso_bkg->SetName("tcMet_relIso_background");
  tcmet_reliso_bkg->Draw("box");
  c1->SaveAs((TString)tcmet_reliso_bkg->GetName()+".png");

  //guess to start with: a : x1=20, x2=51 (oflw), y3=0, y4=0.1
  double x1 = 20.;
  double x2 = 51.;
  double x3 = 0.;
  double x4 = 15.;
  double y1 = 0.2;
  double y2 = 1.;
  double y3 = 0.;
  double y4 = 0.1;
  cout << abcd(tcmet_reliso_bkg, x1, x2, x3, x4, y1, y2, y3, y4) << endl;

  //e
  TH2F* tcmet_reliso_sig_e = h1->get2dHistogram(sigsSources, "lep_tcMet_relIso", "", "e", 1, "_sig"); 
  //tcmet_reliso_sig_e->SetName("tcMet_relIso_signal_e");
  tcmet_reliso_sig_e->Draw("box");
  c1->SaveAs((TString)tcmet_reliso_sig_e->GetName()+".png");

  TH2F* tcmet_reliso_bkg_e = h1->get2dHistogram(bkgSources, "lep_tcMet_relIso", "", "e", 1, "_bkg");
  //tcmet_reliso_bkg->SetName("tcMet_relIso_background_e");
  tcmet_reliso_bkg_e->Draw("box");
  c1->SaveAs((TString)tcmet_reliso_bkg_e->GetName()+".png");

  cout << abcd(tcmet_reliso_bkg_e, x1, x2, x3, x4, y1, y2, y3, y4) << endl;

  //m
  TH2F* tcmet_reliso_sig_m = h1->get2dHistogram(sigsSources, "lep_tcMet_relIso", "", "m", 1, "_sig"); 
  //tcmet_reliso_sig_m->SetName("tcMet_relIso_signal_m");
  tcmet_reliso_sig_m->Draw("box");
  c1->SaveAs((TString)tcmet_reliso_sig_m->GetName()+".png");

  TH2F* tcmet_reliso_bkg_m = h1->get2dHistogram(bkgSources, "lep_tcMet_relIso", "", "m", 1, "_bkg");
  //tcmet_reliso_bkg->SetName("tcMet_relIso_background_m");
  tcmet_reliso_bkg_m->Draw("box");
  c1->SaveAs((TString)tcmet_reliso_bkg_m->GetName()+".png");

  cout << abcd(tcmet_reliso_bkg_m, x1, x2, x3, x4, y1, y2, y3, y4) << endl;
}


//currently, returns a = b*c/d, so call it such that this is what you want
double abcd(TH2F* h, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
  //double a = integrateTH2F(h, x1, x2, y3, y4); //this is what you're predicting, dumb-ass
  double b = integrateTH2F(h, x3, x4, y3, y4);
  double c = integrateTH2F(h, x1, x2, y1, y2);
  double d = integrateTH2F(h, x3, x4, y1, y2);
  return b*c/d;
}


double integrateTH2F(TH2F* h, double xlow, double xhgh, double ylow, double yhgh) {
  int binlowx = h->GetXaxis()->FindFixBin( xlow );
  int binhghx = h->GetXaxis()->FindFixBin( xhgh );
  int binlowy = h->GetYaxis()->FindFixBin( ylow );
  int binhghy = h->GetYaxis()->FindFixBin( yhgh );
  if( xhgh == 0. ) binhghx = h->GetNbinsX(); //default is max bin
  if( yhgh == 0. ) binhghy = h->GetNbinsY(); //default is max bin
  double integral = 0.;

  for( int j=binlowy; j < binhghy; j++ ) {
	for( int i=binlowx; i < binhghx; i++ ) { 
	  integral += h->GetBinContent(i,j);
	}
  }
  return integral;
}


//prints only when bin differences are found
//make sure h1, h2 have same NbinsX,Y
void compareTH2F(TH2F* h1, TH2F* h2) {
  //for( int i=0; i < h1->GetNbinsY(); i++ ) {
  for( int i=h1->GetNbinsY(); i >= 0; i-- ) { //start at top and go down
  //for( int i=h1->GetNbinsY(); i >= h1->GetNbinsY()-10; i-- ) { //start at top and go down, first 10 only
	for( int j=0; j < h1->GetNbinsX(); j++ ) {
	  if( h1->GetBinContent(j,i) - h2->GetBinContent(j,i) > 1.0e-5 )
		cout << "bin " << j << "," << i << "   first " << h1->GetBinContent(j,i) << "   second " << h2->GetBinContent(j,i) << endl;
	}
  }
}

void printTH2F(TH2F* h) {

  //for( int i=0; i < h->GetNbinsY(); i++ ) {
  for( int i=0; i < 10; i++ ) { //10 staring at bottom
  //for( int i=h->GetNbinsY(); i >= 0; i-- ) { //start at top and go down
  //for( int i=h->GetNbinsY(); i >= h->GetNbinsY()-10; i-- ) { //just 10 rows
	//for( int j=0; j < h->GetNbinsX(); j++ ) {
	for( int j=0; j < 30; j++ ) { //first 30 in x
	  printf("%4.1f ", h->GetBinContent(j,i) );
	}
	cout << endl;
  }

}



void plotResultsLep(TString hyp)
{

  gROOT->ProcessLine(".L ~/tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle()");

  // sources/ordering for Z stack plots
  std::vector<DataSource> zSources;

  if (hyp == "e") {
	zSources.push_back(     fH_WJET_ALP()   );
	zSources.push_back(     fH_QCD30()	);
	zSources.push_back(     fH_ZEEJET_ALP()   );
	zSources.push_back(     fH_ZTTJET_ALP()   );
	zSources.push_back(     fH_ZMMJET_ALP()   );
	
  } else if (hyp == "m") {
	zSources.push_back(     fH_WJET_ALP()   );
	zSources.push_back(     fH_MU15_SINGLE()      );
	zSources.push_back(     fH_ZMMJET_ALP()   );
	zSources.push_back(     fH_ZTTJET_ALP()   );
	zSources.push_back(     fH_ZEEJET_ALP()   );
  }
  
  // luminorm for 1pb-1
  HistogramUtilities h1("Results.root", 0.001);
  h1.setOrder(zSources);
  
  TLegend *lg_all = h1.getLegend(theSources, "lep_met", "", hyp);
  THStack *st_lep_met = h1.getStack(theSources, "lep_met", "", hyp);
  THStack *st_lep_met_dphi = h1.getStack(theSources, "lep_met_dphi", "", hyp);
  
  TCanvas *c1 = new TCanvas();
  c1->cd();
  st_lep_met->Draw();
  lg_all->Draw();
  Utilities::saveCanvas(c1, "results/lep_met_" + hyp);
  
  TCanvas *c2 = new TCanvas();
  c2->cd();
  st_lep_met_dphi->Draw();
  lg_all->Draw();
  Utilities::saveCanvas(c2, "results/lep_met_dphi_" + hyp);
  
}


void plotResultsDilep(TString hyp)
{

  gROOT->ProcessLine(".L ~/tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle()");

  // sources/ordering for Z stack plots
  std::vector<DataSource> zSources;

  if (hyp == "ee") {
	zSources.push_back(     fH_ZEEJET_ALP()   );
	zSources.push_back(     fH_QCD30()      );
	zSources.push_back(     fH_WJET_ALP()   );
	zSources.push_back(     fH_ZTTJET_ALP()   );
	zSources.push_back(     fH_ZMMJET_ALP()   );

  } else if (hyp == "mm") {
	zSources.push_back(     fH_ZMMJET_ALP()   );
	zSources.push_back(     fH_MU15_SINGLE()      );
	zSources.push_back(     fH_WJET_ALP()   );
	zSources.push_back(     fH_ZTTJET_ALP()   );
	zSources.push_back(     fH_ZEEJET_ALP()   );
  }

  // luminorm for 1pb-1
  //HistogramUtilities h1("Results.root", 0.001);
  HistogramUtilities h1("Results.root");
  h1.setOrder(zSources);

  TLegend *lg_all = h1.getLegend(theSources, "dilep_mass", "", hyp);
  THStack *st_dilep_mass = h1.getStack(theSources, "dilep_mass", "", hyp, 4);
  THStack *st_dilep_met = h1.getStack(theSources, "dilep_met", "", hyp);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  st_dilep_mass->Draw();
  lg_all->Draw();
  Utilities::saveCanvas(c1, "results/dilep_mass_" + hyp);	

  TCanvas *c2 = new TCanvas();
  c2->cd();
  st_dilep_met->Draw();
  lg_all->Draw();
  Utilities::saveCanvas(c2, "results/dilep_met_" + hyp);


}

