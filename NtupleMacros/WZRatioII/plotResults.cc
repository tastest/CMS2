
//#include "TH2F.h"
#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
//#include "Tools/DataSource.h"
#include "Tools/Utilities.h"
#include "Tools/histtools.cc"

#include "TROOT.h"
//#include "TGraphAsymmErrors.h"

using namespace std;

TString rootfilename = "wzratio.root";

vector<DataSource> getSources() {
  vector<DataSource> vSources;
  vSources.push_back(   fH_WEJET_ALP()   );	
  vSources.push_back(   fH_WMJET_ALP()   );	
  vSources.push_back(   fH_WTJET_ALP()   );	
  vSources.push_back(	fH_ZEEJET_ALP()	);
  vSources.push_back(   fH_ZMMJET_ALP() );
  vSources.push_back(   fH_ZTTJET_ALP() );
  vSources.push_back(	fH_QCD30()	);
  vSources.push_back(   fH_QCD80()      );
  vSources.push_back(	fH_QCDEM()	);
  vSources.push_back(	fH_PHOTONJET()	);
  vSources.push_back(	fH_QCDBCTOE()	);
  vSources.push_back(   fH_MU15_SINGLE() );
  vSources.push_back(	fH_TTBAR()	);
  //sources_.push_back(   fH_TTBAR_SINGLE() ); //it's single, but not called that
  return vSources;
}

void plotResults() {

  //HistogramUtilities* h1 = new HistogramUtilities("Results.root");  //normalization is in weight
  HistogramUtilities* h1 = new HistogramUtilities("ABCDResults.root");  //normalization is in weight
  vector<DataSource> vSources = getSources();
  h1->setOrder(vSources);
  
  TLegend *lg_all = h1->getLegend(theSources, "lep_pt", "", "all");

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
  makeStack(h1, lg_all, theSources, "lep_nlep_nod0iso",		"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_nlep_nometiso",	"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_conversions",		"", "e", true);
  //mus
  makeStack(h1, lg_all, theSources, "lep_trckIso",			"", "m", true);
  makeStack(h1, lg_all, theSources, "lep_ecalIso",			"", "m", true);
  //makeStack(h1, lg_all, theSources, "lep_hcalIso",			"", "m", true);
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

  makeStack(h1, lg_all, theSources, "dilep_mass_ll20",		"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_pt_mgen",		"", "all", true);
  //makeStack(h1, lg_all, theSources, "dilep_0_genpt", "", "all", true);
  //makeStack(h1, lg_all, theSources, "dilep_1_genpt", "", "all", true);
  makeStack(h1, lg_all, theSources, "lep_genpt", "", "all", true);

  makeStack(h1, lg_all, theSources, "dilep_mass_ll20",		"", "ee", true);
  makeStack(h1, lg_all, theSources, "dilep_pt_mgen",		"", "ee", true);
  //makeStack(h1, lg_all, theSources, "dilep_0_genpt", "", "ee", true);
  //makeStack(h1, lg_all, theSources, "dilep_1_genpt", "", "ee", true);
  //makeStack(h1, lg_all, theSources, "dilep_0_genpt_mch", "", "ee", true);
  //makeStack(h1, lg_all, theSources, "dilep_1_genpt_mch", "", "ee", true);
  makeStack(h1, lg_all, theSources, "dilep_genpt", "", "ee", true);
  makeStack(h1, lg_all, theSources, "lep_genpt", "", "e", true);
  makeStack(h1, lg_all, theSources, "lep_genpt_mch", "", "e", true);

}
//end plotResults

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

template <class TH> void printHistStats(TH* h, double xlow, double xhgh) {
  int hghbin = 0;
  if( xhgh == -1. ) hghbin = h->GetXaxis()->GetNbins() + 1; //include overflow
  else hghbin = h->GetXaxis()->FindFixBin( xhgh );
  cout << h->GetName() << " \t\t " << h->GetMean() << "  " << h->Integral()
	   << "  " << h->Integral( h->GetXaxis()->FindFixBin( xlow ), hghbin )/h->Integral() << endl;
}

TH2F* getTH2F( HistogramUtilities* h, sources_t theSources, TString var, TString var2, TString hyptyp, Int_t rebin, TString suffix, TString opt) {
  TCanvas* c = new TCanvas();
  TH2F* hnew = new TH2F(*(h->get2dHistogram(theSources, var, var2, hyptyp, rebin, suffix)));
  hnew->Draw(opt);
  c->SaveAs((TString)hnew->GetName()+".png");
  return hnew;
}

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

void makeStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suffix, bool cpylog, double ymin, double ymax, TString name) {
  THStack *st = h->getStack(theSources, title, subtitle, suffix);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}


