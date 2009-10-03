
#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
#include "Tools/Utilities.h"
#include "Tools/histtools.cc"

#include "TROOT.h"
//#include "TGraphAsymmErrors.h"

using namespace std;

//const static sources_t theSources_22X =
const static sources_t theSources =
(1ll << H_QCD30)	|
(1ll << H_MU15_SINGLE)	|
(1ll << H_ZEEJET_ALP) 	|
(1ll << H_ZMMJET_ALP)   |
(1ll << H_ZTTJET_ALP)   |
//(1ll << H_WJET_ALP) |
(1ll << H_WEJET_ALP) |
(1ll << H_WMJET_ALP) |
(1ll << H_WTJET_ALP) |
(1ll << H_TTBAR)     
  ;

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
  vSources.push_back(	fH_QCD30()	);
  vSources.push_back(   fH_QCD80()      );
  vSources.push_back(   fH_MU15_SINGLE() );
  vSources.push_back(	fH_TTBAR()	);
  //sources_.push_back(   fH_TTBAR_SINGLE() ); //it's single, but not called that
  h1->setOrder(vSources);
  
  TLegend *lg_all = h1->getLegend(theSources, "lep_pt", "", "all");

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
  //mus
  makeStack(h1, lg_all, theSources, "lep_trckIso",			"", "m", true);
  makeStack(h1, lg_all, theSources, "lep_ecalIso",			"", "m", true);
  makeStack(h1, lg_all, theSources, "lep_relIso",			"", "m", true);

  makeStack(h1, lg_all, theSources, "lep_nlep",				"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_njet20",			"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_njet30",			"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_conversions",		"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_0_pt",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_1_pt",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_pt",				"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_mass",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_tcmet",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_calomet_muon",	"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_njet20",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_njet30",			"", "all", true);

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

