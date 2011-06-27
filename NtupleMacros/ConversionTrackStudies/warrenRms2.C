
#include "/code/osgcode/cmssoft/cms/slc5_ia32_gcc434/cms/cmssw/CMSSW_3_4_0/src/Validation/RecoParticleFlow/interface/TH2Analyzer.h"
#include "TH2D.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TF1.h"
#include <iostream>

using namespace std;

void doRms(string datatree, string mctree, bool is2tev=false) {

  TChain* mcchain = new TChain("T1");
  mcchain->Add(mctree.c_str());

  TChain* dachain = new TChain("T1");
  dachain->Add(datatree.c_str());

  //const bool doall=true;
  const bool noMC=false;
  const bool doc3c4=true;

  const unsigned int nEventsMC = mcchain->GetEntries();
  const unsigned int nEventsRD = dachain->GetEntries();
  
  if( !nEventsMC ) {
	cout << "n events in mc is zero" << endl;
	exit(1);
  }
  else {
	cout<<"Number of events in the MC: "<<nEventsMC<<endl;
  } 
	
  if( !nEventsRD ) {
	cout << "n events in mc is zero" << endl;
	exit(1);
  }
  else {
	cout<<"Number of events in the data: "<<nEventsRD<<endl;
  } 


  std::string title;
  if( is2tev )
	title = "CMS Preliminary 2009, 2.36 TeV data";
  else
	title = "CMS Preliminary 2009, 900 GeV data";

  const int nbins=12;
  //const double sumEtMax=55.;
  const double sumEtMax=40.;
  const int nbins_MEX=200;
	
  //// tc MC
  TH2D* MC = new TH2D("MC", title.c_str(), nbins, 0., sumEtMax, nbins_MEX, -50., 50.); 
  MC->GetXaxis()->SetTitle("tcSumEt [GeV]");
  MC->GetYaxis()->SetTitle("tcMEX [GeV]");
  mcchain->Draw("tcmet*cos(tcmetphi):tcsumet>>MC" );

  //test canvas
  //TCanvas *ctest = new TCanvas("ctest", "ctest",14,31,700,500);
  //MC->Draw();
	
  TH2D* MCy = new TH2D("MCy", title.c_str(), nbins, 0., sumEtMax, nbins_MEX, -50., 50.); 
  mcchain->Draw("tcmet*sin(tcmetphi):tcsumet>>MCy" );
  MC->Add(MCy,1.0);
	
  TH2Analyzer TH2Ana(MC,1);
  //TH2Analyzer TH2Ana( &MC, 1, 100, 10, false);
  //TH2Analyzer TH2Ana( &MC, 1, nbins, 1, true);
  TH1D* ha = TH2Ana.SigmaGauss();
  ha->SetLineColor(2);
  ha->SetLineWidth(3);

  //// tc data	
  TH2D* RD = new TH2D("RD", title.c_str(), nbins, 0., sumEtMax, nbins_MEX, -50., 50.); 
  RD->GetXaxis()->SetTitle("tcSumEt [GeV]");
  RD->GetYaxis()->SetTitle("tcMEX [GeV]");
  dachain->Draw("tcmet*cos(tcmetphi):tcsumet>>RD" ); 
	
  TH2D* RDy = new TH2D("RDy", title.c_str(), nbins, 0., sumEtMax, nbins_MEX, -50., 50.); 
  dachain->Draw("tcmet*sin(tcmetphi):tcsumet>>RDy" ); 
  RD->Add(RDy,1.0);
  
  //TH2Analyzer TH2AnaRD( &RD, 1, nbins, 1, true);
  TH2Analyzer TH2AnaRD( RD, 1);
  TH1D* haRD = TH2AnaRD.SigmaGauss();
  haRD->SetMarkerStyle(21);
  haRD->SetTitle(title.c_str());
  haRD->GetXaxis()->SetTitle("SumEt [GeV]");
  haRD->GetYaxis()->SetTitle("#sigma(MEX,MEY) [GeV]");
  haRD->GetYaxis()->SetRangeUser(0.0,6.0); 
	
	
  //// calo MC
  TH2D* MC_calo = new TH2D("MC_calo", title.c_str(), nbins, 0., sumEtMax, nbins_MEX, -50., 50.); 
  MC_calo->GetXaxis()->SetTitle("caloSumEt [GeV]");
  MC_calo->GetYaxis()->SetTitle("caloMEX [GeV]");
  mcchain->Draw("met*cos(metphi):sumet>>MC_calo" ); 
	
  TH2D* MC_caloy = new TH2D("MC_caloy", title.c_str(), nbins, 0., sumEtMax, nbins_MEX, -50., 50.); 
  mcchain->Draw("met*sin(metphi):sumet>>MC_caloy" ); 
  MC_calo->Add(MC_caloy,1.0);
	
  //TH2Analyzer TH2Ana_calo( &MC_calo, 1, nbins, 1, true);
  TH2Analyzer TH2Ana_calo( MC_calo, 1);
  TH1D* ha_calo = TH2Ana_calo.SigmaGauss();
  ha_calo->SetLineColor(3);
  ha_calo->SetLineWidth(3);

  //// calo data
  TH2D* RD_calo = new TH2D("RD_calo", title.c_str(), nbins, 0., sumEtMax, nbins_MEX, -50., 50.); 
  RD_calo->GetXaxis()->SetTitle("caloSumEt [GeV]");
  RD_calo->GetYaxis()->SetTitle("caloMEX [GeV]");
  dachain->Draw("met*cos(metphi):sumet>>RD_calo" );
	
  TH2D* RD_caloy = new TH2D("RD_caloy", title.c_str(), nbins, 0., sumEtMax, nbins_MEX, -50., 50.); 
  dachain->Draw("met*sin(metphi):sumet>>RD_caloy" );
  RD_calo->Add(RD_caloy,1.0);
	
  //TH2Analyzer TH2AnaRD_calo( &RD_calo, 1, nbins, 1, true);
  TH2Analyzer TH2AnaRD_calo( RD_calo, 1);
  TH1D* haRD_calo = TH2AnaRD_calo.SigmaGauss();
  haRD_calo->SetMarkerStyle(21);
  haRD_calo->SetMarkerColor(4);
  haRD_calo->SetLineColor(4);
	
  TCanvas *c1 = new TCanvas("c1", "c1",14,31,700,500);
  c1->SetFillColor(0);
  gStyle->SetOptStat(0);
	
  // draw all 4 above
  haRD->Draw("E1");
  haRD->GetXaxis()->SetTitleSize(0.05);
  if (!noMC) ha->Draw("E1same");
  if (!noMC) ha_calo->Draw("E1same");
  haRD_calo->Draw("E1same");
			
			
  TLegend *lc1 = new TLegend(0.5502874,0.7669492,0.716954,0.8771186,NULL,"brNDC");
  lc1->SetFillColor(10);
  lc1->AddEntry(haRD,"data, TC","p");
  if (!noMC) lc1->AddEntry(ha,"MC, TC","l");
  lc1->AddEntry(haRD_calo,"data, Calo","p");
  if (!noMC) lc1->AddEntry(ha_calo,"MC, Calo","l");
  lc1->Draw();
  //c1->Update();
  
  std::string eps_c1 = "sigmaMEX";
  if( is2tev )
	eps_c1 += "_2tev";
  //gPad->SaveAs( (eps_c1 + ".eps").c_str() );
  //c1->SaveAs( (eps_c1 + ".eps").c_str() );
					
  ///////////////////////
  // above is sigmagauss, below is rms. Forget rms.
  ///////////////////////
  /*					
  TCanvas *c2 = new TCanvas("c2", "c2",14,31,700,500);
  c2->SetFillColor(0);
  gStyle->SetOptStat(0);
					
  //const double sumEtMax_c2=55.;
  const double sumEtMax_c2=40.;
  const int nbins_MEX_c2=200;
  const int nbins_c2=12;
					
  //// tc
  TH2D* MC_c2 = new TH2D("MC_c2", title.c_str(), nbins_c2, 0., sumEtMax_c2, nbins_MEX_c2, -50., 50.); 
  MC_c2->GetXaxis()->SetTitle("tcSumEt [GeV]");
  MC_c2->GetYaxis()->SetTitle("tcMEX [GeV]");
  mcchain->Draw("tcmet*cos(tcmetphi):tcsumet>>MC_c2" );
					
  TH2D* MCy_c2 = new TH2D("MCy_c2", title.c_str(), nbins_c2, 0., sumEtMax_c2, nbins_MEX_c2, -50., 50.); 
  mcchain->Draw("tcmet*sin(tcmetphi):tcsumet>>MCy_c2" );
  MC_c2->Add(MCy_c2,1.0);
					
  TH2Analyzer TH2Ana_c2(MC_c2,1);
  //TH2Analyzer TH2Ana( &MC, 1, 100, 10, false);
  //TH2Analyzer TH2Ana( &MC, 1, nbins, 1, true);
  TH1D* ha_c2=TH2Ana_c2.RMS();
  ha_c2->SetLineColor(2);
  ha_c2->SetLineWidth(3);
					
  TH2D* RD_c2 = new TH2D("RD_c2", title.c_str(), nbins_c2, 0., sumEtMax_c2, nbins_MEX_c2, -50., 50.); 
  RD_c2->GetXaxis()->SetTitle("tcSumEt [GeV]");
  RD_c2->GetYaxis()->SetTitle("tcMEX [GeV]");
  dachain->Draw("tcmet*cos(tcmetphi):tcsumet>>RD_c2" ); 
					
  TH2D* RDy_c2 = new TH2D("RDy_c2", title.c_str(), nbins_c2, 0., sumEtMax_c2, nbins_MEX_c2, -50., 50.); 
  dachain->Draw("tcmet*sin(tcmetphi):tcsumet>>RDy_c2" ); 
  RD_c2->Add(RDy_c2,1.0);
  
  //TH2Analyzer TH2AnaRD( &RD, 1, nbins, 1, true);
  TH2Analyzer TH2AnaRD_c2( RD_c2, 1);
  TH1D* haRD_c2=TH2AnaRD_c2.RMS();
  haRD_c2->SetMarkerStyle(21);
  haRD_c2->SetTitle(title.c_str());
  haRD_c2->GetXaxis()->SetTitle("SumEt [GeV]");
  haRD_c2->GetYaxis()->SetTitle("RMS(MEX,MEY) [GeV]");
  haRD_c2->GetYaxis()->SetRangeUser(0.0,6.0); 
  
  //// calo
  TH2D* MC_calo_c2 = new TH2D("MC_calo_c2", title.c_str(), nbins_c2, 0., sumEtMax_c2, nbins_MEX_c2, -50., 50.); 
  MC_calo_c2->GetXaxis()->SetTitle("caloSumEt [GeV]");
  MC_calo_c2->GetYaxis()->SetTitle("caloMEX [GeV]");
  mcchain->Draw("CaloMet.px():CaloMet.sumEt()>>MC_calo_c2" ); 
  
  TH2D* MC_caloy_c2 = new TH2D("MC_caloy_c2", title.c_str(), nbins_c2, 0., sumEtMax_c2, nbins_MEX_c2, -50., 50.); 
  mcchain->Draw("CaloMet.py():CaloMet.sumEt()>>MC_caloy_c2" ); 
  MC_calo_c2->Add(MC_caloy_c2,1.0);
  
  //TH2Analyzer TH2Ana_calo( &MC_calo, 1, nbins, 1, true);
  TH2Analyzer TH2Ana_calo_c2( MC_calo_c2, 1);
  TH1D* ha_calo_c2=TH2Ana_calo_c2.RMS();
  ha_calo_c2->SetLineColor(3);
  ha_calo_c2->SetLineWidth(3);
					
  TH2D* RD_calo_c2 = new TH2D("RD_calo_c2", title.c_str(), nbins_c2, 0., sumEtMax_c2, nbins_MEX_c2, -50., 50.); 
  RD_calo_c2->GetXaxis()->SetTitle("caloSumEt [GeV]");
  RD_calo_c2->GetYaxis()->SetTitle("caloMEX [GeV]");
  dachain->Draw("CaloMet.px():CaloMet.sumEt()>>RD_calo_c2" );
					
  TH2D* RD_caloy_c2 = new TH2D("RD_caloy_c2", title.c_str(), nbins_c2, 0., sumEtMax_c2, nbins_MEX_c2, -50., 50.); 
  dachain->Draw("CaloMet.py():CaloMet.sumEt()>>RD_caloy_c2" );
  RD_calo_c2->Add(RD_caloy_c2,1.0);
					
  //TH2Analyzer TH2AnaRD_calo( &RD_calo, 1, nbins, 1, true);
  TH2Analyzer TH2AnaRD_calo_c2( RD_calo_c2, 1);
  TH1D* haRD_calo_c2=TH2AnaRD_calo_c2.RMS();
  haRD_calo_c2->SetMarkerStyle(21);
  haRD_calo_c2->SetMarkerColor(4);
  haRD_calo_c2->SetLineColor(4);
  
  // draw !
  haRD_c2->Draw("E1");
  if (!noMC) ha_c2->Draw("E1same");
  if (!noMC) ha_calo_c2->Draw("E1same");
  haRD_calo_c2->Draw("E1same");
  
  
  TLegend *lc2 = new TLegend(0.5502874,0.7669492,0.716954,0.8771186,NULL,"brNDC");
  lc2->SetFillColor(10);
  lc2->AddEntry(haRD_c2,"data, TC","p");
  if (!noMC) lc2->AddEntry(ha_c2,"MC, TC","l");
  lc2->AddEntry(haRD_calo_c2,"data, Calo","p");
  if (!noMC) lc2->AddEntry(ha_calo_c2,"MC, Calo","l");
  lc2->Draw();
  c2->Update();
									
  const std::string eps_c2 = "rmsMEX";
  gPad->SaveAs( (eps_c2 + ".eps").c_str() );

  //commented out the rms half
  */									
									
  if (doc3c4) {
	///////////////////////
	// c3
	///////////////////////
										
	//TCanvas *c3 = new TCanvas("c3", "c3",14,31,700,500);
	TCanvas *c3 = new TCanvas("c3", "c3");
	c3->SetFillColor(0);
	gStyle->SetOptStat(0);

	//old is 12 bins to 40, or 30 to 100. new is 21 to 70.
	//const double sumEtMax_c3=110.;
	const double sumEtMax_c3=70.;
	const int nbins_MEX_c3=200;
	const int nbins_c3=21;
	const int tccolor = 896;
				
	//// tc MC
	TH2D* MC_c3 = new TH2D("MC_c3", title.c_str(), nbins_c3, 0., sumEtMax_c3, nbins_MEX_c3, -50., 50.); 
	MC_c3->GetXaxis()->SetTitle("tcSumEt [GeV]");
	MC_c3->GetYaxis()->SetTitle("tcMEX [GeV]");
	mcchain->Draw("tcmet*cos(tcmetphi):tcsumet>>MC_c3" );
										
	TH2D* MCy_c3 = new TH2D("MCy_c3", title.c_str(), nbins_c3, 0., sumEtMax_c3, nbins_MEX_c3, -50., 50.); 
	mcchain->Draw("tcmet*sin(tcmetphi):tcsumet>>MCy_c3" );
	MC_c3->Add(MCy_c3,1.0);
										
	TH2Analyzer TH2Ana_c3(MC_c3,1);
	//TH2Analyzer TH2Ana( &MC, 1, 100, 10, false);
	//TH2Analyzer TH2Ana( &MC, 1, nbins, 1, true);
	TH1D* ha_c3 = TH2Ana_c3.SigmaGauss();
	ha_c3->SetLineColor(tccolor);
	ha_c3->SetFillColor(0);
	ha_c3->SetLineWidth(3);
	ha_c3->SetMarkerStyle(0);
	ha_c3->SetMarkerColor(tccolor);
	ha_c3->SetTitle("");
	
	//// tc data
	TH2D* RD_c3 = new TH2D("RD_c3", title.c_str(), nbins_c3, 0., sumEtMax_c3, nbins_MEX_c3, -50., 50.); 
	RD_c3->GetXaxis()->SetTitle("tcSumEt [GeV]");
	RD_c3->GetYaxis()->SetTitle("tcMEX [GeV]");
	dachain->Draw("tcmet*cos(tcmetphi):tcsumet>>RD_c3" ); 
										
	TH2D* RDy_c3 = new TH2D("RDy_c3", title.c_str(), nbins_c3, 0., sumEtMax_c3, nbins_MEX_c3, -50., 50.); 
	dachain->Draw("tcmet*sin(tcmetphi):tcsumet>>RDy_c3" ); 
	RD_c3->Add(RDy_c3,1.0);
										
	//TH2Analyzer TH2AnaRD( &RD, 1, nbins, 1, true);
	TH2Analyzer TH2AnaRD_c3( RD_c3, 1);
	TH1D* haRD_c3 = TH2AnaRD_c3.SigmaGauss();
	haRD_c3->SetFillColor(0);

	//add fitting
	//TF1 *fit = new TF1("fit","sqrt(pow([0],2)+pow([1],2)*(x-[3])+pow([2]*(x-[3]),2))", 0, 60); //in note
	//TF1 *fit = new TF1("fit","[0]+sqrt(pow([1],2)*(x-[3])+pow([2]*(x-[3]),2))", 0, 60);
	//haRD_c3->Fit("fit","R");
	//
	//cout << "f(x)=#sqrt{A^{2}+B^{2}(x-D)+C^{2}(x-D)^{2}}" << endl;
	//string sA = "A = "; sA=sA+Form("%.3f",fit->GetParameter(0))+"\\pm"+Form("%.3f",fit->GetParError(0))+" GeV";
	//string sB = "B = "; sB=sB+Form("%.3f",fit->GetParameter(1))+"\\pm"+Form("%.3f",fit->GetParError(1))+" (GeV)^{1/2}";
	//string sC = "C = "; sC=sC+Form("%.3f",fit->GetParameter(2))+"\\pm"+Form("%.3f",fit->GetParError(2));
	//string sD = "D = "; sD=sD+Form("%.3f",fit->GetParameter(3))+"\\pm"+Form("%.3f",fit->GetParError(3))+ " GeV";
	//cout << sA << endl
	//	 << sB << endl
	//	 << sC << endl
	//	 << sD << endl;
	
	// draw
	haRD_c3->SetMarkerStyle(21);
	//haRD_c3->SetTitle(title.c_str());
	ha_c3->SetTitle("");
	string ylabel = "#sigma(tc#slash{E}_{x,y}) [GeV]";
	string sumetlabel = "tc#scale[0.7]{#sum}E_{T} [GeV]";
	ha_c3->GetXaxis()->SetTitle(sumetlabel.c_str());
	ha_c3->GetYaxis()->SetTitle(ylabel.c_str());
	//ha_c3->GetXaxis()->SetTitle("tcMET components [GeV]");
	ha_c3->GetXaxis()->SetTitleSize(0.06);
	ha_c3->GetXaxis()->SetTitleOffset(0.9);
	ha_c3->GetXaxis()->SetLabelOffset(0.005);
	//ha_c3->GetYaxis()->SetTitle("Events/GeV");
	ha_c3->GetYaxis()->SetTitleOffset(0.8);
	ha_c3->GetYaxis()->SetLabelOffset(0.005);
	ha_c3->SetTickLength(0.03,"XYZ");
	//if (!noMC)
	ha_c3->GetYaxis()->SetRangeUser(0.0,6.0); 
	if( is2tev )
	  ha_c3->GetXaxis()->SetRangeUser(0.0,59.9);
	else
	  ha_c3->GetXaxis()->SetRangeUser(0.0,66.6);

	ha_c3->Draw(); //mc
	//haRD_c3->GetXaxis()->SetTitleSize(0.05);
	//haRD_c3->GetXaxis()->SetTitleOffset(0.95);
	haRD_c3->Draw("same"); //data

	string cmsd2 = "#splitline{CMS Preliminary 2009}{#sqrt{s} = 2360 GeV}";
	string cmsd9 = "#splitline{CMS Preliminary 2009}{#sqrt{s} = 900 GeV}";
	TLatex latex;
	latex.SetTextAlign(12);
	latex.SetTextSize(0.04);
	latex.SetTextFont(62);
	latex.SetNDC();
	double x2 = 0.16; //for the latex--this is the middle of the left edge of the tlatex
	double y2 = 0.88;

	//latex.DrawLatex(0.3, 0.85, (cms).c_str());
	if( is2tev )
	  latex.DrawLatex(x2, y2, (cmsd2).c_str());
	//latex.DrawLatex(0.3, 0.79, (cmsd2).c_str());
	else
	  latex.DrawLatex(x2, y2, (cmsd9).c_str());
	//latex.DrawLatex(0.3, 0.79, (cmsd9).c_str());
										
	TLegend &legend2 = *new TLegend(x2, y2-0.16, x2+0.3, y2-0.06); //y2-0.16 = 0.72  y2-0.06 = 0.82
	legend2.SetTextFont(42); //same as TDR style
	legend2.SetFillColor(0);
	//legend2.SetBorderSize(0); //was in code, but not plots in pas
	legend2.SetShadowColor(0);

	legend2.AddEntry(haRD_c3, "Data");
	if (!noMC) legend2.AddEntry(ha_c3, "Simulation");
	legend2.Draw();
	//c3->Update();
										
	//std::string eps_c3 = "sigmaMEX_tc_fit";
	std::string eps_c3 = "sigmaMEX_tc";
	if( is2tev )
	  eps_c3 += "_2tev";
	//gPad->SaveAs( (eps_c3 + ".eps").c_str() );
	c3->SaveAs( (eps_c3 + ".eps").c_str() );
										
	///////////////////////
	// c4 is apparently the same as c3, or not needed anyway
	///////////////////////
	/*
										
	TCanvas *c4 = new TCanvas("c4", "c4",14,31,700,500);
	c4->SetFillColor(0);
	gStyle->SetOptStat(0);
	
	//const double sumEtMax_c4=110.;
	const double sumEtMax_c4=90.;
	const int nbins_MEX_c4=200;
	const int nbins_c4=12;
										
	//// tc
	TH2D* MC_c4 = new TH2D("MC_c4", title.c_str(), nbins_c4, 0., sumEtMax_c4, nbins_MEX_c4, -50., 50.); 
	MC_c4->GetXaxis()->SetTitle("tcSumEt [GeV]");
	MC_c4->GetYaxis()->SetTitle("tcMEX [GeV]");
	mcchain->Draw("tcmet*cos(tcmetphi):tcsumet>>MC_c4" );
										
	TH2D* MCy_c4 = new TH2D("MCy_c4", title.c_str(), nbins_c4, 0., sumEtMax_c4, nbins_MEX_c4, -50., 50.); 
	mcchain->Draw("tcmet*sin(tcmetphi):tcsumet>>MCy_c4" );
	MC_c4->Add(MCy_c4,1.0);
										
	TH2Analyzer TH2Ana_c4(MC_c4,1);
	//TH2Analyzer TH2Ana( &MC, 1, 100, 10, false);
	//TH2Analyzer TH2Ana( &MC, 1, nbins, 1, true);
	TH1D* ha_c4=TH2Ana_c4.RMS();
	ha_c4->SetLineColor(2);
	ha_c4->SetLineWidth(3);
										
	TH2D* RD_c4 = new TH2D("RD_c4", title.c_str(), nbins_c4, 0., sumEtMax_c4, nbins_MEX_c4, -50., 50.); 
	RD_c4->GetXaxis()->SetTitle("tcSumEt [GeV]");
	RD_c4->GetYaxis()->SetTitle("tcMEX [GeV]");
	dachain->Draw("tcmet*cos(tcmetphi):tcsumet>>RD_c4" ); 
										
	TH2D* RDy_c4 = new TH2D("RDy_c4", title.c_str(), nbins_c4, 0., sumEtMax_c4, nbins_MEX_c4, -50., 50.); 
	dachain->Draw("tcmet*sin(tcmetphi):tcsumet>>RDy_c4" ); 
	RD_c4->Add(RDy_c4,1.0);
										
	//TH2Analyzer TH2AnaRD( &RD, 1, nbins, 1, true);
	TH2Analyzer TH2AnaRD_c4( RD_c4, 1);
	TH1D* haRD_c4=TH2AnaRD_c4.RMS();
	haRD_c4->SetMarkerStyle(21);
	haRD_c4->SetTitle(title.c_str());
	haRD_c4->GetXaxis()->SetTitle("SumEt [GeV]");
	haRD_c4->GetYaxis()->SetTitle("RMS(MEX,MEY) [GeV]");
	haRD_c4->GetYaxis()->SetRangeUser(0.0,7.0); 
										
	// draw !
	haRD_c4->Draw("E1");
	if (!noMC) ha_c4->Draw("E1same");
										
	TLegend *lc4 = new TLegend(0.5502874,0.7669492,0.716954,0.8771186,NULL,"brNDC");
	lc4->SetFillColor(10);
	lc4->AddEntry(haRD_c4,"data, TC","p");
	if (!noMC) lc4->AddEntry(ha_c4,"MC, TC","l");
	lc4->Draw();
	c4->Update();
										
	const std::string eps_c4 = "rmsMEX_tc";
	gPad->SaveAs( (eps_c4 + ".eps").c_str() );
	
	////////////////////////////////////////
	//this was apparently redundant
	*/

  } // doc3c4
}
