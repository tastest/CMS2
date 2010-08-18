#include "TSystem.h"
#include "TChain.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TColor.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TEllipse.h"
#include "TLatex.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <strstream>
//#include <iomanip.h>
#include "TStopwatch.h"

#include "xsecLoop.h"

#include <vector>
//#include <multimap>
//#include <map>

using namespace std;

void xsecPlot( TCanvas * myCanv, TChain* myChain, TString drawThese, const string suffix, const string dir ) {
  std::vector<double> lumiVsRun;
  std::vector<int> numMuVsRun;
  std::vector<int> numElVsRun;
  std::vector<int> numJetVsRun;
  std::vector<int> numclMETVsRun;
  std::vector<int> numtcMETVsRun;
  std::vector<int> numpfMETVsRun;
  int runRange = 160000;
  for(int run = 0; run<runRange; ++run) {
    numMuVsRun.push_back(0);
    numElVsRun.push_back(0);
    numJetVsRun.push_back(0);
    numclMETVsRun.push_back(0);
    numtcMETVsRun.push_back(0);
    numpfMETVsRun.push_back(0);
    lumiVsRun.push_back(0.);
  }
   
  //    Bool_t drawMu = true;
//   Bool_t drawEl = !drawMu;
  Bool_t drawRunningMean = true;
   
  xsecLoop t(myChain);

  if( drawThese.Contains("mus") ) {
    std::cout<<"Counting mus"<<std::endl;
    t.Loop(numMuVsRun, drawThese); 
  }
  else if( drawThese.Contains("els") ) {
    std::cout<<"Counting els"<<std::endl;
    t.Loop(numElVsRun, drawThese); 
  }   
  else if( drawThese.Contains("jets") ) {
    std::cout<<"Counting Jets"<<std::endl;
    t.Loop(numJetVsRun, drawThese); 
  }   
  else if( drawThese.Contains("clmet") ) {
    std::cout<<"Counting clMet>25 events"<<std::endl;
    t.Loop(numclMETVsRun, drawThese); 
  }   
  else if( drawThese.Contains("tcmet") ) {
    std::cout<<"Counting tcMet>25 events"<<std::endl;
    t.Loop(numtcMETVsRun, drawThese); 
  }   
  else if( drawThese.Contains("pfmet") ) {
    std::cout<<"Counting pfMet>25 events"<<std::endl;
    t.Loop(numpfMETVsRun, drawThese); 
  }   

  else {
    std::cout<<"BAD drawable in xsecPlot: "<<drawThese <<" is not forseen. Bad things will happen - better fix this!"<<endl;
  }

  //    for(int run = 0; run<runRange; ++run){
  //      if(numMuVsRun[run] > 0) {
  //        cout<<"FINAL Mu run "<<run<<" nmu "<<numMuVsRun[run]<<std::endl;;
  //      }
  //    }
   
  t.readFile("lumi_by_run.txt", lumiVsRun);
   
  // calculate the running mean, here for +- 12 runs
  std::vector<float> nMuRoller; 
  std::vector<float> nElRoller; 
  std::vector<float> nJetRoller; 
  std::vector<float> nclMETRoller; 
  std::vector<float> ntcMETRoller; 
  std::vector<float> npfMETRoller; 
  std::vector<float> nLumiRoller;
#ifndef __CINT__
  if(drawRunningMean) {
    nMuRoller     = rollingThunder(numMuVsRun, lumiVsRun, 12);
    nElRoller     = rollingThunder(numElVsRun, lumiVsRun, 12);
    nJetRoller    = rollingThunder(numJetVsRun, lumiVsRun, 12);
    nclMETRoller  = rollingThunder(numclMETVsRun, lumiVsRun, 12);
    ntcMETRoller  = rollingThunder(numtcMETVsRun, lumiVsRun, 12);
    npfMETRoller  = rollingThunder(numpfMETVsRun, lumiVsRun, 12);
    nLumiRoller   = rollingThunder(lumiVsRun, lumiVsRun, 12);
  }
#endif
  TGraphErrors * aGraph = new TGraphErrors();
  TGraphErrors * anotherGraph = new TGraphErrors();
  int graphpointnr = 0;

  std::vector<int> numLepVsRun;
  std::vector<float> nLepRoller; 
  if( drawThese.Contains("mus") ) {
    numLepVsRun = numMuVsRun;
    nLepRoller = nMuRoller; 
  }
  else if( drawThese.Contains("els") ) {
    numLepVsRun = numElVsRun;
    nLepRoller = nElRoller; 
  }
  else if( drawThese.Contains("jets") ) {
    numLepVsRun = numJetVsRun;
    nLepRoller = nJetRoller; 
  }
  else if( drawThese.Contains("clmet") ) {
    numLepVsRun = numclMETVsRun;
    nLepRoller = nclMETRoller; 
  }
  else if( drawThese.Contains("tcmet") ) {
    numLepVsRun = numtcMETVsRun;
    nLepRoller = ntcMETRoller; 
  }
  else if( drawThese.Contains("pfmet") ) {
    numLepVsRun = numpfMETVsRun;
    nLepRoller = npfMETRoller; 
  }


  for(int run = 0; run<runRange; ++run){
    if( lumiVsRun[run] != 0 && numLepVsRun[run] != 0 ) {
      aGraph->SetPoint( graphpointnr, run, numLepVsRun[run]/lumiVsRun[run] );
      aGraph->SetPointError( graphpointnr, 0, sqrt((float)numLepVsRun[run])/lumiVsRun[run] );
       
      anotherGraph->SetPoint( graphpointnr, run+0.4, nLepRoller[run]/nLumiRoller[run] );
      anotherGraph->SetPointError( graphpointnr, 0, sqrt((float)nLepRoller[run])/nLumiRoller[run] );
       
      graphpointnr +=1;
    }
  }
   
  aGraph->SetFillColor(kWhite);
  if( drawThese.Contains("mus") ) {
    aGraph->SetTitle("Muons per run per lumi");
    aGraph->GetYaxis()->SetTitle("Muon rate, #sigma_{GoodMu Vis p_{T} > 10 GeV} (#mub)");
  }
  else if( drawThese.Contains("els") ) {
    aGraph->SetTitle("Electrons per run per lumi ");
    aGraph->GetYaxis()->SetTitle("Ele rate, #sigma_{GoodEl Vis p_{T} > 10 GeV} (#mub)");
  }
  else if( drawThese.Contains("jets") ) {
    aGraph->SetTitle("Jets per run per lumi ");
    aGraph->GetYaxis()->SetTitle("Jet rate, #sigma_{Jet Vis p_{T}^{corr} > 30 GeV} (#mub)");
  }
  else if( drawThese.Contains("clmet") ) {
    aGraph->SetTitle("clMET > 25 GeV per run per lumi ");
    aGraph->GetYaxis()->SetTitle("high clMET rate, #sigma_{clMET Vis > 25 GeV} (#mub)");
  }
  else if( drawThese.Contains("tcmet") ) {
    aGraph->SetTitle("tcMET > 25 GeV per run per lumi ");
    aGraph->GetYaxis()->SetTitle("high tcMET rate, #sigma_{tcMET Vis > 25 GeV} (#mub)");
  }
  else if( drawThese.Contains("pfmet") ) {
    aGraph->SetTitle("pfMET > 25 GeV per run per lumi ");
    aGraph->GetYaxis()->SetTitle("high pfMET rate, #sigma_{pfMET Vis > 25 GeV} (#mub)");
  }
  aGraph->GetYaxis()->SetTitleOffset(1.2);
  aGraph->GetXaxis()->SetTitle("run number"); //6  
   
  aGraph->SetMarkerColor(4); 
  aGraph->SetMarkerStyle(21);
  aGraph->GetXaxis()->SetNoExponent();
  aGraph->Draw("AP");

  TLatex *tex = new TLatex(0.2,0.75,"#Box blue: Run-by-run Lepton Rate");
  tex->SetNDC();
  tex->SetTextColor(4);
  tex->SetLineWidth(2);
  tex->Draw();
  TLatex *tex2 = new TLatex(0.2,0.68,"#Delta green: 25 Run Running Mean");
  tex2->SetNDC();
  tex2->SetTextColor(8);
  tex2->SetLineWidth(2);
  tex2->Draw();
  TLatex *tex3 = new TLatex(0.2,0.61,"Orange highlight: expect>0.8 events, see none");
  tex3->SetNDC();
  tex3->SetTextColor(kOrange);
  tex3->SetLineWidth(2);
  tex3->Draw();
  TLatex *tex4 = new TLatex(0.2,0.54,"Red highlight: >3#sigma too many or few events");
  tex4->SetNDC();
  tex4->SetTextColor(kRed);
  tex4->SetLineWidth(2);
  tex4->Draw();


  TF1 * pol0Func = new TF1("pol0Func","pol0");
  aGraph->Fit("pol0Func");
   
  double fitMean = pol0Func->GetParameter(0);
  cout<<fitMean<<endl;
   
  anotherGraph->SetMarkerColor(kGreen); 
  anotherGraph->SetMarkerStyle(22);
  if(drawRunningMean) {
    anotherGraph->GetXaxis()->SetNoExponent();
    anotherGraph->Draw("P");
  }
   
  //   graphpointnr = 0;
  for(int run = 0; run<runRange; ++run){
    if( lumiVsRun[run] != 0 && numLepVsRun[run] != 0 ) {
      if( TMath::Abs( ( (numLepVsRun[run]/lumiVsRun[run]) - fitMean) / (sqrt((float)numLepVsRun[run])/lumiVsRun[run]) ) > 3. ) {
        cout<<"EXCESS/DEFICIT Run: "<<run<<" Lep ratio "<<numLepVsRun[run]/lumiVsRun[run]<<" normPull "<< ( (numLepVsRun[run]/lumiVsRun[run]) - fitMean) / ( sqrt( (float)numLepVsRun[run] )/lumiVsRun[run] ) <<" nLep "<<numLepVsRun[run]<<" lumi "<<lumiVsRun[run]<<std::endl;
         
        TEllipse *ellipse = new TEllipse(run, numLepVsRun[run]/lumiVsRun[run]  ,10, ( sqrt( (float)numLepVsRun[run] )/lumiVsRun[run] )  ,0,360,0);
        ellipse->SetFillStyle(0);
        ellipse->SetLineColor(kRed);
        ellipse->SetLineWidth(3);
        ellipse->Draw();
      }
      //       graphpointnr +=1;
    }
    if(lumiVsRun[run] != 0 && numLepVsRun[run] == 0 && (fitMean*lumiVsRun[run]) > 0.8 ) {
      TEllipse *ellipse = new TEllipse(run, fitMean  ,1, ( sqrt( fitMean*lumiVsRun[run] )/lumiVsRun[run] ) ,0,360,0);
      ellipse->SetFillStyle(0);
      ellipse->SetLineColor(kOrange);
      ellipse->SetLineWidth(3);
      ellipse->Draw();
      cout<<"NOLEP Run: "<<run<<" nLep expected "<<fitMean*lumiVsRun[run]<<" lumi "<<lumiVsRun[run]<<std::endl;
    }
  }
  if( drawThese.Contains("mus") ) {
    myCanv->SaveAs((dir+"MuXsecNew"+suffix).c_str());
  }
  else if( drawThese.Contains("els") ) {
    myCanv->SaveAs((dir+"ElXsecNew"+suffix).c_str());
  }
  else if( drawThese.Contains("jets") ) {
    myCanv->SaveAs((dir+"JetXsecNew"+suffix).c_str());
  }
  else if( drawThese.Contains("clmet") ) {
    myCanv->SaveAs((dir+"clMETXsecNew"+suffix).c_str());
  }
  else if( drawThese.Contains("tcmet") ) {
    myCanv->SaveAs((dir+"tcMETXsecNew"+suffix).c_str());
  }
  else if( drawThese.Contains("pfmet") ) {
    myCanv->SaveAs((dir+"tcMETXsecNew"+suffix).c_str());
  }
}


void plotMu(void) {
  //  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->SetStyle("Plain");
  gStyle->SetHistFillColor(kRed);
  gStyle->SetHistMinimumZero();
  gStyle->SetOptStat("nemoui");

  const string suffix = ".png";
  const string dir = "plots/";

  int   N_nbins = 2;
  float N_min = 0.5;
  float N_max = 2.5;

  int pt_nbins = 40;
  float pt_min = 10;
  float pt_max = 50;

  int eta_nbins = 25;
  float eta_min = -2.5;
  float eta_max = 2.5;

  int phi_nbins = 30;
  float phi_min = -3.15;
  float phi_max = 3.15;

  int d0_nbins = 50;
  float d0_min = -.025;
  float d0_max = .025;

  int Iso_nbins = 100;
  float Iso_min = 0;
  float Iso_max = 1;

  int   met_nbins = 75;
  float met_min   = 0;
  float met_max   = 75;

  int jet_pt_nbins  = 70;
  float jet_pt_min  = 30;
  float jet_pt_max  = 100;

  int jet_eta_nbins = 25;
  float jet_eta_min = -2.4;
  float jet_eta_max = 2.4;

  int Njet_nbins  = 6;
  int Njet_min    = 0;
  int Njet_max    = 6; 

  // muon quantities - before
  TH1F* mu1_N        = new TH1F("mu1_Nmus", "mu1_Nmus", N_nbins, N_min, N_max );
  TH1F* mu1_Niso     = new TH1F("mu1_Nmusiso", "mu1_Nmus_isolated", N_nbins, N_min, N_max );
  TH1F* mu1_pt       = new TH1F("mu1_pt", "mu1_pt", pt_nbins, pt_min, pt_max );
  TH1F* mu1_eta      = new TH1F("mu1_eta", "mu1_eta", eta_nbins, eta_min, eta_max );
  TH1F* mu1_phi      = new TH1F("mu1_phi", "mu1_phi", phi_nbins, phi_min, phi_max );
  TH1F* mu1_d0corr   = new TH1F("mu1_d0corr", "mu1_d0corr", d0_nbins, d0_min, d0_max );
  TH1F* mu1_Iso      = new TH1F("mu1_Iso", "mu1_Iso", Iso_nbins, Iso_min, Iso_max );
  TH2F* mu1_phi_d0   = new TH2F("mu1_phi_d0", "mu1_phi_d0", phi_nbins, phi_min, phi_max, d0_nbins, d0_min, d0_max );

  // met & met phi - before
  TH1F* mu1_clmet        = new TH1F("mu1_clmet", "mu1_clmet", met_nbins, met_min, met_max );
  TH1F* mu1_pfmet        = new TH1F("mu1_pfmet", "mu1_pfmet", met_nbins, met_min, met_max );
  TH1F* mu1_tcmet        = new TH1F("mu1_tcmet", "mu1_tcmet", met_nbins, met_min, met_max );
  TH1F* mu1_clmetphi     = new TH1F("mu1_clmetphi", "mu1_clmetphi", phi_nbins, phi_min, phi_max );
  TH1F* mu1_pfmetphi     = new TH1F("mu1_pfmetphi", "mu1_pfmetphi", phi_nbins, phi_min, phi_max );
  TH1F* mu1_tcmetphi     = new TH1F("mu1_tcmetphi", "mu1_tcmetphi", phi_nbins, phi_min, phi_max );

  // jets - before

  TH1F* mu1_Njets     = new TH1F("mu1_Njets", "mu1_Ncalojets", Njet_nbins, Njet_min, Njet_max );
  TH1F* mu1_Npfjets   = new TH1F("mu1_Npfjets", "mu1_Npfjets", Njet_nbins, Njet_min, Njet_max );
  TH1F* mu1_Ntrkjets  = new TH1F("mu1_Ntrkjets", "mu1_Ntrkjets", Njet_nbins, Njet_min, Njet_max );

  TH1F* mu1_jets_pt   	= new TH1F("mu1_jets_pt", "mu1_calojets_pt", jet_pt_nbins, jet_pt_min, jet_pt_max );
  TH1F* mu1_pfjets_pt   	= new TH1F("mu1_pfjets_pt", "mu1_pfjets_pt", jet_pt_nbins, jet_pt_min, jet_pt_max );
  TH1F* mu1_trkjets_pt   	= new TH1F("mu1_trkjets_pt", "mu1_trkjets_pt", jet_pt_nbins, jet_pt_min, jet_pt_max );

  TH1F* mu1_jets_eta   	= new TH1F("mu1_jets_eta", "mu1_calojets_eta", jet_eta_nbins, jet_eta_min, jet_eta_max );
  TH1F* mu1_pfjets_eta   	= new TH1F("mu1_pfjets_eta", "mu1_pfjets_eta", jet_eta_nbins, jet_eta_min, jet_eta_max );
  TH1F* mu1_trkjets_eta   = new TH1F("mu1_trkjets_eta", "mu1_trkjets_eta", jet_eta_nbins, jet_eta_min, jet_eta_max );

  TH1F* mu1_jets_phi   	= new TH1F("mu1_jets_phi", "mu1_calojets_phi", phi_nbins, phi_min, phi_max );
  TH1F* mu1_pfjets_phi   	= new TH1F("mu1_pfjets_phi", "mu1_pfjets_phi", phi_nbins, phi_min, phi_max );
  TH1F* mu1_trkjets_phi   = new TH1F("mu1_trkjets_phi", "mu1_trkjets_phi", phi_nbins, phi_min, phi_max );

  // muons - after
  TH1F* mu2_N        = new TH1F("mu2_Nmus", "mu2_Nmus", N_nbins, N_min, N_max );
  TH1F* mu2_Niso     = new TH1F("mu2_Nmusiso", "mu2_Nmus_isolated", N_nbins, N_min, N_max );
  TH1F* mu2_pt       = new TH1F("mu2_pt", "mu2_pt", pt_nbins, pt_min, pt_max );
  TH1F* mu2_eta      = new TH1F("mu2_eta", "mu2_eta", eta_nbins, eta_min, eta_max );
  TH1F* mu2_phi      = new TH1F("mu2_phi", "mu2_phi", phi_nbins, phi_min, phi_max );
  TH1F* mu2_d0corr   = new TH1F("mu2_d0corr", "mu2_d0corr", d0_nbins, d0_min, d0_max );
  TH1F* mu2_Iso      = new TH1F("mu2_Iso", "mu2_Iso", Iso_nbins, Iso_min, Iso_max );
  TH2F* mu2_phi_d0   = new TH2F("mu2_phi_d0", "mu2_phi_d0", phi_nbins, phi_min, phi_max, d0_nbins, d0_min, d0_max );

  // met & met phi - after
  TH1F* mu2_clmet        = new TH1F("mu2_clmet", "mu2_clmet", met_nbins, met_min, met_max );
  TH1F* mu2_pfmet        = new TH1F("mu2_pfmet", "mu2_pfmet", met_nbins, met_min, met_max );
  TH1F* mu2_tcmet        = new TH1F("mu2_tcmet", "mu2_tcmet", met_nbins, met_min, met_max );
  TH1F* mu2_clmetphi     = new TH1F("mu2_clmetphi", "mu2_clmetphi", phi_nbins, phi_min, phi_max );
  TH1F* mu2_pfmetphi     = new TH1F("mu2_pfmetphi", "mu2_pfmetphi", phi_nbins, phi_min, phi_max );
  TH1F* mu2_tcmetphi     = new TH1F("mu2_tcmetphi", "mu2_tcmetphi", phi_nbins, phi_min, phi_max );

  // jets - after

  TH1F* mu2_Njets     = new TH1F("mu2_Njets", "mu2_Ncalojets", Njet_nbins, Njet_min, Njet_max );
  TH1F* mu2_Npfjets   = new TH1F("mu2_Npfjets", "mu2_Npfjets", Njet_nbins, Njet_min, Njet_max );
  TH1F* mu2_Ntrkjets  = new TH1F("mu2_Ntrkjets", "mu2_Ntrkjets", Njet_nbins, Njet_min, Njet_max );

  TH1F* mu2_jets_pt   	= new TH1F("mu2_jets_pt", "mu2_calojets_pt", jet_pt_nbins, jet_pt_min, jet_pt_max );
  TH1F* mu2_pfjets_pt   	= new TH1F("mu2_pfjets_pt", "mu2_pfjets_pt", jet_pt_nbins, jet_pt_min, jet_pt_max );
  TH1F* mu2_trkjets_pt   	= new TH1F("mu2_trkjets_pt", "mu2_trkjets_pt", jet_pt_nbins, jet_pt_min, jet_pt_max );
  TH1F* mu2_jets_eta   	= new TH1F("mu2_jets_eta", "mu2_calojets_eta", jet_eta_nbins, jet_eta_min, jet_eta_max );
  TH1F* mu2_pfjets_eta   	= new TH1F("mu2_pfjets_eta", "mu2_pfjets_eta", jet_eta_nbins, jet_eta_min, jet_eta_max );
  TH1F* mu2_trkjets_eta   = new TH1F("mu2_trkjets_eta", "mu2_trkjets_eta", jet_eta_nbins, jet_eta_min, jet_eta_max );
  TH1F* mu2_jets_phi   	= new TH1F("mu2_jets_phi", "mu2_calojets_phi", phi_nbins, phi_min, phi_max );
  TH1F* mu2_pfjets_phi   	= new TH1F("mu2_pfjets_phi", "mu2_pfjets_phi", phi_nbins, phi_min, phi_max );
  TH1F* mu2_trkjets_phi   = new TH1F("mu2_trkjets_phi", "mu2_trkjets_phi", phi_nbins, phi_min, phi_max );

  // Chains
  TChain *chain1 = new TChain("tree");
  TChain *chain2 = new TChain("tree");
  chain1->Add("validate_mus_before.root");
  if( 42 == 42 ) {
    chain2->Add("validate_mus_after.root");
  }
  else{
    std::cout<<"ALARM! USING non standard \"after\" chain!!"<<std::endl;
    std::cout<<"ALARM! USING non standard \"after\" chain!!"<<std::endl;
    std::cout<<"ALARM! USING non standard \"after\" chain!!"<<std::endl;
    //    chain2->Add("validate_mus_after_dilepskim.root");
    chain2->Add("validate_mus_after_ptGt20.root");
  }

  if( 42 == 42 ) { 

    // Fill Before
    TCanvas *ctemp = new TCanvas();
    chain1->Draw("nmu >> mu1_Nmus");
    chain1->Draw("Sum$(musid) >> mu1_Nmusiso"); //isolated bc musid is isolated
    chain1->Draw("musp4.pt() >> mu1_pt", "musid");
    chain1->Draw("musp4.eta() >> mu1_eta", "musid");
    chain1->Draw("musp4.phi() >> mu1_phi", "musid");
    chain1->Draw("musd0corr >> mu1_d0corr", "musid");
    chain1->Draw("musiso >> mu1_Iso");
    chain1->Draw("musd0corr:musp4.phi() >> mu1_phi_d0", "musid", "BOX");

    chain1->Draw("clmet >> mu1_clmet");
    chain1->Draw("pfmet >> mu1_pfmet");
    chain1->Draw("tcmet >> mu1_tcmet");

    chain1->Draw("clmetphi >> mu1_clmetphi", "clmet>10");
    chain1->Draw("pfmetphi >> mu1_pfmetphi", "pfmet>10");
    chain1->Draw("tcmetphi >> mu1_tcmetphi", "tcmet>10");

    chain1->Draw("jets@.size() >> mu1_Njets");
    chain1->Draw("pfjets@.size() >> mu1_Npfjets");
    chain1->Draw("trkjets@.size() >> mu1_Ntrkjets");

    chain1->Draw("jets.pt() >> mu1_jets_pt");
    chain1->Draw("pfjets.pt() >> mu1_pfjets_pt");
    chain1->Draw("trkjets.pt() >> mu1_trkjets_pt");
    chain1->Draw("jets.eta() >> mu1_jets_eta");
    chain1->Draw("pfjets.eta() >> mu1_pfjets_eta");
    chain1->Draw("trkjets.eta() >> mu1_trkjets_eta");
    chain1->Draw("jets.phi() >> mu1_jets_phi");
    chain1->Draw("pfjets.phi() >> mu1_pfjets_phi");
    chain1->Draw("trkjets.phi() >> mu1_trkjets_phi");

    // Fill After
    chain2->Draw("nmu >> mu2_Nmus");
    chain2->Draw("Sum$(musid) >> mu2_Nmusiso");
    chain2->Draw("musp4.pt() >> mu2_pt", "musid");
    chain2->Draw("musp4.eta() >> mu2_eta", "musid");
    chain2->Draw("musp4.phi() >> mu2_phi", "musid");
    chain2->Draw("musd0corr >> mu2_d0corr", "musid");
    chain2->Draw("musiso >> mu2_Iso");
    chain2->Draw("musd0corr:musp4.phi() >> mu2_phi_d0", "musid", "BOX PROJ");

    chain2->Draw("clmet >> mu2_clmet");
    chain2->Draw("pfmet >> mu2_pfmet");
    chain2->Draw("tcmet >> mu2_tcmet");

    chain2->Draw("clmetphi >> mu2_clmetphi", "clmet>10");
    chain2->Draw("pfmetphi >> mu2_pfmetphi", "pfmet>10");
    chain2->Draw("tcmetphi >> mu2_tcmetphi", "tcmet>10");

    chain2->Draw("jets@.size() >> mu2_Njets");
    chain2->Draw("pfjets@.size() >> mu2_Npfjets");
    chain2->Draw("trkjets@.size() >> mu2_Ntrkjets");

    chain2->Draw("jets.pt() >> mu2_jets_pt");
    chain2->Draw("pfjets.pt() >> mu2_pfjets_pt");
    chain2->Draw("trkjets.pt() >> mu2_trkjets_pt");
    chain2->Draw("jets.eta() >> mu2_jets_eta");
    chain2->Draw("pfjets.eta() >> mu2_pfjets_eta");
    chain2->Draw("trkjets.eta() >> mu2_trkjets_eta");
    chain2->Draw("jets.phi() >> mu2_jets_phi");
    chain2->Draw("pfjets.phi() >> mu2_pfjets_phi");
    chain2->Draw("trkjets.phi() >> mu2_trkjets_phi");

    delete ctemp;

    TCanvas *c1 = new TCanvas();
    c1->SetWindowSize(1100,850);
    c1->Divide(3,2);
    c1->cd(1)->SetLogy();
    mu1_pt->Draw();
    c1->cd(2);
    mu1_eta->Draw();
    c1->cd(3);
    mu1_phi->Draw();
    c1->cd(4)->SetLogy();
    mu2_pt->Draw();
    c1->cd(5);
    mu2_eta->Draw();
    c1->cd(6);
    mu2_phi->Draw();
    c1->SaveAs((dir+"compare_mu_kinematics"+suffix).c_str());

    TCanvas *c2 = new TCanvas();
    c2->SetWindowSize(1450,850);//wider bc 4 plots
    c2->Divide(4,2);
    c2->cd(1)->SetLogy();
    mu1_N->Draw();
    c2->cd(2)->SetLogy();
    mu1_Niso->Draw();
    c2->cd(3);
    mu1_d0corr->Draw();
    c2->cd(4);
    mu1_Iso->Draw();
    c2->cd(5)->SetLogy();
    mu2_N->Draw();
    c2->cd(6)->SetLogy();
    mu2_Niso->Draw();
    c2->cd(7);
    mu2_d0corr->Draw();
    c2->cd(8);
    mu2_Iso->Draw();
    c2->SaveAs((dir+"compare_mu_N_do_iso"+suffix).c_str());

    TCanvas *c3 = new TCanvas();
    c3->SetWindowSize(1100,900);
    c3->Divide(2,2);
    c3->cd(1);
    mu1_phi_d0->Draw("BOX");
    c3->cd(2);
    mu2_phi_d0->Draw("BOX");
    c3->cd(3);
    mu1_phi_d0->ProfileX()->Draw();
    c3->cd(4);
    mu2_phi_d0->ProfileX()->Draw();
    c3->SaveAs((dir+"compare_mu_2dkin"+suffix).c_str());

    TCanvas *c4 = new TCanvas();
    c4->SetWindowSize(1100,850);
    c4->Divide(3,2);
    c4->cd(1)->SetLogy();
    mu1_clmet->Draw();
    c4->cd(2)->SetLogy();
    mu1_pfmet->Draw();
    c4->cd(3)->SetLogy();
    mu1_tcmet->Draw();
    c4->cd(4)->SetLogy();
    mu2_clmet->Draw();
    c4->cd(5)->SetLogy();
    mu2_pfmet->Draw();
    c4->cd(6)->SetLogy();
    mu2_tcmet->Draw();
    c4->SaveAs((dir+"compare_mu_met"+suffix).c_str());

    TCanvas *c5 = new TCanvas();
    c5->SetWindowSize(1100,850);
    c5->Divide(3,2);
    c5->cd(1);
    mu1_clmetphi->Draw();
    c5->cd(2);
    mu1_pfmetphi->Draw();
    c5->cd(3);
    mu1_tcmetphi->Draw();
    c5->cd(4);
    mu2_clmetphi->Draw();
    c5->cd(5);
    mu2_pfmetphi->Draw();
    c5->cd(6);
    mu2_tcmetphi->Draw();
    c5->SaveAs((dir+"compare_mu_metphi"+suffix).c_str());

    TCanvas *c = new TCanvas();
    c->SetWindowSize(1100,850);
    c->Divide(3,2);
    c->cd(1);
    mu1_Njets->Draw();
    c->cd(2);
    mu1_Npfjets->Draw();
    c->cd(3);
    mu1_Ntrkjets->Draw();
    c->cd(4);
    mu2_Njets->Draw();
    c->cd(5);
    mu2_Npfjets->Draw();
    c->cd(6);
    mu2_Ntrkjets->Draw();
    c->SaveAs((dir+"compare_mu_jet"+suffix).c_str());

    TCanvas *c6 = new TCanvas();
    c6->SetWindowSize(1100,850);
    c6->Divide(3,2);
    c6->cd(1);
    mu1_jets_pt->Draw();
    c6->cd(2);
    mu1_jets_eta->Draw();
    c6->cd(3);
    mu1_jets_phi->Draw();
    c6->cd(4);
    mu2_jets_pt->Draw();
    c6->cd(5);
    mu2_jets_eta->Draw();
    c6->cd(6);
    mu2_jets_phi->Draw();
    c6->SaveAs((dir+"compare_mu_cljetkin"+suffix).c_str());

    TCanvas *c7 = new TCanvas();
    c7->SetWindowSize(1100,850);
    c7->Divide(3,2);
    c7->cd(1);
    mu1_pfjets_pt->Draw();
    c7->cd(2);
    mu1_pfjets_eta->Draw();
    c7->cd(3);
    mu1_pfjets_phi->Draw();
    c7->cd(4);
    mu2_pfjets_pt->Draw();
    c7->cd(5);
    mu2_pfjets_eta->Draw();
    c7->cd(6);
    mu2_pfjets_phi->Draw();
    c7->SaveAs((dir+"compare_mu_pfjetkin"+suffix).c_str());

    TCanvas *c8 = new TCanvas();
    c8->SetWindowSize(1100,850);
    c8->Divide(3,2);
    c8->cd(1);
    mu1_trkjets_pt->Draw();
    c8->cd(2);
    mu1_trkjets_eta->Draw();
    c8->cd(3);
    mu1_trkjets_phi->Draw();
    c8->cd(4);
    mu2_trkjets_pt->Draw();
    c8->cd(5);
    mu2_trkjets_eta->Draw();
    c8->cd(6);
    mu2_trkjets_phi->Draw();
    c8->SaveAs((dir+"compare_mu_tkjetkin"+suffix).c_str());

  }

  //
  // Now for the "lepton cross section" plots:
  //
  //plot reference chain (currently just for lepton)
  TCanvas *c9 = new TCanvas();
  xsecPlot(c9, chain1, "mus", suffix, dir);
  
  //plot new chain
  TCanvas *c10 = new TCanvas();
  xsecPlot(c10, chain2, "mus", suffix, dir);
  
  // temporarily try out jets here
  TCanvas *c11 = new TCanvas();
  xsecPlot(c11, chain2, "jets", suffix, dir);
  
  // temporarily try out clmet here
  TCanvas *c12 = new TCanvas();
  xsecPlot(c12, chain2, "clmet", suffix, dir);
  
  // temporarily try out tcmet here
  TCanvas *c13 = new TCanvas();
  xsecPlot(c13, chain2, "tcmet", suffix, dir);
  
  // temporarily try out pfmet here
  TCanvas *c14 = new TCanvas();
  xsecPlot(c14, chain2, "pfmet", suffix, dir);

}
