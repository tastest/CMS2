// Usage:
//   gROOT->LoadMacro("monitor.C+");
//   TChain* chain1 = new TChain("Events","TTbar-madgraph");
//   chain1->Add("/data/tmp/dmytro/cms2-V01-02-01/TTJets-madgraph/*.root");
//   TChain* chain2 = new TChain("Events","WW-2l");
//   chain2->Add("/data/tmp/dmytro/cms2-V01-02-01/WW_2l-Pythia/merged_ntuple.root");
//   drawAll(chain1,chain2);
//
// Up to 9 samples can be shown at once
//
#include "TH1F.h"
#include "TTree.h"
#include <vector>
#include <string>
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TChain.h"
#include <iostream>
#include "TProfile.h"

struct variable{
  const char* title;
  const char* varexp;
  const char* selection;
  unsigned int nbins;
  float varmin;
  float varmax;
  const char* varexp2;
  unsigned int nbins2;
  float varmin2;
  float varmax2;
  
  variable( const char* ititle, const char* ivarexp, const char* iselection,
	    unsigned int inbins, float ivarmin, float ivarmax,
	    const char* ivarexp2=0, unsigned int inbins2=0, float ivarmin2=0, float ivarmax2=0):
    title(ititle), varexp(ivarexp), selection(iselection), 
    nbins(inbins), varmin(ivarmin), varmax(ivarmax),
    varexp2(ivarexp2), nbins2(inbins2), varmin2(ivarmin2), varmax2(ivarmax2)
  {}
};


void drawAll( TTree* tree1, 
	      TTree* tree2 = 0, 
	      TTree* tree3 = 0, 
	      TTree* tree4 = 0, 
	      TTree* tree5 = 0, 
	      TTree* tree6 = 0, 
	      TTree* tree7 = 0, 
	      TTree* tree8 = 0,
	      TTree* tree9 = 0 )
{
  std::string psname = "monitor.ps";

  // Variables to monitor
  std::vector<variable> vars;
  vars.push_back( variable("Electron d0 corrected for beam spot (pt>20)", "els_d0corr", "els_p4.pt()>20", 400, -0.1, 0.1 ) );
  vars.push_back( variable("Electron d0 corrected for beam spot (pt>20) vs phi", "els_d0corr", "els_p4.pt()>20", 100, -3.14, 3.14, "els_p4.phi()",0,-.1,.1 ) );
  vars.push_back( variable("Electron MC match PDG id (no ID,pt>20)", "abs(els_mc_id)", "els_p4.pt()>20", 350, 0, 350 ) );
  vars.push_back( variable("Electron MC match PDG id (robust ID,pt>20)", "abs(els_mc_id)", "els_egamma_robustId&&els_p4.pt()>20", 350, 0, 350 ) );
  vars.push_back( variable("Electron MC match PDG id (loose ID,pt>20)", "abs(els_mc_id)", "els_egamma_looseId&&els_p4.pt()>20", 350, 0, 350 ) );
  vars.push_back( variable("Electron MC match PDG id (tight ID,pt>20)", "abs(els_mc_id)", "els_egamma_tightId&&els_p4.pt()>20", 350, 0, 350 ) );
  vars.push_back( variable("Electron trk iso (pt>20,zeros are not shown)", "els_tkIso", "els_p4.pt()>20", 100, 0.001, 20 ) );
  vars.push_back( variable("Electron PAT trk iso (pt>20,zeros are not shown)", "els_pat_trackIso", "els_p4.pt()>20", 100, 0.001, 20 ) );
  vars.push_back( variable("Electron ECAL Jurassic basic cluster iso (pt>20,zeros are not shown)", "els_ecalJuraIso", "els_p4.pt()>20", 100, 0.001, 20 ) );
  vars.push_back( variable("Electron PAT ECAL iso (pt>20,zeros are not shown)", "els_pat_ecalIso", "els_p4.pt()>20", 100, 0.001, 20 ) );
  vars.push_back( variable("Electron HCAL tower cone iso (pt>20,zeros are not shown)", "els_hcalConeIso", "els_p4.pt()>20", 100, 0.001, 20 ) );
  vars.push_back( variable("Electron PAT HCAL iso (pt>20,zeros are not shown)", "els_pat_hcalIso", "els_p4.pt()>20", 100, 0.001, 20 ) );
  
  vars.push_back( variable("Muon d0 corrected for beam spot (pt>20)", "mus_d0corr", "mus_p4.pt()>20", 400, -0.1, 0.1 ) );
  vars.push_back( variable("Muon MC match PDG id (no ID,pt>20)", "abs(mus_mc_id)", "mus_p4.pt()>20", 350, 0, 350 ) );
  vars.push_back( variable("Muon type (true MC muon,pt>20)", "mus_type", "abs(mus_mc_id)==13&&mus_p4.pt()>20", 20, 0, 20 ) );
  vars.push_back( variable("Muon type (not MC muon,pt>20)", "mus_type", "abs(mus_mc_id)!=13&&mus_p4.pt()>20", 20, 0, 20 ) );
  vars.push_back( variable("Global muon normalized chi2 (pt>20)", "mus_gfit_chi2/mus_gfit_ndof","mus_type&0x8&&mus_p4.pt()>20", 100, 0, 20 ) );
  vars.push_back( variable("Number of valid hits in inner track (global and tracker muons,pt>20)", "mus_validHits","mus_type&14&&mus_p4.pt()>20", 50, 0, 50 ) );
  vars.push_back( variable("Muon trk R0.3 iso (pt>20,zeros are not shown)", "mus_iso03_sumPt", "mus_p4.pt()>20", 100, 0.001, 20 ) );
  vars.push_back( variable("Muon PAT trk iso (pt>20,zeros are not shown)", "mus_pat_trackIso", "mus_p4.pt()>20", 100, 0.001, 20 ) );
  vars.push_back( variable("Muon ECAL R0.3 iso (pt>20,zeros are not shown)", "mus_iso03_emEt", "mus_p4.pt()>20", 100, 0.001, 20 ) );
  vars.push_back( variable("Muon PAT ECAL iso (pt>20,zeros are not shown)", "mus_pat_ecalIso", "mus_p4.pt()>20", 100, 0.001, 20 ) );
  vars.push_back( variable("Muon HCAL R0.3 iso (pt>20,zeros are not shown)", "mus_iso03_hadEt", "mus_p4.pt()>20", 100, 0.001, 20 ) );
  vars.push_back( variable("Muon PAT HCAL iso (pt>20,zeros are not shown)", "mus_pat_hcalIso", "mus_p4.pt()>20", 100, 0.001, 20 ) );

  vars.push_back( variable("Raw uncorrected MET", "evt_met", "", 100, 0, 100 ) );
  vars.push_back( variable("tcMET", "evt_tcmet", "", 100, 0, 100 ) );
  vars.push_back( variable("Dilepton mass (no cuts)", "hyp_p4.M()", "", 150, 0, 150 ) );

  vars.push_back( variable("Hypothesis jet count", "hyp_njets", "", 10, 0, 10 ) );

  std::vector<TTree*> trees;
  trees.push_back(tree1);
  if (tree2) trees.push_back(tree2);
  if (tree3) trees.push_back(tree3);
  if (tree4) trees.push_back(tree4);
  if (tree5) trees.push_back(tree5);
  if (tree6) trees.push_back(tree6);
  if (tree7) trees.push_back(tree7);
  if (tree8) trees.push_back(tree8);
  if (tree9) trees.push_back(tree9);
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("nemruosk");
  gStyle->SetPalette(1);
  TCanvas* c1 = new TCanvas("c1","c1", 1100,850);

  for ( unsigned int ivar = 0; ivar < vars.size(); ++ivar ){
    std::cout << "Processing var: " << vars[ivar].varexp << std::endl;
    c1->Clear();

    switch ( trees.size() ){
    case 1:
      break;
    case 2: 
      c1->Divide(2,1);
      break;
    case 3: 
    case 4: 
      c1->Divide(2,2);
      break;
    case 5: 
    case 6: 
      c1->Divide(3,2);
    default:
      c1->Divide(3,3);
    }
  
    for ( unsigned int i = 0; i < trees.size(); ++i ){
      c1->cd(i+1);
      char hname[1024];
      sprintf(hname,"h%d",i+10*ivar);
      char htitle[1024];
      sprintf(htitle,"[%s] %s",trees[i]->GetTitle(),vars[ivar].title);
      if ( vars[ivar].varexp2 == 0 ) {
	TH1F* h = new TH1F(hname,htitle,vars[ivar].nbins,vars[ivar].varmin,vars[ivar].varmax);
	h->GetXaxis()->SetTitle(vars[ivar].varexp);
	h->SetFillColor(kBlue);
	char drawCommand[1024];
	sprintf(drawCommand,"%s>>%s",vars[ivar].varexp,hname);
	trees[i]->Draw(drawCommand, vars[ivar].selection);
      } else {
	if ( vars[ivar].nbins2 == 0 ){
	  TProfile* p = new TProfile(hname,htitle,vars[ivar].nbins,vars[ivar].varmin,vars[ivar].varmax,
				     vars[ivar].varmin2,vars[ivar].varmax2);
	  p->GetXaxis()->SetTitle(vars[ivar].varexp2);
	  p->GetYaxis()->SetTitle(vars[ivar].varexp);
	  p->SetLineColor(kBlue);
	  p->SetMinimum(vars[ivar].varmin2);
	  p->SetMaximum(vars[ivar].varmax2);
	  char drawCommand[1024];
	  sprintf(drawCommand,"%s:%s>>%s",vars[ivar].varexp,vars[ivar].varexp2,hname);
	  trees[i]->Draw(drawCommand, vars[ivar].selection);
	}
      }
    }
    if ( ivar == 0 )
      c1->Print((psname+"(").c_str(),"ps");
    else
      c1->Print(psname.c_str(),"ps");
  }
  c1->Print((psname+"]").c_str(),"ps");
}

void drawAll( const char* ntupleSet1, 
	      const char* ntupleSet2 = 0, 
	      const char* ntupleSet3 = 0, 
	      const char* ntupleSet4 = 0, 
	      const char* ntupleSet5 = 0, 
	      const char* ntupleSet6 = 0, 
	      const char* ntupleSet7 = 0, 
	      const char* ntupleSet8 = 0, 
	      const char* ntupleSet9 = 0)
{
  TChain* chain1 = new TChain("Events");
  chain1->Add( ntupleSet1 );
  if ( ntupleSet2 == 0 ) {
    drawAll(chain1);
    return;
  }
  
  TChain* chain2 = new TChain("Events");
  chain2->Add( ntupleSet2 );
  if ( ntupleSet3 == 0 ) {
    drawAll(chain1,chain2);
    return;
  }
  
  TChain* chain3 = new TChain("Events");
  chain3->Add( ntupleSet3 );
  if ( ntupleSet4 == 0 ) {
    drawAll(chain1,chain2,chain3);
    return;
  }
  
  TChain* chain4 = new TChain("Events");
  chain4->Add( ntupleSet4 );
  if ( ntupleSet5 == 0 ) {
    drawAll(chain1,chain2,chain3,chain4);
    return;
  }
  
  TChain* chain5 = new TChain("Events");
  chain5->Add( ntupleSet5 );
  if ( ntupleSet6 == 0 ) {
    drawAll(chain1,chain2,chain3,chain4,chain5);
    return;
  }
  
  TChain* chain6 = new TChain("Events");
  chain6->Add( ntupleSet6 );
  if ( ntupleSet7 == 0 ) {
    drawAll(chain1,chain2,chain3,chain4,chain5,chain6);
    return;
  }
  
  TChain* chain7 = new TChain("Events");
  chain7->Add( ntupleSet7 );
  if ( ntupleSet8 == 0 ) {
    drawAll(chain1,chain2,chain3,chain4,chain5,chain6,chain7);
    return;
  }
  
  TChain* chain8 = new TChain("Events");
  chain8->Add( ntupleSet8 );
  if ( ntupleSet9 == 0 ) {
    drawAll(chain1,chain2,chain3,chain4,chain5,chain6,chain7,chain8);
    return;
  }
  
  TChain* chain9 = new TChain("Events");
  chain9->Add( ntupleSet9 );
  drawAll(chain1,chain2,chain3,chain4,chain5,chain6,chain7,chain8,chain9);
} 
