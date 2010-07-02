#include <iostream>
#include <vector>

#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"

//#include "CMS2.h"
#include "CMS2.cc"

//#include "../CORE/selections.cc"
//#include "../CORE/utilities.cc"
//#include "../Tools/tools.cc"
#include "../CORE/electronSelections.cc"

using namespace tas;
//CMS2 cms2;

TH1F* book1DHist(const char* name, const char* title, unsigned int nbins, float low, float high, const char* xtitle, const char* ytitle, int color) {
  // return histogram instance with called Sumw2
  TH1F *hist = new TH1F(name,title,nbins,low,high);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->Sumw2();
  hist->SetFillColor(color);
  hist->SetLineColor(color);
   
  return hist;   
}

TH2F* book2DHist(const char* name, const char* title, 
                 unsigned int xnbins, float xlow, float xhigh, 
                 unsigned int ynbins, float ylow, float yhigh, 
                 const char* xtitle, const char* ytitle, const char* ztitle, 
                 int color) {
  // return histogram instance with called Sumw2
  TH2F *hist = new TH2F(name,title,xnbins,xlow,xhigh,ynbins,ylow,yhigh);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->SetZTitle(ztitle);
  hist->Sumw2(); 
  hist->SetFillColor(color);
  hist->SetStats(kFALSE);
  //hist->SetLineColor(color);
   
  return hist;   
}

int ScanChain( TChain* chain) {

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=chain->GetEntries();;
  unsigned int nEventsTotal = 0;

  // book histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  unsigned int nBinsPt 	= 55;
  //unsigned int nBinsPt 	= 5;
  float lowBinsPt 	= 0.;
  float highBinsPt 	= 110.;

  unsigned int nBinsEta = 52;
  //unsigned int nBinsEta = 5;
  float lowBinsEta      = -2.6;
  float highBinsEta     =  2.6;
  
  TH1F *els_pt_sim 			= book1DHist("els_pt_sim", 
						     "true Electron: p_{T}^{true}",
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_sim 			= book1DHist("els_eta_sim", 
						     "true Electron: #eta^{true}",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_recosim 			= book1DHist("els_pt_recosim", 
						     "true Electron: p_{T}^{true}",
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_recosim 		= book1DHist("els_eta_recosim", 
						     "true Electron: #eta^{true}",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_reco 			= book1DHist("els_pt_reco", 
						     "reco Electron: p_{T}^{true}",
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_reco 			= book1DHist("els_eta_reco", 
						     "reco Electron: #eta^{true}",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_recosim_incorCharge 	= book1DHist("els_pt_recosim_incorCharge", 
						     "true Electron with incorrect reconstructed Charge: p_{T}^{true}",
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_recosim_incorCharge 	= book1DHist("els_eta_recosim_incorCharge", 
						     "true Electron with incorrect reconstructed Charge: #eta^{true}",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_recosim_corCharge 	= book1DHist("els_pt_recosim_corCharge", 
						     "true Electron with correct reconstructed Charge: p_{T}^{true}",
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_recosim_corCharge 	= book1DHist("els_eta_recosim_corCharge", 
						     "true Electron with correct reconstructed Charge: #eta^{true}",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_reco_corCharge 		= book1DHist("els_pt_reco_corCharge", 
						     "reco Electron with correct reconstructed Charge: p_{T}^{true}",
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_reco_corCharge 		= book1DHist("els_eta_reco_corCharge", 
						     "reco Electron with correct reconstructed Charge: #eta^{true}",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);

  TH1F *els_trkId         		= book1DHist("els_trkId", 
						     "Track ID of associated to reco Electron",
						     30,
						     -1050.,
						     450.,
						     "Track ID",
						     "Electrons",2);
  // 2D histograms 
  TH2F *els_2d_eta_Pt_corCharge		= book2DHist("els_2d_eta_Pt_corCharge",
						     "2D histogram of eta vs Pt corCharge",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "#eta^{true}",
						     "p_{T}^{true} [GeV]",
						     " ",2);

  TH2F *els_2d_eta_Pt_incorCharge	= book2DHist("els_2d_eta_Pt_incorCharge",
						     "2D histogram of eta vs Pt incorCharge",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "#eta^{true}",
						     "p_{T}^{true} [GeV]",
						     " ",2);

  TH2F *els_2d_eta_Pt_ratio 		= book2DHist("els_2d_eta_Pt_ratio",
						     "2D histogram of eta vs Pt ratio",
						     nBinsEta,
						     lowBinsEta,
						     highBinsEta,
						     nBinsPt,
						     lowBinsPt,
						     highBinsPt,
						     "#eta^{true}",
						     "p_{T}^{true} [GeV]",
						     //"MisID rate",2);
						     " ",2);

  int num_beforecut=0;
  int num_ELEID_CAND02=0;
  int num_ELEID_EXTRA=0;
  int num_ELENOTCONV_DISTDCOT002=0;
  int num_ELENOTCONV_HITPATTERN=0;
  int num_ELECHARGE_NOFLIP=0;

  // individual cuts in ELEID_CAND02
  int num_ELEID_CAND02_ELESEED_ECAL=0;
  int num_ELEID_CAND02_ELEETA_250=0;
  int num_ELEID_CAND02_ELENOMUON_010=0;
  int num_ELEID_CAND02_ELEID_CAND02=0;
  int num_ELEID_CAND02_ELEISO_REL010=0;
  int num_ELEID_CAND02_ELENOTCONV_DISTDCOT002=0;
  int num_ELEID_CAND02_ELEIP_200=0;


  int num_cand02=0;
  int num_extra=0;
  int num_hitpattern=0;
  int num_partnertrack=0;
  int num_chargeflip=0;
  int num_cand02_ecal=0;
  int num_cand02_eta250=0;
  int num_cand02_noMuon=0;
  int num_cand02_cand02=0;
  int num_cand02_impact=0;
  int num_cand02_relsusy010=0;


  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    // nEvents = 100;
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;
     
      if ( nEventsTotal%10000 == 0 ) {
	std::cout << "Event: " << nEventsTotal << endl;
      }

      // loop over true electrons
      for ( unsigned int els = 0;
	    els < genps_p4().size();
	    ++els ) {

	// pt cut
	if( (genps_p4()[els].pt()) < 10 ) continue;

	// check that electron is final state electron
	if ( genps_status()[els] != 1 ) continue;
	
	// check for true electron
	if ( TMath::Abs(genps_id()[els]) != 11 ) continue;

	// fill true histrograms
	els_pt_sim->Fill(genps_p4()[els].pt());
	els_eta_sim->Fill(genps_p4()[els].eta());

      }

// down to here, no problem.

      // loop over reco electrons
      for (unsigned int els = 0; 
           els < els_p4().size(); 
	   ++els) {

	//
	// cuts
	//

	// pt
	if( (els_p4().at(els).Pt()) < 10 ) continue;

	//if (!pass_electronSelection(els, electronSelection_oldss)) 
        //   continue;

        // ttbar_v1
        //if (!pass_electronSelection(els, electronSelection_ttbarV1))
        // continue;

	//if (!pass_electronSelection(els, electronSelection_ELECHARGE_NOFLIP))
	//continue; 

	// 3 flip
	//if (!pass_electronSelection(els, electronSelection_ELECHARGE_NOFLIP3AGREE))
	//continue; 

	//if (pass_electronSelection(els, electronSelection_ELECHARGE_NOFLIP))
        // ttbar_v1_chargeflip
        //if (!pass_electronSelection(els, electronSelection_ttbarV1_chargeflip))
        //   continue;



	// new SS
        //if (!pass_electronSelection(els, electronSelection_ttbarV1_noiso)) 
        //   continue;
        //if (!pass_electronSelection(els, electronSelection_ttbarV1_iso))
        //   continue;

	// new cuts debug
///*

	num_beforecut++;

	//if (!pass_electronSelection(els, electronSelection_ELEID_CAND02))
	//continue; 
	//if (pass_electronSelection(els, electronSelection_ELEID_CAND02))
	// in 1.56, mimic 1.26
	if (!pass_electronSelection(els, electronSelection_ELEID_CAND02_old))
	continue;
	//num_ELEID_CAND02++;

	if (!pass_electronSelection(els, electronSelection_ELEID_EXTRA))
	continue; 
	//if (pass_electronSelection(els, electronSelection_ELEID_EXTRA))
	//num_ELEID_EXTRA++;

	if (!pass_electronSelection(els, electronSelection_ELENOTCONV_DISTDCOT002))
	continue; 
	//if (pass_electronSelection(els, electronSelection_ELENOTCONV_DISTDCOT002))
	//num_ELENOTCONV_DISTDCOT002++;

	if (!pass_electronSelection(els, electronSelection_ELENOTCONV_HITPATTERN))
	continue; 
	//if (pass_electronSelection(els, electronSelection_ELENOTCONV_HITPATTERN))
	//num_ELENOTCONV_HITPATTERN++;

	//if (!pass_electronSelection(els, electronSelection_ELECHARGE_NOFLIP))
	//continue; 
	//if (pass_electronSelection(els, electronSelection_ELECHARGE_NOFLIP))
	//num_ELECHARGE_NOFLIP++;

	//if (!pass_electronSelection(els, electronSelection_ELECHARGE_NOFLIP3AGREE))
	//continue; 

	// print out numbers
	if (pass_electronSelection(els, electronSelection_ELEID_CAND02_old_ELESEED_ECAL))
	num_ELEID_CAND02_ELESEED_ECAL++;

	if (pass_electronSelection(els, electronSelection_ELEID_CAND02_old_ELEETA_250))
	num_ELEID_CAND02_ELEETA_250++;

	if (pass_electronSelection(els, electronSelection_ELEID_CAND02_old_ELENOMUON_010))
	num_ELEID_CAND02_ELENOMUON_010++;

	if (pass_electronSelection(els, electronSelection_ELEID_CAND02_old_ELEID_CAND02))
	num_ELEID_CAND02_ELEID_CAND02++;

	if (pass_electronSelection(els, electronSelection_ELEID_CAND02_old_ELEISO_REL010))
	num_ELEID_CAND02_ELEISO_REL010++;

	if (pass_electronSelection(els, electronSelection_ELEID_CAND02_old_ELENOTCONV_DISTDCOT002))
	num_ELEID_CAND02_ELENOTCONV_DISTDCOT002++;

	if (pass_electronSelection(els, electronSelection_ELEID_CAND02_old_ELEIP_200))
	num_ELEID_CAND02_ELEIP_200++;
//*/


	// double check Derek's work
	//if (!pass_electronSelection(els, electronSelection_ss))
	//continue; 
	 
	//if (!pass_electronSelection(els, electronSelection_ELECHARGE_NOFLIP))
	//continue; 

	//if (!pass_electronSelection(els, electronSelection_ELECHARGE_NOFLIP3AGREE))
	//continue; 


	// old cuts std::cout<<"debug"<<endl;

/*
        // eleId_v3(1.26) 
	num_beforecut++;

	if (!electronSelection_cand02(els)) continue;
	//if (electronSelection_cand02(els)) 
	//  num_cand02++;

	if (!electronId_extra(els) ) continue;
	//if (electronId_extra(els) ) 
	//  num_extra++;

	if(isFromConversionPartnerTrack(els)) continue;
	//if(!isFromConversionPartnerTrack(els)) 
	//  num_hitpattern++;

        // need to be commented out for comparison with 336p4
	if(isFromConversionHitPattern(els)) continue; 
	//if(!isFromConversionHitPattern(els))  
	//  num_partnertrack++;

	//if ( isChargeFlip(els) ) continue; 
	//if(!isChargeFlip(els) ) 
	//  num_chargeflip++;


	// debugg
	if (electronSelection_cand02_ecal(els)) 
	  num_cand02_ecal++;

	if (electronSelection_cand02_eta250(els)) 
	  num_cand02_eta250++;

	if (electronSelection_cand02_noMuon(els)) 
	  num_cand02_noMuon++;

	if (electronSelection_cand02_cand02(els)) 
	  num_cand02_cand02++;

	if (electronSelection_cand02_impact(els)) 
	  num_cand02_impact++;

	if (electronSelection_cand02_relsusy010(els)) 
	  num_cand02_relsusy010++;

*/


	// check how many electrons don't have an associated track
	els_trkId->Fill(els_trkidx().at(els));
	// tmp charge variable
	//double charge = els_charge().at(els); //
	double charge = els_trk_charge().at(els);
	
	// fill reco
	els_pt_reco->Fill(els_p4().at(els).Pt());
	els_eta_reco->Fill(els_p4().at(els).eta());

	// fill reco_corCharge
	if ( (charge == -1 && els_mc_id().at(els) == 11) || (charge == 1 && els_mc_id().at(els) == -11) ) {
	  els_pt_reco_corCharge->Fill(els_p4().at(els).Pt());
	  els_eta_reco_corCharge->Fill(els_p4().at(els).eta());
	}

	// exclude reco which has no true electron match
	if ( abs(els_mc_id()[els]) != 11 ) continue;

	// fill recosim
	els_pt_recosim->Fill(els_mc_p4().at(els).Pt());
	els_eta_recosim->Fill(els_mc_p4().at(els).eta());

	// correct charge identified
	if ( (charge == -1 && els_mc_id().at(els) == 11) || (charge == 1 && els_mc_id().at(els) == -11) ) {

	  els_pt_recosim_corCharge->Fill(els_mc_p4().at(els).Pt());
	  els_eta_recosim_corCharge->Fill(els_mc_p4().at(els).eta());
	  els_2d_eta_Pt_corCharge->Fill(els_mc_p4().at(els).eta(), els_mc_p4().at(els).Pt()); //2d


	  // incorrect charge identified
	} else {

	  els_pt_recosim_incorCharge->Fill(els_mc_p4().at(els).Pt());
	  els_eta_recosim_incorCharge->Fill(els_mc_p4().at(els).eta());
	  els_2d_eta_Pt_incorCharge->Fill(els_mc_p4().at(els).eta(), els_mc_p4().at(els).Pt()); //2d
	}

      }
      
    }
  }

  // 2d ratio
  Float_t deno=0, num=0;
  for(Int_t i=0; i<nBinsEta; i++)
  {
    for(Int_t j=0; j<nBinsPt; j++)
    {
       deno = els_2d_eta_Pt_corCharge->GetBinContent(i+1, j+1)	
            +els_2d_eta_Pt_incorCharge->GetBinContent(i+1, j+1);
       num = els_2d_eta_Pt_incorCharge->GetBinContent(i+1, j+1);

       els_2d_eta_Pt_ratio->SetBinContent(i+1, j+1, num/deno);
    }
  } 

/*
  cout << " ------------------------------------------------------------- " << endl;
  cout << "num_beforecut             : " << num_beforecut << endl;
  cout << "num_ELEID_CAND02          : " << num_ELEID_CAND02 << endl;
  cout << "num_ELEID_EXTRA           : " << num_ELEID_EXTRA << endl;
  cout << "num_ELENOTCONV_DISTDCOT002: " << num_ELENOTCONV_DISTDCOT002 << endl;
  cout << "num_ELENOTCONV_HITPATTERN : " << num_ELENOTCONV_HITPATTERN  << endl;
  cout << "num_ELECHARGE_NOFLIP      : " << num_ELECHARGE_NOFLIP << endl;
  cout << "num_ELEID_CAND02_ELESEED_ECAL          : " << num_ELEID_CAND02_ELESEED_ECAL << endl;
  cout << "num_ELEID_CAND02_ELEETA_250            : " << num_ELEID_CAND02_ELEETA_250 << endl;
  cout << "num_ELEID_CAND02_ELENOMUON_010         : " << num_ELEID_CAND02_ELENOMUON_010 << endl;
  cout << "num_ELEID_CAND02_ELEID_CAND02          : " << num_ELEID_CAND02_ELEID_CAND02 << endl;
  cout << "num_ELEID_CAND02_ELEISO_REL010         : " << num_ELEID_CAND02_ELEISO_REL010 << endl;
  cout << "num_ELEID_CAND02_ELENOTCONV_DISTDCOT002: " 
       << num_ELEID_CAND02_ELENOTCONV_DISTDCOT002 << endl;
  cout << "num_ELEID_CAND02_ELEIP_200             : " << num_ELEID_CAND02_ELEIP_200 << endl; 
  cout << " ------------------------------------------------------------- " << endl;
*/
/*
  cout << " ----------------------------------------------- " << endl;
  cout << "num_beforecut          : " << num_beforecut << endl;
  cout << "num_cand02             : " << num_cand02 << endl; 
  cout << "num_extra              : " << num_extra << endl; 
  cout << "num_hitpattern         : " << num_hitpattern << endl; 
  cout << "num_partnertrack       : " << num_partnertrack << endl; 
  cout << "num_chargeflip         : " << num_chargeflip << endl; 
  cout << "num_cand02_ecal          : " << num_cand02_ecal << endl;
  cout << "num_cand02_eta250        : " << num_cand02_eta250  << endl;
  cout << "num_cand02_noMuon        : " << num_cand02_noMuon  << endl;
  cout << "num_cand02_cand02        : " << num_cand02_cand02  << endl; 
  cout << "num_cand02_impact        : " << num_cand02_impact  << endl;
  cout << "num_cand02_relsusy010    : " << num_cand02_relsusy010 << endl;
  cout << " ----------------------------------------------- " << endl;
*/





  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  return 0;
}
