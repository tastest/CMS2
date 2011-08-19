#include <iostream>
#include <vector>

#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

//#include "CMS2.h"
#include "CMS2.cc"

//#include "../CORE/selections.cc"
//#include "../CORE/utilities.cc"
//#include "../Tools/tools.cc"
#include "../../CORE/electronSelections.cc"
//#include "../Tools/goodrun.cc"
#include "../../fakeRate/goodrun.cc"


using namespace tas;
//CMS2 cms2;




// histogram definition
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
                 int color) 
{
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

int ScanChain( TChain* chain) 
{
 	set_goodrun_file("./goodruns.txt");
	TObjArray *listOfFiles = chain->GetListOfFiles();

	unsigned int nEventsChain=chain->GetEntries();;
	unsigned int nEventsTotal = 0;

	// book histograms
	TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
	rootdir->cd();

	unsigned int nBinsPt 	= 55;
	float lowBinsPt 	= 0.;
	float highBinsPt 	= 110.;

	unsigned int nBinsEta = 52;
	float lowBinsEta      = -2.6;
	float highBinsEta     =  2.6;
  
	// 1D histograms example
	TH1F *evt_pfmet_noflip			= book1DHist("evt_pfmet_noflip", 
							     "",
							     20,
							     0,
							     100,
							     "pfMET [GeV]",
							     "",2);

	TH1F *evt_pfmet_3exist			= book1DHist("evt_pfmet_3exist", 
							     "",
							     20,
							     0,
							     100,
							     "pfMET [GeV]",
							     "",2);

	TH1F *evt_pfmet_3same			= book1DHist("evt_pfmet_3same", 
							     "",
							     20,
							     0,
							     100,
							     "pfMET [GeV]",
							     "",2);

	TH1F *els_MT_noflip			= book1DHist("els_MT_nofilp", 
							     "",
							     20,
							     0,
							     150,
							     "MT [GeV]",
							     "",2);

	TH1F *els_MT_3exist			= book1DHist("els_MT_3exist", 
							     "",
							     20,
							     0,
							     150,
							     "MT [GeV]",
							     "",2);

	TH1F *els_MT_3same			= book1DHist("els_MT_3same", 
							     "",
							     20,
							     0,
							     150,
							     "MT [GeV]",
							     "",2);
	// 2D histograms
	TH2F *els_2d_eta_Pt_3exist		= book2DHist("els_2d_eta_Pt_3exist",
							     " ",
							     nBinsEta,
							     lowBinsEta,
							     highBinsEta,
							     nBinsPt,
							     lowBinsPt,
							     highBinsPt,
							     "#eta",
							     "p_{T} [GeV]",
							     " ",2);

	TH2F *els_2d_eta_Pt_3same		= book2DHist("els_2d_eta_Pt_3same",
							     " ",
							     nBinsEta,
							     lowBinsEta,
							     highBinsEta,
							     nBinsPt,
							     lowBinsPt,
							     highBinsPt,
							     "#eta",
							     "p_{T} [GeV]",
							     " ",2);


	// count # events
	unsigned int num_event=0;
	unsigned int num_evt_pfmet=0;
	unsigned int num_ss=0;
	unsigned int num_pt=0;
	unsigned int num_spike=0;
	unsigned int num_elepass=0;

	unsigned int num_ss_cand=0;
	unsigned int num_pt_cand=0;
	unsigned int num_spike_cand=0;

	unsigned int num_noflip=0;
	unsigned int num_3charge_exist=0;
	unsigned int num_3charge_same=0;

	//
	// file loop
	//
	TIter fileIter(listOfFiles);
	TFile *currentFile = 0;
	while ( currentFile = (TFile*)fileIter.Next() ) 
	{
		TFile f(currentFile->GetTitle());
		TTree *tree = (TTree*)f.Get("Events");
		cms2.Init(tree);

	
		//    
		// Event Loop
		//   
		unsigned int nEvents = tree->GetEntries();
		for( unsigned int event = 0; event < nEvents; ++event) 
		{
			cms2.GetEntry(event);

			// goodrun 
			//if( goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) )
			//{ 
				//cout << " evt_run : " << cms2.evt_run()  
				//    << " evt_lumiBlock : " << cms2.evt_lumiBlock() << endl;
				//cout << " ** goodrun ** " << endl; 
			//	continue;
			//}

			++nEventsTotal;
   
			if ( nEventsTotal%100000 == 0 ) 
			{
				std::cout << "Event: " << nEventsTotal << endl;
			}

			// pfMET cut
			//if( cms2.evt_pfmet() < 25 || cms2.evt_pfmet() > 80 ) 
			//continue;
			num_evt_pfmet++;	


			//
			// loop over gsf electrons
			//
			unsigned int  num_pass_inevt=0;
			unsigned int  index_pass_cand=0;
			for (unsigned int els = 0; els < cms2.els_p4().size(); ++els) 
			{

				//
				// cuts
				//
				
				// electron ID : SS
			        if(!pass_electronSelection(els, electronSelection_ss, true)) 
				continue;	

				num_ss++;

				// pT > 20 GeV
				if (cms2.els_p4()[els].pt() < 20) 
				continue;	
				
				num_pt++;

				// spike rejection
        			if(!pass_electronSelection(els, electronSelection_spike_rejection)) 
				continue;			

				num_spike++;
				num_pass_inevt++;
				index_pass_cand=els;

			}
			
			// choose an event with one eletron candidate
			if(num_pass_inevt != 1) 
			continue;

			// calculate MT
			// MT^2 = ET^2 - pT^2
			float MT2;
			TLorentzVector p4_ele(cms2.els_p4()[index_pass_cand].px(),
					      cms2.els_p4()[index_pass_cand].py(),
					      cms2.els_p4()[index_pass_cand].pz(),
					      cms2.els_p4()[index_pass_cand].e());
			TLorentzVector pfMET(cms2.evt_pfmet()*cos(cms2.evt_pfmetPhi()),
					      cms2.evt_pfmet()*sin(cms2.evt_pfmetPhi()),
					      0, cms2.evt_pfmet());
			MT2 = (pfMET.Pt()+p4_ele.Et())*(pfMET.Pt()+p4_ele.Et())
			     -(pfMET.Px()+p4_ele.Px())*(pfMET.Px()+p4_ele.Px())	
			     -(pfMET.Py()+p4_ele.Py())*(pfMET.Py()+p4_ele.Py());	
			//cout << " MT : " << sqrt(MT2) << endl;

			evt_pfmet_noflip->Fill(cms2.evt_pfmet());
			els_MT_noflip->Fill(sqrt(MT2));
			num_noflip++;

			
			// track id exists
			if (cms2.els_trkidx().at(index_pass_cand) < 0) 
			continue;		
		
			// 3 charges exist						
			evt_pfmet_3exist->Fill(cms2.evt_pfmet());
			els_MT_3exist->Fill(sqrt(MT2));
			num_3charge_exist++;
			els_2d_eta_Pt_3exist->Fill(cms2.els_p4()[index_pass_cand].eta(),
						   cms2.els_p4()[index_pass_cand].pt()); 

			// 3 charges same						
			if(!(cms2.els_trk_charge().at(index_pass_cand)
				== cms2.trks_charge().at(cms2.els_trkidx().at(index_pass_cand))
			   && cms2.trks_charge().at(cms2.els_trkidx().at(index_pass_cand)) 
		                == cms2.els_sccharge().at(index_pass_cand))) 
			continue;		

			evt_pfmet_3same->Fill(cms2.evt_pfmet());
			els_MT_3same->Fill(sqrt(MT2));
			num_3charge_same++;
			els_2d_eta_Pt_3same->Fill(cms2.els_p4()[index_pass_cand].eta(),
					  cms2.els_p4()[index_pass_cand].pt());
		

			//cout << cms2.evt_run() << "  " << cms2.evt_lumiBlock() << endl;
			//cout << " --------------------------- " << endl;


		} // event loop

	} // file loop

	cout << " -------------------------------------------------- " << endl;
	cout << "nEventsTotal           : " << nEventsTotal << endl;
	cout << "num_evt_pfmet          : " << num_evt_pfmet << endl;
	cout << "num_ss                 : " << num_ss << endl;
	cout << "num_pt                 : " << num_pt << endl;
	cout << "num_spike(pass)        : " << num_spike << endl;
	cout << " -------------------------------------------------- " << endl;
	cout << "num_noflip             : " << num_noflip << endl;
	cout << "num_3charge_exist      : " << num_3charge_exist << endl;
	cout << "num_3charge_same       : " << num_3charge_same << endl;
	cout << " -------------------------------------------------- " << endl;

	if ( nEventsChain != nEventsTotal ) 
	{
		std::cout << "ERROR: number of events from files is "
                          << "not equal to total number of events " << std::endl;
	}


	return 0;
}
