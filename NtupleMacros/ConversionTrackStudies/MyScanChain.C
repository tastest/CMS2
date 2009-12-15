/* Usage:
   root [0] .L ScanChain.C++
   root [1] TFile *_file0 = TFile::Open("merged_ntuple.root")
   root [2] TChain *chain = new TChain("Events")
   root [3] chain->Add("merged_ntuple.root")

   There are several places where one may create CMS2 cms2
   It can be done here (in a doAll.C script), i.e.:

   root [4] CMS2 cms2 

   It can be done in the source as is done below, or it can be
   ascertained by including CORE/CMS2.cc as is commented out
   below.  They are all the same, and everything will work so
   long as it is created somewhere globally.

   root [5] ScanChain(chain)
 */
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "CMS2.h"
#include "myselections.cc"

CMS2 cms2;

/*
#include "CORE/CMS2.cc"
#include "CORE/selections.cc"
#include "CORE/utilities.cc"
 */

using namespace tas;

int ScanChain(bool isData, TChain* chain, int nEvents = -1, std::string skimFilePrefix="") {

	TObjArray *listOfFiles = chain->GetListOfFiles();

	unsigned int nEventsChain=0;
	if(nEvents==-1) 
		nEvents = chain->GetEntries();
	nEventsChain = nEvents;
	unsigned int nEventsTotal = 0;
	TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

	TH1F *h1_pseudo_ecalEta_ee = new TH1F("h1_pseudo_ecalEta_ee", "ecalEta;ecalEta", 50, -5, 5);
	TH1F *h1_pseudo_ecalPhi_ee = new TH1F("h1_pseudo_ecalPhi_ee", "ecalPhi;ecalPhi", 50, -5, 5);

	TH1F *h1_pseudo_ecalEta_eb = new TH1F("h1_pseudo_ecalEta_eb", "ecalEta;ecalEta", 50, -5, 5);
	TH1F *h1_pseudo_ecalPhi_eb = new TH1F("h1_pseudo_ecalPhi_eb", "ecalPhi;ecalPhi", 50, -5, 5);

	Int_t ecalIsoN = 90;
	Float_t ecalIsoMin = -1.0;
	Float_t ecalIsoMax = 2.0;

        Int_t hcalIsoN = 70;
        Float_t hcalIsoMin = -0.5;
        Float_t hcalIsoMax = 3.0;

        Int_t ecalIsoHitN = 100;
        Float_t ecalIsoHitMin = -0.5;
        Float_t ecalIsoHitMax = 0.5;

	Int_t tkIsoN = 40;
	Float_t tkIsoMin = -1.0;
	Float_t tkIsoMax = 3.0;

	TH1F *h1_pseudo_ecalIso03_E015_eb = new TH1F("h1_pseudo_ecalIso03_E015_eb", "ecalIso03;ecalIso03_E015", ecalIsoN, ecalIsoMin, ecalIsoMax);
	TH1F *h1_pseudo_ecalIso03_E015_ee = new TH1F("h1_pseudo_ecalIso03_E015_ee", "ecalIso03;ecalIso03_E015", ecalIsoN, ecalIsoMin, ecalIsoMax);

        TH1F *h1_pseudo_ecalIso03_E010_eb = new TH1F("h1_pseudo_ecalIso03_E010_eb", "ecalIso03;ecalIso03_E010", ecalIsoN, ecalIsoMin, ecalIsoMax);
        TH1F *h1_pseudo_ecalIso03_E010_ee = new TH1F("h1_pseudo_ecalIso03_E010_ee", "ecalIso03;ecalIso03_E010", ecalIsoN, ecalIsoMin, ecalIsoMax);

	TH1F *h1_pseudo_ecalIso03_recHitEt_eb = new TH1F("h1_pseudo_ecalIso03_recHitEt_eb", "ecalIso03;ecalIso03_recHitEt", ecalIsoHitN, ecalIsoHitMin, ecalIsoHitMax);
	TH1F *h1_pseudo_ecalIso03_recHitEt_ee = new TH1F("h1_pseudo_ecalIso03_recHitEt_ee", "ecalIso03;ecalIso03_recHitEt", ecalIsoHitN, ecalIsoHitMin, ecalIsoHitMax);
        TH1F *h1_pseudo_ecalIso03_recHitE_eb = new TH1F("h1_pseudo_ecalIso03_recHitE_eb", "ecalIso03;ecalIso03_recHitE", ecalIsoHitN, ecalIsoHitMin, ecalIsoHitMax);
        TH1F *h1_pseudo_ecalIso03_recHitE_ee = new TH1F("h1_pseudo_ecalIso03_recHitE_ee", "ecalIso03;ecalIso03_recHitE", ecalIsoHitN, ecalIsoHitMin, ecalIsoHitMax);

	TH1F *h1_pseudo_ecalIso03_recHitN_eb = new TH1F("h1_pseudo_ecalIso03_recHitN_eb", "ecalIso03;ecalIso03_recHitN", 50, 0, 50);
	TH1F *h1_pseudo_ecalIso03_recHitN_ee = new TH1F("h1_pseudo_ecalIso03_recHitN_ee", "ecalIso03;ecalIso03_recHitN", 50, 0, 50);

	TH1F *h1_pseudo_ecalIso03_eb = new TH1F("h1_pseudo_ecalIso03_eb", "ecalIso03;ecalIso03", ecalIsoN, ecalIsoMin, ecalIsoMax);
	TH1F *h1_pseudo_ecalIso03_ebp = new TH1F("h1_pseudo_ecalIso03_ebp", "ecalIso03;ecalIso03", ecalIsoN, ecalIsoMin, ecalIsoMax);
	TH1F *h1_pseudo_ecalIso03_ebm = new TH1F("h1_pseudo_ecalIso03_ebm", "ecalIso03;ecalIso03", ecalIsoN, ecalIsoMin, ecalIsoMax);

	TH1F *h1_pseudo_ecalIso03_qual_eb = new TH1F("h1_pseudo_ecalIso03_qual_eb", "ecalIso03;ecalIso03", ecalIsoN, ecalIsoMin, ecalIsoMax);
	TH1F *h1_pseudo_ecalIso03_fakeSR_eb = new TH1F("h1_pseudo_ecalIso03_fakeSR_eb", "ecalIso03;ecalIso03", ecalIsoN, ecalIsoMin, ecalIsoMax);
	TH1F *h1_pseudo_hcalD1Iso03_eb = new TH1F("h1_pseudo_hcalD1Iso03_eb", "hcalD1Iso03;hcalD1Iso03", hcalIsoN, hcalIsoMin, hcalIsoMax);
	TH1F *h1_pseudo_hcalD2Iso03_eb = new TH1F("h1_pseudo_hcalD2Iso03_eb", "hcalD2Iso03;hcalD2Iso03", hcalIsoN, hcalIsoMin, hcalIsoMax);
	TH1F *h1_pseudo_tkIso03_eb = new TH1F("h1_pseudo_tkIso03_eb", "tkIso03;tkIso03", tkIsoN, tkIsoMin, tkIsoMax);

	TH1F *h1_pseudo_ecalIso03_ee = new TH1F("h1_pseudo_ecalIso03_ee", "ecalIso03;ecalIso03", ecalIsoN, ecalIsoMin, ecalIsoMax);
	TH1F *h1_pseudo_ecalIso03_qual_ee = new TH1F("h1_pseudo_ecalIso03_qual_ee", "ecalIso03;ecalIso03", ecalIsoN, ecalIsoMin, ecalIsoMax);
	TH1F *h1_pseudo_ecalIso03_eep = new TH1F("h1_pseudo_ecalIso03_eep", "ecalIso03;ecalIso03", ecalIsoN, ecalIsoMin, ecalIsoMax);
	TH1F *h1_pseudo_ecalIso03_eem = new TH1F("h1_pseudo_ecalIso03_eem", "ecalIso03;ecalIso03", ecalIsoN, ecalIsoMin, ecalIsoMax);

	TH1F *h1_pseudo_ecalIso03_fakeSR_ee = new TH1F("h1_pseudo_ecalIso03_fakeSR_ee", "ecalIso03;ecalIso03", ecalIsoN, ecalIsoMin, ecalIsoMax);
	TH1F *h1_pseudo_hcalD1Iso03_ee = new TH1F("h1_pseudo_hcalD1Iso03_ee", "hcalD1Iso03;hcalD1Iso03", hcalIsoN, hcalIsoMin, hcalIsoMax);
	TH1F *h1_pseudo_hcalD2Iso03_ee = new TH1F("h1_pseudo_hcalD2Iso03_ee", "hcalD2Iso03;hcalD2Iso03", hcalIsoN, hcalIsoMin, hcalIsoMax);
	TH1F *h1_pseudo_tkIso03_ee = new TH1F("h1_pseudo_tkIso03_ee", "tkIso03;tkIso03", tkIsoN, tkIsoMin, tkIsoMax);

	TH1F *h1_pseudo_dRClosestTower_ee = new TH1F("h1_pseudo_dRClosestTower_ee", "dRClosestTower;dRClosestTower", 50, 0, 5);
	TH1F *h1_pseudo_dRClosestTower_eb = new TH1F("h1_pseudo_dRClosestTower_eb", "dRClosestTower;dRClosestTower", 50, 0, 5);

	TH1F *h1_l1_techbits2_pass = new TH1F("h1_l1_techbits2_pass", "l1_techbits2", 32, -0.5, 31.5);
	TH1F *h1_l1_techbits2_total = new TH1F("h1_l1_techbits2_total", "l1_techbits2", 32, -0.5, 31.5);

	TH1F *h1_l1_techbits2_track_pass = new TH1F("h1_l1_techbits2_track_pass", "l1_techbits2", 32, -0.5, 31.5);
	TH1F *h1_l1_techbits2_track_total = new TH1F("h1_l1_techbits2_track_total", "l1_techbits2", 32, -0.5, 31.5);

	TH1F *h1_l1_techbits2_vtxs_pass = new TH1F("h1_l1_techbits2_vtxs_pass", "l1_techbits2", 32, -0.5, 31.5);
	TH1F *h1_l1_techbits2_vtxs_total = new TH1F("h1_l1_techbits2_vtxs_total", "l1_techbits2", 32, -0.5, 31.5);

	// file loop

	TIter fileIter(listOfFiles);
	TFile *currentFile = 0;
	while ( currentFile = (TFile*)fileIter.Next() ) {
		TFile f(currentFile->GetTitle());
		TTree *tree = (TTree*)f.Get("Events");
		cms2.Init(tree);

		//Event Loop
		unsigned int nEvents = tree->GetEntries();
		for( unsigned int event = 0; event < nEvents; ++event) {
			cms2.GetEntry(event);
			++nEventsTotal;

			// how many good tracks
			//
			int nGoodTracksVz = 0;
			int nHighPurityTracks = 0;
			for (size_t t = 0; t < cms2.trks_ndof().size(); ++t) {

				// count high purity tracks
				if (isTrackQuality(t, (1<<highPurity))) nHighPurityTracks ++;

				// cut on track ndof
				if (cms2.trks_ndof()[t] < 8) continue;
				if (cms2.vtxs_position().size() > 0)
					if (fabs(cms2.vtxs_position()[0].z() - cms2.trks_vertex_p4()[t].z()) < 10.0) continue;
					else if (fabs(cms2.trks_vertex_p4()[t].z()) < 10.0) continue;
					nGoodTracksVz ++;
			}

			// count good vertexs
			//
			int nGoodVertex = 0;
			for (size_t v = 0; v < cms2.vtxs_position().size(); ++v) {

				if (cms2.vtxs_isFake()[v]) continue;
				if (fabs(cms2.vtxs_position()[v].z()) > 10.0
					|| fabs(cms2.vtxs_position()[v].x()) > 0.50
                                        || fabs(cms2.vtxs_position()[v].y()) > 0.50) continue;
				nGoodVertex ++;
			}


			// trigger info
			for (Int_t bin = 1; bin <= 32; ++bin)
			{

				// denominator all
				h1_l1_techbits2_total->Fill(bin, 1);
				if (cms2.l1_techbits2() & 1<<(bin-1)) {
					h1_l1_techbits2_pass->Fill(bin, 1);
				}

				// denominator > 1 good track
				if (nGoodTracksVz > 1) {
					h1_l1_techbits2_track_total->Fill(bin, 1);
					if (cms2.l1_techbits2() & 1<<(bin-1)) {
						h1_l1_techbits2_track_pass->Fill(bin, 1);
					}
				}

				// denominator 1 good vertex
				if (cms2.evt_nvtxs() == 1 && !cms2.vtxs_isFake()[0] ) {
					h1_l1_techbits2_vtxs_total->Fill(bin, 1);
					if (cms2.l1_techbits2() & 1<<(bin-1)) {
						h1_l1_techbits2_vtxs_pass->Fill(bin, 1);
					}					
				}


			}

			// require bit 40 or 41 passed
			//
			if (!(cms2.l1_techbits2() & (1<<8) || cms2.l1_techbits2() & (1<<9))) continue;

			// require bits 36-39 DIDN't pass ???
			//
			if (cms2.l1_techbits2() & (1<<7) || cms2.l1_techbits2() & (1<<6) || 
				cms2.l1_techbits2() & (1<<5) || cms2.l1_techbits2() & (1<<4)) continue;

			// require bit zero for beams 
			//
                        if (isData && !(cms2.l1_techbits1() & (1<<0))) continue;

			// pixel digi requirement to get rid of monster events
			//
			if (nHighPurityTracks < 10 || cms2.trks_ndof().size() > 100) continue;
			if (cms2.trks_ndof().size()/float(nHighPurityTracks) < 0.20) continue;
			

			// fill iso histograms
			//

			// loop on vector of random conesEt in each event
			//

			std::vector<std::vector<float> > &conesEt = cms2.pseudo_ecalIso03_recHitEt();
                        std::vector<std::vector<float> > &conesE = cms2.pseudo_ecalIso03_recHitE();

			for (size_t d = 0; d < cms2.pseudo_ecalEta().size(); ++d)
			{

				if (fabs(cms2.pseudo_ecalEta()[d]) > 1.8) {
					h1_pseudo_ecalEta_ee->Fill(cms2.pseudo_ecalEta()[d]);
					h1_pseudo_ecalPhi_ee->Fill(cms2.pseudo_ecalPhi()[d]);
					h1_pseudo_ecalIso03_ee->Fill(cms2.pseudo_ecalIso03()[d]);

					h1_pseudo_ecalIso03_recHitN_ee->Fill(conesEt[d].size());
                                        float EtSum010 = 0.0;
                                        float EtSum015 = 0.0;
                                        for (size_t h = 0; h < conesEt[d].size(); ++h) {
                                                if (fabs(conesE[d][h]) > 0.10) EtSum010 += conesEt[d][h];
                                                if (fabs(conesE[d][h]) > 0.15) EtSum015 += conesEt[d][h];
                                                h1_pseudo_ecalIso03_recHitEt_ee->Fill(conesEt[d][h]);
                                                h1_pseudo_ecalIso03_recHitE_ee->Fill(conesE[d][h]);
                                        }
                                        h1_pseudo_ecalIso03_E015_ee->Fill(EtSum015);
                                        h1_pseudo_ecalIso03_E010_ee->Fill(EtSum010);

					if (cms2.pseudo_ecalEta()[d] > 1.8) h1_pseudo_ecalIso03_eep->Fill(cms2.pseudo_ecalIso03()[d]);
					if (cms2.pseudo_ecalEta()[d] < -1.8) h1_pseudo_ecalIso03_eem->Fill(cms2.pseudo_ecalIso03()[d]);
					h1_pseudo_hcalD1Iso03_ee->Fill(cms2.pseudo_hcalD1Iso03()[d]);
					h1_pseudo_hcalD2Iso03_ee->Fill(cms2.pseudo_hcalD2Iso03()[d]);
					h1_pseudo_tkIso03_ee->Fill(cms2.pseudo_tkIso03()[d]);
					h1_pseudo_dRClosestTower_ee->Fill(cms2.pseudo_dRClosestTower()[d]);
					if (nGoodTracksVz > 1) h1_pseudo_ecalIso03_qual_ee->Fill(cms2.pseudo_ecalIso03()[d]);
					if (cms2.pseudo_dRClosestTower()[d] > 1.0) h1_pseudo_ecalIso03_fakeSR_ee->Fill(cms2.pseudo_ecalIso03()[d]);
				}

				if (fabs(cms2.pseudo_ecalEta()[d]) < 1.2) {
					h1_pseudo_ecalEta_eb->Fill(cms2.pseudo_ecalEta()[d]);
					h1_pseudo_ecalPhi_eb->Fill(cms2.pseudo_ecalPhi()[d]);
					h1_pseudo_ecalIso03_eb->Fill(cms2.pseudo_ecalIso03()[d]);

                                        h1_pseudo_ecalIso03_recHitN_eb->Fill(conesEt[d].size());
                                        float EtSum010 = 0.0;
                                        float EtSum015 = 0.0;
                                        for (size_t h = 0; h < conesEt[d].size(); ++h) {
                                                if (fabs(conesE[d][h]) > 0.10) EtSum010 += conesEt[d][h];
						if (fabs(conesE[d][h]) > 0.15) EtSum015 += conesEt[d][h];
	                                       	h1_pseudo_ecalIso03_recHitEt_eb->Fill(conesEt[d][h]);
						h1_pseudo_ecalIso03_recHitE_eb->Fill(conesE[d][h]);
					}
                                        h1_pseudo_ecalIso03_E015_eb->Fill(EtSum015);
                                        h1_pseudo_ecalIso03_E010_eb->Fill(EtSum010);


					if (cms2.pseudo_ecalEta()[d] > 0.3) h1_pseudo_ecalIso03_ebp->Fill(cms2.pseudo_ecalIso03()[d]);
					if (cms2.pseudo_ecalEta()[d] < -0.3) h1_pseudo_ecalIso03_ebm->Fill(cms2.pseudo_ecalIso03()[d]);
					h1_pseudo_hcalD1Iso03_eb->Fill(cms2.pseudo_hcalD1Iso03()[d]);
					h1_pseudo_hcalD2Iso03_eb->Fill(cms2.pseudo_hcalD2Iso03()[d]);
					h1_pseudo_tkIso03_eb->Fill(cms2.pseudo_tkIso03()[d]);
					h1_pseudo_dRClosestTower_eb->Fill(cms2.pseudo_dRClosestTower()[d]);
					if (nGoodTracksVz > 1) h1_pseudo_ecalIso03_qual_eb->Fill(cms2.pseudo_ecalIso03()[d]);
					if (cms2.pseudo_dRClosestTower()[d] > 1.0) h1_pseudo_ecalIso03_fakeSR_eb->Fill(cms2.pseudo_ecalIso03()[d]);
				}


			} // end loop on random directions

		} // end loop on files
	} // end loop on events

	if ( nEventsChain != nEventsTotal ) {
		std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
	}

	return 0;
}
