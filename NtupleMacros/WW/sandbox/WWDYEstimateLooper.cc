#include "../CORE/selections.h"
#include "../CORE/CMS2.h"
#include "WWLooper.h"

#include "../CORE/utilities.h"

#include "TFile.h"
#include "TTree.h"

WWDYEstimateLooper::WWDYEstimateLooper (Sample s, cuts_t c, const char *fname)
     : WWLooperBase(s, c, fname),

	// histograms
	dh1_mll_ (sample, "dh1_mll_", 200, 0.0, 400.0),
	dh1_hyp_met_ (sample, "dh1_hyp_met_", 400, 0.0, 400.0),
	dh1_n_trkjet_ (sample, "dh1_n_trkjet_", 10, -0.5, 9.5),
        dh1_n_trkjet_met15_ (sample, "dh1_n_trkjet_met15_", 10, -0.5, 9.5),

        dh1_mll_0j_met15_ (sample, "dh1_mll_0j_met15_", 400, 0.0, 400.0),
        dh1_mll_0j_met25_ (sample, "dh1_mll_0j_met25_", 400, 0.0, 400.0),
        dh1_mll_0j_met35_ (sample, "dh1_mll_0j_met35_", 400, 0.0, 400.0),
        dh1_mll_0j_met45_ (sample, "dh1_mll_0j_met45_", 400, 0.0, 400.0),
        dh1_mll_0j_met55_ (sample, "dh1_mll_0j_met55_", 400, 0.0, 400.0),
        dh1_mll_0j_met65_ (sample, "dh1_mll_0j_met65_", 400, 0.0, 400.0),
        dh1_mll_0j_met75_ (sample, "dh1_mll_0j_met75_", 400, 0.0, 400.0),
        dh1_mll_0j_met85_ (sample, "dh1_mll_0j_met85_", 400, 0.0, 400.0),

        dh1_mll_1j_met15_ (sample, "dh1_mll_1j_met15_", 400, 0.0, 400.0),
        dh1_mll_1j_met25_ (sample, "dh1_mll_1j_met25_", 400, 0.0, 400.0),
        dh1_mll_1j_met35_ (sample, "dh1_mll_1j_met35_", 400, 0.0, 400.0),
        dh1_mll_1j_met45_ (sample, "dh1_mll_1j_met45_", 400, 0.0, 400.0),
        dh1_mll_1j_met55_ (sample, "dh1_mll_1j_met55_", 400, 0.0, 400.0),
        dh1_mll_1j_met65_ (sample, "dh1_mll_1j_met65_", 400, 0.0, 400.0),
        dh1_mll_1j_met75_ (sample, "dh1_mll_1j_met75_", 400, 0.0, 400.0),
        dh1_mll_1j_met85_ (sample, "dh1_mll_1j_met85_", 400, 0.0, 400.0),

        dh1_mll_2j_met15_ (sample, "dh1_mll_2j_met15_", 400, 0.0, 400.0),
        dh1_mll_2j_met25_ (sample, "dh1_mll_2j_met25_", 400, 0.0, 400.0),
        dh1_mll_2j_met35_ (sample, "dh1_mll_2j_met35_", 400, 0.0, 400.0),
        dh1_mll_2j_met45_ (sample, "dh1_mll_2j_met45_", 400, 0.0, 400.0),
        dh1_mll_2j_met55_ (sample, "dh1_mll_2j_met55_", 400, 0.0, 400.0),
        dh1_mll_2j_met65_ (sample, "dh1_mll_2j_met65_", 400, 0.0, 400.0),
        dh1_mll_2j_met75_ (sample, "dh1_mll_2j_met75_", 400, 0.0, 400.0),
        dh1_mll_2j_met85_ (sample, "dh1_mll_2j_met85_", 400, 0.0, 400.0),

        dh1_mll_3j_met15_ (sample, "dh1_mll_3j_met15_", 400, 0.0, 400.0),
        dh1_mll_3j_met25_ (sample, "dh1_mll_3j_met25_", 400, 0.0, 400.0),
        dh1_mll_3j_met35_ (sample, "dh1_mll_3j_met35_", 400, 0.0, 400.0),
        dh1_mll_3j_met45_ (sample, "dh1_mll_3j_met45_", 400, 0.0, 400.0),
        dh1_mll_3j_met55_ (sample, "dh1_mll_3j_met55_", 400, 0.0, 400.0),
        dh1_mll_3j_met65_ (sample, "dh1_mll_3j_met65_", 400, 0.0, 400.0),
        dh1_mll_3j_met75_ (sample, "dh1_mll_3j_met75_", 400, 0.0, 400.0),
        dh1_mll_3j_met85_ (sample, "dh1_mll_3j_met85_", 400, 0.0, 400.0),

        dh1_mll_4j_met15_ (sample, "dh1_mll_4j_met15_", 400, 0.0, 400.0),
        dh1_mll_4j_met25_ (sample, "dh1_mll_4j_met25_", 400, 0.0, 400.0),
        dh1_mll_4j_met35_ (sample, "dh1_mll_4j_met35_", 400, 0.0, 400.0),
        dh1_mll_4j_met45_ (sample, "dh1_mll_4j_met45_", 400, 0.0, 400.0),
        dh1_mll_4j_met55_ (sample, "dh1_mll_4j_met55_", 400, 0.0, 400.0),
        dh1_mll_4j_met65_ (sample, "dh1_mll_4j_met65_", 400, 0.0, 400.0),
        dh1_mll_4j_met75_ (sample, "dh1_mll_4j_met75_", 400, 0.0, 400.0),
        dh1_mll_4j_met85_ (sample, "dh1_mll_4j_met85_", 400, 0.0, 400.0)

{
	std::cout << "[WWDYEstimateLooper::WWDYEstimateLooper]" << std::endl;
	std::string outFileName = "DYEstimateLooper_" + sample.name + ".root";
        std::cout << "making root file " << outFileName << std::endl;

	outFile_ = new TFile(outFileName.c_str(), "RECREATE");
	outFile_->cd();

   	outTree_ = new TTree("T1", "Tree");

   	// event variables
   	outTree_->Branch("sampleId", &sampleId_, "sampleId/I");
        outTree_->Branch("weight", &weight_, "weight/I");

}

WWDYEstimateLooper::~WWDYEstimateLooper ()
{
        std::cout << "[WWDYEstimateLooper::~WWDYEstimateLooper]" << std::endl;
	std::cout << "This looper was analyzing... " << sample.process << "(" << sample.name.c_str() << ")" << std::endl;
	std::cout << "Saving root file " << std::endl;

	outFile_->cd();
   	outTree_->Write();
   	outFile_->Close();
   	delete outFile_;
}


void WWDYEstimateLooper::Dilep (int i_hyp)
{
     cuts_t cuts_passed = DilepSelect(i_hyp);
//      printf("cuts passed: %x, need: %x\n", cuts_passed, ww_old_baseline_cuts);

     // eventually, we want to do the N - 1 plots 
     // ...

     FillHists(i_hyp);

     // for now, we shoot first, ask questions later
     if ((cuts_passed & cuts) != cuts)
          return;

     // if we get an extra muon, print it out
     //int tag_muon_idx = tagMuonIdx(i_hyp);

     enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     double weight =  cms2.evt_scale1fb() * sample.kFactor;

/*
     if (tag_muon_idx != -1 && myType == DILEPTON_EMU)
          printf("tag muon: %d (out of %u), pt = %.1f (%.1f), iso = %f\n",
                 tag_muon_idx, cms2.mus_p4().size(), cms2.mus_p4()[tag_muon_idx].pt(),
                 tagMuonPt(i_hyp), tagMuonRelIso(i_hyp));
*/


	// get the met with F. Golf correction applied
	const TVector3 trkCorr = correctMETforTracks();
  	TVector3 hyp_met;
  	hyp_met.SetPtEtaPhi(cms2.hyp_met()[i_hyp], 0, cms2.hyp_metPhi()[i_hyp]);
  	hyp_met += trkCorr;

	// messing around with DileptonHist
	dh1_mll_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);

	if (nTrkJets(i_hyp) == 0 ){
        	if (hyp_met.Pt() > 15.0) dh1_mll_0j_met15_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
        	if (hyp_met.Pt() > 25.0) dh1_mll_0j_met25_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
        	if (hyp_met.Pt() > 35.0) dh1_mll_0j_met35_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
        	if (hyp_met.Pt() > 45.0) dh1_mll_0j_met45_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
        	if (hyp_met.Pt() > 55.0) dh1_mll_0j_met55_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
        	if (hyp_met.Pt() > 65.0) dh1_mll_0j_met65_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 75.0) dh1_mll_0j_met75_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 85.0) dh1_mll_0j_met85_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
	}

        if (nTrkJets(i_hyp) == 1 ){
                if (hyp_met.Pt() > 15.0) dh1_mll_1j_met15_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 25.0) dh1_mll_1j_met25_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 35.0) dh1_mll_1j_met35_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 45.0) dh1_mll_1j_met45_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 55.0) dh1_mll_1j_met55_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 65.0) dh1_mll_1j_met65_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 75.0) dh1_mll_1j_met75_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 85.0) dh1_mll_1j_met85_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
        }

        if (nTrkJets(i_hyp) == 2 ){
                if (hyp_met.Pt() > 15.0) dh1_mll_2j_met15_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 25.0) dh1_mll_2j_met25_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 35.0) dh1_mll_2j_met35_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 45.0) dh1_mll_2j_met45_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 55.0) dh1_mll_2j_met55_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 65.0) dh1_mll_2j_met65_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 75.0) dh1_mll_2j_met75_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 85.0) dh1_mll_2j_met85_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
        }

        if (nTrkJets(i_hyp) == 3 ){
                if (hyp_met.Pt() > 15.0) dh1_mll_3j_met15_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 25.0) dh1_mll_3j_met25_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 35.0) dh1_mll_3j_met35_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 45.0) dh1_mll_3j_met45_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 55.0) dh1_mll_3j_met55_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 65.0) dh1_mll_3j_met65_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 75.0) dh1_mll_3j_met75_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 85.0) dh1_mll_3j_met85_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
        }

        if (nTrkJets(i_hyp) == 4 ){
                if (hyp_met.Pt() > 15.0) dh1_mll_4j_met15_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 25.0) dh1_mll_4j_met25_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 35.0) dh1_mll_4j_met35_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 45.0) dh1_mll_4j_met45_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 55.0) dh1_mll_4j_met55_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 65.0) dh1_mll_4j_met65_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 75.0) dh1_mll_4j_met75_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
                if (hyp_met.Pt() > 85.0) dh1_mll_4j_met85_.Fill(myType, cms2.hyp_p4()[i_hyp].M(), weight);
        }


	dh1_hyp_met_.Fill(myType, hyp_met.Pt(), weight);

	dh1_n_trkjet_.Fill(myType, nTrkJets(i_hyp), weight);
	if (hyp_met.Pt() > 15.0) dh1_n_trkjet_met15_.Fill(myType, nTrkJets(i_hyp), weight);

	// messing around with trees
	sampleId_ = sample.process;
	weight_ =  weight;
     	outTree_->Fill();

}

