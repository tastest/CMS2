#include <math.h>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"

#include "Math/LorentzVector.h"

Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname)
{
     // zero out the candidate counters (don't comment this out)
     memset(cands_passing_	, 0, sizeof(cands_passing_       ));
     memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
     memset(cands_count_		, 0, sizeof(cands_count_         ));
}

void Looper::BookHistos ()
{

     // set up root file
     std::string outFileName = "Looper_" + sample_.name + ".root";
     std::cout << "[Looper::BookHistos] making root file " << outFileName << std::endl;

     outFile_ = new TFile(outFileName.c_str(), "RECREATE");
     outFile_->cd();
     outTree_ = new TTree("T1", "Tree");
  
     // book the branches
     outTree_->Branch("sampleId", &sampleId_, "sampleId/I");
     outTree_->Branch("evt_weight", &evt_weight_, "evt_weight/F");

	// event properties
     outTree_->Branch("z_pt", &z_pt_, "z_pt/F");
     outTree_->Branch("z_p", &z_p_, "z_p/F");
     outTree_->Branch("evt_tcmet", &evt_tcmet_, "evt_tcmet/F");

	// electron branches
     outTree_->Branch("ele_count", &ele_count_, "ele_count/I");

     outTree_->Branch("ele_pt", &ele_pt_, "ele_pt[ele_count]/F");
     outTree_->Branch("ele_p", &ele_pt_, "ele_p[ele_count]/F");
     outTree_->Branch("ele_eta", &ele_eta_, "ele_eta[ele_count]/F");
     //outTree_->Branch("ele_etaDet", &ele_etaDet_, "ele_etaDet[ele_count]/F");

     outTree_->Branch("ele_hOverE", &ele_hOverE_, "ele_hOverE[ele_count]/F");
     outTree_->Branch("ele_dPhiIn", &ele_dPhiIn_, "ele_dPhiIn[ele_count]/F");
     outTree_->Branch("ele_dEtaIn", &ele_dEtaIn_, "ele_dEtaIn[ele_count]/F");
     outTree_->Branch("ele_sigmaIEtaIEta", &ele_sigmaIEtaIEta_, "ele_sigmaIEtaIEta[ele_count]/F");
     outTree_->Branch("ele_tkIso_", &ele_tkIso_, "ele_tkIso[ele_count]/F");

     outTree_->Branch("ele_matchMC", &ele_matchMC_, "ele_matchMC[ele_count]/I");

}

bool Looper::FilterEvent()
{ 
  return false; 
}

cuts_t Looper::EventSelect ()
{
     //------------------------------------------------------------
     // In an event-based analysis, you would make your cuts here
     //------------------------------------------------------------

     cuts_t ret = 0;
     return ret;
}

void Looper::FillEventHistos ()
{
     //------------------------------------------------------------
     // In an event-based analysis, you would fill your histos here
     //------------------------------------------------------------

	// fill the sampleId and weight branch
	sampleId_ = sample_.process;
	evt_weight_ = cms2.evt_scale1fb() * sample_.kFactor;

	// reset counters
	ele_count_ = 0;

	// event properties
	// get the z momentum
        int zidx = -1;
  	for (size_t i = 0; i < cms2.genps_id().size(); ++i)
	{
        	if( cms2.genps_id()[i] == 23 ) 
		{
          		zidx = i;
        	}
	}

	z_p_ = -1;
	z_pt_ = -1;
	if (zidx != -1)
	{
		z_p_ = cms2.genps_p4()[zidx].P();
		z_pt_ = cms2.genps_p4()[zidx].Pt();
	}

	// met
	evt_tcmet_ = cms2.evt_tcmet();

	// loop on electrons
	for (size_t i = 0; i < cms2.evt_nels(); ++i)
	{
		if (cms2.els_p4()[i].Pt() > 5.0) 
		{
			// Incriment counter and add to tree
			ele_pt_[ele_count_] = cms2.els_p4()[i].Pt();
			ele_p_[ele_count_] = cms2.els_p4()[i].P();
			ele_eta_[ele_count_] = cms2.els_p4()[i].Eta();
			//ele_etaDet_[ele_count_] = cms2.els_etaSC()[i];

			ele_hOverE_[ele_count_] = cms2.els_hOverE()[i];
			ele_sigmaIEtaIEta_[ele_count_] = cms2.els_sigmaIEtaIEta()[i];
                        ele_dEtaIn_[ele_count_] = cms2.els_dEtaIn()[i];
                        ele_dPhiIn_[ele_count_] = cms2.els_dPhiIn()[i];
                        ele_tkIso_[ele_count_] = cms2.els_tkIso()[i];


			// MC matching (temporary fudge)
			Int_t matchMC = 0;
			if (abs(cms2.els_mc3_id()[i]) == 11) matchMC = 1;
			ele_matchMC_[ele_count_] = matchMC;

                        ele_count_ ++;

		}
	}

        // when everything is set
        // fill this entry in the tree
	outTree_->Fill();

}

void Looper::End ()
{
     //------------------------------------------------------------
     //Example status message at the end of a looper; edit for your
     //application
     //------------------------------------------------------------

     int ret = fprintf(logfile_, 
		       "Sample %10s: Total candidate count (ee em mm all): %8u %8u %8u %8u."
		       " Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n",   
		       sample_.name.c_str(),
		       CandsCount(DILEPTON_EE), CandsCount(DILEPTON_EMU), CandsCount(DILEPTON_MUMU), CandsCount(DILEPTON_ALL), 
		       CandsPassing(DILEPTON_EE)  , RMS(DILEPTON_EE),  
		       CandsPassing(DILEPTON_EMU) , RMS(DILEPTON_EMU),  
		       CandsPassing(DILEPTON_MUMU), RMS(DILEPTON_MUMU), 
		       CandsPassing(DILEPTON_ALL) , RMS(DILEPTON_ALL));
     if (ret < 0)
	  perror("writing to log file");

      // Save tree
      std::cout << "[Looper::End] Saving root file " << std::endl;
      outFile_->cd();
      outTree_->Write();
      outFile_->Close();
      delete outFile_; 

}

// matching
float Looper::match(unsigned int leptonIndex)
{


}


