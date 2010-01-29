
//
// ttbar -> ll
// Dave "the one but not the only" Evans 
//

// C++ includes
#include <iostream>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"

#include "Math/LorentzVector.h"


// CMS2 includes
#include "CMS2.h"
#include "../../CORE/electronSelections.h"
#include "../../CORE/selections.h"

#include "../../Tools/DileptonHypType.h"

//
// Namespaces
//
using namespace tas;

double dRbetweenVectors(const LorentzVector &vec1,
                        const LorentzVector &vec2 ){

  double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  double deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
} 

//   jets_p4   
/*
std::vector<LorentzVector> getCorCaloJets(int i_hyp) {
  std::vector<LorentzVector> calo_jets;
  calo_jets.clear();

  for (unsigned int jj=0; jj < cms2.jets_p4().size(); ++jj) {
    if ((dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.jets_p4()[jj]) < 0.4)||
        (dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.jets_p4()[jj]) < 0.4)
        ) continue;
    if (cms2.jets_p4()[jj].pt()*cms2.jets_cor()[jj] < 30) continue;
    if (fabs(cms2.jets_p4()[jj].Eta()) > 2.4) continue;
    //fkw July21 2009 if (cms2.jets_emFrac()[jj] < 0.1) continue;
    calo_jets.push_back(cms2.jets_p4()[jj]*cms2.jets_cor()[jj]);
  }

  if (calo_jets.size() > 1) {
       sort(calo_jets.begin(), calo_jets.end(),  comparePt);
  }
  return calo_jets;
}
*/

enum DileptonHypType hyp_typeToHypType (int hyp_type)
{
     switch (hyp_type) {
     case 0:
          return DILEPTON_MUMU;
     case 1: case 2:
          return DILEPTON_EMU;
     case 3:
          return DILEPTON_EE;
     default:
          assert(hyp_type < 4);
     }
     return DILEPTON_ALL;
}

void Fill(TH1F** hist, unsigned int hyp, float val, float weight)
{
	hist[hyp]->Fill(val, weight);
	hist[DILEPTON_ALL]->Fill(val, weight);
}

void FormatHist(TH1F** hist, std::string sampleName, std::string name, Int_t n, Float_t min, Float_t max)
{       
	// loop on EB, EE
	for (unsigned int i = 0; i < 4; ++i)
	{
		std::string str = dilepton_hypo_names[i];
		std::string title = name + "_" + str;
		hist[i] = new TH1F(Form("%s_%s_%s", sampleName.c_str(), name.c_str(), str.c_str()),
				title.c_str(), n, min, max);
		hist[i]->GetXaxis()->SetTitle(name.c_str());
	}
}    

//
// Main function
//
int ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents = -1, std::string skimFilePrefix="") {

	TObjArray *listOfFiles = chain->GetListOfFiles();

	unsigned int nEventsChain=0;
	if(nEvents==-1) 
		nEvents = chain->GetEntries();
	nEventsChain = nEvents;
	unsigned int nEventsTotal = 0;
	TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

	//
	// define counters
	//

	// count the (weighted and unweighted) number of candidates passing our cuts
	double             cands_passing[4];
	double             cands_passing_w2[4];
	unsigned int       cands_count[4];
	unsigned int 	   hyp_count[4];
        memset(cands_passing   , 0, sizeof(cands_passing       ));
        memset(cands_passing_w2        , 0, sizeof(cands_passing_w2    ));
        memset(cands_count             , 0, sizeof(cands_count         ));
        memset(hyp_count             , 0, sizeof(hyp_count         ));

	//
	// define histograms
	//

	TH1F *h1_hyp_lt_ee_pt[4];
        TH1F *h1_hyp_lt_eb_pt[4];
        FormatHist(h1_hyp_lt_ee_pt, sampleName, "hyp_lt_ee_pt", 20, 0.0, 100.0);
        FormatHist(h1_hyp_lt_eb_pt, sampleName, "hyp_lt_eb_pt", 20, 0.0, 100.0);

	// basic selection quantities
	//
	TH1F *h1_hyp_lt_eb_hoe[4];
        TH1F *h1_hyp_lt_eb_sigmaIEtaIEta[4];
        TH1F *h1_hyp_lt_eb_dEtaIn[4];
        TH1F *h1_hyp_lt_eb_dPhiIn[4];
        TH1F *h1_hyp_lt_eb_d0[4];
        TH1F *h1_hyp_lt_eb_E2x5MaxOver5x5[4];
	TH1F *h1_hyp_lt_eb_ecalIso[4];
        TH1F *h1_hyp_lt_eb_hcalIso[4];
        TH1F *h1_hyp_lt_eb_tkIso[4];

        TH1F *h1_hyp_lt_ee_hoe[4];
        TH1F *h1_hyp_lt_ee_sigmaIEtaIEta[4];
        TH1F *h1_hyp_lt_ee_dEtaIn[4];
        TH1F *h1_hyp_lt_ee_dPhiIn[4];
        TH1F *h1_hyp_lt_ee_d0[4];
        TH1F *h1_hyp_lt_ee_E2x5MaxOver5x5[4];
        TH1F *h1_hyp_lt_ee_ecalIso[4];
        TH1F *h1_hyp_lt_ee_hcalIso[4];
        TH1F *h1_hyp_lt_ee_tkIso[4];

        FormatHist(h1_hyp_lt_eb_hoe, sampleName, "hyp_lt_eb_hoe", 50, 0.0, 0.05);
        FormatHist(h1_hyp_lt_ee_hoe, sampleName, "hyp_lt_ee_hoe", 50, 0.0, 0.05);

        FormatHist(h1_hyp_lt_eb_sigmaIEtaIEta, sampleName, "hyp_lt_eb_sigmaIEtaIEta", 50, 0.0, 0.05);
        FormatHist(h1_hyp_lt_ee_sigmaIEtaIEta, sampleName, "hyp_lt_ee_sigmaIEtaIEta", 50, 0.0, 0.05);

        FormatHist(h1_hyp_lt_eb_dEtaIn, sampleName, "hyp_lt_eb_dEtaIn", 100, -0.05, 0.05);
        FormatHist(h1_hyp_lt_ee_dEtaIn, sampleName, "hyp_lt_ee_dEtaIn", 100, -0.05, 0.05);

        FormatHist(h1_hyp_lt_eb_dPhiIn, sampleName, "hyp_lt_eb_dPhiIn", 100, -0.05, 0.05);
        FormatHist(h1_hyp_lt_ee_dPhiIn, sampleName, "hyp_lt_ee_dPhiIn", 100, -0.05, 0.05);

        FormatHist(h1_hyp_lt_eb_d0, sampleName, "hyp_lt_eb_d0", 100, -0.05, 0.05);
        FormatHist(h1_hyp_lt_ee_d0, sampleName, "hyp_lt_ee_d0", 100, -0.05, 0.05);

        FormatHist(h1_hyp_lt_eb_E2x5MaxOver5x5, sampleName, "hyp_lt_eb_E2x5MaxOver5x5", 110, 0.0, 1.1);
        FormatHist(h1_hyp_lt_ee_E2x5MaxOver5x5, sampleName, "hyp_lt_ee_E2x5MaxOver5x5", 110, 0.0, 1.1);

        FormatHist(h1_hyp_lt_eb_ecalIso, sampleName, "hyp_lt_eb_ecalIso", 100, 0, 25);
        FormatHist(h1_hyp_lt_eb_hcalIso, sampleName, "hyp_lt_eb_hcalIso", 100, 0, 25);
        FormatHist(h1_hyp_lt_eb_tkIso, sampleName, "hyp_lt_eb_tkIso", 100, 0, 25);

        FormatHist(h1_hyp_lt_ee_ecalIso, sampleName, "hyp_lt_ee_ecalIso", 100, 0, 25);
        FormatHist(h1_hyp_lt_ee_hcalIso, sampleName, "hyp_lt_ee_hcalIso", 100, 0, 25);
        FormatHist(h1_hyp_lt_ee_tkIso, sampleName, "hyp_lt_ee_tkIso", 100, 0, 25);


        TH1F *h1_hyp_lt_eb_pt_idnew[4];
        TH1F *h1_hyp_lt_ee_pt_idnew[4];
        TH1F *h1_hyp_lt_eb_pt_idold[4];
        TH1F *h1_hyp_lt_ee_pt_idold[4];
        FormatHist(h1_hyp_lt_eb_pt_idnew, sampleName, "hyp_lt_eb_pt_idnew", 20, 0.0, 100.0);
        FormatHist(h1_hyp_lt_ee_pt_idnew, sampleName, "hyp_lt_ee_pt_idnew", 20, 0.0, 100.0);
        FormatHist(h1_hyp_lt_eb_pt_idold, sampleName, "hyp_lt_eb_pt_idold", 20, 0.0, 100.0);
        FormatHist(h1_hyp_lt_ee_pt_idold, sampleName, "hyp_lt_ee_pt_idold", 20, 0.0, 100.0);

        TH1F *h1_hyp_lt_eb_pt_isonew[4];
        TH1F *h1_hyp_lt_ee_pt_isonew[4];
        TH1F *h1_hyp_lt_eb_pt_isoold[4];
        TH1F *h1_hyp_lt_ee_pt_isoold[4];
        FormatHist(h1_hyp_lt_eb_pt_isonew, sampleName, "hyp_lt_eb_pt_isonew", 20, 0.0, 100.0);
        FormatHist(h1_hyp_lt_ee_pt_isonew, sampleName, "hyp_lt_ee_pt_isonew", 20, 0.0, 100.0);
        FormatHist(h1_hyp_lt_eb_pt_isoold, sampleName, "hyp_lt_eb_pt_isoold", 20, 0.0, 100.0);
        FormatHist(h1_hyp_lt_ee_pt_isoold, sampleName, "hyp_lt_ee_pt_isoold", 20, 0.0, 100.0);

	// file loop
	//
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

			// print out event being processed
			if (nEventsTotal % 1000 == 0)
				std::cout << "Event: " << nEventsTotal << std::endl;

			// work out event weight
		        float weight = cms2.evt_scale1fb();

			//
			// loop on hypothesis
			//
			std::vector<unsigned int> hyp_index_selected;
			for (size_t h = 0; h < cms2.hyp_type().size(); ++h) {

                                DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[h]);

				// apply truth match
				if (abs(cms2.hyp_lt_id()[h]) == 11 && !(abs(cms2.hyp_lt_mc_id()[h]) == 11));

				// apply part of electron denominator first here
				if (abs(cms2.hyp_lt_id()[h]) == 11)
					if (!electron20Eta2p4(cms2.hyp_lt_index()[h])) continue;
                                if (abs(cms2.hyp_ll_id()[h]) == 11)
                                        if (!electron20Eta2p4(cms2.hyp_ll_index()[h])) continue;

				//
				//
			
				//
				// fill denominator histograms
				//
                                if (abs(cms2.hyp_lt_p4()[h].eta()) > 1.5 && abs(cms2.hyp_lt_id()[h]) == 11) {
					Fill(h1_hyp_lt_ee_pt, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                                        Fill(h1_hyp_lt_ee_hoe, hypType, cms2.els_hOverE()[cms2.hyp_lt_index()[h]], weight);
                                        Fill(h1_hyp_lt_ee_d0, hypType, cms2.els_d0corr()[cms2.hyp_lt_index()[h]], weight);
                                        Fill(h1_hyp_lt_ee_dPhiIn, hypType, cms2.els_dPhiIn()[cms2.hyp_lt_index()[h]], weight);
                                        Fill(h1_hyp_lt_ee_dEtaIn, hypType, cms2.els_dEtaIn()[cms2.hyp_lt_index()[h]], weight);
                                        Fill(h1_hyp_lt_ee_sigmaIEtaIEta, hypType, cms2.els_sigmaIEtaIEta()[cms2.hyp_lt_index()[h]], weight);
					float E2x5MaxOver5x5 = cms2.els_e2x5Max()[cms2.hyp_lt_index()[h]] / cms2.els_e5x5()[cms2.hyp_lt_index()[h]];
                                        Fill(h1_hyp_lt_ee_E2x5MaxOver5x5, hypType, E2x5MaxOver5x5, weight);
                                        Fill(h1_hyp_lt_ee_ecalIso, hypType, cms2.els_ecalIso()[h], weight);
                                        Fill(h1_hyp_lt_ee_hcalIso, hypType, cms2.els_hcalIso()[h], weight);
                                        Fill(h1_hyp_lt_ee_tkIso, hypType, cms2.els_tkIso()[h], weight);
				}

                                if (abs(cms2.hyp_lt_p4()[h].eta()) < 1.5 && abs(cms2.hyp_lt_id()[h]) == 11) {
					Fill(h1_hyp_lt_eb_pt, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                                        Fill(h1_hyp_lt_eb_hoe, hypType, cms2.els_hOverE()[cms2.hyp_lt_index()[h]], weight);
                                        Fill(h1_hyp_lt_eb_d0, hypType, cms2.els_d0corr()[cms2.hyp_lt_index()[h]], weight);
                                        Fill(h1_hyp_lt_eb_dPhiIn, hypType, cms2.els_dPhiIn()[cms2.hyp_lt_index()[h]], weight);
                                        Fill(h1_hyp_lt_eb_dEtaIn, hypType, cms2.els_dEtaIn()[cms2.hyp_lt_index()[h]], weight);
                                        Fill(h1_hyp_lt_eb_sigmaIEtaIEta, hypType, cms2.els_sigmaIEtaIEta()[cms2.hyp_lt_index()[h]], weight);
                                        float E2x5MaxOver5x5 = cms2.els_e2x5Max()[cms2.hyp_lt_index()[h]] / cms2.els_e5x5()[cms2.hyp_lt_index()[h]];
                                        Fill(h1_hyp_lt_eb_E2x5MaxOver5x5, hypType, E2x5MaxOver5x5, weight);
					Fill(h1_hyp_lt_eb_ecalIso, hypType, cms2.els_ecalIso()[h], weight);
                                        Fill(h1_hyp_lt_eb_hcalIso, hypType, cms2.els_hcalIso()[h], weight);
                                        Fill(h1_hyp_lt_eb_tkIso, hypType, cms2.els_tkIso()[h], weight);
				}


				//
				//
				//

				// apply lepton id to lt
				bool ltPassOld = false;
				bool ltPassNew = false;
				bool relSusyIso = false;
				bool relSusyIso_cand0 = false;
				if (abs(cms2.hyp_lt_id()[h]) == 11) {
					if (cms2.els_egamma_looseId().at(cms2.hyp_lt_index()[h])) ltPassOld = true;
					if (electronId_cand01(cms2.hyp_lt_index()[h])) ltPassNew = true;
					if (electronIsolation_relsusy(cms2.hyp_lt_index()[h], true) < 0.10) relSusyIso = true;
                                        if (electronIsolation_relsusy_cand0(cms2.hyp_lt_index()[h], true) < 0.10) relSusyIso_cand0 = true;
				}

                                //
                                // fill numerator histograms
                                //
				if (abs(cms2.hyp_lt_p4()[h].eta()) > 1.5) {
                               		if (ltPassOld) Fill(h1_hyp_lt_ee_pt_idold, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                               		if (ltPassNew) Fill(h1_hyp_lt_ee_pt_idnew, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
					if (relSusyIso) Fill(h1_hyp_lt_ee_pt_isoold, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                                        if (relSusyIso_cand0) Fill(h1_hyp_lt_ee_pt_isonew, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
				}
				if (abs(cms2.hyp_lt_p4()[h].eta()) < 1.5)  {
                                       	if (ltPassOld) Fill(h1_hyp_lt_eb_pt_idold, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                               		if (ltPassNew) Fill(h1_hyp_lt_eb_pt_idnew, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                                        if (relSusyIso) Fill(h1_hyp_lt_eb_pt_isoold, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                                        if (relSusyIso_cand0) Fill(h1_hyp_lt_eb_pt_isonew, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
				}
				//
				//
				//

				// apply lepton id to ll
                                if (!looseLeptonSelectionNoIsoTTDil08(cms2.hyp_lt_id()[h], cms2.hyp_lt_index()[h])) continue;
                                if (!looseLeptonSelectionNoIsoTTDil08(cms2.hyp_ll_id()[h], cms2.hyp_ll_index()[h])) continue;

				// apply isolation to ll
                                if (!passLeptonIsolationTTDil08(cms2.hyp_ll_id()[h], cms2.hyp_ll_index()[h])) continue;
                                if (!passLeptonIsolationTTDil08(cms2.hyp_lt_id()[h], cms2.hyp_lt_index()[h])) continue;

				// opposite charge
				if (cms2.hyp_lt_charge()[h] * cms2.hyp_ll_charge()[h] > 0) continue;

				// z mass window
				if (inZmassWindow(cms2.hyp_p4()[h].M())) continue;
			
				//
				// store indices of hypothesese that pass this preselection
				//	
				hyp_index_selected.push_back(h);

			} // end loop on hypothesis

			//
			// perform hypothesis disambiguation
			//

/*
			int strasbourgDilType = 0;
			if (!hyp_index_selected.size()) continue;
			int hyp = eventDilIndexByWeightTTDil08(hyp_index_selected, strasbourgDilType, false, false);
		
			//
			// make requirements of the selected hypothesis
			//

			// trigger
			if (!passTriggersMu9orLisoE15(cms2.hyp_type()[hyp])) continue;

			// met
			if (!passMet_OF20_SF30(hyp, true)) continue;

			//
			// If we got to here then classify the event according to nJets
			//	

			std::vector<LorentzVector> corCaloJets = getCorCaloJets(hyp);

			if (corCaloJets.size() >= 2) {

				DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[hyp]);

				cands_passing[hypType] += weight;
				cands_passing_w2[hypType] += weight * weight;
				cands_count[hypType] ++;
                                cands_passing[DILEPTON_ALL] += weight;
                                cands_passing_w2[DILEPTON_ALL] += weight * weight;
                                cands_count[DILEPTON_ALL] ++;

			}
*/


		} // end loop on files

	} // end loop on events

	if ( nEventsChain != nEventsTotal ) {
		std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
	}

	std::cout << "sampleName \t";
        for (unsigned int i = 0; i < 4; ++i) {
                std::string str = dilepton_hypo_names[i];
		std::cout << dilepton_hypo_names[i] << "\t";
	}
	std::cout << std::endl;
        std::cout << sampleName << "\t";
	for (unsigned int i = 0; i < 4; ++i) {
		std::cout << cands_passing[i] << " $\pm$ " << sqrt(cands_passing_w2[i]) << "\t";
	}
	std::cout << std::endl;

	return 0;
}

