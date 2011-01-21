#include "babymakercommon.h"
#include "emubabymaker.h" 

#include <algorithm>
#include <iostream>
#include <string>

#include "TChain.h"
#include "TCollection.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TChainElement.h"

#include "../CORE/CMS2.h"
#include "../CORE/mcSelections.h"
#include "../CORE/electronSelections.h"
#include "../CORE/electronSelectionsParameters.h"
#include "../CORE/metSelections.h"
#include "../CORE/muonSelections.h"
#include "../CORE/trackSelections.h"

void emubabymaker::ScanChain (const char *inputFilename, const char *babyFilename, int nEvents)
{

    TChain *chain = new TChain("Events");
    chain->Add(inputFilename);
    TObjArray *listOfFiles = chain->GetListOfFiles();

    unsigned int nEventsChain=0;
    if (nEvents == -1) 
        nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;

    // make a baby ntuple
    MakeBabyNtuple(babyFilename);

    // file loop
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;
    while ((currentFile = (TFile*)fileIter.Next()))
    {
        TFile f(currentFile->GetTitle());
        TTree *tree = (TTree*)f.Get("Events");
        cms2.Init(tree);

        //Event Loop
        unsigned int nEvents = tree->GetEntries();
        for (unsigned int event = 0; event < nEvents; ++event)
        {
            cms2.GetEntry(event);
            ++nEventsTotal;
            // Progress feedback to the user
            if (nEventsTotal%1000 == 0)
            {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1))
                {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                            "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
                    fflush(stdout);
                }
            }

            // muon stuff
            for (unsigned mui = 0; mui < cms2.mus_p4().size(); ++mui)
            {
                // global and tracker muons only
                if (!(cms2.mus_type()[mui] & 6))
                    continue;
                // pt > 20
                if (cms2.mus_p4()[mui].pt() <= 20.)
                    continue;

                // initialize baby quantities
                InitBabyNtuple();

                // event stuff
                strcpy(dataset_, cms2.evt_dataset().Data());
                run_        = cms2.evt_run();
                ls_         = cms2.evt_lumiBlock();
                evt_        = cms2.evt_event();
                isdata_     = cms2.evt_isRealData();
                nvtx_ = 0;
                for (unsigned int i = 0; i < cms2.vtxs_isFake().size(); i++) {
                    if (!cms2.vtxs_isFake()[i] && cms2.vtxs_isValid()[i])
                        ++nvtx_;
                }										

                if (!isdata_) {
                    int nlep  = leptonGenpCount_lepTauDecays(ngenels_, ngenmus_, ngentaus_);
                    scale1fb_ = cms2.evt_scale1fb();
                    pthat_    = cms2.genps_pthat();
                }

                pfmet_      = cms2.evt_pfmet();
                tcmet_      = cms2.evt_pf_tcmet();
                calotcmet_  = cms2.evt_tcmet();
                ntrks_      = cms2.trks_trk_p4().size();

                float thePFMetPhi     = cms2.evt_pfmetPhi();
                float theTCMetPhi     = cms2.evt_pf_tcmetPhi();
                float theCaloTCMetPhi = cms2.evt_tcmetPhi();

                float metx  = tcmet_ * cos(theTCMetPhi);
                float mety  = tcmet_ * sin(theTCMetPhi);
                float cmetx = calotcmet_ * cos(theCaloTCMetPhi);
                float cmety = calotcmet_ * sin(theCaloTCMetPhi);
                for (unsigned int muj = 0; muj < cms2.mus_p4().size(); ++muj) {
                    if (cms2.mus_p4()[muj].Pt() <= 10.0) continue;
                    if (!wasMetCorrectedForThisMuon(muj, usingTcMet) && muonIdNotIsolated(muj, NominalTTbarV2)) {
                        fixMetForThisMuon(muj, metx, mety, usingTcMet);
                        fixMetForThisMuon(muj, cmetx, cmety, usingTcMet);
                    }
                }
                tcmet_ = sqrt(metx * metx + mety * mety);
                theTCMetPhi = atan2(mety, metx);
                calotcmet_ = sqrt(cmetx * cmetx + cmety * cmety);
                theCaloTCMetPhi = atan2(cmety, cmetx);

                // loop over muons and electrons to get ngoodlep
                ngoodlep_ = 0;
                ngoodmus_ = 0;
                ngoodels_ = 0;
                for (unsigned muii = 0; muii < cms2.mus_p4().size(); ++muii) {
                    if (cms2.mus_p4()[muii].pt() <= 20.) continue;
                    if (muonId(muii, NominalTTbarV2)) {
                        ++ngoodlep_;
                        ++ngoodmus_;
                    }
                }
                for (unsigned eli = 0; eli < cms2.els_p4().size(); ++eli) {
                    if (cms2.els_p4()[eli].pt() > 20. && pass_electronSelection(eli, electronSelection_ttbarV2)) {
                        ++ngoodlep_;
                        ++ngoodels_;
                    }
                }

                eormu_    = 13*cms2.mus_charge()[mui]*-1; // *-1 to follow pdg conventions
                type_     = cms2.mus_type()[mui];
                pt_       = cms2.mus_p4()[mui].pt();
                eta_      = cms2.mus_p4()[mui].eta();
                phi_      = cms2.mus_p4()[mui].phi();
                iso_      = muonIsoValue(mui);
                d0corr_   = cms2.mus_d0corr()[mui];

                int trkidx = cms2.mus_trkidx()[mui];
                d0vtx_     = cms2.trks_d0vtx()[trkidx];

                dphipfmet_ = deltaPhi(thePFMetPhi, cms2.mus_p4()[mui].phi());
                dphitcmet_ = deltaPhi(theTCMetPhi, cms2.mus_p4()[mui].phi());

                // clean jets for _this_ hyp lepton
                std::vector<unsigned int> theJetIndices;
                sumjetpt_ = 0.;
                for (unsigned int jeti = 0; jeti < cms2.pfjets_p4().size(); ++jeti) {
                    LorentzVector vjet = cms2.pfjets_p4()[jeti];

                    bool jetIsLep = false;
                    for (unsigned int muj = 0; muj < cms2.mus_p4().size(); ++muj) {
                        if (cms2.mus_p4()[muj].pt() <= 20.0) continue;
                        LorentzVector vlep = cms2.mus_p4()[muj];
                        if (dRbetweenVectors(vjet, vlep) < 0.4 && muonId(muj, NominalTTbarV2))
                            jetIsLep = true;
                    }
                    for (unsigned int elj = 0; elj < cms2.els_p4().size(); ++elj) {
                        LorentzVector vlep = cms2.els_p4()[elj];
                        if (dRbetweenVectors(vjet, vlep) < 0.4 && cms2.els_p4()[elj].pt() > 20. && pass_electronSelection(elj, electronSelection_ttbarV2))
                            jetIsLep = true;
                    }
                    if (jetIsLep) continue;

                    if (cms2.pfjets_p4()[jeti].pt() > 30. && fabs(cms2.pfjets_p4()[jeti].eta()) < 2.5 && isGoodPFJet(jeti)) {
						 theJetIndices.push_back(jeti);
						 sumjetpt_ += vjet.Pt();  
					}
                }
                std::sort(theJetIndices.begin(), theJetIndices.end(), sortByPFJetPt);

                njets_        = theJetIndices.size();
                jet1pt_       = theJetIndices.size() > 0 ? cms2.pfjets_p4()[theJetIndices[0]].pt()  : -999999.;
                jet1eta_      = theJetIndices.size() > 0 ? cms2.pfjets_p4()[theJetIndices[0]].eta() : -999999.;
                jet1phi_      = theJetIndices.size() > 0 ? cms2.pfjets_p4()[theJetIndices[0]].phi() : -999999.;
                jet2pt_       = theJetIndices.size() > 1 ? cms2.pfjets_p4()[theJetIndices[1]].pt()  : -999999.;
                jet2eta_      = theJetIndices.size() > 1 ? cms2.pfjets_p4()[theJetIndices[1]].eta() : -999999.;
                jet2phi_      = theJetIndices.size() > 1 ? cms2.pfjets_p4()[theJetIndices[1]].phi() : -999999.;
                jet3pt_       = theJetIndices.size() > 2 ? cms2.pfjets_p4()[theJetIndices[2]].pt()  : -999999.;
                jet3eta_      = theJetIndices.size() > 2 ? cms2.pfjets_p4()[theJetIndices[2]].eta() : -999999.;
                jet3phi_      = theJetIndices.size() > 2 ? cms2.pfjets_p4()[theJetIndices[2]].phi() : -999999.;

                LorentzVector dijetP4;
                jetmass_ = theJetIndices.size() > 1 ? sqrt((cms2.pfjets_p4()[theJetIndices[0]]+cms2.pfjets_p4()[theJetIndices[1]]).M2()) : -999999.; 

                // comment out initializes no longer needed now that call to InitBabyNtuple() moved inside loop
                double mindphipfmet = 999999.;
                double mindphitcmet = 999999.;
                neffbtags_   = 0;
                npurbtags_   = 0;
                ntceffbtags_ = 0;
                ntcpurbtags_ = 0;
                jet1isBtag_  = 0;
                jet2isBtag_  = 0;
                jet3isBtag_  = 0;
                for (unsigned int jeti = 0; jeti < theJetIndices.size(); ++jeti)
                {
                    if (cms2.pfjets_simpleSecondaryVertexHighEffBJetTag()[theJetIndices[jeti]] > 1.74)
                    {
                        ++neffbtags_;

                        if (jeti == 0)
                            jet1isBtag_ = 1;
                        else if (jeti == 1)
                            jet2isBtag_ = 1;
                        else if (jeti == 2)
                            jet3isBtag_ = 1;
                    }
                    if (cms2.pfjets_simpleSecondaryVertexHighPurBJetTags()[theJetIndices[jeti]] > 2.)
                    {
                        ++npurbtags_;

                        if (jeti == 0)
                            jet1isBtag_ = 1;
                        else if (jeti == 1)
                            jet2isBtag_ = 1;
                        else if (jeti == 2)
                            jet3isBtag_ = 1;
                    }
					if (cms2.pfjets_trackCountingHighEffBJetTag()[jeti] > 1.7)
						 ++ntceffbtags_;
					if (cms2.pfjets_trackCountingHighPurBJetTag()[jeti] > 1.19)
						 ++ntcpurbtags_;

                    float currdphipfmet = deltaPhi(thePFMetPhi, cms2.pfjets_p4()[theJetIndices[jeti]].phi());
                    if (currdphipfmet < mindphipfmet)
                        mindphipfmet = currdphipfmet;

                    float currdphitcmet = deltaPhi(theTCMetPhi, cms2.pfjets_p4()[theJetIndices[jeti]].phi());
                    if (currdphitcmet < mindphitcmet)
                        mindphitcmet = currdphitcmet;
                }

                dphipfmetjet_ = mindphipfmet;
                dphitcmetjet_ = mindphitcmet;

                float mindrjet = 999999.;
                for (unsigned int jeti = 0; jeti < theJetIndices.size(); ++jeti) {
                    float deta = cms2.mus_p4()[mui].eta()-cms2.pfjets_p4()[theJetIndices[jeti]].eta();
                    float dphi = deltaPhi(cms2.mus_p4()[mui].phi(), cms2.pfjets_p4()[theJetIndices[jeti]].phi());
                    float currdrjet = sqrt(deta*deta+dphi*dphi);
                    if (currdrjet < mindrjet)
                        mindrjet = currdrjet;
                }

                drjet_           = mindrjet;
                mt_              = sqrt(2.*pt_*pfmet_*(1.-cos(dphipfmet_)));
                tcmt_            = sqrt(2.*pt_*tcmet_*(1.-cos(dphitcmet_)));
                mu_muonidfull_   = muonId(mui, NominalTTbarV2);
                mu_muonid_       = muonIdNotIsolated(mui, NominalTTbarV2);
                mu_muonidfullV1_ = muonId(mui, NominalTTbar);
                mu_muonidV1_     = muonIdNotIsolated(mui, NominalTTbar);
                mu_goodmask_     = cms2.mus_goodmask()[mui];
                mu_gfitchi2_     = cms2.mus_gfit_chi2()[mui] < -9000. ? -999999. : cms2.mus_gfit_chi2()[mui]/cms2.mus_gfit_ndof()[mui];
                mu_cosmic_       = isCosmics(mui);
                mu_siHits_       = cms2.mus_validHits()[mui];
                mu_saHits_       = cms2.mus_gfit_validSTAHits()[mui];
                mu_emVetoDep_    = cms2.mus_iso_ecalvetoDep()[mui];
                mu_hadVetoDep_   = cms2.mus_iso_hcalvetoDep()[mui];
				if (cms2.mus_pfmusidx()[mui] > -1)
					 mu_isPFmuon_ = 1;

                // look for this guy in a same-flavor hyp and report
                // hyp_mass; if found in >1 hyp report one with high-
                // est mass
                for (unsigned int hypi = 0; hypi < cms2.hyp_p4().size(); ++hypi)
                {
                    // mumu
                    if (cms2.hyp_type()[hypi] == 0)
                        if ((unsigned int)cms2.hyp_lt_index()[hypi] == mui || (unsigned int)cms2.hyp_ll_index()[hypi] == mui)
                        {
                            float mass = cms2.hyp_p4()[hypi].mass2() > 0 ? cms2.hyp_p4()[hypi].mass() : TMath::Sqrt(-1 * cms2.hyp_p4()[hypi].mass2());
                            if (mass > sf_mass_)
                                sf_mass_ = mass;
                        }
                }

                FillBabyNtuple();
            }

            // electron stuff
            for (unsigned eli = 0; eli < cms2.els_p4().size(); ++eli)
            {
                // pt > 20
                if (cms2.els_p4()[eli].pt() <= 20.)
                    continue;

                // initialize baby quantities
                InitBabyNtuple();

                // event stuff
                strcpy(dataset_, cms2.evt_dataset().Data());
                run_        = cms2.evt_run();
                ls_         = cms2.evt_lumiBlock();
                evt_        = cms2.evt_event();
                isdata_     = cms2.evt_isRealData();
                nvtx_ = 0;
                for (unsigned int i = 0; i < cms2.vtxs_isFake().size(); i++) {
                    if (!cms2.vtxs_isFake()[i] && cms2.vtxs_isValid()[i])
                        ++nvtx_;
                }										

                if (!isdata_) {
                    int nlep  = leptonGenpCount_lepTauDecays(ngenels_, ngenmus_, ngentaus_);
                    scale1fb_ = cms2.evt_scale1fb();
                    pthat_    = cms2.genps_pthat();
                }

                pfmet_      = cms2.evt_pfmet();
                tcmet_      = cms2.evt_pf_tcmet();
                calotcmet_  = cms2.evt_tcmet();
                ntrks_      = cms2.trks_trk_p4().size();

                float thePFMetPhi     = cms2.evt_pfmetPhi();
                float theTCMetPhi     = cms2.evt_pf_tcmetPhi();
                float theCaloTCMetPhi = cms2.evt_tcmetPhi();

                float metx  = tcmet_ * cos(theTCMetPhi);
                float mety  = tcmet_ * sin(theTCMetPhi);
                float cmetx = calotcmet_ * cos(theCaloTCMetPhi);
                float cmety = calotcmet_ * sin(theCaloTCMetPhi);
                for (unsigned int muj = 0; muj < cms2.mus_p4().size(); ++muj) {
                    if (cms2.mus_p4()[muj].Pt() <= 10.0) continue;
                    if (!wasMetCorrectedForThisMuon(muj, usingTcMet) && muonIdNotIsolated(muj, NominalTTbarV2)) {
                        fixMetForThisMuon(muj, metx, mety, usingTcMet);
                        fixMetForThisMuon(muj, cmetx, cmety, usingTcMet);
                    }
                }
                tcmet_ = sqrt(metx * metx + mety * mety);
                theTCMetPhi = atan2(mety, metx);
                calotcmet_ = sqrt(cmetx * cmetx + cmety * cmety);
                theCaloTCMetPhi = atan2(cmety, cmetx);

                // loop over muons and electrons to get ngoodlep
                ngoodlep_ = 0;
                ngoodmus_ = 0;
                ngoodels_ = 0;
                for (unsigned m = 0; m < cms2.mus_p4().size(); ++m) {
                    if (cms2.mus_p4()[m].pt() <= 20.) continue;
                    if (muonId(m, NominalTTbarV2)) {
                        ++ngoodlep_;
                        ++ngoodmus_;
                    }
                }
                for (unsigned e = 0; e < cms2.els_p4().size(); ++e) {
                    if (cms2.els_p4()[e].pt() > 20. && pass_electronSelection(e, electronSelection_ttbarV2)) {
                        ++ngoodlep_;
                        ++ngoodels_;
                    }
                }

                eormu_    = 11*cms2.els_charge()[eli]*-1; // *-1 to follow pdg conventions
                type_     = cms2.els_type()[eli];
                pt_       = cms2.els_p4()[eli].pt();
                eta_      = cms2.els_p4()[eli].eta();
                phi_      = cms2.els_p4()[eli].phi();
                iso_      = electronIsolation_rel(eli, true);
                d0corr_   = cms2.els_d0corr()[eli];

                int trkidx= cms2.els_trkidx()[eli];
                if (trkidx >= 0)
                    d0vtx_= cms2.trks_d0vtx()[trkidx];

                dphipfmet_= deltaPhi(thePFMetPhi, cms2.els_p4()[eli].phi());
                dphitcmet_= deltaPhi(theTCMetPhi, cms2.els_p4()[eli].phi());

                // clean jets for _this_ hyp lepton
                std::vector<unsigned int> theJetIndices;
                sumjetpt_ = 0.0;
                for (unsigned int jeti = 0; jeti < cms2.pfjets_p4().size(); ++jeti) {
                    LorentzVector vjet = cms2.pfjets_p4()[jeti];

                    bool jetIsLep = false;
                    for (unsigned int muj = 0; muj < cms2.mus_p4().size(); ++muj) {
                        LorentzVector vlep = cms2.mus_p4()[muj];
                        if (dRbetweenVectors(vjet, vlep) < 0.4 && cms2.mus_p4()[muj].pt() > 20. && muonId(muj, NominalTTbarV2))
                            jetIsLep = true;
                    }
                    for (unsigned int elj = 0; elj < cms2.els_p4().size(); ++elj) {
                        LorentzVector vlep = cms2.els_p4()[elj];
                        if (dRbetweenVectors(vjet, vlep) < 0.4 && cms2.els_p4()[elj].pt() > 20. && pass_electronSelection(elj, electronSelection_ttbarV2))
                            jetIsLep = true;
                    }
                    if (jetIsLep) continue;

                    if (cms2.pfjets_p4()[jeti].pt() > 30. && fabs(cms2.pfjets_p4()[jeti].eta()) < 2.5 && isGoodPFJet(jeti)) {
                        theJetIndices.push_back(jeti);
						sumjetpt_ += vjet.Pt();  
                    }
                }
                std::sort(theJetIndices.begin(), theJetIndices.end(), sortByPFJetPt);

                njets_        = theJetIndices.size();
                jet1pt_       = theJetIndices.size() > 0 ? cms2.pfjets_p4()[theJetIndices[0]].pt()  : -999999.;
                jet1eta_      = theJetIndices.size() > 0 ? cms2.pfjets_p4()[theJetIndices[0]].eta() : -999999.;
                jet1phi_      = theJetIndices.size() > 0 ? cms2.pfjets_p4()[theJetIndices[0]].phi() : -999999.;
                jet2pt_       = theJetIndices.size() > 1 ? cms2.pfjets_p4()[theJetIndices[1]].pt()  : -999999.;
                jet2eta_      = theJetIndices.size() > 1 ? cms2.pfjets_p4()[theJetIndices[1]].eta() : -999999.;
                jet2phi_      = theJetIndices.size() > 1 ? cms2.pfjets_p4()[theJetIndices[1]].phi() : -999999.;
                jet3pt_       = theJetIndices.size() > 2 ? cms2.pfjets_p4()[theJetIndices[2]].pt()  : -999999.;
                jet3eta_      = theJetIndices.size() > 2 ? cms2.pfjets_p4()[theJetIndices[2]].eta() : -999999.;
                jet3phi_      = theJetIndices.size() > 2 ? cms2.pfjets_p4()[theJetIndices[2]].phi() : -999999.;

                LorentzVector dijetP4;
                jetmass_ = theJetIndices.size() > 1 ? sqrt((cms2.pfjets_p4()[theJetIndices[0]]+cms2.pfjets_p4()[theJetIndices[1]]).M2()) : -999999.; 

                // comment out initializes no longer needed now that call to InitBabyNtuple() moved inside loop
                double mindphipfmet = 999999.;
                double mindphitcmet = 999999.;
                neffbtags_   = 0;
                npurbtags_   = 0;
                ntceffbtags_ = 0;
                ntcpurbtags_ = 0;
                jet1isBtag_  = 0;
                jet2isBtag_  = 0;
                jet3isBtag_  = 0;
                for (unsigned int jeti = 0; jeti < theJetIndices.size(); ++jeti)
                {
                    if (cms2.pfjets_simpleSecondaryVertexHighEffBJetTag()[theJetIndices[jeti]] > 1.74)
                    {
                        ++neffbtags_;

                        if (jeti == 0)
                            jet1isBtag_ = 1;
                        else if (jeti == 1)
                            jet2isBtag_ = 1;
                        else if (jeti == 2)
                            jet3isBtag_ = 1;
                    }
                    if (cms2.pfjets_simpleSecondaryVertexHighPurBJetTags()[theJetIndices[jeti]] > 2.)
                    {
                        ++npurbtags_;

                        if (jeti == 0)
                            jet1isBtag_ = 1;
                        else if (jeti == 1)
                            jet2isBtag_ = 1;
                        else if (jeti == 2)
                            jet3isBtag_ = 1;
                    }
					if (cms2.pfjets_trackCountingHighEffBJetTag()[jeti] > 1.7)
						 ++ntceffbtags_;
					if (cms2.pfjets_trackCountingHighPurBJetTag()[jeti] > 1.19)
						 ++ntcpurbtags_;

                    float currdphipfmet = deltaPhi(thePFMetPhi, cms2.pfjets_p4()[theJetIndices[jeti]].phi());
                    if (currdphipfmet < mindphipfmet)
                        mindphipfmet = currdphipfmet;

                    float currdphitcmet = deltaPhi(theTCMetPhi, cms2.pfjets_p4()[theJetIndices[jeti]].phi());
                    if (currdphitcmet < mindphitcmet)
                        mindphitcmet = currdphitcmet;
                }

                dphipfmetjet_ = mindphipfmet;
                dphitcmetjet_ = mindphitcmet;

                float mindrjet = 999999.;
                for (unsigned int jeti = 0; jeti < theJetIndices.size(); ++jeti)
                {
                    float deta = cms2.els_p4()[eli].eta()-cms2.pfjets_p4()[theJetIndices[jeti]].eta();
                    float dphi = deltaPhi(cms2.els_p4()[eli].phi(), cms2.pfjets_p4()[theJetIndices[jeti]].phi());
                    float currdrjet = sqrt(deta*deta+dphi*dphi);
                    if (currdrjet < mindrjet)
                        mindrjet = currdrjet;
                }

                drjet_        = mindrjet;
                mt_           = sqrt(2.*pt_*pfmet_*(1.-cos(dphipfmet_)));
                tcmt_         = sqrt(2.*pt_*tcmet_*(1.-cos(dphitcmet_)));
                e_cand01full_ = pass_electronSelection(eli, electronSelection_ttbar);
                e_cand01_     = electronId_cand(eli, CAND_01);
                e_vbtf90full_ = pass_electronSelection(eli, electronSelection_ttbarV2);
                electronIdComponent_t answer_vbtf90 = electronId_VBTF(eli, VBTF_35X_90);
                e_vbtf90_     = (answer_vbtf90 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                electronIdComponent_t answer_vbtf85 = electronId_VBTF(eli, VBTF_35X_85);
                e_vbtf85_     = (answer_vbtf85 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                electronIdComponent_t answer_vbtf80 = electronId_VBTF(eli, VBTF_35X_80);
                e_vbtf80_     = (answer_vbtf80 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                electronIdComponent_t answer_vbtf70 = electronId_VBTF(eli, VBTF_35X_70);
                e_vbtf70_     = (answer_vbtf70 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                e_scet_       = cms2.els_eSC()[eli] / cosh(cms2.els_etaSC()[eli]);
                e_eopin_      = cms2.els_eOverPIn()[eli];
                e_hoe_        = cms2.els_hOverE()[eli];
                e_dphiin_     = cms2.els_dPhiIn()[eli];
                e_detain_     = cms2.els_dEtaIn()[eli];
                e_e25Me55_    = cms2.els_e2x5Max()[eli]/cms2.els_e5x5()[eli];
                e_sigieie_    = cms2.els_sigmaIEtaIEta()[eli];
                e_eMe55_      = cms2.els_eMax()[eli]/cms2.els_e5x5()[eli];
                e_nmHits_     = cms2.els_exp_innerlayers()[eli];
                e_dcot_       = cms2.els_conv_dcot()[eli];
                e_dist_       = cms2.els_conv_dist()[eli];
                e_drmu_       = cms2.els_closestMuon()[eli] < 0 ? -999999. : cms2.els_musdr()[eli];
                e_isspike_    = isSpikeElectron(eli);
                e_scCharge_   = cms2.els_sccharge()[eli];
                e_gsfCharge_  = cms2.els_trk_charge()[eli];
                e_ctfCharge_  = cms2.els_trkidx()[eli] > -1 ? cms2.trks_charge()[cms2.els_trkidx()[eli]] : -999999;
                // look for this guy in a same-flavor hyp and report
                // hyp_mass; if found in >1 hyp report one with high-
                // est mass
                for(unsigned int hypi = 0; hypi < cms2.hyp_p4().size(); ++hypi)
                {
                    // ee
                    if (cms2.hyp_type()[hypi] == 3)
                        if ((unsigned int)cms2.hyp_lt_index()[hypi] == eli || (unsigned int)cms2.hyp_ll_index()[hypi] == eli)
                        {
                            float mass = cms2.hyp_p4()[hypi].mass2() > 0 ? cms2.hyp_p4()[hypi].mass() : TMath::Sqrt(-1 * cms2.hyp_p4()[hypi].mass2());
                            if (mass > sf_mass_)
                                sf_mass_ = mass;
                        }
                }

                FillBabyNtuple();
            }
        }

    }

    if (nEventsChain != nEventsTotal)
    {
        std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    }

    CloseBabyNtuple();

}

void emubabymaker::InitBabyNtuple ()
{
    // event stuff
    memset(dataset_, '\0', 200);
    run_          = -999999;
    ls_           = -999999;
    evt_          = -999999;
    isdata_       = 1;
    scale1fb_     = -999999.;
    pthat_        = -999999.;
    nvtx_         = -999999;

    pfmet_        = -999999.;
    tcmet_        = -999999.;
    calotcmet_    = -999999.;	 
    ntrks_        = -999999;
    njets_        = -999999;
    jet1pt_       = -999999.;
    jet2pt_       = -999999.;
    jet3pt_       = -999999.;
    sumjetpt_     = -999999.;
    jet1eta_      = -999999.;
    jet2eta_      = -999999.;
    jet3eta_      = -999999.;
    jet1phi_      = -999999.;
    jet2phi_      = -999999.;
    jet3phi_      = -999999.;
    jetmass_      = -999999.;
    jet1isBtag_   = 0;
    jet2isBtag_   = 0;
    jet3isBtag_   = 0;
    dphipfmetjet_ = -999999.;
    dphitcmetjet_ = -999999.;
    neffbtags_    = -999999;
    npurbtags_    = -999999;
    ntceffbtags_  = -999999;
    ntcpurbtags_  = -999999;

    // lepton stuff
    ngoodlep_     = -999999;
    ngoodmus_     = -999999;
    ngoodels_     = -999999;
    ngenels_      = -999999;
    ngenmus_      = -999999;
    ngentaus_     = -999999;
    eormu_        = -999999;
    type_         = -999999;
    pt_           = -999999.;
    eta_          = -999999.;
    phi_          = -999999.;
    iso_          = -999999.;
    d0corr_       = -999999.;
    d0vtx_        = -999999.;
    dphipfmet_    = -999999.;
    dphitcmet_    = -999999.;
    drjet_        = -999999.;
    mt_           = -999999.;
    tcmt_         = -999999.;
    sf_mass_      = -999999.;

    // muon stuff
    mu_muonidfull_  = 0;
    mu_muonid_      = 0;
    mu_muonidfullV1_= 0;
    mu_muonidV1_    = 0;
    mu_goodmask_    = -999999;
    mu_gfitchi2_    = -999999.;
    mu_cosmic_      = 0;
    mu_siHits_      = -999999;
    mu_saHits_      = -999999;
    mu_emVetoDep_   = -999999.;
    mu_hadVetoDep_  = -999999.;
	mu_isPFmuon_    = 0;

    // electron stuff
    e_cand01full_      = 0;
    e_cand01_          = 0;
    e_vbtf90full_      = 0;
    e_vbtf90_          = 0;
    e_vbtf85_          = 0;
    e_vbtf80_          = 0;
    e_vbtf70_          = 0;
    e_scet_            = -999999.;
    e_eopin_           = -999999.;
    e_hoe_             = -999999.;
    e_dphiin_          = -999999.;
    e_detain_          = -999999.;
    e_e25Me55_         = -999999.;
    e_sigieie_         = -999999.;
    e_eMe55_           = -999999.;
    e_nmHits_          = -999999;
    e_dcot_            = -999999.;
    e_dist_            = -999999.;
    e_drmu_            = -999999.;
    e_isspike_         = 0;
    e_ctfCharge_       = -999999;
    e_gsfCharge_       = -999999;
    e_scCharge_        = -999999;
}

void emubabymaker::MakeBabyNtuple(const char *babyFilename)
{
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();

    babyFile_ = new TFile(Form("%s", babyFilename), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");

    // event stuff
    babyTree_->Branch("dataset",      &dataset_,      "dataset[200]/C");
    babyTree_->Branch("run",          &run_,          "run/I"         );
    babyTree_->Branch("ls",           &ls_,           "ls/I"          );
    babyTree_->Branch("evt",          &evt_,          "evt/I"         );
    babyTree_->Branch("nvtx",         &nvtx_,         "nvtx/I"        );
    babyTree_->Branch("isdata",       &isdata_,       "isdata/I"      );
    babyTree_->Branch("scale1fb",     &scale1fb_,     "scale1fb/F"    );
    babyTree_->Branch("pthat",        &pthat_,        "pthat/F"       );
    babyTree_->Branch("pfmet",        &pfmet_,        "pfmet/F"       );
    babyTree_->Branch("tcmet",        &tcmet_,        "tcmet/F"       );
    babyTree_->Branch("calotcmet",    &calotcmet_,    "calotcmet/F"   );
    babyTree_->Branch("ntrks",        &ntrks_,        "ntrks/I"       );
    babyTree_->Branch("njets",        &njets_,        "njets/I"       ); // uncorrected pt > 20
    babyTree_->Branch("jet1pt",       &jet1pt_,       "jet1pt/F"      );
    babyTree_->Branch("jet2pt",       &jet2pt_,       "jet2pt/F"      );
    babyTree_->Branch("jet3pt",       &jet3pt_,       "jet3pt/F"      );
    babyTree_->Branch("sumjetpt",     &sumjetpt_,     "sumjetpt/F"    );    
    babyTree_->Branch("jet1eta",      &jet1eta_,      "jet1eta/F"     );
    babyTree_->Branch("jet2eta",      &jet2eta_,      "jet2eta/F"     );
    babyTree_->Branch("jet3eta",      &jet3eta_,      "jet3eta/F"     );
    babyTree_->Branch("jet1phi",      &jet1phi_,      "jet1phi/F"     );
    babyTree_->Branch("jet2phi",      &jet2phi_,      "jet2phi/F"     );
    babyTree_->Branch("jet3phi",      &jet3phi_,      "jet3phi/F"     );
    babyTree_->Branch("jetmass",      &jetmass_,      "jetmass/F"     );
    babyTree_->Branch("jet1isBtag",   &jet1isBtag_,   "jet1isBtag/O"  );
    babyTree_->Branch("jet2isBtag",   &jet2isBtag_,   "jet2isBtag/O"  );
    babyTree_->Branch("jet3isBtag",   &jet3isBtag_,   "jet3isBtag/O"  );
    babyTree_->Branch("dphipfmetjet", &dphipfmetjet_, "dphipfmetjet/F");
    babyTree_->Branch("dphitcmetjet", &dphitcmetjet_, "dphitcmetjet/F");
    babyTree_->Branch("neffbtags",    &neffbtags_,    "neffbtags/I"   );
    babyTree_->Branch("npurbtags",    &npurbtags_,    "npurbtags/I"   );
    babyTree_->Branch("ntceffbtags",  &ntceffbtags_,  "ntceffbtags/I" );
    babyTree_->Branch("ntcpurbtags",  &ntcpurbtags_,  "ntcpurbtags/I" );

    // lepton stuff
    babyTree_->Branch("ngoodlep",  &ngoodlep_,  "ngoodlep/I" );
    babyTree_->Branch("ngoodmus",  &ngoodmus_,  "ngoodmus/I" );
    babyTree_->Branch("ngoodels",  &ngoodels_,  "ngoodels/I" );
    babyTree_->Branch("ngenels",   &ngenels_,   "ngenels/I" );
    babyTree_->Branch("ngenmus",   &ngenmus_,   "ngenmus/I" );
    babyTree_->Branch("ngentaus",  &ngentaus_,  "ngentaus/I" );
    babyTree_->Branch("eormu",     &eormu_,     "eormu/I"    );
    babyTree_->Branch("type",      &type_,      "type/I"     );
    babyTree_->Branch("pt",        &pt_,        "pt/F"       );
    babyTree_->Branch("eta",       &eta_,       "eta/F"      );
    babyTree_->Branch("phi",       &phi_,       "phi/F"      );
    babyTree_->Branch("iso",       &iso_,       "iso/F"      );
    babyTree_->Branch("d0corr",    &d0corr_,    "d0corr/F"   );
    babyTree_->Branch("d0vtx",     &d0vtx_,     "d0vtx/F"    );
    babyTree_->Branch("dphipfmet", &dphipfmet_, "dphipfmet/F");
    babyTree_->Branch("dphitcmet", &dphitcmet_, "dphitcmet/F");
    babyTree_->Branch("drjet",     &drjet_,     "drjet/F"    );
    babyTree_->Branch("mt",        &mt_,        "mt/F"       );
    babyTree_->Branch("tcmt",      &tcmt_,      "tcmt/F"     );
    babyTree_->Branch("sf_mass",   &sf_mass_,   "sf_mass/F"  );

    // muon stuff
    babyTree_->Branch("mu_muonidfull",   &mu_muonidfull_,   "mu_muonidfull/O"  );
    babyTree_->Branch("mu_muonid",       &mu_muonid_,       "mu_muonid/O"      );
    babyTree_->Branch("mu_muonidfullV1", &mu_muonidfullV1_, "mu_muonidfullV1/O");
    babyTree_->Branch("mu_muonidV1",     &mu_muonidV1_,     "mu_muonidV1/O"    );
    babyTree_->Branch("mu_goodmask",     &mu_goodmask_,     "mu_goodmask/I"    );
    babyTree_->Branch("mu_gfitchi2",     &mu_gfitchi2_,     "mu_gfitchi2/F"    );
    babyTree_->Branch("mu_cosmic",       &mu_cosmic_,       "mu_cosmic/O"      );
    babyTree_->Branch("mu_siHits",       &mu_siHits_,       "mu_siHits/I"      );
    babyTree_->Branch("mu_saHits",       &mu_saHits_,       "mu_saHits/I"      );
    babyTree_->Branch("mu_emVetoDep",    &mu_emVetoDep_,    "mu_emVetoDep/F"   );
    babyTree_->Branch("mu_hadVetoDep",   &mu_hadVetoDep_,   "mu_hadVetoDep/F"  );
	babyTree_->Branch("mu_isPFmuon",     &mu_isPFmuon_,     "mu_isPFmuon/O"    );

    // electron stuff
    babyTree_->Branch("e_cand01full",      &e_cand01full_,      "e_cand01full/O"     );
    babyTree_->Branch("e_cand01",          &e_cand01_,          "e_cand01/O"         );
    babyTree_->Branch("e_vbtf90full",      &e_vbtf90full_,      "e_vbtf90full/O"     );
    babyTree_->Branch("e_vbtf90",          &e_vbtf90_,          "e_vbtf90/O"         );
    babyTree_->Branch("e_vbtf85",          &e_vbtf85_,          "e_vbtf85/O"         );
    babyTree_->Branch("e_vbtf80",          &e_vbtf80_,          "e_vbtf80/O"         );
    babyTree_->Branch("e_vbtf70",          &e_vbtf70_,          "e_vbtf70/O"         );
    babyTree_->Branch("e_scet",            &e_scet_,            "e_scet/F"           );
    babyTree_->Branch("e_eopin",           &e_eopin_,           "e_eopin/F"          );
    babyTree_->Branch("e_hoe",             &e_hoe_,             "e_hoe/F"            );
    babyTree_->Branch("e_dphiin",          &e_dphiin_,          "e_dphiin/F"         );
    babyTree_->Branch("e_detain",          &e_detain_,          "e_detain/F"         );
    babyTree_->Branch("e_e25Me55",         &e_e25Me55_,         "e_e25Me55/F"        );
    babyTree_->Branch("e_sigieie",         &e_sigieie_,         "e_sigieie/F"        );
    babyTree_->Branch("e_eMe55",           &e_eMe55_,           "e_eMe55/F"          ); // for spikes
    babyTree_->Branch("e_nmHits",          &e_nmHits_,          "e_nmHits/I"         ); // els_exp_innerlayers
    babyTree_->Branch("e_dcot",            &e_dcot_,            "e_dcot/F"           ); // els_conv_dcot
    babyTree_->Branch("e_dist",            &e_dist_,            "e_dist/F"           ); // els_conv_dist
    babyTree_->Branch("e_drmu",            &e_drmu_,            "e_drmu/F"           );
    babyTree_->Branch("e_isspike",         &e_isspike_,         "e_isspike/O"        ); // isSpikeElectron
    babyTree_->Branch("e_ctfCharge",       &e_ctfCharge_,       "e_ctfCharge/I"      );
    babyTree_->Branch("e_gsfCharge",       &e_gsfCharge_,       "e_gsfCharge/I"      );
    babyTree_->Branch("e_scCharge",        &e_scCharge_,        "e_scCharge/I"       );
}

void emubabymaker::FillBabyNtuple()
{
    babyTree_->Fill();
}

void emubabymaker::CloseBabyNtuple()
{
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();
}

emubabymaker::~emubabymaker() {
    std::cout << "emubabymaker::~emubabymaker" << std::endl;
    delete babyFile_;
    delete babyTree_;
};


