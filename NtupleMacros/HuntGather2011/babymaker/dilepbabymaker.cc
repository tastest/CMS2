#include "babymakercommon.h"
#include "dilepbabymaker.h" 
#include "readTriggerList.h"

#include <algorithm>
#include <iostream>

#include "TChain.h"
#include "TCollection.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TTreeCache.h"

// CMS2 CORE includes 
#include "CORE/MT2/MT2.h"
#include "CORE/CMS2.h"
#include "CORE/mcSelections.h"
#include "CORE/electronSelections.h"
#include "CORE/electronSelectionsParameters.h"
#include "CORE/metSelections.h"
#include "CORE/muonSelections.h"
#include "CORE/trackSelections.h"
#include "CORE/ssSelections.h"
#include "CORE/eventSelections.h"
#include "CORE/triggerUtils.h"

// tools include is for
// the duplicate cleaning
#include "Tools/tools.cc"
#include "Tools/ElectronIDMVA.cc"

void dilepbabymaker::NewRun()
{
    // if the run has changed
    // then update the list of trigger to check
    triggers_e_ = get_trigger_names(cms2.evt_run(), "e");
    triggers_m_ = get_trigger_names(cms2.evt_run(), "m");
    triggers_ee_ = get_trigger_names(cms2.evt_run(), "ee");
    triggers_mm_ = get_trigger_names(cms2.evt_run(), "mm");
    triggers_em_ = get_trigger_names(cms2.evt_run(), "em");
    triggers_ehad_ee_ = get_trigger_names(cms2.evt_run(), "eehad");
    triggers_mhad_mm_ = get_trigger_names(cms2.evt_run(), "mmhad");
    triggers_mhad_em_ = get_trigger_names(cms2.evt_run(), "emhad");
    //std::cout << "[dilepbabymaker::NewRun] " << cms2.evt_run() << std::endl;
    //for (unsigned int i = 0; i < triggers_m_.size(); ++i) std::cout << triggers_m_[i] << std::endl;

}

void dilepbabymaker::ScanChain (const char *inputFilename, const char *babyFilename, int nEvents)
{

    // set trigger file
    std::string cms2_location = "..";
    char* ppath = getenv("CMS2_LOCATION");
    if (ppath != NULL) {
        cms2_location = ppath;
        cms2_location += "/NtupleMacros/HuntGather2011";
    }
    set_trigger_file((cms2_location + "/runlists/trigger.txt").c_str());

    TChain *chain = new TChain("Events");
    chain->Add(inputFilename);
    TObjArray *listOfFiles = chain->GetListOfFiles();    

    unsigned int nEventsChain=0;
    if (nEvents==-1) 
        nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;

    // make a baby ntuple
    MakeBabyNtuple(babyFilename);

    // make a random number generator
    TRandom3 rndm;

    // file loop
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;
    while ((currentFile = (TFile*)fileIter.Next()))
    {
        TFile* f = TFile::Open(currentFile->GetTitle());
        TTree *tree = (TTree*)f->Get("Events");
        TTreeCache::SetLearnEntries(10);
        tree->SetCacheSize(128*1024*1024);
        cms2.Init(tree);

        //Event Loop
        unsigned int thisRun = 0;
        unsigned int nEvents = tree->GetEntries();
        for(unsigned int event = 0; event < nEvents; ++event)
        {
            tree->LoadTree(event);
            cms2.GetEntry(event);
            ++nEventsTotal;

            // call progress counted
            CMS2::progress(nEventsTotal, nEventsChain );
            //if (event == 0) PrintTriggers();

            // every time we find a new run,
            // update trigger monitoring from configuration
            if (cms2.evt_run() != thisRun) {
                NewRun();
                thisRun = cms2.evt_run();
            } 

            // dilepton hypothesis stuff
            for (unsigned hypi = 0; hypi < cms2.hyp_p4().size(); ++hypi)
            {

                //
                // select a basic hypothesis (5 GeV cut on muon pT, 10 GeV cut on electron pT)
                //
                int index1 = cms2.hyp_lt_index()[hypi];
                int index2 = cms2.hyp_ll_index()[hypi];

                if (abs(cms2.hyp_lt_id()[hypi]) == 13 && !(cms2.mus_type()[index1] & 6))
                    continue;
                if (abs(cms2.hyp_ll_id()[hypi]) == 13 && !(cms2.mus_type()[index2] & 6))
                    continue;

                //
                // initialize baby quantities
                //

                InitBabyNtuple();

                //
                // Fill event level information
                //

                SetEventLevelInfo();
                rndm_ = rndm.Uniform();

                //
                // Fill mc id of hypothesis leptons if this is MC
                //

                int lt_id = cms2.hyp_lt_id()[hypi];
                int ll_id = cms2.hyp_ll_id()[hypi];

                if (!isdata_) {
                    int lt_idx = cms2.hyp_lt_index()[hypi];
                    int ll_idx = cms2.hyp_ll_index()[hypi];
                    mcid1_ = abs(lt_id) == 13 ? cms2.mus_mc_id()[lt_idx] : cms2.els_mc_id()[lt_idx];
                    mcmotherid1_ = abs(lt_id) == 13 ? cms2.mus_mc_motherid()[lt_idx] : cms2.els_mc_motherid()[lt_idx];
                    mc3id1_ = abs(lt_id) == 13 ? cms2.mus_mc3_id()[lt_idx] : cms2.els_mc3_id()[lt_idx];
                    mc3motherid1_ = abs(lt_id) == 13 ? cms2.mus_mc3_motherid()[lt_idx] : cms2.els_mc3_motherid()[lt_idx];
                    lep1isFromW_ = leptonIsFromW(lt_idx, lt_id, true);
                    mcpt1_ = abs(lt_id) == 13 ? cms2.mus_mc_p4()[lt_idx].pt() : cms2.els_mc_p4()[lt_idx].pt();
                    mcmotherpt1_ = abs(lt_id) == 13 ? cms2.mus_mc_motherp4()[lt_idx].pt() : cms2.els_mc_motherp4()[lt_idx].pt();
                    mc3pt1_ = abs(lt_id) == 13 ? cms2.genps_p4()[cms2.mus_mc3idx()[lt_idx]].pt() : cms2.genps_p4()[cms2.els_mc3idx()[lt_idx]].pt();
                    mc3motherpt1_ = abs(lt_id) == 13 ? cms2.genps_p4()[cms2.mus_mc3_motheridx()[lt_idx]].pt() : cms2.genps_p4()[cms2.els_mc3_motheridx()[lt_idx]].pt();


                    mcid2_ = abs(ll_id) == 13 ? cms2.mus_mc_id()[ll_idx] : cms2.els_mc_id()[ll_idx];
                    mcmotherid2_ = abs(ll_id) == 13 ? cms2.mus_mc_motherid()[ll_idx] : cms2.els_mc_motherid()[ll_idx];
                    mc3id2_ = abs(ll_id) == 13 ? cms2.mus_mc3_id()[ll_idx] : cms2.els_mc3_id()[ll_idx];
                    mc3motherid2_ = abs(ll_id) == 13 ? cms2.mus_mc3_motherid()[ll_idx] : cms2.els_mc3_motherid()[ll_idx];
                    lep2isFromW_ = leptonIsFromW(ll_idx, ll_id, true);
                    mcpt2_ = abs(ll_id) == 13 ? cms2.mus_mc_p4()[ll_idx].pt() : cms2.els_mc_p4()[ll_idx].pt();
                    mcmotherpt2_ = abs(ll_id) == 13 ? cms2.mus_mc_motherp4()[ll_idx].pt() : cms2.els_mc_motherp4()[ll_idx].pt();
                    mc3pt2_ = abs(ll_id) == 13 ? cms2.genps_p4()[cms2.mus_mc3idx()[ll_idx]].pt() : cms2.genps_p4()[cms2.els_mc3idx()[ll_idx]].pt();
                    mc3motherpt2_ = abs(ll_id) == 13 ? cms2.genps_p4()[cms2.mus_mc3_motheridx()[ll_idx]].pt() : cms2.genps_p4()[cms2.els_mc3_motheridx()[ll_idx]].pt();
                }

                //
                // Fill trigger information
                //

                if (isdata_) {

                    // do electrons (single)
                    if (abs(lt_id) == 11) {
                        PassTriggerGroup(triggers_e_, pass_trg_single_e1_);
                        PassTriggerGroup(triggers_e_, cms2.hyp_lt_p4()[hypi], match_trg_single_e1_);
                    }
                    if (abs(ll_id) == 11) {
                        PassTriggerGroup(triggers_e_, pass_trg_single_e2_);
                        PassTriggerGroup(triggers_e_, cms2.hyp_ll_p4()[hypi], match_trg_single_e2_);
                    }

                    // do muons (single)
                    if (abs(lt_id) == 13) {
                        PassTriggerGroup(triggers_m_, pass_trg_single_mu1_);
                        PassTriggerGroup(triggers_m_, cms2.hyp_lt_p4()[hypi], match_trg_single_mu1_);
                    }
                    if (abs(ll_id) == 13) {
                        PassTriggerGroup(triggers_m_, pass_trg_single_mu2_);
                        PassTriggerGroup(triggers_m_, cms2.hyp_ll_p4()[hypi], match_trg_single_mu2_);
                    }

                    // do double electron
                    if (abs(lt_id) == 11 && abs(ll_id) == 11) {
                        PassTriggerGroup(triggers_ee_, pass_trg_double_e1_);
                        PassTriggerGroup(triggers_ee_, cms2.hyp_lt_p4()[hypi], match_trg_double_e1_);

                        PassTriggerGroup(triggers_ee_, pass_trg_double_e2_);
                        PassTriggerGroup(triggers_ee_, cms2.hyp_ll_p4()[hypi], match_trg_double_e2_);

                        // for now, don't require matching of leptons and trigger legs due to bug in CMS2
                        // however, once the bug is fixed, revert back to matching
                        PassTriggerGroup(triggers_ehad_ee_, pass_trg_had_double_e1_);
                        PassTriggerGroup(triggers_ehad_ee_, cms2.hyp_lt_p4()[hypi], match_trg_had_double_e1_);

                        PassTriggerGroup(triggers_ehad_ee_, pass_trg_had_double_e2_);
                        PassTriggerGroup(triggers_ehad_ee_, cms2.hyp_ll_p4()[hypi], match_trg_had_double_e2_);
                    }

                    // do double muon
                    if (abs(lt_id) == 13 && abs(ll_id) == 13) {
                        PassTriggerGroup(triggers_mm_, pass_trg_double_mu1_);
                        PassTriggerGroup(triggers_mm_, cms2.hyp_lt_p4()[hypi], match_trg_double_mu1_);

                        PassTriggerGroup(triggers_mm_, pass_trg_double_mu2_);
                        PassTriggerGroup(triggers_mm_, cms2.hyp_ll_p4()[hypi], match_trg_double_mu2_);

                        // for now, don't require matching of leptons and trigger legs due to bug in CMS2
                        // however, once the bug is fixed, revert back to matching
                        PassTriggerGroup(triggers_mhad_mm_, pass_trg_had_double_mu1_);
                        PassTriggerGroup(triggers_mhad_mm_, cms2.hyp_lt_p4()[hypi], match_trg_had_double_mu1_);

                        PassTriggerGroup(triggers_mhad_mm_, pass_trg_had_double_mu2_);
                        PassTriggerGroup(triggers_mhad_mm_, cms2.hyp_ll_p4()[hypi], match_trg_had_double_mu2_);

                    }

                    // do cross trigger (e-mu)
                    if (abs(lt_id) != abs(ll_id)) {
                        PassTriggerGroup(triggers_em_, pass_trg_cross_emu_);
                        PassTriggerGroup(triggers_mhad_em_, pass_trg_had_cross_emu_);
                    }



                    // integrated luminosity per luminosity section
                    if (isdata_)
                        intLumiPerLS_ = cms2.ls_lumiSectionLength() * cms2.ls_avgInsRecLumi();

                } else {
                    pass_trg_single_e1_      = 1;
                    pass_trg_single_e2_      = 1;
                    pass_trg_single_mu1_     = 1;
                    pass_trg_single_mu2_     = 1;
                    pass_trg_double_e1_      = 1;
                    pass_trg_double_e2_      = 1;
                    pass_trg_double_mu1_     = 1;
                    pass_trg_double_mu2_     = 1;
                    pass_trg_cross_emu_      = 1;
                    pass_trg_had_double_e1_  = 1;
                    pass_trg_had_double_e2_  = 1;
                    pass_trg_had_double_mu1_ = 1;
                    pass_trg_had_double_mu2_ = 1;
                    pass_trg_had_cross_emu_  = 1;
                    match_trg_single_e1_      = 1;
                    match_trg_single_e2_      = 1;
                    match_trg_single_mu1_     = 1;
                    match_trg_single_mu2_     = 1;
                    match_trg_double_e1_      = 1;
                    match_trg_double_e2_      = 1;
                    match_trg_double_mu1_     = 1;
                    match_trg_double_mu2_     = 1;
                    match_trg_cross_emu_      = 1;
                    match_trg_had_double_e1_  = 1;
                    match_trg_had_double_e2_  = 1;
                    match_trg_had_double_mu1_ = 1;
                    match_trg_had_double_mu2_ = 1;
                    match_trg_had_cross_emu_  = 1;
                }

                // 
                // Fill MET information
                //

                pfmet_      = cms2.evt_pfmet();

                float theTCMetPhi     = 0.0;
                if (cms2.evt_CMS2tag() != "V03-06-14" && cms2.evt_dataset() != "/Mu/Run2010B-Nov4ReReco_v1/RECO") {
                    tcmet_          = cms2.evt_pf_tcmet();
                    theTCMetPhi     = cms2.evt_pf_tcmetPhi();
                } else {
                    tcmet_          = cms2.evt_tcmet();
                    theTCMetPhi     = cms2.evt_tcmetPhi();
                }
                calotcmet_  = cms2.evt_tcmet();

                float thePFMetPhi     = cms2.evt_pfmetPhi();
                float theCaloTCMetPhi = cms2.evt_tcmetPhi();
                pfmetphi_ = thePFMetPhi;
                tcmetphi_ = theCaloTCMetPhi;

                float metx   = tcmet_ * cos(theTCMetPhi);
                float mety   = tcmet_ * sin(theTCMetPhi);
                float cmetx  = calotcmet_ * cos(theCaloTCMetPhi);
                float cmety  = calotcmet_ * sin(theCaloTCMetPhi);
                float pfmetx = pfmet_ * cos(thePFMetPhi);
                float pfmety = pfmet_ * sin(thePFMetPhi);
                for (unsigned int muj = 0; muj < cms2.mus_p4().size(); ++muj) {
                    if (cms2.mus_p4()[muj].Pt() <= 10.0) continue;
                    if (!wasMetCorrectedForThisMuon(muj, usingTcMet) && muonIdNotIsolated(muj, OSGeneric_v3)) {
                        fixMetForThisMuon(muj, metx, mety, usingTcMet);
                        fixMetForThisMuon(muj, cmetx, cmety, usingTcMet);
                    }
                }
                tcmet_          = sqrt(metx * metx + mety * mety);
                theTCMetPhi     = atan2(mety, metx);
                calotcmet_      = sqrt(cmetx * cmetx + cmety * cmety);
                theCaloTCMetPhi = atan2(cmety, cmetx);

                //
                // Fill lepton counters
                //

                //
                // note that this is also used to define what leptons the jets
                // will be cleaned for
                //

                ngoodlep_   = 0;
                ngoodmus_   = 0;
                ngoodels_   = 0;
                ngoodlepSS_ = 0;
                ngoodmusSS_ = 0;
                ngoodelsSS_ = 0;

                std::vector<unsigned int> goodElectronIndicesTTBarV2;
                std::vector<unsigned int> goodMuonIndicesTTBarV2;
                std::vector<unsigned int> goodElectronIndicesSSV2;
                std::vector<unsigned int> goodMuonIndicesSSV2;
                std::vector<unsigned int> goodElectronIndicesSShighPt;

                // muon loop
                for (unsigned muii = 0; muii < cms2.mus_p4().size(); ++muii) {
                    // for SS
                    if (cms2.mus_p4()[muii].pt() > 5. && muonId(muii, NominalSSv4) && fabs(cms2.mus_p4()[muii].eta()) < 2.4) {
                        goodMuonIndicesSSV2.push_back(muii);
                        ++ngoodlepSS_;
                        ++ngoodmusSS_;
                    }
                    // for TTBarV2
                    if (cms2.mus_p4()[muii].pt() > 10. && muonId(muii, OSGeneric_v3)) {
                        goodMuonIndicesTTBarV2.push_back(muii);
                        ++ngoodlep_;
                        ++ngoodmus_;
                    }
                }

                // electron loop
                for (unsigned eli = 0; eli < cms2.els_p4().size(); ++eli) {
                    cuts_t cuts_passed = electronSelection(eli);
                    // for SS
                    if (cms2.els_p4()[eli].pt() > 10. && pass_electronSelectionCompareMask(cuts_passed, electronSelection_ssV5) && fabs(cms2.els_p4()[eli].eta()) < 2.4) {
                        goodElectronIndicesSSV2.push_back(eli);
                        ++ngoodlepSS_;
                        ++ngoodelsSS_;
                        if (cms2.els_ecalIso()[eli]/cms2.els_p4()[eli].pt() < 0.2) {
                            goodElectronIndicesSShighPt.push_back(eli);
                        }
                    }
                    // for TTBarV2
                    if(cms2.els_p4()[eli].pt() > 10. && pass_electronSelectionCompareMask(cuts_passed, electronSelection_el_OSV3)) {
                        goodElectronIndicesTTBarV2.push_back(eli);
                        ++ngoodlep_;
                        ++ngoodels_;
                    }
                }

                //
                // hypothesis stuff
                //

                hyp_type_    = cms2.hyp_type()[hypi];
                pt1_         = cms2.hyp_lt_p4()[hypi].pt();
                pt2_         = cms2.hyp_ll_p4()[hypi].pt();
                eta1_        = cms2.hyp_lt_p4()[hypi].eta();
                eta2_        = cms2.hyp_ll_p4()[hypi].eta();
                phi1_        = cms2.hyp_lt_p4()[hypi].phi();
                phi2_        = cms2.hyp_ll_p4()[hypi].phi();
                d0corr1_     = cms2.hyp_lt_d0corr()[hypi];
                d0corr2_     = cms2.hyp_ll_d0corr()[hypi];		 
                eormu1_      = cms2.hyp_lt_id()[hypi];
                eormu2_      = cms2.hyp_ll_id()[hypi];
                dphipfmet1_  = deltaPhi(cms2.hyp_lt_p4()[hypi].phi(), thePFMetPhi);
                dphipfmet2_  = deltaPhi(cms2.hyp_ll_p4()[hypi].phi(), thePFMetPhi);
                dphitcmet1_  = deltaPhi(cms2.hyp_lt_p4()[hypi].phi(), theTCMetPhi);
                dphitcmet2_  = deltaPhi(cms2.hyp_ll_p4()[hypi].phi(), theTCMetPhi);
                mass_        = cms2.hyp_p4()[hypi].mass2() > 0 ? cms2.hyp_p4()[hypi].mass() : TMath::Sqrt(-1 * cms2.hyp_p4()[hypi].mass2());
                dilpt_       = cms2.hyp_p4()[hypi].pt();		 
                dileta_      = cms2.hyp_p4()[hypi].eta();
                dilphi_      = cms2.hyp_p4()[hypi].phi();
                deltaphi_    = deltaPhi(phi1_, phi2_);

                // calculate projected pfmet
                float mindphi_pfmet = min(dphipfmet1_, dphipfmet2_);
                if (mindphi_pfmet < TMath::Pi()/2.)
                    proj_pfmet_ = pfmet_ * sin(mindphi_pfmet);
                else
                    proj_pfmet_ = pfmet_;
                // calculate projected pfmet
                float mindphi_tcmet = min(dphitcmet1_, dphitcmet2_);
                if (mindphi_tcmet < TMath::Pi()/2.)
                    proj_tcmet_ = tcmet_ * sin(mindphi_tcmet);
                else
                    proj_tcmet_ = tcmet_;				

                // initialize meff to 0
                pfmeff_ = 0.;
                tcmeff_ = 0.;

                // add MET to meff
                pfmeff_ += pfmet_;
                tcmeff_ += tcmet_;

                // add lepton pt to meff
                pfmeff_ += (pt1_ + pt2_);
                tcmeff_ += (pt1_ + pt2_);

                //
                // check if leptons are from the same vertex
                //
                int vertex_index = firstGoodVertex();
                lepsFromSameVtx_ = hypsFromFirstGoodVertex(hypi);

                //
                // clean jets for ALL leptons
                //

                VofP4 theJets;
                std::vector<unsigned int> theJetIndices;
                std::vector<unsigned int> theBtagJetIndices;
                sumjetpt_         = 0.;
                sumjetpt25_       = 0.;
                sumjetptSS30_     = 0.;
                sumjetptSS_       = 0.;
                sumjetptSShighPt_ = 0.;
                njets25_          = 0;
                njetsSS30_        = 0;
                njetsSS_          = 0;
                njetsSShighPt_    = 0;
                njets_            = 0;

                for (unsigned int jeti = 0; jeti < cms2.pfjets_p4().size(); ++jeti)
                {

                    LorentzVector vjet = cms2.pfjets_p4()[jeti]*cms2.pfjets_corL1FastL2L3()[jeti];
                    bool jetIsLepTTBarV2 = false;
                    bool jetIsLepSSV2 = false;
                    bool jetIsLepSShighPt = false;
                    
                    // we previously weren't performing jet-lepton overlap removal with the hypothesis leptons
                    // but we should have been; fixing this now
                    if (dRbetweenVectors(vjet, cms2.hyp_lt_p4()[hypi]) < 0.4)
                        continue;
                    if (dRbetweenVectors(vjet, cms2.hyp_ll_p4()[hypi]) < 0.4)
                        continue;

                    // check if jet is a lepton that passes the ttbar selection
                    for (unsigned int muj = 0; muj < goodMuonIndicesTTBarV2.size(); ++muj) {
                        LorentzVector vlep = cms2.mus_p4()[goodMuonIndicesTTBarV2[muj]];
                        if (!jetIsLepTTBarV2 && dRbetweenVectors(vjet, vlep) < 0.4)
                            jetIsLepTTBarV2 = true;
                    }
                    for (unsigned int elj = 0; elj < goodElectronIndicesTTBarV2.size(); ++elj) {
                        LorentzVector vlep = cms2.els_p4()[goodElectronIndicesTTBarV2[elj]];
                        if (!jetIsLepTTBarV2 && dRbetweenVectors(vjet, vlep) < 0.4)
                            jetIsLepTTBarV2 = true;
                    }

                    // 30 GeV TTBarV2 cleaning
                    if (!jetIsLepTTBarV2 && vjet.Pt() > 30. && fabs(cms2.pfjets_p4()[jeti].eta()) < 2.5 && isGoodPFJet(jeti)) {
                        sumjetpt_ += vjet.Pt();
                        njets_++;
                    }
                    // 25 GeV TTBarV2 cleaning
                    if (!jetIsLepTTBarV2 && vjet.Pt() > 25. && fabs(cms2.pfjets_p4()[jeti].eta()) < 2.5 && isGoodPFJet(jeti)) {
                        sumjetpt25_ += vjet.Pt();
                        njets25_ ++;
                    }

                    // check if jet is a lepton that passes the ss selection
                    for (unsigned int muj = 0; muj < goodMuonIndicesSSV2.size(); ++muj) {
                        LorentzVector vlep = cms2.mus_p4()[goodMuonIndicesSSV2[muj]];
                        if (!jetIsLepSSV2 && dRbetweenVectors(vjet, vlep) < 0.4)
                            jetIsLepSSV2 = true;
                    }
                    for (unsigned int elj = 0; elj < goodElectronIndicesSSV2.size(); ++elj) {
                        LorentzVector vlep = cms2.els_p4()[goodElectronIndicesSSV2[elj]];
                        if (!jetIsLepSSV2 && dRbetweenVectors(vjet, vlep) < 0.4)
                            jetIsLepSSV2 = true;
                    }
                    if (!jetIsLepSSV2 && fabs(cms2.pfjets_p4()[jeti].eta()) < 2.5 && isGoodPFJet(jeti)) {
                        if (vjet.Pt() > 20.)
                            theBtagJetIndices.push_back(jeti);
                        if (vjet.Pt() > 30.) {
                            theJets.push_back(cms2.pfjets_p4()[jeti]);
                            theJetIndices.push_back(jeti);
                            sumjetptSS30_ += vjet.Pt();
                        }
                    }
                    if (!jetIsLepSSV2 && vjet.Pt() > 40. && fabs(cms2.pfjets_p4()[jeti].eta()) < 2.5 && isGoodPFJet(jeti)) {
                        sumjetptSS_ += vjet.Pt();
                        njetsSS_++;
                    }


                    // check if jet is a lepton that passes the ss selection with offline isolation mockup
                    for (unsigned int muj = 0; muj < goodMuonIndicesSSV2.size(); ++muj) {
                        LorentzVector vlep = cms2.mus_p4()[goodMuonIndicesSSV2[muj]];
                        if (!jetIsLepSShighPt && dRbetweenVectors(vjet, vlep) < 0.4)
                            jetIsLepSShighPt = true;
                    }
                    for (unsigned int elj = 0; elj < goodElectronIndicesSShighPt.size(); ++elj) {
                        LorentzVector vlep = cms2.els_p4()[goodElectronIndicesSShighPt[elj]];
                        if (!jetIsLepSShighPt && dRbetweenVectors(vjet, vlep) < 0.4)
                            jetIsLepSShighPt = true;
                    }
                    if (!jetIsLepSShighPt && vjet.Pt() > 40. && fabs(cms2.pfjets_p4()[jeti].eta()) < 2.5 && isGoodPFJet(jeti)) {
                        sumjetptSShighPt_ += vjet.Pt();
                        njetsSShighPt_++;
                    }

                }
                std::sort(theJets.begin(), theJets.end(), sortByPt);
                std::sort(theJetIndices.begin(), theJetIndices.end(), sortByPFJetPt);

                //
                // jet kinematic quantities 
                // for TTBarV2 selection cleaning
                // 

                njetsSS30_    = theJetIndices.size();
                jet1pt_       = theJetIndices.size() > 0 ? cms2.pfjets_p4()[theJetIndices[0]].pt()*cms2.pfjets_corL1FastL2L3()[theJetIndices[0]] : -999999.;
                jet1eta_      = theJetIndices.size() > 0 ? cms2.pfjets_p4()[theJetIndices[0]].eta() : -999999.;
                jet1phi_      = theJetIndices.size() > 0 ? cms2.pfjets_p4()[theJetIndices[0]].phi() : -999999.;
                jet1cor_      = theJetIndices.size() > 0 ? cms2.pfjets_corL1FastL2L3()[theJetIndices[0]] : -999999.;
                jet2pt_       = theJetIndices.size() > 1 ? cms2.pfjets_p4()[theJetIndices[1]].pt()*cms2.pfjets_corL1FastL2L3()[theJetIndices[1]] : -999999.;
                jet2eta_      = theJetIndices.size() > 1 ? cms2.pfjets_p4()[theJetIndices[1]].eta() : -999999.;
                jet2phi_      = theJetIndices.size() > 1 ? cms2.pfjets_p4()[theJetIndices[1]].phi() : -999999.;
                jet2cor_      = theJetIndices.size() > 1 ? cms2.pfjets_corL1FastL2L3()[theJetIndices[1]] : -999999.;
                jet3pt_       = theJetIndices.size() > 2 ? cms2.pfjets_p4()[theJetIndices[2]].pt()*cms2.pfjets_corL1FastL2L3()[theJetIndices[2]] : -999999.;
                jet3eta_      = theJetIndices.size() > 2 ? cms2.pfjets_p4()[theJetIndices[2]].eta() : -999999.;
                jet3phi_      = theJetIndices.size() > 2 ? cms2.pfjets_p4()[theJetIndices[2]].phi() : -999999.;
                jet3cor_      = theJetIndices.size() > 2 ? cms2.pfjets_corL1FastL2L3()[theJetIndices[2]] : -999999.;

                LorentzVector dijetP4;
                jetmass_ = theJetIndices.size() > 1 ? sqrt((cms2.pfjets_p4()[theJetIndices[0]]*cms2.pfjets_corL1FastL2L3()[theJetIndices[0]] 
                                                            + cms2.pfjets_p4()[theJetIndices[1]]*cms2.pfjets_corL1FastL2L3()[theJetIndices[1]]).M2()) : -999999.; 

                mlljj_ = theJetIndices.size() > 1 ? sqrt((cms2.hyp_p4()[hypi] + cms2.pfjets_p4()[theJetIndices[0]]*cms2.pfjets_corL1FastL2L3()[theJetIndices[0]]
                                                          + cms2.pfjets_p4()[theJetIndices[1]]*cms2.pfjets_corL1FastL2L3()[theJetIndices[1]]).M2()) : -999999.;

                mllj_ = theJetIndices.size() > 0 ? sqrt((cms2.hyp_p4()[hypi] 
                                                         + cms2.pfjets_p4()[theJetIndices[0]]*cms2.pfjets_corL1FastL2L3()[theJetIndices[0]]).M2()) : -999999.;

                // rescaled jet1 taking into account met
                if (theJetIndices.size() == 1) {
                    float dphi_jet_met = cms2.pfjets_p4()[theJetIndices[0]].phi() - thePFMetPhi;
                    jet1fmetcor_ = (cms2.pfjets_p4()[theJetIndices[0]].pt() + pfmet_*cos(dphi_jet_met)) / cms2.pfjets_p4()[theJetIndices[0]].pt();
                    mlljmetcor_ = sqrt((cms2.hyp_p4()[hypi] + cms2.pfjets_p4()[theJetIndices[0]]*jet1fmetcor_).M2());

                    // sanity check
                    float metx_prime = pfmet_ * cos(thePFMetPhi) - (cms2.pfjets_p4()[theJetIndices[0]].px()*jet1fmetcor_ - cms2.pfjets_p4()[theJetIndices[0]].px());
                    float mety_prime = pfmet_ * sin(thePFMetPhi) - (cms2.pfjets_p4()[theJetIndices[0]].py()*jet1fmetcor_ - cms2.pfjets_p4()[theJetIndices[0]].py());
                    pfmetcor_ = sqrt(metx_prime*metx_prime + mety_prime*mety_prime);
                    pfmetphicor_ = atan2(mety_prime, metx_prime);

                }

                double mindphipfmet = 999999.;
                double mindphitcmet = 999999.;
                ntchelbtags_  = 0;
                nssvhembtags_ = 0;
                nssvhetbtags_ = 0;
                nssvhptbtags_ = 0;
                jet1isBtag_   = 0;
                jet2isBtag_   = 0;
                jet3isBtag_   = 0;

                //
                // loop on the jets
                // BTag WP are documented in:
                // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
                //

                for (unsigned int jeti = 0; jeti < theJetIndices.size(); ++jeti)
                {

                    // TCHEL
                    if (cms2.pfjets_trackCountingHighEffBJetTag()[jeti] > 1.7) {
                        ++ntchelbtags_;
                        if (jeti == 0) jet1isBtag_ |= 0x1;
                        if (jeti == 1) jet2isBtag_ |= 0x1;
                        if (jeti == 2) jet3isBtag_ |= 0x1;                        
                    }

                    // SSVHEM
                    if (cms2.pfjets_simpleSecondaryVertexHighEffBJetTag()[theJetIndices[jeti]] > 1.74)
                    {
                        ++nssvhembtags_;
                        if (jeti == 0) jet1isBtag_ |= 0x2;
                        if (jeti == 1) jet2isBtag_ |= 0x2;
                        if (jeti == 2) jet3isBtag_ |= 0x2;
                    }

                    // SSVHET
                    if (cms2.pfjets_simpleSecondaryVertexHighEffBJetTag()[theJetIndices[jeti]] > 3.05)
                    {
                        ++nssvhetbtags_;
                        if (jeti == 0) jet1isBtag_ |= 0x4;
                        if (jeti == 1) jet2isBtag_ |= 0x4;
                        if (jeti == 2) jet3isBtag_ |= 0x4;
                    }

                    // SSVHPT
                    if (cms2.pfjets_simpleSecondaryVertexHighPurBJetTags()[theJetIndices[jeti]] > 2.)
                    {
                        ++nssvhptbtags_;
                        if (jeti == 0) jet1isBtag_ |= 0x8;
                        if (jeti == 1) jet2isBtag_ |= 0x8;
                        if (jeti == 2) jet3isBtag_ |= 0x8;

                    }

                    // TCHEM
                    if (cms2.pfjets_trackCountingHighEffBJetTag()[jeti] > 3.3) {
                        if (jeti == 0) jet1isBtag_ |= 0x16;
                        if (jeti == 1) jet2isBtag_ |= 0x16;
                        if (jeti == 2) jet3isBtag_ |= 0x16;                        
                    }

                    float currdphipfmet = deltaPhi(thePFMetPhi, cms2.pfjets_p4()[theJetIndices[jeti]].phi());
                    if (currdphipfmet < mindphipfmet)
                        mindphipfmet = currdphipfmet;

                    float currdphitcmet = deltaPhi(theTCMetPhi, cms2.pfjets_p4()[theJetIndices[jeti]].phi());
                    if (currdphitcmet < mindphitcmet)
                        mindphitcmet = currdphitcmet;

                    // add jet pt to meff
                    pfmeff_ += cms2.pfjets_p4()[theJetIndices[jeti]].pt() * cms2.pfjets_cor()[theJetIndices[jeti]];
                    tcmeff_ += cms2.pfjets_p4()[theJetIndices[jeti]].pt() * cms2.pfjets_cor()[theJetIndices[jeti]];
                }

                // initialize btag counts
                ntchembtags20_ = 0;
                ntchembtags30_ = 0;
                ntchembtags40_ = 0;

                ntchelbtags20_ = 0;
                ntchelbtags30_ = 0;
                ntchelbtags40_ = 0;

                nssvhembtags20_ = 0;
                nssvhembtags30_ = 0;
                nssvhembtags40_ = 0;
                
                for (unsigned int jeti = 0; jeti < theBtagJetIndices.size(); jeti++) {
                    
                    float pt = cms2.pfjets_p4()[theBtagJetIndices[jeti]].pt() * cms2.pfjets_corL1FastL2L3()[theBtagJetIndices[jeti]];

                    // TCHEM
                    if (cms2.pfjets_trackCountingHighEffBJetTag()[theBtagJetIndices[jeti]] > 3.3) {
                        if (pt > 20.)
                            ++ntchembtags20_;
                        if (pt > 30.)
                            ++ntchembtags30_;
                        if (pt > 40.)
                            ++ntchembtags40_;
                    }
                    
                    // TCHEL
                    if (cms2.pfjets_trackCountingHighEffBJetTag()[theBtagJetIndices[jeti]] > 1.7) {
                        if (pt > 20.)
                            ++ntchelbtags20_;
                        if (pt > 30.)
                            ++ntchelbtags30_;
                        if (pt > 40.)
                            ++ntchelbtags40_;
                    }

                    // SSVHEM
                    if (cms2.pfjets_simpleSecondaryVertexHighEffBJetTag()[theBtagJetIndices[jeti]] > 1.74) {
                        if (pt > 20.)
                            ++nssvhembtags20_;
                        if (pt > 30.)
                            ++nssvhembtags30_;
                        if (pt > 40.)
                            ++nssvhembtags40_;
                    }
                }

                dphipfmetjet_ = mindphipfmet;
                dphitcmetjet_ = mindphitcmet;

                // now, find jet closest to each hyp lepton
                float mindrjet1 = 999999.;
                float mindrjet2 = 999999.;
                for (unsigned int jeti = 0; jeti < theJetIndices.size(); ++jeti)
                {
                    // for tight hyp lepton
                    float deta1 = cms2.hyp_lt_p4()[hypi].eta()-cms2.pfjets_p4()[theJetIndices[jeti]].eta();
                    float dphi1 = deltaPhi(cms2.hyp_lt_p4()[hypi].phi(), cms2.pfjets_p4()[theJetIndices[jeti]].phi());
                    float currdrjet1 = sqrt(deta1*deta1+dphi1*dphi1);
                    if (currdrjet1 < mindrjet1)
                        mindrjet1 = currdrjet1;

                    // for loose hyp lepton
                    float deta2 = cms2.hyp_ll_p4()[hypi].eta()-cms2.pfjets_p4()[theJetIndices[jeti]].eta();
                    float dphi2 = deltaPhi(cms2.hyp_ll_p4()[hypi].phi(), cms2.pfjets_p4()[theJetIndices[jeti]].phi());
                    float currdrjet2 = sqrt(deta2*deta2+dphi2*dphi2);
                    if (currdrjet2 < mindrjet2)
                        mindrjet2 = currdrjet2;
                }

                drjet1_ = mindrjet1;
                drjet2_ = mindrjet2;

                // mt2
                mt2_ = MT2(pfmet_, thePFMetPhi, cms2.hyp_ll_p4()[hypi], cms2.hyp_lt_p4()[hypi]);
                // mt2j
                if (theJets.size() > 1)
                    mt2j_ = MT2J(pfmet_, thePFMetPhi, cms2.hyp_ll_p4()[hypi], cms2.hyp_lt_p4()[hypi], theJets);

                // MT for Higgs

                float enell = TMath::Sqrt(dilpt_*dilpt_ + mass_*mass_);
                float pfenenn = TMath::Sqrt(pfmet_*pfmet_  + mass_*mass_);
                float pfenex  = cms2.hyp_p4()[hypi].Px() + pfmetx;
                float pfeney  = cms2.hyp_p4()[hypi].Py() + pfmety;
                float tcenenn = TMath::Sqrt(tcmet_*tcmet_  + mass_*mass_);
                float tcenex  = cms2.hyp_p4()[hypi].Px() + metx;
                float tceney  = cms2.hyp_p4()[hypi].Py() + mety;
                pfmth_ = sqrt(2.0*mass_*mass_ + 2.0*(enell*pfenenn - pfenex*pfenex - pfeney*pfeney));
                tcmth_ = sqrt(2.0*mass_*mass_ + 2.0*(enell*tcenenn - tcenex*tcenex - tceney*tceney));

                // check if hypothesis makes an extra Z
//                extraZveto_ = makesExtraZ(hypi, false, false);
                extraZveto_ = samesign::overlapsOtherNNHypInZ(hypi);

                //
                // Now fill the detailed electron/muon specific variables
                // for each leg of the hypothesis
                //

                //
                // First do LT
                //

                // if LT is a mu, fill mu info
                if (abs(cms2.hyp_lt_id()[hypi]) == 13) {
                    iso1_   = muonIsoValue(index1);
                    ntiso1_ = muonIsoValue(index1, false); 
                    type1_  = cms2.mus_type()[index1];
                    mu1_numSSv4_      = muonId(index1, NominalSSv4);
                    mu1_numSSv4noIso_ = muonIdNotIsolated(index1, NominalSSv4);
                    mu1_foSSv4_       = muonId(index1, muonSelectionFO_ssV4);
                    mu1_foSSv4noIso_  = muonIdNotIsolated(index1, muonSelectionFO_ssV4);
                    mu1_numSSv3_      = muonId(index1, NominalSSv3);
                    mu1_foSSv3_       = muonId(index1, muonSelectionFO_ssV3);
                    mu1_muonidfull_   = muonId(index1, OSGeneric_v3);
                    mu1_muonid_       = muonIdNotIsolated(index1, OSGeneric_v3);
                    mu1_muonidfullV1_ = muonId(index1, OSGeneric_v3);
                    mu1_muonidV1_     = muonIdNotIsolated(index1, OSGeneric_v3);
                    mu1_goodmask_     = cms2.mus_goodmask()[index1];
                    mu1_gfitchi2_     = cms2.mus_gfit_chi2()[index1] < -9000. ? -999999. : cms2.mus_gfit_chi2()[index1]/cms2.mus_gfit_ndof()[index1];
                    mu1_cosmic_       = isCosmics(index1);
                    mu1_siHits_       = cms2.mus_validHits()[index1];
                    mu1_saHits_       = cms2.mus_gfit_validSTAHits()[index1];
                    mu1_emVetoDep_    = cms2.mus_iso_ecalvetoDep()[index1];
                    mu1_hadVetoDep_   = cms2.mus_iso_hcalvetoDep()[index1];
                    mu1_relPtErr_     = cms2.mus_ptErr()[index1]/cms2.mus_p4()[index1].pt();
                    mu1_trkpt_        = cms2.mus_trk_p4()[index1].pt();
                    if (cms2.mus_pfmusidx()[index1] > -1) mu1_isPFmuon_ = 1;
                    int trkidx1 = cms2.mus_trkidx()[index1];
                    int vtxidx = firstGoodVertex();
                    if (trkidx1 >= 0 && vtxidx >=0)
                        d0vtx1_ = trks_d0_pv(trkidx1, vtxidx, true).first;
                    trkIso1_ = cms2.mus_iso03_sumPt()[index1];
                    ecalIso1_ = cms2.mus_iso03_emEt()[index1];
                    hcalIso1_ = cms2.mus_iso03_hadEt()[index1];
                    ecalIso1ps_ = cms2.mus_iso03_emEt()[index1];
                    lep1_pfiso0p3_ = muonIsoValuePF(index1, 0, 0.3);
                    lep1_pfiso0p4_ = muonIsoValuePF(index1, 0, 0.4);
                }
                // if LT is an ele, fill ele info
                if (abs(cms2.hyp_lt_id()[hypi]) == 11) {
                    iso1_   = electronIsolation_rel(index1, true);
                    ntiso1_ = electronIsolation_rel_v1(index1, true);
                    type1_  = cms2.els_type()[index1];
                    electronIdComponent_t answer_vbtf90 = electronId_VBTF(index1, VBTF_35X_90);
                    e1_vbtf90_      = (answer_vbtf90 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf85 = electronId_VBTF(index1, VBTF_35X_85);
                    e1_vbtf85_      = (answer_vbtf85 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf80 = electronId_VBTF(index1, VBTF_35X_80);
                    e1_vbtf80_      = (answer_vbtf80 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf70 = electronId_VBTF(index1, VBTF_35X_70);
                    e1_vbtf70_      = (answer_vbtf70 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);

                    cuts_t electron_selection = electronSelection(index1);
                    e1_vbtf90full_   = pass_electronSelectionCompareMask(electron_selection, electronSelection_el_OSV3);
                    e1_smurfV3_      = pass_electronSelectionCompareMask(electron_selection, electronSelection_smurfV3_id);
                    e1_numSSv3_      = pass_electronSelectionCompareMask(electron_selection, electronSelection_ssV3);
                    e1_foSSv3_       = pass_electronSelectionCompareMask(electron_selection, electronSelectionFOV3_ssVBTF80_v3);
                    e1_numSSv4_      = pass_electronSelectionCompareMask(electron_selection, electronSelection_ssV4);
                    e1_foSSv4_       = pass_electronSelectionCompareMask(electron_selection, electronSelectionFOV4_ssVBTF80_v3);
                    e1_numSSv5_      = pass_electronSelectionCompareMask(electron_selection, electronSelection_ssV5);
                    e1_foSSv5_       = pass_electronSelectionCompareMask(electron_selection, electronSelectionFOV5_ssVBTF80_v3);
                    e1_numSSv5noIso_ = pass_electronSelectionCompareMask(electron_selection, electronSelection_ssV5_noIso);
                    e1_numSSv5_noConvCuts_      = pass_electronSelectionCompareMask(electron_selection, electronSelection_ssV5_noConvCuts);
                    e1_foSSv5_noConvCuts_       = pass_electronSelectionCompareMask(electron_selection, electronSelectionFOV5_ssVBTF80_noConvCuts_v3);
                    e1_numSSv5noIso_noConvCuts_ = pass_electronSelectionCompareMask(electron_selection, electronSelection_ssV5_noIso_noConvCuts);

                    e1_scet_        = cms2.els_eSC()[index1] / cosh(cms2.els_etaSC()[index1]);
                    e1_sceta_       = cms2.els_etaSC()[index1];
                    e1_eopin_       = cms2.els_eOverPIn()[index1];
                    e1_hoe_         = cms2.els_hOverE()[index1];
                    e1_dphiin_      = cms2.els_dPhiIn()[index1];
                    e1_detain_      = cms2.els_dEtaIn()[index1];
                    e1_e25Me55_     = cms2.els_e2x5Max()[index1] / cms2.els_e5x5()[index1];
//                    e1_e15Me55_     = cms2.els_e1x5Max()[index1] / cms2.els_e5x5()[index1];
                    e1_sigieie_     = cms2.els_sigmaIEtaIEta()[index1];
                    e1_eMe55_       = cms2.els_eMax()[index1] / cms2.els_e5x5()[index1];
                    e1_nmHits_      = cms2.els_exp_innerlayers()[index1];
                    e1_dcot_        = cms2.els_conv_dcot()[index1];
                    e1_dist_        = cms2.els_conv_dist()[index1];
                    e1_dcot_old_    = cms2.els_conv_old_dcot()[index1];
                    e1_dist_old_    = cms2.els_conv_old_dist()[index1];
                    e1_drmu_        = cms2.els_closestMuon()[index1] < 0 ? -999999. : cms2.els_musdr()[index1];
                    if (cms2.els_closestMuon()[index1] < 0)
                        e1_drmuSS_ = -999999.;
                    else if (!muonId(cms2.els_closestMuon()[index1], NominalSSv4))
                        e1_drmuSS_ = -999999.;
                    else
                        e1_drmuSS_ = cms2.els_musdr()[index1];
                    e1_isspike_     = isSpikeElectron(index1);
                    e1_scCharge_   = cms2.els_sccharge()[index1];
                    e1_gsfCharge_  = cms2.els_trk_charge()[index1];
                    e1_ctfCharge_  = cms2.els_trkidx()[index1] > -1 ? cms2.trks_charge()[cms2.els_trkidx()[index1]] : -999999;
                    int trkidx1 = cms2.els_gsftrkidx()[index1];
                    int vtxidx = firstGoodVertex();
                    if (trkidx1 >= 0 && vtxidx >= 0)
                        d0vtx1_ = gsftrks_d0_pv(trkidx1, vtxidx, true).first;
                    e1_fbrem_ = cms2.els_fbrem()[index1];
                    trkIso1_ = cms2.els_tkIso()[index1];
                    hcalIso1_ = cms2.els_hcalIso()[index1];
                    ecalIso1_ = cms2.els_ecalIso()[index1];
                    e1_mitConv_ = !isFromConversionMIT(index1);
                    ecalIso1ps_ = fabs(cms2.els_p4()[index1].eta()) < 1.479 ? max(cms2.els_ecalIso()[index1] - 1., 0.) : cms2.els_ecalIso()[index1];
                    lep1_pfiso0p3_ = electronIsoValuePF(index1, 0, 0.3);
                    lep1_pfiso0p4_ = electronIsoValuePF(index1, 0, 0.4);

                    e1_eSeedOverPout_         = cms2.els_eSeedOverPOut().at(index1);
                    e1_eSeedOverPin_          = cms2.els_eSeedOverPIn().at(index1);
                    e1_sigipip_               = cms2.els_sigmaIPhiIPhi().at(index1);
                    e1_nbrem_                 = cms2.els_nSeed().at(index1);
                    e1_ip3d_                  = cms2.els_ubIp3d().at(index1);
                    if (cms2.els_ubIp3derr().at(index1) == 0.)                        
                        e1_ip3dSig_ = 0.;
                    else
                        e1_ip3dSig_ = e1_ip3d_ / cms2.els_ubIp3derr().at(index1);
                    e1_lh_                    = cms2.els_lh().at(index1);
                    TVector3 pIn(cms2.els_trk_p4().at(index1).px(), cms2.els_trk_p4().at(index1).py(), cms2.els_trk_p4().at(index1).pz());
                    e1_oneOverEminusOneOverP_ = (1. / pIn.Mag()) * ( (1/e1_eopin_) - 1.);
                    e1_mva_                   = electronIdMVA->MVAValue(index1, 0);
                }

                //
                // Now do LL...
                //

                // if LL is a mu, fill mu info
                if (abs(cms2.hyp_ll_id()[hypi]) == 13) {
                    iso2_   = muonIsoValue(index2);
                    ntiso2_ = muonIsoValue(index2, false);
                    type2_  = cms2.mus_type()[index2];
                    mu2_numSSv4_      = muonId(index2, NominalSSv4);
                    mu2_numSSv4noIso_ = muonIdNotIsolated(index2, NominalSSv4);
                    mu2_foSSv4_       = muonId(index2, muonSelectionFO_ssV4);
                    mu2_foSSv4noIso_  = muonIdNotIsolated(index2, muonSelectionFO_ssV4);
                    mu2_numSSv3_      = muonId(index2, NominalSSv3);
                    mu2_foSSv3_       = muonId(index2, muonSelectionFO_ssV3);
                    mu2_muonidfull_   = muonId(index2, OSGeneric_v3);
                    mu2_muonid_       = muonIdNotIsolated(index2, OSGeneric_v3);
                    mu2_muonidfullV1_ = muonId(index2, OSGeneric_v3);
                    mu2_muonidV1_     = muonIdNotIsolated(index2, OSGeneric_v3);
                    mu2_goodmask_     = cms2.mus_goodmask()[index2];
                    mu2_gfitchi2_     = cms2.mus_gfit_chi2()[index2] < -9000. ? -999999. : cms2.mus_gfit_chi2()[index2]/cms2.mus_gfit_ndof()[index2];
                    mu2_cosmic_       = isCosmics(index2);
                    mu2_siHits_       = cms2.mus_validHits()[index2];
                    mu2_saHits_       = cms2.mus_gfit_validSTAHits()[index2];
                    mu2_emVetoDep_    = cms2.mus_iso_ecalvetoDep()[index2];
                    mu2_hadVetoDep_   = cms2.mus_iso_hcalvetoDep()[index2];
                    mu2_relPtErr_     = cms2.mus_ptErr()[index2]/cms2.mus_p4()[index2].pt();
                    mu2_trkpt_        = cms2.mus_trk_p4()[index2].pt();
                    if (cms2.mus_pfmusidx()[index2] > -1) mu2_isPFmuon_ = 1;
                    int trkidx2 = cms2.mus_trkidx()[index2];
                    int vtxidx = firstGoodVertex();
                    if (trkidx2 >= 0 && vtxidx >= 0)
                        d0vtx2_ = trks_d0_pv(trkidx2, vtxidx, true).first;
                    trkIso2_ = cms2.mus_iso03_sumPt()[index2];
                    ecalIso2_ = cms2.mus_iso03_emEt()[index2];
                    hcalIso2_ = cms2.mus_iso03_hadEt()[index2];
                    ecalIso2ps_ = cms2.mus_iso03_emEt()[index2];
                    lep2_pfiso0p3_ = muonIsoValuePF(index2, 0, 0.3);
                    lep2_pfiso0p4_ = muonIsoValuePF(index2, 0, 0.4);
                }
                // if LL is an ele, fill ele info
                if (abs(cms2.hyp_ll_id()[hypi]) == 11) {
                    iso2_   = electronIsolation_rel(index2, true);
                    ntiso2_ = electronIsolation_rel_v1(index2, true);
                    type2_  = cms2.els_type()[index2];
                    electronIdComponent_t answer_vbtf90   = electronId_VBTF(index2, VBTF_35X_90);
                    e2_vbtf90_      = (answer_vbtf90 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf85   = electronId_VBTF(index2, VBTF_35X_85);
                    e2_vbtf85_      = (answer_vbtf85 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf80   = electronId_VBTF(index2, VBTF_35X_80);
                    e2_vbtf80_      = (answer_vbtf80 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf70   = electronId_VBTF(index2, VBTF_35X_70);
                    e2_vbtf70_      = (answer_vbtf70 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);

                    cuts_t electron_selection = electronSelection(index2);
                    e2_vbtf90full_ = pass_electronSelectionCompareMask(electron_selection, electronSelection_el_OSV3);
                    e2_smurfV3_ = pass_electronSelectionCompareMask(electron_selection, electronSelection_smurfV3_id);
                    e2_numSSv3_ = pass_electronSelectionCompareMask(electron_selection, electronSelection_ssV3);
                    e2_foSSv3_  = pass_electronSelectionCompareMask(electron_selection, electronSelectionFOV3_ssVBTF80_v3);
                    e2_numSSv4_ = pass_electronSelectionCompareMask(electron_selection, electronSelection_ssV4);
                    e2_foSSv4_  = pass_electronSelectionCompareMask(electron_selection, electronSelectionFOV4_ssVBTF80_v3);
                    e2_numSSv5_ = pass_electronSelectionCompareMask(electron_selection, electronSelection_ssV5);
                    e2_foSSv5_  = pass_electronSelectionCompareMask(electron_selection, electronSelectionFOV5_ssVBTF80_v3);
                    e2_numSSv5noIso_ = pass_electronSelectionCompareMask(electron_selection, electronSelection_ssV5_noIso);
                    e2_numSSv5_noConvCuts_      = pass_electronSelectionCompareMask(electron_selection, electronSelection_ssV5_noConvCuts);
                    e2_foSSv5_noConvCuts_       = pass_electronSelectionCompareMask(electron_selection, electronSelectionFOV5_ssVBTF80_noConvCuts_v3);
                    e2_numSSv5noIso_noConvCuts_ = pass_electronSelectionCompareMask(electron_selection, electronSelection_ssV5_noIso_noConvCuts);

                    e2_scet_        = cms2.els_eSC()[index2] / cosh(cms2.els_etaSC()[index2]);
                    e2_sceta_       = cms2.els_etaSC()[index2];
                    e2_eopin_       = cms2.els_eOverPIn()[index2];
                    e2_hoe_         = cms2.els_hOverE()[index2];
                    e2_dphiin_      = cms2.els_dPhiIn()[index2];
                    e2_detain_      = cms2.els_dEtaIn()[index2]; 
                    e2_e25Me55_     = cms2.els_e2x5Max()[index2] / cms2.els_e5x5()[index2];
//                    e2_e15Me55_     = cms2.els_e1x5Max()[index2] / cms2.els_e5x5()[index2];
                    e2_sigieie_     = cms2.els_sigmaIEtaIEta()[index2];
                    e2_eMe55_       = cms2.els_eMax()[index2] / cms2.els_e5x5()[index2];
                    e2_nmHits_      = cms2.els_exp_innerlayers()[index2];
                    e2_dcot_        = cms2.els_conv_dcot()[index2];
                    e2_dist_        = cms2.els_conv_dist()[index2];
                    e2_dcot_old_    = cms2.els_conv_old_dcot()[index2];
                    e2_dist_old_    = cms2.els_conv_old_dist()[index2];
                    e2_drmu_        = cms2.els_closestMuon()[index2] < 0 ? -999999. : cms2.els_musdr()[index2];
                    if (cms2.els_closestMuon()[index2] < 0)
                        e2_drmuSS_ = -999999.;
                    else if (!muonId(cms2.els_closestMuon()[index2], NominalSSv4))
                        e2_drmuSS_ = -999999.;
                    else
                        e2_drmuSS_ = cms2.els_musdr()[index2];
                    e2_isspike_     = isSpikeElectron(index2);
                    e2_scCharge_    = cms2.els_sccharge()[index2];
                    e2_gsfCharge_   = cms2.els_trk_charge()[index2]; 
                    e2_ctfCharge_   = cms2.els_trkidx()[index2] > -1 ? cms2.trks_charge()[cms2.els_trkidx()[index2]] : -999999;
                    int trkidx2 = cms2.els_gsftrkidx()[index2];
                    int vtxidx = firstGoodVertex();
                    if (trkidx2 >= 0 && vtxidx >= 0)
                        d0vtx2_ = gsftrks_d0_pv(trkidx2, vtxidx, true).first;
                    e2_fbrem_ = cms2.els_fbrem()[index2];
                    trkIso2_ = cms2.els_tkIso()[index2];
                    hcalIso2_ = cms2.els_hcalIso()[index2];
                    ecalIso2_ = cms2.els_ecalIso()[index2];
                    e2_mitConv_ = !isFromConversionMIT(index2);
                    ecalIso2ps_ = fabs(cms2.els_p4()[index2].eta()) < 1.479 ? max(cms2.els_ecalIso()[index2] - 1., 0.) : cms2.els_ecalIso()[index2];
                    lep2_pfiso0p3_ = electronIsoValuePF(index2, 0, 0.3);
                    lep2_pfiso0p4_ = electronIsoValuePF(index2, 0, 0.4);

                    e2_eSeedOverPout_         = cms2.els_eSeedOverPOut().at(index2);
                    e2_eSeedOverPin_          = cms2.els_eSeedOverPIn().at(index2);
                    e2_sigipip_               = cms2.els_sigmaIPhiIPhi().at(index2);
                    e2_nbrem_                 = cms2.els_nSeed().at(index2);
                    e2_ip3d_                  = cms2.els_ubIp3d().at(index2);
                    if (cms2.els_ubIp3derr().at(index2) == 0.)                        
                        e2_ip3dSig_ = 0.;
                    else
                        e2_ip3dSig_ = e2_ip3d_ / cms2.els_ubIp3derr().at(index2);
                    e2_lh_                    = cms2.els_lh().at(index2);
                    TVector3 pIn(cms2.els_trk_p4().at(index2).px(), cms2.els_trk_p4().at(index2).py(), cms2.els_trk_p4().at(index2).pz());
                    e2_oneOverEminusOneOverP_ = (1. / pIn.Mag()) * ( (1/e2_eopin_) - 1.);
                    e2_mva_                   = electronIdMVA->MVAValue(index2, 0);
                }

                //
                // Now write everything to the baby
                //

                FillBabyNtuple();
            }
        }

        f->Close();

    }

    if (nEventsChain != nEventsTotal)
        std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;


    CloseBabyNtuple();
}

void dilepbabymaker::InitBabyNtuple ()
{
    // event stuff
    memset(dataset_, '\0', 200);

    rndm_             = -999999.;
    run_              = -999999;
    ls_               = -999999;
    evt_              = 0;
    isdata_           = 1;
    evt_clean082010_  = -999999;
    evt_clean102010_  = -999999;    
    evt_clean042011_  = -999999;    
    nvtx_             = -999999;
    scale1fb_         = -999999.;
    pthat_            = -999999.;
    hyp_type_         = -999999;
    ngenjets_         = -999999;
    genmet_           = -999999.;
    gensumjetpt_      = -999999.;
    pfmet_            = -999999.;
    tcmet_            = -999999.;
    pfmetphi_         = -999999.;
    tcmetphi_         = -999999.;
    calotcmet_        = -999999.;
    proj_pfmet_       = -999999.;
    proj_tcmet_       = -999999.;
    ntrks_            = -999999;
    njets_            = -999999;
    njets25_          = -999999;
    njetsSS30_        = -999999;
    njetsSS_          = -999999;
    njetsSShighPt_    = -999999;
    jet1pt_           = -999999.;
    jet2pt_           = -999999.;
    jet3pt_           = -999999.;
    sumjetpt_         = -999999.;
    sumjetpt25_       = -999999.;
    sumjetptSS30_     = -999999.;
    sumjetptSS_       = -999999.;
    sumjetptSShighPt_ = -999999.;
    jet1cor_          = -9999999.;
    jet2cor_          = -9999999.;
    jet3cor_          = -9999999.;
    jet1eta_          = -999999.;
    jet2eta_          = -999999.;
    jet3eta_          = -999999.;
    jet1phi_          = -999999.;
    jet2phi_          = -999999.;
    jet3phi_          = -999999.;
    jetmass_          = -999999.;
    mllj_             = -999999.;
    mlljj_            = -999999.;
    pfmth_            = -999999.;
    tcmth_            = -999999.;
    jet1isBtag_       = 0;
    jet2isBtag_       = 0;
    jet3isBtag_       = 0;
    dphipfmetjet_     = -999999.;
    dphitcmetjet_     = -999999.;
    deltaphi_         = -999999.;

    jet1fmetcor_      = -9999999.;
    mlljmetcor_       = -9999999.;
    pfmetcor_         = -9999999.;
    pfmetphicor_      = -9999999.;

    ntchelbtags_      = -999999;
    nssvhembtags_     = -999999;
    nssvhetbtags_     = -999999;
    nssvhptbtags_     = -999999;
    ntchelbtags20_    = -999999;
    ntchelbtags30_    = -999999;
    ntchelbtags40_    = -999999;
    ntchembtags20_    = -999999;
    ntchembtags30_    = -999999;
    ntchembtags40_    = -999999;
    nssvhembtags20_   = -999999;
    nssvhembtags30_   = -999999;
    nssvhembtags40_   = -999999;

    pfmeff_           = -999999.;
    tcmeff_           = -999999.;

    intLumiPerLS_     = -999999.;

    // lepton stuff
    ngoodlep_         = -999999;
    ngoodmus_         = -999999;
    ngoodels_         = -999999;
    ngoodlepSS_       = -999999;
    ngoodmusSS_       = -999999;
    ngoodelsSS_       = -999999;
    ngenels_          = -999999;
    ngenmus_          = -999999;
    ngentaus_         = -999999;

    dilpt_            = -999999.;
    dileta_           = -999999.;
    dilphi_           = -999999.;
    mass_             = -999999.;
    eormu1_           = -999999;
    type1_            = -999999;
    pt1_              = -999999.;
    eta1_             = -999999.;
    phi1_             = -999999.;
    iso1_             = -999999.;
    ntiso1_           = -999999.;
    d0corr1_          = -999999.;
    d0vtx1_           = -999999.;
    dphipfmet1_       = -999999.;
    dphitcmet1_       = -999999.;
    drjet1_           = -999999.;
    mcid1_            = -999999;
    mcmotherid1_      = -999999;
    mcpt1_            = -999999;
    mcmotherpt1_      = -999999.;
    mc3id1_           = -999999;
    mc3motherid1_     = -999999;
    mc3pt1_           = -999999;
    mc3motherpt1_     = -999999.;
    eormu2_           = -999999;
    type2_            = -999999;
    pt2_              = -999999.;
    eta2_             = -999999.;
    phi2_             = -999999.;
    iso2_             = -999999.;
    ntiso2_           = -999999.;
    d0corr2_          = -999999.;
    d0vtx2_           = -999999.;
    dphipfmet2_       = -999999.;
    dphitcmet2_       = -999999.;
    drjet2_           = -999999.;
    mcid2_            = -999999;
    mcmotherid2_      = -999999;
    mcpt2_            = -999999.;
    mcmotherpt2_      = -999999.;
    mc3id2_           = -999999;
    mc3motherid2_     = -999999;
    mc3pt2_           = -999999.;
    mc3motherpt2_     = -999999.;
    mt2_              = -999999.;
    mt2j_             = -999999.;
    extraZveto_       = 0;
    trkIso1_          = -999999.;
    ecalIso1_         = -999999.;
    hcalIso1_         = -999999.;
    trkIso2_          = -999999.;
    ecalIso2_         = -999999.;
    hcalIso2_         = -999999.;
    ecalIso1ps_       = -999999.;
    ecalIso2ps_       = -999999.;
    rho_              = -999999.;
    lepsFromSameVtx_  = 0;
    lep1isFromW_      = -999999;
    lep2isFromW_      = -999999;
    lep1_pfiso0p3_    = -999999.;
    lep1_pfiso0p4_    = -999999.;
    lep2_pfiso0p3_    = -999999.;
    lep2_pfiso0p4_    = -999999.;

    // muon stuff
    mu1_numSSv4_      = 0;
    mu1_numSSv4noIso_      = 0;
    mu1_foSSv4_       = 0;
    mu1_foSSv4noIso_       = 0;
    mu1_numSSv3_      = 0;
    mu1_foSSv3_       = 0;
    mu1_muonidfull_   = 0;
    mu1_muonid_       = 0;
    mu1_muonidfullV1_ = 0;
    mu1_muonidV1_     = 0;
    mu1_goodmask_     = -999999;
    mu1_gfitchi2_     = -999999.;
    mu1_cosmic_       = 0;
    mu1_siHits_       = -999999;
    mu1_saHits_       = -999999;
    mu1_emVetoDep_    = -999999.;
    mu1_hadVetoDep_   = -999999.;
    mu1_isPFmuon_     = 0;
    mu1_relPtErr_     = -999999.;
    mu1_trkpt_        = -999999.;

    mu2_numSSv4_      = 0;
    mu2_numSSv4noIso_      = 0;
    mu2_foSSv4_       = 0;
    mu2_foSSv4noIso_       = 0;
    mu2_numSSv3_      = 0;
    mu2_foSSv3_       = 0;
    mu2_muonidfull_   = 0;
    mu2_muonid_       = 0;
    mu2_muonidfullV1_ = 0;
    mu2_muonidV1_     = 0;
    mu2_goodmask_     = -999999;
    mu2_gfitchi2_     = -999999.;
    mu2_cosmic_       = 0;
    mu2_siHits_       = -999999;
    mu2_saHits_       = -999999;
    mu2_emVetoDep_    = -999999.;
    mu2_hadVetoDep_   = -999999.;
    mu2_isPFmuon_     = 0;
    mu2_relPtErr_     = -999999.;
    mu2_trkpt_        = -999999.;

    // electron stuff
    e1_numSSv5_       = 0;
    e1_numSSv5noIso_  = 0;
    e1_foSSv5_        = 0;
    e1_numSSv5_noConvCuts_       = 0;
    e1_numSSv5noIso_noConvCuts_  = 0;
    e1_foSSv5_noConvCuts_        = 0;
    e1_numSSv4_       = 0;
    e1_foSSv4_        = 0;
    e1_numSSv3_       = 0;
    e1_foSSv3_        = 0;
    e1_vbtf90full_    = 0;
    e1_vbtf90_        = 0;
    e1_vbtf85_        = 0;
    e1_vbtf80_        = 0;
    e1_vbtf70_        = 0;
    e1_smurfV3_       = 0;
    e1_scet_          = -999999.;
    e1_sceta_         = -999999.;
    e1_eopin_         = -999999.;
    e1_hoe_           = -999999.;
    e1_dphiin_        = -999999.;
    e1_detain_        = -999999.;
    e1_e25Me55_       = -999999.;
    e1_e15Me55_       = -999999.;
    e1_sigieie_       = -999999.;
    e1_eMe55_         = -999999.;
    e1_nmHits_        = -999999;
    e1_dcot_          = -999999.;
    e1_dist_          = -999999.;
    e1_dcot_old_      = -999999.;
    e1_dist_old_      = -999999.;
    e1_drmu_          = -999999.;
    e1_drmuSS_        = -999999.;
    e1_isspike_       = 0;
    e1_ctfCharge_     = -999999;
    e1_gsfCharge_     = -999999;
    e1_scCharge_      = -999999;
    e1_fbrem_         = -999999.;
    e1_mitConv_       = 0;
    e1_eSeedOverPout_         = -999999.;
    e1_eSeedOverPin_          = -999999.;
    e1_sigipip_               = -999999.;
    e1_nbrem_                 = -999999.;
    e1_ip3d_                  = -999999.;
    e1_ip3dSig_               = -999999.;
    e1_lh_                    = -999999.;
    e1_oneOverEminusOneOverP_ = -999999.;
    e1_mva_                   = -999999.;

    e2_numSSv5_       = 0;
    e2_numSSv5noIso_  = 0;
    e2_foSSv5_        = 0;
    e2_numSSv5_noConvCuts_       = 0;
    e2_numSSv5noIso_noConvCuts_  = 0;
    e2_foSSv5_noConvCuts_        = 0;
    e2_numSSv4_       = 0;
    e2_foSSv4_        = 0;
    e2_numSSv3_       = 0;
    e2_foSSv3_        = 0;
    e2_vbtf90full_    = 0;
    e2_vbtf90_        = 0;
    e2_vbtf85_        = 0;
    e2_vbtf80_        = 0;
    e2_vbtf70_        = 0;
    e2_smurfV3_       = 0;
    e2_scet_          = -999999.;
    e2_sceta_         = -999999.;
    e2_eopin_         = -999999.;
    e2_hoe_           = -999999.;
    e2_dphiin_        = -999999.;
    e2_detain_        = -999999.;
    e2_e25Me55_       = -999999.;
    e2_e15Me55_       = -999999.;
    e2_sigieie_       = -999999.;
    e2_eMe55_         = -999999.;
    e2_nmHits_        = -999999;
    e2_dcot_          = -999999.;
    e2_dist_          = -999999.;
    e2_dcot_old_      = -999999.;
    e2_dist_old_      = -999999.;
    e2_drmu_          = -999999.;
    e2_drmuSS_        = -999999.;
    e2_isspike_       = 0;
    e2_ctfCharge_     = -999999;
    e2_gsfCharge_     = -999999;
    e2_scCharge_      = -999999;
    e2_fbrem_                 = -999999.;
    e2_mitConv_               = 0;
    e2_eSeedOverPout_         = -999999.;
    e2_eSeedOverPin_          = -999999.;
    e2_sigipip_               = -999999.;
    e2_nbrem_                 = -999999.;
    e2_ip3d_                  = -999999.;
    e2_ip3dSig_               = -999999.;
    e2_lh_                    = -999999.;
    e2_oneOverEminusOneOverP_ = -999999.;
    e2_mva_                   = -999999.;

    pass_trg_single_e1_       = 0;
    pass_trg_single_e2_       = 0;
    pass_trg_single_mu1_      = 0;
    pass_trg_single_mu2_      = 0;
    pass_trg_double_e1_       = 0;
    pass_trg_double_e2_       = 0;
    pass_trg_double_mu1_      = 0;
    pass_trg_double_mu2_      = 0;
    pass_trg_cross_emu_       = 0;

    pass_trg_had_double_e1_   = 0;
    pass_trg_had_double_e2_   = 0;
    pass_trg_had_double_mu1_  = 0;
    pass_trg_had_double_mu2_  = 0;
    pass_trg_had_cross_emu_   = 0;

    match_trg_single_e1_      = 0;
    match_trg_single_e2_      = 0;
    match_trg_single_mu1_     = 0;
    match_trg_single_mu2_     = 0;
    match_trg_double_e1_      = 0;
    match_trg_double_e2_      = 0;
    match_trg_double_mu1_     = 0;
    match_trg_double_mu2_     = 0;
    match_trg_cross_emu_      = 0;

    match_trg_had_double_e1_  = 0;
    match_trg_had_double_e2_  = 0;
    match_trg_had_double_mu1_ = 0;
    match_trg_had_double_mu2_ = 0;
    match_trg_had_cross_emu_  = 0;
}

void dilepbabymaker::MakeBabyNtuple(const char *babyFilename)
{
    // TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    // rootdir->cd();

    babyFile_ = TFile::Open(Form("%s", babyFilename), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");
    babyTree_->SetDirectory(0);

    // event stuff
    babyTree_->Branch("rndm",          &rndm_,         "rndm/F"         );
    babyTree_->Branch("dataset",      &dataset_,     "dataset[200]/C");
    babyTree_->Branch("run",          &run_,         "run/I"         );
    babyTree_->Branch("ls",           &ls_,          "ls/I"          );
    babyTree_->Branch("evt",          &evt_,         "evt/i"         );
    babyTree_->Branch("nvtx",         &nvtx_,        "nvtx/I"        );
    babyTree_->Branch("isdata",       &isdata_,      "isdata/I"      );
    babyTree_->Branch("evt_clean082010",       &evt_clean082010_,      "evt_clean082010/I"      );
    babyTree_->Branch("evt_clean102010",       &evt_clean102010_,      "evt_clean102010/I"      );
    babyTree_->Branch("evt_clean042011",       &evt_clean042011_,      "evt_clean042011/I"      );
    babyTree_->Branch("scale1fb",     &scale1fb_,    "scale1fb/F"    );
    babyTree_->Branch("pthat",        &pthat_,       "pthat/F"       );
    babyTree_->Branch("hyp_type",     &hyp_type_,    "hyp_type/I"    );
    babyTree_->Branch("ngenjets", &ngenjets_, "gen_jets/I");
    babyTree_->Branch("genmet", &genmet_, "genmet/F");
    babyTree_->Branch("gensumjetpt", &gensumjetpt_, "gensumjetpt/F");
    babyTree_->Branch("pfmet",        &pfmet_,       "pfmet/F"       );
    babyTree_->Branch("tcmet",        &tcmet_,       "tcmet/F"       );
    babyTree_->Branch("pfmetphi",        &pfmetphi_,       "pfmetphi/F"       );
    babyTree_->Branch("tcmetphi",        &tcmetphi_,       "tcmetphi/F"       );
    babyTree_->Branch("calotcmet",    &calotcmet_,   "calotcmet/F"   );
    babyTree_->Branch("proj_pfmet",   &proj_pfmet_,  "proj_pfmet/F"  );
    babyTree_->Branch("proj_tcmet",   &proj_tcmet_,  "proj_tcmet/F"  );
    babyTree_->Branch("ntrks",        &ntrks_,       "ntrks/I"       );
    babyTree_->Branch("njets",        &njets_,       "njets/I"       );
    babyTree_->Branch("njets25",        &njets25_,       "njets25/I"       ); 
    babyTree_->Branch("njetsSS30",      &njetsSS30_,     "njetsSS30/I"     );
    babyTree_->Branch("njetsSS",      &njetsSS_,     "njetsSS/I"     );
    babyTree_->Branch("njetsSShighPt",      &njetsSShighPt_,     "njetsSShighPt/I"     );
    babyTree_->Branch("jet1cor", &jet1cor_, "jet1cor/F");
    babyTree_->Branch("jet2cor", &jet2cor_, "jet2cor/F");
    babyTree_->Branch("jet3cor", &jet3cor_, "jet3cor/F");
    babyTree_->Branch("jet1pt",       &jet1pt_,      "jet1pt/F"      );
    babyTree_->Branch("jet2pt",       &jet2pt_,      "jet2pt/F"      );
    babyTree_->Branch("jet3pt",       &jet3pt_,      "jet3pt/F"      );
    babyTree_->Branch("sumjetpt",     &sumjetpt_,    "sumjetpt/F"    );
    babyTree_->Branch("sumjetpt25",     &sumjetpt25_,    "sumjetpt25/F"    );
    babyTree_->Branch("sumjetptSS30",   &sumjetptSS30_,  "sumjetptSS30/F"  );
    babyTree_->Branch("sumjetptSS",   &sumjetptSS_,  "sumjetptSS/F"  );
    babyTree_->Branch("sumjetptSShighPt",   &sumjetptSShighPt_,  "sumjetptSShighPt/F"  );
    babyTree_->Branch("jet1eta",      &jet1eta_,     "jet1eta/F"     );
    babyTree_->Branch("jet2eta",      &jet2eta_,     "jet2eta/F"     );
    babyTree_->Branch("jet3eta",      &jet3eta_,     "jet3eta/F"     );
    babyTree_->Branch("jet1phi",      &jet1phi_,     "jet1phi/F"     );
    babyTree_->Branch("jet2phi",      &jet2phi_,     "jet2phi/F"     );
    babyTree_->Branch("jet3phi",      &jet3phi_,     "jet3phi/F"     );
    babyTree_->Branch("jetmass",      &jetmass_,     "jetmass/F"     );
    babyTree_->Branch("pfmth",          &pfmth_,         "pfmth/F"         );
    babyTree_->Branch("tcmth",          &tcmth_,         "tcmth/F"         );
    babyTree_->Branch("jet1isBtag",   &jet1isBtag_,  "jet1isBtag/I"  );
    babyTree_->Branch("jet2isBtag",   &jet2isBtag_,  "jet2isBtag/I"  );
    babyTree_->Branch("jet3isBtag",   &jet3isBtag_,  "jet3isBtag/I"  );
    babyTree_->Branch("dphipfmetjet", &dphipfmetjet_,"dphipfmetjet/F");
    babyTree_->Branch("dphitcmetjet", &dphitcmetjet_,"dphitcmetjet/F");
    babyTree_->Branch("deltaphi",     &deltaphi_,    "deltaphi/F"    );

    babyTree_->Branch("jet1fmetcor", &jet1fmetcor_, "jet1fmetcor/F");
    babyTree_->Branch("mlljmetcor", &mlljmetcor_, "mlljmetcor/F");
    babyTree_->Branch("pfmetcor", &pfmetcor_, "pfmetcor/F");
    babyTree_->Branch("pfmetphicor", &pfmetphicor_, "pfmetphicor/F");


    babyTree_->Branch("ntchelbtags"   ,    &ntchelbtags_   ,      "ntchelbtags/I"   );
    babyTree_->Branch("nssvhembtags"  ,   &nssvhembtags_   ,     "nssvhembtags/I"   );
    babyTree_->Branch("nssvhetbtags"  ,   &nssvhetbtags_   ,     "nssvhetbtags/I"   );
    babyTree_->Branch("nssvhptbtags"  ,   &nssvhptbtags_   ,     "nssvhptbtags/I"   );
    babyTree_->Branch("ntchelbtags20" ,    &ntchelbtags20_ ,      "ntchelbtags20/I" );
    babyTree_->Branch("ntchelbtags30" ,    &ntchelbtags30_ ,      "ntchelbtags30/I" );
    babyTree_->Branch("ntchelbtags40" ,    &ntchelbtags40_ ,      "ntchelbtags40/I" );
    babyTree_->Branch("ntchembtags20" ,    &ntchembtags20_ ,      "ntchembtags20/I" );
    babyTree_->Branch("ntchembtags30" ,    &ntchembtags30_ ,      "ntchembtags30/I" );
    babyTree_->Branch("ntchembtags40" ,    &ntchembtags40_ ,      "ntchembtags40/I" );
    babyTree_->Branch("nssvhembtags20"  ,   &nssvhembtags20_   ,     "nssvhembtags20/I"   );
    babyTree_->Branch("nssvhembtags30"  ,   &nssvhembtags30_   ,     "nssvhembtags30/I"   );
    babyTree_->Branch("nssvhembtags40"  ,   &nssvhembtags40_   ,     "nssvhembtags40/I"   );

    babyTree_->Branch("pfmeff",       &pfmeff_,      "pfmeff/F"      );
    babyTree_->Branch("tcmeff",       &tcmeff_,      "tcmeff/F"      );
    babyTree_->Branch("intLumiPerLS", &intLumiPerLS_, "intLumiPerLS/F");

    // lepton stuff
    babyTree_->Branch("ngoodlep",   &ngoodlep_,   "ngoodlep/I"  );
    babyTree_->Branch("ngoodmus",   &ngoodmus_,   "ngoodmus/I"  );
    babyTree_->Branch("ngoodels",   &ngoodels_,   "ngoodels/I"  );
    babyTree_->Branch("ngoodlepSS", &ngoodlepSS_, "ngoodlepSS/I"  );
    babyTree_->Branch("ngoodmusSS", &ngoodmusSS_, "ngoodmusSS/I"  );
    babyTree_->Branch("ngoodelsSS", &ngoodelsSS_, "ngoodelsSS/I"  );
    babyTree_->Branch("dilpt",      &dilpt_,      "dilpt/F"     );
    babyTree_->Branch("dileta",      &dileta_,      "dileta/F"     );
    babyTree_->Branch("dilphi",     &dilphi_,     "dilphi/F"    );
    babyTree_->Branch("mass",       &mass_,       "mass/F"      );
    babyTree_->Branch("mlljj",       &mlljj_,       "mlljj/F"      );
    babyTree_->Branch("mllj",       &mllj_,       "mllj/F"      );
    babyTree_->Branch("eormu1",     &eormu1_,     "eormu1/I"    );
    babyTree_->Branch("type1",      &type1_,      "type1/I"     );
    babyTree_->Branch("ngenels",    &ngenels_,    "ngenels/I"   );
    babyTree_->Branch("ngenmus",    &ngenmus_,    "ngenmus/I"   );
    babyTree_->Branch("ngentaus",   &ngentaus_,   "ngentaus/I"  );
    babyTree_->Branch("pt1",        &pt1_,        "pt1/F"       );
    babyTree_->Branch("eta1",       &eta1_,       "eta1/F"      );
    babyTree_->Branch("phi1",       &phi1_,       "phi1/F"      );
    babyTree_->Branch("iso1",       &iso1_,       "iso1/F"      );
    babyTree_->Branch("ntiso1",     &ntiso1_,     "ntiso1/F"    );
    babyTree_->Branch("d0corr1",    &d0corr1_,    "d0corr1/F"   );
    babyTree_->Branch("d0vtx1",     &d0vtx1_,     "d0vtx1/F"    );
    babyTree_->Branch("dphipfmet1", &dphipfmet1_, "dphipfmet1/F");
    babyTree_->Branch("dphitcmet1", &dphitcmet1_, "dphitcmet1/F");
    babyTree_->Branch("drjet1",     &drjet1_,     "drjet1/F"    );
    babyTree_->Branch("mcid1",      &mcid1_,      "mcid1/I"     );
    babyTree_->Branch("mcmotherid1",&mcmotherid1_,"mcmotherid1/I");
    babyTree_->Branch("mcpt1",      &mcpt1_,      "mcpt1/F"     );
    babyTree_->Branch("mcmotherpt1",&mcmotherpt1_,"mcmotherpt1/F");
    babyTree_->Branch("mc3id1",      &mc3id1_,      "mc3id1/I"     );
    babyTree_->Branch("mc3motherid1",&mc3motherid1_,"mc3motherid1/I");
    babyTree_->Branch("mc3pt1",      &mc3pt1_,      "mc3pt1/F"     );
    babyTree_->Branch("mc3motherpt1",&mc3motherpt1_,"mc3motherpt1/F");
    babyTree_->Branch("eormu2",     &eormu2_,     "eormu2/I"    );
    babyTree_->Branch("type2",      &type2_,      "type2/I"     );
    babyTree_->Branch("pt2",        &pt2_,        "pt2/F"       );
    babyTree_->Branch("eta2",       &eta2_,       "eta2/F"      );
    babyTree_->Branch("phi2",       &phi2_,       "phi2/F"      );
    babyTree_->Branch("iso2",       &iso2_,       "iso2/F"      );
    babyTree_->Branch("ntiso2",     &ntiso2_,     "ntiso2/F"    );
    babyTree_->Branch("d0corr2",    &d0corr2_,    "d0corr2/F"   );
    babyTree_->Branch("d0vtx2",     &d0vtx2_,     "d0vtx2/F"    );
    babyTree_->Branch("dphipfmet2", &dphipfmet2_, "dphipfmet2/F");
    babyTree_->Branch("dphitcmet2", &dphitcmet2_, "dphitcmet2/F");
    babyTree_->Branch("drjet2",     &drjet2_,     "drjet2/F"    );
    babyTree_->Branch("mcid2",      &mcid2_,      "mcid2/I"     );
    babyTree_->Branch("mcmotherid2",&mcmotherid2_,"mcmotherid2/I");
    babyTree_->Branch("mcpt2",      &mcpt2_,      "mcpt2/F"     );
    babyTree_->Branch("mcmotherpt2",&mcmotherpt2_,"mcmotherpt2/F");
    babyTree_->Branch("mc3id2",      &mc3id2_,      "mc3id2/I"     );
    babyTree_->Branch("mc3motherid2",&mc3motherid2_,"mc3motherid2/I");
    babyTree_->Branch("mc3pt2",      &mc3pt2_,      "mc3pt2/F"     );
    babyTree_->Branch("mc3motherpt2",&mc3motherpt2_,"mc3motherpt2/F");
    babyTree_->Branch("mt2",        &mt2_,        "mt2/F"       );
    babyTree_->Branch("mt2j",       &mt2j_,       "mt2j/F"      );
    babyTree_->Branch("extraZveto", &extraZveto_, "extraZveto/O");
    babyTree_->Branch("trkIso1", &trkIso1_, "trkIso1/F");
    babyTree_->Branch("ecalIso1", &ecalIso1_, "ecalIso1/F");
    babyTree_->Branch("hcalIso1", &hcalIso1_, "hcalIso1/F");
    babyTree_->Branch("trkIso2", &trkIso2_, "trkIso2/F");
    babyTree_->Branch("ecalIso2", &ecalIso2_, "ecalIso2/F");
    babyTree_->Branch("hcalIso2", &hcalIso2_, "hcalIso2/F");
    babyTree_->Branch("ecalIso1ps", &ecalIso1ps_, "ecalIso1ps/F");
    babyTree_->Branch("ecalIso2ps", &ecalIso2ps_, "ecalIso2ps/F");
    babyTree_->Branch("rho", &rho_, "rho/F");
    babyTree_->Branch("lepsFromSameVtx", &lepsFromSameVtx_, "lepsFromSameVtx/O");
    babyTree_->Branch("lep1isFromW", &lep1isFromW_, "lep1isFromW/I");
    babyTree_->Branch("lep2isFromW", &lep2isFromW_, "lep2isFromW/I");
    babyTree_->Branch("lep1_pfiso0p3", &lep1_pfiso0p3_, "lep1_pfiso0p3/F");
    babyTree_->Branch("lep1_pfiso0p4", &lep1_pfiso0p4_, "lep1_pfiso0p4/F");
    babyTree_->Branch("lep2_pfiso0p3", &lep2_pfiso0p3_, "lep2_pfiso0p3/F");
    babyTree_->Branch("lep2_pfiso0p4", &lep2_pfiso0p4_, "lep2_pfiso0p4/F");

    // Muon stuff
    babyTree_->Branch("mu1_numSSv4",      &mu1_numSSv4_,      "mu1_numSSv4/O"     );
    babyTree_->Branch("mu1_numSSv4noIso",      &mu1_numSSv4noIso_,      "mu1_numSSv4noIso/O"     );
    babyTree_->Branch("mu1_foSSv4",       &mu1_foSSv4_,       "mu1_foSSv4/O"      );
    babyTree_->Branch("mu1_foSSv4noIso",       &mu1_foSSv4noIso_,       "mu1_foSSv4noIso/O"      );
    babyTree_->Branch("mu1_numSSv3",      &mu1_numSSv3_,      "mu1_numSSv3/O"     );
    babyTree_->Branch("mu1_foSSv3",       &mu1_foSSv3_,       "mu1_foSSv3/O"      );
    babyTree_->Branch("mu1_muonidfull",   &mu1_muonidfull_,   "mu1_muonidfull/O"  );
    babyTree_->Branch("mu1_muonid",       &mu1_muonid_,       "mu1_muonid/O"      );
    babyTree_->Branch("mu1_muonidfullV1", &mu1_muonidfullV1_, "mu1_muonidfullV1/O");
    babyTree_->Branch("mu1_muonidV1",     &mu1_muonidV1_,     "mu1_muonidV1/O"    );
    babyTree_->Branch("mu1_goodmask",     &mu1_goodmask_,     "mu1_goodmask/I"    );
    babyTree_->Branch("mu1_gfitchi2",     &mu1_gfitchi2_,     "mu1_gfitchi2/F"    );
    babyTree_->Branch("mu1_cosmic",       &mu1_cosmic_,       "mu1_cosmic/O"      );
    babyTree_->Branch("mu1_siHits",       &mu1_siHits_,       "mu1_siHits/I"      );
    babyTree_->Branch("mu1_saHits",       &mu1_saHits_,       "mu1_saHits/I"      );
    babyTree_->Branch("mu1_emVetoDep",    &mu1_emVetoDep_,    "mu1_emVetoDep/F"   );
    babyTree_->Branch("mu1_hadVetoDep",   &mu1_hadVetoDep_,   "mu1_hadVetoDep/F"  );
    babyTree_->Branch("mu1_isPFmuon",     &mu1_isPFmuon_,     "mu1_isPFmuon/O"    );
    babyTree_->Branch("mu1_relPtErr",     &mu1_relPtErr_,     "mu1_relPtErr/F"    );
    babyTree_->Branch("mu1_trkpt",     &mu1_trkpt_,     "mu1_trkpt/F"    );

    babyTree_->Branch("mu2_numSSv4",      &mu2_numSSv4_,      "mu2_numSSv4/O"     );
    babyTree_->Branch("mu2_numSSv4noIso",      &mu2_numSSv4noIso_,      "mu2_numSSv4noIso/O"     );
    babyTree_->Branch("mu2_foSSv4",       &mu2_foSSv4_,       "mu2_foSSv4/O"      );
    babyTree_->Branch("mu2_foSSv4noIso",       &mu2_foSSv4noIso_,       "mu2_foSSv4noIso/O"      );
    babyTree_->Branch("mu2_numSSv3",      &mu2_numSSv3_,      "mu2_numSSv3/O"     );
    babyTree_->Branch("mu2_foSSv3",       &mu2_foSSv3_,       "mu2_foSSv3/O"      );
    babyTree_->Branch("mu2_muonidfull",   &mu2_muonidfull_,   "mu2_muonidfull/O"  );
    babyTree_->Branch("mu2_muonid",       &mu2_muonid_,       "mu2_muonid/O"      );
    babyTree_->Branch("mu2_muonidfullV1", &mu2_muonidfullV1_, "mu2_muonidfullV1/O");
    babyTree_->Branch("mu2_muonidV1",     &mu2_muonidV1_,     "mu2_muonidV1/O"    );
    babyTree_->Branch("mu2_goodmask",     &mu2_goodmask_,     "mu2_goodmask/I"    );
    babyTree_->Branch("mu2_gfitchi2",     &mu2_gfitchi2_,     "mu2_gfitchi2/F"    );
    babyTree_->Branch("mu2_cosmic",       &mu2_cosmic_,       "mu2_cosmic/O"      );
    babyTree_->Branch("mu2_siHits",       &mu2_siHits_,       "mu2_siHits/I"      );
    babyTree_->Branch("mu2_saHits",       &mu2_saHits_,       "mu2_saHits/I"      );
    babyTree_->Branch("mu2_emVetoDep",    &mu2_emVetoDep_,    "mu2_emVetoDep/F"   );
    babyTree_->Branch("mu2_hadVetoDep",   &mu2_hadVetoDep_,   "mu2_hadVetoDep/F"  );
    babyTree_->Branch("mu2_isPFmuon",     &mu2_isPFmuon_,     "mu2_isPFmuon/O"    );
    babyTree_->Branch("mu2_relPtErr",     &mu2_relPtErr_,     "mu2_relPtErr/F"    );
    babyTree_->Branch("mu2_trkpt",     &mu2_trkpt_,     "mu2_trkpt/F"    );

    // electron stuff

    babyTree_->Branch("e1_numSSv5",    &e1_numSSv5_,    "e1_numSSv5/O"   );
    babyTree_->Branch("e1_numSSv5noIso",    &e1_numSSv5noIso_,    "e1_numSSv5noIso/O"   );
    babyTree_->Branch("e1_foSSv5",     &e1_foSSv5_,     "e1_foSSv5/O"    );
    babyTree_->Branch("e1_numSSv5_noConvCuts",    &e1_numSSv5_noConvCuts_,    "e1_numSSv5_noConvCuts/O"   );
    babyTree_->Branch("e1_numSSv5noIso_noConvCuts",    &e1_numSSv5noIso_noConvCuts_,    "e1_numSSv5noIso_noConvCuts/O"   );
    babyTree_->Branch("e1_foSSv5_noConvCuts",     &e1_foSSv5_noConvCuts_,     "e1_foSSv5_noConvCuts/O"    );
    babyTree_->Branch("e1_numSSv4",    &e1_numSSv4_,    "e1_numSSv4/O"   );
    babyTree_->Branch("e1_foSSv4",     &e1_foSSv4_,     "e1_foSSv4/O"    );
    babyTree_->Branch("e1_numSSv3",    &e1_numSSv3_,    "e1_numSSv3/O"   );
    babyTree_->Branch("e1_foSSv3",     &e1_foSSv3_,     "e1_foSSv3/O"    );
    babyTree_->Branch("e1_vbtf90full", &e1_vbtf90full_, "e1_vbtf90full/O");
    babyTree_->Branch("e1_vbtf90",     &e1_vbtf90_,     "e1_vbtf90/O"    );
    babyTree_->Branch("e1_vbtf85",     &e1_vbtf85_,     "e1_vbtf85/O"    );
    babyTree_->Branch("e1_vbtf80",     &e1_vbtf80_,     "e1_vbtf80/O"    );
    babyTree_->Branch("e1_vbtf70",     &e1_vbtf70_,     "e1_vbtf70/O"    );
    babyTree_->Branch("e1_smurfV3",    &e1_smurfV3_,    "e1_smurfV3/O"   );
    babyTree_->Branch("e1_scet",       &e1_scet_,       "e1_scet/F"      );
    babyTree_->Branch("e1_sceta",       &e1_sceta_,       "e1_sceta/F"      );
    babyTree_->Branch("e1_eopin",      &e1_eopin_,      "e1_eopin/F"     );
    babyTree_->Branch("e1_hoe",        &e1_hoe_,        "e1_hoe/F"       );
    babyTree_->Branch("e1_dphiin",     &e1_dphiin_,     "e1_dphiin/F"    );
    babyTree_->Branch("e1_detain",     &e1_detain_,     "e1_detain/F"    );
    babyTree_->Branch("e1_e25Me55",    &e1_e25Me55_,    "e1_e25Me55/F"   );
    babyTree_->Branch("e1_e15Me55",    &e1_e15Me55_,    "e1_e15Me55/F"   );
    babyTree_->Branch("e1_sigieie",    &e1_sigieie_,    "e1_sigieie/F"   );
    babyTree_->Branch("e1_eMe55",      &e1_eMe55_,      "e1_eMe55/F"     ); // for spikes
    babyTree_->Branch("e1_nmHits",     &e1_nmHits_,     "e1_nmHits/I"    );
    babyTree_->Branch("e1_dcot",       &e1_dcot_,       "e1_dcot/F"      );
    babyTree_->Branch("e1_dist",       &e1_dist_,       "e1_dist/F"      );
    babyTree_->Branch("e1_dcot_old",   &e1_dcot_old_,   "e1_dcot_old/F"      );
    babyTree_->Branch("e1_dist_old",   &e1_dist_old_,   "e1_dist_old/F"      );
    babyTree_->Branch("e1_drmu",       &e1_drmu_,       "e1_drmu/F"      );
    babyTree_->Branch("e1_drmuSS",     &e1_drmuSS_,     "e1_drmuSS/F"    );
    babyTree_->Branch("e1_isspike",    &e1_isspike_,    "e1_isspike/O"   );
    babyTree_->Branch("e1_ctfCharge",  &e1_ctfCharge_,  "e1_ctfCharge/I" );
    babyTree_->Branch("e1_gsfCharge",  &e1_gsfCharge_,  "e1_gsfCharge/I" );
    babyTree_->Branch("e1_scCharge",   &e1_scCharge_,   "e1_scCharge/I"  );
    babyTree_->Branch("e1_fbrem",      &e1_fbrem_,      "e1_fbrem/F"     );
    babyTree_->Branch("e1_mitConv",    &e1_mitConv_,    "e1_mitConv/O"   );
    babyTree_->Branch("e1_eSeedOverPout", &e1_eSeedOverPout_, "e1_eSeedOverPout/F");
    babyTree_->Branch("e1_eSeedOverPin", &e1_eSeedOverPin_, "e1_eSeedOverPin/F");
    babyTree_->Branch("e1_sigipip", &e1_sigipip_, "e1_sigipip/F");
    babyTree_->Branch("e1_nbrem", &e1_nbrem_, "e1_nbrem/F");
    babyTree_->Branch("e1_ip3d", &e1_ip3d_, "e1_ip3d/F");
    babyTree_->Branch("e1_ip3dSig", &e1_ip3dSig_, "e1_ip3dSig/F");
    babyTree_->Branch("e1_lh", &e1_lh_, "e1_lh/F");
    babyTree_->Branch("e1_oneOverEminusOneOverP", &e1_oneOverEminusOneOverP_, "e1_oneOverEminusOneOverP/F");
    babyTree_->Branch("e1_mva", &e1_mva_, "e1_mva/F");
    babyTree_->Branch("e2_numSSv5",    &e2_numSSv5_,    "e2_numSSv5/O"   );
    babyTree_->Branch("e2_numSSv5noIso",    &e2_numSSv5noIso_,    "e2_numSSv5noIso/O"   );
    babyTree_->Branch("e2_foSSv5",     &e2_foSSv5_,     "e2_foSSv5/O"    );
    babyTree_->Branch("e2_numSSv5_noConvCuts",    &e2_numSSv5_noConvCuts_,    "e2_numSSv5_noConvCuts/O"   );
    babyTree_->Branch("e2_numSSv5noIso_noConvCuts",    &e2_numSSv5noIso_noConvCuts_,    "e2_numSSv5noIso_noConvCuts/O"   );
    babyTree_->Branch("e2_foSSv5_noConvCuts",     &e2_foSSv5_noConvCuts_,     "e2_foSSv5_noConvCuts/O"    );
    babyTree_->Branch("e2_numSSv4",    &e2_numSSv4_,    "e2_numSSv4/O"   );
    babyTree_->Branch("e2_foSSv4",     &e2_foSSv4_,     "e2_foSSv4/O"    );
    babyTree_->Branch("e2_numSSv3",    &e2_numSSv3_,    "e2_numSSv3/O"   );
    babyTree_->Branch("e2_foSSv3",     &e2_foSSv3_,     "e2_foSSv3/O"    );
    babyTree_->Branch("e2_vbtf90full", &e2_vbtf90full_, "e2_vbtf90full/O");
    babyTree_->Branch("e2_vbtf90",     &e2_vbtf90_,     "e2_vbtf90/O"    );
    babyTree_->Branch("e2_vbtf85",     &e2_vbtf85_,     "e2_vbtf85/O"    );
    babyTree_->Branch("e2_vbtf80",     &e2_vbtf80_,     "e2_vbtf80/O"    );
    babyTree_->Branch("e2_vbtf70",     &e2_vbtf70_,     "e2_vbtf70/O"    );
    babyTree_->Branch("e2_smurfV3",    &e2_smurfV3_,    "e2_smurfV3/O"   );
    babyTree_->Branch("e2_scet",       &e2_scet_,       "e2_scet/F"      );
    babyTree_->Branch("e2_sceta",      &e2_sceta_,      "e2_sceta/F"     );
    babyTree_->Branch("e2_eopin",      &e2_eopin_,      "e2_eopin/F"     );
    babyTree_->Branch("e2_hoe",        &e2_hoe_,        "e2_hoe/F"       );
    babyTree_->Branch("e2_dphiin",     &e2_dphiin_,     "e2_dphiin/F"    );
    babyTree_->Branch("e2_detain",     &e2_detain_,     "e2_detain/F"    );
    babyTree_->Branch("e2_e25Me55",    &e2_e25Me55_,    "e2_e25Me55/F"   );
    babyTree_->Branch("e2_e15Me55",    &e2_e15Me55_,    "e2_e15Me55/F"   );
    babyTree_->Branch("e2_sigieie",    &e2_sigieie_,    "e2_sigieie/F"   );
    babyTree_->Branch("e2_eMe55",      &e2_eMe55_,      "e2_eMe55/F"     ); // for spikes
    babyTree_->Branch("e2_nmHits",     &e2_nmHits_,     "e2_nmHits/I"    );
    babyTree_->Branch("e2_dcot",       &e2_dcot_,       "e2_dcot/F"      );
    babyTree_->Branch("e2_dist",       &e2_dist_,       "e2_dist/F"      );
    babyTree_->Branch("e2_dcot_old",   &e2_dcot_old_,   "e2_dcot_old/F"  );
    babyTree_->Branch("e2_dist_old",   &e2_dist_old_,   "e2_dist_old/F"  );
    babyTree_->Branch("e2_drmu",       &e2_drmu_,       "e2_drmu/F"      );
    babyTree_->Branch("e2_drmuSS",     &e2_drmuSS_,     "e2_drmuSS/F"    );
    babyTree_->Branch("e2_isspike",    &e2_isspike_,    "e2_isspike/O"   );
    babyTree_->Branch("e2_ctfCharge",  &e2_ctfCharge_,  "e2_ctfCharge/I" );
    babyTree_->Branch("e2_gsfCharge",  &e2_gsfCharge_,  "e2_gsfCharge/I" );
    babyTree_->Branch("e2_scCharge",   &e2_scCharge_,   "e2_scCharge/I"  );
    babyTree_->Branch("e2_fbrem",      &e2_fbrem_,      "e2_fbrem/F"     );
    babyTree_->Branch("e2_mitConv",    &e2_mitConv_,    "e2_mitConv/O"   );
    babyTree_->Branch("e2_eSeedOverPout", &e2_eSeedOverPout_, "e2_eSeedOverPout/F");
    babyTree_->Branch("e2_eSeedOverPin", &e2_eSeedOverPin_, "e2_eSeedOverPin/F");
    babyTree_->Branch("e2_sigipip", &e2_sigipip_, "e2_sigipip/F");
    babyTree_->Branch("e2_nbrem", &e2_nbrem_, "e2_nbrem/F");
    babyTree_->Branch("e2_ip3d", &e2_ip3d_, "e2_ip3d/F");
    babyTree_->Branch("e2_ip3dSig", &e2_ip3dSig_, "e2_ip3dSig/F");
    babyTree_->Branch("e2_lh", &e2_lh_, "e2_lh/F");
    babyTree_->Branch("e2_oneOverEminusOneOverP", &e2_oneOverEminusOneOverP_, "e2_oneOverEminusOneOverP/F");
    babyTree_->Branch("e2_mva", &e2_mva_, "e2_mva/F");

    // trigger stuff
    babyTree_->Branch("pass_trg_single_mu1",   &pass_trg_single_mu1_,   "pass_trg_single_mu1/I"  );
    babyTree_->Branch("pass_trg_single_mu2",   &pass_trg_single_mu2_,   "pass_trg_single_mu2/I"  );

    babyTree_->Branch("pass_trg_single_e1",   &pass_trg_single_e1_,   "pass_trg_single_e1/I"  );
    babyTree_->Branch("pass_trg_single_e2",   &pass_trg_single_e2_,   "pass_trg_single_e2/I"  );

    babyTree_->Branch("pass_trg_double_mu1",   &pass_trg_double_mu1_,   "pass_trg_double_mu1/I"  );
    babyTree_->Branch("pass_trg_double_mu2",   &pass_trg_double_mu2_,   "pass_trg_double_mu2/I"  );

    babyTree_->Branch("pass_trg_double_e1",   &pass_trg_double_e1_,   "pass_trg_double_e1/I"  );
    babyTree_->Branch("pass_trg_double_e2",   &pass_trg_double_e2_,   "pass_trg_double_e2/I"  );

    babyTree_->Branch("pass_trg_cross_emu",   &pass_trg_cross_emu_,   "pass_trg_cross_emu/I"  );

    babyTree_->Branch("pass_trg_had_double_e1" , &pass_trg_had_double_e1_ , "pass_trg_had_double_e1/I" );
    babyTree_->Branch("pass_trg_had_double_e2" , &pass_trg_had_double_e2_ , "pass_trg_had_double_e2/I" );

    babyTree_->Branch("pass_trg_had_double_mu1", &pass_trg_had_double_mu1_, "pass_trg_had_double_mu1/I");
    babyTree_->Branch("pass_trg_had_double_mu2", &pass_trg_had_double_mu2_, "pass_trg_had_double_mu2/I");
    babyTree_->Branch("pass_trg_had_cross_emu", &pass_trg_had_cross_emu_, "pass_trg_had_cross_emu/I");

    // trigger stuff
    babyTree_->Branch("match_trg_single_mu1",   &match_trg_single_mu1_,   "match_trg_single_mu1/I"  );
    babyTree_->Branch("match_trg_single_mu2",   &match_trg_single_mu2_,   "match_trg_single_mu2/I"  );

    babyTree_->Branch("match_trg_single_e1",   &match_trg_single_e1_,   "match_trg_single_e1/I"  );
    babyTree_->Branch("match_trg_single_e2",   &match_trg_single_e2_,   "match_trg_single_e2/I"  );

    babyTree_->Branch("match_trg_double_mu1",   &match_trg_double_mu1_,   "match_trg_double_mu1/I"  );
    babyTree_->Branch("match_trg_double_mu2",   &match_trg_double_mu2_,   "match_trg_double_mu2/I"  );

    babyTree_->Branch("match_trg_double_e1",   &match_trg_double_e1_,   "match_trg_double_e1/I"  );
    babyTree_->Branch("match_trg_double_e2",   &match_trg_double_e2_,   "match_trg_double_e2/I"  );

    babyTree_->Branch("match_trg_cross_emu",   &match_trg_cross_emu_,   "match_trg_cross_emu/I"  );

    babyTree_->Branch("match_trg_had_double_e1" , &match_trg_had_double_e1_ , "match_trg_had_double_e1/I" );
    babyTree_->Branch("match_trg_had_double_e2" , &match_trg_had_double_e2_ , "match_trg_had_double_e2/I" );

    babyTree_->Branch("match_trg_had_double_mu1", &match_trg_had_double_mu1_, "match_trg_had_double_mu1/I");
    babyTree_->Branch("match_trg_had_double_mu2", &match_trg_had_double_mu2_, "match_trg_had_double_mu2/I");
    babyTree_->Branch("match_trg_had_cross_emu", &match_trg_had_cross_emu_, "match_trg_had_cross_emu/I");
}

bool dilepbabymaker::PassTriggerGroup(const std::vector<std::pair<std::string, unsigned int> > &triggers, const LorentzVector &obj)
{
    for (unsigned int i = 0; i < triggers.size(); ++i) {
        if (passUnprescaledHLTTrigger(triggers[i].first.c_str(), obj)) return true;
    }
    return false;
}

void dilepbabymaker::PassTriggerGroup(const std::vector<std::pair<std::string, unsigned int> > &triggers, const LorentzVector &obj, Int_t &mask)
{
    for (unsigned int i = 0; i < triggers.size(); ++i) {
        if (passUnprescaledHLTTrigger(triggers[i].first.c_str(), obj)) mask |= (1<<triggers[i].second);
    }    

}  

void dilepbabymaker::PassTriggerGroup(const std::vector<std::pair<std::string, unsigned int> > &triggers, Int_t &mask)
{
    for (unsigned int i = 0; i < triggers.size(); ++i) {
        if (passUnprescaledHLTTrigger(triggers[i].first.c_str())) mask |= (1<<triggers[i].second);
    }

}

void dilepbabymaker::SetEventLevelInfo ()
{

    strcpy(dataset_, cms2.evt_dataset().Data());
    run_        = cms2.evt_run();
    ls_         = cms2.evt_lumiBlock();
    evt_        = cms2.evt_event();
    isdata_     = cms2.evt_isRealData();

    evt_clean082010_ = cleaning_standardAugust2010(isdata_);
    evt_clean102010_ = cleaning_standardOctober2010();
    evt_clean042011_ = cleaning_standardApril2011();

    nvtx_ = 0;
    for (unsigned int i=0; i<cms2.vtxs_isFake().size(); i++) {
        // isGoodVertex is defined in CORE/eventSelections.cc
        // !isFake, ndof >= 4, rho <= 2, Z <= 24.
        if (isGoodVertex(i))
            ++nvtx_;
    }

    if (!isdata_) {
        int nlep  = leptonGenpCount_lepTauDecays(ngenels_, ngenmus_, ngentaus_);
        scale1fb_ = cms2.evt_scale1fb();
        pthat_    = cms2.genps_qScale();
        genmet_   = cms2.gen_met();
        unsigned int ngjets = 0;
        float ght = 0.;
        for (unsigned int gidx = 0; gidx < cms2.genps_p4().size(); gidx++) {
            if (cms2.genps_status().at(gidx) != 3)
                continue;
            if (fabs(cms2.genps_p4().at(gidx).eta()) > 2.5)
                continue;
            if (cms2.genps_p4().at(gidx).pt() < 40.)
                continue;
            if ((abs(cms2.genps_id().at(gidx)) < 1 || abs(cms2.genps_id().at(gidx)) > 5) && abs(cms2.genps_id().at(gidx)) != 21)
                continue;

            ++ngjets;
            ght += cms2.genps_p4().at(gidx).pt();
        }
        ngenjets_ = ngjets;
        gensumjetpt_ = ght;
    }

    ntrks_      = cms2.trks_trk_p4().size();

    rho_        = cms2.evt_rho();
}

void dilepbabymaker::FillBabyNtuple()
{
    babyTree_->Fill();
}

void dilepbabymaker::CloseBabyNtuple()
{
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();
    delete babyFile_;
}

dilepbabymaker::dilepbabymaker()
{
    electronIdMVA = new ElectronIDMVA();
    std::string cms2_location = "../..";
    char* ppath = getenv("CMS2_LOCATION");
    if (ppath != NULL) {
        cms2_location = ppath;
        cms2_location += "/NtupleMacros";
    }
    electronIdMVA->Initialize("BDTG method", 2,
                             cms2_location + "/Tools/EgammaAnalysisTools/data/Subdet0LowPt_WithIPInfo_BDTG.weights.xml",
                             cms2_location + "/Tools/EgammaAnalysisTools/data/Subdet1LowPt_WithIPInfo_BDTG.weights.xml",
                             cms2_location + "/Tools/EgammaAnalysisTools/data/Subdet2LowPt_WithIPInfo_BDTG.weights.xml",
                             cms2_location + "/Tools/EgammaAnalysisTools/data/Subdet0HighPt_WithIPInfo_BDTG.weights.xml",
                             cms2_location + "/Tools/EgammaAnalysisTools/data/Subdet1HighPt_WithIPInfo_BDTG.weights.xml",
                             cms2_location + "/Tools/EgammaAnalysisTools/data/Subdet2HighPt_WithIPInfo_BDTG.weights.xml");
}
