#include "babymaker.h" 
#include <algorithm>
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TROOT.h"

#include "CORE/CMS2.cc"
#include "CORE/electronSelections.cc"
#include "CORE/muonSelections.cc"

#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;
float deltaPhi (float, float);
bool sortByPt (const LorentzVector &, const LorentzVector &);

void babymaker::ScanChain (TChain *chain, char *babyFilename, int nEvents)
{
    if (chain->GetEntries() < 1)
    {
        std::cout<<"babymaker: no entries in chain, exiting..."<<std::endl;
        return;
    }
    TObjArray *listOfFiles = chain->GetListOfFiles();

    unsigned int nEventsChain=0;
    if (nEvents==-1) 
        nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;

    // make a baby ntuple
    MakeBabyNtuple(babyFilename);
    //BookHists();

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
        for(unsigned int event = 0; event < nEvents; ++event)
        {
            cms2.GetEntry(event);
            ++nEventsTotal;
            // Progress feedback to the user
            if(nEventsTotal%1000 == 0)
            {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1))
                {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                            "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
                    fflush(stdout);
                }
            }

            //
            // N.B. BABY NTUPLE IS FILLED
            // FOR EACH MUON/ELECTRON
            //
            InitBabyNtuple();

            // event stuff
            run_        = cms2.evt_run();
            ls_         = cms2.evt_lumiBlock();
            evt_        = cms2.evt_event();
            pfmet_      = cms2.evt_pfmet();

            float theMetPhi = cms2.evt_pfmetPhi();

            VofP4 theJets;
            for(unsigned int jeti = 0; jeti < cms2.pfjets_p4().size(); ++jeti)
            {
                if (cms2.pfjets_p4()[jeti].pt() > 20.)
                    theJets.push_back(cms2.pfjets_p4()[jeti]);
            }
            std::sort(theJets.begin(), theJets.end(), sortByPt);

            double mindphimet = 999999.;
            for(unsigned int jeti = 0; jeti < theJets.size(); ++jeti)
            {
                float currdphimet = deltaPhi(theMetPhi, theJets[jeti].phi());
                if (currdphimet < mindphimet)
                    mindphimet = currdphimet;
            }

            njets_      = theJets.size();
            jet1pt_     = theJets.size() ? theJets[0].pt() : -999999.;
            dphimetjet_ = mindphimet;

            // muon stuff
            for(unsigned mui = 0; mui < cms2.mus_p4().size(); ++mui)
            {
                // global and tracker muons only
                if (! (cms2.mus_type()[mui]&6))
                    continue;
                // pt > 20
                if (cms2.mus_p4()[mui].pt() <= 20)
                    continue;

                eormu_    = 13;
                type_     = cms2.mus_type()[mui];
                pt_       = cms2.mus_p4()[mui].pt();
                iso_      = muonIsoValue(mui);
                d0corr_   = cms2.mus_d0corr()[mui];
                dphimet_  = deltaPhi(theMetPhi, cms2.mus_p4()[mui].phi());

                float mindrjet = 999999.;
                for(unsigned int jeti = 0; jeti < theJets.size(); ++jeti)
                {
                    float deta = cms2.mus_p4()[mui].eta()-theJets[jeti].eta();
                    float dphi = deltaPhi(cms2.mus_p4()[mui].phi(), theJets[jeti].phi());
                    float currdrjet = sqrt(deta*deta+dphi*dphi);
                    if (currdrjet < mindrjet)
                        mindrjet = currdrjet;
                }

                drjet_       = mindrjet;
                mu_muonid_   = muonId(mui);
                mu_goodmask_ = cms2.mus_goodmask()[mui];
                mu_gfitchi2_ = cms2.mus_gfit_chi2()[mui];

                FillBabyNtuple();
            }

            // electron stuff
            for(unsigned eli = 0; eli < cms2.els_p4().size(); ++eli)
            {
                // pt > 20
                if (cms2.els_p4()[eli].pt() <= 20)
                    continue;

                eormu_   = 11;
                type_    = cms2.els_type()[eli];
                pt_      = cms2.els_p4()[eli].pt();
                iso_     = electronIsolation_relsusy_cand0(eli, true);
                d0corr_  = cms2.els_d0corr()[eli];
                dphimet_ = deltaPhi(theMetPhi, cms2.els_p4()[eli].phi());

                float mindrjet = 999999.;
                for(unsigned int jeti = 0; jeti < theJets.size(); ++jeti)
                {
                    float deta = cms2.els_p4()[eli].eta()-theJets[jeti].eta();
                    float dphi = deltaPhi(cms2.els_p4()[eli].phi(), theJets[jeti].phi());
                    float currdrjet = sqrt(deta*deta+dphi*dphi);
                    if (currdrjet < mindrjet)
                        mindrjet = currdrjet;
                }

                drjet_    = mindrjet;
                e_cand01_ = electronId_cand01(eli);
                e_eopin_  = cms2.els_eOverPIn()[eli];
                e_hoe_    = cms2.els_hOverE()[eli];
                e_dphiin_ = cms2.els_dPhiIn()[eli];
                e_detain_ = cms2.els_dEtaIn()[eli];
                e_eMe55_  = cms2.els_eMax()[eli]/cms2.els_e5x5()[eli];
                e_nmHits_ = cms2.els_exp_innerlayers()[eli];

                FillBabyNtuple();
            }
        }

        if (nEventsChain != nEventsTotal)
        {
            std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
        }
    }

    CloseBabyNtuple();
}

void babymaker::InitBabyNtuple ()
{
    // event stuff
    run_          = -999999;
    ls_           = -999999;
    evt_          = -999999;
    pfmet_        = -999999.;
    njets_        = -999999;
    jet1pt_       = -999999.;
    dphimetjet_   = -999999.;
    eormu_        = -999999;

    // lepton stuff
    type_         = -999999;
    pt_           = -999999.;
    iso_          = -999999.;
    d0corr_       = -999999.;
    dphimet_      = -999999.;
    drjet_        = -999999.;

    // muon stuff
    mu_muonid_   = -999999;
    mu_goodmask_ = -999999;
    mu_gfitchi2_ = -999999.;

    // electron stuff
    e_cand01_   = -999999;
    e_eopin_    = -999999.;
    e_hoe_      = -999999.;
    e_dphiin_   = -999999.;
    e_detain_   = -999999.;
    e_eMe55_    = -999999.;
    e_nmHits_   = -999999;
}

void babymaker::MakeBabyNtuple(char *babyFilename)
{
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();

    babyFile_ = new TFile(Form("%s", babyFilename), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");

    // event stuff
    babyTree_->Branch("run",        &run_,       "run/I"       );
    babyTree_->Branch("ls",         &ls_,        "ls/I"        );
    babyTree_->Branch("evt",        &evt_,       "evt/I"       );
    babyTree_->Branch("pfmet",      &pfmet_,     "pfmet/F"     );
    babyTree_->Branch("njets",      &njets_,     "njets/I"     ); // uncorrected pt > 20
    babyTree_->Branch("jet1pt",     &jet1pt_,    "jet1pt/F"    );
    babyTree_->Branch("dphimetjet", &dphimetjet_,"dphimetjet/F");
    babyTree_->Branch("eormu",      &eormu_,     "eormu/I"     );

    // lepton stuff
    babyTree_->Branch("type",    &type_,    "type/I"   );
    babyTree_->Branch("pt",      &pt_,      "pt/F"     );
    babyTree_->Branch("iso",     &iso_,     "iso/F"    );
    babyTree_->Branch("d0corr",  &d0corr_,  "d0corr/F" );
    babyTree_->Branch("dphimet", &dphimet_, "dphimet/F");
    babyTree_->Branch("drjet",   &drjet_,   "drjet/F"  );

    // muon stuff
    babyTree_->Branch("mu_muonid",   &mu_muonid_,   "mu_muonid/O"  );
    babyTree_->Branch("mu_goodmask", &mu_goodmask_, "mu_goodmask/I");
    babyTree_->Branch("mu_gfitchi2", &mu_gfitchi2_, "mu_gfitchi2/F");

    // eectron stuff
    babyTree_->Branch("e_cand01", &e_cand01_, "e_cand01/O");
    babyTree_->Branch("e_eopin",  &e_eopin_,  "e_eopin/F" );
    babyTree_->Branch("e_hoe",    &e_hoe_,    "e_hoe/F"   );
    babyTree_->Branch("e_dphiin", &e_dphiin_, "e_dphiin/F");
    babyTree_->Branch("e_detain", &e_detain_, "e_detain/F");
    babyTree_->Branch("e_eMe55",  &e_eMe55_,  "e_eMe55/F" ); // for spikes
    babyTree_->Branch("e_nmHits", &e_nmHits_, "e_nmHits/I"); // for conversions
}

void babymaker::FillBabyNtuple()
{
    babyTree_->Fill();
}

void babymaker::CloseBabyNtuple()
{
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();
}

float deltaPhi (float phi1, float phi2)
{
    float dphi = phi1-phi2;
    while (dphi < 0)
        dphi += TMath::Pi();
    while (dphi > TMath::Pi())
        dphi -= TMath::Pi();
    return dphi;
}

bool sortByPt (const LorentzVector &vec1, const LorentzVector &vec2)
{
    return vec1.pt() > vec2.pt();
}
