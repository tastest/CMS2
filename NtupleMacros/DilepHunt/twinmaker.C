#include "twinmaker.h" 
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

#include "CORE/CMS2.h"
#include "CORE/electronSelections.cc"
#include "CORE/muonSelections.cc"
#include "CORE/metSelections.cc"
#include "CORE/trackSelections.cc"

#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;
float deltaPhi (float, float);
bool sortByPt (const LorentzVector &, const LorentzVector &);
bool isGoodElectron(const int);

void twinmaker::ScanChain (const char *inputFilename, const char *twinFilename, int nEvents)
{
    TChain *chain = new TChain("Events");
    chain->Add(inputFilename);
    TObjArray *listOfFiles = chain->GetListOfFiles();    

    unsigned int nEventsChain=0;
    if (nEvents==-1) 
        nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;

    // make a twin ntuple
    MakeTwinNtuple(twinFilename);

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
            // N.B. TWIN NTUPLE IS FILLED
            // FOR EACH MUON/ELECTRON PAIR
            //
            InitTwinNtuple();

            // event stuff
            run_        = cms2.evt_run();
            ls_         = cms2.evt_lumiBlock();
            evt_        = cms2.evt_event();
            pfmet_      = cms2.evt_pfmet();

            metStruct tcMET = correctedTCMET();
            float tcmet    = tcMET.met;
            float tcmetPhi = tcMET.metphi;

            // now, need to loop over calotwers and correct tcMET for HF spikes
            for(unsigned int i = 0; i < cms2.twrs_emEnergy().size(); ++i)
            {
                // if tower is in the HF, check if it is a spike
                if (fabs(cms2.twrs_eta().at(i)) < 3.)
                    continue;

                double alpha = cms2.twrs_emEnergy().at(i) / (cms2.twrs_emEnergy().at(i) + cms2.twrs_hadEnergy().at(i)); // alpha = (L-S)/(L+S)

                // if tower is "spiking" correct MET by adding components of tower to components of MET
                if ( (alpha <= -0.8 || alpha >= 0.99) && (cms2.twrs_emEt().at(i) + cms2.twrs_hadEt().at(i)) > 5. )
                {

                    const float towerET  = cms2.twrs_emEt().at(i) + cms2.twrs_hadEt().at(i);			      
                    const float hfspikeX = towerET * cos(cms2.twrs_phi().at(i));
                    const float hfspikeY = towerET * sin(cms2.twrs_phi().at(i));

                    const float tcmetCorrX = tcmet * cos(tcmetPhi) + hfspikeX;
                    const float tcmetCorrY = tcmet * sin(tcmetPhi) + hfspikeY;
                    tcmet    = sqrt( tcmetCorrX * tcmetCorrX + tcmetCorrY * tcmetCorrY );
                    tcmetPhi = atan2(tcmetCorrY, tcmetCorrX);
                }
            }

            tcmet_ = tcmet;

            float thePFMetPhi = cms2.evt_pfmetPhi();
            float theTCMetPhi = tcmetPhi;

            VofP4 theJets;
            for(unsigned int jeti = 0; jeti < cms2.pfjets_p4().size(); ++jeti)
            {
                if (cms2.pfjets_p4()[jeti].pt() > 20.)
                    theJets.push_back(cms2.pfjets_p4()[jeti]);
            }
            std::sort(theJets.begin(), theJets.end(), sortByPt);

            njets_      = theJets.size();
            jet1pt_     = theJets.size() ? theJets[0].pt() : -999999.;

            double mindphipfmet = 999999.;
            double mindphitcmet = 999999.;
            for(unsigned int jeti = 0; jeti < theJets.size(); ++jeti)
            {
                float currdphipfmet = deltaPhi(thePFMetPhi, theJets[jeti].phi());
                if (currdphipfmet < mindphipfmet)
                    mindphipfmet = currdphipfmet;

                float currdphitcmet = deltaPhi(theTCMetPhi, theJets[jeti].phi());
                if (currdphitcmet < mindphitcmet)
                    mindphitcmet = currdphitcmet;
            }

            dphipfmetjet_ = mindphipfmet;
            dphitcmetjet_ = mindphitcmet;

            // dilepton hypothesis stuff
            for(unsigned hypi = 0; hypi < cms2.hyp_p4().size(); ++hypi)
            {
                // require 5/5 hypothesis
                if (min(cms2.hyp_lt_p4()[hypi].pt(), cms2.hyp_ll_p4()[hypi].pt()) < 5.)
                    continue;

                int index1 = cms2.hyp_lt_index()[hypi];
                int index2 = cms2.hyp_ll_index()[hypi];

                if (abs(cms2.hyp_lt_id()[hypi]) == 13 && !(cms2.mus_type()[index1] & 6))
                    continue;
                if (abs(cms2.hyp_ll_id()[hypi]) == 13 && !(cms2.mus_type()[index2] & 6))
                    continue;

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
                deltaphi_    = deltaPhi(phi1_, phi2_);
                ntrks_       = cms2.trks_trk_p4().size();

                // now, find jet closest to each hyp lepton
                float mindrjet1 = 999999.;
                float mindrjet2 = 999999.;
                for(unsigned int jeti = 0; jeti < theJets.size(); ++jeti)
                {
                    // for tight hyp lepton
                    float deta1 = cms2.hyp_lt_p4()[hypi].eta()-theJets[jeti].eta();
                    float dphi1 = deltaPhi(cms2.hyp_lt_p4()[hypi].phi(), theJets[jeti].phi());
                    float currdrjet1 = sqrt(deta1*deta1+dphi1*dphi1);
                    if (currdrjet1 < mindrjet1)
                        mindrjet1 = currdrjet1;

                    // for loose hyp lepton
                    float deta2 = cms2.hyp_ll_p4()[hypi].eta()-theJets[jeti].eta();
                    float dphi2 = deltaPhi(cms2.hyp_ll_p4()[hypi].phi(), theJets[jeti].phi());
                    float currdrjet2 = sqrt(deta2*deta2+dphi2*dphi2);
                    if (currdrjet2 < mindrjet2)
                        mindrjet2 = currdrjet2;
                }

                drjet1_ = mindrjet1;
                drjet2_ = mindrjet2;

                // figure out hypothesis type and fill appropriate iso, type variables
                if(hyp_type_ == 0)
                {
                    iso1_   = muonIsoValue(index1);
                    iso2_   = muonIsoValue(index2);
                    type1_  = cms2.mus_type()[index1];
                    type2_  = cms2.mus_type()[index2];

                    mu1_muonid_    = muonIdNotIsolated(index1, NominalTTbar); 
                    mu1_goodmask_  = cms2.mus_goodmask()[index1];
                    mu1_gfitchi2_  = cms2.mus_gfit_chi2()[index1] < -9000. ? -999999. : cms2.mus_gfit_chi2()[index1]/cms2.mus_gfit_ndof()[index1];

                    mu2_muonid_    = muonIdNotIsolated(index2, NominalTTbar);
                    mu2_goodmask_  = cms2.mus_goodmask()[index2];
                    mu2_gfitchi2_  = cms2.mus_gfit_chi2()[index2] < -9000. ? -999999. : cms2.mus_gfit_chi2()[index2]/cms2.mus_gfit_ndof()[index2];

                    int trkidx1 = cms2.mus_trkidx()[index1];
                    int trkidx2 = cms2.mus_trkidx()[index2];
                    d0vtx1_ = cms2.trks_d0vtx()[trkidx1];
                    d0vtx2_ = cms2.trks_d0vtx()[trkidx2];
                }
                else if(hyp_type_ == 1)
                {
                    iso1_   = muonIsoValue(index1);
                    iso2_   = electronIsolation_relsusy_cand1(index2, true);
                    type1_  = cms2.mus_type()[index1];
                    type2_  = cms2.els_type()[index2];

                    mu1_muonid_    = muonIdNotIsolated(index1, NominalTTbar); 
                    mu1_goodmask_  = cms2.mus_goodmask()[index1];
                    mu1_gfitchi2_  = cms2.mus_gfit_chi2()[index1] < -9000. ? -999999. : cms2.mus_gfit_chi2()[index1]/cms2.mus_gfit_ndof()[index1];

                    int trkidx1 = cms2.mus_trkidx()[index1];
                    int trkidx2 = cms2.els_trkidx()[index2];
                    d0vtx1_ = cms2.trks_d0vtx()[trkidx1];
                    if(trkidx2 >= 0)
                        d0vtx2_ = cms2.trks_d0vtx()[trkidx2];

                    e2_cand01_  = isGoodElectron(index2);
                    e2_eopin_   = cms2.els_eOverPIn()[index2];
                    e2_hoe_     = cms2.els_hOverE()[index2];
                    e2_dphiin_  = cms2.els_dPhiIn()[index2];
                    e2_detain_  = cms2.els_dEtaIn()[index2];
                    e2_eMe55_   = cms2.els_eMax()[index2] / cms2.els_e5x5()[index2];
                    e2_nmHits_  = cms2.els_exp_innerlayers()[index2];
                    e2_dcot_    = cms2.els_conv_dcot()[index2];
                    e2_dist_    = cms2.els_conv_dist()[index2];
                    e2_drmu_    = cms2.els_closestMuon()[index2] < 0 ? -999999. : cms2.els_musdr()[index2];
                }
                else if(hyp_type_ == 2)
                {
                    iso1_   = electronIsolation_relsusy_cand1(index1, true);
                    iso2_   = muonIsoValue(index2);
                    type1_  = cms2.els_type()[index1];
                    type2_  = cms2.mus_type()[index2];

                    e1_cand01_  = isGoodElectron(index1);
                    e1_eopin_   = cms2.els_eOverPIn()[index1];
                    e1_hoe_     = cms2.els_hOverE()[index1];
                    e1_dphiin_  = cms2.els_dPhiIn()[index1];
                    e1_detain_  = cms2.els_dEtaIn()[index1];
                    e1_eMe55_   = cms2.els_eMax()[index1] / cms2.els_e5x5()[index1];
                    e1_nmHits_  = cms2.els_exp_innerlayers()[index1];
                    e1_dcot_    = cms2.els_conv_dcot()[index1];
                    e1_dist_    = cms2.els_conv_dist()[index1];
                    e1_drmu_    = cms2.els_closestMuon()[index1] < 0 ? -999999. : cms2.els_musdr()[index1];

                    mu2_muonid_    = muonIdNotIsolated(index2, NominalTTbar); 
                    mu2_goodmask_  = cms2.mus_goodmask()[index2];
                    mu2_gfitchi2_  = cms2.mus_gfit_chi2()[index2] < -9000. ? -999999. : cms2.mus_gfit_chi2()[index2]/cms2.mus_gfit_ndof()[index2];

                    int trkidx1 = cms2.els_trkidx()[index1];
                    int trkidx2 = cms2.mus_trkidx()[index2];
                    if(trkidx1 >= 0)
                        d0vtx1_ = cms2.trks_d0vtx()[index1];
                    d0vtx2_ = cms2.trks_d0vtx()[trkidx2];
                }
                else if(hyp_type_ == 3)
                {
                    iso1_   = electronIsolation_relsusy_cand1(index1, true);
                    iso2_   = electronIsolation_relsusy_cand1(index2, true);
                    type1_  = cms2.els_type()[index1];
                    type2_  = cms2.els_type()[index2];

                    e1_cand01_  = isGoodElectron(index1);
                    e1_eopin_   = cms2.els_eOverPIn()[index1];
                    e1_hoe_     = cms2.els_hOverE()[index1];
                    e1_dphiin_  = cms2.els_dPhiIn()[index1];
                    e1_detain_  = cms2.els_dEtaIn()[index1];
                    e1_eMe55_   = cms2.els_eMax()[index1] / cms2.els_e5x5()[index1];
                    e1_nmHits_  = cms2.els_exp_innerlayers()[index1];
                    e1_dcot_    = cms2.els_conv_dcot()[index1];
                    e1_dist_    = cms2.els_conv_dist()[index1];
                    e1_drmu_    = cms2.els_closestMuon()[index1] < 0 ? -999999. : cms2.els_musdr()[index1];

                    e2_cand01_  = isGoodElectron(index2);
                    e2_eopin_   = cms2.els_eOverPIn()[index2];
                    e2_hoe_     = cms2.els_hOverE()[index2];
                    e2_dphiin_  = cms2.els_dPhiIn()[index2];
                    e2_detain_  = cms2.els_dEtaIn()[index2];
                    e2_eMe55_   = cms2.els_eMax()[index2] / cms2.els_e5x5()[index2];
                    e2_nmHits_  = cms2.els_exp_innerlayers()[index2];
                    e2_dcot_    = cms2.els_conv_dcot()[index2];
                    e2_dist_    = cms2.els_conv_dist()[index2];
                    e2_drmu_    = cms2.els_closestMuon()[index2] < 0 ? -999999. : cms2.els_musdr()[index2];

                    int trkidx1 = cms2.els_trkidx()[index1];
                    int trkidx2 = cms2.els_trkidx()[index2];
                    if(trkidx1 >= 0)
                        d0vtx1_ = cms2.trks_d0vtx()[trkidx1];
                    if(trkidx2 >= 0)
                        d0vtx2_ = cms2.trks_d0vtx()[trkidx2];
                }

                FillTwinNtuple();
            }
        }

        if (nEventsChain != nEventsTotal)
        {
            std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
        }
    }

    CloseTwinNtuple();
}

void twinmaker::InitTwinNtuple ()
{
    // event stuff
    run_          = -999999;
    ls_           = -999999;
    evt_          = -999999;
    njets_        = -999999;
    hyp_type_     = -999999;
    pfmet_        = -999999.;
    tcmet_        = -999999.;
    jet1pt_       = -999999.;
    dphipfmetjet_ = -999999.;
    dphitcmetjet_ = -999999.;
    dilpt_        = -999999.;
    deltaphi_     = -999999.;
    ntrks_        = -999999;

    // lepton stuff
    eormu1_        = -999999;
    eormu2_        = -999999;
    type1_         = -999999;
    type2_         = -999999;
    pt1_           = -999999.;
    pt2_           = -999999.;
    eta1_          = -999999.;
    eta2_          = -999999.;
    phi1_          = -999999.;
    phi2_          = -999999.;
    iso1_          = -999999.;
    iso2_          = -999999.;
    d0corr1_       = -999999.;
    d0corr2_       = -999999.;
    d0vtx1_        = -999999.;
    d0vtx2_        = -999999.;
    dphipfmet1_    = -999999.;
    dphipfmet2_    = -999999.;
    dphitcmet1_    = -999999.;
    dphitcmet2_    = -999999.;
    drjet1_        = -999999.;
    drjet2_        = -999999.;
    mass_          = -999999.;

    // muon stuff
    mu1_muonid_   = 0;
    mu2_muonid_   = 0;
    mu1_goodmask_ = -999999;
    mu2_goodmask_ = -999999;
    mu1_gfitchi2_ = -999999.;
    mu2_gfitchi2_ = -999999.;

    // electron stuff
    e1_cand01_   = 0;
    e2_cand01_   = 0;
    e1_eopin_    = -999999.;
    e2_eopin_    = -999999.;
    e1_hoe_      = -999999.;
    e2_hoe_      = -999999.;
    e1_dphiin_   = -999999.;
    e2_dphiin_   = -999999.;
    e1_detain_   = -999999.;
    e2_detain_   = -999999.;
    e1_eMe55_    = -999999.;
    e2_eMe55_    = -999999.;
    e1_nmHits_   = -999999;
    e2_nmHits_   = -999999;
    e1_dcot_     = -999999.;
    e2_dcot_     = -999999.;
    e1_dist_     = -999999.;
    e2_dist_     = -999999.;
    e1_drmu_     = -999999.;
    e2_drmu_     = -999999.;
}

void twinmaker::MakeTwinNtuple(const char *twinFilename)
{
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();

    twinFile_ = new TFile(Form("%s", twinFilename), "RECREATE");
    twinFile_->cd();
    twinTree_ = new TTree("tree", "A Twin Ntuple");

    // event stuff
    twinTree_->Branch("run",          &run_,         "run/I"         );
    twinTree_->Branch("ls",           &ls_,          "ls/I"          );
    twinTree_->Branch("evt",          &evt_,         "evt/I"         );
    twinTree_->Branch("njets",        &njets_,       "njets/I"       ); // uncorrected pt > 20
    twinTree_->Branch("hyp_type",     &hyp_type_,    "hyp_type/I"    );
    twinTree_->Branch("pfmet",        &pfmet_,       "pfmet/F"       );
    twinTree_->Branch("tcmet",        &tcmet_,       "tcmet/F"       );
    twinTree_->Branch("jet1pt",       &jet1pt_,      "jet1pt/F"      );
    twinTree_->Branch("dphipfmetjet", &dphipfmetjet_,"dphipfmetjet/F");
    twinTree_->Branch("dphitcmetjet", &dphitcmetjet_,"dphitcmetjet/F");
    twinTree_->Branch("dilpt",        &dilpt_,       "dilpt/F"       );
    twinTree_->Branch("deltaphi",     &deltaphi_,    "deltaphi/F"    );
    twinTree_->Branch("ntrks",        &ntrks_,       "ntrks/I"       );

    // lepton stuff
    twinTree_->Branch("eormu1",     &eormu1_,     "eormu1/I"    );
    twinTree_->Branch("eormu2",     &eormu2_,     "eormu2/I"    );
    twinTree_->Branch("type1",      &type1_,      "type1/I"     );
    twinTree_->Branch("type2",      &type2_,      "type2/I"     );
    twinTree_->Branch("pt1",        &pt1_,        "pt1/F"       );
    twinTree_->Branch("pt2",        &pt2_,        "pt2/F"       );
    twinTree_->Branch("eta1",       &eta1_,       "eta1/F"      );
    twinTree_->Branch("eta2",       &eta2_,       "eta2/F"      );
    twinTree_->Branch("phi1",       &phi1_,       "phi1/F"      );
    twinTree_->Branch("phi2",       &phi2_,       "phi2/F"      );
    twinTree_->Branch("iso1",       &iso1_,       "iso1/F"      );
    twinTree_->Branch("iso2",       &iso2_,       "iso2/F"      );
    twinTree_->Branch("d0corr1",    &d0corr1_,    "d0corr1/F"   );
    twinTree_->Branch("d0corr2",    &d0corr2_,    "d0corr2/F"   );
    twinTree_->Branch("d0vtx1",     &d0vtx1_,     "d0vtx1/F"    );
    twinTree_->Branch("d0vtx2",     &d0vtx2_,     "d0vtx2/F"    );
    twinTree_->Branch("dphipfmet1", &dphipfmet1_, "dphipfmet1/F");
    twinTree_->Branch("dphipfmet2", &dphipfmet2_, "dphipfmet2/F");
    twinTree_->Branch("dphitcmet1", &dphitcmet1_, "dphitcmet1/F");
    twinTree_->Branch("dphitcmet2", &dphitcmet2_, "dphitcmet2/F");
    twinTree_->Branch("drjet1",     &drjet1_,     "drjet1/F"    );
    twinTree_->Branch("drjet2",     &drjet2_,     "drjet2/F"    );
    twinTree_->Branch("mass",       &mass_,       "mass/F"      );

    // muon stuff
    twinTree_->Branch("mu1_muonid",   &mu1_muonid_,   "mu1_muonid/O"  );
    twinTree_->Branch("mu2_muonid",   &mu2_muonid_,   "mu2_muonid/O"  );
    twinTree_->Branch("mu1_goodmask", &mu1_goodmask_, "mu1_goodmask/I");
    twinTree_->Branch("mu2_goodmask", &mu2_goodmask_, "mu2_goodmask/I");
    twinTree_->Branch("mu1_gfitchi2", &mu1_gfitchi2_, "mu1_gfitchi2/F");
    twinTree_->Branch("mu2_gfitchi2", &mu2_gfitchi2_, "mu2_gfitchi2/F");

    // eectron stuff
    twinTree_->Branch("e1_cand01", &e1_cand01_, "e1_cand01/O");
    twinTree_->Branch("e2_cand01", &e2_cand01_, "e2_cand01/O");
    twinTree_->Branch("e1_eopin",  &e1_eopin_,  "e1_eopin/F" );
    twinTree_->Branch("e2_eopin",  &e2_eopin_,  "e2_eopin/F" );
    twinTree_->Branch("e1_hoe",    &e1_hoe_,    "e1_hoe/F"   );
    twinTree_->Branch("e2_hoe",    &e2_hoe_,    "e2_hoe/F"   );
    twinTree_->Branch("e1_dphiin", &e1_dphiin_, "e1_dphiin/F");
    twinTree_->Branch("e2_dphiin", &e2_dphiin_, "e2_dphiin/F");
    twinTree_->Branch("e1_detain", &e1_detain_, "e1_detain/F");
    twinTree_->Branch("e2_detain", &e2_detain_, "e2_detain/F");
    twinTree_->Branch("e1_eMe55",  &e1_eMe55_,  "e1_eMe55/F" ); // for spikes
    twinTree_->Branch("e2_eMe55",  &e2_eMe55_,  "e2_eMe55/F" ); // for spikes
    twinTree_->Branch("e1_nmHits", &e1_nmHits_, "e1_nmHits/I");
    twinTree_->Branch("e2_nmHits", &e2_nmHits_, "e2_nmHits/I");
    twinTree_->Branch("e1_dcot",   &e1_dcot_,   "e1_dcot/F"  );
    twinTree_->Branch("e2_dcot",   &e2_dcot_,   "e2_dcot/F"  );
    twinTree_->Branch("e1_dist",   &e1_dist_,   "e1_dist/F"  );
    twinTree_->Branch("e2_dist",   &e2_dist_,   "e2_dist/F"  );
    twinTree_->Branch("e1_drmu",   &e1_drmu_,   "e1_drmu/F"  );
    twinTree_->Branch("e2_drmu",   &e2_drmu_,   "e2_drmu/F"  );
}

void twinmaker::FillTwinNtuple()
{
    twinTree_->Fill();
}

void twinmaker::CloseTwinNtuple()
{
    twinFile_->cd();
    twinTree_->Write();
    twinFile_->Close();
}

float deltaPhi (float phi1, float phi2)
{
    float dphi = phi1-phi2;
    if (dphi < 0.) dphi = -1.*dphi;
    if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
    return dphi;
}

bool sortByPt (const LorentzVector &vec1, const LorentzVector &vec2)
{
    return vec1.pt() > vec2.pt();
}

bool isGoodElectron(const int index)
{
    if (! electronId_noMuon(index))
        return false;
    if (! electronId_cand01(index))
        return false;
    //if (! electronImpact_cand01(index))
    //    return false;
    //if (isFromConversionPartnerTrack(index))
    //    return false;
    // Note that this is not currently
    // in electronSelection_cand01
    //if (isFromConversionHitPattern(index))
    //    return false;

    return true;
}
