
//
// Dave "the one but not the only" Evans 
//

#include "MyScanChain.h"

// ROOT includes
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TDirectory.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

// CMS2 includes
#include "../../CORE/CMS2.h"

#include "../../CORE/eventSelections.h"
#include "../../CORE/electronSelections.h"
#include "../../CORE/muonSelections.h"
#include "../../CORE/jetSelections.h"
#include "../../CORE/trackSelections.h"
#include "../../CORE/ttbarSelections.h"

#include "../../Tools/DileptonHypType.h"
#include "../../Tools/tools.cc"

#include "../../Tools/goodrun.cc"
#include "../../CORE/utilities.h"
#include "../../CORE/triggerUtils.h"
#include "../../CORE/mcSelections.h"
#include "../../Tools/EgammaAnalysisTools/include/LikelihoodUtil.h"

#include <set>
#include <algorithm>

//
// Namespaces
//
using namespace tas;

//
// functions
//

void MyScanChain::InitBaby()
{
    run_ = cms2.evt_run();
    ls_  = cms2.evt_lumiBlock();
    evt_ = cms2.evt_event();
    type_ = -999;
    weight_ = -999.;

    reco_ptlt_  = -999.;
    reco_etalt_ = -999.;
    reco_philt_ = -999.;
    reco_isolt_ = -999.;
    reco_typelt_ = -999;
    reco_mctypelt_ = -999;
    reco_algolt_ = -999;
    reco_vbtf95lt_ = -999;
    reco_vbtf90lt_ = -999;
    reco_vbtf80lt_ = -999;
    reco_vbtf70lt_ = -999;

    reco_ptll_  = -999.;
    reco_etall_ = -999.;
    reco_phill_ = -999.;
    reco_isoll_ = -999.;
    reco_typell_ = -999;
    reco_mctypell_ = -999;
    reco_algoll_ = -999;
    reco_vbtf95ll_ = -999;
    reco_vbtf90ll_ = -999;
    reco_vbtf80ll_ = -999;
    reco_vbtf70ll_ = -999;

    reco_mdil_ = -999.;
    reco_vdilpt_ = -999.;
    reco_vdilphi_ = -999.;

    // met and jets
    reco_tcmet_ = -999.;
    reco_pfmet_ = -999.;
    reco_npfjets_ = -999;

    // which event selections 
    // passed
    reco_lhlt_ = -999.9;
    reco_lhll_ = -999.9;
    reco_pfmvalt_ = -999.9;
    reco_pfmvall_ = -999.9;

}

//
// Dilepton analysis
//

void MyScanChain::AnalyseDilepton(const float &weight)
{

    //
    // loop on hypothesis
    // and pick one if multiple
    //

    // ttbar pass5
    // do not apply aplignment correction in endcap
    // do not remove detain cuts in endcap
    bool applyAlignmentCorrection = false;
    bool removedEtaCutInEndcap = false;

    // list of selected genericsearch hyps
    std::vector<unsigned int> selectedHyps;

    for (size_t h = 0; h < cms2.hyp_type().size(); ++h) {

        //
        // select good hypothesis
        //

        //
        // pt
        //

        if (!(cms2.hyp_ll_p4()[h].Pt() > 10.0 && cms2.hyp_lt_p4()[h].Pt() > 10.0)) continue;

        //
        // hyps from same vertex
        //

        if (!hypsFromSameVtx(h)) continue;

        //
        // leptons must have opposite charge
        //

        if (cms2.hyp_lt_id()[h] * cms2.hyp_ll_id()[h] > 0) continue;

        //
        // check which triggers this hyp passed
        //

        // record the event decision
        // for the generic search trigger path
        //pass_trig_gs_ = 1;

        /*
        // single object triggers
        pass_trig_single_e_ = passEGTrigger(!isData_);
        pass_trig_single_mu_ = passMuTrigger(!isData_);

        // pass5 requires at least one muon with eta < 2.1
        if (cms2.hyp_type()[h] == 0) {
        if (!(fabs(cms2.hyp_lt_p4()[h].Eta()) > 2.1 && fabs(cms2.hyp_ll_p4()[h].Eta()) > 2.1)
        && pass_trig_single_mu_) pass_trig_gs_ = 1;
        }

        // pass5 requires for emu fiducial cuts on the e and the m
        // em
        if ((cms2.hyp_type()[h] == 1 
        && fabs(cms2.hyp_lt_p4()[h].Eta()) < 2.1 
        && cms2.els_eSC()[cms2.hyp_ll_index()[h]]/cosh(cms2.els_etaSC()[cms2.hyp_ll_index()[h]]) > 17.0
        && pass_trig_single_mu_ || pass_trig_single_e_)) pass_trig_gs_ = 1;
        // me
        if ((cms2.hyp_type()[h] == 2
        && fabs(cms2.hyp_ll_p4()[h].Eta()) < 2.1 
        && cms2.els_eSC()[cms2.hyp_lt_index()[h]]/cosh(cms2.els_etaSC()[cms2.hyp_lt_index()[h]]) > 17.0
        && pass_trig_single_mu_ || pass_trig_single_e_)) pass_trig_gs_ = 1;

        // electrons
        if (cms2.hyp_type()[h] == 3 && pass_trig_single_e_) pass_trig_gs_ = 1; 

        //
        // begin by checking generic search selection
        //
         */

        //if(pass_trig_gs_) {
        selectedHyps.push_back(h);
        //} 

    } // end loop on hyps

    if (selectedHyps.size() == 0) return;

    //
    // now examine the hypotheses selected for this event
    //

    unsigned int h = 0;

    // case 1:
    // - there is exactly one genericSearch dilepton hypothesis
    // -- do not need to do anything
    if (selectedHyps.size() == 1) h = selectedHyps[0];

    //
    // case 2:
    // - there are two or more genericSearch dilepton hypotheses
    // -- if there exactly one same flavour, pick it
    // -- if there are multiple same flavour, pick closest to Z mass

    if (selectedHyps.size() > 1) {
        float closestToZMass = 9999999.99;
        unsigned int closestToZ = 0;
        std::vector<unsigned int> selectedSameFlavorHyps;
        for (size_t i = 0; i < selectedHyps.size(); ++i) {
            unsigned int idx = selectedHyps[i];
            if (cms2.hyp_type()[idx] == 0 || cms2.hyp_type()[idx] == 3) {
                selectedSameFlavorHyps.push_back(idx);
                float mass = sqrt(cms2.hyp_p4()[selectedHyps[i]].M2());
                if (fabs(mass - 91.0) < closestToZMass) {
                    closestToZMass = mass;
                    closestToZ = idx;
                }
            }
        } //end of loop over selectedHyps

        // set the chosen hypothesis to cloestToZ
        // which is the index in the hyps block of the same flavor 
        // hyp with the mass closest to the Z.  This will be the one 
        // which is closest to the Z out of multiple choices if there are
        // multiple choices, otherwise it is just the only choice
        //Note: if you have 2 or more hyps then you have 3 or more leptons.
        //      if you have 3 or more leptons then you always have at least one same flavor hypothesis.
        //      In an eemu event, we pick the ee as our preferred hypothesis, irrespective of charge.
        //      Once we picked a hypothesis, we record the charge of that hypothesis a few lines below.
        //      This means that sometimes, the only indication that you had a ss hypothesis in the event is that
        //      tri = true, as ss was not set given that we picked the os dilepton hypothesis.
        h = closestToZ;

    }

    //
    // now compute what cuts passed and 
    // other quantities for plotting
    //

    //
    // record charge flip
    //

    //
    // hyp type
    //

    DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[h]);

    //
    // PF jets
    //
    reco_npfjets_ = 0;//nJets(h, JETS_TYPE_PF_CORR, JETS_CLEAN_HYP_E_MU, 0.4, 30.0, 2.5);
    //
    // invariant mass
    //

    LorentzVector vec_dilepton = cms2.hyp_lt_p4()[h] + cms2.hyp_ll_p4()[h];
    reco_mdil_ = sqrt(vec_dilepton.M2());
    reco_vdilpt_ = vec_dilepton.Pt();
    reco_vdilphi_ = vec_dilepton.Phi();

    //
    // set the baby variables for any of the three analyses
    // which have not yet been set
    //

    electronIdComponent_t answer_vbtf = 0;

    // kinematics
    reco_ptlt_ = cms2.hyp_lt_p4()[h].Pt();
    reco_etalt_ = cms2.hyp_lt_p4()[h].Eta();
    reco_philt_ = cms2.hyp_lt_p4()[h].Phi();
    reco_typelt_ = abs(cms2.hyp_lt_id()[h]);
    reco_mctypelt_ = abs(cms2.hyp_lt_mc_id()[h]);
    if (reco_typelt_ == 11) {
        reco_algolt_ = cms2.els_type()[cms2.hyp_lt_index()[h]];
        reco_pfmvalt_ = cms2.els_mva()[cms2.hyp_lt_index()[h]];
        reco_lhlt_ = likelihoodUtil_->getValue(cms2.hyp_lt_index()[h]);
        reco_isolt_ = electronIsolation_rel(cms2.hyp_lt_index()[h], true);

        cuts_t cuts_passed = electronSelection(cms2.hyp_lt_index()[h], applyAlignmentCorrection, removedEtaCutInEndcap);
        reco_vbtf95lt_ = pass_electronSelectionCompareMask(cuts_passed, (1<<ELEID_VBTF_35X_95));
        reco_vbtf90lt_ = pass_electronSelectionCompareMask(cuts_passed, (1<<ELEID_VBTF_35X_90));
        reco_vbtf80lt_ = pass_electronSelectionCompareMask(cuts_passed, (1<<ELEID_VBTF_35X_80));
        reco_vbtf70lt_ = pass_electronSelectionCompareMask(cuts_passed, (1<<ELEID_VBTF_35X_70));

    } if (reco_typelt_ == 13) {
        reco_isolt_ = muonIsoValue(cms2.hyp_lt_index()[h]);
    }


    reco_ptll_ = cms2.hyp_ll_p4()[h].Pt();
    reco_etall_ = cms2.hyp_ll_p4()[h].Eta();
    reco_phill_ = cms2.hyp_ll_p4()[h].Phi();
    reco_typell_ = abs(cms2.hyp_ll_id()[h]);
    reco_mctypell_ = abs(cms2.hyp_ll_mc_id()[h]);
    if (reco_typell_ == 11) {
        reco_algoll_ = cms2.els_type()[cms2.hyp_ll_index()[h]];
        reco_pfmvall_ = cms2.els_mva()[cms2.hyp_ll_index()[h]];
        reco_lhll_ = likelihoodUtil_->getValue(cms2.hyp_ll_index()[h]);
        reco_isoll_ = electronIsolation_rel(cms2.hyp_ll_index()[h], true);

        cuts_t cuts_passed = electronSelection(cms2.hyp_ll_index()[h], applyAlignmentCorrection, removedEtaCutInEndcap);
        reco_vbtf95ll_ = pass_electronSelectionCompareMask(cuts_passed, (1<<ELEID_VBTF_35X_95));
        reco_vbtf90ll_ = pass_electronSelectionCompareMask(cuts_passed, (1<<ELEID_VBTF_35X_90));
        reco_vbtf80ll_ = pass_electronSelectionCompareMask(cuts_passed, (1<<ELEID_VBTF_35X_80));
        reco_vbtf70ll_ = pass_electronSelectionCompareMask(cuts_passed, (1<<ELEID_VBTF_35X_70));

    } if (reco_typell_ == 13) {
        reco_isoll_ = muonIsoValue(cms2.hyp_ll_index()[h]);
    }

    // event properties
    reco_tcmet_ = cms2.evt_tcmet();
    reco_pfmet_ = cms2.evt_pfmet();
    type_ = hypType;

    //
    // fill the baby tree
    //

    weight_ = weight;
    babyTree_->Fill();

    return;

}

//
// set good run list
//

void MyScanChain::setGoodRunList(std::string fname) 
{
    set_goodrun_file(fname.c_str());
}

//
// Main function
//

int MyScanChain::ScanChain(bool isData, std::string sampleName, TChain *chain, float kFactor, int nEvents, std::string skimFilePrefix) {

    std::cout << "scanning " << sampleName << std::endl;

    //
    // Tree related
    // 

    babyFile_ = new TFile(Form("baby_%s.root", sampleName.c_str()), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");

    babyTree_->Branch("run",                &run_,          "run/I"         );
    babyTree_->Branch("ls",                 &ls_,           "ls/I"          );
    babyTree_->Branch("evt",                &evt_,          "evt/I"         );
    babyTree_->Branch("type",               &type_,         "type/I"        );
    babyTree_->Branch("weight",             &weight_,       "weight/F"         );

    babyTree_->Branch("reco_ptlt",          &reco_ptlt_,    "reco_ptlt/F"       );
    babyTree_->Branch("reco_etalt",         &reco_etalt_,   "reco_etalt/F"      );
    babyTree_->Branch("reco_philt",         &reco_philt_,   "reco_philt/F"      );
    babyTree_->Branch("reco_isolt",         &reco_isolt_,   "reco_isolt/F"      );
    babyTree_->Branch("reco_typelt",          &reco_typelt_,    "reco_typelt/I"       );
    babyTree_->Branch("reco_mctypelt",          &reco_mctypelt_,    "reco_mctypelt/I"       );
    babyTree_->Branch("reco_algolt",          &reco_algolt_,    "reco_algolt/I"       );
    babyTree_->Branch("reco_vbtf95lt",          &reco_vbtf95lt_,    "reco_vbtf95lt/I"       );
    babyTree_->Branch("reco_vbtf90lt",          &reco_vbtf90lt_,    "reco_vbtf90lt/I"       );
    babyTree_->Branch("reco_vbtf80lt",          &reco_vbtf80lt_,    "reco_vbtf80lt/I"       );
    babyTree_->Branch("reco_vbtf70lt",          &reco_vbtf70lt_,    "reco_vbtf70lt/I"       );

    babyTree_->Branch("reco_ptll",          &reco_ptll_,            "reco_ptll/F"          );
    babyTree_->Branch("reco_etall",         &reco_etall_,           "reco_etall/F"         );
    babyTree_->Branch("reco_phill",         &reco_phill_,           "reco_phill/F"         );
    babyTree_->Branch("reco_isoll",         &reco_isoll_,   "reco_isoll/F"      );    
    babyTree_->Branch("reco_typell",          &reco_typell_,    "reco_typell/I"       );
    babyTree_->Branch("reco_mctypell",          &reco_mctypell_,    "reco_mctypell/I"       );
    babyTree_->Branch("reco_algoll",          &reco_algoll_,    "reco_algoll/I"       );
    babyTree_->Branch("reco_vbtf95ll",          &reco_vbtf95ll_,    "reco_vbtf95ll/I"       );
    babyTree_->Branch("reco_vbtf90ll",          &reco_vbtf90ll_,    "reco_vbtf90ll/I"       );
    babyTree_->Branch("reco_vbtf80ll",          &reco_vbtf80ll_,    "reco_vbtf80ll/I"       );
    babyTree_->Branch("reco_vbtf70ll",          &reco_vbtf70ll_,    "reco_vbtf70ll/I"       );

    babyTree_->Branch("reco_tcmet",         &reco_tcmet_,           "reco_tcmet/F"          );
    babyTree_->Branch("reco_pfmet",         &reco_pfmet_,           "reco_pfmet/F"          );

    babyTree_->Branch("reco_mdil",          &reco_mdil_,            "reco_mdil/F"           );
    babyTree_->Branch("reco_vdilpt",        &reco_vdilpt_,          "reco_vdilpt/F"         );
    babyTree_->Branch("reco_vdilphi",        &reco_vdilphi_,          "reco_vdilphi/F"         );
    babyTree_->Branch("reco_npfjets",       &reco_npfjets_,         "reco_npfjets/I"        );

    babyTree_->Branch("reco_lhlt",          &reco_lhlt_,            "reco_lhlt/F"           );
    babyTree_->Branch("reco_lhll",          &reco_lhll_,            "reco_lhll/F"           );
    babyTree_->Branch("reco_pfmvalt",          &reco_pfmvalt_,            "reco_pfmvalt/F"           );
    babyTree_->Branch("reco_pfmvall",          &reco_pfmvall_,            "reco_pfmvall/F"           );

    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    if (rootdir == 0){
        std::cout<<"Head directory root: not found. Try Rint: ..."<<std::endl;
        rootdir = gROOT->GetDirectory("Rint:");
        if (rootdir){
            std::cout<<"OK: Got Rint:"<<std::endl;
        } else {
            std::cout<<"ERROR: no root: or Rint: found. Histograms will likely be lost"<<std::endl;
        }
    } 
    rootdir->cd(); 

    //
    // set data flag
    //

    isData_ = isData;

    //
    // set up likelihood id utility
    //

    likelihoodUtil_ = new LikelihoodUtil("../../Tools/EgammaAnalysisTools/PDFs/pdfs_MC.root");

    //
    // file loop
    //

    unsigned int nEventsChain=0;
    if(nEvents == -1) nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;
    int i_permille_old = 0;

    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
        TFile *f = TFile::Open(currentFile->GetTitle());
        TTree *tree = (TTree*)f->Get("Events");
        cms2.Init(tree);

        //Event Loop
        ULong64_t nEvents = tree->GetEntries();
        for(ULong64_t event = 0; event < nEvents; ++event) {
            cms2.GetEntry(event);
            ++nEventsTotal;

            // Progress feedback to the user
            int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
            if (i_permille != i_permille_old) {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1)) {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                            "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
                    fflush(stdout);
                }
                i_permille_old = i_permille;
            }

            // work out event weight
            float weight = 1.0;
            if (!isData) {
                weight = cms2.evt_scale1fb() * kFactor;
            }

            //
            // get the MC source correct
            //

            if (!isData) {
                int nels, nmus, ntaus;
                int nlep = leptonGenpCount_lepTauDecays(nels, nmus, ntaus);
                if (sampleName == "dyee"     &&  nels != 2) continue;
                if (sampleName == "dymm"     &&  nmus != 2) continue;
                if (sampleName == "dytt"     &&  ntaus != 2) continue;
            }


            // if we are looking at a drell yan sample then it is made up from diffeent mass regions
            // - these need to be stitched together to avoid double counting
            // - additionally, the individual pieces have different k factors which must be set

            if ((sampleName == "dyee" || sampleName == "dymm" || sampleName == "dytt")) {

                // if we are not looking at the madgraph then we must have one of the pythia pieces
                // - if these have an event with mass above the start of the madgraph sample then skip

                if (TString(cms2.evt_dataset()).Contains("madgraph") == false) {
                    bool skipEvent = false;
                    for (unsigned int i = 0; i < genps_p4().size(); i++){
                        if(abs(cms2.genps_id()[i]) == 23 && cms2.genps_p4()[i].M() > 50.) {
                            skipEvent = true;
                            break;
                        }
                    }

                    // skip the event if necessary
                    if (skipEvent) continue;

                    // if the event was not skipped then set the k factor
                    // for the appropriate sample

                    if (TString(evt_dataset()).Contains("M10to20") == true) { //10 < mll < 20 
                        weight = weight*3457./2659.;
                    } else { // 20 < mll < 50
                        weight = weight * 1666./1300.;
                    }

                    // otherwise we are looking at the madgraph sample
                    // this has a mass > 50.0 and we always include it

                } else {
                    // set the weight for madgraph
                    weight = weight*3048./2400.;
                }

            }

            //dumpDocLines();

            //
            // do good run check
            //

            if (isData) {
                if (!goodrun(cms2.evt_run(), cms2.evt_lumiBlock())) continue;
            }

            //
            // do duplicate check
            //

            if (isData) {
                DorkyEventIdentifier id(cms2.evt_run(), cms2.evt_event(), cms2.evt_lumiBlock());
                if (is_duplicate(id)) {
                    continue;
                }
            }

            //
            // init the baby variables
            //

            InitBaby();
            weight_ = weight;

            //
            // standard event cleaning
            //

            if (!cleaning_standardAugust2010(isData)) continue;

            //
            // Do dilepton study
            //

            AnalyseDilepton(weight);


        } // end loop on files

        delete f;

    } // end loop on events

    if ( nEventsChain != nEventsTotal ) {
        std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    }

    //
    // close the baby file
    //

    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();

    //
    // delete the likelihood util
    //
    std::cout << "deleting lh" << std::endl;
    delete likelihoodUtil_;

    //
    // make sure we're back in the right root dir
    //

    rootdir = gROOT->GetDirectory("root:");
    if (rootdir) rootdir->cd();
    else{
        std::cout<<"Cant find root: . Current dir is "<<gDirectory->GetName()<<std::endl;
        rootdir = gROOT->GetDirectory("Rint:");
        if (rootdir){
            std::cout<<"OK, got Rint: "<<std::endl;
            rootdir->cd();
        } else {
            std::cout<<"Cant find Rint: either . Current dir is "<<gDirectory->GetName()<<std::endl;
        }
    }

    return 0;
}

