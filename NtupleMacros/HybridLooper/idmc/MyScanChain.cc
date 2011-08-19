
//
// ttbar -> ll
// Dave "the one but not the only" Evans 
//

#include "MyScanChain.h"

// ROOT includes
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"

#include "Math/LorentzVector.h"

// CMS2 includes
#include "../../CORE/CMS2.h"
#include "../../CORE/muonSelections.h"
#include "../../CORE/jetSelections.h"
#include "../../CORE/mcSelections.h"
#include "../../CORE/utilities.h"

#include "../../CORE/electronSelectionsParameters.h"
#include "../../CORE/electronSelections.h"

//
// Namespaces
//
using namespace tas;

//
//
//

enum ElectronSelection {
    ELEPASS_PT10,
    ELEPASS_PT20,
    ELEPASS_PT10NOT20,
};


static const char det_names[][128] = { "EB", "EE"};

//
// for hyps
//

enum DileptonHypType MyScanChain::hyp_typeToHypType (int hyp_type)
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

void MyScanChain::Fill(TH1F** hist, const unsigned int hyp, const float &val, const float &weight)
{
    hist[hyp]->Fill(val, weight);
    hist[DILEPTON_ALL]->Fill(val, weight);
}

void MyScanChain::Fill2D(TH2F** hist, const unsigned int hyp, const float &valx, const float &valy, const float &weight)
{   
    hist[hyp]->Fill(valx, valy, weight); 
    hist[DILEPTON_ALL]->Fill(valx, valy, weight);
}


void MyScanChain::FormatHist(TH1F** hist, std::string sampleName, std::string name, int n, float min, float max)
{       
    // loop on EB, EE
    for (unsigned int i = 0; i < 4; ++i)
    {
        std::string str = dilepton_hypo_names[i];
        std::string title = name + "_" + str;
        hist[i] = new TH1F(Form("%s_%s_%s", sampleName.c_str(), name.c_str(), str.c_str()),
                title.c_str(), n, min, max);
        hist[i]->GetXaxis()->SetTitle(name.c_str());
        hist[i]->Sumw2();
    }
}    

void MyScanChain::FormatHist2D(TH2F** hist, std::string sampleName, std::string name, int nx, float minx, float maxx, int ny, float miny, float maxy)
{       
    // loop on EB, EE 
    for (unsigned int i = 0; i < 4; ++i)
    {   
        std::string str = dilepton_hypo_names[i];
        std::string title = name + "_" + str;
        hist[i] = new TH2F(Form("%s_%s_%s", sampleName.c_str(), name.c_str(), str.c_str()),
                title.c_str(), nx, minx, maxx, ny, miny, maxy);
        hist[i]->GetXaxis()->SetTitle(name.c_str());
        hist[i]->Sumw2();
    }
}

bool MyScanChain::CheckCutsNM1(cuts_t apply, cuts_t remove, cuts_t passed)
{
    if ((passed & (apply & (~remove))) == (apply & (~remove))) return true;
    return false;
}

bool MyScanChain::CheckCuts(cuts_t apply, cuts_t passed)
{
    if ((apply & passed) == apply) return true;
    return false;
}


void MyScanChain::FormatAllEleIdHistograms(std::string sampleName)
{

    for (unsigned int i = 0; i < 2; ++i) {
        std::string detname = det_names[i];

        FormatHist(h1_hyp_pt_[i], sampleName, "h1_hyp_pt_" + detname, 200, 0.0, 200.0);              
        FormatHist(h1_hyp_reliso_[i], sampleName, "h1_hyp_reliso_" + detname, 100, 0, 1.0);
        FormatHist(h1_hyp_pdgid_[i], sampleName, "h1_hyp_pdgid_" + detname, 10, -0.5, 9.5);

        FormatHist(h1_hyp_cand01_pt_[i], sampleName, "h1_hyp_cand01_pt_" + detname, 200, 0.0, 200.0);
        FormatHist(h1_hyp_cand01_pdgid_[i], sampleName, "h1_hyp_cand01_pdgid_" + detname, 10, -0.5, 9.5);
        FormatHist(h1_hyp_cand01_reliso_[i], sampleName, "h1_hyp_cand01_reliso_" + detname, 100, 0.0, 1.0);

        FormatHist(h1_hyp_distdcot002_pt_[i], sampleName, "h1_hyp_distdcot002_pt_" + detname, 200, 0.0, 200.0);
        FormatHist(h1_hyp_distdcot002_pdgid_[i], sampleName, "h1_hyp_distdcot002_pdgid_" + detname, 10, -0.5, 9.5);

        FormatHist(h1_hyp_hitpattern_pt_[i], sampleName, "h1_hyp_hitpattern_pt_" + detname, 200, 0.0, 200.0);
        FormatHist(h1_hyp_hitpattern_pdgid_[i], sampleName, "h1_hyp_hitpattern_pdgid_" + detname, 10, -0.5, 9.5);

        FormatHist(h1_hyp_convboth_pt_[i], sampleName, "h1_hyp_convboth_pt_" + detname, 200, 0.0, 200.0);
        FormatHist(h1_hyp_convboth_pdgid_[i], sampleName, "h1_hyp_convboth_pdgid_" + detname, 10, -0.5, 9.5);

        // check reimplementation of sani id in the looper
        FormatHist(h1_hyp_idstudy_classExpLooseRecompId_[i], sampleName, "h1_hyp_idstudy_classExpLooseRecompId_" + detname, 5, -0.5, 4.5);
        FormatHist(h1_hyp_idstudy_classExpLooseRecompIso_[i], sampleName, "h1_hyp_idstudy_classExpLooseRecompIso_" + detname, 5, -0.5, 4.5);
        FormatHist(h1_hyp_idstudy_classExpTightRecompId_[i], sampleName, "h1_hyp_idstudy_classExpTightRecompId_" + detname, 5, -0.5, 4.5);
        FormatHist(h1_hyp_idstudy_classExpTightRecompIso_[i], sampleName, "h1_hyp_idstudy_classExpTightRecompIso_" + detname, 5, -0.5, 4.5);

        FormatHist(h1_hyp_idstudy_classExpSaniLooseId_[i], sampleName, "h1_hyp_idstudy_classExpSaniLooseId_" + detname, 5, -0.5, 4.5);
        FormatHist(h1_hyp_idstudy_classExpSaniLooseIso_[i], sampleName, "h1_hyp_idstudy_classExpSaniLooseIso_" + detname, 5, -0.5, 4.5);
        FormatHist(h1_hyp_idstudy_classExpSaniTightId_[i], sampleName, "h1_hyp_idstudy_classExpSaniTightId_" + detname, 5, -0.5, 4.5);
        FormatHist(h1_hyp_idstudy_classExpSaniTightIso_[i], sampleName, "h1_hyp_idstudy_classExpSaniTightIso_" + detname, 5, -0.5, 4.5);


    }

}

void MyScanChain::FillAllEleIdHistograms(const unsigned int index, const float &weight, const TString &sampleName, const unsigned int hyp)
{

    // apply truth match behavior if ttbar
    if (sampleName == "ttbar") {
        if(!((abs(cms2.els_mc_id()[index]) == 11) && abs(cms2.els_mc_motherid()[index]) == 24) ) return;
    }
    // apply truth match behavior if wjets
    if (sampleName == "wjets") {
        if(((abs(cms2.els_mc_id()[index]) == 11) && abs(cms2.els_mc_motherid()[index]) == 24) ) return;
    }

    //
    // work out which cuts passed
    //

    DileptonHypType hypType = DILEPTON_EE;


    // find detector 
    unsigned int det = 0;
    if (fabs(cms2.els_etaSC()[index]) > 1.479) det = 1;

    //
    // controllable cuts
    //

    cuts_t general_cuts_passed = 0;
    if (cms2.els_p4()[index].Pt() > 10.0) general_cuts_passed |= (1<<ELEPASS_PT10);
    if (cms2.els_p4()[index].Pt() > 20.0) general_cuts_passed |= (1<<ELEPASS_PT20);
    if (cms2.els_p4()[index].Pt() > 10.0 && cms2.els_p4()[index].Pt() < 20.0) general_cuts_passed |= (1<<ELEPASS_PT10NOT20);
    if (!((general_cuts_passed & configured_cuts_) == configured_cuts_)) return;

    //
    // compute common variables
    //

    float E2x5MaxOver5x5 = cms2.els_e2x5Max()[index] / cms2.els_e5x5()[index];
    float iso_rel = electronIsolation_rel(index, true);
    int pdgidCatagory = 4;
    int mcId = abs(cms2.els_mc_id()[index]);
    int motherId = abs(cms2.els_mc_motherid()[index]);
    if ((mcId == 11) && (motherId == 24 || motherId == 23)) pdgidCatagory = 0;
    else pdgidCatagory = elFakeMCCategory(index);

    //
    // get the standard cut bits for this electron
    //

    cuts_t electron_cuts_passed = electronSelection(index);

    //
    // basic denominator must pass
    //
    static const cuts_t denominator = (1ll<<ELEIP_200) | (1ll<<ELEETA_250) | (1ll<<ELENOMUON_010) | (1ll<<ELESEED_ECAL);

    if (CheckCuts(denominator, electron_cuts_passed)) {

        //
        // fill histograms
        //

        Fill(h1_hyp_pt_[det], hypType, cms2.els_p4()[index].Pt(), weight);
        Fill(h1_hyp_reliso_[det], hypType, iso_rel, weight);    
        Fill(h1_hyp_pdgid_[det], hypType, pdgidCatagory, weight);

        if (CheckCuts((1ll<<ELEID_CAND01), electron_cuts_passed)) {
            Fill(h1_hyp_cand01_pt_[det], hypType, cms2.els_p4()[index].Pt(), weight);
            Fill(h1_hyp_cand01_reliso_[det], hypType, iso_rel, weight);
            Fill(h1_hyp_cand01_pdgid_[det], hypType, pdgidCatagory, weight);
        }

        if (CheckCuts((1ll<<ELENOTCONV_DISTDCOT002), electron_cuts_passed)) {
            Fill(h1_hyp_distdcot002_pt_[det], hypType, cms2.els_p4()[index].Pt(), weight);
            Fill(h1_hyp_distdcot002_pdgid_[det], hypType, pdgidCatagory, weight);
        }

        if (CheckCuts((1ll<<ELENOTCONV_HITPATTERN), electron_cuts_passed)) {
            Fill(h1_hyp_hitpattern_pt_[det], hypType, cms2.els_p4()[index].Pt(), weight);
            Fill(h1_hyp_hitpattern_pdgid_[det], hypType, pdgidCatagory, weight);
        }

        if (CheckCuts(((1ll<<ELENOTCONV_DISTDCOT002) | (1ll<<ELENOTCONV_HITPATTERN)), electron_cuts_passed)) {
            Fill(h1_hyp_convboth_pt_[det], hypType, cms2.els_p4()[index].Pt(), weight);
            Fill(h1_hyp_convboth_pdgid_[det], hypType, pdgidCatagory, weight);
        }

        //
        // validate sani ID
        // 

        int answerLoose = electronId_CIC(index, 4, CIC_LOOSE);
        int answerTight = electronId_CIC(index, 4, CIC_TIGHT);

        // check that class based id passed independently of isolation decision
        Fill(h1_hyp_idstudy_classExpLooseRecompId_[det], hypType, bool(answerLoose & (1<<ELEID_ID)), 1.0);
        Fill(h1_hyp_idstudy_classExpLooseRecompIso_[det], hypType, bool(answerLoose & (1<<ELEID_ISO)), 1.0);
        Fill(h1_hyp_idstudy_classExpTightRecompId_[det], hypType, bool(answerTight & (1<<ELEID_ID)), 1.0);
        Fill(h1_hyp_idstudy_classExpTightRecompIso_[det], hypType, bool(answerTight & (1<<ELEID_ISO)), 1.0);

        Fill(h1_hyp_idstudy_classExpSaniLooseId_[det], hypType, bool(int(cms2.els_egamma_looseId()[index]) & (1<<0)), 1.0);
        Fill(h1_hyp_idstudy_classExpSaniLooseIso_[det], hypType, bool(int(cms2.els_egamma_looseId()[index]) & (1<<1)), 1.0);
        Fill(h1_hyp_idstudy_classExpSaniTightId_[det], hypType, bool(int(cms2.els_egamma_tightId()[index]) & (1<<0)), 1.0);
        Fill(h1_hyp_idstudy_classExpSaniTightIso_[det], hypType, bool(int(cms2.els_egamma_tightId()[index]) & (1<<1)), 1.0);



    }

}

//
// Main function
//
int MyScanChain::ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents, std::string skimFilePrefix) {

    std::cout << "scanning " << sampleName << std::endl;

    //
    //
    //
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

    //
    // format histograms
    //
    FormatAllEleIdHistograms(sampleName);

    // file loop
    //

    unsigned int nEventsChain=0;
    if(nEvents == -1) nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;
    int i_permille_old = 0;

    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;
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
            int i_permille = (int)floor(1001 * nEventsTotal / float(nEventsChain));
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
            float weight = cms2.evt_scale1fb()*0.01;

            //
            // Fill event level electron histograms
            //
            //FillAllEleIdHistograms(weight, sampleName);

            //
            // loop on hypothesis
            //

            for (size_t h = 0; h < cms2.hyp_type().size(); ++h) {
                //
                // fill basic electron ID histograms
                // depending on the hypothesis
                //


                // apply truth match behavior if ttbar
                //    if (sampleName == "ttbar") {
                //        if(!((abs(cms2.els_mc_id()[index]) == 11) && abs(cms2.els_mc_motherid()[index]) == 24) ) return;
                //    }
                // apply truth match behavior if wjets
                //    if (sampleName == "wjets") {
                //        if(((abs(cms2.els_mc_id()[index]) == 11) && abs(cms2.els_mc_motherid()[index]) == 24) ) return;
                //    }

                DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[h]);
                if (hypType == DILEPTON_EMU && sampleName == "ttbar") {
                    if(abs(cms2.hyp_ll_id()[h]) == 13) {
                        if (leptonIsFromW(cms2.hyp_ll_index()[h], cms2.hyp_ll_id()[h]) > 0 )
                            FillAllEleIdHistograms(cms2.hyp_lt_index()[h], weight, sampleName, h);
                    }                if(abs(cms2.hyp_lt_id()[h]) == 13) {
                        if (leptonIsFromW(cms2.hyp_lt_index()[h], cms2.hyp_lt_id()[h]) > 0 )
                            FillAllEleIdHistograms(cms2.hyp_ll_index()[h], weight, sampleName, h);
                    }
                }       

                if (hypType == DILEPTON_EMU && sampleName == "wjets") {
                    if(abs(cms2.hyp_ll_id()[h]) == 13) {
                        if (leptonIsFromW(cms2.hyp_ll_index()[h], cms2.hyp_ll_id()[h]) > 0 )
                            FillAllEleIdHistograms(cms2.hyp_lt_index()[h], weight, sampleName, h);
                    }
                    if(abs(cms2.hyp_lt_id()[h]) == 13) {
                        if (leptonIsFromW(cms2.hyp_lt_index()[h], cms2.hyp_lt_id()[h]) > 0 )
                            FillAllEleIdHistograms(cms2.hyp_ll_index()[h], weight, sampleName, h);
                    }

                }

                if (hypType == DILEPTON_EMU && sampleName == "QCDpt30") {
                    if(abs(cms2.hyp_lt_id()[h]) == 11) {
                        FillAllEleIdHistograms(cms2.hyp_lt_index()[h], weight, sampleName, h);
                    }
                }

                if (sampleName == "eleidval") {
                    FillAllEleIdHistograms(cms2.hyp_lt_index()[h], weight, sampleName, h);
                }                

            } // end loop on hypothesis



        } // end loop on events

        delete f;
    } // end loop on files

    if ( nEventsChain != nEventsTotal ) {
        std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    }

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

