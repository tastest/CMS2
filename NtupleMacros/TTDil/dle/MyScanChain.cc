
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
#include "TProfile.h"

#include "Math/LorentzVector.h"

// CMS2 includes
#include "../../CORE/CMS2.h"

#include "../../CORE/eventSelections.h"
#include "../../CORE/electronSelections.h"
#include "../../CORE/muonSelections.h"
#include "../../CORE/jetSelections.h"
#include "../../CORE/metSelections.h"
#include "../../CORE/ttbarSelections.h"

#include "../../Tools/tools.cc"
#include "../../Tools/goodrun.cc"
#include "../../CORE/utilities.h"

//
// Namespaces
//
using namespace tas;

//
// definitions...
//

enum {
    PASS_TOPLEPTON_PT10,
    PASS_TOPLEPTON_PTLOW,
    PASS_TOPLEPTON_PT20,
    PASS_TOPLEPTON_ASYMPT,
    PASS_TOPLEPTON_PRESELECTION,
    PASS_TOPLEPTON_ID,
    PASS_TOPLEPTON_ISO,
    PASS_TOPLEPTON_OS,
    PASS_TOPLEPTON_MET,
    PASS_TOPLEPTON_PMET,
    PASS_TOPLEPTON_ZVETO,
    PASS_TOPLEPTON_TRIGGER,
};

enum {
    TOPELEID_CAND01,
    TOPELEID_CICSUPERTIGHT,
    TOPELEID_CICSUPERTIGHTREL,
    TOPELEID_CICTIGHT,
    TOPELEID_CICTIGHTREL,
    TOPELEID_CICLOOSE,
    TOPELEID_CICLOOSEREL,
    TOPELEID_CICVERYLOOSE,
    TOPELEID_CICVERYLOOSEREL,
    TOPELEID_VBTF_TOP95,
    TOPELEID_VBTF_TOP90,
    TOPELEID_VBTF_TOP85,
    TOPELEID_VBTF_TOP80,
    TOPELEID_VBTF_TOP70,
    TOPELEID_VBTF_TOP95REL,
    TOPELEID_VBTF_TOP90REL,
    TOPELEID_VBTF_TOP85REL,
    TOPELEID_VBTF_TOP80REL,
    TOPELEID_VBTF_TOP70REL,
    TOPELEID_VBTF35X_TOP90,
    TOPELEID_VBTF35X_TOP85,
    TOPELEID_VBTF35X_TOP80,
    TOPELEID_VBTF35X_TOP70,
    TOPELEID_VBTF35X_TOP90REL,
    TOPELEID_VBTF35X_TOP85REL,
    TOPELEID_VBTF35X_TOP80REL,
    TOPELEID_VBTF35X_TOP70REL,
};

//
// for jets
//

static const char jetbin_names[][128] = { "0j", "1j", "2j", "allj"};

//
// functions
//

float MyScanChain::leptonIsolation(const unsigned int index, const int id)
{
    if (abs(id) == 11) {
        return electronIsolation_rel(index, true);
    }
    else if (abs(id) == 13) {
        return muonIsoValue(index);
    }
    return -999.9;

}



unsigned int MyScanChain::leptonSelectTopNotTopDown(const int id, const unsigned int lepIdx)
{           

    unsigned int cuts_passed = 0;

    //          
    // The basic pre selection for electrons
    //
    //---------------------------------------------------------
    static const unsigned int preSelection = 
        (1ll<<ELEIP_400) |
        (1ll<<ELEETA_250) |
        (1ll<<ELENOMUON_010) |
        (1ll<<ELESEED_ECAL);
    //---------------------------------------------------------


    // electrons
    //
    if (abs(id) == 11) {
        // preselection
        bool passPreSelection = pass_electronSelection(lepIdx, preSelection);
        if (passPreSelection) cuts_passed |= (1ll<<PASS_TOPLEPTON_PRESELECTION);

        // pt
        if (cms2.els_p4()[lepIdx].Pt() > 10.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PT10);
        if (cms2.els_p4()[lepIdx].Pt() > 10.0 && cms2.els_p4()[lepIdx].Pt() < 20.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PTLOW);
        if (cms2.els_p4()[lepIdx].Pt() > 20.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PT20);

    }

    // muons
    //
    if (abs(id) == 13) {

        // preselection
        if (fabs(cms2.mus_d0corr()[lepIdx]) < 0.02 
                && fabs(cms2.mus_p4()[lepIdx].Eta()) < 2.50) cuts_passed |= (1ll<<PASS_TOPLEPTON_PRESELECTION);

        // pt
        if (cms2.mus_p4()[lepIdx].Pt() > 20.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PT20);
        if (cms2.mus_p4()[lepIdx].Pt() > 10.0 && cms2.mus_p4()[lepIdx].Pt() < 20.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PTLOW);
        if (cms2.mus_p4()[lepIdx].Pt() > 10.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PT10);

        if (muonIdNotIsolated(lepIdx, NominalTTbar)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
        if (muonIsoValue(lepIdx) < 0.15) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
    }

    return cuts_passed;

}


unsigned int MyScanChain::leptonSelect(const int id, const unsigned int lepIdx)
{

    unsigned int cuts_passed = 0;

    //
    // The basic TTBar selection
    //
    //---------------------------------------------------------
    static const unsigned int basicSelection =
        (1ll<<ELEIP_400) |
        (1ll<<ELENOTCONV_HITPATTERN) |
        (1ll<<ELEETA_250) |
        (1ll<<ELENOMUON_010) |
        (1ll<<ELESEED_ECAL);
    //---------------------------------------------------------


    // electrons
    //
    if (abs(id) == 11) {
        // test basic part of electron selection before trying
        // different electron ids
        bool passBasicSelection = pass_electronSelection(lepIdx, basicSelection);
        if (passBasicSelection)
        {
            // pt
            if (cms2.els_p4()[lepIdx].Pt() > 10.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PT10);
            if (cms2.els_p4()[lepIdx].Pt() > 10.0 && cms2.els_p4()[lepIdx].Pt() < 20.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PTLOW);
            if (cms2.els_p4()[lepIdx].Pt() > 20.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PT20);


            // cand01
            if (electronId_ == TOPELEID_CAND01) {
                if (electronId_cand(lepIdx, CAND_01)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }


            // vbtf95_top
            if (electronId_ == TOPELEID_VBTF_TOP95 || electronId_ == TOPELEID_VBTF_TOP95REL) {
                unsigned int answer_vbtf = electronId_VBTF(lepIdx, VBTF_TOP95);
                if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_VBTF_TOP95 && (answer_vbtf & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_VBTF_TOP95REL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }

            // vbtf90_top
            if (electronId_ == TOPELEID_VBTF_TOP90 || electronId_ == TOPELEID_VBTF_TOP90REL) {
                unsigned int answer_vbtf = electronId_VBTF(lepIdx, VBTF_TOP90);
                if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_VBTF_TOP90 && (answer_vbtf & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_VBTF_TOP90REL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }
            // vbtf85_top
            if (electronId_ == TOPELEID_VBTF_TOP85 || electronId_ == TOPELEID_VBTF_TOP85REL) {
                unsigned int answer_vbtf = electronId_VBTF(lepIdx, VBTF_TOP85);
                if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_VBTF_TOP85 && (answer_vbtf & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_VBTF_TOP85REL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }
            // vbtf80_top
            if (electronId_ == TOPELEID_VBTF_TOP80 || electronId_ == TOPELEID_VBTF_TOP80REL) {
                unsigned int answer_vbtf = electronId_VBTF(lepIdx, VBTF_TOP80);
                if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_VBTF_TOP80 && (answer_vbtf & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_VBTF_TOP80REL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }
            // vbtf70_top
            if (electronId_ == TOPELEID_VBTF_TOP70 || electronId_ == TOPELEID_VBTF_TOP70REL) {
                unsigned int answer_vbtf = electronId_VBTF(lepIdx, VBTF_TOP70);
                if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_VBTF_TOP70 && (answer_vbtf & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_VBTF_TOP70REL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }

            //
            //
            //
            // vbtf35x_90_top
            if (electronId_ == TOPELEID_VBTF35X_TOP90 || electronId_ == TOPELEID_VBTF35X_TOP90REL) {
                unsigned int answer_vbtf35x_ = electronId_VBTF(lepIdx, VBTF_35X_90);
                if ((answer_vbtf35x_ & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_VBTF35X_TOP90 && (answer_vbtf35x_ & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_VBTF35X_TOP90REL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }
            // vbtf35x_85_top
            if (electronId_ == TOPELEID_VBTF35X_TOP85 || electronId_ == TOPELEID_VBTF35X_TOP85REL) {
                unsigned int answer_vbtf35x_ = electronId_VBTF(lepIdx, VBTF_35X_85);
                if ((answer_vbtf35x_ & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_VBTF35X_TOP85 && (answer_vbtf35x_ & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_VBTF35X_TOP85REL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }
            // vbtf35x_80_top
            if (electronId_ == TOPELEID_VBTF35X_TOP80 || electronId_ == TOPELEID_VBTF35X_TOP80REL) {
                unsigned int answer_vbtf35x_ = electronId_VBTF(lepIdx, VBTF_35X_80);
                if ((answer_vbtf35x_ & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_VBTF35X_TOP80 && (answer_vbtf35x_ & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_VBTF35X_TOP80REL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }
            // vbtf35x_70_top
            if (electronId_ == TOPELEID_VBTF35X_TOP70 || electronId_ == TOPELEID_VBTF35X_TOP70REL) {
                unsigned int answer_vbtf35x_ = electronId_VBTF(lepIdx, VBTF_35X_70);
                if ((answer_vbtf35x_ & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_VBTF35X_TOP70 && (answer_vbtf35x_ & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_VBTF35X_TOP70REL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }

            //
            //
            //
            // CIC super tight
            if (electronId_ == TOPELEID_CICSUPERTIGHT || electronId_ == TOPELEID_CICSUPERTIGHTREL) {
                unsigned int answer_cic = electronId_CIC(lepIdx, CIC_SUPERTIGHT);
                if ((answer_cic & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_CICSUPERTIGHT && (answer_cic & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_CICSUPERTIGHTREL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }

            // CIC tight
            if (electronId_ == TOPELEID_CICTIGHT || electronId_ == TOPELEID_CICTIGHTREL) {
                unsigned int answer_cic = electronId_CIC(lepIdx, CIC_TIGHT);
                if ((answer_cic & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_CICTIGHT && (answer_cic & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_CICTIGHTREL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }

            // CIC loose
            if (electronId_ == TOPELEID_CICLOOSE || electronId_ == TOPELEID_CICLOOSEREL) {
                unsigned int answer_cic = electronId_CIC(lepIdx, CIC_LOOSE);
                if ((answer_cic & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_CICLOOSE && (answer_cic & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_CICLOOSEREL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }

            // CIC very loose
            if (electronId_ == TOPELEID_CICVERYLOOSE || electronId_ == TOPELEID_CICVERYLOOSEREL) {
                unsigned int answer_cic = electronId_CIC(lepIdx, CIC_VERYLOOSE);
                if ((answer_cic & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_CICVERYLOOSE && (answer_cic & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_CICVERYLOOSEREL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }



        }
    }
    // muons
    //
    if (abs(id) == 13) {

        // pt
        if (cms2.mus_p4()[lepIdx].Pt() > 20.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PT20);
        if (cms2.mus_p4()[lepIdx].Pt() > 10.0 && cms2.mus_p4()[lepIdx].Pt() < 20.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PTLOW);
        if (cms2.mus_p4()[lepIdx].Pt() > 10.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PT10);

        if (muonIdNotIsolated(lepIdx, NominalTTbar)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
        if (muonIsoValue(lepIdx) < 0.15) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
    }

    return cuts_passed;

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

void MyScanChain::FormatAllAnaHistograms(std::string sampleName)
{

    // profiles for dEtaIn analysis
    Int_t nBinsProf = 16;
    Float_t pi = acos(-1);
    p1_detain_eeplus_ = new TProfile("p1_detain_eeplus", "p1_detain_eeplus", nBinsProf, -1*pi, pi);
    p1_detain_corrected_eeplus_ = new TProfile("p1_detain_corrected_eeplus", "p1_detain_corrected_eeplus", nBinsProf, -1*pi, pi);
    p1_detain_correctedfit_eeplus_ = new TProfile("p1_detain_correctedfit_eeplus", "p1_detain_correctedfit_eeplus", nBinsProf, -1*pi, pi);

    p1_detain_eeminus_ = new TProfile("p1_detain_eeminus", "p1_detain_eeminus", nBinsProf, -1*pi, pi);
    p1_detain_corrected_eeminus_ = new TProfile("p1_detain_corrected_eeminus", "p1_detain_corrected_eeminus", nBinsProf, -1*pi, pi);
    p1_detain_correctedfit_eeminus_ = new TProfile("p1_detain_correctedfit_eeminus", "p1_detain_correctedfit_eeminus", nBinsProf, -1*pi, pi);

    h1_dx_eeminus_ = new TH1F("h1_dx_eeminus", "h1_dx_eeminus", 40, -10, 10);
    h1_dy_eeminus_ = new TH1F("h1_dy_eeminus", "h1_dy_eeminus", 40, -10, 10);
    h1_dx_eeplus_ = new TH1F("h1_dx_eeplus", "h1_dx_eeplus", 40, -10, 10);
    h1_dy_eeplus_ = new TH1F("h1_dy_eeplus", "h1_dy_eeplus", 40, -10, 10);

    // spike investigations
    h2_spike_scatteret_ = new TH2F("h2_spike_scatteret", "h2_spike_scatteret;R4;E_T (SC)", 150, 0, 1.5, 200, 0, 200);
    h2_spike_scatteret_eid_ = new TH2F("h2_spike_scatteret_eid", "h2_spike_scatteret_eid;R4;E_T (SC)", 150, 0, 1.5, 200, 0, 200);
    h2_spike_scattermet_ = new TH2F("h2_spike_scattermet", "h2_spike_scattermet;R4;tcMET", 150, 0, 1.5, 200, 0, 200);
    h2_spike_scattermet_fixed_ = new TH2F("h2_spike_scattermet_fixed", "h2_spike_scattermet_fixed;R4;tcMET", 150, 0, 1.5, 200, 0, 200);

    // event level
    FormatHist(h1_hyp_njets_, sampleName, "hyp_njets", 10, -0.5, 9.5);
    FormatHist(h1_hyp_presel_njets_, sampleName, "hyp_presel_njets", 10, -0.5, 9.5);
    FormatHist(h1_hyp_presel_looseiso_njets_, sampleName, "hyp_presel_looseiso_njets", 10, -0.5, 9.5);

    FormatHist(h1_hyp_presel_looseiso_lt_passEle10_, sampleName, "hyp_presel_looseiso_lt_passEle10", 2, -0.5, 1.5);
    FormatHist(h1_hyp_presel_looseiso_ll_pt_reducedBias_, sampleName, "hyp_presel_looseiso_ll_pt_reducedBias", 40, 0.0, 200.0);
    FormatHist(h1_hyp_presel_looseiso_ll_pt_reducedBias_passEle10_, sampleName, "hyp_presel_looseiso_ll_pt_reducedBias_passEle10", 40, 0.0, 200.0);
    FormatHist(h1_hyp_presel_looseiso_ll_eta_reducedBias_, sampleName, "hyp_presel_looseiso_ll_eta_reducedBias", 30, -3.0, 3.0);
    FormatHist(h1_hyp_presel_looseiso_ll_eta_reducedBias_passEle10_, sampleName, "hyp_presel_looseiso_ll_eta_reducedBias_passEle10", 30, -3.0, 3.0);

    // f(njets)
    for (unsigned int j = 0; j < 4; ++j) {
        std::string jetbin = jetbin_names[j];
        FormatHist(h1_hyp_tcmet_[j], sampleName, "hyp_tcmet_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_lt_pt_[j], sampleName, "hyp_lt_pt_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_ll_pt_[j], sampleName, "hyp_ll_pt_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_lt_d0corr_[j], sampleName, "hyp_lt_d0corr_" + jetbin, 100, -0.1, 0.1);
        FormatHist(h1_hyp_ll_d0corr_[j], sampleName, "hyp_ll_d0corr_" + jetbin, 100, -0.1, 0.1);

        // test dependence of muon on if it is global+tracker or just any muon
        FormatHist(h1_hyp_allmuon_d0corr_[j], sampleName, "hyp_allmuon_d0corr_" + jetbin, 100, -0.1, 0.1);
        FormatHist(h1_hyp_glbtrkmuon_d0corr_[j], sampleName, "hyp_glbtrkmuon_d0corr_" + jetbin, 100, -0.1, 0.1);

        FormatHist(h1_hyp_presel_mll_[j], sampleName, "hyp_presel_mll_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_presel_pt_[j], sampleName, "hyp_presel_pt_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_presel_tcmet_[j], sampleName, "hyp_presel_tcmet_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_presel_lt_iso_[j], sampleName, "hyp_presel_lt_iso_" + jetbin, 20, 0.0, 2.0);
        FormatHist(h1_hyp_presel_ll_iso_[j], sampleName, "hyp_presel_ll_iso_" + jetbin, 20, 0.0, 2.0);
        FormatHist(h1_hyp_presel_lt_pt_[j], sampleName, "hyp_presel_lt_pt_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_presel_ll_pt_[j], sampleName, "hyp_presel_ll_pt_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_presel_lt_d0corr_[j], sampleName, "hyp_presel_lt_d0corr_" + jetbin, 100, -0.1, 0.1);
        FormatHist(h1_hyp_presel_ll_d0corr_[j], sampleName, "hyp_presel_ll_d0corr_" + jetbin, 100, -0.1, 0.1);

        FormatHist(h1_hyp_presel_looseiso_mll_[j], sampleName, "hyp_presel_looseiso_mll_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_presel_looseiso_pt_[j], sampleName, "hyp_presel_looseiso_pt_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_presel_looseiso_tcmet_[j], sampleName, "hyp_presel_looseiso_tcmet_" + jetbin, 40, 0.0, 200.0);

    }   

}

void MyScanChain::FormatAllDYEstHistograms(std::string sampleName)
{

    for (unsigned int j = 0; j < 4; ++j) {
        std::string jetbin = jetbin_names[j];
        FormatHist(h1_dyest_mll_met_[j], sampleName, "dyest_mll_met_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_dyest_mll_nomet_[j], sampleName, "dyest_mll_nomet_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_dyest_met_in_[j], sampleName, "dyest_met_in_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_dyest_met_out_[j], sampleName, "dyest_met_out_" + jetbin, 40, 0.0, 200.0);
    }

}

void MyScanChain::FillAllDYEstHistograms(const unsigned int h, const float &weight, const unsigned int njet)
{

    // which hypothesis type
    DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[h]); 

    // which jet bin to fill
    unsigned int jetbin = 0;
    if (njet == 1) jetbin = 1;
    if (njet >= 2) jetbin = 2;

    // fill the mass histogram
    float mass = cms2.hyp_p4()[h].mass();
    Fill(h1_dyest_mll_nomet_[jetbin], hypType, mass, weight);
    Fill(h1_dyest_mll_nomet_[3], hypType, mass, weight);
    if (passMetAsIs_OF20_SF30(cms2.evt_pfmet(), hypType)) {
        Fill(h1_dyest_mll_met_[jetbin], hypType, mass, weight);
        Fill(h1_dyest_mll_met_[3], hypType, mass, weight);
    }
    // fill the met histograms for "in" and "out" regions
    float mymet = cms2.evt_pfmet();
    if (inZmassWindow(mass)) {
        Fill(h1_dyest_met_in_[jetbin], hypType, mymet, weight);
        Fill(h1_dyest_met_in_[3], hypType, mymet, weight);
    }
    else {
        Fill(h1_dyest_met_out_[jetbin], hypType, mymet, weight);
        Fill(h1_dyest_met_out_[3], hypType, mymet, weight);
    }
}

//
// Spike Check
//

void MyScanChain::SpikeCheck(const float &weight)
{

    // require sc et > 10 GeV
    // and ecal seeded bit set
    static const cuts_t selection =
        (1ll<<ELESCET_010) |
        (1ll<<ELESEED_ECAL);

    // loop on electrons
    for (size_t i = 0; i < cms2.evt_nels(); ++i) {

        // restrict to barrel with valid sc index
        if (fabs(cms2.els_etaSC()[i]) > 1.5) continue;
        int scidx = cms2.els_scindex()[i];
        if (scidx == -1) continue;

        // check basic selection passed
        cuts_t cuts_passed = electronSelection(i); 
        if (!pass_electronSelectionCompareMask(cuts_passed, selection)) continue;

        // get the r4 variable
        const float r4 = (cms2.scs_e1x3()[scidx] + cms2.scs_e3x1()[scidx] 
                            - 2*cms2.scs_eMax()[scidx])/cms2.scs_eMax()[scidx];

        // fill histograms
        h2_spike_scattermet_->Fill(r4, cms2.evt_tcmet());
        h2_spike_scatteret_->Fill(r4, cms2.els_eSC()[i]/cosh(cms2.els_etaSC()[i]));
        if (pass_electronSelectionCompareMask(cuts_passed, (1ll<<ELENOSPIKE_SWISS005))) {
            h2_spike_scattermet_fixed_->Fill(r4, cms2.evt_tcmet());
        }
        if (pass_electronSelectionCompareMask(cuts_passed, (1ll<<ELEID_VBTF_35X_90))) {
            h2_spike_scatteret_eid_->Fill(r4, cms2.els_eSC()[i]/cosh(cms2.els_etaSC()[i]));
        }



    }

}

//
// Debugging selections
//

void MyScanChain::Debug(const float &weight)
{

    for (size_t i = 0; i < cms2.evt_nels(); ++i) {

        // get result
        bool pass = false;
        cuts_t cuts_passed = electronSelection(i);
        if ((cuts_passed & electronSelection_ttbarV1) == electronSelection_ttbarV1) pass = true;

        // print debug
    
        std::cout << bool(cuts_passed & (1ll<<ELEID_VBTF_35X_90)) << "\t";;
        std::cout << bool(cuts_passed & (1ll<<ELEIP_400)) << "\t";
        std::cout << bool(cuts_passed & (1ll<<ELEISO_REL015)) << "\t";
        std::cout << bool(cuts_passed & (1ll<<ELENOMUON_010)) << "\t";
        std::cout << bool(cuts_passed & (1ll<<ELENOTCONV_HITPATTERN)) << "\t";
        std::cout << bool(cuts_passed & (1ll<<ELENOTCONV_DISTDCOT002)) << "\t";
        std::cout << bool(cuts_passed & (1ll<<ELESCET_010)) << "\t";
        std::cout << bool(cuts_passed & (1ll<<ELEPT_010)) << "\t";
        std::cout << bool(cuts_passed & (1ll<<ELEETA_250)) << "\t";
        std::cout << bool(cuts_passed & (1ll<<ELESEED_ECAL)) << std::endl;

    }


}

//
// Early data studies
//

void MyScanChain::TopNotTopDown(const float &weight)
{

    //
    // loop on hypothesis
    //

    for (size_t h = 0; h < cms2.hyp_type().size(); ++h) {

        //
        // find out what cuts passed
        //

        // lepton level cuts
        unsigned int cuts_passed_ll = leptonSelectTopNotTopDown(cms2.hyp_ll_id()[h], cms2.hyp_ll_index()[h]);
        unsigned int cuts_passed_lt = leptonSelectTopNotTopDown(cms2.hyp_lt_id()[h], cms2.hyp_lt_index()[h]);

        // hypothesis level cuts
        unsigned int cuts_passed_hyp = 0;
        if (cms2.hyp_lt_id()[h] * cms2.hyp_ll_id()[h] < 0) cuts_passed_hyp |= (1ll<<PASS_TOPLEPTON_OS);

        //
        // fill plots before selecting good hyps
        //

        // determine which has highest pt

        // ll meets pt low and lt meets pt high
        if (    ((cuts_passed_ll & (1ll<<PASS_TOPLEPTON_PT10)) == (1ll<<PASS_TOPLEPTON_PT10) &&
                    (cuts_passed_lt & (1ll<<PASS_TOPLEPTON_PT20)) == (1ll<<PASS_TOPLEPTON_PT20)) ||
                // lt meets pt low and ll meets pt high
                ((cuts_passed_lt & (1ll<<PASS_TOPLEPTON_PT10)) == (1ll<<PASS_TOPLEPTON_PT10) &&
                 (cuts_passed_ll & (1ll<<PASS_TOPLEPTON_PT20)) == (1ll<<PASS_TOPLEPTON_PT20)) ) {
            cuts_passed_hyp |= (1ll<<PASS_TOPLEPTON_ASYMPT);
        }

        float lt_iso = 0.0;
        // get hyp type
        DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[h]);
        // get njets
        unsigned int nJetsFound = nJets(h, JETS_TYPE_PF_CORR, JETS_CLEAN_HYP_E_MU, 0.4, 30.0, 2.4);

        // VERY BASIC SELECTION
        // one of the leptons must be above 20 GeV
        // the other can be as low as 10 GeV
        if (    (cuts_passed_hyp & (1ll<<PASS_TOPLEPTON_OS)) == (1ll<<PASS_TOPLEPTON_OS) &&
                (cuts_passed_hyp & (1ll<<PASS_TOPLEPTON_ASYMPT)) == (1ll<<PASS_TOPLEPTON_ASYMPT)) {

            // if pt(lt) > pt(ll)
            if (cms2.hyp_lt_p4()[h].Pt() > cms2.hyp_ll_p4()[h].Pt()) {
                Fill(h1_hyp_lt_pt_[3], hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                Fill(h1_hyp_ll_pt_[3], hypType, cms2.hyp_ll_p4()[h].Pt(), weight);
                Fill(h1_hyp_lt_d0corr_[3], hypType, cms2.hyp_lt_d0corr()[h], weight);
                Fill(h1_hyp_ll_d0corr_[3], hypType, cms2.hyp_ll_d0corr()[h], weight);

            }
            // if pt(ll) >= pt(lt) 
            else if (cms2.hyp_ll_p4()[h].Pt() >= cms2.hyp_lt_p4()[h].Pt()) {
                Fill(h1_hyp_lt_pt_[3], hypType, cms2.hyp_ll_p4()[h].Pt(), weight);
                Fill(h1_hyp_ll_pt_[3], hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                Fill(h1_hyp_lt_d0corr_[3], hypType, cms2.hyp_ll_d0corr()[h], weight);
                Fill(h1_hyp_ll_d0corr_[3], hypType, cms2.hyp_lt_d0corr()[h], weight);

            }

            Fill(h1_hyp_tcmet_[3], hypType, cms2.evt_tcmet(), weight);

            // investigate the agreement between d0 in data and mc
            // for all muons and for when the muon is global and tracker
            if (abs(cms2.hyp_lt_id()[h]) == 13) {
                Fill(h1_hyp_allmuon_d0corr_[3], hypType, cms2.hyp_lt_d0corr()[h], weight);
                if ( (((cms2.mus_type().at(hyp_lt_index()[h])) & (1<<1)) != 0) 
                        && (((cms2.mus_type().at(hyp_lt_index()[h])) & (1<<2)) != 0) ) {
                    Fill(h1_hyp_glbtrkmuon_d0corr_[3], hypType, cms2.hyp_lt_d0corr()[h], weight);
                }
            }
            if (abs(cms2.hyp_ll_id()[h]) == 13) {
                Fill(h1_hyp_allmuon_d0corr_[3], hypType, cms2.hyp_ll_d0corr()[h], weight);
                if ( (((cms2.mus_type().at(hyp_ll_index()[h])) & (1<<1)) != 0) 
                        && (((cms2.mus_type().at(hyp_ll_index()[h])) & (1<<2)) != 0) ) {
                    Fill(h1_hyp_glbtrkmuon_d0corr_[3], hypType, cms2.hyp_ll_d0corr()[h], weight);
                }
            }

        }

        // VERY BASIC SELECTION
        // + PRESELECTION
        bool lt_is_highpt = true;
        if (    (cuts_passed_lt & (1ll<<PASS_TOPLEPTON_PRESELECTION)) == (1ll<<PASS_TOPLEPTON_PRESELECTION) &&
                (cuts_passed_ll & (1ll<<PASS_TOPLEPTON_PRESELECTION)) == (1ll<<PASS_TOPLEPTON_PRESELECTION) &&
                (cuts_passed_hyp & (1ll<<PASS_TOPLEPTON_OS)) == (1ll<<PASS_TOPLEPTON_OS) &&
                (cuts_passed_hyp & (1ll<<PASS_TOPLEPTON_ASYMPT)) == (1ll<<PASS_TOPLEPTON_ASYMPT)) {

            // if pt(lt) > pt(ll)
            if (cms2.hyp_lt_p4()[h].Pt() > cms2.hyp_ll_p4()[h].Pt()) {
                lt_iso = leptonIsolation(cms2.hyp_lt_index()[h], cms2.hyp_lt_id()[h]);
                Fill(h1_hyp_presel_lt_pt_[3], hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                Fill(h1_hyp_presel_ll_pt_[3], hypType, cms2.hyp_ll_p4()[h].Pt(), weight);
                Fill(h1_hyp_presel_lt_iso_[3], hypType, lt_iso, weight);
                Fill(h1_hyp_presel_ll_iso_[3], hypType, leptonIsolation(cms2.hyp_ll_index()[h], cms2.hyp_ll_id()[h]), weight);
                Fill(h1_hyp_presel_lt_d0corr_[3], hypType, cms2.hyp_lt_d0corr()[h], weight);
                Fill(h1_hyp_presel_ll_d0corr_[3], hypType, cms2.hyp_ll_d0corr()[h], weight);

            }
            // if pt(ll) >= pt(lt) 
            else if (cms2.hyp_ll_p4()[h].Pt() >= cms2.hyp_lt_p4()[h].Pt()) {
                lt_is_highpt = false;
                lt_iso = leptonIsolation(cms2.hyp_ll_index()[h], cms2.hyp_ll_id()[h]);
                Fill(h1_hyp_presel_lt_pt_[3], hypType, cms2.hyp_ll_p4()[h].Pt(), weight);
                Fill(h1_hyp_presel_ll_pt_[3], hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                Fill(h1_hyp_presel_lt_iso_[3], hypType, leptonIsolation(cms2.hyp_ll_index()[h], cms2.hyp_ll_id()[h]), weight);
                Fill(h1_hyp_presel_ll_iso_[3], hypType, lt_iso, weight);
                Fill(h1_hyp_presel_lt_d0corr_[3], hypType, cms2.hyp_ll_d0corr()[h], weight);
                Fill(h1_hyp_presel_ll_d0corr_[3], hypType, cms2.hyp_lt_d0corr()[h], weight);

            }

            Fill(h1_hyp_presel_njets_, hypType, nJetsFound, weight);
            Fill(h1_hyp_presel_mll_[3], hypType, cms2.hyp_p4()[h].mass(), weight);
            Fill(h1_hyp_presel_pt_[3], hypType, cms2.hyp_p4()[h].Pt(), weight);
            Fill(h1_hyp_presel_tcmet_[3], hypType, cms2.evt_tcmet(), weight);

        }


        // VERY BASIC SELECTION
        // + PRESELECTION
        // + LOOSE ISO (<0.4 on the high pT leg)
        if (    lt_iso < 0.4 &&
                (cuts_passed_lt & (1ll<<PASS_TOPLEPTON_PRESELECTION)) == (1ll<<PASS_TOPLEPTON_PRESELECTION) &&
                (cuts_passed_ll & (1ll<<PASS_TOPLEPTON_PRESELECTION)) == (1ll<<PASS_TOPLEPTON_PRESELECTION) &&
                (cuts_passed_hyp & (1ll<<PASS_TOPLEPTON_OS)) == (1ll<<PASS_TOPLEPTON_OS) &&
                (cuts_passed_hyp & (1ll<<PASS_TOPLEPTON_ASYMPT)) == (1ll<<PASS_TOPLEPTON_ASYMPT)) {


            Fill(h1_hyp_presel_looseiso_njets_, hypType, nJetsFound, weight);
            Fill(h1_hyp_presel_looseiso_mll_[3], hypType, cms2.hyp_p4()[h].mass(), weight);
            Fill(h1_hyp_presel_looseiso_pt_[3], hypType, cms2.hyp_p4()[h].Pt(), weight);
            Fill(h1_hyp_presel_looseiso_tcmet_[3], hypType, cms2.evt_tcmet(), weight);

            // something on the trigger (8e29?)
            int trigIndx;
            const std::vector<LorentzVector> *trigObjs = 0;
            TString trigName = "HLT_Ele10_LW_L1R";
            vector<TString>::const_iterator begin_it = hlt_trigNames().begin();
            vector<TString>::const_iterator end_it = hlt_trigNames().end();
            vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
            if(found_it != end_it) {
                trigIndx = found_it - begin_it;
                trigObjs = &cms2.hlt_trigObjs_p4()[trigIndx];
            }
            else {
                std::cout << "Cannot find Trigger " << trigName << std::endl;
            }

            // lt pass
            bool lt_passEle10 = false;
            bool ll_passEle10 = false;
            for (unsigned int i = 0; i < trigObjs->size(); ++i) {
                // is there a trigger lorentz vector within dR < 0.1
                // of the candidate?
                float dR_lt = 999.99;
                float dR_ll = 999.99;
                if (lt_is_highpt) {
                    dR_lt = dRbetweenVectors(cms2.hyp_lt_p4()[h], (*trigObjs)[i]);
                    dR_ll = dRbetweenVectors(cms2.hyp_ll_p4()[h], (*trigObjs)[i]);
                } else {
                    dR_lt = dRbetweenVectors(cms2.hyp_ll_p4()[h], (*trigObjs)[i]);
                    dR_ll = dRbetweenVectors(cms2.hyp_lt_p4()[h], (*trigObjs)[i]);
                }
            
                // if so, then it passed
                if (dR_lt < 0.1) lt_passEle10 = true;
                if (dR_ll < 0.1) ll_passEle10 = true;
            }

            Fill(h1_hyp_presel_looseiso_lt_passEle10_, hypType, lt_passEle10, weight);

            if (lt_passEle10) {
                if (lt_is_highpt) {
                    Fill(h1_hyp_presel_looseiso_ll_pt_reducedBias_, hypType, cms2.hyp_ll_p4()[h].Pt(), weight);
                    Fill(h1_hyp_presel_looseiso_ll_eta_reducedBias_, hypType, cms2.hyp_ll_p4()[h].Eta(), weight);
                    if (ll_passEle10) {
                        Fill(h1_hyp_presel_looseiso_ll_pt_reducedBias_passEle10_, hypType, cms2.hyp_ll_p4()[h].Pt(), weight);
                        Fill(h1_hyp_presel_looseiso_ll_eta_reducedBias_passEle10_, hypType, cms2.hyp_ll_p4()[h].Eta(), weight);                    
                    }
                } else {
                    Fill(h1_hyp_presel_looseiso_ll_pt_reducedBias_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                    Fill(h1_hyp_presel_looseiso_ll_eta_reducedBias_, hypType, cms2.hyp_lt_p4()[h].Eta(), weight);
                    if (ll_passEle10) {
                        Fill(h1_hyp_presel_looseiso_ll_pt_reducedBias_passEle10_, hypType, cms2.hyp_lt_p4()[h].Pt(), weight);
                        Fill(h1_hyp_presel_looseiso_ll_eta_reducedBias_passEle10_, hypType, cms2.hyp_lt_p4()[h].Eta(), weight);
                    }
                }
            }

        }


    } // end loop on hyps

}

//
// check deta and correction
//
void MyScanChain::SuperClusterCheck(const float &weight)
{

    for (size_t i = 0; i < cms2.evt_nels(); ++i) {

        if (fabs(cms2.els_etaSC()[i]) > 2.1) continue;
        if (cms2.els_eSC()[i]/cosh(cms2.els_etaSC()[i]) < 10.0) continue;    
        if (cms2.els_scindex()[i] == -1) continue;

        // plot dEtaIn as a function of phi
        // for EE+ and EE-
        // and for the "corrected" dEta too

        // get corrected dEtaIn
        LorentzVector initial_pos = cms2.scs_pos_p4()[cms2.els_scindex()[i]];

        float eta_trackin_at_calo = cms2.els_dEtaIn()[i] + cms2.els_etaSC()[cms2.els_scindex()[i]];
        float theta_trackin_at_calo = 2*atan(exp(-1*eta_trackin_at_calo));
        float phi_trackin_at_calo = cms2.els_dPhiIn()[i] + cms2.els_phiSC()[cms2.els_scindex()[i]];
        
        //float R = 322.5;
        float R = initial_pos.z();
        if (cms2.els_etaSC()[i] < 0) R *= -1;
        float x_trackin_at_calo = R*tan(theta_trackin_at_calo)*cos(phi_trackin_at_calo);
        float y_trackin_at_calo = R*tan(theta_trackin_at_calo)*sin(phi_trackin_at_calo);
        float x_cluster = initial_pos.x();
        float y_cluster = initial_pos.y();
     
        if (fabs(cms2.els_etaSC()[i]) > 1.479) { 
            std::cout << "x, y, z: " << initial_pos.x() << ", " << initial_pos.y() << ", " << initial_pos.z() << std::endl; 
            std::cout << "\t calculated x, y: " << x_trackin_at_calo << "\t" << y_trackin_at_calo << std::endl;
        }

        if (cms2.els_sigmaIEtaIEta()[i] < 0.03) {
            // EE+
            if (cms2.els_etaSC()[i] > 1.479) {
                h1_dx_eeplus_->Fill(x_cluster - x_trackin_at_calo);
                h1_dy_eeplus_->Fill(y_cluster - y_trackin_at_calo);
            }
            // EE-
            if (cms2.els_etaSC()[i] < 1.479) {
                h1_dx_eeminus_->Fill(x_cluster - x_trackin_at_calo);
                h1_dy_eeminus_->Fill(y_cluster - y_trackin_at_calo);
            }
        }

            if (cms2.els_etaSC()[i] > 1.479) {
                LorentzVector corrected_pos = LorentzVector(initial_pos.x() + 0.5, initial_pos.y() - 0.8, initial_pos.z(), 0.0);
                p1_detain_eeplus_->Fill(cms2.els_phiSC()[i], cms2.els_dEtaIn()[i]);
                p1_detain_corrected_eeplus_->Fill(cms2.els_phiSC()[i], cms2.els_dEtaIn()[i] - (initial_pos.eta() - corrected_pos.eta()));
            }
            if (cms2.els_etaSC()[i] < 1.479) {
                LorentzVector corrected_pos = LorentzVector(initial_pos.x(), initial_pos.y() - 0.8, initial_pos.z(), 0.0);
                p1_detain_eeminus_->Fill(cms2.els_phiSC()[i], cms2.els_dEtaIn()[i]);
                p1_detain_corrected_eeminus_->Fill(cms2.els_phiSC()[i], cms2.els_dEtaIn()[i] - (initial_pos.eta() - corrected_pos.eta()));
            }


    }   

}

//
// Main function
//

int MyScanChain::ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents, std::string skimFilePrefix) {

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

    FormatAllAnaHistograms(sampleName);
    FormatAllDYEstHistograms(sampleName);

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
                weight = cms2.evt_scale1fb()*0.01;
            }

            //
            // do good run check
            //

            if (isData) {
                if (!goodrun(cms2.evt_run(), cms2.evt_lumiBlock())) continue;
            }

            //
            // standard event cleaning
            //

            if (!cleaning_standard(isData)) continue;

            //
            // Do early data study
            //

            //TopNotTopDown(weight);
            //SuperClusterCheck(weight);
            
            //
            // Run debugging function
            // 
            //Debug(weight);

            //
            // Run the spike check
            //
            SpikeCheck(weight);

            //
            // define cuts
            //
            
               const unsigned int selection_lt = (1ll<<PASS_TOPLEPTON_PT20) | 
               (1ll<<PASS_TOPLEPTON_ISO) | (1ll<<PASS_TOPLEPTON_ID);
               const unsigned int selection_ll = selection_lt;
               const unsigned int selection_hyp = (1ll<<PASS_TOPLEPTON_OS);

            //
            // loop on hypothesis
            //

            std::vector<unsigned int> hyp_index_selected;
            hyp_index_selected.clear();

            for (size_t h = 0; h < cms2.hyp_type().size(); ++h) {

            //
            // find out what cuts passed
            //

            // find what cuts lt and ll passed
            unsigned int cuts_passed_ll = leptonSelect(cms2.hyp_ll_id()[h], cms2.hyp_ll_index()[h]);
            unsigned int cuts_passed_lt = leptonSelect(cms2.hyp_lt_id()[h], cms2.hyp_lt_index()[h]);
            // hypothesis level cuts
            unsigned int cuts_passed_hyp = 0;
            if (cms2.hyp_lt_id()[h] * cms2.hyp_ll_id()[h] < 0) cuts_passed_hyp |= (1ll<<PASS_TOPLEPTON_OS);

            //
            // apply selection
            //
            if (!((cuts_passed_lt & selection_lt) == selection_lt)) continue;
            if (!((cuts_passed_ll & selection_ll) == selection_ll)) continue;
            if (!((cuts_passed_hyp & selection_hyp) == selection_hyp)) continue;

            //
            // store indices of hypothesese that pass this preselection
            //	
            hyp_index_selected.push_back(h);

            } // end loop on hypothesis

            //
            // require a hypothesis was found
            //

            if (hyp_index_selected.size() == 0) continue;
            int hyp = 0;

            //
            // perform hypothesis disambiguation and get hyp type for selected hyp
            //

            hyp = eventDilIndexByWeightTTDil10(hyp_index_selected);
            DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[hyp]);

            //
            // make requirements of the selected hypothesis
            //

            unsigned int cuts_passed_event = 0;

            // trigger
            // FIXME
            //if (isData) cuts_passed_event |= (1ll<<PASS_TOPLEPTON_TRIGGER);
            //if (!isData && passTriggersMu9orLisoE15(cms2.hyp_type()[hyp])) cuts_passed_event |= (1ll<<PASS_TOPLEPTON_TRIGGER);
            cuts_passed_event |= (1ll<<PASS_TOPLEPTON_TRIGGER);

            // met
            if (passMetAsIs_OF20_SF30(cms2.evt_pfmet(), hypType)) cuts_passed_event |= (1ll<<PASS_TOPLEPTON_MET);

            // projected met
            float pMet = projectedMET(cms2.evt_pfmet(), cms2.evt_pfmetPhi(), hyp);
            if (passMetAsIs_OF20_SF30(pMet, hypType)) cuts_passed_event |= (1ll<<PASS_TOPLEPTON_PMET);

            // z mass window
            if (cms2.hyp_type()[hyp] == 0 || cms2.hyp_type()[hyp] == 3) {
                if (!inZmassWindow(cms2.hyp_p4()[hyp].mass())) cuts_passed_event |= (1ll<<PASS_TOPLEPTON_ZVETO);
            }
            else {
                cuts_passed_event |= (1ll<<PASS_TOPLEPTON_ZVETO);
            }

            //
            // classify the event according to nJets
            //

            unsigned int nJetsFound = nJets(hyp, JETS_TYPE_PF_CORR, JETS_CLEAN_HYP_E_MU, 0.4, 30.0, 2.4);

            //
            // estimate the DY background before the z veto and MET cuts are applied
            //

            if ((cuts_passed_event & (1ll<<PASS_TOPLEPTON_TRIGGER)) == (1ll<<PASS_TOPLEPTON_TRIGGER)) {
                FillAllDYEstHistograms(hyp, weight, nJetsFound);
            }

            //
            // fill analysis results histograms
            //

            // the full standard event selection
            const unsigned int selection_event_full = (1ll<<PASS_TOPLEPTON_TRIGGER) | (1ll<<PASS_TOPLEPTON_MET) | 
                (1ll<<PASS_TOPLEPTON_PMET) | (1ll<<PASS_TOPLEPTON_ZVETO); 

            // if the full standard event selection passed
            // then fill final analysis histograms
            if ((cuts_passed_event & selection_event_full) == selection_event_full) {
                Fill(h1_hyp_njets_, hypType, nJetsFound, weight);
            }
            


        } // end loop on files

    } // end loop on events

    if ( nEventsChain != nEventsTotal ) {
        std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    }

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

