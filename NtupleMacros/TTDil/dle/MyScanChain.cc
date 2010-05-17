
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

#include "../../CORE/eventSelections.h"
#include "../../CORE/electronSelections.h"
#include "../../CORE/muonSelections.h"
#include "../../CORE/jetSelections.h"
#include "../../CORE/metSelections.h"
#include "../../CORE/ttbarSelections.h"

#include "../../Tools/tools.cc"

//
// Namespaces
//
using namespace tas;

//
// definitions...
//

enum {

    PASS_TOPLEPTON_PT7,
    PASS_TOPLEPTON_PTLOW,
    PASS_TOPLEPTON_PTHIGH,
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
            if (cms2.els_p4()[lepIdx].Pt() > 7.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PT7);
            if (cms2.els_p4()[lepIdx].Pt() > 10.0 && cms2.els_p4()[lepIdx].Pt() < 20.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PTLOW);
            if (cms2.els_p4()[lepIdx].Pt() > 20.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PTHIGH);


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
                unsigned int answer_vbtf35x_ = electronId_VBTF(lepIdx, VBTF_35X_TOP90);
                if ((answer_vbtf35x_ & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_VBTF35X_TOP90 && (answer_vbtf35x_ & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_VBTF35X_TOP90REL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }
            // vbtf35x_85_top
            if (electronId_ == TOPELEID_VBTF35X_TOP85 || electronId_ == TOPELEID_VBTF35X_TOP85REL) {
                unsigned int answer_vbtf35x_ = electronId_VBTF(lepIdx, VBTF_35X_TOP85);
                if ((answer_vbtf35x_ & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_VBTF35X_TOP85 && (answer_vbtf35x_ & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_VBTF35X_TOP85REL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }
            // vbtf35x_80_top
            if (electronId_ == TOPELEID_VBTF35X_TOP80 || electronId_ == TOPELEID_VBTF35X_TOP80REL) {
                unsigned int answer_vbtf35x_ = electronId_VBTF(lepIdx, VBTF_35X_TOP80);
                if ((answer_vbtf35x_ & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ID);
                if (electronId_ == TOPELEID_VBTF35X_TOP80 && (answer_vbtf35x_ & (1ll<<ELEID_ISO)) == (1ll<<ELEID_ISO)) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
                if (electronId_ == TOPELEID_VBTF35X_TOP80REL && electronIsolation_rel(lepIdx, true) < 0.1) cuts_passed |= (1ll<<PASS_TOPLEPTON_ISO);
            }
            // vbtf35x_70_top
            if (electronId_ == TOPELEID_VBTF35X_TOP70 || electronId_ == TOPELEID_VBTF35X_TOP70REL) {
                unsigned int answer_vbtf35x_ = electronId_VBTF(lepIdx, VBTF_35X_TOP70);
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
        if (cms2.mus_p4()[lepIdx].Pt() > 20.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PTHIGH);
        if (cms2.mus_p4()[lepIdx].Pt() > 10.0 && cms2.mus_p4()[lepIdx].Pt() < 20.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PTLOW);
        if (cms2.mus_p4()[lepIdx].Pt() > 7.0) cuts_passed |= (1ll<<PASS_TOPLEPTON_PT7);

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

    // event level
    FormatHist(h1_hyp_njets_, sampleName, "hyp_njets", 10, -0.5, 9.5);

    // f(njets)
    for (unsigned int j = 0; j < 4; ++j) {
        std::string jetbin = jetbin_names[j];
        FormatHist(h1_hyp_mll_[j], sampleName, "hyp_mll_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_pt_[j], sampleName, "hyp_pt_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_tcmet_[j], sampleName, "hyp_tcmet_" + jetbin, 40, 0.0, 200.0);
        FormatHist(h1_hyp_pfmet_[j], sampleName, "hyp_pfmet_" + jetbin, 40, 0.0, 200.0);
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
            // standard event cleaning
            //

            if (!cleaning_standard(isData)) continue;

            //
            // define cuts
            //

            const unsigned int selection_lt = (1ll<<PASS_TOPLEPTON_PTHIGH) | (1ll<<PASS_TOPLEPTON_ISO) | (1ll<<PASS_TOPLEPTON_ID);
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
                // fill plots before selecting good hyps
                //

                // some early data studies
                // w.r.t. pt > 7 lt and ll
                if ( ((cuts_passed_ll & (1ll<<PASS_TOPLEPTON_PT7)) == (1ll<<PASS_TOPLEPTON_PT7)) &&
                        ((cuts_passed_lt & (1ll<<PASS_TOPLEPTON_PT7)) == (1ll<<PASS_TOPLEPTON_PT7)) ) {

                    // get hyp type
                    DileptonHypType hypType = hyp_typeToHypType(cms2.hyp_type()[h]);

                    // mll
                    float mass = cms2.hyp_p4()[h].mass();
                    Fill(h1_hyp_mll_[3], hypType, cms2.hyp_p4()[h].mass(), weight);
                    Fill(h1_hyp_mll_[3], hypType, cms2.hyp_p4()[h].Pt(), weight);
                    Fill(h1_hyp_tcmet_[3], hypType, cms2.evt_tcmet(), weight);
                    Fill(h1_hyp_pfmet_[3], hypType, cms2.evt_pfmet(), weight);

                }

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

