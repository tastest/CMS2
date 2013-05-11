#define __CMS2_SLIM__ // this is now assumed to be default since moving forward all data is slimmed
#include "myBabyMaker.h"

// C++ includes
#include <iostream>
#include <fstream>
#include <set>
#include <exception>
#include <string>

// ROOT includes
#include "TSystem.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TChainElement.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Math/VectorUtil.h"
#include "TObjArray.h"
#include "TString.h"
#include "TVector2.h"
#include "TDatabasePDG.h"
#include "TBenchmark.h"

// TAS includes
// This is for those using a makefile 
// for linking CORE and Tools as a standalone librabry 
// you need to define __NON_ROOT_BUILD__ in your build script
// (i.e. g++ ... -D__NON_ROOT_BUILD__ ... )
#ifdef __NON_ROOT_BUILD__  
#include "CMS2.h"
#include "utilities.h"
#include "electronSelections.h"
#include "electronSelectionsParameters.h"
#include "eventSelections.h"
#include "jetSelections.h"
#include "metSelections.h"
#include "MITConversionUtilities.h"
#include "muonSelections.h"
#include "trackSelections.h"
#include "triggerUtils.h"
#include "goodrun.h"
#include "mcSelections.h"
#include "susySelections.h"
#include "ssSelections.h"
#include "ttvSelections.h"
#include "jetcorr/FactorizedJetCorrector.h"
#else
// for compiling in ACLiC (.L myBabyMaker.c++ method)
// since the source files are included
#include "CMS2.cc"  
#ifndef __CINT__
#include "../CORE/utilities.cc"
#include "../CORE/electronSelections.cc"
#include "../CORE/electronSelectionsParameters.cc"
#include "../CORE/eventSelections.cc"
#include "../CORE/jetSelections.cc"
#include "../CORE/metSelections.cc"
#include "../CORE/MITConversionUtilities.cc"
#include "../CORE/muonSelections.cc"
#include "../CORE/trackSelections.cc"
#include "../CORE/triggerUtils.cc"
#include "../Tools/goodrun.cc"
#include "../CORE/mcSelections.cc"
#include "../CORE/ssSelections.cc"
#include "../CORE/susySelections.cc"
#include "../CORE/jetcorr/FactorizedJetCorrector.h"
#include "../CORE/ttvSelections.cc"
#endif // __CINT__
#endif // __NON_ROOT_BUILD__

// namespaces
using namespace std;
using namespace tas;

#ifndef __CINT__
bool header1 = false;
bool header2 = false;

void PrintTriggerDebugHeader(string outfileName)
{
    int width  = 7;
    ofstream outfile( Form("triggerStudy/%s", outfileName.c_str() ), ios::app );
    outfile  
        <<  setw(width) << "itrg"
        <<  setw(width) << "id"
        <<  setw(width) << "match"
        <<  setw(width) << "matchId"
        <<  setw(width) << "dr"
        <<  setw(width) << "lep pt"
        <<  setw(width) << "trg pt"
        <<  setw(width) << "lep eta"
        <<  setw(width) << "trg eta"
        <<  setw(width) << "lep phi"
        <<  setw(width) << "trg phi" 
        << "\t"         << "trigString" << endl << endl;
    outfile.close();
}

void PrintTriggerDebugLine
(
    int itrg, 
    int id, 
    bool match, 
    bool matchId, 
    double dr, 
    const LorentzVector& lepton_p4, 
    const LorentzVector& p4tr, 
    const string& trigString, 
    int nTrig, 
    const string& outfileName
    )
{
    ofstream outfile( Form("triggerStudy/%s", outfileName.c_str() ), ios::app );

    int precis = 2;
    int width  = 7;
    outfile.setf( ios::fixed, ios::floatfield );
    outfile << setprecision(precis) << setw(width) << setfill(' ') << itrg
            << setprecision(precis) << setw(width) << setfill(' ') << id
            << setprecision(precis) << setw(width) << setfill(' ') << match
            << setprecision(precis) << setw(width) << setfill(' ') << matchId
            << setprecision(precis) << setw(width) << setfill(' ') << dr
            << setprecision(precis) << setw(width) << setfill(' ') << lepton_p4.pt()
            << setprecision(precis) << setw(width) << setfill(' ') << p4tr.pt()
            << setprecision(precis) << setw(width) << setfill(' ') << lepton_p4.eta()
            << setprecision(precis) << setw(width) << setfill(' ') << p4tr.eta()
            << setprecision(precis) << setw(width) << setfill(' ') << lepton_p4.phi()
            << setprecision(precis) << setw(width) << setfill(' ') << p4tr.phi()
            << "\t" << trigString << endl;
    if( itrg == nTrig-1 ) outfile << endl;
    outfile.close();
}

bool found_ele8 = false;
bool found_ele8_CaloIdL_TrkIdVL = false;
bool found_ele8_CaloIdL_CaloIsoVL = false;
bool found_ele17_CaloIdL_CaloIsoVL = false;
bool found_ele8_CaloIdL_CaloIsoVL_Jet40 = false;

// function for dR matching offline letpon to trigger object 
pair<int, float> TriggerMatch( LorentzVector lepton_p4, const char* trigString, double dR_cut = 0.4, int pid = 11 )
{
    float dR_min = numeric_limits<float>::max();
    dR_min = 99.0;
    int nTrig = nHLTObjects(trigString);
    if (nTrig > 0) {
        bool match   = false;
        bool matchId = false;

        for (int itrg=0; itrg<nTrig; itrg++) 
        {
            LorentzVector p4tr = p4HLTObject( trigString, itrg );
            int id             = idHLTObject( trigString, itrg );
            double dr = ROOT::Math::VectorUtil::DeltaR( lepton_p4, p4tr);
            if ( dr < dR_cut ){
                match = true;
                if( abs(id) == abs(pid) ) matchId = true;
            }
            if (dr < dR_min) dR_min = dr;


            //////////////////////////
            // Debug Mixed Triggers //
            //////////////////////////
            vector<string> triggers;
            triggers.push_back("HLT_Ele8");
            triggers.push_back("HLT_Ele8_CaloIdL_TrkIdVL");
            triggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL");
            triggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL");
            triggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40");
            for (unsigned int i=0; i < triggers.size(); i++) {
                if (
                    ( strcmp( trigString, Form( "%s_%s", triggers.at(i).c_str(), "v1" ) ) == 0 ) || 
                    ( strcmp( trigString, Form( "%s_%s", triggers.at(i).c_str(), "v2" ) ) == 0 )
                    ) {

                    //
                    if(
                        ( strcmp( trigString, "HLT_Ele8_v1") == 0 ) || 
                        ( strcmp( trigString, "HLT_Ele8_v2") == 0 )
                        ){
                        if(found_ele8 == false){
                            found_ele8 = true;
                            PrintTriggerDebugHeader( Form("%s.txt",triggers.at(i).c_str()) );
                        }
                    }

                    //
                    if(
                        ( strcmp( trigString, "HLT_Ele8_CaloIdL_TrkIdVL_v1") == 0 ) || 
                        ( strcmp( trigString, "HLT_Ele8_CaloIdL_TrkIdVL_v2") == 0 )
                        ){
                        if(found_ele8_CaloIdL_TrkIdVL == false){
                            found_ele8_CaloIdL_TrkIdVL = true;
                            PrintTriggerDebugHeader( Form("%s.txt",triggers.at(i).c_str()) );
                        }
                    }

                    //
                    if(
                        ( strcmp( trigString, "HLT_Ele8_CaloIdL_CaloIsoVL_v1") == 0 ) || 
                        ( strcmp( trigString, "HLT_Ele8_CaloIdL_CaloIsoVL_v2") == 0 )
                        ){
                        if(found_ele8_CaloIdL_CaloIsoVL == false){
                            found_ele8_CaloIdL_CaloIsoVL = true;
                            PrintTriggerDebugHeader( Form("%s.txt",triggers.at(i).c_str()) );
                        }
                    }

                    //
                    if(
                        ( strcmp( trigString, "HLT_Ele17_CaloIdL_CaloIsoVL_v1") == 0 ) ||
                        ( strcmp( trigString, "HLT_Ele17_CaloIdL_CaloIsoVL_v2") == 0 )
                        ){
                        if(found_ele17_CaloIdL_CaloIsoVL == false){
                            found_ele17_CaloIdL_CaloIsoVL = true;
                            PrintTriggerDebugHeader( Form("%s.txt",triggers.at(i).c_str()) );
                        }
                    }

                    //
                    if(
                        ( strcmp( trigString, "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1") == 0 ) || 
                        ( strcmp( trigString, "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2") == 0 )
                        ){
                        if(found_ele8_CaloIdL_CaloIsoVL_Jet40 == false){
                            found_ele8_CaloIdL_CaloIsoVL_Jet40 = true;
                            PrintTriggerDebugHeader( Form("%s.txt",triggers.at(i).c_str()) );
                        }
                    }

                    PrintTriggerDebugLine( itrg, id, match, matchId, dr, lepton_p4, p4tr, trigString, nTrig, Form("%s.txt", triggers.at(i).c_str() ) );
                }
            }

        } // end loop on triggers
        if(matchId){
            nTrig = 3;
        }
        else {
            if (match) {
                nTrig = 2;
            } 
            else {
                nTrig = 1;
            }
        }
    }
    pair<int, float> answer;
    answer.first  = nTrig;
    answer.second = dR_min;
    return answer;
}

// struct for trigger matching
struct triggerMatchStruct 
{
    triggerMatchStruct(int NumHLTObjs, float deltaR, int vers, int hltps = -1, int l1ps = -1) 
        : nHLTObjects_(NumHLTObjs)
        , dR_(deltaR)
        , version_(vers)
        , hltps_(hltps)
        , l1ps_(l1ps)
    {}

    int nHLTObjects_;
    float dR_;
    int version_;
    int hltps_;
    int l1ps_;  // not used yet
};

// wrapper around TriggerMatch that takes a TRegExp for matching a class of triggers
triggerMatchStruct MatchTriggerClass(LorentzVector lepton_p4, TPMERegexp& regexp, int pid = 11, double dR_cut = 0.4)
{
    std::pair<int, float> triggerMatchValues = make_pair (0, 99.);
    triggerMatchStruct triggerMatchInfo = triggerMatchStruct(triggerMatchValues.first, triggerMatchValues.second, -1, -1);
    
    unsigned int loopCounts = 0;
    for (unsigned int tidx = 0; tidx < cms2.hlt_trigNames().size(); tidx++) {
        if (regexp.Match(cms2.hlt_trigNames().at(tidx)) == 0)
            continue;

        ++loopCounts;

        // get lepton-trigger matching information
        triggerMatchValues = TriggerMatch(lepton_p4, cms2.hlt_trigNames().at(tidx).Data(), dR_cut, pid);

        int version = -1;
        TString tversion = regexp[1];
        if (tversion.IsDigit())
            version = tversion.Atoi();
        int hltprescale = HLT_prescale(cms2.hlt_trigNames().at(tidx).Data());

        triggerMatchInfo = triggerMatchStruct(triggerMatchValues.first, triggerMatchValues.second, version, hltprescale);
    }

    //assert (loopCounts < 2);
    // check that we counted at least two matches
    if (loopCounts >= 2)
    {
        throw std::logic_error(Form("[FR baby maker]: MatchTriggerClass - looper cout is greater than two. loopCounts = %u", loopCounts));
    }
    return triggerMatchInfo;
}

//////////////////////////////
// THIS NEEDS TO BE IN CORE //
//////////////////////////////

struct DorkyEventIdentifier
{
    // this is a workaround for not having unique event id's in MC
    unsigned long int run, event,lumi;
    bool operator < (const DorkyEventIdentifier &) const;
    bool operator == (const DorkyEventIdentifier &) const;
};

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
    if (run != other.run)
        return run < other.run;
    if (event != other.event)
        return event < other.event;
    if(lumi != other.lumi)
        return lumi < other.lumi;
    return false;
}

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
    if (run != other.run)
        return false;
    if (event != other.event)
        return false;
    return true;
}

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id)
{
    std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
        already_seen.insert(id);
    return !ret.second;
}

// transverse mass
float Mt( LorentzVector p4, float met, float met_phi )
{
    return sqrt( 2*met*( p4.pt() - ( p4.Px()*cos(met_phi) + p4.Py()*sin(met_phi) ) ) );
}

#endif // __CINT__

// set good run list
void myBabyMaker::SetGoodRunList(const char* fileName, bool goodRunIsJson)
{
    if (!std::string(fileName).empty())
    {
        if (goodRunIsJson)
            set_goodrun_file_json(fileName);
        else
            set_goodrun_file(fileName);

        goodrun_is_json = goodRunIsJson;
    }
}

// lepton effecitve area
// calculate Effective area (updated to value from Egamma)
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaEARhoCorrection
// Topic revision: r12 - 28-Nov-2012
Float_t EffectiveArea(float eta, float cone_size, int eormu, bool use_tight)
{
    float etaAbs = fabs(eta);
    float eff_area = 0.0;

    if (abs(eormu) == 11)
    {
        if (abs(cone_size - 0.3) < 0.01)
        {
            eff_area = fastJetEffArea03_v2(etaAbs); 
        }
        else if (abs(cone_size - 0.4) < 0.01)
        {
            eff_area = fastJetEffArea04_v2(etaAbs); 
        }
    }
    else if (abs(eormu) == 13)
    {
        if (abs(cone_size - 0.3) < 0.01)
        {
            if (use_tight)
            {
                if (etaAbs < 1.0)      eff_area = 0.207;
                else if (etaAbs < 1.5) eff_area = 0.183;
                else if (etaAbs < 2.0) eff_area = 0.177;
                else if (etaAbs < 2.2) eff_area = 0.271;
                else if (etaAbs < 2.3) eff_area = 0.348;
                else if (etaAbs < 2.4) eff_area = 0.246;
            }
            else
            {
                if (etaAbs < 1.0)      eff_area = 0.382;
                else if (etaAbs < 1.5) eff_area = 0.317;
                else if (etaAbs < 2.0) eff_area = 0.242;
                else if (etaAbs < 2.2) eff_area = 0.326;
                else if (etaAbs < 2.3) eff_area = 0.462;
                else if (etaAbs < 2.4) eff_area = 0.372;
            }
        }
        else if (abs(cone_size - 0.4) < 0.01)
        {
            if (use_tight)
            {
                if (etaAbs < 1.0)      eff_area = 0.340;
                else if (etaAbs < 1.5) eff_area = 0.310;
                else if (etaAbs < 2.0) eff_area = 0.315;
                else if (etaAbs < 2.2) eff_area = 0.415;
                else if (etaAbs < 2.3) eff_area = 0.658;
                else if (etaAbs < 2.4) eff_area = 0.405;
            }
            else
            {
                if (etaAbs < 1.0)      eff_area = 0.674;
                else if (etaAbs < 1.5) eff_area = 0.565;
                else if (etaAbs < 2.0) eff_area = 0.442;
                else if (etaAbs < 2.2) eff_area = 0.515;
                else if (etaAbs < 2.3) eff_area = 0.821;
                else if (etaAbs < 2.4) eff_area = 0.660;
            }
        }

    }

    // done
    return eff_area;
}

// muons only
Float_t EffectiveArea_nh(float eta, float cone_size, bool use_tight)
{
    float etaAbs = fabs(eta);
    float eff_area = 0.0;

    if (abs(cone_size - 0.3) < 0.01)
    {
        if (use_tight)
        {
            if (etaAbs < 1.0)      eff_area = 0.093;
            else if (etaAbs < 1.5) eff_area = 0.116;
            else if (etaAbs < 2.0) eff_area = 0.144;
            else if (etaAbs < 2.2) eff_area = 0.101;
            else if (etaAbs < 2.3) eff_area = 0.105;
            else if (etaAbs < 2.4) eff_area = 0.178;
        }

        else
        {
            if (etaAbs < 1.0)      eff_area = 0.107;
            else if (etaAbs < 1.5) eff_area = 0.141;
            else if (etaAbs < 2.0) eff_area = 0.159;
            else if (etaAbs < 2.2) eff_area = 0.102;
            else if (etaAbs < 2.3) eff_area = 0.096;
            else if (etaAbs < 2.4) eff_area = 0.104;
        }
    }
    else if (abs(cone_size - 0.4) < 0.01)
    {
        if (use_tight)
        {
            if (etaAbs < 1.0)      eff_area = 0.140;
            else if (etaAbs < 1.5) eff_area = 0.204;
            else if (etaAbs < 2.0) eff_area = 0.224;
            else if (etaAbs < 2.2) eff_area = 0.229;
            else if (etaAbs < 2.3) eff_area = 0.322;
            else if (etaAbs < 2.4) eff_area = 0.178;
        }
        else
        {
            if (etaAbs < 1.0)      eff_area = 0.166;
            else if (etaAbs < 1.5) eff_area = 0.259;
            else if (etaAbs < 2.0) eff_area = 0.247;
            else if (etaAbs < 2.2) eff_area = 0.220;
            else if (etaAbs < 2.3) eff_area = 0.340;
            else if (etaAbs < 2.4) eff_area = 0.216;
        }
    }

    // done
    return eff_area;
}

// muons only
Float_t EffectiveArea_em(float eta, float cone_size, bool use_tight)
{
    float etaAbs = fabs(eta);
    float eff_area = 0.0;

    if (abs(cone_size - 0.3) < 0.01)
    {
        if (use_tight)
        {
            if (etaAbs < 1.0)      eff_area = 0.118;
            else if (etaAbs < 1.5) eff_area = 0.053;
            else if (etaAbs < 2.0) eff_area = 0.015;
            else if (etaAbs < 2.2) eff_area = 0.112;
            else if (etaAbs < 2.3) eff_area = 0.302;
            else if (etaAbs < 2.4) eff_area = 0.251;
        }
        else
        {
            if (etaAbs < 1.0)      eff_area = 0.274;
            else if (etaAbs < 1.5) eff_area = 0.161;
            else if (etaAbs < 2.0) eff_area = 0.079;
            else if (etaAbs < 2.2) eff_area = 0.168;
            else if (etaAbs < 2.3) eff_area = 0.359;
            else if (etaAbs < 2.4) eff_area = 0.294;
        }
    }
    else if (abs(cone_size - 0.4) < 0.01)
    {
        if (use_tight)
        {
            if (etaAbs < 1.0)      eff_area = 0.200;
            else if (etaAbs < 1.5) eff_area = 0.109;
            else if (etaAbs < 2.0) eff_area = 0.087;
            else if (etaAbs < 2.2) eff_area = 0.184;
            else if (etaAbs < 2.3) eff_area = 0.425;
            else if (etaAbs < 2.4) eff_area = 0.350;
        }
        else
        {
            if (etaAbs < 1.0)      eff_area = 0.504;
            else if (etaAbs < 1.5) eff_area = 0.306;
            else if (etaAbs < 2.0) eff_area = 0.198;
            else if (etaAbs < 2.2) eff_area = 0.287;
            else if (etaAbs < 2.3) eff_area = 0.525;
            else if (etaAbs < 2.4) eff_area = 0.488;
        }
    }

    // done
    return eff_area;
}



//------------------------------------------
// Initialize baby ntuple variables
//------------------------------------------
void myBabyMaker::InitBabyNtuple() 
{
    /////////////////////////// 
    // Event Information     //
    ///////////////////////////
    
    // Basic Event Information
    run_          = -1;
    ls_           = -1;
    evt_          = 0;
    weight_       = 1.0;
    dataset_      = "";
    filename_     = "";
    is_real_data_ = false;
  
    // Pileup - PUSummaryInfoMaker
    pu_nPUvertices_ = -1;
    pu_nPUtrueint_  = -1.0;
  
    // Pileup - VertexMaker
    evt_nvtxs_ = -1;

    // event level variables
    nFOels_ = 0;
    nFOmus_ = 0;
    ngsfs_  = 0;
    nmus_   = 0;
    nvetomus_ = 0;
    nvetoels_ = 0;

    /////////////////////////// 
    // End Event Information //
    ///////////////////////////


    //////////////////////////// 
    // Lepton Information     //
    ////////////////////////////

    id_       = -1;
    pt_       = -999.;
    eta_      = -999.;
    sceta_    = -999.;
    phi_      = -999.;
    scet_     = -999.;
    pfmet_    = -999.;
    pfmetphi_ = -999.;
    hoe_      = -999.;
    d0_       = -999.;
    dz_       = -999.;
    ip3d_     = -999.;
    d0err_    = -999.;
    dzerr_    = -999.;
    ip3derr_  = -999.;

    lp4_.SetCoordinates(0,0,0,0);     // 4-vector of the lepton
    foel_p4_.SetCoordinates(0,0,0,0); // 4-vector of the highest pt additional FO in the event
    fomu_p4_.SetCoordinates(0,0,0,0); // 4-vector of the highest pt additional FO in the event
    foel_id_   = -999;
    fomu_id_   = -999;
    foel_mass_ = -999.0;
    fomu_mass_ = -999.0;

    iso_                = -999.;
    iso_nps_            = -999.;
    trck_iso_           = -999.;
    ecal_iso_           = -999.;
    ecal_iso_nps_       = -999.;
    hcal_iso_           = -999.;
    pfiso03_            = -999.;
    ch_pfiso03_         = -999.;
    nh_pfiso03_         = -999.;
    em_pfiso03_         = -999.;
    pfiso03_bv_         = -999.;
    ch_pfiso03_bv_      = -999.;
    nh_pfiso03_bv_      = -999.;
    pfiso04_            = -999.;
    ch_pfiso04_         = -999.;
    nh_pfiso04_         = -999.;
    em_pfiso04_         = -999.;
    pfiso04_bv_         = -999.;
    ch_pfiso04_bv_      = -999.;
    nh_pfiso04_bv_      = -999.;
    em_pfiso04_bv_      = -999.;
    radiso_et1p0_       = -999.;
    ch_radiso_et1p0_    = -999.;
    nh_radiso_et1p0_    = -999.;
    em_radiso_et1p0_    = -999.;
    radiso_et0p5_       = -999.;
    ch_radiso_et0p5_    = -999.;
    nh_radiso_et0p5_    = -999.;
    em_radiso_et0p5_    = -999.;
    radiso_et1p0_bv_    = -999.;
    ch_radiso_et1p0_bv_ = -999.;
    nh_radiso_et1p0_bv_ = -999.;
    em_radiso_et1p0_bv_ = -999.;
    radiso_et0p5_bv_    = -999.;
    ch_radiso_et0p5_bv_ = -999.;
    nh_radiso_et0p5_bv_ = -999.;
    em_radiso_et0p5_bv_ = -999.;
    pfpupt03_           = -999.;
    pfpupt04_           = -999.;
    cpfiso03_rho_       = -999.;
    cpfiso04_rho_       = -999.;
    cpfiso03_db_        = -999.;

    closestMuon_      = false;
    el_id_sieie_      = -999.;
    el_id_detain_     = -999.;
    el_id_dphiin_     = -999.;
    el_id_smurfV5_    = false;
    el_id_vbtf80_     = false;
    el_id_vbtf90_     = false;
    convHitPattern_   = false;
    convPartnerTrack_ = false;
    convMIT_          = false;

    el_effarea03_          = -999.;
    el_effarea04_          = -999.;
    mu_effarea03_          = -999.;
    mu_nh_effarea03_       = -999.;
    mu_em_effarea03_       = -999.;
    mu_effarea03_tight_    = -999.;
    mu_nh_effarea03_tight_ = -999.;
    mu_em_effarea03_tight_ = -999.;
    mu_effarea04_          = -999.;
    mu_nh_effarea04_       = -999.;
    mu_em_effarea04_       = -999.;
    mu_effarea04_tight_    = -999.;
    mu_nh_effarea04_tight_ = -999.;
    mu_em_effarea04_tight_ = -999.;


    // Z mass variables
    mz_fo_gsf_       = -999.;
    mz_gsf_iso_      = -999.;
    mz_fo_ctf_       = -999.;
    mz_ctf_iso_      = -999.;
    mupsilon_fo_mu_  = -999.;
    mupsilon_mu_iso_ = -999.;

    d0PV_wwV1_       = -999.;
    dzPV_wwV1_       = -999.;

    mu_isCosmic_ = false;

    mu_ecal_veto_dep_ = -999.;
    mu_hcal_veto_dep_ = -999.;
    mu_nchi2_         = -999.;

    mt_                   = -999;
    pfmt_                 = -999;
    q3_                   = false;
    els_exp_innerlayers_  = 999;
    mcid_                 = 0;
    mcmotherid_           = 0;
      
    // HT
#ifndef __CMS2_SLIM__
    ht_calo_          = -999;           
    ht_calo_L2L3_     = -999;      
#endif
    ht_pf_            = -999;            
    ht_pf_L2L3_       = -999;        
    ht_pf_L1FastL2L3_ = -999;

    // MC truth information
    mc3id_         = -999;
    mc3pt_         = -999.;
    mc3dr_         = -999.;
    mc3p4_.SetCoordinates(0,0,0,0);
    leptonIsFromW_ = -999;

    //////////////////////////// 
    // End Lepton Information //
    ////////////////////////////



    //////////////////////////////////////////////////////
    // Fake Rate Numerator & Denominator Selections     //
    //////////////////////////////////////////////////////

    //////////
    // 2012 //
    //////////

    // SS

    // Electrons
    num_el_ssV7_       = false;
    num_el_ssV7_noIso_ = false;
    v1_el_ssV7_        = false;
    v2_el_ssV7_        = false;
    v3_el_ssV7_        = false;
    
    // Muons
    num_mu_ssV5_       = false;        
    num_mu_ssV5_noIso_ = false;  
    fo_mu_ssV5_        = false;         
    fo_mu_ssV5_noIso_  = false;   

    // TTZ

    // Electrons
    num_el_TTZcuttightv1_       = false;
    num_el_TTZcuttightv1_noIso_ = false;
    fo_el_TTZcuttightv1_        = false;
    fo_el_TTZcuttightv1_noIso_  = false;

    num_el_TTZcutloosev1_       = false;
    num_el_TTZcutloosev1_noIso_ = false;
    fo_el_TTZcutloosev1_        = false;
    fo_el_TTZcutloosev1_noIso_  = false;

    num_el_TTZMVAtightv1_       = false;
    num_el_TTZMVAtightv1_noIso_ = false;
    fo_el_TTZMVAtightv1_        = false;
    fo_el_TTZMVAtightv1_noIso_  = false;

    num_el_TTZMVAloosev1_       = false;
    num_el_TTZMVAloosev1_noIso_ = false;
    fo_el_TTZMVAloosev1_        = false;
    fo_el_TTZMVAloosev1_noIso_  = false;

    // Muons
    num_mu_TTZtightv1_       = false;
    num_mu_TTZtightv1_noIso_ = false;
    fo_mu_TTZtightv1_        = false;
    fo_mu_TTZtightv1_noIso_  = false;

    num_mu_TTZloosev1_       = false;
    num_mu_TTZloosev1_noIso_ = false;
    fo_mu_TTZloosev1_        = false;
    fo_mu_TTZloosev1_noIso_  = false;



    //////////
    // 2011 //
    //////////

    // SS

    // Electrons
    num_el_ssV6_       = false;
    num_el_ssV6_noIso_ = false;
    v1_el_ssV6_        = false;
    v2_el_ssV6_        = false;
    v3_el_ssV6_        = false;

    // Muons
    numNomSSv4_      = false;
    numNomSSv4noIso_ = false; 
    fo_mussV4_04_    = false;
    fo_mussV4_noIso_ = false;

    // WW, HWW

    // Electrons
    num_el_smurfV6_ = false;
    num_el_smurfV6lh_ = false;
    v1_el_smurfV1_  = false;
    v2_el_smurfV1_  = false;
    v3_el_smurfV1_  = false;
    v4_el_smurfV1_  = false;

    // Muons
    num_mu_smurfV6_ = false;
    fo_mu_smurf_04_ = false;
    fo_mu_smurf_10_ = false;


    // OS
    num_el_OSV2_  = false;
    num_mu_OSGV2_ = false;
    num_mu_OSZV2_ = false;
    fo_el_OSV2_   = false;
    fo_mu_OSGV2_  = false;

    // OS
    num_el_OSV3_  = false;
    num_mu_OSGV3_ = false;
    fo_el_OSV3_   = false;
    fo_mu_OSGV3_  = false;
    
    //////////////////////////////////////////////////////
    // End Fake Rate Numerator & Denominator Selections //
    //////////////////////////////////////////////////////

    ///////////////////////  
    // 2012 Triggers     //
    ///////////////////////

    // Electrons
    ele17_CaloIdL_CaloIsoVL_vstar_                                     = 0;
    ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                    = 0;
    ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_              = 0;
    ele8_CaloIdL_CaloIsoVL_vstar_                                      = 0;
    ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                     = 0;
    ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_               = 0;
    ele8_CaloIdT_TrkIdVL_vstar_                                        = 0;
    ele8_CaloIdT_TrkIdVL_Jet30_vstar_                                  = 0;
    ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_       = 0;
    ele27_WP80_vstar_                                                  = 0;

    ele17_CaloIdL_CaloIsoVL_version_                                   = - 1;
    ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_                  = - 1;
    ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_            = - 1;
    ele8_CaloIdL_CaloIsoVL_version_                                    = - 1;
    ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_                   = - 1;
    ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_             = - 1;
    ele8_CaloIdT_TrkIdVL_version_                                      = - 1;
    ele8_CaloIdT_TrkIdVL_Jet30_version_                                = - 1;
    ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_version_     = - 1;
    ele27_WP80_version_                                                = - 1;

    dr_ele17_CaloIdL_CaloIsoVL_vstar_                                  = 99.0;
    dr_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                 = 99.0;
    dr_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_           = 99.0;
    dr_ele8_CaloIdL_CaloIsoVL_vstar_                                   = 99.0;
    dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                  = 99.0;
    dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_            = 99.0;
    dr_ele8_CaloIdT_TrkIdVL_vstar_                                     = 99.0;
    dr_ele8_CaloIdT_TrkIdVL_Jet30_vstar_                               = 99.0;
    dr_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_    = 99.0;
    dr_ele27_WP80_vstar_                                               = 99.0;

    hltps_ele17_CaloIdL_CaloIsoVL_vstar_                               = - 1;
    hltps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_              = - 1;
    hltps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_        = - 1;
    hltps_ele8_CaloIdL_CaloIsoVL_vstar_                                = - 1;
    hltps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_               = - 1;
    hltps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_         = - 1;
    hltps_ele8_CaloIdT_TrkIdVL_vstar_                                  = - 1;
    hltps_ele8_CaloIdT_TrkIdVL_Jet30_vstar_                            = - 1;
    hltps_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_ = - 1;
    hltps_ele27_WP80_vstar_                                            = - 1;

#ifndef __CMS2_SLIM__
    l1ps_ele17_CaloIdL_CaloIsoVL_vstar_                                = - 1;
    l1ps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_               = - 1;
    l1ps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_         = - 1;
    l1ps_ele8_CaloIdL_CaloIsoVL_vstar_                                 = - 1;
    l1ps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                = - 1;
    l1ps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_          = - 1;
    l1ps_ele8_CaloIdT_TrkIdVL_vstar_                                   = - 1;
    l1ps_ele8_CaloIdT_TrkIdVL_Jet30_vstar_                             = - 1;
    l1ps_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_  = - 1;
    l1ps_ele27_WP80_vstar_                                             = - 1;
#endif

    // Muons
    mu5_vstar_                  = 0;
    mu8_vstar_                  = 0;
    mu12_vstar_                 = 0;
    mu17_vstar_                 = 0;
    mu15_eta2p1_vstar_          = 0;
    mu24_eta2p1_vstar_          = 0;
    mu30_eta2p1_vstar_          = 0;
    isoMu20_eta2p1_vstar_       = 0;
    isoMu24_eta2p1_vstar_       = 0;
    isoMu30_eta2p1_vstar_       = 0;
    relIso1p0Mu17_vstar_        = 0;
    relIso1p0Mu20_vstar_        = 0;
    relIso1p0Mu5_vstar_         = 0;

    mu5_version_                = -1;
    mu8_version_                = -1;
    mu12_version_               = -1;
    mu17_version_               = -1;
    mu15_eta2p1_version_        = -1;
    mu24_eta2p1_version_        = -1;
    mu30_eta2p1_version_        = -1;
    isoMu20_eta2p1_version_     = -1;
    isoMu24_eta2p1_version_     = -1;
    isoMu30_eta2p1_version_     = -1;
    relIso1p0Mu17_version_      = -1;
    relIso1p0Mu20_version_      = -1;
    relIso1p0Mu5_version_       = -1;

    dr_mu8_vstar_               = 99.0;
    dr_mu5_vstar_               = 99.0;
    dr_mu12_vstar_              = 99.0;
    dr_mu17_vstar_              = 99.0;
    dr_mu15_eta2p1_vstar_       = 99.0;
    dr_mu24_eta2p1_vstar_       = 99.0;
    dr_mu30_eta2p1_vstar_       = 99.0;
    dr_isoMu20_eta2p1_vstar_    = 99.0;
    dr_isoMu24_eta2p1_vstar_    = 99.0;
    dr_isoMu30_eta2p1_vstar_    = 99.0;
    dr_relIso1p0Mu17_vstar_     = 99.0;
    dr_relIso1p0Mu20_vstar_     = 99.0;
    dr_relIso1p0Mu5_vstar_      = 99.0;

    hltps_mu8_vstar_            = -1;
    hltps_mu5_vstar_            = -1;
    hltps_mu12_vstar_           = -1;
    hltps_mu17_vstar_           = -1;
    hltps_mu15_eta2p1_vstar_    = -1;
    hltps_mu24_eta2p1_vstar_    = -1;
    hltps_mu30_eta2p1_vstar_    = -1;
    hltps_isoMu20_eta2p1_vstar_ = -1;
    hltps_isoMu24_eta2p1_vstar_ = -1;
    hltps_isoMu30_eta2p1_vstar_ = -1;
    hltps_relIso1p0Mu17_vstar_  = -1;
    hltps_relIso1p0Mu20_vstar_  = -1;
    hltps_relIso1p0Mu5_vstar_   = -1;

#ifndef __CMS2_SLIM__
    l1ps_mu8_vstar_             = -1;
    l1ps_mu5_vstar_             = -1;
    l1ps_mu12_vstar_            = -1;
    l1ps_mu17_vstar_            = -1;
    l1ps_mu15_eta2p1_vstar_     = -1;
    l1ps_mu24_eta2p1_vstar_     = -1;
    l1ps_mu30_eta2p1_vstar_     = -1;
    l1ps_isoMu20_eta2p1_vstar_  = -1;
    l1ps_isoMu24_eta2p1_vstar_  = -1;
    l1ps_isoMu30_eta2p1_vstar_  = -1;
    l1ps_relIso1p0Mu17_vstar_   = -1;
    l1ps_relIso1p0Mu5_vstar_    = -1;
#endif

    ///////////////////////  
    // End 2012 Triggers //
    ///////////////////////


   
  
    ///////////////////////  
    // 2011 Triggers     //
    ///////////////////////

    // Electrons
    ele8_vstar_                                             = 0;
    ele8_CaloIdL_TrkIdVL_vstar_                             = 0;
    ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                     = 0;
    ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_          = 0;
    photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_    = 0;

    ele8_version_                                           = -1;
    ele8_CaloIdL_TrkIdVL_version_                           = -1;
    ele8_CaloIdL_CaloIsoVL_Jet40_version_                   = -1;
    ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_version_        = -1;
    photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_version_  = -1;

    dr_ele8_vstar_                                          = 99.0; 
    dr_ele8_CaloIdL_TrkIdVL_vstar_                          = 99.0; 
    dr_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                  = 99.0; 
    dr_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_       = 99.0;
    dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_ = 99.0; 

    hltps_ele8_vstar_                                          = -1; 
    hltps_ele8_CaloIdL_TrkIdVL_vstar_                          = -1; 
    hltps_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                  = -1; 
    hltps_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_       = -1;
    hltps_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_ = -1; 

    // Muons
    mu3_vstar_          = 0;
    mu15_vstar_         = 0;  
    mu20_vstar_         = 0;  
    mu24_vstar_         = 0;  
    mu30_vstar_         = 0;  
    mu8_Jet40_vstar_    = 0;    

    mu3_version_        = -1;
    mu15_version_       = -1;  
    mu20_version_       = -1;  
    mu24_version_       = -1;  
    mu30_version_       = -1;  
    mu8_Jet40_version_  = -1;    

    dr_mu3_vstar_       = 99.0;
    dr_mu15_vstar_      = 99.0; 
    dr_mu20_vstar_      = 99.0; 
    dr_mu24_vstar_      = 99.0; 
    dr_mu30_vstar_      = 99.0;
    dr_mu8_Jet40_vstar_ = 99.0;

    hltps_mu3_vstar_       = -1;
    hltps_mu15_vstar_      = -1; 
    hltps_mu20_vstar_      = -1; 
    hltps_mu24_vstar_      = -1; 
    hltps_mu30_vstar_      = -1;
    hltps_mu8_Jet40_vstar_ = -1;

    ///////////////////////  
    // End 2011 Triggers //
    ///////////////////////

  
    //////////////
    // Jets     //
    //////////////

#ifndef __CMS2_SLIM__
    // Calo Jets
    ptj1_       = 0.;
    ptj1_b2b_   = -999.;
    dphij1_b2b_ = -999.;
    nj1_        = 0;
#endif

    // PF Jets
    ptpfj1_       = 0.;
    ptpfj1_b2b_   = -999.;
    dphipfj1_b2b_ = -999.;
    npfj1_        = 0;

    // PF L2L3 Corrected jets
    ptpfcj1_        = 0.;
    ptpfcj1_b2b_    = -999.;
    dphipfcj1_b2b_  = -999.;
    npfcj1_         = 0;
    btagpfc_        = false;

    // PF L1FastL2L3 Corrected jets
    emfpfcL1Fj1_       = -999.;
    ptpfcL1Fj1_        = 0.;
    dphipfcL1Fj1_      = -999.;
    ptpfcL1Fj1_b2b_    = -999.;
    dphipfcL1Fj1_b2b_  = -999.;
    npfcL1Fj1_         = 0;
    npfc30L1Fj1_       = 0;
    npfc40L1Fj1_       = 0;
    nbpfc40L1Fj1_      = 0;
    btagpfcL1F_        = false;
    npfc50L1Fj1_eth_   = 0;
    npfc65L1Fj1_eth_   = 0;

    // PF L1FastL2L3Residual Corrected jets
    emfpfcL1Fj1res_       = -999.;
    ptpfcL1Fj1res_        = 0.;
    dphipfcL1Fj1res_      = -999.;
    ptpfcL1Fj1res_b2b_    = -999.;
    dphipfcL1Fj1res_b2b_  = -999.;
    npfcL1Fj1res_         = 0;
    npfc30L1Fj1res_       = 0;
    npfc40L1Fj1res_       = 0;
    nbpfc40L1Fj1res_      = 0;
    btagpfcL1Fres_        = false;
    npfc50L1Fj1res_eth_   = 0;
    npfc65L1Fj1res_eth_   = 0;

    // btag PF L1FastL2L3 Corrected jets
    ptbtagpfcL1Fj1_        = 0.;
    dphibtagpfcL1Fj1_      = -999.;
    
    // btag PF L1FastL2L3Residual Corrected jets
    ptbtagpfcL1Fj1res_        = 0.;
    dphibtagpfcL1Fj1res_      = -999.;
    
#ifndef __CMS2_SLIM__
    // CALO Jet info
    nbjet_  = 0;
    dRbNear_ = 99.;
    dRbFar_ = -99.;
#endif
    nbpfcjet_  = 0;
    dRbpfcNear_ = 99.;
    dRbpfcFar_ = -99.;

    rho_ = -999.;

    //////////////
    // End Jets //
    //////////////

}

// Book the baby ntuple
void myBabyMaker::MakeBabyNtuple(const char *babyFilename)
{
    babyFile_ = TFile::Open(Form("%s", babyFilename), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");

    /////////////////////////// 
    // Event Information     //
    ///////////////////////////

    babyTree_->Branch("run"            , &run_                );
    babyTree_->Branch("ls"             , &ls_                 );
    babyTree_->Branch("evt"            , &evt_                );
    babyTree_->Branch("weight"         , &weight_             );
    babyTree_->Branch("is_real_data"   , &is_real_data_       );
    babyTree_->Branch("dataset"        , "TString", &dataset_ );
    babyTree_->Branch("filename"       , "TString", &filename_);
  
    // Pileup
    babyTree_->Branch("pu_nPUvertices" , &pu_nPUvertices_ );
    babyTree_->Branch("pu_nPUtrueint"  , &pu_nPUtrueint_ );
    babyTree_->Branch("evt_nvtxs"      , &evt_nvtxs_      );

    // event level variables
    babyTree_->Branch("nFOels"         , &nFOels_         );
    babyTree_->Branch("nFOmus"         , &nFOmus_         );
    babyTree_->Branch("ngsfs"          , &ngsfs_          );
    babyTree_->Branch("nmus"           , &nmus_           );
    babyTree_->Branch("nvetoels"       , &nvetoels_       );
    babyTree_->Branch("nvetomus"       , &nvetomus_       );

    /////////////////////////// 
    // End Event Information //
    ///////////////////////////
        
        
        
    //////////////////////////// 
    // Lepton Information     //
    ////////////////////////////

    babyTree_->Branch("lp4"    , "LorentzVector", &lp4_    ); 
    babyTree_->Branch("mc3p4"  , "LorentzVector", &mc3p4_  ); 
    babyTree_->Branch("foel_p4", "LorentzVector", &foel_p4_); 
    babyTree_->Branch("fomu_p4", "LorentzVector", &fomu_p4_); 

    babyTree_->Branch("foel_id"               , &foel_id_               ); 
    babyTree_->Branch("fomu_id"               , &fomu_id_               ); 
    babyTree_->Branch("foel_mass"             , &foel_mass_             ); 
    babyTree_->Branch("fomu_mass"             , &fomu_mass_             ); 
    babyTree_->Branch("pt"                    , &pt_                    ); 
    babyTree_->Branch("eta"                   , &eta_                   ); 
    babyTree_->Branch("sceta"                 , &sceta_                 ); 
    babyTree_->Branch("phi"                   , &phi_                   ); 
    babyTree_->Branch("scet"                  , &scet_                  ); 
    babyTree_->Branch("d0"                    , &d0_                    ); 
    babyTree_->Branch("dz"                    , &dz_                    ); 
    babyTree_->Branch("ip3d"                  , &ip3d_                  ); 
    babyTree_->Branch("d0err"                 , &d0err_                 ); 
    babyTree_->Branch("dzerr"                 , &dzerr_                 ); 
    babyTree_->Branch("ip3derr"               , &ip3derr_               ); 
    babyTree_->Branch("hoe"                   , &hoe_                   ); 
    babyTree_->Branch("pfmet"                 , &pfmet_                 ); 
    babyTree_->Branch("pfmetphi"              , &pfmetphi_              ); 
    babyTree_->Branch("iso"                   , &iso_                   ); 
    babyTree_->Branch("iso_nps"               , &iso_nps_               ); 
    babyTree_->Branch("trck_iso"              , &trck_iso_              ); 
    babyTree_->Branch("ecal_iso"              , &ecal_iso_              ); 
    babyTree_->Branch("ecal_iso_nps"          , &ecal_iso_nps_          ); 
    babyTree_->Branch("hcal_iso"              , &hcal_iso_              ); 
    babyTree_->Branch("pfiso03"               , &pfiso03_               ); 
    babyTree_->Branch("ch_pfiso03"            , &ch_pfiso03_            ); 
    babyTree_->Branch("nh_pfiso03"            , &nh_pfiso03_            ); 
    babyTree_->Branch("em_pfiso03"            , &em_pfiso03_            ); 
    //babyTree_->Branch("pfiso03_bv"            , &pfiso03_bv_            ); 
    //babyTree_->Branch("ch_pfiso03_bv"         , &ch_pfiso03_bv_         ); 
    //babyTree_->Branch("nh_pfiso03_bv"         , &nh_pfiso03_bv_         ); 
    //babyTree_->Branch("em_pfiso03_bv"         , &em_pfiso03_bv_         ); 
    babyTree_->Branch("pfiso04"               , &pfiso04_               ); 
    babyTree_->Branch("ch_pfiso04"            , &ch_pfiso04_            ); 
    babyTree_->Branch("nh_pfiso04"            , &nh_pfiso04_            ); 
    babyTree_->Branch("em_pfiso04"            , &em_pfiso04_            ); 
    //babyTree_->Branch("pfiso04_bv"            , &pfiso04_bv_            ); 
    //babyTree_->Branch("ch_pfiso04_bv"         , &ch_pfiso04_bv_         ); 
    //babyTree_->Branch("nh_pfiso04_bv"         , &nh_pfiso04_bv_         ); 
    //babyTree_->Branch("em_pfiso04_bv"         , &em_pfiso04_bv_         ); 
    //babyTree_->Branch("radiso_et1p0"          , &radiso_et1p0_          ); 
    //babyTree_->Branch("ch_radiso_et1p0"       , &ch_radiso_et1p0_       ); 
    //babyTree_->Branch("nh_radiso_et1p0"       , &nh_radiso_et1p0_       ); 
    //babyTree_->Branch("em_radiso_et1p0"       , &em_radiso_et1p0_       ); 
    //babyTree_->Branch("radiso_et0p5"          , &radiso_et0p5_          ); 
    //babyTree_->Branch("ch_radiso_et0p5"       , &ch_radiso_et0p5_       ); 
    //babyTree_->Branch("nh_radiso_et0p5"       , &nh_radiso_et0p5_       ); 
    //babyTree_->Branch("em_radiso_et0p5"       , &em_radiso_et0p5_       ); 
    //babyTree_->Branch("radiso_et1p0_bv"       , &radiso_et1p0_bv_       ); 
    //babyTree_->Branch("ch_radiso_et1p0_bv"    , &ch_radiso_et1p0_bv_    ); 
    //babyTree_->Branch("nh_radiso_et1p0_bv"    , &nh_radiso_et1p0_bv_    ); 
    //babyTree_->Branch("em_radiso_et1p0_bv"    , &em_radiso_et1p0_bv_    ); 
    //babyTree_->Branch("radiso_et0p5_bv"       , &radiso_et0p5_bv_       ); 
    //babyTree_->Branch("ch_radiso_et0p5_bv"    , &ch_radiso_et0p5_bv_    ); 
    //babyTree_->Branch("nh_radiso_et0p5_bv"    , &nh_radiso_et0p5_bv_    ); 
    //babyTree_->Branch("em_radiso_et0p5_bv"    , &em_radiso_et0p5_bv_    ); 
    babyTree_->Branch("pfpupt03"              , &pfpupt03_              ); 
    babyTree_->Branch("pfpupt04"              , &pfpupt04_              ); 
    babyTree_->Branch("cpfiso03_rho"          , &cpfiso03_rho_          ); 
    babyTree_->Branch("cpfiso04_rho"          , &cpfiso04_rho_          ); 
    babyTree_->Branch("cpfiso03_db"           , &cpfiso03_db_           ); 
    babyTree_->Branch("id"                    , &id_                    ); 
    babyTree_->Branch("closestMuon"           , &closestMuon_           ); 
    babyTree_->Branch("el_id_sieie"           , &el_id_sieie_           );
    babyTree_->Branch("el_id_detain"          , &el_id_detain_          );
    babyTree_->Branch("el_id_dphiin"          , &el_id_dphiin_          );
    babyTree_->Branch("el_id_smurfV5"         , &el_id_smurfV5_         ); 
    babyTree_->Branch("el_id_vbtf80"          , &el_id_vbtf80_          ); 
    babyTree_->Branch("el_id_vbtf90"          , &el_id_vbtf90_          ); 
    babyTree_->Branch("el_effarea03"          , &el_effarea03_          ); 
    babyTree_->Branch("el_effarea04"          , &el_effarea04_          ); 
    babyTree_->Branch("mu_effarea03"          , &mu_effarea03_          ); 
    babyTree_->Branch("mu_nh_effarea03"       , &mu_nh_effarea03_       ); 
    babyTree_->Branch("mu_em_effarea03"       , &mu_em_effarea03_       ); 
    babyTree_->Branch("mu_effarea03_tight"    , &mu_effarea03_tight_    ); 
    babyTree_->Branch("mu_nh_effarea03_tight" , &mu_nh_effarea03_tight_ ); 
    babyTree_->Branch("mu_em_effarea03_tight" , &mu_em_effarea03_tight_ ); 
    babyTree_->Branch("mu_effarea04"          , &mu_effarea04_          ); 
    babyTree_->Branch("mu_nh_effarea04"       , &mu_nh_effarea04_       ); 
    babyTree_->Branch("mu_em_effarea04"       , &mu_em_effarea04_       ); 
    babyTree_->Branch("mu_effarea04_tight"    , &mu_effarea04_tight_    ); 
    babyTree_->Branch("mu_nh_effarea04_tight" , &mu_nh_effarea04_tight_ ); 
    babyTree_->Branch("mu_em_effarea04_tight" , &mu_em_effarea04_tight_ ); 
    babyTree_->Branch("conv0MissHits"         , &conv0MissHits_         ); 
    babyTree_->Branch("convHitPattern"        , &convHitPattern_        ); 
    babyTree_->Branch("convPartnerTrack"      , &convPartnerTrack_      ); 
    babyTree_->Branch("convMIT"               , &convMIT_               ); 
    babyTree_->Branch("mt"                    , &mt_                    ); 
    babyTree_->Branch("pfmt"                  , &pfmt_                  ); 
    babyTree_->Branch("q3"                    , &q3_                    ); 
    babyTree_->Branch("els_exp_innerlayers"   , &els_exp_innerlayers_   ); 
    babyTree_->Branch("d0PV_wwV1"             , &d0PV_wwV1_             ); 
    babyTree_->Branch("dzPV_wwV1"             , &dzPV_wwV1_             ); 
#ifndef __CMS2_SLIM__
    babyTree_->Branch("ht_calo"               , &ht_calo_               ); 
    babyTree_->Branch("ht_calo_L2L3"          , &ht_calo_L2L3_          ); 
#endif
    babyTree_->Branch("ht_pf"                 , &ht_pf_                 ); 
    babyTree_->Branch("ht_pf_L2L3"            , &ht_pf_L2L3_            ); 
    babyTree_->Branch("ht_pf_L1FastL2L3"      , &ht_pf_L1FastL2L3_      ); 
    babyTree_->Branch("mcid"                  , &mcid_                  ); 
    babyTree_->Branch("mcmotherid"            , &mcmotherid_            ); 
    babyTree_->Branch("mc3id"                 , &mc3id_                 ); 
    babyTree_->Branch("mc3pt"                 , &mc3pt_                 ); 
    babyTree_->Branch("mc3dr"                 , &mc3dr_                 ); 
    babyTree_->Branch("leptonIsFromW"         , &leptonIsFromW_         ); 
    babyTree_->Branch("mu_isCosmic"           , &mu_isCosmic_           ); 
    babyTree_->Branch("mu_ecal_veto_dep"      , &mu_ecal_veto_dep_      );
    babyTree_->Branch("mu_hcal_veto_dep"      , &mu_hcal_veto_dep_      );
    babyTree_->Branch("mu_nchi2"              , &mu_nchi2_              );

    // Z mass variables
    babyTree_->Branch("mz_fo_gsf"      , &mz_fo_gsf_      );
    babyTree_->Branch("mz_gsf_iso"     , &mz_gsf_iso_     );
    babyTree_->Branch("mz_fo_ctf"      , &mz_fo_ctf_      );
    babyTree_->Branch("mz_ctf_iso"     , &mz_ctf_iso_     );
    babyTree_->Branch("mupsilon_fo_mu" , &mupsilon_fo_mu_ );
    babyTree_->Branch("mupsilon_mu_iso", &mupsilon_mu_iso_);

    //////////////////////////// 
    // End Lepton Information //
    ////////////////////////////

    //////////////////////////////////////////////////////
    // Fake Rate Numerator & Denominator Selections     //
    //////////////////////////////////////////////////////

    //////////
    // 2012 //
    //////////

    // SS
    // Electrons
    babyTree_->Branch("num_el_ssV7"       , &num_el_ssV7_       );
    babyTree_->Branch("num_el_ssV7_noIso" , &num_el_ssV7_noIso_ );
    babyTree_->Branch("v1_el_ssV7"        , &v1_el_ssV7_        );
    babyTree_->Branch("v2_el_ssV7"        , &v2_el_ssV7_        );
    babyTree_->Branch("v3_el_ssV7"        , &v3_el_ssV7_        );
    
    // Muons
    babyTree_->Branch("num_mu_ssV5"       , &num_mu_ssV5_       );
    babyTree_->Branch("num_mu_ssV5_noIso" , &num_mu_ssV5_noIso_ );
    babyTree_->Branch("fo_mu_ssV5"        , &fo_mu_ssV5_        );
    babyTree_->Branch("fo_mu_ssV5_noIso"  , &fo_mu_ssV5_noIso_  );

    // TTZ
    // Electrons
    babyTree_->Branch("num_el_TTZcuttightv1",       &num_el_TTZcuttightv1_);
    babyTree_->Branch("num_el_TTZcuttightv1_noIso", &num_el_TTZcuttightv1_noIso_);
    babyTree_->Branch("fo_el_TTZcuttightv1",        &fo_el_TTZcuttightv1_);
    babyTree_->Branch("fo_el_TTZcuttightv1_noIso",  &fo_el_TTZcuttightv1_noIso_);

    babyTree_->Branch("num_el_TTZcutloosev1",       &num_el_TTZcutloosev1_);
    babyTree_->Branch("num_el_TTZcutloosev1_noIso", &num_el_TTZcutloosev1_noIso_);
    babyTree_->Branch("fo_el_TTZcutloosev1",        &fo_el_TTZcutloosev1_);
    babyTree_->Branch("fo_el_TTZcutloosev1_noIso",  &fo_el_TTZcutloosev1_noIso_);

    babyTree_->Branch("num_el_TTZMVAtightv1",       &num_el_TTZMVAtightv1_);
    babyTree_->Branch("num_el_TTZMVAtightv1_noIso", &num_el_TTZMVAtightv1_noIso_);
    babyTree_->Branch("fo_el_TTZMVAtightv1",        &fo_el_TTZMVAtightv1_);
    babyTree_->Branch("fo_el_TTZMVAtightv1_noIso",  &fo_el_TTZMVAtightv1_noIso_);

    babyTree_->Branch("num_el_TTZMVAloosev1",       &num_el_TTZMVAloosev1_);
    babyTree_->Branch("num_el_TTZMVAloosev1_noIso", &num_el_TTZMVAloosev1_noIso_);
    babyTree_->Branch("fo_el_TTZMVAloosev1",        &fo_el_TTZMVAloosev1_);
    babyTree_->Branch("fo_el_TTZMVAloosev1_noIso",  &fo_el_TTZMVAloosev1_noIso_);

    // Muons
    babyTree_->Branch("num_mu_TTZtightv1",       &num_mu_TTZtightv1_);
    babyTree_->Branch("num_mu_TTZtightv1_noIso", &num_mu_TTZtightv1_noIso_);
    babyTree_->Branch("fo_mu_TTZtightv1",        &fo_mu_TTZtightv1_);
    babyTree_->Branch("fo_mu_TTZtightv1_noIso",  &fo_mu_TTZtightv1_noIso_);

    babyTree_->Branch("num_mu_TTZloosev1",       &num_mu_TTZloosev1_);
    babyTree_->Branch("num_mu_TTZloosev1_noIso", &num_mu_TTZloosev1_noIso_);
    babyTree_->Branch("fo_mu_TTZloosev1",        &fo_mu_TTZloosev1_);
    babyTree_->Branch("fo_mu_TTZloosev1_noIso",  &fo_mu_TTZloosev1_noIso_);



    //////////
    // 2011 //
    //////////

    // SS

    // Electrons
    babyTree_->Branch("num_el_ssV6"       , &num_el_ssV6_       );
    babyTree_->Branch("num_el_ssV6_noIso" , &num_el_ssV6_noIso_ );
    babyTree_->Branch("v1_el_ssV6"        , &v1_el_ssV6_        );
    babyTree_->Branch("v2_el_ssV6"        , &v2_el_ssV6_        );
    babyTree_->Branch("v3_el_ssV6"        , &v3_el_ssV6_        );

    // Muons
    babyTree_->Branch("numNomSSv4"      , &numNomSSv4_      );
    babyTree_->Branch("numNomSSv4noIso" , &numNomSSv4noIso_ );
    babyTree_->Branch("fo_mussV4_04"    , &fo_mussV4_04_    );
    babyTree_->Branch("fo_mussV4_noIso" , &fo_mussV4_noIso_ );

    // WW, HWW

    // Electrons
    babyTree_->Branch("num_el_smurfV6"   , &num_el_smurfV6_   );
    babyTree_->Branch("num_el_smurfV6lh" , &num_el_smurfV6lh_ );
    babyTree_->Branch("v1_el_smurfV1"    , &v1_el_smurfV1_    );
    babyTree_->Branch("v2_el_smurfV1"    , &v2_el_smurfV1_    );
    babyTree_->Branch("v3_el_smurfV1"    , &v3_el_smurfV1_    );
    babyTree_->Branch("v4_el_smurfV1"    , &v4_el_smurfV1_    );

    // Muons
    babyTree_->Branch("num_mu_smurfV6",  &num_mu_smurfV6_ );
    babyTree_->Branch("fo_mu_smurf_04",  &fo_mu_smurf_04_ );
    babyTree_->Branch("fo_mu_smurf_10",  &fo_mu_smurf_10_ );
  
    // OS
    babyTree_->Branch("num_el_OSV2"  , &num_el_OSV2_  );
    babyTree_->Branch("num_mu_OSGV2" , &num_mu_OSGV2_ );
    babyTree_->Branch("num_mu_OSZV2" , &num_mu_OSZV2_ );
    babyTree_->Branch("fo_el_OSV2"   , &fo_el_OSV2_   );
    babyTree_->Branch("fo_mu_OSGV2"  , &fo_mu_OSGV2_  );

    babyTree_->Branch("num_el_OSV3"  , &num_el_OSV3_  );
    babyTree_->Branch("num_mu_OSGV3" , &num_mu_OSGV3_ );
    babyTree_->Branch("fo_el_OSV3"   , &fo_el_OSV3_   );
    babyTree_->Branch("fo_mu_OSGV3"  , &fo_mu_OSGV3_  );


    //////////////////////////////////////////////////////
    // End Fake Rate Numerator & Denominator Selections //
    //////////////////////////////////////////////////////

    ///////////////////////  
    // Triggers          //
    ///////////////////////

    // Electrons
    babyTree_->Branch("ele8_vstar"                                                        , &ele8_vstar_                                                        );
    babyTree_->Branch("ele8_CaloIdL_TrkIdVL_vstar"                                        , &ele8_CaloIdL_TrkIdVL_vstar_                                        );
    babyTree_->Branch("ele8_CaloIdL_CaloIsoVL_Jet40_vstar"                                , &ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                                );
    babyTree_->Branch("ele8_CaloIdL_CaloIsoVL_vstar"                                      , &ele8_CaloIdL_CaloIsoVL_vstar_                                      );
    babyTree_->Branch("ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar"                     , &ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_                     );
    babyTree_->Branch("ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar"                     , &ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                     );
    babyTree_->Branch("ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar"               , &ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_               );
    babyTree_->Branch("ele8_CaloIdT_TrkIdVL_vstar"                                        , &ele8_CaloIdT_TrkIdVL_vstar_                                        );
    babyTree_->Branch("ele8_CaloIdT_TrkIdVL_Jet30_vstar"                                  , &ele8_CaloIdT_TrkIdVL_Jet30_vstar_                                  );
    babyTree_->Branch("ele17_CaloIdL_CaloIsoVL_vstar"                                     , &ele17_CaloIdL_CaloIsoVL_vstar_                                     );
    babyTree_->Branch("ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar"                    , &ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                    );
    babyTree_->Branch("ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar"              , &ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_              );
    babyTree_->Branch("photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar"               , &photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_               );
    babyTree_->Branch("ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar"       , &ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_       );
    babyTree_->Branch("ele27_WP80_vstar"                                                  , &ele27_WP80_vstar_                                                  );

    babyTree_->Branch("ele8_version"                                                      , &ele8_version_                                                      );
    babyTree_->Branch("ele8_CaloIdL_TrkIdVL_version"                                      , &ele8_CaloIdL_TrkIdVL_version_                                      );
    babyTree_->Branch("ele8_CaloIdL_CaloIsoVL_Jet40_version"                              , &ele8_CaloIdL_CaloIsoVL_Jet40_version_                              );
    babyTree_->Branch("ele8_CaloIdL_CaloIsoVL_version"                                    , &ele8_CaloIdL_CaloIsoVL_version_                                    );
    babyTree_->Branch("ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_version"                   , &ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_version_                   );
    babyTree_->Branch("ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version"                   , &ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_                   );
    babyTree_->Branch("ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version"             , &ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_             );
    babyTree_->Branch("ele8_CaloIdT_TrkIdVL_version"                                      , &ele8_CaloIdT_TrkIdVL_version_                                      );
    babyTree_->Branch("ele8_CaloIdT_TrkIdVL_Jet30_version"                                , &ele8_CaloIdT_TrkIdVL_Jet30_version_                                );
    babyTree_->Branch("ele17_CaloIdL_CaloIsoVL_version"                                   , &ele17_CaloIdL_CaloIsoVL_version_                                   );
    babyTree_->Branch("ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version"                  , &ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_                  );
    babyTree_->Branch("ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version"            , &ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_            );
    babyTree_->Branch("photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_version"             , &photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_version_             );
    babyTree_->Branch("ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_version"     , &ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_version_     );
    babyTree_->Branch("ele27_WP80_version"                                                , &ele27_WP80_version_                                                );

    babyTree_->Branch("dr_ele8_vstar"                                                     , &dr_ele8_vstar_                                                     );
    babyTree_->Branch("dr_ele8_CaloIdL_TrkIdVL_vstar"                                     , &dr_ele8_CaloIdL_TrkIdVL_vstar_                                     );
    babyTree_->Branch("dr_ele8_CaloIdL_CaloIsoVL_Jet40_vstar"                             , &dr_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                             );
    babyTree_->Branch("dr_ele8_CaloIdL_CaloIsoVL_vstar"                                   , &dr_ele8_CaloIdL_CaloIsoVL_vstar_                                   );
    babyTree_->Branch("dr_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar"                  , &dr_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_                  );
    babyTree_->Branch("dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar"                  , &dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                  );
    babyTree_->Branch("dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar"            , &dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_            );
    babyTree_->Branch("dr_ele8_CaloIdT_TrkIdVL_vstar"                                     , &dr_ele8_CaloIdT_TrkIdVL_vstar_                                     );
    babyTree_->Branch("dr_ele8_CaloIdT_TrkIdVL_Jet30_vstar"                               , &dr_ele8_CaloIdT_TrkIdVL_Jet30_vstar_                               );
    babyTree_->Branch("dr_ele17_CaloIdL_CaloIsoVL_vstar"                                  , &dr_ele17_CaloIdL_CaloIsoVL_vstar_                                  );
    babyTree_->Branch("dr_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar"                 , &dr_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                 );
    babyTree_->Branch("dr_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar"           , &dr_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_           );
    babyTree_->Branch("dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar"            , &dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_            );
    babyTree_->Branch("dr_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar"    , &dr_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_    );
    babyTree_->Branch("dr_ele27_WP80_vstar"                                               , &dr_ele27_WP80_vstar_                                               );

    babyTree_->Branch("hltps_ele8_vstar"                                                  , &hltps_ele8_vstar_                                                  );
    babyTree_->Branch("hltps_ele8_CaloIdL_TrkIdVL_vstar"                                  , &hltps_ele8_CaloIdL_TrkIdVL_vstar_                                  );
    babyTree_->Branch("hltps_ele8_CaloIdL_CaloIsoVL_Jet40_vstar"                          , &hltps_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                          );
    babyTree_->Branch("hltps_ele8_CaloIdL_CaloIsoVL_vstar"                                , &hltps_ele8_CaloIdL_CaloIsoVL_vstar_                                );
    babyTree_->Branch("hltps_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar"               , &hltps_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_               );
    babyTree_->Branch("hltps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar"               , &hltps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_               );
    babyTree_->Branch("hltps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar"         , &hltps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_         );
    babyTree_->Branch("hltps_ele8_CaloIdT_TrkIdVL_vstar"                                  , &hltps_ele8_CaloIdT_TrkIdVL_vstar_                                  );
    babyTree_->Branch("hltps_ele8_CaloIdT_TrkIdVL_Jet30_vstar"                            , &hltps_ele8_CaloIdT_TrkIdVL_Jet30_vstar_                            );
    babyTree_->Branch("hltps_ele17_CaloIdL_CaloIsoVL_vstar"                               , &hltps_ele17_CaloIdL_CaloIsoVL_vstar_                               );
    babyTree_->Branch("hltps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar"              , &hltps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_              );
    babyTree_->Branch("hltps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar"        , &hltps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_        );
    babyTree_->Branch("hltps_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar"         , &hltps_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_         );
    babyTree_->Branch("hltps_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar" , &hltps_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_ );
    babyTree_->Branch("hltps_ele27_WP80_vstar"                                            , &hltps_ele27_WP80_vstar_                                            );

#ifndef __CMS2_SLIM__
    babyTree_->Branch("l1ps_ele8_CaloIdL_CaloIsoVL_vstar"                                 , &l1ps_ele8_CaloIdL_CaloIsoVL_vstar_                                 );
    babyTree_->Branch("l1ps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar"                , &l1ps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                );
    babyTree_->Branch("l1ps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar"          , &l1ps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_          );
    babyTree_->Branch("l1ps_ele8_CaloIdT_TrkIdVL_vstar"                                   , &l1ps_ele8_CaloIdT_TrkIdVL_vstar_                                   );
    babyTree_->Branch("l1ps_ele8_CaloIdT_TrkIdVL_Jet30_vstar"                             , &l1ps_ele8_CaloIdT_TrkIdVL_Jet30_vstar_                             );
    babyTree_->Branch("l1ps_ele17_CaloIdL_CaloIsoVL_vstar"                                , &l1ps_ele17_CaloIdL_CaloIsoVL_vstar_                                );
    babyTree_->Branch("l1ps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar"               , &l1ps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_               );
    babyTree_->Branch("l1ps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar"         , &l1ps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_         );
    babyTree_->Branch("l1ps_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar"  , &l1ps_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_  );
    babyTree_->Branch("l1ps_ele27_WP80_vstar"                                             , &l1ps_ele27_WP80_vstar_                                             );
#endif

    // Muons
    babyTree_->Branch("mu3_vstar"                  , &mu3_vstar_                  ); 
    babyTree_->Branch("mu5_vstar"                  , &mu5_vstar_                  ); 
    babyTree_->Branch("mu8_vstar"                  , &mu8_vstar_                  ); 
    babyTree_->Branch("mu12_vstar"                 , &mu12_vstar_                 ); 
    babyTree_->Branch("mu15_vstar"                 , &mu15_vstar_                 ); 
    babyTree_->Branch("mu17_vstar"                 , &mu17_vstar_                 ); 
    babyTree_->Branch("mu20_vstar"                 , &mu20_vstar_                 ); 
    babyTree_->Branch("mu24_vstar"                 , &mu24_vstar_                 ); 
    babyTree_->Branch("mu30_vstar"                 , &mu30_vstar_                 ); 
    babyTree_->Branch("mu15_eta2p1_vstar"          , &mu15_eta2p1_vstar_          ); 
    babyTree_->Branch("mu24_eta2p1_vstar"          , &mu24_eta2p1_vstar_          ); 
    babyTree_->Branch("mu30_eta2p1_vstar"          , &mu30_eta2p1_vstar_          ); 
    babyTree_->Branch("mu8_Jet40_vstar"            , &mu8_Jet40_vstar_            ); 
    babyTree_->Branch("isoMu20_eta2p1_vstar"       , &isoMu20_eta2p1_vstar_       ); 
    babyTree_->Branch("isoMu24_eta2p1_vstar"       , &isoMu24_eta2p1_vstar_       ); 
    babyTree_->Branch("isoMu30_eta2p1_vstar"       , &isoMu30_eta2p1_vstar_       ); 
    babyTree_->Branch("relIso1p0Mu17_vstar"        , &relIso1p0Mu17_vstar_        ); 
    babyTree_->Branch("relIso1p0Mu20_vstar"        , &relIso1p0Mu20_vstar_        ); 
    babyTree_->Branch("relIso1p0Mu5_vstar"         , &relIso1p0Mu5_vstar_         ); 

    babyTree_->Branch("mu3_version"                , &mu3_version_                ); 
    babyTree_->Branch("mu5_version"                , &mu5_version_                ); 
    babyTree_->Branch("mu8_version"                , &mu8_version_                ); 
    babyTree_->Branch("mu12_version"               , &mu12_version_               ); 
    babyTree_->Branch("mu15_version"               , &mu15_version_               ); 
    babyTree_->Branch("mu17_version"               , &mu17_version_               ); 
    babyTree_->Branch("mu20_version"               , &mu20_version_               ); 
    babyTree_->Branch("mu24_version"               , &mu24_version_               ); 
    babyTree_->Branch("mu30_version"               , &mu30_version_               ); 
    babyTree_->Branch("mu15_eta2p1_version"        , &mu15_eta2p1_version_        ); 
    babyTree_->Branch("mu24_eta2p1_version"        , &mu24_eta2p1_version_        ); 
    babyTree_->Branch("mu30_eta2p1_version"        , &mu30_eta2p1_version_        ); 
    babyTree_->Branch("mu8_Jet40_version"          , &mu8_Jet40_version_          ); 
    babyTree_->Branch("isoMu20_eta2p1_version"     , &isoMu20_eta2p1_version_     ); 
    babyTree_->Branch("isoMu24_eta2p1_version"     , &isoMu24_eta2p1_version_     ); 
    babyTree_->Branch("isoMu30_eta2p1_version"     , &isoMu30_eta2p1_version_     ); 
    babyTree_->Branch("relIso1p0Mu17_version"      , &relIso1p0Mu17_version_      ); 
    babyTree_->Branch("relIso1p0Mu20_version"      , &relIso1p0Mu20_version_      ); 
    babyTree_->Branch("relIso1p0Mu5_version"       , &relIso1p0Mu5_version_       ); 

    babyTree_->Branch("dr_mu3_vstar"               , &dr_mu3_vstar_               ); 
    babyTree_->Branch("dr_mu5_vstar"               , &dr_mu5_vstar_               ); 
    babyTree_->Branch("dr_mu8_vstar"               , &dr_mu8_vstar_               ); 
    babyTree_->Branch("dr_mu12_vstar"              , &dr_mu12_vstar_              ); 
    babyTree_->Branch("dr_mu15_vstar"              , &dr_mu15_vstar_              ); 
    babyTree_->Branch("dr_mu17_vstar"              , &dr_mu17_vstar_              ); 
    babyTree_->Branch("dr_mu20_vstar"              , &dr_mu20_vstar_              ); 
    babyTree_->Branch("dr_mu24_vstar"              , &dr_mu24_vstar_              ); 
    babyTree_->Branch("dr_mu30_vstar"              , &dr_mu30_vstar_              ); 
    babyTree_->Branch("dr_mu15_eta2p1_vstar"       , &dr_mu15_eta2p1_vstar_       ); 
    babyTree_->Branch("dr_mu24_eta2p1_vstar"       , &dr_mu24_eta2p1_vstar_       ); 
    babyTree_->Branch("dr_mu30_eta2p1_vstar"       , &dr_mu30_eta2p1_vstar_       ); 
    babyTree_->Branch("dr_mu8_Jet40_vstar"         , &dr_mu8_Jet40_vstar_         ); 
    babyTree_->Branch("dr_isoMu20_eta2p1_vstar"    , &dr_isoMu20_eta2p1_vstar_    ); 
    babyTree_->Branch("dr_isoMu24_eta2p1_vstar"    , &dr_isoMu24_eta2p1_vstar_    ); 
    babyTree_->Branch("dr_isoMu30_eta2p1_vstar"    , &dr_isoMu30_eta2p1_vstar_    ); 
    babyTree_->Branch("dr_relIso1p0Mu17_vstar"     , &dr_relIso1p0Mu17_vstar_     ); 
    babyTree_->Branch("dr_relIso1p0Mu20_vstar"     , &dr_relIso1p0Mu20_vstar_     ); 
    babyTree_->Branch("dr_relIso1p0Mu5_vstar"      , &dr_relIso1p0Mu5_vstar_      ); 

    babyTree_->Branch("hltps_mu3_vstar"            , &hltps_mu3_vstar_            ); 
    babyTree_->Branch("hltps_mu5_vstar"            , &hltps_mu5_vstar_            ); 
    babyTree_->Branch("hltps_mu8_vstar"            , &hltps_mu8_vstar_            ); 
    babyTree_->Branch("hltps_mu12_vstar"           , &hltps_mu12_vstar_           ); 
    babyTree_->Branch("hltps_mu15_vstar"           , &hltps_mu15_vstar_           ); 
    babyTree_->Branch("hltps_mu17_vstar"           , &hltps_mu17_vstar_           ); 
    babyTree_->Branch("hltps_mu20_vstar"           , &hltps_mu20_vstar_           ); 
    babyTree_->Branch("hltps_mu24_vstar"           , &hltps_mu24_vstar_           ); 
    babyTree_->Branch("hltps_mu30_vstar"           , &hltps_mu30_vstar_           ); 
    babyTree_->Branch("hltps_mu15_eta2p1_vstar"    , &hltps_mu15_eta2p1_vstar_    ); 
    babyTree_->Branch("hltps_mu24_eta2p1_vstar"    , &hltps_mu24_eta2p1_vstar_    ); 
    babyTree_->Branch("hltps_mu30_eta2p1_vstar"    , &hltps_mu30_eta2p1_vstar_    ); 
    babyTree_->Branch("hltps_mu8_Jet40_vstar"      , &hltps_mu8_Jet40_vstar_      ); 
    babyTree_->Branch("hltps_isoMu20_eta2p1_vstar" , &hltps_isoMu20_eta2p1_vstar_ ); 
    babyTree_->Branch("hltps_isoMu24_eta2p1_vstar" , &hltps_isoMu24_eta2p1_vstar_ ); 
    babyTree_->Branch("hltps_isoMu30_eta2p1_vstar" , &hltps_isoMu30_eta2p1_vstar_ ); 
    babyTree_->Branch("hltps_relIso1p0Mu17_vstar"  , &hltps_relIso1p0Mu17_vstar_  ); 
    babyTree_->Branch("hltps_relIso1p0Mu20_vstar"  , &hltps_relIso1p0Mu20_vstar_  ); 
    babyTree_->Branch("hltps_relIso1p0Mu5_vstar"   , &hltps_relIso1p0Mu5_vstar_   ); 

#ifndef __CMS2_SLIM__
    babyTree_->Branch("l1ps_mu5_vstar"             , &l1ps_mu5_vstar_             ); 
    babyTree_->Branch("l1ps_mu8_vstar"             , &l1ps_mu8_vstar_             ); 
    babyTree_->Branch("l1ps_mu12_vstar"            , &l1ps_mu12_vstar_            ); 
    babyTree_->Branch("l1ps_mu17_vstar"            , &l1ps_mu17_vstar_            ); 
    babyTree_->Branch("l1ps_mu15_eta2p1_vstar"     , &l1ps_mu15_eta2p1_vstar_     ); 
    babyTree_->Branch("l1ps_mu24_eta2p1_vstar"     , &l1ps_mu24_eta2p1_vstar_     ); 
    babyTree_->Branch("l1ps_mu30_eta2p1_vstar"     , &l1ps_mu30_eta2p1_vstar_     ); 
    babyTree_->Branch("l1ps_isoMu20_eta2p1_vstar"  , &l1ps_isoMu20_eta2p1_vstar_  ); 
    babyTree_->Branch("l1ps_isoMu24_eta2p1_vstar"  , &l1ps_isoMu24_eta2p1_vstar_  ); 
    babyTree_->Branch("l1ps_isoMu30_eta2p1_vstar"  , &l1ps_isoMu30_eta2p1_vstar_  ); 
    babyTree_->Branch("l1ps_relIso1p0Mu17_vstar"   , &l1ps_relIso1p0Mu17_vstar_   ); 
    babyTree_->Branch("l1ps_relIso1p0Mu5_vstar"    , &l1ps_relIso1p0Mu5_vstar_    ); 
#endif

    ///////////////////////  
    // End 2011 Triggers //
    ///////////////////////
        
    //////////////
    // Jets     //
    //////////////

    // Information to do offline jet trigger selection
#ifndef __CMS2_SLIM__
    babyTree_->Branch("ptj1"         , &ptj1_         );
    babyTree_->Branch("nj1"          , &nj1_          );
    babyTree_->Branch("ptj1_b2b"     , &ptj1_b2b_     );
    babyTree_->Branch("dphij1_b2b"   , &dphij1_b2b_   );
#endif
    babyTree_->Branch("ptpfj1"       , &ptpfj1_       );
    babyTree_->Branch("npfj1"        , &npfj1_        );
    babyTree_->Branch("ptpfj1_b2b"   , &ptpfj1_b2b_   );
    babyTree_->Branch("dphipfj1_b2b" , &dphipfj1_b2b_ );
      
    // PF L2L3 Corrected jets
    babyTree_->Branch("ptpfcj1"      , &ptpfcj1_      );
    babyTree_->Branch("npfcj1"       , &npfcj1_       );
    babyTree_->Branch("ptpfcj1_b2b"  , &ptpfcj1_b2b_  );
    babyTree_->Branch("dphipfcj1_b2b", &dphipfcj1_b2b_);
    babyTree_->Branch("btagpfc"      , &btagpfc_      );
      
    // PF L1FastL2L3 Corrected jets         
    babyTree_->Branch("emfpfcL1Fj1"     , &emfpfcL1Fj1_      );
    babyTree_->Branch("ptpfcL1Fj1"      , &ptpfcL1Fj1_       );       
    babyTree_->Branch("dphipfcL1Fj1"    , &dphipfcL1Fj1_     );       
    babyTree_->Branch("npfcL1Fj1"       , &npfcL1Fj1_        );
    babyTree_->Branch("npfc30L1Fj1"     , &npfc30L1Fj1_      );
    babyTree_->Branch("npfc40L1Fj1"     , &npfc40L1Fj1_      );
    babyTree_->Branch("nbpfc40L1Fj1"    , &nbpfc40L1Fj1_     );
    babyTree_->Branch("ptpfcL1Fj1_b2b"  , &ptpfcL1Fj1_b2b_   );       
    babyTree_->Branch("dphipfcL1Fj1_b2b", &dphipfcL1Fj1_b2b_ );     
    babyTree_->Branch("btagpfcL1F"      , &btagpfcL1F_       );
    babyTree_->Branch("npfc50L1Fj1_eth" , &npfc50L1Fj1_eth_  );
    babyTree_->Branch("npfc65L1Fj1_eth" , &npfc65L1Fj1_eth_  );

    // PF L1FastL2L3Residual Corrected jets         
    babyTree_->Branch("emfpfcL1Fj1res"     , &emfpfcL1Fj1res_      );
    babyTree_->Branch("ptpfcL1Fj1res"      , &ptpfcL1Fj1res_       );       
    babyTree_->Branch("dphipfcL1Fj1res"    , &dphipfcL1Fj1res_     );       
    babyTree_->Branch("npfcL1Fj1res"       , &npfcL1Fj1res_        );
    babyTree_->Branch("npfc30L1Fj1res"     , &npfc30L1Fj1res_      );
    babyTree_->Branch("npfc40L1Fj1res"     , &npfc40L1Fj1res_      );
    babyTree_->Branch("nbpfc40L1Fj1res"    , &nbpfc40L1Fj1res_     );
    babyTree_->Branch("ptpfcL1Fj1res_b2b"  , &ptpfcL1Fj1res_b2b_   );       
    babyTree_->Branch("dphipfcL1Fj1res_b2b", &dphipfcL1Fj1res_b2b_ );     
    babyTree_->Branch("btagpfcL1Fres"      , &btagpfcL1Fres_       );
    babyTree_->Branch("npfc50L1Fj1res_eth" , &npfc50L1Fj1res_eth_  );
    babyTree_->Branch("npfc65L1Fj1res_eth" , &npfc65L1Fj1res_eth_  );

    // B-tagged PF L1FastL2L3 Corrected jets         
    babyTree_->Branch("ptbtagpfcL1Fj1"   , &ptbtagpfcL1Fj1_   );
    babyTree_->Branch("dphibtagpfcL1Fj1" , &dphibtagpfcL1Fj1_ );
    
    // B-tagged PF L1FastL2L3Residual Corrected jets         
    babyTree_->Branch("ptbtagpfcL1Fj1res"   , &ptbtagpfcL1Fj1res_   );
    babyTree_->Branch("dphibtagpfcL1Fj1res" , &dphibtagpfcL1Fj1res_ );
    
#ifndef __CMS2_SLIM__
    // B Tagging
    babyTree_->Branch("nbjet"        , &nbjet_        );
    babyTree_->Branch("dRNear"       , &dRbNear_      );
    babyTree_->Branch("dRFar"        , &dRbFar_       );
#endif
      
    babyTree_->Branch("nbpfcjet"     , &nbpfcjet_     );
    babyTree_->Branch("dRpfcNear"    , &dRbpfcNear_   );
    babyTree_->Branch("dRpfcFar"     , &dRbpfcFar_    );

    babyTree_->Branch("rho", &rho_);
    
    //////////////
    // End Jets //
    //////////////
}

// Fill the baby
void myBabyMaker::FillBabyNtuple()
{ 
    babyTree_->Fill(); 
}

// Close the baby
void myBabyMaker::CloseBabyNtuple()
{
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();
}

// constructor
myBabyMaker::myBabyMaker () 
    : nEvents_(-1)
    , verbose_(false)
    , ele8_regexp                                                 ("HLT_Ele8_v(\\d+)"                                                 , "o")
    , ele8_CaloIdL_TrkIdVL_regexp                                 ("HLT_Ele8_CaloIdL_TrkIdVL_v(\\d+)"                                 , "o")
    , ele8_CaloIdL_CaloIsoVL_regexp                               ("HLT_Ele8_CaloIdL_CaloIsoVL_v(\\d+)"                               , "o")
    , ele8_CaloIdL_CaloIsoVL_Jet40_regexp                         ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v(\\d+)"                         , "o")
    , ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_regexp              ("HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v(\\d+)"              , "o")
    , ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_regexp              ("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v(\\d+)"              , "o")
    , ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_regexp        ("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v(\\d+)"        , "o")
    , ele8_CaloIdT_TrkIdVL_regexp                                 ("HLT_Ele8_CaloIdT_TrkIdVL_v(\\d+)"                                 , "o")
    , ele8_CaloIdT_TrkIdVL_Jet30_regexp                           ("HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v(\\d+)"                           , "o")
    , ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_regexp        ("HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v(\\d+)"              , "o")
    , ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_regexp        ("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v(\\d+)"              , "o")
    , ele17_CaloIdL_CaloIsoVL_regexp                              ("HLT_Ele17_CaloIdL_CaloIsoVL_v(\\d+)"                              , "o")
    , ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_regexp             ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v(\\d+)"             , "o")
    , ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_regexp       ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v(\\d+)"       , "o")
    , ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_rexexp("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_v(\\d+)", "o")
    , ele27_WP80_rexexp                                           ("HLT_Ele27_WP80_v(\\d+)"                                           , "o")
    , photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_regexp        ("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v(\\d+)"        , "o")
    , mu3_regexp            ("HLT_Mu3_v(\\d+)"           , "o")
    , mu5_regexp            ("HLT_Mu5_v(\\d+)"           , "o")          
    , mu8_regexp            ("HLT_Mu8_v(\\d+)"           , "o")      
    , mu12_regexp           ("HLT_Mu12_v(\\d+)"          , "o")     
    , mu15_regexp           ("HLT_Mu15_v(\\d+)"          , "o")     
    , mu17_regexp           ("HLT_Mu17_v(\\d+)"          , "o")     
    , mu20_regexp           ("HLT_Mu20_v(\\d+)"          , "o")     
    , mu24_regexp           ("HLT_Mu24_v(\\d+)"          , "o")     
    , mu30_regexp           ("HLT_Mu30_v(\\d+)"          , "o")     
    , mu15_eta2p1_regexp    ("HLT_Mu15_eta2p1_v(\\d+)"   , "o")     
    , mu24_eta2p1_regexp    ("HLT_Mu24_eta2p1_v(\\d+)"   , "o")     
    , mu30_eta2p1_regexp    ("HLT_Mu30_eta2p1_v(\\d+)"   , "o")     
    , mu8_Jet40_regexp      ("HLT_Mu8_Jet40_v(\\d+)"     , "o")
    , isoMu20_eta2p1_regexp ("HLT_IsoMu20_eta2p1_v(\\d+)", "o")
    , isoMu24_eta2p1_regexp ("HLT_IsoMu24_eta2p1_v(\\d+)", "o")
    , isoMu30_eta2p1_regexp ("HLT_IsoMu30_eta2p1_v(\\d+)", "o")
    , relIso1p0Mu17_regexp  ("HLT_RelIso1p0Mu17_v(\\d+)" , "o")
    , relIso1p0Mu20_regexp  ("HLT_RelIso1p0Mu20_v(\\d+)" , "o")
    , relIso1p0Mu5_regexp   ("HLT_RelIso1p0Mu5_v(\\d+)"  , "o")
{
}

//-----------------------------------
// Looper code starts here
// eormu=-1 do both e and mu
//      =11 do electrons
//      =13 do muons
//-----------------------------------
void myBabyMaker::ScanChain(TChain* chain, const char *babyFilename, int eormu, bool applyFOfilter, const std::string& jetcorrPath)
{
    try
    {
        already_seen.clear();

        // Make a baby ntuple
        MakeBabyNtuple(babyFilename);

        // Jet Corrections
        std::vector<std::string> jetcorr_pf_L2L3_filenames;
        //if (isData) {        
        string data_pf_l2 = jetcorrPath;
        data_pf_l2.append("/GR_R_52_V7_L2Relative_AK5PF.txt");
        string data_pf_l3 = jetcorrPath;
        data_pf_l3.append("/GR_R_52_V7_L3Absolute_AK5PF.txt");
        jetcorr_pf_L2L3_filenames.push_back(data_pf_l2.c_str());
        jetcorr_pf_L2L3_filenames.push_back(data_pf_l3.c_str());
        //}
        //else {
        //    string mc_pf_l2 = jetcorrPath;
        //    mc_pf_l2.append("/START41_V0_AK5PF_L2Relative.txt");
        //    string mc_pf_l3 = jetcorrPath;
        //    mc_pf_l3.append("/START41_V0_AK5PF_L3Absolute.txt");
        //    jetcorr_pf_L2L3_filenames.push_back(mc_pf_l2.c_str());
        //    jetcorr_pf_L2L3_filenames.push_back(mc_pf_l3.c_str());
        //}    

        std::cout << "making jet corrector with the following files: " << std::endl;
        for (unsigned int idx = 0; idx < jetcorr_pf_L2L3_filenames.size(); idx++)
            std::cout << jetcorr_pf_L2L3_filenames.at(idx) << std::endl;
        FactorizedJetCorrector *jet_pf_L2L3corrector = makeJetCorrector(jetcorr_pf_L2L3_filenames);

        //// set up on-the-fly L1FastJetL2L3 JEC
        //std::vector<std::string> jetcorr_pf_L1FastJetL2L3_filenames;
        //jetcorr_pf_L1FastJetL2L3_filenames.push_back(Form("%s/GR_R_52_V7_L1FastJet_AK5PF.txt"   , jetcorrPath.empty() ? "." : jetcorrPath.c_str()));
        //jetcorr_pf_L1FastJetL2L3_filenames.push_back(Form("%s/GR_R_52_V7_L2Relative_AK5PF.txt"  , jetcorrPath.empty() ? "." : jetcorrPath.c_str()));
        //jetcorr_pf_L1FastJetL2L3_filenames.push_back(Form("%s/GR_R_52_V7_L3Absolute_AK5PF.txt"  , jetcorrPath.empty() ? "." : jetcorrPath.c_str()));
        //FactorizedJetCorrector* jet_pf_L1FastJetL2L3_corrector = makeJetCorrector(jetcorr_pf_L1FastJetL2L3_filenames); 

        //// set up on-the-fly L1FastJetL2L3 residual JEC
        //std::vector<std::string> jetcorr_pf_L1FastJetL2L3Residual_filenames;
        //jetcorr_pf_L1FastJetL2L3Residual_filenames.push_back(Form("%s/GR_R_52_V7_L1FastJet_AK5PF.txt"   , jetcorrPath.empty() ? "." : jetcorrPath.c_str()));
        //jetcorr_pf_L1FastJetL2L3Residual_filenames.push_back(Form("%s/GR_R_52_V7_L2Relative_AK5PF.txt"  , jetcorrPath.empty() ? "." : jetcorrPath.c_str()));
        //jetcorr_pf_L1FastJetL2L3Residual_filenames.push_back(Form("%s/GR_R_52_V7_L3Absolute_AK5PF.txt"  , jetcorrPath.empty() ? "." : jetcorrPath.c_str()));
        //jetcorr_pf_L1FastJetL2L3Residual_filenames.push_back(Form("%s/GR_R_52_V7_L2L3Residual_AK5PF.txt", jetcorrPath.empty() ? "." : jetcorrPath.c_str()));
        //FactorizedJetCorrector* jet_pf_L1FastJetL2L3Residual_corrector = makeJetCorrector(jetcorr_pf_L1FastJetL2L3Residual_filenames); 

        // The deltaR requirement between objects and jets to remove the jet trigger dependence
        float deltaRCut   = 1.0;
        float deltaPhiCut = 2.5;

        //--------------------------
        // File and Event Loop
        //---------------------------

        // benchmark
        TBenchmark bmark;
        bmark.Start("benchmark");

        int i_permilleOld = 0;
        unsigned int nEventsTotal = 0;
        unsigned int nEventsChain = 0;
        int nEvents = nEvents_; 
        if (nEvents==-1){
            nEventsChain = chain->GetEntries();
        } else {
            nEventsChain = nEvents;
        }
        TObjArray *listOfFiles = chain->GetListOfFiles();
        TIter fileIter(listOfFiles);
        bool finish_looping = false;

        std::cout << "looping on " << nEventsChain << " out of " << chain->GetEntries() << " events..." << std::endl;
        std::cout << "nEventTotal = " << nEventsTotal << endl;
        std::cout << "nEventChain = " << nEventsChain << endl;

        while(TChainElement *currentFile = (TChainElement*)fileIter.Next())
        {
            if (finish_looping) {
                break;
            }

            TString filename = currentFile->GetTitle();
            if (verbose_)
            {
                cout << filename << endl;
            }

            TFile* f = TFile::Open(filename.Data());
            TTree* tree = (TTree*)f->Get("Events");
            cms2.Init(tree);

            unsigned int nEntries = tree->GetEntries();
            unsigned int nGoodEvents(0);
            unsigned int nLoop = nEntries;
            unsigned int z;

            // Event Loop
            for( z = 0; z < nLoop; z++)
            { 
                cms2.GetEntry(z);

                if (nEventsTotal >= nEventsChain) {
                    finish_looping = true;
                    break;
                }

                bool isData = evt_isRealData();

                if(isData){
                    // Good  Runs
                    if (goodrun_is_json) {
                        if(!goodrun_json(evt_run(), evt_lumiBlock())) continue;   
                    }
                    else {
                        if(!goodrun(evt_run(), evt_lumiBlock())) continue;   
                    }

                    // check for duplicated
                    DorkyEventIdentifier id = {evt_run(), evt_event(), evt_lumiBlock()};
                    if (is_duplicate(id) ) { 
                        cout << "\t! ERROR: found duplicate." << endl;
                        continue;
                    }
                }

                // looper progress
                ++nEventsTotal;
                ++nGoodEvents;
                int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
                if (i_permille != i_permilleOld) {
                    printf("  \015\033[32m ---> \033[1m\033[31m%4.1f%%" "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
                    fflush(stdout);
                    i_permilleOld = i_permille;
                }

                // Event cleaning (careful, it requires technical bits)
                //if (!cleaning_BPTX(isData))   continue;
                //if (!cleaning_beamHalo())   continue;
                //if (!cleaning_goodVertexAugust2010()) continue;
                //if (!cleaning_goodTracks()) continue;
                if (!cleaning_standardApril2011()) continue;

                // Loop over jets and see what is btagged
                // Medium operating point from https://twiki.cern.ch/twiki/bin/view/CMS/BTagPerformanceOP

                // #ifndef __CMS2_SLIM__
                //          int this_nbjet = 0;
                //             vector<unsigned int> bindex;
                //             for (unsigned int iJet = 0; iJet < jets_p4().size(); iJet++) {
                //                 if (jets_p4().at(iJet).pt() < 15.) continue;
                //                 if (cms2.jets_combinedSecondaryVertexBJetTag().at(iJet) < 0.679) continue;
                //                 this_nbjet++;
                //                 bindex.push_back(iJet);
                //             }
                // #endif

                // PF Jets
                int this_nbpfjet = 0;
                //FactorizedJetCorrector* jet_pf_corrector = evt_isRealData() ? jet_pf_L1FastJetL2L3Residual_corrector : jet_pf_L1FastJetL2L3_corrector; 
                vector<unsigned int> bpfindex;
                for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                    if ( !passesPFJetID(iJet)) continue;
                    LorentzVector jp4 = pfjets_p4().at(iJet);
                    //jet_pf_corrector->setRho(cms2.evt_ww_rho_vor());
                    //jet_pf_corrector->setJetA(cms2.pfjets_area().at(iJet));
                    //jet_pf_corrector->setJetPt(cms2.pfjets_p4().at(iJet).pt());
                    //jet_pf_corrector->setJetEta(cms2.pfjets_p4().at(iJet).eta()); 
                    //float jet_cor = jetCorrection(jp4, jet_pf_corrector);
                    float jet_cor = evt_isRealData() ? cms2.pfjets_corL1FastL2L3residual().at(iJet) : cms2.pfjets_corL1FastL2L3().at(iJet);
                    LorentzVector jp4cor = jp4 * jet_cor;
                    if (jp4cor.pt() < 15) continue;
                    if (cms2.pfjets_combinedSecondaryVertexBJetTag().at(iJet) < 0.679) continue;
                    this_nbpfjet++;
                    bpfindex.push_back(iJet);
                }

                // Electrons
                if (eormu == -1 || eormu==11) {
                    for (unsigned int iLep = 0 ; iLep < els_p4().size(); iLep++) {

                        // Apply a pt cut (Changed it from 5 GeV to 10 GeV...Claudio 10 July 2010)
                        if ( els_p4().at(iLep).pt() < 10.) continue;

                        // Initialize baby ntuple
                        InitBabyNtuple();

                        //////////////////////////////////////////////////////
                        // Fake Rate Numerator & Denominator Selections     //
                        //////////////////////////////////////////////////////

                        // store number of electron FOs in event (use SS FO definition)
                        nFOels_ = 0;
                        ngsfs_ = 0;
                        for (unsigned int iel = 0; iel < cms2.els_p4().size(); iel++)
                        {
                            if (iel == iLep)
                                continue;

                            if (cms2.els_p4().at(iel).pt() < 10.)
                                continue;

                            if (pass_electronSelection(iel, electronSelectionFOV7_v3, false, false))
                                ++ngsfs_;

                            if (samesign::isDenominatorLepton(11, iel)) {
                                ++nFOels_;
                                if (cms2.els_p4().at(iel).pt() > foel_p4_.pt() && iel != iLep) {
                                    foel_p4_ = cms2.els_p4().at(iel);
                                    foel_id_ = 11*cms2.els_charge().at(iel);
                                }
                                continue;
                            }
                            if (samesign::isDenominatorLepton(11, iel)) {
                                ++nFOels_;
                                if (cms2.els_p4().at(iel).pt() > foel_p4_.pt() && iel != iLep) {
                                    foel_p4_ = cms2.els_p4().at(iel);
                                    foel_id_ = 11*cms2.els_charge().at(iel);
                                }
                                continue;
                            }
                        }

                        // store number of muon FOs in event (use SS FO definition)
                        nFOmus_ = 0;
                        nmus_ = 0;
                        for (unsigned int imu = 0; imu < cms2.mus_p4().size(); imu++)
                        {
                            if (cms2.mus_p4().at(imu).pt() < 10.)
                                continue;

                            if (muonIdNotIsolated(imu, muonSelectionFO_ssV5))
                                ++nmus_;

                            if (samesign::isDenominatorLepton(13, imu)) {
                                ++nFOmus_;
                                if (cms2.mus_p4().at(imu).pt() > fomu_p4_.pt()) {
                                    fomu_p4_ = cms2.mus_p4().at(imu);
                                    fomu_id_ = 13*cms2.mus_charge().at(imu);
                                }
                                continue;
                            }
                            if (samesign::isDenominatorLepton(13, imu)) {
                                ++nFOmus_;
                                if (cms2.mus_p4().at(imu).pt() > fomu_p4_.pt()) {
                                    fomu_p4_ = cms2.mus_p4().at(imu);
                                    fomu_id_ = 13*cms2.mus_charge().at(imu);
                                }
                                continue;
                            }
                        }

                        // store number of "veto" electrons in event
                        nvetoels_ = 0;
                        for (unsigned int iel = 0; iel < cms2.els_p4().size(); iel++)
                        {
                            if (iel == iLep) // skip the current electron
                                continue;

                            if (cms2.els_p4().at(iel).pt() < 5.0f)
                                continue;

                            if (fabs(cms2.els_p4().at(iel).eta()) > 2.4f)
                                continue;

                            const float iso = electronIsoValuePF2012_FastJetEffArea_v3(iel, /*conesize=*/0.3, /*vtx=*/-999, /*52X iso=*/false);
                            if (iso > 1.0)
                                continue;

                            // if we get here, then count it
                            nvetoels_++;
                        }

                        // store number of "veto" muons in event
                        nvetomus_ = 0;
                        for (unsigned int imu = 0; imu < cms2.mus_p4().size(); imu++)
                        {
                            if (cms2.mus_p4().at(imu).pt() < 5.0f)
                                continue;

                            if (fabs(cms2.mus_p4().at(imu).eta()) > 2.4f)
                                continue;

                            const bool is_global      = ((cms2.mus_type().at(imu) & (1<<1)) != 0);
                            const bool is_tracker     = ((cms2.mus_type().at(imu) & (1<<2)) != 0);
                            const bool is_pfmu        = ((cms2.mus_type().at(imu) & (1<<5)) != 0);
                            const bool passes_mu_type = ((is_global or is_tracker) and is_pfmu);
                            if (not passes_mu_type)
                                continue;

                            const float iso = muonIsoValuePF2012_deltaBeta(imu); 
                            if (iso > 1.0)
                                continue;

                            // if we get here, then count it
                            nvetomus_++;
                        }

                        //////////
                        // 2012 //
                        //////////

                        // SS
                        num_el_ssV7_       = pass_electronSelection(iLep, electronSelection_ssV7       );
                        num_el_ssV7_noIso_ = pass_electronSelection(iLep, electronSelection_ssV7_noIso );
                        v1_el_ssV7_        = pass_electronSelection(iLep, electronSelectionFOV7_v1     );
                        v2_el_ssV7_        = pass_electronSelection(iLep, electronSelectionFOV7_v2     );
                        v3_el_ssV7_        = pass_electronSelection(iLep, electronSelectionFOV7_v3     );

                        // TTZ

                        num_el_TTZcuttightv1_       = ttv::isNumeratorLepton(11, iLep, ttv::LeptonType::TIGHT);
                        num_el_TTZcuttightv1_noIso_ = ttv::isGoodLepton(11, iLep, ttv::LeptonType::TIGHT);
                        fo_el_TTZcuttightv1_        = (ttv::isDenominatorLepton(11, iLep, ttv::LeptonType::TIGHT) && electronIsoValuePF2012_FastJetEffArea_v2( iLep ) < 0.6);
                        fo_el_TTZcuttightv1_noIso_  = ttv::isDenominatorLepton(11, iLep, ttv::LeptonType::TIGHT);

                        num_el_TTZcutloosev1_       = ttv::isNumeratorLepton(11, iLep, ttv::LeptonType::LOOSE);
                        num_el_TTZcutloosev1_noIso_ = ttv::isGoodLepton(11, iLep, ttv::LeptonType::LOOSE);
                        fo_el_TTZcutloosev1_        = (ttv::isDenominatorLepton(11, iLep, ttv::LeptonType::LOOSE) && electronIsoValuePF2012_FastJetEffArea_v2( iLep ) < 0.6);
                        fo_el_TTZcutloosev1_noIso_  = ttv::isDenominatorLepton(11, iLep, ttv::LeptonType::LOOSE);

                        num_el_TTZMVAtightv1_       = false;
                        num_el_TTZMVAtightv1_noIso_ = false;
                        fo_el_TTZMVAtightv1_        = false;
                        fo_el_TTZMVAtightv1_noIso_  = false;

                        num_el_TTZMVAloosev1_       = false;
                        num_el_TTZMVAloosev1_noIso_ = false;
                        fo_el_TTZMVAloosev1_        = false;
                        fo_el_TTZMVAloosev1_noIso_  = false;

                        //////////
                        // 2011 //
                        //////////

                        // SS
                        num_el_ssV6_       = pass_electronSelection( iLep, electronSelection_ssV6            );
                        num_el_ssV6_noIso_ = pass_electronSelection( iLep, electronSelection_ssV6_noIso      );
                        v1_el_ssV6_        = pass_electronSelection( iLep, electronSelectionFOV6_ssVBTF80_v1 );
                        v2_el_ssV6_        = pass_electronSelection( iLep, electronSelectionFOV6_ssVBTF80_v2 );
                        v3_el_ssV6_        = pass_electronSelection( iLep, electronSelectionFOV6_ssVBTF80_v3 );

                        // WW
                        num_el_smurfV6_ = pass_electronSelection( iLep, electronSelection_smurfV6          );
                        v1_el_smurfV1_  = pass_electronSelection( iLep, electronSelectionFO_el_smurf_v1    );
                        v2_el_smurfV1_  = pass_electronSelection( iLep, electronSelectionFO_el_smurf_v2    );
                        v3_el_smurfV1_  = pass_electronSelection( iLep, electronSelectionFO_el_smurf_v3    );
                        v4_el_smurfV1_  = pass_electronSelection( iLep, electronSelectionFO_el_smurf_v4    );

                        //OS
                        num_el_OSV2_   = pass_electronSelection( iLep, electronSelection_el_OSV2          );
                        fo_el_OSV2_    = pass_electronSelection( iLep, electronSelection_el_OSV2_FO       );
                        num_el_OSV3_   = pass_electronSelection( iLep, electronSelection_el_OSV3          );
                        fo_el_OSV3_    = pass_electronSelection( iLep, electronSelection_el_OSV3_FO       );

                        ////////////////////////////////////////////////////////////
                        // Skip this muon if it fails the loosest denominator.    //
                        // Ignore this OR if applyFOfilter is set to false.       //
                        ////////////////////////////////////////////////////////////
                        if (applyFOfilter) {  
                            if (
                                !v1_el_ssV7_          && !v2_el_ssV7_          && !v3_el_ssV7_    &&                    // SS 2012
                                !fo_el_TTZMVAtightv1_ && !fo_el_TTZMVAloosev1_ &&                                       // TTZ MVA 2012
                                !fo_el_TTZcuttightv1_ && !fo_el_TTZcutloosev1_ &&                                       // TTZ cut 2012
                                !v1_el_ssV6_          && !v2_el_ssV6_          && !v3_el_ssV6_    &&                    // SS 2011
                                !fo_el_OSV2_          && !fo_el_OSV3_          &&                                       // OS 2011
                                !v1_el_smurfV1_       && !v1_el_smurfV1_       && !v3_el_smurfV1_ && !v4_el_smurfV1_    // WW 2011
                                ) 
                                continue;
                        }

                        //////////////////////////////////////////////////////
                        // End Fake Rate Numerator & Denominator Selections //
                        //////////////////////////////////////////////////////


                        ////////////////////////////////////////////////////////////////////////
                        // NEED TO THINK ABOUT THIS... Z'S ARE VETOED BASED ON TOP SELECTIONS //
                        ////////////////////////////////////////////////////////////////////////

                        // If it is above 20 GeV see if we can make a 
                        // Z with another pt>20 FO.  Will use the v1 FO since 
                        // these are the loosest
                        bool isaZ = false;
                        if (els_p4().at(iLep).pt() > 20.) {
                            for (unsigned int jEl = 0 ; jEl < els_p4().size(); jEl++) {
                                if (iLep == jEl)                             continue;
                                if (els_p4().at(jEl).pt() < 20.)            continue;
                                if ( ! pass_electronSelection( jEl, electronSelection_el_OSV3_FO ) ) continue;
                                if ( ! fo_el_OSV3_ ) continue;
                                LorentzVector w = els_p4().at(iLep) + els_p4().at(jEl);
                                if (abs(w.mass()-91.) > 20.) continue;
                                isaZ = true;
                            }
                        }
                        if (isaZ) continue;

                        ////////////////////////////////////////////////////////////////////////
                        // STORE SOME Z MASS VARIABLES //
                        ////////////////////////////////////////////////////////////////////////
                        mz_fo_gsf_  = -999.;
                        mz_gsf_iso_ = -999.;
                        LorentzVector p4fo = cms2.els_p4().at(iLep);
                        for (unsigned int iel = 0; iel < cms2.els_p4().size(); iel++) {
                            if (iel == iLep) continue;

                            if (fabs(cms2.els_p4().at(iel).eta()) > 2.5)
                                continue;

                            if (cms2.els_p4().at(iel).pt() < 10.)
                                continue;

                            LorentzVector zp4 = p4fo + cms2.els_p4().at(iel);
                            float zcandmass = sqrt(fabs(zp4.mass2()));
                            if ( fabs(zcandmass - 91.) > fabs(mz_fo_gsf_ - 91.) )
                                continue;

                            mz_fo_gsf_  = zcandmass;
                            mz_gsf_iso_ = electronIsolation_rel_v1(iel, true);
                        }

                        mz_fo_ctf_  = -999.;
                        mz_ctf_iso_ = -999.;
                        for (int ictf = 0; ictf < static_cast<int>(cms2.trks_trk_p4().size()); ictf++) {
                            if (ictf == cms2.els_trkidx().at(iLep)) continue;

                            if (fabs(cms2.trks_trk_p4().at(ictf).eta()) > 2.5)
                                continue;

                            if (cms2.trks_trk_p4().at(ictf).pt() < 10.)
                                continue;

                            LorentzVector zp4 = p4fo + cms2.trks_trk_p4().at(ictf);
                            float zcandmass = sqrt(fabs(zp4.mass2()));
                            if ( fabs(zcandmass - 91.) > fabs(mz_fo_ctf_ - 91.) )
                                continue;

                            mz_fo_ctf_  = zcandmass;
                            mz_ctf_iso_ = ctfIsoValuePF(ictf, associateTrackToVertex(ictf));
                        }


                        /////////////////////////// 
                        // Event Information     //
                        ///////////////////////////

                        // Load the electron and event quantities
                        run_          = evt_run();
                        ls_           = evt_lumiBlock();
                        evt_          = evt_event();
                        weight_       = isData ? 1.0 : evt_scale1fb();
                        filename_     = f->GetName();
                        dataset_      = evt_dataset().front();
                        is_real_data_ = evt_isRealData();

                        if(!isData){
                            // Pileup - PUSummaryInfoMaker                        
                            for (unsigned int vidx = 0; vidx < cms2.puInfo_nPUvertices().size(); vidx++) {
                                if (cms2.puInfo_bunchCrossing().at(vidx) != 0)
                                    continue;
                                pu_nPUvertices_ = cms2.puInfo_nPUvertices().at(vidx);
                                pu_nPUtrueint_  = cms2.puInfo_trueNumInteractions().at(vidx);
                            }

                        }

                        // Pileup - VertexMaker
                        bool first_good_vertex_found         = false;
                        unsigned int first_good_vertex_index = 0;
                        for (unsigned int vidx = 0; vidx < cms2.vtxs_position().size(); vidx++)
                        {
                            if (!isGoodVertex(vidx))
                            {
                                continue;
                            }
                            if (!first_good_vertex_found)
                            {
                                first_good_vertex_found = true;
                                first_good_vertex_index = vidx;
                            }
                            ++evt_nvtxs_;
                        }

                        /////////////////////////// 
                        // End Event Information //
                        ///////////////////////////



                        //////////////////////////// 
                        // Lepton Information     //
                        ////////////////////////////

                        // Basic Quantities
                        lp4_       = cms2.els_p4().at(iLep);
                        pt_        = els_p4().at(iLep).pt();
                        eta_       = els_p4().at(iLep).eta();
                        sceta_     = els_etaSC().at(iLep);
                        phi_       = els_p4().at(iLep).phi();
                        scet_      = els_eSC().at(iLep) / cosh( els_etaSC().at(iLep));
                        hoe_       = els_hOverE().at(iLep);
                        id_        = 11*els_charge().at(iLep);

                        // ip (2d and 3d)
                        const int elgsftkid = cms2.els_gsftrkidx().at(iLep);
                        const int eltkid    = cms2.els_trkidx().at(iLep);
                        const int ivtx      = firstGoodVertex();
                        if (ivtx >= 0 && eltkid >= 0) 
                        {
                            d0_        = elgsftkid>=0 ? gsftrks_d0_pv(elgsftkid,ivtx).first  : trks_d0_pv(eltkid,ivtx).first;
                            d0err_     = elgsftkid>=0 ? gsftrks_d0_pv(elgsftkid,ivtx).second : trks_d0_pv(eltkid,ivtx).second;
                            dz_        = elgsftkid>=0 ? gsftrks_dz_pv(elgsftkid,ivtx).first  : trks_dz_pv(eltkid,ivtx).first;
                            dzerr_     = elgsftkid>=0 ? gsftrks_dz_pv(elgsftkid,ivtx).second : trks_dz_pv(eltkid,ivtx).second;
                        }
                        else
                        {
                            d0_        = cms2.els_d0().at(iLep);
                            d0err_     = cms2.els_d0Err().at(iLep);
                            dz_        = cms2.els_z0().at(iLep);
                            dzerr_     = cms2.els_z0Err().at(iLep);
                        }
                        ip3d_      = els_ip3d().at(iLep);;
                        ip3derr_   = els_ip3derr().at(iLep);;


                        pfmet_     = evt_pfmet();
                        pfmetphi_  = evt_pfmetPhi();
                        foel_mass_ = sqrt(fabs((lp4_ + foel_p4_).mass2()));
                        fomu_mass_ = sqrt(fabs((lp4_ + fomu_p4_).mass2()));

                        // Isolation
                        iso_          = electronIsolation_rel_v1      (iLep, /*use_calo_iso=*/true ); 
                        iso_nps_      = electronIsolation_rel_v1      (iLep, /*use_calo_iso=*/true ); 
                        trck_iso_     = electronIsolation_rel_v1      (iLep, /*use_calo_iso=*/false) * els_p4().at(iLep).pt(); 
                        ecal_iso_     = electronIsolation_ECAL_rel_v1 (iLep, /*use EBps=*/true     ) * els_p4().at(iLep).pt(); 
                        ecal_iso_nps_ = electronIsolation_ECAL_rel_v1 (iLep, /*use EBps=*/false    ) * els_p4().at(iLep).pt(); 
                        hcal_iso_     = electronIsolation_HCAL_rel    (iLep                        ) * els_p4().at(iLep).pt(); 

                        // PF Isolation
                        ch_pfiso03_ = cms2.els_iso03_pf2012ext_ch().at(iLep);
                        em_pfiso03_ = cms2.els_iso03_pf2012ext_em().at(iLep);
                        nh_pfiso03_ = cms2.els_iso03_pf2012ext_nh().at(iLep);
                        pfiso03_ = (ch_pfiso03_ + em_pfiso03_ + nh_pfiso03_)/els_p4().at(iLep).pt(); 

                        ch_pfiso04_ = cms2.els_iso04_pf2012ext_ch().at(iLep);
                        em_pfiso04_ = cms2.els_iso04_pf2012ext_em().at(iLep);
                        nh_pfiso04_ = cms2.els_iso04_pf2012ext_nh().at(iLep);
                        pfiso04_ = (ch_pfiso04_ + em_pfiso04_ + nh_pfiso04_)/els_p4().at(iLep).pt(); 

                        // correct isolaion (using fastjet effective area)
                        cpfiso03_rho_ = electronIsoValuePF2012_FastJetEffArea_v3(iLep, /*conesize=*/0.3, /*vtx=*/-999, /*52X iso=*/false);
                        cpfiso04_rho_ = electronIsoValuePF2012_FastJetEffArea_v3(iLep, /*conesize=*/0.4, /*vtx=*/-999, /*52X iso=*/false);

                        // mc information
                        if (!isData) {
                            mcid_       = els_mc_id().at(iLep);
                            mcmotherid_ = els_mc_motherid().at(iLep);
                            int status3_index = mc3idx_eormu(11, iLep);
                            if (status3_index >= 0)
                            {
                                mc3id_ = cms2.genps_id().at(status3_index);
                                mc3pt_ = cms2.genps_p4().at(status3_index).pt();                            
                                mc3p4_ = cms2.genps_p4().at(status3_index);                            
                            }
                            mc3dr_ = mc3dr_eormu(11, iLep);            
                            leptonIsFromW_ = leptonIsFromW(iLep, -11 * cms2.els_charge().at(iLep), true);
                        }

                        // ID
                        el_id_sieie_   = cms2.els_sigmaIEtaIEta().at(iLep);
                        el_id_detain_  = cms2.els_dEtaIn().at(iLep);
                        el_id_dphiin_  = cms2.els_dPhiIn().at(iLep);
                        el_id_smurfV5_ = pass_electronSelection( iLep, electronSelection_smurfV5_id );
                        el_id_vbtf80_  = electronId_VBTF(iLep, VBTF_35X_80, false, false);
                        el_id_vbtf90_  = electronId_VBTF(iLep, VBTF_35X_90, false, false);
                        if( els_closestMuon().at(iLep) == -1 )
                            closestMuon_ = true;

                        // electron ID effective area
                        el_effarea03_ = EffectiveArea(eta_, /*cone=*/0.3, /*eormu=*/11, /*use_tight=*/false);
                        el_effarea04_ = EffectiveArea(eta_, /*cone=*/0.4, /*eormu=*/11, /*use_tight=*/false);

                        // PV
                        d0PV_wwV1_ = electron_d0PV_wwV1(iLep);
                        dzPV_wwV1_ = electron_dzPV_wwV1(iLep);

                        // W transverse mass
                        mt_   = Mt( els_p4().at(iLep), pfmet_, pfmetphi_ );
                        pfmt_ = Mt( els_p4().at(iLep), pfmet_, pfmetphi_ );

                        // Do the 3 electron charges agree?
                        int iCTF = els_trkidx().at(iLep);
                        if( iCTF >= 0 ){
                            int qCTF = trks_charge().at( iCTF );
                            int qGSF = els_trk_charge().at(iLep);
                            int qPIX = els_sccharge().at(iLep);
                            if( qCTF == qGSF && qCTF == qPIX && qGSF == qPIX ) q3_ = true;
                        }

                        // Missing hits info
                        els_exp_innerlayers_ = els_exp_innerlayers().at(iLep);

                        // Conversion Rejection  
                        convHitPattern_   = isFromConversionHitPattern(iLep);
                        convPartnerTrack_ = isFromConversionPartnerTrack(iLep);
                        convMIT_          = isFromConversionMIT(iLep);
                        if( els_exp_innerlayers().at(iLep) == 0 ) conv0MissHits_ = true;

                        // HT
#ifndef __CMS2_SLIM__
                        ht_calo_           = (float) sumPt (iLep, JETS_TYPE_CALO_UNCORR  , JETS_CLEAN_SINGLE_E );
                        ht_calo_L2L3_      = (float) sumPt (iLep, JETS_TYPE_CALO_CORR    , JETS_CLEAN_SINGLE_E );
#endif
                        ht_pf_             = (float) sumPt (iLep, JETS_TYPE_PF_UNCORR    , JETS_CLEAN_SINGLE_E );
                        ht_pf_L2L3_        = (float) sumPt (iLep, JETS_TYPE_PF_CORR      , JETS_CLEAN_SINGLE_E );
                        ht_pf_L1FastL2L3_  = (float) sumPt (iLep, JETS_TYPE_PF_FAST_CORR , JETS_CLEAN_SINGLE_E );

                        //////////////////////////// 
                        // End Lepton Information //
                        ////////////////////////////

                        ///////////////////////  
                        // 2012 Triggers     //
                        ///////////////////////

                        // Electrons
                        triggerMatchStruct struct_ele8_CaloIdL_CaloIsoVL_vstar                                = MatchTriggerClass(els_p4().at(iLep), ele8_CaloIdL_CaloIsoVL_regexp                               );
                        triggerMatchStruct struct_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar               = MatchTriggerClass(els_p4().at(iLep), ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_regexp              );
                        triggerMatchStruct struct_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar         = MatchTriggerClass(els_p4().at(iLep), ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_regexp        );
                        triggerMatchStruct struct_ele8_CaloIdT_TrkIdVL_vstar                                  = MatchTriggerClass(els_p4().at(iLep), ele8_CaloIdT_TrkIdVL_regexp                                 );
                        triggerMatchStruct struct_ele8_CaloIdT_TrkIdVL_Jet30_vstar                            = MatchTriggerClass(els_p4().at(iLep), ele8_CaloIdT_TrkIdVL_Jet30_regexp                           );
                        triggerMatchStruct struct_ele17_CaloIdL_CaloIsoVL_vstar                               = MatchTriggerClass(els_p4().at(iLep), ele17_CaloIdL_CaloIsoVL_regexp                              );
                        triggerMatchStruct struct_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar              = MatchTriggerClass(els_p4().at(iLep), ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_regexp             );
                        triggerMatchStruct struct_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar        = MatchTriggerClass(els_p4().at(iLep), ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_regexp       );
                        triggerMatchStruct struct_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar = MatchTriggerClass(els_p4().at(iLep), ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_rexexp);
                        triggerMatchStruct struct_ele27_WP80_vstar                                            = MatchTriggerClass(els_p4().at(iLep), ele27_WP80_rexexp                                           );

                        ele8_CaloIdL_CaloIsoVL_vstar_                                      = struct_ele8_CaloIdL_CaloIsoVL_vstar.nHLTObjects_;
                        ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                     = struct_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar.nHLTObjects_;
                        ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_               = struct_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar.nHLTObjects_;
                        ele8_CaloIdT_TrkIdVL_vstar_                                        = struct_ele8_CaloIdT_TrkIdVL_vstar.nHLTObjects_;
                        ele8_CaloIdT_TrkIdVL_Jet30_vstar_                                  = struct_ele8_CaloIdT_TrkIdVL_Jet30_vstar.nHLTObjects_;
                        ele17_CaloIdL_CaloIsoVL_vstar_                                     = struct_ele17_CaloIdL_CaloIsoVL_vstar.nHLTObjects_;
                        ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                    = struct_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar.nHLTObjects_;
                        ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_              = struct_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar.nHLTObjects_;
                        ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_       = struct_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar.nHLTObjects_;
                        ele27_WP80_vstar_                                                  = struct_ele27_WP80_vstar.nHLTObjects_;

                        ele8_CaloIdL_CaloIsoVL_version_                                    = struct_ele8_CaloIdL_CaloIsoVL_vstar.version_;
                        ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_                   = struct_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar.version_;
                        ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_             = struct_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar.version_;
                        ele8_CaloIdT_TrkIdVL_version_                                      = struct_ele8_CaloIdT_TrkIdVL_vstar.version_;
                        ele8_CaloIdT_TrkIdVL_Jet30_version_                                = struct_ele8_CaloIdT_TrkIdVL_Jet30_vstar.version_;
                        ele17_CaloIdL_CaloIsoVL_version_                                   = struct_ele17_CaloIdL_CaloIsoVL_vstar.version_;
                        ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_version_                  = struct_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar.version_;
                        ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_version_            = struct_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar.version_;
                        ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_version_     = struct_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar.version_;
                        ele27_WP80_version_                                                = struct_ele27_WP80_vstar.version_;

                        dr_ele8_CaloIdL_CaloIsoVL_vstar_                                   = struct_ele8_CaloIdL_CaloIsoVL_vstar.dR_;
                        dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                  = struct_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar.dR_;
                        dr_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_            = struct_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar.dR_;
                        dr_ele8_CaloIdT_TrkIdVL_vstar_                                     = struct_ele8_CaloIdT_TrkIdVL_vstar.dR_;
                        dr_ele8_CaloIdT_TrkIdVL_Jet30_vstar_                               = struct_ele8_CaloIdT_TrkIdVL_Jet30_vstar.dR_;
                        dr_ele17_CaloIdL_CaloIsoVL_vstar_                                  = struct_ele17_CaloIdL_CaloIsoVL_vstar.dR_;
                        dr_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                 = struct_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar.dR_;
                        dr_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_           = struct_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar.dR_;
                        dr_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_    = struct_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar.dR_;
                        dr_ele27_WP80_vstar_                                               = struct_ele27_WP80_vstar.dR_;

                        hltps_ele8_CaloIdL_CaloIsoVL_vstar_                                = struct_ele8_CaloIdL_CaloIsoVL_vstar.hltps_;
                        hltps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_               = struct_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar.hltps_;
                        hltps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_         = struct_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar.hltps_;
                        hltps_ele8_CaloIdT_TrkIdVL_vstar_                                  = struct_ele8_CaloIdT_TrkIdVL_vstar.hltps_;
                        hltps_ele8_CaloIdT_TrkIdVL_Jet30_vstar_                            = struct_ele8_CaloIdT_TrkIdVL_Jet30_vstar.hltps_;
                        hltps_ele17_CaloIdL_CaloIsoVL_vstar_                               = struct_ele17_CaloIdL_CaloIsoVL_vstar.hltps_;
                        hltps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_              = struct_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar.hltps_;
                        hltps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_        = struct_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar.hltps_;
                        hltps_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_ = struct_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar.hltps_;
                        hltps_ele27_WP80_vstar_                                            = struct_ele27_WP80_vstar.hltps_;

                        // These are hardcoded to the value in Dima's table:
                        // http://dmytro.web.cern.ch/dmytro/trigger/triggerEvolution_all.html
#ifndef __CMS2_SLIM__
                        l1ps_ele8_CaloIdL_CaloIsoVL_vstar_                                    = L1_prescale("L1_SingleEG5");
                        l1ps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                   = L1_prescale("L1_SingleEG7");
                        l1ps_ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_             = L1_prescale("L1_SingleEG7");
                        l1ps_ele8_CaloIdT_TrkIdVL_vstar_                                      = L1_prescale("L1_SingleEG5");
                        l1ps_ele8_CaloIdT_TrkIdVL_Jet30_vstar_                                = L1_prescale("L1_SingleEG5");
                        l1ps_ele17_CaloIdL_CaloIsoVL_vstar_                                   = L1_prescale("L1_SingleEG12");
                        l1ps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_vstar_                  = L1_prescale("L1_SingleEG12");
                        l1ps_ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_vstar_            = L1_prescale("L1_SingleEG12");
                        l1ps_ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_vstar_     = L1_prescale("L1_SingleEG20");
                        l1ps_ele27_WP80_vstar_                                                = L1_prescale("L1_SingleEG20");
#endif


                        ///////////////////////  
                        // end 2012 Triggers //
                        ///////////////////////

                        ///////////////////////  
                        // 2011 Triggers     //
                        ///////////////////////


                        // Electrons
                        triggerMatchStruct struct_ele8_vstar                                          = MatchTriggerClass(els_p4().at(iLep), ele8_regexp                                         );
                        triggerMatchStruct struct_ele8_CaloIdL_TrkIdVL_vstar                          = MatchTriggerClass(els_p4().at(iLep), ele8_CaloIdL_TrkIdVL_regexp                         );
                        triggerMatchStruct struct_ele8_CaloIdL_CaloIsoVL_Jet40_vstar                  = MatchTriggerClass(els_p4().at(iLep), ele8_CaloIdL_CaloIsoVL_Jet40_regexp                 );
                        triggerMatchStruct struct_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar       = MatchTriggerClass(els_p4().at(iLep), ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_regexp);
                        triggerMatchStruct struct_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar = MatchTriggerClass(els_p4().at(iLep), photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_regexp);

                        ele8_vstar_                                             = struct_ele8_vstar.nHLTObjects_;
                        ele8_CaloIdL_TrkIdVL_vstar_                             = struct_ele8_CaloIdL_TrkIdVL_vstar.nHLTObjects_; 
                        ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                     = struct_ele8_CaloIdL_CaloIsoVL_Jet40_vstar.nHLTObjects_;
                        ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_          = struct_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar.nHLTObjects_;
                        photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_    = struct_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar.nHLTObjects_;

                        ele8_version_                                           = struct_ele8_vstar.version_;
                        ele8_CaloIdL_TrkIdVL_version_                           = struct_ele8_CaloIdL_TrkIdVL_vstar.version_; 
                        ele8_CaloIdL_CaloIsoVL_Jet40_version_                   = struct_ele8_CaloIdL_CaloIsoVL_Jet40_vstar.version_;
                        ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_version_        = struct_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar.version_;
                        photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_version_  = struct_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar.nHLTObjects_;

                        dr_ele8_vstar_                                          = struct_ele8_vstar.dR_;
                        dr_ele8_CaloIdL_TrkIdVL_vstar_                          = struct_ele8_CaloIdL_TrkIdVL_vstar.dR_;
                        dr_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                  = struct_ele8_CaloIdL_CaloIsoVL_Jet40_vstar.dR_;
                        dr_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_       = struct_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar.dR_;
                        dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_ = struct_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar.dR_; 

                        hltps_ele8_vstar_                                          = struct_ele8_vstar.hltps_;
                        hltps_ele8_CaloIdL_TrkIdVL_vstar_                          = struct_ele8_CaloIdL_TrkIdVL_vstar.hltps_;
                        hltps_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                  = struct_ele8_CaloIdL_CaloIsoVL_Jet40_vstar.hltps_;
                        hltps_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar_       = struct_ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_vstar.hltps_;
                        hltps_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_ = struct_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar.hltps_; 

                        ///////////////////////  
                        // end 2011 Triggers //
                        ///////////////////////

                        //////////////
                        // Jets     //
                        //////////////

                        // Calo Jets
                        // Find the highest Pt jet separated by at least dRcut from this lepton and fill the jet Pt
                        // #ifndef __CMS2_SLIM__
                        //                     ptj1_       = -999.0;
                        //                     ptj1_b2b_   = -999.0;
                        //                     dphij1_b2b_ = -999.0;
                        //                     nj1_        = 0;
                        //                     for (unsigned int iJet = 0; iJet < jets_p4().size(); iJet++) {
                        //                         double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jets_p4().at(iJet) );
                        //                         if( dr > deltaRCut && jets_p4().at(iJet).pt() > 10 ) nj1_++;
                        //                         if ( dr > deltaRCut && jets_p4().at(iJet).pt() > ptj1_ ){
                        //                             ptj1_ = jets_p4().at(iJet).pt();

                        //                             // back to back in phi
                        //                             float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iLep), jets_p4().at(iJet) ) );
                        //                             if( dphi > deltaPhiCut && jets_p4().at(iJet).pt() > ptj1_b2b_ ){ 
                        //                                 ptj1_b2b_   = jets_p4().at(iJet).pt();
                        //                                 dphij1_b2b_ = dphi;
                        //                             }
                        //                         }
                        //                     }
                        // #endif

                        // PF Jets
                        // Find the highest Pt pfjet separated by at least dRcut from this lepton and fill the pfjet Pt
                        ptpfj1_       = -999.0;
                        ptpfj1_b2b_   = -999.0;
                        dphipfj1_b2b_ = -999.0;
                        npfj1_        = 0;
                        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                            double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), pfjets_p4().at(iJet) );
                            if( dr > deltaRCut && pfjets_p4().at(iJet).pt() > 10 ) npfj1_++;
                            if ( dr > deltaRCut && pfjets_p4().at(iJet).pt() > ptpfj1_ ){
                                ptpfj1_ = pfjets_p4().at(iJet).pt();

                                // back to back in phi
                                float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iLep), pfjets_p4().at(iJet) ) );
                                if( dphi > deltaPhiCut && pfjets_p4().at(iJet).pt() > ptpfj1_b2b_ ){ 
                                    ptpfj1_b2b_   = pfjets_p4().at(iJet).pt();
                                    dphipfj1_b2b_ = dphi;
                                }
                            }
                        }

                        // L2L3 PF Jets
                        // Find the highest Pt PF L2L3 corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                        ptpfcj1_       = -999.0; 
                        ptpfcj1_b2b_   = -999.0;
                        dphipfcj1_b2b_ = -999.0;
                        npfcj1_        = 0;
                        btagpfc_       = false;
                        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                            if ( !passesPFJetID(iJet)) continue;
                            LorentzVector jp4 = pfjets_p4().at(iJet);
                            float jet_cor = jetCorrection(jp4, jet_pf_L2L3corrector);
                            //float jet_cor = pfjets_corL2L3().at(iJet);
                            LorentzVector jp4cor = jp4 * jet_cor;
                            if (jp4cor.pt() > 15 && pfjets_combinedSecondaryVertexBJetTag().at(iJet) > 0.679) btagpfc_ = true;
                            double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jp4cor );
                            if( dr > deltaRCut && jp4cor.pt() > 10 ) npfcj1_++;
                            if ( dr > deltaRCut && jp4cor.pt() > ptpfcj1_ ){
                                ptpfcj1_ = jp4cor.pt();

                                // back to back in phi
                                float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iLep), jp4cor ) );
                                if( dphi > deltaPhiCut && jp4cor.pt() > ptpfcj1_b2b_ ){
                                    ptpfcj1_b2b_   = jp4cor.pt();
                                    dphipfcj1_b2b_ = dphi;
                                } 
                            }
                        }

                        // L1FastL2L3 PF Jets
                        // Find the highest Pt PF L1FastL2L3 corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                        emfpfcL1Fj1_      = -999.0;
                        ptpfcL1Fj1_       = -999.0;
                        dphipfcL1Fj1_     = -999.0;
                        ptpfcL1Fj1_b2b_   = -999.0;
                        dphipfcL1Fj1_b2b_ = -999.0;
                        npfcL1Fj1_        = 0;
                        npfc30L1Fj1_      = 0;
                        npfc40L1Fj1_      = 0;
                        nbpfc40L1Fj1_     = 0;
                        npfc50L1Fj1_eth_  = 0;
                        npfc65L1Fj1_eth_  = 0;
                        btagpfcL1F_       = false;
                        rho_ = cms2.evt_rho();
                        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                            if ( !passesPFJetID(iJet)) continue;
                            LorentzVector jp4 = pfjets_p4().at(iJet);
                            float jet_cor = cms2.pfjets_corL1FastL2L3().at(iJet);
                            //jet_pf_L1FastJetL2L3_corrector->setRho(cms2.evt_ww_rho_vor());
                            //jet_pf_L1FastJetL2L3_corrector->setJetA(cms2.pfjets_area().at(iJet));
                            //jet_pf_L1FastJetL2L3_corrector->setJetPt(cms2.pfjets_p4().at(iJet).pt());
                            //jet_pf_L1FastJetL2L3_corrector->setJetEta(cms2.pfjets_p4().at(iJet).eta()); 
                            //float jet_cor = jetCorrection(jp4, jet_pf_L1FastJetL2L3_corrector);
                            LorentzVector jp4cor = jp4 * jet_cor;
                            if (jp4cor.pt() > 15 && pfjets_combinedSecondaryVertexBJetTag().at(iJet) > 0.679) btagpfcL1F_ = true;
                            double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jp4cor );
                            if( dr > deltaRCut && jp4cor.pt() > 10 ) npfcL1Fj1_++;
                            if( dr > deltaRCut && jp4cor.pt() > 30 ) npfc30L1Fj1_++;
                            if( dr > deltaRCut && jp4cor.pt() > 40 ) npfc40L1Fj1_++;
                            if( dr > deltaRCut && jp4cor.pt() > 40 && pfjets_combinedSecondaryVertexBJetTag().at(iJet) > 0.679) nbpfc40L1Fj1_++;
                            if (dr > 0.4       && jp4cor.pt() > 50 ) npfc50L1Fj1_eth_++;
                            if (dr > 0.4       && jp4cor.pt() > 65 ) npfc65L1Fj1_eth_++;
                            if ( dr > deltaRCut && jp4cor.pt() > ptpfcL1Fj1_ ){
                                emfpfcL1Fj1_ = (cms2.pfjets_chargedEmE().at(iJet) + cms2.pfjets_neutralEmE().at(iJet)) / pfjets_p4().at(iJet).E();
                                ptpfcL1Fj1_ = jp4cor.pt();
                                float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iLep), jp4cor ) );
                                dphipfcL1Fj1_ = dphi;

                                // back to back in phi
                                if( dphi > deltaPhiCut && jp4cor.pt() > ptpfcL1Fj1_b2b_ ){
                                    ptpfcL1Fj1_b2b_   = jp4cor.pt();
                                    dphipfcL1Fj1_b2b_ = dphi;
                                }
                            }
                        }

                        // L1FastL2L3Residual PF Jets
                        // Find the highest Pt PF L1FastL2L3Residual corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                        emfpfcL1Fj1res_      = -999.0;
                        ptpfcL1Fj1res_       = -999.0;
                        dphipfcL1Fj1res_      = -999.0;
                        ptpfcL1Fj1res_b2b_   = -999.0;
                        dphipfcL1Fj1res_b2b_ = -999.0;
                        npfcL1Fj1res_        = 0;
                        npfc30L1Fj1res_      = 0;
                        npfc40L1Fj1res_      = 0;
                        nbpfc40L1Fj1res_     = 0;
                        npfc50L1Fj1res_eth_  = 0;
                        npfc65L1Fj1res_eth_  = 0;
                        btagpfcL1Fres_       = false;
                        rho_ = cms2.evt_rho();
                        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                            if ( !passesPFJetID(iJet)) continue;
                            LorentzVector jp4 = pfjets_p4().at(iJet);
                            float jet_cor = cms2.pfjets_corL1FastL2L3residual().at(iJet);
                            //jet_pf_L1FastJetL2L3Residual_corrector->setRho(cms2.evt_ww_rho_vor());
                            //jet_pf_L1FastJetL2L3Residual_corrector->setJetA(cms2.pfjets_area().at(iJet));
                            //jet_pf_L1FastJetL2L3Residual_corrector->setJetPt(cms2.pfjets_p4().at(iJet).pt());
                            //jet_pf_L1FastJetL2L3Residual_corrector->setJetEta(cms2.pfjets_p4().at(iJet).eta()); 
                            //float jet_cor = jetCorrection(jp4, jet_pf_L1FastJetL2L3Residual_corrector);
                            LorentzVector jp4cor = jp4 * jet_cor;
                            if (jp4cor.pt() > 15 && pfjets_combinedSecondaryVertexBJetTag().at(iJet) > 0.679) btagpfcL1Fres_ = true;
                            double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jp4cor );
                            if( dr > deltaRCut && jp4cor.pt() > 10 ) npfcL1Fj1res_++;
                            if( dr > deltaRCut && jp4cor.pt() > 30 ) npfc30L1Fj1res_++;
                            if( dr > deltaRCut && jp4cor.pt() > 40 ) npfc40L1Fj1res_++;
                            if( dr > deltaRCut && jp4cor.pt() > 40 && pfjets_combinedSecondaryVertexBJetTag().at(iJet) > 0.679) nbpfc40L1Fj1res_++;
                            if (dr > 0.4       && jp4cor.pt() > 50 ) npfc50L1Fj1res_eth_++;
                            if (dr > 0.4       && jp4cor.pt() > 65 ) npfc65L1Fj1res_eth_++;
                            if ( dr > deltaRCut && jp4cor.pt() > ptpfcL1Fj1res_ ){
                                emfpfcL1Fj1res_ = (cms2.pfjets_chargedEmE().at(iJet) + cms2.pfjets_neutralEmE().at(iJet)) / pfjets_p4().at(iJet).E();
                                ptpfcL1Fj1res_ = jp4cor.pt();
                                float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iLep), jp4cor ) );
                                dphipfcL1Fj1res_ = dphi;

                                // back to back in phi
                                if( dphi > deltaPhiCut && jp4cor.pt() > ptpfcL1Fj1res_b2b_ ){
                                    ptpfcL1Fj1res_b2b_   = jp4cor.pt();
                                    dphipfcL1Fj1res_b2b_ = dphi;
                                }
                            }
                        }

                        // *** Doing B-tagging correctly ***
                        // B-tagged L1FastL2L3 PF Jets
                        // Find the highest Pt B-tagged PF L1FastL2L3 corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                        ptbtagpfcL1Fj1_       = -999.0;
                        dphibtagpfcL1Fj1_       = -999.0;
                        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                            if ( !passesPFJetID(iJet)) continue;
                            LorentzVector jp4 = pfjets_p4().at(iJet);
                            float jet_cor = cms2.pfjets_corL1FastL2L3().at(iJet);
                            //jet_pf_L1FastJetL2L3_corrector->setRho(cms2.evt_ww_rho_vor());
                            //jet_pf_L1FastJetL2L3_corrector->setJetA(cms2.pfjets_area().at(iJet));
                            //jet_pf_L1FastJetL2L3_corrector->setJetPt(cms2.pfjets_p4().at(iJet).pt());
                            //jet_pf_L1FastJetL2L3_corrector->setJetEta(cms2.pfjets_p4().at(iJet).eta()); 
                            //float jet_cor = jetCorrection(jp4, jet_pf_L1FastJetL2L3_corrector);
                            LorentzVector jp4cor = jp4 * jet_cor;
                            if (pfjets_combinedSecondaryVertexBJetTag().at(iJet) < 0.679) continue;
                            double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jp4cor );
                            if ( dr > deltaRCut && jp4cor.pt() > ptbtagpfcL1Fj1_ ){
                                ptbtagpfcL1Fj1_ = jp4cor.pt();
                                float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iLep), jp4cor ) );
                                dphibtagpfcL1Fj1_ = dphi; 
                            }
                        }

                        // *** Doing B-tagging correctly ***
                        // B-tagged L1FastL2L3Residual PF Jets
                        // Find the highest Pt B-tagged PF L1FastL2L3Residual corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                        ptbtagpfcL1Fj1res_       = -999.0;
                        dphibtagpfcL1Fj1res_       = -999.0;
                        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                            if ( !passesPFJetID(iJet)) continue;
                            LorentzVector jp4 = pfjets_p4().at(iJet);
                            float jet_cor = cms2.pfjets_corL1FastL2L3residual().at(iJet);
                            //jet_pf_L1FastJetL2L3Residual_corrector->setRho(cms2.evt_ww_rho_vor());
                            //jet_pf_L1FastJetL2L3Residual_corrector->setJetA(cms2.pfjets_area().at(iJet));
                            //jet_pf_L1FastJetL2L3Residual_corrector->setJetPt(cms2.pfjets_p4().at(iJet).pt());
                            //jet_pf_L1FastJetL2L3Residual_corrector->setJetEta(cms2.pfjets_p4().at(iJet).eta()); 
                            //float jet_cor = jetCorrection(jp4, jet_pf_L1FastJetL2L3Residual_corrector);
                            LorentzVector jp4cor = jp4 * jet_cor;
                            if (pfjets_combinedSecondaryVertexBJetTag().at(iJet) < 0.679) continue;
                            double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jp4cor );
                            if ( dr > deltaRCut && jp4cor.pt() > ptbtagpfcL1Fj1res_ ){
                                ptbtagpfcL1Fj1res_ = jp4cor.pt();
                                float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iLep), jp4cor ) );
                                dphibtagpfcL1Fj1res_ = dphi; 
                            }
                        }
                        //////////////
                        // End Jets //
                        //////////////


                        ///////////////////
                        // B Tagging     //
                        ///////////////////

                        // The btag information
                        // #ifndef __CMS2_SLIM__
                        //                     nbjet_ = this_nbjet;
                        //                     dRbNear_ = 99.;
                        //                     dRbFar_  = -99.;
                        //                     for (int ii=0; ii<nbjet_; ii++) {
                        //                         unsigned int iJet = bindex.at(ii);
                        //                         float dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jets_p4().at(iJet));
                        //                         if (dr < dRbNear_) dRbNear_ = dr;
                        //                         if (dr > dRbFar_)   dRbFar_  = dr;
                        //                     }
                        // #endif

                        // btag info for corrected pfjet
                        nbpfcjet_ = this_nbpfjet;
                        dRbpfcNear_ = 99.;
                        dRbpfcFar_  = -99.;
                        for (int ii=0; ii<nbpfcjet_; ii++) {
                            unsigned int iJet = bpfindex.at(ii);
                            LorentzVector jp4 = pfjets_p4().at(iJet);
                            float jet_cor = jetCorrection(jp4, jet_pf_L2L3corrector);
                            //float jet_corr = pfjets_corL2L3().at(iJet);
                            LorentzVector jp4cor = jp4 * jet_cor;
                            float dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jp4cor);
                            if (dr < dRbpfcNear_) dRbpfcNear_ = dr;
                            if (dr > dRbpfcFar_)   dRbpfcFar_  = dr;
                        }

                        ///////////////////
                        // End B Tagging //
                        ///////////////////



                        // Time to fill the baby for the electrons
                        FillBabyNtuple();

                    } // closes loop over electrons
                } // closes if statements about whether we want to fill electrons


                // Muons
                if (eormu == -1 || eormu==13) {
                    for ( unsigned int iLep = 0; iLep < mus_p4().size(); iLep++) {

                        // Apply a pt cut
                        if ( mus_p4().at(iLep).pt() < 5.0) continue;

                        ////////////////////////////////////////////////////////////////////////
                        // NEED TO THINK ABOUT THIS... Z'S ARE VETOED BASED ON TOP SELECTIONS //
                        ////////////////////////////////////////////////////////////////////////

                        // If it is above 20 GeV see if we can make a 
                        // Z with another pt>20 FO.  
                        bool isaZ = false;
                        if (mus_p4().at(iLep).pt() > 20.) {
                            for (unsigned int jMu = 0 ; jMu < mus_p4().size(); jMu++) {
                                if (iLep == jMu)                             continue;
                                if (mus_p4().at(jMu).pt() < 20.)            continue;
                                if ( ! muonId( jMu,  OSGeneric_v3_FO) ) continue;
                                if ( ! muonId( iLep, OSGeneric_v3_FO) ) continue;
                                LorentzVector w = mus_p4().at(iLep) + mus_p4().at(jMu);
                                if (abs(w.mass()-91.) > 20.) continue;
                                isaZ = true;
                            }
                        }
                        if (isaZ) continue;

                        // Initialize baby ntuple
                        InitBabyNtuple();

                        // store number of electron FOs in event (use SS FO definition)
                        nFOels_ = 0;
                        ngsfs_ = 0;
                        for (unsigned int iel = 0; iel < cms2.els_p4().size(); iel++) {
                            if (cms2.els_p4().at(iel).pt() < 10.)
                                continue;

                            if (pass_electronSelection(iel, electronSelectionFOV7_v3, false, false))
                                ++ngsfs_;

                            if (samesign::isDenominatorLepton(11, iel)) {
                                ++nFOels_;
                                if (cms2.els_p4().at(iel).pt() > foel_p4_.pt()) {
                                    foel_p4_ = cms2.els_p4().at(iel);
                                    foel_id_ = 11*cms2.els_charge().at(iel);
                                }
                                continue;
                            }
                            if (samesign::isDenominatorLepton(11, iel)) {
                                ++nFOels_;
                                if (cms2.els_p4().at(iel).pt() > foel_p4_.pt()) {
                                    foel_p4_ = cms2.els_p4().at(iel);
                                    foel_id_ = 11*cms2.els_charge().at(iel);
                                }
                                continue;
                            }
                        }

                        // store number of muon FOs in event (use SS FO definition)
                        nFOmus_ = 0;
                        nmus_ = 0;
                        for (unsigned int imu = 0; imu < cms2.mus_p4().size(); imu++) {
                            if (imu == iLep)
                                continue;

                            if (cms2.mus_p4().at(imu).pt() < 10.)
                                continue;

                            if (muonIdNotIsolated(imu, muonSelectionFO_ssV5))
                                ++nmus_;

                            if (samesign::isDenominatorLepton(13, imu)) {
                                ++nFOmus_;
                                if (cms2.mus_p4().at(imu).pt() > fomu_p4_.pt() && imu != iLep) {
                                    fomu_p4_ = cms2.mus_p4().at(imu);
                                    fomu_id_ = 13*cms2.mus_charge().at(imu);
                                }
                                continue;
                            }
                            if (samesign::isDenominatorLepton(13, imu)) {
                                ++nFOmus_;
                                if (cms2.mus_p4().at(imu).pt() > fomu_p4_.pt() && imu != iLep) {
                                    fomu_p4_ = cms2.mus_p4().at(imu);
                                    fomu_id_ = 13*cms2.mus_charge().at(imu);
                                }
                                continue;
                            }
                        }

                        // store number of "veto" electrons in event
                        nvetoels_ = 0;
                        for (unsigned int iel = 0; iel < cms2.els_p4().size(); iel++)
                        {
                            if (cms2.els_p4().at(iel).pt() < 5.0f)
                                continue;

                            if (fabs(cms2.els_p4().at(iel).eta()) > 2.4f)
                                continue;

                            const float iso = electronIsoValuePF2012_FastJetEffArea_v3(iel, /*conesize=*/0.3, /*vtx=*/-999, /*52X iso=*/false);
                            if (iso > 1.0)
                                continue;

                            // if we get here, then count it
                            nvetoels_++;
                        }

                        // store number of "veto" muons in event
                        nvetomus_ = 0;
                        for (unsigned int imu = 0; imu < cms2.mus_p4().size(); imu++)
                        {
                            if (imu == iLep) // skip the current muon
                                continue;

                            if (cms2.mus_p4().at(imu).pt() < 5.0f)
                                continue;

                            if (fabs(cms2.mus_p4().at(imu).eta()) > 2.4f)
                                continue;

                            const bool is_global      = ((cms2.mus_type().at(imu) & (1<<1)) != 0);
                            const bool is_tracker     = ((cms2.mus_type().at(imu) & (1<<2)) != 0);
                            const bool is_pfmu        = ((cms2.mus_type().at(imu) & (1<<5)) != 0);
                            const bool passes_mu_type = ((is_global or is_tracker) and is_pfmu);
                            if (not passes_mu_type)
                                continue;

                            const float iso = muonIsoValuePF2012_deltaBeta(imu); 
                            if (iso > 1.0)
                                continue;

                            // if we get here, then count it
                            nvetomus_++;
                        }

                        ////////////////////////////////////////////////////////////////////////
                        // STORE SOME Z MASS VARIABLES //
                        ////////////////////////////////////////////////////////////////////////

                        mz_fo_ctf_  = -999.;
                        mz_ctf_iso_ = -999.;
                        mupsilon_fo_mu_ = -999.;
                        mupsilon_mu_iso_ = -999.;
                        LorentzVector p4fo = cms2.mus_p4().at(iLep);
                        for (unsigned int imu = 0; imu < cms2.mus_p4().size(); imu++) {
                            if (imu == iLep) continue;

                            if (fabs(cms2.mus_p4().at(imu).eta()) > 2.5)
                                continue;

                            if (cms2.mus_p4().at(imu).pt() < 10.)
                                continue;

                            LorentzVector zp4 = p4fo + cms2.mus_p4().at(imu);
                            float zcandmass = sqrt(fabs(zp4.mass2()));
                            if ( fabs(zcandmass - 91.) < fabs(mz_fo_ctf_ - 91.) ) {
                                mz_fo_ctf_  = zcandmass;
                                mz_ctf_iso_ = muonIsoValue(imu, false);
                            }
                            if ( fabs(zcandmass - 9.5) < fabs(mupsilon_fo_mu_ - 9.5) ) {
                                mupsilon_fo_mu_  = zcandmass;
                                mupsilon_mu_iso_ = muonIsoValue(imu, false); 
                            }
                        }


                        /////////////////////////// 
                        // Event Information     //
                        ///////////////////////////

                        // Load the electron and event quantities
                        run_          = evt_run();
                        ls_           = evt_lumiBlock();
                        evt_          = evt_event();
                        weight_       = isData ? 1.0 : evt_scale1fb();
                        filename_     = f->GetName();
                        dataset_      = evt_dataset().front();
                        is_real_data_ = evt_isRealData();

                        if(!isData){
                            // Pileup - PUSummaryInfoMaker
                            for (unsigned int vidx = 0; vidx < cms2.puInfo_nPUvertices().size(); vidx++) {
                                if (cms2.puInfo_bunchCrossing().at(vidx) != 0)
                                    continue;
                                pu_nPUvertices_ = cms2.puInfo_nPUvertices().at(vidx);
                                pu_nPUtrueint_  = cms2.puInfo_trueNumInteractions().at(vidx);
                            }
                        } 

                        // Pileup - VertexMaker
                        bool first_good_vertex_found         = false;
                        unsigned int first_good_vertex_index = 0;
                        for (unsigned int vidx = 0; vidx < cms2.vtxs_position().size(); vidx++)
                        {
                            if (!isGoodVertex(vidx))
                            {
                                continue;
                            }
                            if (!first_good_vertex_found)
                            {
                                first_good_vertex_found = true;
                                first_good_vertex_index = vidx;
                            }
                            ++evt_nvtxs_;
                        }

                        /////////////////////////// 
                        // End Event Information //
                        ///////////////////////////



                        //////////////////////////// 
                        // Lepton Information     //
                        ////////////////////////////

                        // Basic Quantities
                        lp4_       = cms2.mus_p4().at(iLep);
                        pt_        = mus_p4().at(iLep).pt();
                        eta_       = mus_p4().at(iLep).eta();
                        phi_       = mus_p4().at(iLep).phi();
                        id_        = 13*mus_charge().at(iLep);
                        pfmet_     = evt_pfmet();
                        pfmetphi_  = evt_pfmetPhi();
                        foel_mass_ = sqrt(fabs((lp4_ + foel_p4_).mass2()));
                        fomu_mass_ = sqrt(fabs((lp4_ + fomu_p4_).mass2()));

                        // ip (2d and 3d)
                        const int mutkid = cms2.mus_trkidx().at(iLep);
                        const int ivtx   = firstGoodVertex();
                        if (ivtx >= 0 && mutkid >= 0) 
                        {
                            d0_        = trks_d0_pv(mutkid,ivtx).first;
                            d0err_     = trks_d0_pv(mutkid,ivtx).second;
                            dz_        = trks_dz_pv(mutkid,ivtx).first;
                            dzerr_     = trks_dz_pv(mutkid,ivtx).second;
                        }
                        else
                        {
                            d0_        = cms2.mus_d0().at(iLep);
                            d0err_     = cms2.mus_d0Err().at(iLep);
                            dz_        = cms2.mus_z0().at(iLep);
                            dzerr_     = cms2.mus_z0Err().at(iLep);
                        }
                        ip3d_      = mus_ip3d().at(iLep);;
                        ip3derr_   = mus_ip3derr().at(iLep);;

                        // Isolation
                        iso_      = muonIsoValue     (iLep, /*truncated=*/false);
                        trck_iso_ = muonIsoValue_TRK (iLep, /*truncated=*/false) * mus_p4().at(iLep).pt();
                        ecal_iso_ = muonIsoValue_ECAL(iLep, /*truncated=*/false) * mus_p4().at(iLep).pt();
                        hcal_iso_ = muonIsoValue_HCAL(iLep, /*truncated=*/false) * mus_p4().at(iLep).pt();

                        // PF Isolation 03
                        ch_pfiso03_ = mus_isoR03_pf_ChargedHadronPt().at(iLep);
                        nh_pfiso03_ = mus_isoR03_pf_NeutralHadronEt().at(iLep);
                        em_pfiso03_ = mus_isoR03_pf_PhotonEt().at(iLep);
                        pfiso03_    = (ch_pfiso03_ + em_pfiso03_ + nh_pfiso03_)/mus_p4().at(iLep).pt(); 

                        // PF Isolation 04
                        ch_pfiso04_ = mus_isoR04_pf_ChargedHadronPt().at(iLep);
                        nh_pfiso04_ = mus_isoR04_pf_NeutralHadronEt().at(iLep);
                        em_pfiso04_ = mus_isoR04_pf_PhotonEt().at(iLep);
                        pfiso04_    = (ch_pfiso04_ + em_pfiso04_ + nh_pfiso04_)/mus_p4().at(iLep).pt(); 

                        // Radial Isolation 03
                        //radiso_et1p0_ = muonRadialIsolation(iLep, ch_radiso_et1p0_, nh_radiso_et1p0_, em_radiso_et1p0_, /*neutral_et_threshold=*/1.0, /*cone size=*/0.3, verbose_); 
                        //radiso_et0p5_ = muonRadialIsolation(iLep, ch_radiso_et0p5_, nh_radiso_et0p5_, em_radiso_et0p5_, /*neutral_et_threshold=*/0.5, /*cone size=*/0.3, verbose_); 

                        // PF Pile UP Sim pT
                        pfpupt03_ = mus_isoR03_pf_PUPt().at(iLep);
                        pfpupt04_ = mus_isoR04_pf_PUPt().at(iLep);

                        // correct isolaion (for SS2012)
                        cpfiso03_db_ = muonIsoValuePF2012_deltaBeta(iLep); 

                        // mc information
                        if (!isData) {
                            mcid_       = mus_mc_id().at(iLep);
                            mcmotherid_ = mus_mc_motherid().at(iLep);
                            int status3_index = mc3idx_eormu(13, iLep);
                            if (status3_index >= 0)
                            {
                                mc3id_ = cms2.genps_id().at(status3_index);
                                mc3pt_ = cms2.genps_p4().at(status3_index).pt();                            
                                mc3p4_ = cms2.genps_p4().at(status3_index);                            
                            }
                            mc3dr_ = mc3dr_eormu(13, iLep);            
                            leptonIsFromW_ = leptonIsFromW(iLep, -13 * cms2.mus_charge().at(iLep), true);
                        }

                        // muon effective area
                        // 2012 working point effective id (take from https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=188494)
                        mu_effarea03_          = EffectiveArea   (eta_, /*cone=*/0.3, /*eormu=*/13, /*use_tight=*/false);
                        mu_nh_effarea03_       = EffectiveArea_nh(eta_, /*cone=*/0.3, /*use_tight=*/false);
                        mu_em_effarea03_       = EffectiveArea_em(eta_, /*cone=*/0.3, /*use_tight=*/false);
                        mu_effarea03_tight_    = EffectiveArea   (eta_, /*cone=*/0.3, /*eormu=*/13, /*use_tight=*/true);
                        mu_nh_effarea03_tight_ = EffectiveArea_nh(eta_, /*cone=*/0.3, /*use_tight=*/true);
                        mu_em_effarea03_tight_ = EffectiveArea_em(eta_, /*cone=*/0.3, /*use_tight=*/true);
                        mu_effarea04_          = EffectiveArea   (eta_, /*cone=*/0.4, /*eormu=*/13, /*use_tight=*/false);
                        mu_nh_effarea04_       = EffectiveArea_nh(eta_, /*cone=*/0.4, /*use_tight=*/false);
                        mu_em_effarea04_       = EffectiveArea_em(eta_, /*cone=*/0.4, /*use_tight=*/false);
                        mu_effarea04_tight_    = EffectiveArea   (eta_, /*cone=*/0.4, /*eormu=*/13, /*use_tight=*/true);
                        mu_nh_effarea04_tight_ = EffectiveArea_nh(eta_, /*cone=*/0.4, /*use_tight=*/true);
                        mu_em_effarea04_tight_ = EffectiveArea_em(eta_, /*cone=*/0.4, /*use_tight=*/true);

                        // cosmics
                        mu_isCosmic_           = isCosmics(iLep);

                        // some muon quantities
                        mu_ecal_veto_dep_ = cms2.mus_iso_ecalvetoDep().at(iLep);
                        mu_hcal_veto_dep_ = cms2.mus_iso_hcalvetoDep().at(iLep);
                        mu_nchi2_         = cms2.mus_chi2().at(iLep) / cms2.mus_ndof().at(iLep);

                        // W transverse mass
                        mt_   = Mt( mus_p4().at(iLep), pfmet_, pfmetphi_ );
                        pfmt_ = Mt( mus_p4().at(iLep), pfmet_, pfmetphi_ );

                        // HT
#ifndef __CMS2_SLIM__
                        ht_calo_           = (float) sumPt (iLep, JETS_TYPE_CALO_UNCORR  , JETS_CLEAN_SINGLE_MU );
                        ht_calo_L2L3_      = (float) sumPt (iLep, JETS_TYPE_CALO_CORR    , JETS_CLEAN_SINGLE_MU );
#endif
                        ht_pf_             = (float) sumPt (iLep, JETS_TYPE_PF_UNCORR    , JETS_CLEAN_SINGLE_MU );
                        ht_pf_L2L3_        = (float) sumPt (iLep, JETS_TYPE_PF_CORR      , JETS_CLEAN_SINGLE_MU );

                        //////////////////////////// 
                        // End Lepton Information //
                        ////////////////////////////

                        //////////////////////////////////////////////////////
                        // Fake Rate Numerator & Denominator Selections     //
                        //////////////////////////////////////////////////////

                        //////////
                        // 2012 //
                        //////////

                        // SS
                        num_mu_ssV5_       = muonId(iLep, NominalSSv5);
                        num_mu_ssV5_noIso_ = muonIdNotIsolated(iLep, NominalSSv5);
                        fo_mu_ssV5_        = muonId(iLep, muonSelectionFO_ssV5);
                        fo_mu_ssV5_noIso_  = muonIdNotIsolated(iLep, muonSelectionFO_ssV5);

                        // TTZ
                        num_mu_TTZtightv1_       = ttv::isNumeratorLepton(13, iLep, ttv::LeptonType::TIGHT);
                        num_mu_TTZtightv1_noIso_ = ttv::isGoodLepton(13, iLep, ttv::LeptonType::TIGHT);
                        fo_mu_TTZtightv1_        = muonId(iLep, NominalTTZ_tightFO_v1);
                        fo_mu_TTZtightv1_noIso_  = ttv::isDenominatorLepton(13, iLep, ttv::LeptonType::TIGHT);

                        num_mu_TTZloosev1_       = ttv::isNumeratorLepton(13, iLep, ttv::LeptonType::LOOSE);
                        num_mu_TTZloosev1_noIso_ = ttv::isGoodLepton(13, iLep, ttv::LeptonType::LOOSE);
                        fo_mu_TTZloosev1_        = muonId(iLep, NominalTTZ_looseFO_v1);
                        fo_mu_TTZloosev1_noIso_  = ttv::isDenominatorLepton(13, iLep, ttv::LeptonType::LOOSE);



                        //////////
                        // 2011 //
                        //////////

                        // SS
                        numNomSSv4_      = muonId(iLep, NominalSSv4          );
                        fo_mussV4_04_    = muonId(iLep, muonSelectionFO_ssV4 );
                        numNomSSv4noIso_ = muonIdNotIsolated(iLep, NominalSSv4);
                        fo_mussV4_noIso_ = muonIdNotIsolated(iLep, muonSelectionFO_ssV4 );

                        //OS
                        num_mu_OSGV3_     = muonId(iLep, OSGeneric_v3);
                        fo_mu_OSGV3_      = muonId(iLep, OSGeneric_v3_FO);

                        // WW
                        num_mu_smurfV6_    = muonId(iLep, NominalSmurfV6                    );
                        fo_mu_smurf_04_    = muonId(iLep, muonSelectionFO_mu_smurf_04       );
                        fo_mu_smurf_10_    = muonId(iLep, muonSelectionFO_mu_smurf_10       );

                        ////////////////////////////////////////////////////////////
                        // Skip this muon if it fails the loosest denominator.    //
                        // Ignore this OR if applyFOfilter is set to false.       //
                        ////////////////////////////////////////////////////////////
                        if (applyFOfilter) {
                            if (
                                !fo_mussV4_04_     && !fo_mu_ssV5_       &&                    // SS
                                !fo_mu_TTZtightv1_ && !fo_mu_TTZloosev1_ &&                 // TTZ 2012
                                !fo_mu_OSGV2_      && !fo_mu_OSGV3_      &&                    // OS
                                !fo_mu_smurf_04_   && !fo_mu_smurf_10_                       // WW
                                )
                                continue;
                        }

                        //////////////////////////////////////////////////////
                        // End Fake Rate Numerator & Denominator Selections //
                        //////////////////////////////////////////////////////

                        ///////////////////////  
                        // 2012 Triggers     //
                        ///////////////////////

                        // Muons
                        triggerMatchStruct struct_mu5_vstar             = MatchTriggerClass(mus_p4().at(iLep), mu5_regexp            , 13);
                        triggerMatchStruct struct_mu8_vstar             = MatchTriggerClass(mus_p4().at(iLep), mu8_regexp            , 13);
                        triggerMatchStruct struct_mu12_vstar            = MatchTriggerClass(mus_p4().at(iLep), mu12_regexp           , 13);
                        triggerMatchStruct struct_mu17_vstar            = MatchTriggerClass(mus_p4().at(iLep), mu17_regexp           , 13);
                        triggerMatchStruct struct_mu15_eta2p1_vstar     = MatchTriggerClass(mus_p4().at(iLep), mu15_eta2p1_regexp    , 13);
                        triggerMatchStruct struct_mu24_eta2p1_vstar     = MatchTriggerClass(mus_p4().at(iLep), mu24_eta2p1_regexp    , 13);
                        triggerMatchStruct struct_mu30_eta2p1_vstar     = MatchTriggerClass(mus_p4().at(iLep), mu30_eta2p1_regexp    , 13);
                        triggerMatchStruct struct_isoMu20_eta2p1_vstar  = MatchTriggerClass(mus_p4().at(iLep), isoMu20_eta2p1_regexp , 13);
                        triggerMatchStruct struct_isoMu24_eta2p1_vstar  = MatchTriggerClass(mus_p4().at(iLep), isoMu24_eta2p1_regexp , 13);
                        triggerMatchStruct struct_isoMu30_eta2p1_vstar  = MatchTriggerClass(mus_p4().at(iLep), isoMu30_eta2p1_regexp , 13);
                        triggerMatchStruct struct_relIso1p0Mu17_vstar   = MatchTriggerClass(mus_p4().at(iLep), relIso1p0Mu17_regexp  , 13);
                        triggerMatchStruct struct_relIso1p0Mu20_vstar   = MatchTriggerClass(mus_p4().at(iLep), relIso1p0Mu20_regexp  , 13);
                        triggerMatchStruct struct_relIso1p0Mu5_vstar    = MatchTriggerClass(mus_p4().at(iLep), relIso1p0Mu5_regexp   , 13);

                        mu5_vstar_                  = struct_mu5_vstar.nHLTObjects_;
                        mu8_vstar_                  = struct_mu8_vstar.nHLTObjects_;
                        mu12_vstar_                 = struct_mu12_vstar.nHLTObjects_;
                        mu17_vstar_                 = struct_mu17_vstar.nHLTObjects_;
                        mu15_eta2p1_vstar_          = struct_mu15_eta2p1_vstar.nHLTObjects_;
                        mu24_eta2p1_vstar_          = struct_mu24_eta2p1_vstar.nHLTObjects_;
                        mu30_eta2p1_vstar_          = struct_mu30_eta2p1_vstar.nHLTObjects_;
                        isoMu20_eta2p1_vstar_       = struct_isoMu20_eta2p1_vstar.nHLTObjects_;
                        isoMu24_eta2p1_vstar_       = struct_isoMu24_eta2p1_vstar.nHLTObjects_;
                        isoMu30_eta2p1_vstar_       = struct_isoMu30_eta2p1_vstar.nHLTObjects_;
                        relIso1p0Mu17_vstar_        = struct_relIso1p0Mu17_vstar.nHLTObjects_;
                        relIso1p0Mu20_vstar_        = struct_relIso1p0Mu20_vstar.nHLTObjects_;
                        relIso1p0Mu5_vstar_         = struct_relIso1p0Mu5_vstar.nHLTObjects_;

                        mu5_version_                = struct_mu5_vstar.version_;
                        mu8_version_                = struct_mu8_vstar.version_;
                        mu12_version_               = struct_mu12_vstar.version_;
                        mu17_version_               = struct_mu17_vstar.version_;
                        mu15_eta2p1_version_        = struct_mu15_eta2p1_vstar.version_;
                        mu24_eta2p1_version_        = struct_mu24_eta2p1_vstar.version_;
                        mu30_eta2p1_version_        = struct_mu30_eta2p1_vstar.version_;
                        isoMu20_eta2p1_version_     = struct_isoMu20_eta2p1_vstar.version_;
                        isoMu24_eta2p1_version_     = struct_isoMu24_eta2p1_vstar.version_;
                        isoMu30_eta2p1_version_     = struct_isoMu30_eta2p1_vstar.version_;
                        relIso1p0Mu17_version_      = struct_relIso1p0Mu17_vstar.version_;
                        relIso1p0Mu20_version_      = struct_relIso1p0Mu20_vstar.version_;
                        relIso1p0Mu5_version_       = struct_relIso1p0Mu5_vstar.version_;

                        dr_mu5_vstar_               = struct_mu5_vstar.dR_;
                        dr_mu8_vstar_               = struct_mu8_vstar.dR_;
                        dr_mu12_vstar_              = struct_mu12_vstar.dR_;
                        dr_mu17_vstar_              = struct_mu17_vstar.dR_;
                        dr_mu15_eta2p1_vstar_       = struct_mu15_eta2p1_vstar.dR_;
                        dr_mu24_eta2p1_vstar_       = struct_mu24_eta2p1_vstar.dR_;
                        dr_mu30_eta2p1_vstar_       = struct_mu30_eta2p1_vstar.dR_;
                        dr_isoMu20_eta2p1_vstar_    = struct_isoMu20_eta2p1_vstar.dR_;
                        dr_isoMu24_eta2p1_vstar_    = struct_isoMu24_eta2p1_vstar.dR_;
                        dr_isoMu30_eta2p1_vstar_    = struct_isoMu30_eta2p1_vstar.dR_;
                        dr_relIso1p0Mu17_vstar_     = struct_relIso1p0Mu17_vstar.dR_;
                        dr_relIso1p0Mu20_vstar_     = struct_relIso1p0Mu20_vstar.dR_;
                        dr_relIso1p0Mu5_vstar_      = struct_relIso1p0Mu5_vstar.dR_;

                        hltps_mu5_vstar_            = struct_mu5_vstar.hltps_;
                        hltps_mu8_vstar_            = struct_mu8_vstar.hltps_;
                        hltps_mu12_vstar_           = struct_mu12_vstar.hltps_;
                        hltps_mu17_vstar_           = struct_mu17_vstar.hltps_;
                        hltps_mu15_eta2p1_vstar_    = struct_mu15_eta2p1_vstar.hltps_;
                        hltps_mu24_eta2p1_vstar_    = struct_mu24_eta2p1_vstar.hltps_;
                        hltps_mu30_eta2p1_vstar_    = struct_mu30_eta2p1_vstar.hltps_;
                        hltps_isoMu20_eta2p1_vstar_ = struct_isoMu20_eta2p1_vstar.hltps_;
                        hltps_isoMu24_eta2p1_vstar_ = struct_isoMu24_eta2p1_vstar.hltps_;
                        hltps_isoMu30_eta2p1_vstar_ = struct_isoMu30_eta2p1_vstar.hltps_;
                        hltps_relIso1p0Mu17_vstar_  = struct_relIso1p0Mu17_vstar.hltps_;
                        hltps_relIso1p0Mu20_vstar_  = struct_relIso1p0Mu20_vstar.hltps_;
                        hltps_relIso1p0Mu5_vstar_   = struct_relIso1p0Mu5_vstar.hltps_;

                        // These are hardcoded to the value in Dima's table:
                        // http://dmytro.web.cern.ch/dmytro/trigger/triggerEvolution_all.html 
#ifndef __CMS2_SLIM__
                        l1ps_mu5_vstar_            = L1_prescale("L1_SingleMu3"         );
                        l1ps_mu8_vstar_            = L1_prescale("L1_SingleMu3"         );
                        l1ps_mu12_vstar_           = L1_prescale("L1_SingleMu7"         );
                        l1ps_mu17_vstar_           = L1_prescale("L1_SingleMu12"        );
                        l1ps_mu15_eta2p1_vstar_    = L1_prescale("L1_SingleMu7"         );
                        l1ps_mu24_eta2p1_vstar_    = L1_prescale("L1_SingleMu16_Eta2p1" );
                        l1ps_mu30_eta2p1_vstar_    = L1_prescale("L1_SingleMu16_Eta2p1" );
                        l1ps_isoMu20_eta2p1_vstar_ = L1_prescale("L1_SingleMu16_Eta2p1" );
                        l1ps_isoMu24_eta2p1_vstar_ = L1_prescale("L1_SingleMu16_Eta2p1" );
                        l1ps_isoMu30_eta2p1_vstar_ = L1_prescale("L1_SingleMu16_Eta2p1" );
                        l1ps_relIso1p0Mu17_vstar_  = L1_prescale("L1_SingleMu12"        );
                        l1ps_relIso1p0Mu5_vstar_   = L1_prescale("L1_SingleMu3"         );
#endif


                        ///////////////////////  
                        // End 2012 Triggers //
                        ///////////////////////

                        ///////////////////////  
                        // 2011 Triggers     //
                        ///////////////////////

                        // Muons
                        triggerMatchStruct struct_mu3_vstar       = MatchTriggerClass(mus_p4().at(iLep), mu3_regexp      , 13);
                        triggerMatchStruct struct_mu15_vstar      = MatchTriggerClass(mus_p4().at(iLep), mu15_regexp     , 13);
                        triggerMatchStruct struct_mu20_vstar      = MatchTriggerClass(mus_p4().at(iLep), mu20_regexp     , 13);
                        triggerMatchStruct struct_mu24_vstar      = MatchTriggerClass(mus_p4().at(iLep), mu24_regexp     , 13);
                        triggerMatchStruct struct_mu30_vstar      = MatchTriggerClass(mus_p4().at(iLep), mu30_regexp     , 13);
                        triggerMatchStruct struct_mu8_Jet40_vstar = MatchTriggerClass(mus_p4().at(iLep), mu8_Jet40_regexp, 13);

                        mu3_vstar_          = struct_mu3_vstar.nHLTObjects_;
                        mu15_vstar_         = struct_mu15_vstar.nHLTObjects_;
                        mu20_vstar_         = struct_mu20_vstar.nHLTObjects_;
                        mu24_vstar_         = struct_mu24_vstar.nHLTObjects_;
                        mu30_vstar_         = struct_mu30_vstar.nHLTObjects_;
                        mu8_Jet40_vstar_    = struct_mu8_Jet40_vstar.nHLTObjects_;

                        mu3_version_        = struct_mu3_vstar.version_;
                        mu15_version_       = struct_mu15_vstar.version_;
                        mu20_version_       = struct_mu20_vstar.version_;
                        mu24_version_       = struct_mu24_vstar.version_;
                        mu30_version_       = struct_mu30_vstar.version_;
                        mu8_Jet40_version_  = struct_mu8_Jet40_vstar.version_;

                        dr_mu3_vstar_       = struct_mu3_vstar.dR_;
                        dr_mu15_vstar_      = struct_mu15_vstar.dR_;
                        dr_mu20_vstar_      = struct_mu20_vstar.dR_;
                        dr_mu24_vstar_      = struct_mu24_vstar.dR_;
                        dr_mu30_vstar_      = struct_mu30_vstar.dR_;
                        dr_mu8_Jet40_vstar_ = struct_mu8_Jet40_vstar.dR_;

                        hltps_mu3_vstar_       = struct_mu3_vstar.hltps_;
                        hltps_mu15_vstar_      = struct_mu15_vstar.hltps_;
                        hltps_mu20_vstar_      = struct_mu20_vstar.hltps_;
                        hltps_mu24_vstar_      = struct_mu24_vstar.hltps_;
                        hltps_mu30_vstar_      = struct_mu30_vstar.hltps_;
                        hltps_mu8_Jet40_vstar_ = struct_mu8_Jet40_vstar.hltps_;

                        ///////////////////////  
                        // End 2011 Triggers //
                        ///////////////////////

                        //////////////
                        // Jets     //
                        //////////////

                        // Calo Jets
                        // Find the highest Pt jet separated by at least dRcut from this lepton and fill the jet Pt
                        // #ifndef __CMS2_SLIM__
                        //                     ptj1_       = -999.0;
                        //                     ptj1_b2b_   = -999.0;
                        //                     dphij1_b2b_ = -999.0;
                        //                     nj1_        = 0;
                        //                     for (unsigned int iJet = 0; iJet < jets_p4().size(); iJet++) {
                        //                         double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jets_p4().at(iJet) );
                        //                         if( dr > deltaRCut && jets_p4().at(iJet).pt() > 10 ) nj1_++;
                        //                         if ( dr > deltaRCut && jets_p4().at(iJet).pt() > ptj1_ ){
                        //                             ptj1_ = jets_p4().at(iJet).pt();

                        //                             // back to back in phi
                        //                             float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iLep), jets_p4().at(iJet) ) );
                        //                             if( dphi > deltaPhiCut && jets_p4().at(iJet).pt() > ptj1_b2b_ ){        
                        //                                 ptj1_b2b_   = jets_p4().at(iJet).pt();
                        //                                 dphij1_b2b_ = dphi;
                        //                             }
                        //                         }
                        //                     }
                        // #endif

                        // PF Jets
                        // Find the highest Pt pfjet separated by at least dRcut from this lepton and fill the pfjet Pt
                        ptpfj1_       = -999.0;
                        ptpfj1_b2b_   = -999.0;
                        dphipfj1_b2b_ = -999.0;
                        npfj1_        = 0;
                        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                            double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), pfjets_p4().at(iJet) );
                            if( dr > deltaRCut && pfjets_p4().at(iJet).pt() > 10 ) npfj1_++;
                            if ( dr > deltaRCut && pfjets_p4().at(iJet).pt() > ptpfj1_ ){
                                ptpfj1_ = pfjets_p4().at(iJet).pt();

                                // back to back in phi
                                float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iLep), pfjets_p4().at(iJet) ) );
                                if( dphi > deltaPhiCut && pfjets_p4().at(iJet).pt() > ptpfj1_b2b_ ){        
                                    ptpfj1_b2b_   = pfjets_p4().at(iJet).pt();
                                    dphipfj1_b2b_ = dphi;
                                }
                            }
                        }

                        // L2L3 PF Jets
                        // Find the highest Pt PF corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                        ptpfcj1_       = -999.0;
                        ptpfcj1_b2b_   = -999.0;
                        dphipfcj1_b2b_ = -999.0;
                        npfcj1_        = 0;
                        btagpfc_       = false;
                        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                            // JetID
                            if ( !passesPFJetID(iJet)) continue;
                            LorentzVector jp4 = pfjets_p4().at(iJet);
                            float jet_cor = jetCorrection(jp4, jet_pf_L2L3corrector);
                            //float jet_cor = pfjets_corL2L3().at(iJet);
                            LorentzVector jp4cor = jp4 * jet_cor;
                            if (jp4cor.pt() > 15 && pfjets_combinedSecondaryVertexBJetTag().at(iJet) > 0.679 ) btagpfc_ = true;
                            double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jp4cor );
                            if( dr > deltaRCut && jp4cor.pt() > 10 ) npfcj1_++;
                            if ( dr > deltaRCut && jp4cor.pt() > ptpfcj1_ ){
                                ptpfcj1_ = jp4cor.pt();

                                // back to back in phi
                                float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iLep), jp4cor ) );
                                if( dphi > deltaPhiCut && jp4cor.pt() > ptpfcj1_b2b_ ){
                                    ptpfcj1_b2b_   = jp4cor.pt();
                                    dphipfcj1_b2b_ = dphi;
                                }
                            }
                        }

                        // L1FastL2L3 PF Jets
                        // Find the highest Pt PF corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                        emfpfcL1Fj1_      = -999.0;
                        ptpfcL1Fj1_       = -999.0;
                        dphipfcL1Fj1_     = -999.0;
                        ptpfcL1Fj1_b2b_   = -999.0;
                        dphipfcL1Fj1_b2b_ = -999.0;
                        npfcL1Fj1_        = 0;
                        npfc30L1Fj1_      = 0;
                        npfc40L1Fj1_      = 0;
                        nbpfc40L1Fj1_     = 0;
                        npfc50L1Fj1_eth_  = 0;
                        npfc65L1Fj1_eth_  = 0;
                        btagpfcL1F_       = false;
                        rho_ = cms2.evt_rho();
                        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                            // JetID
                            if ( !passesPFJetID(iJet)) continue;
                            LorentzVector jp4 = pfjets_p4().at(iJet);
                            float jet_cor = cms2.pfjets_corL1FastL2L3().at(iJet);
                            //jet_pf_L1FastJetL2L3_corrector->setRho(cms2.evt_ww_rho_vor());
                            //jet_pf_L1FastJetL2L3_corrector->setJetA(cms2.pfjets_area().at(iJet));
                            //jet_pf_L1FastJetL2L3_corrector->setJetPt(cms2.pfjets_p4().at(iJet).pt());
                            //jet_pf_L1FastJetL2L3_corrector->setJetEta(cms2.pfjets_p4().at(iJet).eta()); 
                            //float jet_cor = jetCorrection(jp4, jet_pf_L1FastJetL2L3_corrector);
                            LorentzVector jp4cor = jp4 * jet_cor;
                            if (jp4cor.pt() > 15 && pfjets_combinedSecondaryVertexBJetTag().at(iJet) > 0.679 ) btagpfcL1F_ = true;
                            double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jp4cor );
                            if( dr > deltaRCut && jp4cor.pt() > 10 ) npfcL1Fj1_++;
                            if( dr > deltaRCut && jp4cor.pt() > 30 ) npfc30L1Fj1_++;
                            if( dr > deltaRCut && jp4cor.pt() > 40 ) npfc40L1Fj1_++;
                            if( dr > deltaRCut && jp4cor.pt() > 40 && pfjets_combinedSecondaryVertexBJetTag().at(iJet) > 0.679) nbpfc40L1Fj1_++;
                            if( dr > 0.4       && jp4cor.pt() > 50 ) npfc50L1Fj1_eth_++;
                            if( dr > 0.4       && jp4cor.pt() > 65 ) npfc65L1Fj1_eth_++;
                            if ( dr > deltaRCut && jp4cor.pt() > ptpfcL1Fj1_ ){
                                emfpfcL1Fj1_ = (cms2.pfjets_chargedEmE().at(iJet) + cms2.pfjets_neutralEmE().at(iJet)) / pfjets_p4().at(iJet).E();
                                ptpfcL1Fj1_ = jp4cor.pt();
                                float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iLep), jp4cor ) );
                                dphipfcL1Fj1_ = dphi;

                                // back to back in phi
                                if( dphi > deltaPhiCut && jp4cor.pt() > ptpfcL1Fj1_b2b_ ){
                                    ptpfcL1Fj1_b2b_   = jp4cor.pt();
                                    dphipfcL1Fj1_b2b_ = dphi;
                                }
                            }
                        }

                        // L1FastL2L3Residual PF Jets
                        // Find the highest Pt PF corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                        emfpfcL1Fj1res_      = -999.0;
                        ptpfcL1Fj1res_       = -999.0;
                        dphipfcL1Fj1res_      = -999.0;
                        ptpfcL1Fj1res_b2b_   = -999.0;
                        dphipfcL1Fj1res_b2b_ = -999.0;
                        npfcL1Fj1res_        = 0;
                        npfc30L1Fj1res_      = 0;
                        npfc40L1Fj1res_      = 0;
                        nbpfc40L1Fj1res_     = 0;
                        npfc50L1Fj1res_eth_  = 0;
                        npfc65L1Fj1res_eth_  = 0;
                        btagpfcL1Fres_       = false;
                        rho_ = cms2.evt_rho();
                        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                            // JetID
                            if ( !passesPFJetID(iJet)) continue;
                            LorentzVector jp4 = pfjets_p4().at(iJet);
                            float jet_cor = cms2.pfjets_corL1FastL2L3residual().at(iJet);
                            //jet_pf_L1FastJetL2L3Residual_corrector->setRho(cms2.evt_ww_rho_vor());
                            //jet_pf_L1FastJetL2L3Residual_corrector->setJetA(cms2.pfjets_area().at(iJet));
                            //jet_pf_L1FastJetL2L3Residual_corrector->setJetPt(cms2.pfjets_p4().at(iJet).pt());
                            //jet_pf_L1FastJetL2L3Residual_corrector->setJetEta(cms2.pfjets_p4().at(iJet).eta()); 
                            //float jet_cor = jetCorrection(jp4, jet_pf_L1FastJetL2L3Residual_corrector);
                            LorentzVector jp4cor = jp4 * jet_cor;
                            if (jp4cor.pt() > 15 && pfjets_combinedSecondaryVertexBJetTag().at(iJet) > 0.679 ) btagpfcL1Fres_ = true;
                            double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jp4cor );
                            if( dr > deltaRCut && jp4cor.pt() > 10 ) npfcL1Fj1res_++;
                            if( dr > deltaRCut && jp4cor.pt() > 30 ) npfc30L1Fj1res_++;
                            if( dr > deltaRCut && jp4cor.pt() > 40 ) npfc40L1Fj1res_++;
                            if( dr > deltaRCut && jp4cor.pt() > 40 && pfjets_combinedSecondaryVertexBJetTag().at(iJet) > 0.679) nbpfc40L1Fj1res_++;
                            if( dr > 0.4       && jp4cor.pt() > 50 ) npfc50L1Fj1res_eth_++;
                            if( dr > 0.4       && jp4cor.pt() > 65 ) npfc65L1Fj1res_eth_++;
                            if ( dr > deltaRCut && jp4cor.pt() > ptpfcL1Fj1res_ ){
                                emfpfcL1Fj1res_ = (cms2.pfjets_chargedEmE().at(iJet) + cms2.pfjets_neutralEmE().at(iJet)) / pfjets_p4().at(iJet).E();
                                ptpfcL1Fj1res_ = jp4cor.pt();
                                float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iLep), jp4cor ) );
                                dphipfcL1Fj1res_ = dphi;

                                // back to back in phi
                                if( dphi > deltaPhiCut && jp4cor.pt() > ptpfcL1Fj1res_b2b_ ){
                                    ptpfcL1Fj1res_b2b_   = jp4cor.pt();
                                    dphipfcL1Fj1res_b2b_ = dphi;
                                }
                            }
                        }

                        //***  Doing B-tagging correctly ***
                        // B-tagged L1FastL2L3 PF Jets
                        // Find the highest Pt B-tagged PF L1FastL2L3 corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                        ptbtagpfcL1Fj1_       = -999.0;
                        dphibtagpfcL1Fj1_       = -999.0;
                        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                            if ( !passesPFJetID(iJet)) continue;
                            LorentzVector jp4 = pfjets_p4().at(iJet);
                            float jet_cor = cms2.pfjets_corL1FastL2L3().at(iJet);
                            //jet_pf_L1FastJetL2L3_corrector->setRho(cms2.evt_ww_rho_vor());
                            //jet_pf_L1FastJetL2L3_corrector->setJetA(cms2.pfjets_area().at(iJet));
                            //jet_pf_L1FastJetL2L3_corrector->setJetPt(cms2.pfjets_p4().at(iJet).pt());
                            //jet_pf_L1FastJetL2L3_corrector->setJetEta(cms2.pfjets_p4().at(iJet).eta()); 
                            //float jet_cor = jetCorrection(jp4, jet_pf_L1FastJetL2L3_corrector);
                            LorentzVector jp4cor = jp4 * jet_cor;
                            if (pfjets_combinedSecondaryVertexBJetTag().at(iJet) < 0.679 ) continue;
                            double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jp4cor );
                            if ( dr > deltaRCut && jp4cor.pt() > ptbtagpfcL1Fj1_ ){
                                ptbtagpfcL1Fj1_ = jp4cor.pt();
                                float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iLep), jp4cor ) );
                                dphibtagpfcL1Fj1_ = dphi; 
                            }
                        }

                        //***  Doing B-tagging correctly ***
                        // B-tagged L1FastL2L3Residual PF Jets
                        // Find the highest Pt B-tagged PF L1FastL2L3Residual corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                        ptbtagpfcL1Fj1res_       = -999.0;
                        dphibtagpfcL1Fj1res_     = -999.0;
                        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                            if ( !passesPFJetID(iJet)) continue;
                            LorentzVector jp4 = pfjets_p4().at(iJet);
                            float jet_cor = cms2.pfjets_corL1FastL2L3residual().at(iJet);
                            //jet_pf_L1FastJetL2L3Residual_corrector->setRho(cms2.evt_ww_rho_vor());
                            //jet_pf_L1FastJetL2L3Residual_corrector->setJetA(cms2.pfjets_area().at(iJet));
                            //jet_pf_L1FastJetL2L3Residual_corrector->setJetPt(cms2.pfjets_p4().at(iJet).pt());
                            //jet_pf_L1FastJetL2L3Residual_corrector->setJetEta(cms2.pfjets_p4().at(iJet).eta()); 
                            //float jet_cor = jetCorrection(jp4, jet_pf_L1FastJetL2L3Residual_corrector);
                            LorentzVector jp4cor = jp4 * jet_cor;
                            if (pfjets_combinedSecondaryVertexBJetTag().at(iJet) < 0.679 ) continue;
                            double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jp4cor );
                            if ( dr > deltaRCut && jp4cor.pt() > ptbtagpfcL1Fj1res_ ){
                                ptbtagpfcL1Fj1res_ = jp4cor.pt();
                                float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iLep), jp4cor ) );
                                dphibtagpfcL1Fj1res_ = dphi; 
                            }
                        }


                        //////////////
                        // End Jets //
                        //////////////

                        ///////////////////
                        // B Tagging     //
                        ///////////////////

                        // #ifndef __CMS2_SLIM__
                        //                     // The btag information
                        //                     nbjet_ = this_nbjet;
                        //                     dRbNear_ =  99.;
                        //                     dRbFar_  = -99.;
                        //                     for (int ii=0; ii<nbjet_; ii++) {
                        //                         unsigned int iJet = bindex.at(ii);
                        //                         float dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jets_p4().at(iJet));
                        //                         if (dr < dRbNear_) dRbNear_ = dr;
                        //                         if (dr > dRbFar_)  dRbFar_  = dr;
                        //                     }
                        // #endif

                        // The btag information for pfjets
                        nbpfcjet_ = this_nbpfjet;
                        dRbpfcNear_ = 99.;
                        dRbpfcFar_  = -99.;
                        //FactorizedJetCorrector* jet_pf_corrector = evt_isRealData() ? jet_pf_L1FastJetL2L3Residual_corrector : jet_pf_L1FastJetL2L3_corrector; 
                        for (int ii=0; ii<nbpfcjet_; ii++) {
                            unsigned int iJet = bpfindex.at(ii);
                            LorentzVector jp4 = pfjets_p4().at(iJet);
                            //jet_pf_corrector->setRho(cms2.evt_ww_rho_vor());
                            //jet_pf_corrector->setJetA(cms2.pfjets_area().at(iJet));
                            //jet_pf_corrector->setJetPt(cms2.pfjets_p4().at(iJet).pt());
                            //jet_pf_corrector->setJetEta(cms2.pfjets_p4().at(iJet).eta()); 
                            //float jet_cor = jetCorrection(jp4, jet_pf_corrector);
                            float jet_cor = evt_isRealData() ? cms2.pfjets_corL1FastL2L3residual().at(iJet) : cms2.pfjets_corL1FastL2L3().at(iJet);
                            LorentzVector jp4cor = jp4 * jet_cor;
                            float dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jp4cor);
                            if (dr < dRbpfcNear_) dRbpfcNear_ = dr;
                            if (dr > dRbpfcFar_)   dRbpfcFar_  = dr;
                        }

                        ///////////////////
                        // End B Tagging //
                        ///////////////////


                        // Time to fill the baby for the muons
                        FillBabyNtuple();

                    }// closes loop over muons
                } // closes if statements about whether we want to fill muons

            }// closes loop over events
            //printf("Good events found: %d out of %d\n",nGoodEvents,nEntries);

            f->Close();
            //delete f;

        }  // closes loop over files

        std::cout << "nEventTotal = " << nEventsTotal << endl;
        std::cout << "nEventChain = " << nEventsChain << endl;

        bmark.Stop("benchmark");
        cout << endl;
        cout << nEventsTotal << " Events Processed" << endl;
        //cout << "# of bad events filtered = " << bad_events << endl; 
        //cout << "# of duplicates filtered = " << duplicates << endl; 
        cout << "------------------------------" << endl;
        cout << "CPU  Time: " << Form("%.01f", bmark.GetCpuTime("benchmark" )) << endl;
        cout << "Real Time: " << Form("%.01f", bmark.GetRealTime("benchmark")) << endl;
        cout << endl;

        CloseBabyNtuple();
        return;

    }
    catch (std::exception& e)
    {
        cout << e.what() << endl;
    }

} // closes myLooper function  

