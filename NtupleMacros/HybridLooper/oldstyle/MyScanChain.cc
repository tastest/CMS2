
//
// Yanyan "glorious" Gao
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
#include "../../CORE/trackSelections.h"
#include "../../CORE/electronSelections.h"
#include "../../Tools/DileptonHypType.h"
#include "../../CORE/eventSelections.h"
#include "../../CORE/muonSelections.h"
#include "../../CORE/metSelections.h"

//
// Namespaces
//
using namespace tas;

//
//
//

enum mu_selection {
    PASS_MU_PT,
    PASS_MU_NOSECOND,
    PASS_MU_ISFIDUCIAL,
    PASS_MU_MET,
    PASS_MU_ISO,
};

enum ele_selection {
    PASS_ELE_PT,
    PASS_ELE_NOSECOND,
    PASS_ELE_ISFIDUCIAL,
    PASS_ELE_MET,
    PASS_ELE_ISO,
    PASS_ELE_JETVETO,
    PASS_ELE_R19,
    PASS_ELE_DPHI,
    PASS_ELE_DETA,
    PASS_ELE_HOE,
    PASS_ELE_LSHAPE,
    PASS_ELE_D0,
    PASS_ELE_DETA_CAND02,
    PASS_ELE_LSHAPE_CAND02,
    PASS_ELE_EXTRA,
    PASS_ELE_CONV,
    PASS_ELE_NOMUON,
    PASS_ELE_CLEANEVENT,
    PASS_ELE_ANTIMET,
};

double dRbetweenVectors(const LorentzVector &vec1, const LorentzVector &vec2 ){

    double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
    double deta = vec1.Eta() - vec2.Eta();
    return sqrt(dphi*dphi + deta*deta);
} 

// to decdie if to fill the EB histo (zero in the array)
// or the EE histo (one in the array)
enum DetType MyScanChain::getSubdet(int eleIndex)
{
    // determine what detector the electron is in
    if (cms2.els_fiduciality()[eleIndex] & (1<<ISEB)) return DET_EB;
    else if (cms2.els_fiduciality()[eleIndex] & (1<<ISEE)) return DET_EE;
    std::cout << "ERROR! not in EB or EE" << std::endl;
    return DET_ALL;
}

void MyScanChain::FillHist(TH1F** hist, DetType det, const float value, const float weight) 
{
    hist[det]->Fill(value, weight);
    hist[DET_ALL]->Fill(value, weight);
}

void MyScanChain::Format2DHist(TH2F** hist, std::string name, Int_t nx, Float_t minx, Float_t maxx,
        Int_t ny, Float_t miny, Float_t maxy)
{
    // loop on EB, EE
    for (unsigned int i = 0; i < 3; ++i)
    {
        std::string det = det_names[i];
        hist[i] = new TH2F(Form("%s_%s_%s", sampleName_.c_str(), name.c_str(), det.c_str()),
                name.c_str(), nx, minx, maxx, ny, miny, maxy);
    }
}

void MyScanChain::FormatHist(TH1F** hist, std::string name, Int_t n, Float_t min, Float_t max)
{
    // loop on EB, EE
    for (unsigned int i = 0; i < 3; ++i)
    {
        std::string det = det_names[i];
        hist[i] = new TH1F(Form("%s_%s_%s", sampleName_.c_str(), name.c_str(), det.c_str()),
                name.c_str(), n, min, max);
        hist[i]->GetXaxis()->SetTitle(name.c_str());
    }
}

void MyScanChain::FormatEffHist(EffMulti** hist, bool lessThan, float thresholdEB, float thresholdEE, std::string name)
{
    // loop on EB, EE
    for (unsigned int i = 0; i < 3; ++i)
    {
        std::string det = det_names[i];
        hist[i] = new EffMulti(lessThan, thresholdEB, thresholdEE, sampleName_, name, det);
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

//
// Analyse Electrons
//

void MyScanChain::AnalyseElectrons(const float &weight) {

    std::cout << "doing electrons" << std::endl;

    // find candidate electron
    // assumes electrons are sorted by pT descending
    int eleIndex = 0;
    int eleSecondIndex = 0;
    bool foundFirst = false;
    bool foundSecond = false;
    for (size_t i = 0; i < cms2.evt_nels(); ++i)
    {
        // must be found by the ecal driven algorithms
        if (! cms2.els_type()[i] & (1<<ISECALDRIVEN)) continue;

        if (foundFirst && !foundSecond) {
            eleSecondIndex = i;
            foundSecond = true;
            break;
        }
        if (!foundFirst) {
            eleIndex = i;
            foundFirst = true;
        }
    }

    // must have found first electron
    if (!foundFirst) return;

    // get subdetector (for histogramming)
    DetType det = getSubdet(eleIndex);

    //
    // work out what cuts this event passes
    //

    // the cuts that this event passes
    cuts_t cuts_passed = 0;

    // pt cut
    if (cms2.els_p4()[eleIndex].Pt() > 20.0) cuts_passed |= (1<<PASS_ELE_PT);

    // not a spike
    float r19 = cms2.els_eMax()[eleIndex]/cms2.els_e5x5()[eleIndex];
    if (r19 < 0.95) cuts_passed |= (1<<PASS_ELE_R19);

    // don't allow events with a second electron above 20.0 GeV
    float secondPt = 0.0;
    if (foundSecond) secondPt = cms2.els_p4()[eleSecondIndex].Pt();
    if (secondPt < 20.0) cuts_passed |= (1<<PASS_ELE_NOSECOND);

    // impose fiducial cuts in Eta
    if (fabs(cms2.els_etaSC()[eleIndex]) < 1.4442
            || (fabs(cms2.els_etaSC()[eleIndex]) > 1.560 && fabs(cms2.els_etaSC()[eleIndex]) < 2.500)) cuts_passed |= (1<<PASS_ELE_ISFIDUCIAL);

    // isolation
    float iso_relsusy = electronIsolation_relsusy_cand1(eleIndex, true);
    if (iso_relsusy < 0.10) cuts_passed |= (1<<PASS_ELE_ISO);

    // met
    //if (cms2.evt_tcmet() > 20.0 ) cuts_passed |= (1<<PASS_ELE_MET);
    if (cms2.evt_tcmet() > 20.0 && cms2.evt_pfmet() > 20) cuts_passed |= (1<<PASS_ELE_MET);
    if (cms2.evt_tcmet() < 15.0 && cms2.evt_pfmet() < 15.0) cuts_passed |= (1<<PASS_ELE_ANTIMET);   

    // ratio of the met to the muon pt
    float tcmetratio = cms2.evt_tcmet() / cms2.els_p4()[eleIndex].Pt();
    float pfmetratio = cms2.evt_pfmet() / cms2.els_p4()[eleIndex].Pt();

    // phi angle between the met and the muon
    float tcmetdphi = acos(cos(cms2.evt_tcmetPhi() - cms2.els_p4()[eleIndex].Phi()));
    float pfmetdphi = acos(cos(cms2.evt_pfmetPhi() - cms2.els_p4()[eleIndex].Phi()));

    // jet veto
    // leading pT JPT jet that is dR > 0.4 from the nearest electron
    float leadingJPT = 0.0;
    int leadingJPTIndex = 0;
    for (size_t j = 0; j < cms2.jpts_p4().size(); ++j)
    {
        if ( TMath::Abs(cms2.jpts_p4()[j].eta()) > 2.5 ) continue;
        if ( TMath::Abs(dRbetweenVectors(cms2.els_p4()[eleIndex], cms2.jpts_p4()[j])) < 0.4) continue;

        // find leading pT JPT
        if (cms2.jpts_p4()[j].Et() > leadingJPT) {
            leadingJPT = cms2.jpts_p4()[j].Et();
            leadingJPTIndex = j;
        }
    }

    if (leadingJPT < 30.0) cuts_passed |= (1<<PASS_ELE_JETVETO);

    //
    // do plotting
    //

    // compute some common variables 
    float e2x5MaxOver5x5 = cms2.els_e2x5Max()[eleIndex]/cms2.els_e5x5()[eleIndex];
    float ppfmet = projectedMETW(cms2.evt_pfmet(), cms2.evt_pfmetPhi(), cms2.els_p4()[eleIndex].Phi());
    float ptcmet = projectedMETW(cms2.evt_tcmet(), cms2.evt_tcmetPhi(), cms2.els_p4()[eleIndex].Phi());
    float pfmetsignificance = cms2.evt_pfmet() / sqrt(cms2.evt_pfsumet());
    float tcmetsignificance = cms2.evt_tcmet() / sqrt(cms2.evt_tcsumet());
    float pftransmass = sqrt( 2.0 * cms2.els_p4()[eleIndex].Pt() * cms2.evt_pfmet() 
            * (1 - cos(cms2.evt_pfmetPhi() - cms2.els_p4()[eleIndex].Phi() )));
    float tctransmass = sqrt( 2.0 * cms2.els_p4()[eleIndex].Pt() * cms2.evt_tcmet() 
            * (1 - cos(cms2.evt_tcmetPhi() - cms2.els_p4()[eleIndex].Phi() )));

    FillHist(h1_ele_pt_, det, cms2.els_p4()[eleIndex].Pt(), weight);
    FillHist(h1_ele_eta_, det, cms2.els_etaSC()[eleIndex], weight);
    FillHist(h1_ele_phi_, det, cms2.els_phiSC()[eleIndex], weight);
    FillHist(h1_ele_tcmet_, det, cms2.evt_tcmet(), weight);
    FillHist(h1_ele_pfmet_, det, cms2.evt_pfmet(), weight);
    FillHist(h1_ele_tcmetdphi_, det, tcmetdphi, weight);
    FillHist(h1_ele_pfmetdphi_, det, pfmetdphi, weight);
    FillHist(h1_ele_tcmetratio_, det, tcmetratio, weight);
    FillHist(h1_ele_pfmetratio_, det, pfmetratio, weight);

    const cuts_t pass_all =     (1<<PASS_ELE_PT) | (1<<PASS_ELE_NOSECOND) | 
        (1<<PASS_ELE_ISFIDUCIAL) | (1<<PASS_ELE_ISO) | 
        (1<<PASS_ELE_MET) | (1<<PASS_ELE_JETVETO) |
        (1<<PASS_ELE_R19);

    const cuts_t pass_all_antiselection = (1<<PASS_ELE_PT) | (1<<PASS_ELE_NOSECOND) |
        (1<<PASS_ELE_ISFIDUCIAL) | (1<<PASS_ELE_ISO) |
        (1<<PASS_ELE_ANTIMET) | (1<<PASS_ELE_JETVETO) |
        (1<<PASS_ELE_R19);

    if (CheckCutsNM1(pass_all, (1<<PASS_ELE_MET), cuts_passed)) {
        FillHist(h1_ele_nm1_tcmet_, det, cms2.evt_tcmet(), weight);
        FillHist(h1_ele_nm1_pfmet_, det, cms2.evt_pfmet(), weight);
        FillHist(h1_ele_nm1_tcmetdphi_, det, tcmetdphi, weight);
        FillHist(h1_ele_nm1_pfmetdphi_, det, pfmetdphi, weight);
        FillHist(h1_ele_nm1_tcmetratio_, det, tcmetratio, weight);
        FillHist(h1_ele_nm1_pfmetratio_, det, pfmetratio, weight);
    }

    if (CheckCutsNM1(pass_all, (1<<PASS_ELE_ISO), cuts_passed)) {
        FillHist(h1_ele_nm1_iso_, det, iso_relsusy, weight);
    }

    if (CheckCutsNM1(pass_all, (1<<PASS_ELE_JETVETO), cuts_passed)) {
        FillHist(h1_ele_nm1_jetveto_, det, leadingJPT, weight);
    }

    if (CheckCutsNM1(pass_all, (1<<PASS_ELE_NOSECOND), cuts_passed)) {
        FillHist(h1_ele_nm1_secondpt_, det, secondPt, weight);
    }

    if (CheckCutsNM1(pass_all, (1<<PASS_ELE_R19), cuts_passed)) {
        FillHist(h1_ele_nm1_r19_, det, r19, weight);
    }

    if (CheckCutsNM1(pass_all, ((1<<PASS_ELE_R19) | (1<<PASS_ELE_MET)), cuts_passed)) {
        FillHist(h1_ele_nm1nor19_tcmet_, det, cms2.evt_tcmet(), weight);
        FillHist(h1_ele_nm1nor19_pfmet_, det, cms2.evt_pfmet(), weight);
        FillHist(h1_ele_nm1nor19_tcmetratio_, det, tcmetratio, weight);
        FillHist(h1_ele_nm1nor19_pfmetratio_, det, pfmetratio, weight);
    }

    // apply all cuts
    if (CheckCuts(pass_all, cuts_passed)) {

        // print out details of event passing
        // the full selection
        if (isData_) {
            _asciifile    << "********************************************************************"
                << std::endl;
            _asciifile    << "ELECTRONS: \t"
                << cms2.evt_run() << "\t\t" 
                << cms2.evt_lumiBlock() << "\t" 
                << cms2.evt_event() << std::endl;
            _asciifile    << "--------------------------------------------------------------------"
                << std::endl;
            _asciifile   << "Pt = "<<cms2.els_p4()[eleIndex].Pt() 
                << "\t tcMet = "<< cms2.evt_tcmet()
                << "\t Projected tcMet = "<< ptcmet
                << "\t Transverse Mass (tcMet) = "<< tctransmass
                <<"\n"
                << "\t pfMet = "<< cms2.evt_pfmet()
                << "\t Transverse Mass (pfMet) = "<< pftransmass
                << "\t Projected pfMet = "<< ppfmet
                << std::endl;
        }

        FillHist(h1_ele_selected_pt_, det, cms2.els_p4()[eleIndex].Pt(), weight);
        FillHist(h1_ele_selected_eta_, det, cms2.els_etaSC()[eleIndex], weight);
        FillHist(h1_ele_selected_phi_, det, cms2.els_phiSC()[eleIndex], weight);
        FillHist(h1_ele_selected_tcmet_, det, cms2.evt_tcmet(), weight);
        FillHist(h1_ele_selected_pfmet_, det, cms2.evt_pfmet(), weight);
        FillHist(h1_ele_selected_tcmetdphi_, det, tcmetdphi, weight);
        FillHist(h1_ele_selected_pfmetdphi_, det, pfmetdphi, weight);
        FillHist(h1_ele_selected_tcmetratio_, det, tcmetratio, weight);
        FillHist(h1_ele_selected_pfmetratio_, det, pfmetratio, weight);        
        FillHist(h1_ele_selected_ptcmet_, det, ptcmet, weight);
        FillHist(h1_ele_selected_ppfmet_, det, ppfmet, weight);
        FillHist(h1_ele_selected_tcmetsignificance_, det, tcmetsignificance, weight);
        FillHist(h1_ele_selected_pfmetsignificance_, det, pfmetsignificance, weight);
        FillHist(h1_ele_selected_tctransmass_, det, tctransmass, weight);
        FillHist(h1_ele_selected_pftransmass_, det, pftransmass, weight);
        FillHist(h1_ele_selected_d0corr_, det, cms2.els_d0corr()[eleIndex], weight);
        FillHist(h1_ele_selected_nmhits_, det, cms2.els_exp_innerlayers()[eleIndex], weight);

        FillHist(h1_ele_selected_sigmaIEtaIEta_, det, cms2.els_sigmaIEtaIEta()[eleIndex], weight);
        FillHist(h1_ele_selected_dEtaIn_, det, cms2.els_dEtaIn()[eleIndex], weight);
        FillHist(h1_ele_selected_dPhiIn_, det, cms2.els_dPhiIn()[eleIndex], weight);
        FillHist(h1_ele_selected_hOverE_, det, cms2.els_hOverE()[eleIndex], weight);
        FillHist(h1_ele_selected_e2x5MaxOver5x5_, det, e2x5MaxOver5x5, weight);
        FillHist(h1_ele_selected_fbrem_, det, cms2.els_fbrem()[eleIndex], weight);
        FillHist(h1_ele_selected_eOverPIn_, det, cms2.els_eOverPIn()[eleIndex], weight);

    }

    // apply all cuts with antiselection on met instead of selection
    // to select background like electrons

    if (CheckCuts(pass_all_antiselection, cuts_passed)) {
        FillHist(h1_ele_antiselected_sigmaIEtaIEta_, det, cms2.els_sigmaIEtaIEta()[eleIndex], weight);
        FillHist(h1_ele_antiselected_dEtaIn_, det, cms2.els_dEtaIn()[eleIndex], weight);
        FillHist(h1_ele_antiselected_dPhiIn_, det, cms2.els_dPhiIn()[eleIndex], weight);
        FillHist(h1_ele_antiselected_hOverE_, det, cms2.els_hOverE()[eleIndex], weight);
        FillHist(h1_ele_antiselected_e2x5MaxOver5x5_, det, e2x5MaxOver5x5, weight);
        FillHist(h1_ele_antiselected_fbrem_, det, cms2.els_fbrem()[eleIndex], weight);
        FillHist(h1_ele_antiselected_eOverPIn_, det, cms2.els_eOverPIn()[eleIndex], weight);

    }

}

//
// Analyse Muons
//

void MyScanChain::AnalyseMuons(const float &weight) {

    std::cout << "doing muons" << std::endl;

    // find candidate muon
    // assumes muons are sorted by pT descending
    int muIndex = 0;
    int muSecondIndex = 0;
    bool foundFirst = false;
    bool foundSecond = false;
    for (size_t i = 0; i < cms2.mus_p4().size(); ++i)
    {
        // must be tracker and global
        // see CORE/muonSelections.cc
        if (((cms2.mus_type()[i]) & (1<<1)) == 0)    continue; // global muon
        if (((cms2.mus_type()[i]) & (1<<2)) == 0)    continue; // tracker muon

        if (foundFirst && !foundSecond) {
            muSecondIndex = i;
            foundSecond = true;
            break;
        }
        if (!foundFirst) {
            muIndex = i;
            foundFirst = true;
        }
    }

    // must have found first muon
    if (!foundFirst) return;

    // get subdetector (for histogramming)
    // not sure what the most sensible eta division
    // is for muons!
    DetType det = DET_EE;
    if (fabs(cms2.mus_p4()[muIndex].Eta()) < 1.479) det = DET_EB;

    //
    // work out what cuts this event passes
    //

    std::cout << "going to work out cuts" << std::endl;

    // the cuts that this event passes
    cuts_t cuts_passed = 0;

    // pt cut
    if (cms2.mus_p4()[muIndex].Pt() > 20.0) cuts_passed |= (1<<PASS_MU_PT);

    // don't allow events with a second muon above 20.0 GeV
    float secondPt = 0.0;
    if (foundSecond) secondPt = cms2.mus_p4()[muSecondIndex].Pt();
    if (secondPt < 20.0) cuts_passed |= (1<<PASS_MU_NOSECOND);

    // impose fiducial cuts in Eta
    if (fabs(cms2.mus_p4()[muIndex].Eta()) < 2.5) cuts_passed |= (1<<PASS_MU_ISFIDUCIAL);

    std::cout << "mark 1" << std::endl;

    // met
    if (cms2.evt_tcmet() > 20.0 && cms2.evt_pfmet() > 20) cuts_passed |= (1<<PASS_MU_MET);

    // ratio of the met to the muon pt
    float tcmetratio = cms2.evt_tcmet() / cms2.mus_p4()[muIndex].Pt();
    float pfmetratio = cms2.evt_pfmet() / cms2.mus_p4()[muIndex].Pt();

    std::cout << "mark 2" << std::endl;

    // phi angle between the met and the muon
    float tcmetdphi = acos(cos(cms2.evt_tcmetPhi() - cms2.mus_p4()[muIndex].Phi()));
    float pfmetdphi = acos(cos(cms2.evt_pfmetPhi() - cms2.mus_p4()[muIndex].Phi()));

    std::cout << "mark 3" << std::endl;

    // muon isolation value
    if (muonIsoValue(muIndex) < 0.1)   cuts_passed |= (1<<PASS_MU_ISO);
    //
    // do plotting
    //

    std::cout << "going to do muon plotting" << std::endl;

    // compute some common variables 
    float ppfmet = projectedMETW(cms2.evt_pfmet(), cms2.evt_pfmetPhi(), cms2.mus_p4()[muIndex].Phi());
    float ptcmet = projectedMETW(cms2.evt_tcmet(), cms2.evt_tcmetPhi(), cms2.mus_p4()[muIndex].Phi());
    float pfmetsignificance = cms2.evt_pfmet() / sqrt(cms2.evt_pfsumet());
    float tcmetsignificance = cms2.evt_tcmet() / sqrt(cms2.evt_tcsumet());
    float pftransmass = sqrt( 2.0 * cms2.mus_p4()[muIndex].Pt() * cms2.evt_pfmet()
            * (1 - cos(cms2.evt_pfmetPhi() - cms2.mus_p4()[muIndex].Phi() )));
    float tctransmass = sqrt( 2.0 * cms2.mus_p4()[muIndex].Pt() * cms2.evt_tcmet()
            * (1 - cos(cms2.evt_tcmetPhi() - cms2.mus_p4()[muIndex].Phi() )));

    FillHist(h1_mu_pt_, det, cms2.mus_p4()[muIndex].Pt(), weight);
    FillHist(h1_mu_eta_, det, cms2.mus_p4()[muIndex].Eta(), weight);
    FillHist(h1_mu_phi_, det, cms2.mus_p4()[muIndex].Phi(), weight);
    FillHist(h1_mu_tcmet_, det, cms2.evt_tcmet(), weight);
    FillHist(h1_mu_pfmet_, det, cms2.evt_pfmet(), weight);
    FillHist(h1_mu_tcmetdphi_, det, tcmetdphi, weight);
    FillHist(h1_mu_pfmetdphi_, det, pfmetdphi, weight);
    FillHist(h1_mu_tcmetratio_, det, tcmetratio, weight);
    FillHist(h1_mu_pfmetratio_, det, pfmetratio, weight);


    const cuts_t pass_all =     (1<<PASS_MU_PT) | (1<<PASS_MU_NOSECOND) |
        (1<<PASS_MU_ISO) | (1<<PASS_MU_ISFIDUCIAL) | (1<<PASS_MU_MET)
        ;

    if (CheckCutsNM1(pass_all, (1<<PASS_MU_NOSECOND), cuts_passed)) {
        FillHist(h1_mu_nm1_secondpt_, det, secondPt, weight);
    }

    if (CheckCutsNM1(pass_all, (1<<PASS_MU_MET), cuts_passed)) {
        FillHist(h1_mu_nm1_tcmet_, det, cms2.evt_tcmet(), weight);
        FillHist(h1_mu_nm1_pfmet_, det, cms2.evt_pfmet(), weight);
    }

    if (CheckCutsNM1(pass_all, (1<<PASS_MU_ISO), cuts_passed)) {
        FillHist(h1_mu_nm1_iso_, det, muonIsoValue(muIndex), weight);
    }

    // apply all cuts
    if (CheckCuts(pass_all, cuts_passed)) {

        // print out details of event passing
        // the full selection
        if (isData_) {
            _asciifile    << "********************************************************************"
                << std::endl;
            _asciifile       << "MUONS: \t\t"
                << cms2.evt_run() << "\t\t"
                << cms2.evt_lumiBlock() << "\t"
                << cms2.evt_event() << std::endl;
            _asciifile    << "--------------------------------------------------------------------"
                << std::endl;
            _asciifile   << "Pt = "<<cms2.mus_p4()[muIndex].Pt() 
                << "\t tcMet = "<< cms2.evt_tcmet()
                << "\t Projected tcMet = "<< ptcmet
                << "\t Transverse Mass (tcMet) = "<< tctransmass
                << "\n"
                << "\t pfMet = "<< cms2.evt_pfmet()
                << "\t Transverse Mass (pfMet) = "<< pftransmass
                << "\t Projected pfMet = "<< ppfmet
                << std::endl;
        }

        FillHist(h1_mu_selected_pt_, det, cms2.mus_p4()[muIndex].Pt(), weight);
        FillHist(h1_mu_selected_eta_, det, cms2.mus_p4()[muIndex].Eta(), weight);
        FillHist(h1_mu_selected_phi_, det, cms2.mus_p4()[muIndex].Phi(), weight);
        FillHist(h1_mu_selected_tcmet_, det, cms2.evt_tcmet(), weight);
        FillHist(h1_mu_selected_pfmet_, det, cms2.evt_pfmet(), weight);
        FillHist(h1_mu_selected_tcmetdphi_, det, tcmetdphi, weight);
        FillHist(h1_mu_selected_pfmetdphi_, det, pfmetdphi, weight);
        FillHist(h1_mu_selected_tcmetratio_, det, tcmetratio, weight);
        FillHist(h1_mu_selected_pfmetratio_, det, pfmetratio, weight);
        FillHist(h1_mu_selected_ptcmet_, det, ptcmet, weight);
        FillHist(h1_mu_selected_ppfmet_, det, ppfmet, weight);
        FillHist(h1_mu_selected_tcmetsignificance_, det, tcmetsignificance, weight);
        FillHist(h1_mu_selected_pfmetsignificance_, det, pfmetsignificance, weight);
        FillHist(h1_mu_selected_tctransmass_, det, tctransmass, weight);
        FillHist(h1_mu_selected_pftransmass_, det, pftransmass, weight);
        FillHist(h1_mu_selected_d0corr_, det, cms2.mus_d0corr()[muIndex], weight);

    }


}

//
// Main function
//
int MyScanChain::ScanChain(bool isData, std::string sampleName, TChain *chain, int nEvents, std::string skimFilePrefix) {

    //
    // define counters
    //

    // count the (weighted and unweighted) number of candidates passing our cuts
    float             cands_passing[3];
    float             cands_passing_w2[3];
    unsigned int       cands_count[3];
    memset(cands_passing   , 0, sizeof(cands_passing       ));
    memset(cands_passing_w2        , 0, sizeof(cands_passing_w2    ));
    memset(cands_count             , 0, sizeof(cands_count         ));

    // set sampleName
    sampleName_ = sampleName;
    isData_ = isData;

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

    //
    // Electrons
    //

    // General (before any cuts)
    FormatHist(h1_ele_pt_, "ele_pt", 100, 0, 100);
    FormatHist(h1_ele_eta_, "ele_eta", 100, -3, 3);
    FormatHist(h1_ele_phi_, "ele_phi", 100, -4, 4);
    FormatHist(h1_ele_tcmet_, "ele_tcmet", 100, 0, 100);
    FormatHist(h1_ele_pfmet_, "ele_pfmet", 100, 0, 100);
    FormatHist(h1_ele_tcmetdphi_, "ele_tcmetdphi", 100, -4, 4);
    FormatHist(h1_ele_pfmetdphi_, "ele_pfmetdphi", 100, -4, 4);
    FormatHist(h1_ele_tcmetratio_, "ele_tcmetratio", 50, 0, 5);
    FormatHist(h1_ele_pfmetratio_, "ele_pfmetratio", 50, 0, 5);

    // N-1
    FormatHist(h1_ele_nm1_tcmet_, "ele_nm1_tcmet", 100, 0, 100);
    FormatHist(h1_ele_nm1_pfmet_, "ele_nm1_pfmet", 100, 0, 100);
    FormatHist(h1_ele_nm1_tcmetdphi_, "ele_nm1_tcmetdphi", 100, -4, 4);
    FormatHist(h1_ele_nm1_pfmetdphi_, "ele_nm1_pfmetdphi", 100, -4, 4);
    FormatHist(h1_ele_nm1_tcmetratio_, "ele_nm1_tcmetratio", 50, 0, 5);
    FormatHist(h1_ele_nm1_pfmetratio_, "ele_nm1_pfmetratio", 50, 0, 5);
    FormatHist(h1_ele_nm1_jetveto_, "ele_nm1_jetveto", 100, 0, 100);
    FormatHist(h1_ele_nm1_iso_, "ele_nm1_iso", 100, 0, 1);
    FormatHist(h1_ele_nm1_secondpt_, "ele_nm1_secondpt", 100, 0, 100);

    FormatHist(h1_ele_nm1_r19_, "ele_nm1_r19", 120, 0, 1.2);
    FormatHist(h1_ele_nm1nor19_tcmet_, "ele_nm1nor19_tcmet", 100, 0, 100);
    FormatHist(h1_ele_nm1nor19_pfmet_, "ele_nm1nor19_pfmet", 100, 0, 100);
    FormatHist(h1_ele_nm1nor19_tcmetratio_, "ele_nm1nor19_tcmetratio", 50, 0, 5);
    FormatHist(h1_ele_nm1nor19_pfmetratio_, "ele_nm1nor19_pfmetratio", 50, 0, 5);

    // after all selections
    FormatHist(h1_ele_selected_pt_, "ele_selected_pt", 100, 0, 100);
    FormatHist(h1_ele_selected_eta_, "ele_selected_eta", 100, -3, 3);
    FormatHist(h1_ele_selected_phi_, "ele_selected_phi", 100, -4, 4);
    FormatHist(h1_ele_selected_tcmet_, "ele_selected_tcmet", 100, 0, 100);
    FormatHist(h1_ele_selected_pfmet_, "ele_selected_pfmet", 100, 0, 100);
    FormatHist(h1_ele_selected_tcmetdphi_, "ele_selected_tcmetdphi", 100, -4, 4);
    FormatHist(h1_ele_selected_pfmetdphi_, "ele_selected_pfmetdphi", 100, -4, 4);
    FormatHist(h1_ele_selected_tcmetratio_, "ele_selected_tcmetratio", 50, 0, 5);
    FormatHist(h1_ele_selected_pfmetratio_, "ele_selected_pfmetratio", 50, 0, 5);
    FormatHist(h1_ele_selected_tcmetsignificance_, "ele_selected_tcmetsignificance", 100, 0, 10.0);
    FormatHist(h1_ele_selected_pfmetsignificance_, "ele_selected_pfmetsignificance", 100, 0, 10.0);
    FormatHist(h1_ele_selected_ptcmet_, "ele_selected_ptcmet", 100, 0, 100);
    FormatHist(h1_ele_selected_ppfmet_, "ele_selected_ppfmet", 100, 0, 100);
    FormatHist(h1_ele_selected_tctransmass_, "ele_selected_tctransmass", 200, 0, 200);
    FormatHist(h1_ele_selected_pftransmass_, "ele_selected_pftransmass", 200, 0, 200);
    FormatHist(h1_ele_selected_d0corr_, "ele_selected_d0corr", 100, -0.2, 0.2);
    FormatHist(h1_ele_selected_nmhits_, "ele_selected_nmhits", 10, -0.5, 9.5);   

    FormatHist(h1_ele_selected_sigmaIEtaIEta_, "ele_selected_sigmaIEtaIEta", 100, 0, 0.1);
    FormatHist(h1_ele_selected_dEtaIn_, "ele_selected_dEtaIn", 100, -0.2, 0.2);
    FormatHist(h1_ele_selected_dPhiIn_, "ele_selected_dPhiIn", 100, -0.5, 0.5);
    FormatHist(h1_ele_selected_hOverE_, "ele_selected_hOverE", 100, 0, 1.0);
    FormatHist(h1_ele_selected_e2x5MaxOver5x5_, "ele_selected_e2x5MaxOver5x5", 110, 0, 1.1);
    FormatHist(h1_ele_selected_fbrem_, "ele_selected_fbrem", 100, 0, 1.0);
    FormatHist(h1_ele_selected_eOverPIn_, "ele_selected_eOverPIn", 100, 0, 5.0);


    // anti-selection on met
    FormatHist(h1_ele_antiselected_sigmaIEtaIEta_, "ele_antiselected_sigmaIEtaIEta", 100, 0, 0.1);    
    FormatHist(h1_ele_antiselected_dEtaIn_, "ele_antiselected_dEtaIn", 100, -0.2, 0.2);
    FormatHist(h1_ele_antiselected_dPhiIn_, "ele_antiselected_dPhiIn", 100, -0.5, 0.5);
    FormatHist(h1_ele_antiselected_hOverE_, "ele_antiselected_hOverE", 100, 0, 1.0);
    FormatHist(h1_ele_antiselected_e2x5MaxOver5x5_, "ele_antiselected_e2x5MaxOver5x5", 110, 0, 1.1);
    FormatHist(h1_ele_antiselected_fbrem_, "ele_antiselected_fbrem", 100, 0, 1.0);
    FormatHist(h1_ele_antiselected_eOverPIn_, "ele_antiselected_eOverPIn", 100, 0, 5.0);


    //
    // Muons
    //

    // General (before any cuts)
    FormatHist(h1_mu_pt_, "mu_pt", 100, 0, 100);
    FormatHist(h1_mu_eta_, "mu_eta", 100, -3, 3);
    FormatHist(h1_mu_phi_, "mu_phi", 100, -4, 4);
    FormatHist(h1_mu_tcmet_, "mu_tcmet", 100, 0, 100);
    FormatHist(h1_mu_pfmet_, "mu_pfmet", 100, 0, 100);
    FormatHist(h1_mu_tcmetdphi_, "mu_tcmetdphi", 100, -4, 4);
    FormatHist(h1_mu_pfmetdphi_, "mu_pfmetdphi", 100, -4, 4);
    FormatHist(h1_mu_tcmetratio_, "mu_tcmetratio", 50, 0, 5);
    FormatHist(h1_mu_pfmetratio_, "mu_pfmetratio", 50, 0, 5);

    // N-1
    FormatHist(h1_mu_nm1_secondpt_, "mu_nm1_secondpt", 100, 0, 100);
    FormatHist(h1_mu_nm1_tcmet_, "mu_nm1_tcmet", 100, 0, 100);
    FormatHist(h1_mu_nm1_pfmet_, "mu_nm1_pfmet", 100, 0, 100);
    FormatHist(h1_mu_nm1_iso_, "mu_nm1_iso", 100, 0, 1);

    // after all selections
    FormatHist(h1_mu_selected_pt_, "mu_selected_pt", 100, 0, 100);
    FormatHist(h1_mu_selected_eta_, "mu_selected_eta", 100, -3, 3);
    FormatHist(h1_mu_selected_phi_, "mu_selected_phi", 100, -4, 4);
    FormatHist(h1_mu_selected_tcmet_, "mu_selected_tcmet", 100, 0, 100);
    FormatHist(h1_mu_selected_pfmet_, "mu_selected_pfmet", 100, 0, 100);
    FormatHist(h1_mu_selected_tcmetdphi_, "mu_selected_tcmetdphi", 100, -4, 4);
    FormatHist(h1_mu_selected_pfmetdphi_, "mu_selected_pfmetdphi", 100, -4, 4);
    FormatHist(h1_mu_selected_tcmetratio_, "mu_selected_tcmetratio", 50, 0, 5);
    FormatHist(h1_mu_selected_pfmetratio_, "mu_selected_pfmetratio", 50, 0, 5);
    FormatHist(h1_mu_selected_tcmetsignificance_, "mu_selected_tcmetsignificance", 100, 0, 10.0);
    FormatHist(h1_mu_selected_pfmetsignificance_, "mu_selected_pfmetsignificance", 100, 0, 10.0);
    FormatHist(h1_mu_selected_ptcmet_, "mu_selected_ptcmet", 100, 0, 100);
    FormatHist(h1_mu_selected_ppfmet_, "mu_selected_ppfmet", 100, 0, 100);
    FormatHist(h1_mu_selected_tctransmass_, "mu_selected_tctransmass", 200, 0, 200);
    FormatHist(h1_mu_selected_pftransmass_, "mu_selected_pftransmass", 200, 0, 200);
    FormatHist(h1_mu_selected_d0corr_, "mu_selected_d0corr", 100, -0.2, 0.2);

    // open an asciifile to store results
    if(isData_) {
        _asciifile.open("whunt_output.txt"); 
        _asciifile << "Type \t\t"<<  "Run # \t\t" << "LumiBlock # \t" << "Event #\t\t" <<std::endl;
    }

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


            std::cout << "Event number: " << event << std::endl;

            // work out event weight
            float weight = 1.0;
            if (!isData) weight = cms2.evt_scale1fb();

            //
            // do event cleaning
            //
            if (!cleaning_standard(isData)) continue;

            //
            // Do the analysis
            //
            AnalyseElectrons(weight);
            std::cout << "done electrons" << std::endl;
            AnalyseMuons(weight);
            std::cout << "done muons" << std::endl;


            //
            // Count... Hmm think this bit needs fixing as it doesn't
            // mean much right now
            //
            /*
               cands_passing[det] += weight;
               cands_passing_w2[det] += weight * weight;
               cands_count[det] ++;

               cands_passing[DET_ALL] += weight;
               cands_passing_w2[DET_ALL] += weight * weight;
               cands_count[DET_ALL] ++;
             */

        } // end loop on files

    } // end loop on events

    if ( nEventsChain != nEventsTotal ) {
        std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    }

    //
    // print table entry
    //

    std::cout.flush();
    std::cout << std::endl;
    for (unsigned int i = 0; i < 3; ++i) {
        std::string str = det_names[i];
        std::cout << " & " << det_names[i] << "\t";
    }
    std::cout << "\\\\ \\hline" << std::endl;
    std::cout << sampleName << "\t";
    for (unsigned int i = 0; i < 3; ++i) {
        std::cout << " & " << cands_passing[i] << " $\\pm$ " << sqrt(cands_passing_w2[i]) << "\t";
    }
    std::cout << "\\\\ \\hline" << std::endl;

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

