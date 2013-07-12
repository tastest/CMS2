#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <set>
#include <algorithm>
// ROOT includes
#include "TSystem.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TChainElement.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TTreeCache.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"
// TAS includes
#include "../CORE/CMS2.h"
#include "../CORE/SimpleFakeRate.h"
#include "../CORE/trackSelections.h"
#include "../CORE/eventSelections.h"

#include "../CORE/muonSelections.h"
#include "../CORE/electronSelections.h"
#include "../CORE/electronSelectionsParameters.h"
#include "../CORE/metSelections.h"
#include "../CORE/mcSelections.h"
#include "../CORE/jetSelections.h"
#include "../CORE/ttbarSelections.h"
#include "../CORE/susySelections.h"
#include "../CORE/mcSUSYkfactor.h"
#include "../CORE/triggerSuperModel.h"
#include "../CORE/triggerUtils.h"
#include "../Tools/goodrun.cc"
#include "../CORE/utilities.cc"
#include "./topAFB_looper.h"
#include "Histograms.cc"
#include "../Tools/vtxreweight.cc"
#include "../Tools/bTagEff_BTV.cc"
#include "../CORE/jetcorr/SimpleJetCorrectionUncertainty.icc"
#include "../CORE/jetcorr/JetCorrectionUncertainty.icc"

// // top mass
// #include "../CORE/topmass/ttdilepsolve.cpp"
// #include "../CORE/topmass/getTopMassEstimate.icc"

using namespace std;
using namespace tas;

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;


double topAFB_looper::TopPtWeight(double topPt)
{
    if ( topPt < 0 ) return 1;
    if (topPt > 400) topPt = 400;

    //double result = (1.4 / 1000000.0) * topPt * topPt - (2.0 / 1000.0) * topPt + 1.2; //old 7TeV fit (l+j and 2l combined)
    double result = exp(0.199 - 0.00166 * topPt); //new 7TeV fit (l+j and 2l combined)
    //note this fit is for data/madgraph, and we are using MC@NLO



    return result;
}

double topAFB_looper::LeptonPtWeight(double leptonPt)
{

  if ( leptonPt < 20 ) return 1;
  if ( (leptonPt >= 20) && (leptonPt < 40) ) return 1.03545;
  if ( (leptonPt >= 40) && (leptonPt < 60) ) return 1.01794;
  if ( (leptonPt >= 60) && (leptonPt < 80) ) return 1.00371;
  if ( (leptonPt >= 80) && (leptonPt < 100) ) return 0.932273;
  if ( (leptonPt >= 100) && (leptonPt < 120) ) return 0.881938;
  if ( (leptonPt >= 120) && (leptonPt < 140) ) return 0.867833;
  if ( (leptonPt >= 140) && (leptonPt < 160) ) return 0.864478;
  if ( (leptonPt >= 160) && (leptonPt < 180) ) return 0.811505;
  if ( (leptonPt >= 180) && (leptonPt < 200) ) return 0.7997;
  if ( (leptonPt >= 200) && (leptonPt < 220) ) return 0.876389;
  if ( (leptonPt >= 220) && (leptonPt < 240) ) return 0.737906;  
  if ( leptonPt >= 240 ) return 1;
  
  return 1;
}

double topAFB_looper::JetPtWeight(double jetPt)
{

  if ( jetPt < 30 ) return 1;
  if ( (jetPt >= 30) && (jetPt < 60) ) return 1.04352;
  if ( (jetPt >= 60) && (jetPt < 90) ) return 1.02321;
  if ( (jetPt >= 90) && (jetPt < 120) ) return 0.952548;
  if ( (jetPt >= 120) && (jetPt < 150) ) return 0.89809;
  if ( (jetPt >= 150) && (jetPt < 180) ) return 0.862834;
  if ( (jetPt >= 180) && (jetPt < 210) ) return 0.857449;
  if ( (jetPt >= 210) && (jetPt < 240) ) return 0.896524;
  if ( (jetPt >= 240) && (jetPt < 270) ) return 0.871748;
  if ( (jetPt >= 270) && (jetPt < 300) ) return 0.747337;
  if ( (jetPt >= 300) && (jetPt < 330) ) return 0.74634;
  if ( (jetPt >= 330) && (jetPt < 360) ) return 0.692822;
  if ( jetPt >= 360 ) return 1;
  
  return 1;
}


bool topAFB_looper::passbTagging(const unsigned int jet_idx, const string jetAlgo, const string bTagDiscriminator)
{
    if (jetAlgo == "jptJets")
    {
        if (bTagDiscriminator == "trackCountingHighEffBJetTag" && jpts_trackCountingHighEffBJetTag()[jet_idx] > 3.3)
            return true;
        else if (bTagDiscriminator == "combinedSecondaryVertexBJetTag" && jpts_combinedSecondaryVertexBJetTag()[jet_idx] > 0.679)
            return true;
        else if (bTagDiscriminator == "simpleSecondaryVertexHighEffBJetTag" && jpts_simpleSecondaryVertexHighEffBJetTag()[jet_idx] > 1.74)
            return true;
        else if (bTagDiscriminator == "simpleSecondaryVertexHighPurBJetTag" && jpts_simpleSecondaryVertexHighPurBJetTags()[jet_idx] > 2)
            return true;
    }


    if (jetAlgo == "caloJets")
    {
        if (bTagDiscriminator == "trackCountingHighEffBJetTag" && jets_trackCountingHighEffBJetTag()[jet_idx] > 3.3)
            return true;
        else if (bTagDiscriminator == "combinedSecondaryVertexBJetTag" && jets_combinedSecondaryVertexBJetTag()[jet_idx] > 0.679)
            return true;
        else if (bTagDiscriminator == "simpleSecondaryVertexHighEffBJetTag" && jets_simpleSecondaryVertexHighEffBJetTag()[jet_idx] > 1.74)
            return true;
        else if (bTagDiscriminator == "simpleSecondaryVertexHighPurBJetTag" && jets_simpleSecondaryVertexHighPurBJetTags()[jet_idx] > 2)
            return true;
    }


    if (jetAlgo == "pfJets")
    {
        if (bTagDiscriminator == "trackCountingHighEffBJetTag" && pfjets_trackCountingHighEffBJetTag()[jet_idx] > 3.3)
            return true;
        else if (bTagDiscriminator == "combinedSecondaryVertexBJetTag" && pfjets_combinedSecondaryVertexBJetTag()[jet_idx] > 0.679)
            return true;
        else if (bTagDiscriminator == "simpleSecondaryVertexHighEffBJetTag" && pfjets_simpleSecondaryVertexHighEffBJetTag()[jet_idx] > 1.74)
            return true;
        else if (bTagDiscriminator == "simpleSecondaryVertexHighPurBJetTag" && pfjets_simpleSecondaryVertexHighPurBJetTags()[jet_idx] > 2)
            return true;
    }

    return false;

}
void topAFB_looper::FillHistograms(const unsigned int hypIdx, const vector<unsigned int> v_jets, const vector<unsigned int> v_jetsNoEtaCut,
                                   const pair<float, float> p_met, const float weight, const string prefix)
{

}

/* THIS NEEDS TO BE IN CORE */

struct DorkyEventIdentifier
{
    // this is a workaround for not having unique event id's in MC
    unsigned long int run, event, lumi;
    bool operator < (const DorkyEventIdentifier &) const;
    bool operator == (const DorkyEventIdentifier &) const;
};

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
    if (run != other.run)
        return run < other.run;
    if (event != other.event)
        return event < other.event;
    if (lumi != other.lumi)
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
    return sqrt( 2 * met * ( p4.pt() - ( p4.Px() * cos(met_phi) + p4.Py() * sin(met_phi) ) ) );
}
/*
void setupJetCorrectors() {

  vector<string> v_jptjetcorrs;
  vector<string> v_pfjetcorrs;

  v_jptjetcorrs.push_back("../CORE/jetcorr/START38_V13_AK5JPT_L2Relative.txt");
  v_jptjetcorrs.push_back("../CORE/jetcorr/START38_V13_AK5JPT_L3Absolute.txt");

  v_pfjetcorrs.push_back("../CORE/jetcorr/START38_V13_AK5PF_L2Relative.txt");
  v_pfjetcorrs.push_back("../CORE/jetcorr/START38_V13_AK5PF_L3Absolute.txt");

  vector<string> v_jptjetcorrs_wResidual = v_jptjetcorrs;
  vector<string> v_pfjetcorrs_wResidual = v_pfjetcorrs;

  v_jptjetcorrs_wResidual.push_back("../CORE/jetcorr/START38_V13_AK5JPT_L2L3Residual.txt");
  v_pfjetcorrs_wResidual.push_back("../CORE/jetcorr/START38_V13_AK5PF_L2L3Residual.txt");


  jptL2L3Corr = makeJetCorrector(v_jptjetcorrs);
  pfL2L3Corr = makeJetCorrector(v_pfjetcorrs);

  jptL2L3ResidualCorr = makeJetCorrector(v_jptjetcorrs_wResidual);
  pfL2L3ResidualCorr = makeJetCorrector(v_pfjetcorrs_wResidual);


}
*/


double topAFB_looper::JERsf(double eta)
{

    eta = fabs(eta);
    double SF = 1.;

    if(eta<0.5) SF = 1.052;
    else if(eta<1.1) SF = 1.057;
    else if(eta<1.7) SF = 1.096;
    else if(eta<2.3) SF = 1.134;
    else SF = 1.288;

    return SF;

}

double topAFB_looper::triggerEff(const int hypIdx, bool scaleTrigSFup, bool scaleTrigSFdown)
{
    LorentzVector lt_p4  = hyp_lt_p4()[hypIdx];
    LorentzVector ll_p4  = hyp_ll_p4()[hypIdx];
    float lt_pt = lt_p4.Pt();
    float ll_pt = ll_p4.Pt();
    float lt_eta = fabs(lt_p4.Eta());
    float ll_eta = fabs(ll_p4.Eta());
    double weight_lt;
    double weight_ll;
    double weight;
    weight_lt = 1.0;
    weight_ll = 1.0;
    weight = 1.0;
    int id_lt = fabs(hyp_lt_id()[hypIdx]);
    int id_ll = fabs(hyp_ll_id()[hypIdx]);
    //reference AN2011/456 v2
    if (id_lt == 11)
    {
        if (lt_eta >= 0 && lt_eta < 1.5)
        {
            if (lt_pt > 20.0 && lt_pt <= 30.0) weight_lt = 0.9849;
            else if (lt_pt > 30.0) weight_lt = 0.9928;
        }
        else if (lt_eta >= 1.5)
        {
            if (lt_pt > 20.0 && lt_pt <= 30.0) weight_lt = 0.9774;
            else if (lt_pt > 30.0) weight_lt = 0.9938;
        }
    }

    if (id_ll == 11)
    {
        if (ll_eta >= 0 && ll_eta < 1.5)
        {
            if (ll_pt > 20.0 && ll_pt <= 30.0) weight_ll = 0.9923;
            else if (ll_pt > 30.0) weight_ll = 0.9948;
        }
        else if (ll_eta >= 1.5)
        {
            if (ll_pt > 20.0 && ll_pt <= 30.0) weight_ll = 0.9953;
            else if (ll_pt > 30.0) weight_ll = 0.9956;
        }
    }

    if (id_lt == 13)
    {
        if (lt_eta > 0 && lt_eta < 0.8)
        {
            if (lt_pt > 20.0 && lt_pt <= 30.0) weight_lt = 0.9648;
            else if (lt_pt > 30.0) weight_lt = 0.9666;
        }
        else if (lt_eta >= 0.8 && lt_eta < 1.2)
        {
            if (lt_pt > 20.0 && lt_pt <= 30.0) weight_lt = 0.9516;
            else if (lt_pt > 30.0) weight_lt = 0.9521;
        }
        else if (lt_eta >= 1.2 && lt_eta < 2.1)
        {
            if (lt_pt > 20.0 && lt_pt <= 30.0) weight_lt = 0.9480;
            else if (lt_pt > 30.0) weight_lt = 0.9485;
        }
        else if (lt_eta >= 2.1)
        {
            if (lt_pt > 20.0 && lt_pt <= 30.0) weight_lt = 0.8757;
            else if (lt_pt > 30.0) weight_lt = 0.8772;
        }
    }

    if (id_ll == 13)
    {
        if (ll_eta > 0 && ll_eta < 0.8)
        {
            if (ll_pt > 20.0 && ll_pt <= 30.0) weight_ll = 0.9655;
            else if (ll_pt > 30.0) weight_ll = 0.9670;
        }
        else if (ll_eta >= 0.8 && ll_eta < 1.2)
        {
            if (ll_pt > 20.0 && ll_pt <= 30.0) weight_ll = 0.9535;
            else if (ll_pt > 30.0) weight_ll = 0.9537;
        }
        else if (ll_eta >= 1.2 && ll_eta < 2.1)
        {
            if (ll_pt > 20.0 && ll_pt <= 30.0) weight_ll = 0.9558;
            else if (ll_pt > 30.0) weight_ll = 0.9530;
        }
        else if (ll_eta >= 2.1)
        {
            if (ll_pt > 20.0 && ll_pt <= 30.0) weight_ll = 0.9031;
            else if (ll_pt > 30.0) weight_ll = 0.8992;
        }
    }
    
    if ( scaleTrigSFup ) {
      if (id_lt == 11) weight_lt*=1.02;
      if (id_ll == 11) weight_ll*=1.02;
      if (id_lt == 13) weight_lt*=0.98;
      if (id_ll == 13) weight_ll*=0.98;
    } else if ( scaleTrigSFdown )    {
      if (id_lt == 11) weight_lt*=0.98;
      if (id_ll == 11) weight_ll*=0.98;
      if (id_lt == 13) weight_lt*=1.02;
      if (id_ll == 13) weight_ll*=1.02;
    }
    
    weight = weight_lt * weight_ll;
    return weight;
}

double topAFB_looper::getBFRWeight(const int hypIdx, vector<LorentzVector> &v_goodNonBtagJets_p4, vector<LorentzVector> &v_goodBtagJets_p4,  bool isData)
{

    LorentzVector lt_p4  = hyp_lt_p4()[hypIdx];
    LorentzVector ll_p4  = hyp_ll_p4()[hypIdx];
    int nbtag_jet = v_goodBtagJets_p4.size();
    if (nbtag_jet > 2) return -9999;

    double mass_ltb, mass_llb, nonb_pt, nonb_eta;

    if (nbtag_jet == 2 )
    {
        //fill only events where >=1 of the btagged jets does not match a gen-level b for MC closure test
        for (unsigned int i = 0; i < 2; i++)
        {
            mass_ltb = (lt_p4 + v_goodBtagJets_p4.at(i)).M();
            mass_llb = (ll_p4 + v_goodBtagJets_p4.at(i)).M();
            if (! (mass_ltb > 170 && mass_llb > 170) ) return -9999.;
        }
    }

    if (nbtag_jet == 1 )
    {
        for (unsigned int i = 0; i < nbtag_jet; i++)
        {
            mass_ltb = (lt_p4 + v_goodBtagJets_p4.at(i)).M();
            mass_llb = (ll_p4 + v_goodBtagJets_p4.at(i)).M();
            if (! (mass_ltb > 170 && mass_llb > 170) ) return -9999.;
        }
    }

    double bjetfr = 1.0;

    double weight = 0.0;
    int nFO = 0;  //FOs are non-b jets that give mass_ltb > 170 && mass_llb > 170
    bool isFO = false;
    vector<double> weight_tmp ;
    for (unsigned int i = 0; i < v_goodNonBtagJets_p4.size(); i++)
    {
        isFO = false;
        mass_ltb = (lt_p4 + v_goodNonBtagJets_p4.at(i)).M();
        mass_llb = (ll_p4 + v_goodNonBtagJets_p4.at(i)).M();
        nonb_pt =  v_goodNonBtagJets_p4.at(i).Pt();
        nonb_eta =  fabs(v_goodNonBtagJets_p4.at(i).Eta());
        //to keep within the range of the function
        if (nonb_pt > 499.0)
        {
            nonb_pt =  499.0;
        }
        if (nonb_eta > 2.3)
        {
            nonb_eta =  2.3;
        }

        bjetfr =  getMisTagRate( nonb_pt, nonb_eta  , "TCHEM");
        //not needed when the mistag events in MC are being weighted by the MisTagSF
        //if(!isData) bjetfr/=getMisTagSF( nonb_pt, nonb_eta  , "TCHEM");

        if (mass_ltb > 170 && mass_llb > 170 )
        {
            nFO++;
            isFO = true;
        }
        if (isFO)
        {
            weight_tmp.push_back(bjetfr / (1 - bjetfr));
        }
    }


    if (nbtag_jet == 0 )
    {
        if (nFO >= 2) weight = weight_tmp.at(0) * weight_tmp.at(1);
        if (nFO >= 3) weight += weight_tmp.at(0) * weight_tmp.at(2) + weight_tmp.at(1) * weight_tmp.at(2);
        if (nFO >= 4 ) weight += weight_tmp.at(0) * weight_tmp.at(3) + weight_tmp.at(1) * weight_tmp.at(3) + weight_tmp.at(2) * weight_tmp.at(3);
        if (nFO >= 5 ) weight += weight_tmp.at(0) * weight_tmp.at(4) + weight_tmp.at(1) * weight_tmp.at(4) + weight_tmp.at(2) * weight_tmp.at(4) + weight_tmp.at(3) * weight_tmp.at(4);
        if (nFO >= 6 ) weight += weight_tmp.at(0) * weight_tmp.at(5) + weight_tmp.at(1) * weight_tmp.at(5) + weight_tmp.at(2) * weight_tmp.at(5) + weight_tmp.at(3) * weight_tmp.at(5) + weight_tmp.at(4) * weight_tmp.at(5);
        if (nFO >= 7 ) weight += weight_tmp.at(0) * weight_tmp.at(6) + weight_tmp.at(1) * weight_tmp.at(6) + weight_tmp.at(2) * weight_tmp.at(6) + weight_tmp.at(3) * weight_tmp.at(6) + weight_tmp.at(4) * weight_tmp.at(6) + weight_tmp.at(5) * weight_tmp.at(6);
    }


    if (nbtag_jet == 1 )
    {
        if (nFO >= 1) weight = weight_tmp.at(0);
        if (nFO >= 2) weight += weight_tmp.at(1);
        if (nFO >= 3) weight += weight_tmp.at(2);
        if (nFO >= 4 ) weight += weight_tmp.at(3);
        if (nFO >= 5 ) weight += weight_tmp.at(4);
        if (nFO >= 6 ) weight += weight_tmp.at(5);
        if (nFO >= 7 ) weight += weight_tmp.at(6);
    }

    //there are no events with >6 FOs (3.23/fb)

    //if (weight !=1 )  cout <<"nbtag_jet, nFO, weight =  " <<nbtag_jet<<" ,  "<<nFO<<" ,  "<<weight <<endl;

    if (nFO > 1 && nbtag_jet == 0 ) return weight;

    if (nFO > 0 && nbtag_jet == 1 ) return weight;

    if (nbtag_jet >= 2 ) return 1.0;

    if (nFO > 1 - nbtag_jet) cout << "something went wrong " << nbtag_jet << " ,  " << nFO << " ,  " << weight << endl;
    return -9999.;
}

// *****************************************************************
//get the FR weight
// *****************************************************************


double topAFB_looper::getFRWeight(const int hypIdx, SimpleFakeRate *mufr, SimpleFakeRate *elfr, FREnum frmode, bool isData)
{

    //std::cout<<"Called topAFB_looper::getFRWeight"<<std::endl;

    bool  estimateQCD   = false;
    bool  estimateWJets = false;

    if ( frmode == e_qcd )
    {
        estimateQCD   = true;
        estimateWJets = false;
    }
    else if ( frmode == e_wjets )
    {
        estimateQCD   = false;
        estimateWJets = true;
    }
    else
    {
        std::cout << "topAFB_looper::getFRWeight: bad FR mode given, fix this!" << std::endl;
        return -9999.;
    }
    if (hyp_type()[hypIdx] == 0)
    {

        bool isGoodMut = false;
        bool isGoodMul = false;
        bool isFOMut   = false;
        bool isFOMul   = false;

        unsigned int iMut = hyp_lt_index()[hypIdx];
        unsigned int iMul = hyp_ll_index()[hypIdx];

        if ( muonId( iMut , OSGeneric_v3 ) )
        {
            isGoodMut = true;
        }
        if ( muonId( iMul , OSGeneric_v3 ) )
        {
            isGoodMul = true;
        }
        if ( muonId( iMut , OSGeneric_v3_FO ) )
        {
            isFOMut = true;
        }
        if ( muonId( iMul , OSGeneric_v3_FO ) )
        {
            isFOMul = true;
        }

        //for both WJets and QCD, we need both to be FOs at least
        if (!isFOMut || !isFOMul)
            return -9999.;

        //if we want to estimate the fakes for QCD, then we ask that
        //both are not num objects, and that both are FO
        if (estimateQCD)
        {

            //if at least one is a Numerator lepton, we return
            if ( isGoodMut || isGoodMul)
                return -9999.;

            double FRMut = mufr->getFR(mus_p4()[iMut].pt(), mus_p4()[iMut].eta());
            double FRMul = mufr->getFR(mus_p4()[iMul].pt(), mus_p4()[iMul].eta());
            return (FRMut / (1 - FRMut)) * (FRMul / (1 - FRMul));
        }
        else if (estimateWJets)
        {

            //need one to be a Numerator lepton, and the other to be FO but not num
            if ( isGoodMut && !isGoodMul && isFOMul)
            {
                double FR = mufr->getFR(mus_p4()[iMul].pt(), mus_p4()[iMul].eta());
                //cout << "mm, FR and FR/(1-FR) " << FR << ", " << FR/(1-FR) << endl;
                return FR / (1 - FR);
            }

            //check the other muon
            if ( isGoodMul && !isGoodMut && isFOMut)
            {
                double FR = mufr->getFR(mus_p4()[iMut].pt(), mus_p4()[iMut].eta());
                //cout << "mm, FR and FR/(1-FR) " << FR << ", " << FR/(1-FR) << endl;
                return FR / (1 - FR);
            }
        }//estimate WJets
        return -9999.;
    }//mumu case


    //now we do the ee case
    if (hyp_type()[hypIdx] == 3)
    {

        unsigned int iElt = hyp_lt_index()[hypIdx];
        unsigned int iEll = hyp_ll_index()[hypIdx];

        bool isGoodElt = false;
        bool isGoodEll = false;
        bool isFOElt   = false;
        bool isFOEll   = false;

        if ( pass_electronSelection( iElt , electronSelection_el_OSV3 ) )
        {
            isGoodElt = true;
        }
        if ( pass_electronSelection( iEll , electronSelection_el_OSV3 ) )
        {
            isGoodEll = true;
        }
        if ( pass_electronSelection( iElt , electronSelection_el_OSV3_FO ) )
        {
            isFOElt   = true;
        }
        if ( pass_electronSelection( iEll , electronSelection_el_OSV3_FO ) )
        {
            isFOEll   = true;
        }

        //for both WJets and QCD, we need both to be FOs at least
        //if both are good, we continue
        if ( !isFOElt || !isFOEll)
            return -9999.;

        if (estimateQCD)
        {

            //if at least one is a Numerator object, then we return -9999.
            if ( isGoodElt || isGoodEll)
                return -9999.;

            double FRElt = elfr->getFR(els_p4()[iElt].pt(), els_p4()[iElt].eta());
            double FREll = elfr->getFR(els_p4()[iEll].pt(), els_p4()[iEll].eta());
            //      sumfr = sumfr +  (FRElt/(1-FRElt))*(FREll/(1-FREll));
            return (FRElt / (1 - FRElt)) * (FREll / (1 - FREll));
        }
        else if (estimateWJets)
        {

            if (isGoodElt && !isGoodEll && isFOEll)
            {
                double FR = elfr->getFR(els_p4()[iEll].pt(), els_p4()[iEll].eta());
                //cout << "ee, FR and FR/(1-FR) " << FR << ", " << FR/(1-FR) << endl;
                return FR / (1 - FR);
            }
            //check the other electron
            if (isGoodEll && !isGoodElt && isFOElt)
            {
                double FR = elfr->getFR(els_p4()[iElt].pt(), els_p4()[iElt].eta());
                //cout << "ee, FR and FR/(1-FR) " << FR << ", " << FR/(1-FR) << endl;
                return FR / (1 - FR);
            }
            return -9999.;
        }//estimateWJets

    }//ee case

    if (hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2)
    {
        int iEl = 0;
        int iMu = 0;
        /*
        if(hyp_type()[hypIdx] == 2) {
          iEl = hyp_lt_index()[hypIdx];
          iMu = hyp_ll_index()[hypIdx];
        }
        if (hyp_type()[hypIdx] == 1) {
          iEl = hyp_ll_index()[hypIdx];
          iMu = hyp_lt_index()[hypIdx];
        }
        */

        if     ( abs(hyp_ll_id()[hypIdx]) == 11 && abs(hyp_lt_id()[hypIdx]) == 13 )
        {
            iEl = hyp_ll_index()[hypIdx];
            iMu = hyp_lt_index()[hypIdx];
        }
        else if ( abs(hyp_ll_id()[hypIdx]) == 13 && abs(hyp_lt_id()[hypIdx]) == 11 )
        {
            iEl = hyp_lt_index()[hypIdx];
            iMu = hyp_ll_index()[hypIdx];
        }
        else
        {
            cout << "ID ll " << hyp_ll_id()[hypIdx] << endl;
            cout << "ID lt " << hyp_lt_id()[hypIdx] << endl;
            cout << "Error in getFRWeight, quitting!" << endl;
            exit(0);
        }

        bool isGoodEl = false;
        bool isFOEl   = false;
        bool isGoodMu = false;
        bool isFOMu   = false;

        if ( pass_electronSelection( iEl , electronSelection_el_OSV3 ) )
        {
            isGoodEl = true;
        }
        if ( muonId( iMu , OSGeneric_v3 ) )
        {
            isGoodMu = true;
        }
        if ( pass_electronSelection( iEl , electronSelection_el_OSV3_FO ) )
        {
            isFOEl = true;
        }
        if ( muonId( iMu , OSGeneric_v3_FO ) )
        {
            isFOMu = true;
        }

        //if either fail FO, return!!!
        if (!isFOMu || !isFOEl)
            return -9999.;

        if (estimateQCD )
        {

            //if at least one is a numerator, then we fail
            if (isGoodMu || isGoodEl)
                return -9999.;

            double FRMu = mufr->getFR(mus_p4()[iMu].pt(), mus_p4()[iMu].eta());
            double FREl = elfr->getFR(els_p4()[iEl].pt(), els_p4()[iEl].eta());
            return FRMu * FREl / (1 - FRMu) / (1 - FREl);
            //if we get here somehow, we're in trouble
            cout << "We have gotten to line: " << __LINE__  << " in the FR code.";
            cout << "We should never get here, so something is wrong with our logic!!" << endl;
        }
        else if (estimateWJets)
        {

            //need one to be a numerator lepton and the other to be a FO
            if (isGoodMu && !isGoodEl && isFOEl)
            {
                double FR = elfr->getFR(els_p4()[iEl].pt(), els_p4()[iEl].eta());
                //cout << "emu, el FR, FR/(1-FR): " << FR << ", " << FR/(1-FR) << endl;
                return FR / (1 - FR);
            }

            if (isGoodEl && !isGoodMu && isFOMu)
            {
                double FR = mufr->getFR(mus_p4()[iMu].pt(), mus_p4()[iMu].eta());
                //cout << "emu, mu FR, FR/(1-FR): " << FR << ", " << FR/(1-FR) << endl;
                return FR / (1 - FR);
            }
            return -9999.;
        }
    } //emu case

    return -9999.;
}


topAFB_looper::topAFB_looper()
{
    applyNoCuts = false;
    getVtxDistOnly = false;
    usePtGt2020 = false;
    usePtGt2010 = false;
    excludePtGt2020 = false;
    applylepIDCuts = false;
    applyFOv1Cuts = false;
    applyFOv2Cuts = false;
    applyFOv3Cuts = false;
    applylepIsoCuts = false;
    applylepLooseIsoCuts = false;
    applyTriggers = false;
    vetoZmass = false;
    requireZmass = false;
    hypDisamb = false;
    useCorMET = false;
    usetcMET = false;
    usetcMET35X                = false;
    usepfMET = false;
    vetoMET = false;
    vetoMET50 = false;
    vetoProjectedMET = false;
    usejptJets = false;
    usecaloJets = false;
    usepfJets = false;
    veto1Jet = false;
    veto2Jets = false;
    requireEcalEls = false;
    useOS = false;
    useSS = false;
    applyTopPtWeighting = false;
    applyLeptonPtWeighting = false;
    applyJetPtWeighting = false;
    useReweightingUncorrelated = false;
    useReweightingLeadingObject = false;
    applyAlignmentCorrection   = false;
    vetoHypMassLt10            = false;
    vetoHypMassLt12            = false;
    scaleJESMETUp              = false;
    scaleJESMETDown            = false;
    scaleJER                   = false;
    scaleLeptonEnergyUp = false;
    scaleLeptonEnergyDown = false;
    scaleBTAGSFup = false;
    scaleBTAGSFdown = false;
    scaleTrigSFup = false;
    scaleTrigSFdown = false;
    noVertexReweighting = false;
    weighttaudecay = false;
    estimateQCD = false;
    estimateWJets = false;
    requireBTag                = false;
    require2BTag               = false;
    sortJetCandidatesbyPt      = false;
    sortJetCandidatesbyDR      = false;
    applyLeptonJetInvMassCut450 = false;
    doBFR  = false;
    requireExact2BTag           = false;
    applyTopSystEta              = false;
    globalJESRescale = 1.;

    jptL2L3Corr = NULL;
    pfL2L3Corr = NULL;
    jptL2L3ResidualCorr = NULL;
    pfL2L3ResidualCorr = NULL;
    d_llsol = new ttdilepsolve;
}


topAFB_looper::~topAFB_looper()
{
    delete babyFile_;
    delete babyTree_;
    delete d_llsol;
}
void topAFB_looper::ScanChain(TChain *chain, vector<TString> v_Cuts, string prefix,
                              bool doFRestimation, float lumi , float kFactor , bool verbose , FREnum frmode, double AMWTmass )
{


  // reset JES scale variable in ScanChain
  globalJESRescale = 1.;
  
    //deal with the cuts
    applyNoCuts = find(v_Cuts.begin(), v_Cuts.end(), "applyNoCuts") != v_Cuts.end();
    getVtxDistOnly = find(v_Cuts.begin(), v_Cuts.end(), "getVtxDistOnly") != v_Cuts.end();
    usePtGt2020 = find(v_Cuts.begin(), v_Cuts.end(), "usePtGt2020") != v_Cuts.end();
    usePtGt2010 = find(v_Cuts.begin(), v_Cuts.end(), "usePtGt2010") != v_Cuts.end();
    excludePtGt2020 = find(v_Cuts.begin(), v_Cuts.end(), "excludePtGt2020") != v_Cuts.end();
    applylepIDCuts = find(v_Cuts.begin(), v_Cuts.end(), "applylepIDCuts") != v_Cuts.end();
    applyFOv1Cuts = find(v_Cuts.begin(), v_Cuts.end(), "applyFOv1Cuts") != v_Cuts.end();
    applyFOv2Cuts = find(v_Cuts.begin(), v_Cuts.end(), "applyFOv2Cuts") != v_Cuts.end();
    applyFOv3Cuts = find(v_Cuts.begin(), v_Cuts.end(), "applyFOv3Cuts") != v_Cuts.end();
    applylepIsoCuts = find(v_Cuts.begin(), v_Cuts.end(), "applylepIsoCuts") != v_Cuts.end();
    applylepLooseIsoCuts = find(v_Cuts.begin(), v_Cuts.end(), "applylepLooseIsoCuts") != v_Cuts.end();
    applyTriggers = find(v_Cuts.begin(), v_Cuts.end(), "applyTriggers") != v_Cuts.end();
    vetoZmass = find(v_Cuts.begin(), v_Cuts.end(), "vetoZmass") != v_Cuts.end();
    requireZmass = find(v_Cuts.begin(), v_Cuts.end(), "requireZmass") != v_Cuts.end();
    hypDisamb = find(v_Cuts.begin(), v_Cuts.end(), "hypDisamb") != v_Cuts.end();
    useCorMET = find(v_Cuts.begin(), v_Cuts.end(), "useCorMET") != v_Cuts.end();
    usetcMET = find(v_Cuts.begin(), v_Cuts.end(), "usetcMET") != v_Cuts.end();
    usetcMET35X = find(v_Cuts.begin(), v_Cuts.end(), "usetcMET35X") != v_Cuts.end();
    usepfMET = find(v_Cuts.begin(), v_Cuts.end(), "usepfMET") != v_Cuts.end();
    vetoMET = find(v_Cuts.begin(), v_Cuts.end(), "vetoMET") != v_Cuts.end();
    vetoMET50 = find(v_Cuts.begin(), v_Cuts.end(), "vetoMET50") != v_Cuts.end();
    vetoProjectedMET = find(v_Cuts.begin(), v_Cuts.end(), "vetoProjectedMET") != v_Cuts.end();
    usecaloJets = find(v_Cuts.begin(), v_Cuts.end(), "usecaloJets") != v_Cuts.end();
    usejptJets = find(v_Cuts.begin(), v_Cuts.end(), "usejptJets") != v_Cuts.end();
    usepfJets = find(v_Cuts.begin(), v_Cuts.end(), "usepfJets") != v_Cuts.end();
    veto1Jet = find(v_Cuts.begin(), v_Cuts.end(), "veto1Jet") != v_Cuts.end();
    veto2Jets = find(v_Cuts.begin(), v_Cuts.end(), "veto2Jets") != v_Cuts.end();
    requireEcalEls = find(v_Cuts.begin(), v_Cuts.end(), "requireEcalEls") != v_Cuts.end();
    useOS = find(v_Cuts.begin(), v_Cuts.end(), "useOS") != v_Cuts.end();
    useSS = find(v_Cuts.begin(), v_Cuts.end(), "useSS") != v_Cuts.end();
    applyTopPtWeighting = find(v_Cuts.begin(), v_Cuts.end(), "applyTopPtWeighting" ) != v_Cuts.end();
    applyLeptonPtWeighting = find(v_Cuts.begin(), v_Cuts.end(), "applyLeptonPtWeighting" ) != v_Cuts.end();
    applyJetPtWeighting = find(v_Cuts.begin(), v_Cuts.end(), "applyJetPtWeighting" ) != v_Cuts.end();
    useReweightingUncorrelated = find(v_Cuts.begin(), v_Cuts.end(), "useReweightingUncorrelated" ) != v_Cuts.end();
    useReweightingLeadingObject = find(v_Cuts.begin(), v_Cuts.end(), "useReweightingLeadingObject" ) != v_Cuts.end();
    applyAlignmentCorrection = find(v_Cuts.begin(), v_Cuts.end(), "applyAlignmentCorrection" ) != v_Cuts.end();
    vetoHypMassLt10               = find(v_Cuts.begin(), v_Cuts.end(), "vetoHypMassLt10"                  ) != v_Cuts.end();
    vetoHypMassLt12               = find(v_Cuts.begin(), v_Cuts.end(), "vetoHypMassLt12"                  ) != v_Cuts.end();
    scaleJESMETUp = find(v_Cuts.begin(), v_Cuts.end(), "scaleJESMETUp"                    ) != v_Cuts.end();
    scaleJESMETDown = find(v_Cuts.begin(), v_Cuts.end(), "scaleJESMETDown"                  ) != v_Cuts.end();
    scaleJER = find(v_Cuts.begin(), v_Cuts.end(), "scaleJER"                  ) != v_Cuts.end();
    scaleLeptonEnergyUp = find(v_Cuts.begin(), v_Cuts.end(), "scaleLeptonEnergyUp") != v_Cuts.end();
    scaleLeptonEnergyDown = find(v_Cuts.begin(), v_Cuts.end(), "scaleLeptonEnergyDown") != v_Cuts.end();
    scaleBTAGSFup = find(v_Cuts.begin(), v_Cuts.end(), "scaleBTAGSFup") != v_Cuts.end();
    scaleBTAGSFdown = find(v_Cuts.begin(), v_Cuts.end(), "scaleBTAGSFdown") != v_Cuts.end();
    scaleTrigSFup = find(v_Cuts.begin(), v_Cuts.end(), "scaleTrigSFup") != v_Cuts.end();
    scaleTrigSFdown = find(v_Cuts.begin(), v_Cuts.end(), "scaleTrigSFdown") != v_Cuts.end();
    noVertexReweighting = find(v_Cuts.begin(), v_Cuts.end(), "noVertexReweighting") != v_Cuts.end();
    weighttaudecay = find(v_Cuts.begin(), v_Cuts.end(), "weighttaudecay") != v_Cuts.end();
    estimateQCD    = find(v_Cuts.begin(), v_Cuts.end(), "estimateQCD") != v_Cuts.end();
    estimateWJets = find(v_Cuts.begin(), v_Cuts.end(), "estimateWJets") != v_Cuts.end();
    requireBTag                   = find(v_Cuts.begin(), v_Cuts.end(), "requireBTag"                      ) != v_Cuts.end();
    require2BTag                   = find(v_Cuts.begin(), v_Cuts.end(), "require2BTag"                      ) != v_Cuts.end();
    sortJetCandidatesbyPt         = find(v_Cuts.begin(), v_Cuts.end(), "sortJetCandidatesbyPt"                      ) != v_Cuts.end();
    sortJetCandidatesbyDR         = find(v_Cuts.begin(), v_Cuts.end(), "sortJetCandidatesbyDR"                      ) != v_Cuts.end();
    applyLeptonJetInvMassCut450   = find(v_Cuts.begin(), v_Cuts.end(), "applyLeptonJetInvMassCut450"                ) != v_Cuts.end();
    generalLeptonVeto   = find(v_Cuts.begin(), v_Cuts.end(), "generalLeptonVeto"                ) != v_Cuts.end();
    applyHTCut   = find(v_Cuts.begin(), v_Cuts.end(), "applyHTCut"                ) != v_Cuts.end();
    matchLeptonJetbyMaxDR   = find(v_Cuts.begin(), v_Cuts.end(), "matchLeptonJetbyMaxDR"                ) != v_Cuts.end();
    BTagAlgTCHE = find(v_Cuts.begin(), v_Cuts.end(), "BTagAlgTCHE"                ) != v_Cuts.end();
    createBabyNtuples =  find(v_Cuts.begin(), v_Cuts.end(), "createBabyNtuples"                ) != v_Cuts.end();
    doBFR =  find(v_Cuts.begin(), v_Cuts.end(), "doBFR"                ) != v_Cuts.end();
    requireExact2BTag    = find(v_Cuts.begin(), v_Cuts.end(), "requireExact2BTag"                      ) != v_Cuts.end();
    applyTopSystEta       = find(v_Cuts.begin(), v_Cuts.end(), "applyTopSystEta"                      ) != v_Cuts.end();
    // top mass

    //std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
    //FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3;

    //jetcorr_filenames_pfL1FastJetL2L3.clear();

    string pfUncertaintyFile;

    if ( TString(prefix).Contains("data") )
    {
        //jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L1FastJet.txt");
        //jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L2Relative.txt");
        //jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L3Absolute.txt");
        //jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_R_42_V23_AK5PF_L2L3Residual.txt");

        pfUncertaintyFile = "../CORE/jetcorr/jetCorrectionFiles/GR_R_42_V23_AK5PF_Uncertainty.txt";
    }
    else
    {
        //jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/DESIGN42_V17_AK5PF_L1FastJet.txt");
        //jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/DESIGN42_V17_AK5PF_L2Relative.txt");
        //jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/DESIGN42_V17_AK5PF_L3Absolute.txt");

        pfUncertaintyFile = "../CORE/jetcorr/jetCorrectionFiles/DESIGN42_V17_AK5PF_Uncertainty.txt";
    }

    //jet_corrector_pfL1FastJetL2L3  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);

    JetCorrectionUncertainty *pfUncertainty   = new JetCorrectionUncertainty( pfUncertaintyFile );

    //if ( scaleJESMETUp) globalJESRescale = 1.075 ;
    //else if (scaleJESMETDown)globalJESRescale = 0.925;



    already_seen.clear();

    // Make a baby ntuple
    if (createBabyNtuples && !applyNoCuts)
    {
        TString babyFilename = prefix + ".root";
        MakeBabyNtuple(babyFilename.Data());
    }
    bool isData = false;
    // Set the JSON file
    if (prefix == "data")
    {
        //set_goodrun_file("Cert_160404-163869_7TeV_PromptReco_Collisions11_JSON_goodruns.txt");
        //set_goodrun_file("Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_goodruns.txt");
        // set_goodrun_file_json("Cert_160404_165970_7TeV_May10ReRecoPlusPromptReco_349pb-1.json");
        // set_goodrun_file("Cert_EPSFINAL_May10ReReco_v2_PromptReco_160404_167913_JSON_goodruns.txt");
        //set_goodrun_file_json("Cert_160404-173692_7TeV_PromptReco_Collisions11_JSON.txt");
        //set_goodrun_file_json("Cert_160404-177515_7TeV_PromptReco-May10ReRecov3-ReReco5Augv2.txt");
        set_goodrun_file_json("Cert_160404-180252_7TeV_May10ReRecoV3_ReReco5AugV3_PromptReco_Collisions11_JSON.txt");
        cout << "DATA!!!" << endl;
        isData = true;
    }

    //---------------------------
    // lepton energy systematic
    // vary electron energy by +- 1 sigma (0.3%) for data
    float leptonEnergyScaleFactor = 1.;
    if ( isData ) {
      if ( scaleLeptonEnergyUp ) {
        leptonEnergyScaleFactor     = 1.003;
      }
      if ( scaleLeptonEnergyDown ) {
        leptonEnergyScaleFactor     = 0.997;
      }
    }
        
    // pileup vertex reweighting, switch off with v_baseCuts.push_back("noVertexReweighting") in doAll.C
    if ( !noVertexReweighting ) {
      if ( prefix == "ttdil" || prefix == "ttotr" )
      {
          cout << "using Fall11 vertex weighting" << endl;
          set_vtxreweight_rootfile("vtxreweight_Fall11MC_PUS6_4p7fb_Zselection.root", false);
      }
      else
      {
          cout << "using Summer11 vertex weighting" << endl;
          set_vtxreweight_rootfile("vtxreweight_Summer11MC_PUS4_4p7fb_Zselection.root", false);
      }
    }

    //set_vtxreweight_rootfile("vtxreweight_Summer11MC_PUS4_4p7fb_Zselection.root",false);
    //set_vtxreweight_rootfile("vtxreweight_Fall11MC_PUS6_4p7fb_Zselection.root",false);
    //set_vtxreweight_rootfile("vtxreweight_Summer11MC_PUS4_3p5fb_Zselection.root",false);
    //set_vtxreweight_rootfile("vtxreweight_Spring11MC_336pb_Zselection.root",false);
    //set_vtxreweight_rootfile("vtxreweight_Summer11MC_1160pb_Zselection.root",false);
    //set_vtxreweight_rootfile_tprime("vtxreweight_Summer11MC_1160pb_Zselection.root",false);

    int nchs = NCHANNELS;
    int nhists = NHISTS;
    bookHistos(prefix.c_str(), nchs, nhists);


    //instantiate SimpleFakeRate class for electrons and muons
    SimpleFakeRate *mufr = 0;
    //this is the default, can change it below if needed
    SimpleFakeRate *elfr = 0;

    if (doFRestimation)
    {
        std::cout << "**************************" << std::endl;
        std::cout << "Running FR application job" << std::endl;
        std::cout << "**************************" << std::endl;

        if (prefix == "data")
        {
            std::cout << "Using data derived FR files" << std::endl;
            mufr = new SimpleFakeRate("fr_os7June2011.root", "fr_mu_OSGV3" );
            elfr = new SimpleFakeRate("fr_os7June2011.root", "fr_el_OSGV3" );
        }
        else
        {
            std::cout << "Using data derived FR files" << std::endl;
            std::cout << "CURRENTLY USING DATA FR FOR MC FIXME!!!!!" << std::endl;
            mufr = new SimpleFakeRate("fr_os7June2011.root", "fr_mu_OSGV3" );
            elfr = new SimpleFakeRate("fr_os7June2011.root", "fr_el_OSGV3" );
        }
    }

    //--------------------------
    // File and Event Loop
    //---------------------------
    TObjArray *listOfFiles = chain->GetListOfFiles();
    unsigned int nEventsChain = 0;
    unsigned int nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;
    double nEvents_noCuts = 0;
    double nEvents_noCuts_dil = 0;
    ULong64_t nEvents_noCuts_novtxweight = 0;
    ULong64_t nEvents_noCuts_novtxweight_dil = 0;
    double nSelectedEvents = 0;
    unsigned int npreSelectedEvents = 0;
    unsigned int npreSelectedEvents_genmatch1 = 0;
    unsigned int npreSelectedEvents_genmatch2 = 0;
    TIter fileIter(listOfFiles);
    map<int, int> m_events;
    while (TChainElement *currentFile = (TChainElement *)fileIter.Next() )
    {
        TString filename = currentFile->GetTitle();

        TFile f(filename.Data());
        TTree *tree = (TTree *)f.Get("Events");
        TTreeCache::SetLearnEntries(10);
        tree->SetCacheSize(128 * 1024 * 1024);
        cms2.Init(tree);
        unsigned int nEvents = tree->GetEntries();
        bool print_evt_weight = true;
        for (unsigned int event = 0; event < nEvents; ++event)
        {
            // Event Loop
            tree->LoadTree(event);
            cms2.GetEntry(event);


            if (print_evt_weight && !isData)
            {
                cout << "Event Weight = " << kFactor *evt_scale1fb() * lumi << endl;
                print_evt_weight = false;
            }
            // skip duplicates
            if ( isData )
            {
                DorkyEventIdentifier id = { evt_run(), evt_event(), evt_lumiBlock() };
                if (is_duplicate(id) )
                {
                    //  cout << "Found duplicate! " << evt_dataset()<< " " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
                    continue;
                }
            }

            int myType = 2;
            unsigned int jetBin = 2;
            int ndavtx = 0;
            double weight = 1.0;
            float lepPlus_costheta_cms , lep_azimuthal_asymmetry , lep_azimuthal_asymmetry_2 , lep_charge_asymmetry , lep_pseudorap_diff , top_costheta_cms;
            float lepMinus_costheta_cms;
            float top_pseudorapiditydiff_cms , top_rapiditydiff_Marco , top_rapiditydiff_cms , top_spin_correlation , ttRapidity ,ttRapidity2, tt_mass , tt_mass_nojetsmear , tt_pT , tt_pT_nojetsmear, massllbb;
            float m_top = -999.0;
            float m_top_nojetsmear = -999.0;
            double mass_ltb, mass_llb;

            vector <float> AMWTweight, AMWTweight_nojetsmear;
            vector <TLorentzVector> top1_p4, top2_p4, top1_nojetsmear_p4, top2_nojetsmear_p4;
            TLorentzVector cms, cms_nojetsmear, lepPlus, lepMinus, jet1, jet2;
            int Nsolns = -999;
            int imaxAMWTweight = -999;
            float maxAMWTweight = -999;
            double sumAMWTweight = -999;
            float aveAMWTweight = -999;
            bool useOnlyMaxWsoln = false;

            float ndavtxweight = 1.;
            if ( !noVertexReweighting ) {
              ndavtxweight = vtxweight(isData, true);              
            }

            //if (prefix == "ttdil" || prefix == "ttotr") ndavtxweight = 1.;  //no vtx weighting for fastsim samples

            //get the channels correct
            int nels = 0;
            int nmus = 0;
            int ntaus = 0;
            int nleps = 0;
            if (!isData)
                nleps = leptonGenpCount_lepTauDecays(nels, nmus, ntaus);
            if (!applyNoCuts)
            {
                if (prefix == "ttdil"    &&  nleps != 2) continue;
                if (prefix == "ttotr"    &&  nleps == 2) continue;
                if (prefix == "DYee"     &&  nels != 2) continue;
                if (prefix == "DYmm"     &&  nmus != 2) continue;
                if (prefix == "DYtautau" &&  ntaus != 2) continue;
            }

            if (prefix == "ttdil"    &&  nleps == 2) ++nEventsTotal;
            // Progress feedback to the user
            if (nEventsTotal % 2000 == 0)
            {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1))
                {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                           "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal / (nEventsChain * 0.01));
                    fflush(stdout);
                }
            }//if(nEventsTotal%20000 == 0) {


            //daughter lepton angle in tau rest frame to check if MC is correctly using the tau polarisation
            double weight_taudecay = 1.;
            if (ntaus > 0)
            {

                TLorentzVector lepPlus_gen(0, 0, 0, 0), lepMinus_gen(0, 0, 0, 0), lepPlus_gen_status1(0, 0, 0, 0), lepMinus_gen_status1(0, 0, 0, 0);
                bool lepPlusIsTau = false;
                bool lepPlusIsTauEle = false;
                bool lepPlusIsTauMuo = false;
                bool lepMinusIsTau = false;
                bool lepMinusIsTauEle = false;
                bool lepMinusIsTauMuo = false;

                for (unsigned int i = 0; i < genps_p4().size(); i++)
                {
                    if (genps_status()[i] == 3 && fabs(genps_id_mother()[i]) == 24 )
                    {

                        if (genps_id()[i] == -15 )
                        {
                            //if(genps_id()[i]==-15 && genps_lepdaughter_id()[i].size()==3 ) {   //cosTheta* distribution is only flat for unpolarised tau decay in the simple case of 3 daughters
                            double lepPlus_genX = 0;
                            double lepPlus_genY = 0;
                            double lepPlus_genZ = 0;
                            double lepPlus_genT = 0;

                            for (unsigned int kk = 0; kk < genps_lepdaughter_id()[i].size(); kk++)
                            {
                                int daughterID = abs(genps_lepdaughter_id()[i][kk]);
                                if ( daughterID == 12)
                                {
                                    lepPlusIsTauEle = true;
                                    lepPlusIsTau = true;
                                    break;
                                }
                                if ( daughterID == 14)
                                {
                                    lepPlusIsTauMuo = true;
                                    lepPlusIsTau = true;
                                    break;
                                }
                                // only neutrinos are definitely from leptonic tau decay (see comment for leptonGenpCount_lepTauDecays)
                            }

                            for (unsigned int kk = 0; kk < genps_lepdaughter_id()[i].size(); kk++)
                            {
                                int daughterID = abs(genps_lepdaughter_id()[i][kk]);
                                if ( daughterID == 11 && lepPlusIsTauEle)
                                {
                                    lepPlus_gen_status1.SetXYZT( genps_lepdaughter_p4()[i][kk].x(), genps_lepdaughter_p4()[i][kk].y(), genps_lepdaughter_p4()[i][kk].z(), genps_lepdaughter_p4()[i][kk].t() );
                                    lepPlusIsTauEle = false;  //so that any secondary daughter electrons don't overwrite the first one
                                }
                                if ( daughterID == 13 && lepPlusIsTauMuo)
                                {
                                    lepPlus_gen_status1.SetXYZT( genps_lepdaughter_p4()[i][kk].x(), genps_lepdaughter_p4()[i][kk].y(), genps_lepdaughter_p4()[i][kk].z(), genps_lepdaughter_p4()[i][kk].t() );
                                    lepPlusIsTauMuo = false;  //so that any secondary daughter muons don't overwrite the first one
                                }
                                lepPlus_genX += genps_lepdaughter_p4()[i][kk].x();
                                lepPlus_genY += genps_lepdaughter_p4()[i][kk].y();
                                lepPlus_genZ += genps_lepdaughter_p4()[i][kk].z();
                                lepPlus_genT += genps_lepdaughter_p4()[i][kk].T();
                                lepPlus_gen.SetXYZT(lepPlus_genX, lepPlus_genY, lepPlus_genZ, lepPlus_genT);
                            }

                        }

                        if (genps_id()[i] == 15 )
                        {
                            //if(genps_id()[i]==15 && genps_lepdaughter_id()[i].size()==3 ) {   //cosTheta* distribution is only flat for unpolarised tau decay in the simple case of 3 daughters
                            double lepMinus_genX = 0;
                            double lepMinus_genY = 0;
                            double lepMinus_genZ = 0;
                            double lepMinus_genT = 0;

                            for (unsigned int kk = 0; kk < genps_lepdaughter_id()[i].size(); kk++)
                            {
                                int daughterID = abs(genps_lepdaughter_id()[i][kk]);
                                if ( daughterID == 12)
                                {
                                    lepMinusIsTauEle = true;
                                    lepMinusIsTau = true;
                                    break;
                                }
                                if ( daughterID == 14)
                                {
                                    lepMinusIsTauMuo = true;
                                    lepMinusIsTau = true;
                                    break;
                                }
                                // only neutrinos are definitely from leptonic tau decay (see comment for leptonGenpCount_lepTauDecays)
                            }

                            for (unsigned int kk = 0; kk < genps_lepdaughter_id()[i].size(); kk++)
                            {
                                int daughterID = abs(genps_lepdaughter_id()[i][kk]);
                                if ( daughterID == 11 && lepMinusIsTauEle)
                                {
                                    lepMinus_gen_status1.SetXYZT( genps_lepdaughter_p4()[i][kk].x(), genps_lepdaughter_p4()[i][kk].y(), genps_lepdaughter_p4()[i][kk].z(), genps_lepdaughter_p4()[i][kk].t() );
                                    lepMinusIsTauEle = false;  //so that any secondary daughter electrons don't overwrite the first one
                                }
                                if ( daughterID == 13 && lepMinusIsTauMuo)
                                {
                                    lepMinus_gen_status1.SetXYZT( genps_lepdaughter_p4()[i][kk].x(), genps_lepdaughter_p4()[i][kk].y(), genps_lepdaughter_p4()[i][kk].z(), genps_lepdaughter_p4()[i][kk].t() );
                                    lepMinusIsTauMuo = false;  //so that any secondary daughter muons don't overwrite the first one
                                }
                                lepMinus_genX += genps_lepdaughter_p4()[i][kk].x();
                                lepMinus_genY += genps_lepdaughter_p4()[i][kk].y();
                                lepMinus_genZ += genps_lepdaughter_p4()[i][kk].z();
                                lepMinus_genT += genps_lepdaughter_p4()[i][kk].T();
                                lepMinus_gen.SetXYZT(lepMinus_genX, lepMinus_genY, lepMinus_genZ, lepMinus_genT);
                            }

                        }

                    }
                }

                bool fillEvent = true;
                if (prefix == "ttdil"    &&  nleps != 2) fillEvent = false;
                if (prefix == "ttotr"    &&  nleps == 2) fillEvent = false;
                if (prefix == "DYee"     &&  nels != 2) fillEvent = false;
                if (prefix == "DYmm"     &&  nmus != 2) fillEvent = false;
                if (prefix == "DYtautau" &&  ntaus != 2) fillEvent = false;

                if (lepPlusIsTau)
                {
                    lepPlus_gen_status1.Boost(-lepPlus_gen.BoostVector());
                    double cosTheta_lepPlus_status1 = lepPlus_gen_status1.Vect().Dot(lepPlus_gen.Vect()) / (lepPlus_gen_status1.Vect().Mag() * lepPlus_gen.Vect().Mag());
                    double EoverEmax_lepPlus = 2.*lepPlus_gen_status1.E() / lepPlus_gen.M();
                    if (EoverEmax_lepPlus > 1.) EoverEmax_lepPlus = 1.;
                    //if(EoverEmax_lepPlus>1.) {cout<<"**********x>1********"<<endl; cout<<EoverEmax_lepPlus<<" "<<lepPlus_gen_status1.E()<<" "<<lepPlus_gen.M()<<endl;}
                    //double weight_lepPlus = 1.+cosTheta_lepPlus_status1/3.;  //approximation: cosTheta dependence varies with x
                    double weight_lepPlus = 1. + cosTheta_lepPlus_status1 * ( 1. - 2.*EoverEmax_lepPlus ) / ( 2.*EoverEmax_lepPlus - 3. );
                    //cout<<"P "<<(weighttaudecay?weight_lepPlus:1.)<<endl;
                    if (fillEvent) fillHistos( hlepPlusCosThetaTau_gen,  cosTheta_lepPlus_status1, weighttaudecay ? weight_lepPlus : 1., myType, jetBin);
                    if (fillEvent) fillHistos( hlepPlusxTau_gen,  EoverEmax_lepPlus, weighttaudecay ? weight_lepPlus : 1., myType, jetBin);
                    weight_taudecay *= weight_lepPlus;
                }
                if (lepMinusIsTau)
                {
                    lepMinus_gen_status1.Boost(-lepMinus_gen.BoostVector());
                    double cosTheta_lepMinus_status1 = lepMinus_gen_status1.Vect().Dot(lepMinus_gen.Vect()) / (lepMinus_gen_status1.Vect().Mag() * lepMinus_gen.Vect().Mag());
                    double EoverEmax_lepMinus = 2.*lepMinus_gen_status1.E() / lepMinus_gen.M();
                    if (EoverEmax_lepMinus > 1.) EoverEmax_lepMinus = 1.;
                    //if(EoverEmax_lepMinus>1.) {cout<<"**********x>1********"<<endl; cout<<EoverEmax_lepMinus<<" "<<lepMinus_gen_status1.E()<<" "<<lepMinus_gen.M()<<endl;}
                    //double weight_lepMinus = 1.+cosTheta_lepMinus_status1/3.;  //approximation: cosTheta dependence varies with x
                    double weight_lepMinus = 1. + cosTheta_lepMinus_status1 * ( 1. - 2.*EoverEmax_lepMinus ) / ( 2.*EoverEmax_lepMinus - 3. );
                    //cout<<"M "<<(weighttaudecay?weight_lepMinus:1.)<<endl;
                    if (fillEvent) fillHistos( hlepMinusCosThetaTau_gen,  cosTheta_lepMinus_status1, weighttaudecay ? weight_lepMinus : 1., myType, jetBin);
                    if (fillEvent) fillHistos( hlepMinusxTau_gen,  EoverEmax_lepMinus, weighttaudecay ? weight_lepMinus : 1., myType, jetBin);
                    weight_taudecay *= weight_lepMinus;
                }

            } //ntaus>0
            //cout<<weight_taudecay<<endl;


            if (!applyNoCuts)
            {
                if ( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;
                //if( isData && !goodrun_json(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;
                //did it pass all the good event cuts?
                if ( !cleaning_standardApril2011() )                            continue;

                //skip events with bad els_conv_dist
                bool skipEvent = false;
                for ( unsigned int iEl = 0 ; iEl < els_conv_dist().size() ; ++iEl )
                {
                    if ( els_conv_dist().at(iEl) != els_conv_dist().at(iEl) )
                    {
                        skipEvent = true;
                    }
                    if ( els_sigmaIEtaIEta().at(iEl) != els_sigmaIEtaIEta().at(iEl) )
                    {
                        skipEvent = true;
                    }
                    if ( els_sigmaIEtaIEtaSC().at(iEl) != els_sigmaIEtaIEtaSC().at(iEl) )
                    {
                        skipEvent = true;
                    }
                }
                if ( skipEvent )
                {
                    continue;
                }


                if (verbose)
                    cout << "Event passed all standard event cleaning cuts" << endl;

                float pthat_cutoff = 30.;
                if (prefix == "qcdpt15" && genps_pthat() > pthat_cutoff)
                    continue;


                if ( isData )
                {
                    weight = 1;
                }
                else
                {
		  weight = kFactor * evt_scale1fb() * lumi * ndavtxweight;
		  //weight = kFactor * evt_scale1fb() * lumi ;
                    //negative weights for MC@NLO
                    if (prefix == "ttdil" || prefix == "ttotr") weight = weight * fabs(genps_weight()) / genps_weight();
                    //tau decay cosTheta* weighting
                    if (weighttaudecay && (prefix == "ttdil" || prefix == "ttotr")  && ntaus > 0) weight *= weight_taudecay;
                }
            }//!applyNoCuts

            //splice together the DY samples - if its madgraph, then we do nothing
            if (TString(prefix).Contains("DY") && TString(evt_dataset()).Contains("madgraph") == false)
            {
                bool doNotContinue = false;
                for (unsigned int i = 0; i < genps_p4().size(); i++)
                {
                    if (abs(genps_id()[i]) == 23 && genps_p4()[i].M() > 50.)
                        doNotContinue = true;
                }
                if (doNotContinue)
                    continue;
            }

            if (applyNoCuts)
            {
                if (!isData) weight = kFactor * evt_scale1fb() * lumi;
                //negative weights for MC@NLO
                if (prefix == "ttdil" || prefix == "ttotr") weight = weight * fabs(genps_weight()) / genps_weight();
                //tau decay cosTheta* weighting
                if (weighttaudecay && (prefix == "ttdil" || prefix == "ttotr")  && ntaus > 0) weight *= weight_taudecay;
                for (size_t v = 0; v < cms2.davtxs_position().size(); ++v)
                {
                    if (isGoodDAVertex(v)) ++ndavtx;
                }
                hnVtx[3][3]->Fill(ndavtx, ndavtxweight);
                nEvents_noCuts_novtxweight += 1;
                if (isData) nEvents_noCuts += 1.;
                else  nEvents_noCuts += ndavtxweight;

                if (prefix == "ttdil"    &&  nleps != 2) continue;
                if (prefix == "ttotr"    &&  nleps == 2) continue;
                if (prefix == "DYee"     &&  nels != 2) continue;
                if (prefix == "DYmm"     &&  nmus != 2) continue;
                if (prefix == "DYtautau" &&  ntaus != 2) continue;

                nEvents_noCuts_novtxweight_dil += 1;
                if (isData) nEvents_noCuts_dil += 1.;
                else  nEvents_noCuts_dil += ndavtxweight;
            }

            vector<unsigned int> v_goodHyps;
            v_goodHyps.clear();
            vector<float> v_weights;
            v_weights.clear();

            VofP4 goodLeptons;

            ngoodlep_ = 0;
            ngoodel_  = 0;
            ngoodmu_  = 0;

            if (!applyNoCuts)
            {

                if ( generalLeptonVeto )
                {

                    for ( unsigned int iel = 0 ; iel < els_p4().size(); ++iel )
                    {
                        if ( els_p4().at(iel).pt() < 20 )                                                 continue;
                        if ( !pass_electronSelection( iel , electronSelection_el_OSV3 , false , false ) ) continue;
                        goodLeptons.push_back( els_p4().at(iel) );
                        ngoodel_++;
                        ngoodlep_++;
                    }

                    for ( unsigned int imu = 0 ; imu < mus_p4().size(); ++imu )
                    {
                        if ( mus_p4().at(imu).pt() < 20 )           continue;
                        if ( !muonId( imu , OSGeneric_v3 ))         continue;
                        goodLeptons.push_back( mus_p4().at(imu) );
                        ngoodmu_++;
                        ngoodlep_++;
                    }

                }

                if (ngoodlep_ > 2) continue; //veto the third lepton

                for (unsigned int ihyp = 0; ihyp < hyp_p4().size(); ++ihyp)
                {
                    //InitBabyNtuple();
                    //int mc_match;
                    //if(prefix == "ttdil" ){
                    //  mc_match=ttbarconstituents(ihyp);
                    //}
                    //if( !isData && (mc_match != 1) ) continue;
                    if ( !passSUSYTrigger2011_v1( isData , hyp_type()[ihyp] , true ) ) continue;
                    //if( !hypsFromFirstGoodDAvertx(ihyp, 1.0) ) continue;
                    int type = hyp_type()[ihyp];
                    VofP4 vjpts_p4;
                    LorentzVector lt_p4  = hyp_lt_p4()[ihyp];
                    LorentzVector ll_p4  = hyp_ll_p4()[ihyp];
                    int id_lt = hyp_lt_id()[ihyp];
                    int id_ll = hyp_ll_id()[ihyp];

                    // lepton energy scale systematic variation
                    if(fabs(id_lt) == 11) lt_p4*=leptonEnergyScaleFactor;
                    if(fabs(id_ll) == 11) ll_p4*=leptonEnergyScaleFactor;

                    int idx_lt = hyp_lt_index()[ihyp];
                    int idx_ll = hyp_ll_index()[ihyp];
                    // opposite charge
                    if (useOS)
                        if (id_lt * id_ll > 0)
                            continue;

                    if (verbose && useOS)
                        cout << "Hyp passes OS cuts" << endl;

                    // same sign
                    if (useSS)
                        if (id_lt * id_ll < 0)
                            continue;

                    if (verbose && useOS)
                        cout << "Hyp passes OS cuts" << endl;
                    /*

                    //if a muon, always require global and tracker
                    if(abs(id_lt)==13) {
                      if (((mus_type()[idx_lt]) & (1<<1)) == 0)    continue; // global muon
                      if (((mus_type()[idx_lt]) & (1<<2)) == 0)    continue; // tracker muon

                    }
                    if(abs(id_ll)==13) {
                      if (((mus_type()[idx_ll]) & (1<<1)) == 0)    continue; // global muon
                      if (((mus_type()[idx_ll]) & (1<<2)) == 0)    continue; // tracker muon
                    }
                    */

                    if (requireEcalEls)
                    {
                        //ask that the electron is ecal driven
                        if (abs(id_lt) == 11)
                        {
                            if (!(els_type()[idx_lt] & (1 << 2)))
                                continue;
                        }
                        if (abs(id_ll) == 11)
                        {
                            if (!(els_type()[idx_ll] & (1 << 2)))
                                continue;
                        }
                    }

                    //regardless of jet bin
                    //cut at Pt > 20, 20
                    if (usePtGt2020)
                    {
                        if (lt_p4.Pt() < 20. || ll_p4.Pt() < 20.)
                            continue;
                        if (verbose)
                            cout << "Hyp Passes Pt > 20, 20 cuts" << endl;
                    }

                    //cut at tight Pt > 20, loose Pt > 10
                    if (usePtGt2010)
                    {
                        if (TMath::Max(lt_p4.Pt(), ll_p4.Pt()) < 20)
                            continue;
                        if (TMath::Min(lt_p4.Pt(), ll_p4.Pt()) < 10)
                            continue;
                    }


                    //only look at events where loose lepton Pt < 20, tight > 20
                    if (excludePtGt2020)
                    {
                        if (lt_p4.Pt() > 20. && ll_p4.Pt() > 20.)
                            continue;
                        if (lt_p4.Pt() < 10 || ll_p4.Pt() < 10.)
                            continue;
                    }


                    if (vetoHypMassLt10)
                    {
                        if (hyp_p4()[ihyp].mass() < 10.)
                            continue;
                        if (verbose)
                            cout << "Passed Hyp Mass > 10 cut" << endl;
                    }

                    if (vetoHypMassLt12)
                    {
                        // if(hyp_p4()[ihyp].mass() < 12.)
                        if (hyp_p4()[ihyp].mass() < 20.)
                            continue;
                        if (verbose)
                            cout << "Passed Hyp Mass > 12 cut" << endl;
                    }
                    float FRweight = 1 ;
                    if (doFRestimation)
                    {

                        //          float FRweight = getFRWeight(hypIdx, elFRversion, mufr, elfr);
                        FRweight = getFRWeight(ihyp, mufr, elfr, frmode, isData);

                        // if getFRWeight returns less then 1 it means the current hyp
                        // does not fulfill the FO selections.
                        if (FRweight < -1.)
                            continue;

                    }//
                    if (applylepIDCuts)
                    {
                        //muon ID
                        if (abs(hyp_ll_id()[ihyp]) == 13  && !( muonId(hyp_ll_index()[ihyp] , OSGeneric_v3 ) ) )   continue;
                        if (abs(hyp_lt_id()[ihyp]) == 13  && !( muonId(hyp_lt_index()[ihyp] , OSGeneric_v3 ) ) )   continue;

                        //OSV1
                        if (abs(hyp_ll_id()[ihyp]) == 11  && !( pass_electronSelection( hyp_ll_index()[ihyp] , electronSelection_el_OSV3  ))) continue;
                        if (abs(hyp_lt_id()[ihyp]) == 11  && !( pass_electronSelection( hyp_lt_index()[ihyp] , electronSelection_el_OSV3  ))) continue;

                    }

                    v_goodHyps.push_back(ihyp);
                    v_weights.push_back(FRweight);
                }
                //
                // perform hypothesis disambiguation
                //
                if (v_goodHyps.size() == 0)
                    continue;

                if (hypDisamb)
                {
                    //returns the index of the best hypothesis in the vector of hypotheses
                    unsigned int goodHyp = selectHypByHighestSumPt(v_goodHyps);
                    vector<unsigned int>::const_iterator goodHyp_it = find(v_goodHyps.begin(), v_goodHyps.end(), goodHyp);
                    if (goodHyp_it == v_goodHyps.end())
                    {
                        cout << "The weight index does not correspond to the index of the best hypothesis!!!! Something is wrong"
                             << "We will quit" << endl;
                        return;
                    }

                    //clear this vector and put in the goodHyp vector in here so we can then save some space
                    //and loop over this vector below. Useful when we're not using the hypDisambiguation
                    //get the index of the goodHyp in the vector of goodHyps
                    unsigned int goodHyp_idx = goodHyp_it - v_goodHyps.begin();
                    //get the weight of the corresponding goodHyp
                    float goodHyp_weight = v_weights[goodHyp_idx];
                    v_goodHyps.clear();
                    v_weights.clear();
                    v_goodHyps.push_back(goodHyp);
                    v_weights.push_back(goodHyp_weight);
                }//if(hypDisamb)

            }//!applyNoCuts


            //now loop over the good hypotheses. If we require hypothesis disambiguation,
            //we will only have one entry in the vector of good hypothesis
            bool hasGoodHyp = false;
            int goodHyps_size = v_goodHyps.size();
            if (applyNoCuts) goodHyps_size = 1;

            for (unsigned int i = 0; i < goodHyps_size; i++)
            {

                unsigned int hypIdx = 999;
                vector<LorentzVector>  v_goodJets_cand_p4;
                vector<unsigned int> v_goodJets;
                int nBtagJets = -1;
                LorentzVector lt_p4;
                LorentzVector ll_p4;
                pair<float, float> p_met; //met and met phi
                float thefirstJet_pt = -999;
                float thesecondJet_pt = -999;
                float theSumBtagJetPt = -999;
                int i_ltbjet = -1;
                int i_llbjet = -1;

                if (!applyNoCuts)
                {
                    hypIdx = v_goodHyps[i];
                    weight = weight * v_weights[i];
                    int type = hyp_type()[hypIdx];
                    if ( !isData )
                    {
                        /*
                        if(type == 0) {       //mm
                          weight = weight*0.90;
                        }
                        else if (type == 1 || type == 2){ //em or me
                          weight = weight*0.95;
                        }
                        */
                        double trigger_weight = triggerEff(hypIdx,scaleTrigSFup,scaleTrigSFdown);
                        weight = weight * trigger_weight;
                    }

                    lt_p4  = hyp_lt_p4()[hypIdx];
                    ll_p4  = hyp_ll_p4()[hypIdx];
                    int id_lt = hyp_lt_id()[hypIdx];
                    int id_ll = hyp_ll_id()[hypIdx];

                    // lepton energy scale systematic variation
                    if(fabs(id_lt) == 11) lt_p4*=leptonEnergyScaleFactor;
                    if(fabs(id_ll) == 11) ll_p4*=leptonEnergyScaleFactor;

                    int idx_lt = hyp_lt_index()[hypIdx];
                    int idx_ll = hyp_ll_index()[hypIdx];
                    // z mass window
                    if (vetoZmass)
                    {
                        if (type == 0 || type == 3)
                        {
                            if (inZmassWindow(hyp_p4()[hypIdx].mass()))
                                continue;
                        }
                        if (verbose)
                            cout << "Hyp is not in Zmass window" << endl;
                    }//vetoZmass

                    if (requireZmass)
                    {
                        if (type == 0 || type == 3)
                        {
                            if (!inZmassWindow(hyp_p4()[hypIdx].mass()))
                                continue;
                        }
                    }//requireZmass

                    //  if(requireZmass && getVtxDistOnly) {
                    //           if(type == 1 || type == 2) continue;
                    //    int ndavtx = 0;

                    //    for (size_t v = 0; v < cms2.davtxs_position().size(); ++v){
                    //      if(isGoodDAVertex(v)) ++ndavtx;
                    //    }
                    //    //hnVtx[myType]                     ->Fill(ndavtx,               1);
                    //    //hnVtx[3]                          ->Fill(ndavtx,               1);
                    //    //      if (cms2.davtxs_position().size() >0 ) nEvents++ ;
                    //    // continue;
                    //         }//fill vtx using Z events

                    //get the jets passing cuts
                    //vector<unsigned int> v_goodJets;
                    vector<unsigned int> v_goodJetsNoEtaCut;
                    vector<LorentzVector> v_jetP4s;
                    float dmetx  = 0.0;
                    float dmety  = 0.0;
                    float jetptx = 0.0;
                    float jetpty = 0.0;


                    if (usecaloJets)
                    {
                        for (unsigned int i = 0; i < jets_p4().size(); i++)
                            v_jetP4s.push_back(jets_p4()[i]*jets_cor()[i]*globalJESRescale); //jets are uncorrected in our ntuples
                    }
                    if (usejptJets)
                    {
                        for (unsigned int i = 0; i < jpts_p4().size(); i++)
                            v_jetP4s.push_back(jpts_p4()[i] * jpts_corL1FastL2L3()[i]*globalJESRescale);
                    }
                    if (usepfJets)
                    {
                        // for (unsigned int i = 0; i < jpts_p4().size(); i++) //this is a bug
                        for (unsigned int i = 0; i < pfjets_p4().size(); i++)
                        {
                            LorentzVector vjet   = pfjets_p4()[i] * pfjets_corL1FastL2L3()[i];
                            pfUncertainty->setJetEta(vjet.eta());
                            pfUncertainty->setJetPt(vjet.pt());   // here you must use the CORRECTED jet pt
                            double unc = pfUncertainty->getUncertainty(true);
                            if ( scaleJESMETUp && isData )globalJESRescale = 1 + unc;
                            if ( scaleJESMETDown && isData )globalJESRescale = 1 - unc;
                            if ( scaleJER && !isData ) {
                                TRandom3 r3(  evt_event()*1000 + i  );
                                double JER_smear_width = sqrt( pow( JERsf(vjet.eta()) , 2 ) - 1. ) * unc;
                                double JER_smear = r3.Gaus(0,JER_smear_width);
                                globalJESRescale = 1. + JER_smear/vjet.pt();
                                if(globalJESRescale < 0.1 ) globalJESRescale = 0.1;
                                if(globalJESRescale > 10. ) globalJESRescale = 10.;
                                //cout<<JER_smear_width<<" "<<JER_smear<<" "<<globalJESRescale<<endl;
                            }
                            v_jetP4s.push_back(pfjets_p4()[i] * pfjets_corL1FastL2L3()[i]*globalJESRescale);

                        }
                    }

                    for (unsigned int j = 0; j < v_jetP4s.size(); ++j)
                    {
                        if (usejptJets && !passesCaloJetID(v_jetP4s.at(j)))
                            continue;
                        if (usepfJets && !passesPFJetID(j))
                            continue;
                        LorentzVector vlt  = hyp_lt_p4()[hypIdx];
                        LorentzVector vll  = hyp_ll_p4()[hypIdx];

                        if ( generalLeptonVeto )
                        {
                            bool rejectJet = false;
                            for ( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ )
                            {
                                if ( ROOT::Math::VectorUtil::DeltaR(v_jetP4s.at(j)  , goodLeptons.at(ilep) ) < 0.4 ) rejectJet = true;
                            }
                            if ( rejectJet ) continue;
                        }

                        if ( ROOT::Math::VectorUtil::DeltaR(v_jetP4s.at(j), vll) < 0.4 )  continue;
                        if ( ROOT::Math::VectorUtil::DeltaR(v_jetP4s.at(j), vlt) < 0.4 )  continue;
                        if ( v_jetP4s.at(j).pt() < JETPTCUT)              continue;
                        v_goodJetsNoEtaCut.push_back(j);
                        if ( fabs(v_jetP4s.at(j).eta() ) > 2.5 )            continue;
                        v_goodJets.push_back(j);

                        if (  v_jetP4s.at(j).pt()  > 10 )
                        {
                            dmetx  += v_jetP4s.at(j).px() - (pfjets_p4()[j] * pfjets_corL1FastL2L3()[j]).px()  ;
                            dmety  += v_jetP4s.at(j).py() - (pfjets_p4()[j] * pfjets_corL1FastL2L3()[j]).py()  ;
                            jetptx += (pfjets_p4()[j] * pfjets_corL1FastL2L3()[j]).px() ;
                            jetpty += (pfjets_p4()[j] * pfjets_corL1FastL2L3()[j]).py() ;
                        }
                    }



                    //if we want to veto on nJets, do it here
                    if (veto1Jet)
                    {
                        if (v_goodJets.size() < 1)
                            continue;
                        if (verbose)
                            cout << "Event passes 1 Jet cut" << endl;
                    }

                    if (veto2Jets)
                    {
                        if (v_goodJets.size() < 2)
                            continue;
                        if (verbose)
                            cout << "Event passes 2 Jet cut" << endl;
                    }


                    ///b-tagging

                    nBtagJets = 0;
                    vector<unsigned int> v_goodBtagJets;
                    vector<unsigned int> v_goodNonBtagJets;
                    TString btag_algo, btag_algo_name ;
                    if (BTagAlgTCHE)
                    {
                        btag_algo = "trackCountingHighEffBJetTag";
                        btag_algo_name = "TCHEM" ;
                    }
                    //else btag_algo = "simpleSecondaryVertexHighEffBJetTag";
                    else
                    {
                        btag_algo = "combinedSecondaryVertexBJetTag";
                        btag_algo_name = "CSVM" ;
                    }

                    for (unsigned int i = 0; i < v_goodJets.size(); i++)
                    {

                        if (usecaloJets)
                        {
                            if (passbTagging(v_goodJets[i],  "caloJets", btag_algo.Data()) )
                            {
                                nBtagJets++;
                                v_goodBtagJets.push_back(v_goodJets[i]);
                            }
                            else v_goodNonBtagJets.push_back(v_goodJets[i]);
                        }
                        else if (usejptJets)
                        {
                            if (passbTagging(v_goodJets[i],  "jptJets", btag_algo.Data()) )
                            {
                                nBtagJets++;
                                v_goodBtagJets.push_back(v_goodJets[i]);
                            }
                            else v_goodNonBtagJets.push_back(v_goodJets[i]);
                        }
                        else if (usepfJets)
                        {
                            if (passbTagging(v_goodJets[i],  "pfJets", btag_algo.Data()) )
                            {
                                nBtagJets++;
                                v_goodBtagJets.push_back(v_goodJets[i]);
                            }
                            else v_goodNonBtagJets.push_back(v_goodJets[i]);
                        }
                    }


                    //b tagging weights now taken care of later
                    /*
                        if( !isData ){
                          if(nBtagJets <2) {       //leave these unweighted for now (needs to be fixed)
                            weight = weight*1.0;
                          }
                          else if (nBtagJets>=2){
                            weight = weight*0.95*0.95; //0.95 is the MC-data scale factor for bjet tagging
                          }
                        }
                    */

                    if (requireBTag && nBtagJets < 1 )
                        continue;

                    if (requireExact2BTag && nBtagJets != 2)
                        continue;

                    if (require2BTag && nBtagJets < 2)
                        continue;


                    // MET cut
                    string metAlgo;
                    if (useCorMET   ) metAlgo  = "CorMET";
                    if (usetcMET    ) metAlgo  = "tcMET";
                    if (usetcMET35X ) metAlgo  = "tcMET35X";
                    if (usepfMET    ) metAlgo  = "pfMET";
                    //pair<float, float> p_met; //met and met phi
                    if (usetcMET || usepfMET || usetcMET35X)
                    {
                        p_met = getMet(metAlgo, hypIdx);

                        if (p_met.first < 0)
                        {
                            cout << "Something is wrong with the MET. Exiting" << endl;
                            return;
                        }

                        //if we are scaling the MET
                        if (isData && ( scaleJESMETDown || scaleJESMETUp ) )
                        {
                            float met = p_met.first;
                            float metPhi = p_met.second;
                            float metx = p_met.first * cos(metPhi);
                            float mety = p_met.first * sin(metPhi);

                            float lepx = hyp_p4()[hypIdx].Px();
                            float lepy = hyp_p4()[hypIdx].Py();

                            //--------------------------------------------------------
                            // calculate unclustered energy x and y components
                            // unclustered energy = -1 X ( MET + jets + leptons )
                            //--------------------------------------------------------

                            float unclustered_x = -1 * ( metx + jetptx );
                            float unclustered_y = -1 * ( mety + jetpty );

                            for ( int ilep = 0 ; ilep < goodLeptons.size() ; ilep++ )
                            {
                                unclustered_x -= goodLeptons.at(ilep).px();
                                unclustered_y -= goodLeptons.at(ilep).py();
                            }

                            //------------------------------------------------------------------------------
                            // now vary jets according to JEC uncertainty, vary unclustered energy by 10%
                            //------------------------------------------------------------------------------
                            float metNewx;
                            float metNewy;
                            if ( scaleJESMETUp && isData )
                            {
                                metNewx = metx - dmetx - 0.1 * unclustered_x;
                                metNewy = mety - dmety - 0.1 * unclustered_y;
                            }
                            if ( scaleJESMETDown && isData )
                            {
                                metNewx = metx + dmetx + 0.1 * unclustered_x;
                                metNewy = mety + dmety + 0.1 * unclustered_y;
                            }
                            /*
                            //hadronic component of MET (well, mostly), scaled
                            float metHx = (metx + lepx)*globalJESRescale;
                            float metHy = (mety + lepy)*globalJESRescale;
                            float metNewx = metHx - lepx;
                            float metNewy = metHy - lepy;
                            */
                            float metNewPhi = atan2(metNewy, metNewx);
                            p_met = make_pair(sqrt(metNewx * metNewx + metNewy * metNewy), metNewPhi);
                        }
                        //the cut is 30 for ee/mm hyps to reject DY
                        //20 for emu
                        if (vetoMET)
                        {
                            if (hyp_type()[hypIdx] == 0 || hyp_type()[hypIdx] == 3)
                            {
                                //if(p_met.first < 30.) continue;
                                if (p_met.first < 40.) continue;
                            }
                            /*
                            if(hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) {
                              if(p_met.first < 30.)   continue;

                            }
                            */
                            if (verbose)
                                cout << "Event passes MET cut" << endl;
                        }

                        //the cut is 50 for all hypotheses
                        if (vetoMET50)
                        {
                            if (p_met.first < 50.) continue;
                            if (verbose)
                                cout << "Event passes MET50 cut" << endl;
                        }

                        if (vetoProjectedMET)
                        {
                            if (hyp_type()[hypIdx] == 0 || hyp_type()[hypIdx] == 3)
                            {
                                if (p_met.first < 30.) continue;
                                // if(p_met.first < 50.) continue;
                            }
                            if (hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2)
                            {
                                if (p_met.first < 20.)   continue;
                                // if(p_met.first < 50.) continue;
                            }
                            if (v_goodJets.size() < 2 && hyp_p4()[hypIdx].M() < 80.)
                            {
                                if (projectedMET(p_met.first, p_met.second, hypIdx) < 10)
                                    continue;
                            }
                        }

                    }
                    else if (useCorMET)
                    {

                        cout << "THIS HAS BEEN COMMENTED OUT FOR NOW. DOES NOTHING" << endl;
                    }//if vetoing on corrected caloMet



                    //make a vector of corrected jets
                    vector<LorentzVector> v_goodJets_p4;
                    vector<LorentzVector> v_goodBtagJets_p4;
                    vector<LorentzVector> v_goodNonBtagJets_p4;

                    for (unsigned int i = 0; i < v_goodJets.size(); i++)
                    {

                        if (usecaloJets)
                            v_goodJets_p4.push_back(jets_p4()[v_goodJets[i]]*jets_cor()[v_goodJets[i]]*globalJESRescale);

                        else if (usejptJets)
                            v_goodJets_p4.push_back(jpts_p4()[v_goodJets[i]]*jpts_corL1FastL2L3()[v_goodJets[i]]*globalJESRescale);

                        else if (usepfJets)
                        {
                            LorentzVector vjet   = pfjets_p4()[v_goodJets[i]] * pfjets_corL1FastL2L3()[v_goodJets[i]];
                            pfUncertainty->setJetEta(vjet.eta());
                            pfUncertainty->setJetPt(vjet.pt());   // here you must use the CORRECTED jet pt
                            double unc = pfUncertainty->getUncertainty(true);
                            if ( scaleJESMETUp && isData )globalJESRescale = 1 + unc;
                            if ( scaleJESMETDown && isData )globalJESRescale = 1 - unc;
                            if ( scaleJER && !isData ) {
                                TRandom3 r3(  evt_event()*1000 + v_goodJets[i]  );
                                double JER_smear_width = sqrt( pow( JERsf(vjet.eta()) , 2 ) - 1. ) * unc;
                                double JER_smear = r3.Gaus(0,JER_smear_width);
                                globalJESRescale = 1. + JER_smear/vjet.pt();
                                if(globalJESRescale < 0.1 ) globalJESRescale = 0.1;
                                if(globalJESRescale > 10. ) globalJESRescale = 10.;
                                //cout<<JER_smear_width<<" "<<JER_smear<<" "<<globalJESRescale<<endl;
                            }
                            v_goodJets_p4.push_back(pfjets_p4()[v_goodJets[i]]*pfjets_corL1FastL2L3()[v_goodJets[i]]*globalJESRescale);
                        }
                    }

                    float theSumJetPt = 0.;
                    for (unsigned int ijet = 0; ijet < v_goodJets_p4.size(); ijet++)
                    {
                        theSumJetPt += v_goodJets_p4.at(ijet).Pt();
                    }

                    //make a vector of corrected Btag jets


                    for (unsigned int i = 0; i < v_goodBtagJets.size(); i++)
                    {

                        if (usecaloJets)
                            v_goodBtagJets_p4.push_back(jets_p4()[v_goodBtagJets[i]]*jets_cor()[v_goodBtagJets[i]]*globalJESRescale);

                        else if (usejptJets)
                            v_goodBtagJets_p4.push_back(jpts_p4()[v_goodBtagJets[i]]*jpts_corL1FastL2L3()[v_goodBtagJets[i]]*globalJESRescale);

                        else if (usepfJets)
                        {
                            LorentzVector vjet   = pfjets_p4()[v_goodBtagJets[i]] * pfjets_corL1FastL2L3()[v_goodBtagJets[i]];
                            pfUncertainty->setJetEta(vjet.eta());
                            pfUncertainty->setJetPt(vjet.pt());   // here you must use the CORRECTED jet pt
                            double unc = pfUncertainty->getUncertainty(true);
                            if ( scaleJESMETUp && isData )globalJESRescale = 1 + unc;
                            if ( scaleJESMETDown && isData )globalJESRescale = 1 - unc;
                            if ( scaleJER && !isData ) {
                                TRandom3 r3(  evt_event()*1000 + v_goodBtagJets[i]  );
                                double JER_smear_width = sqrt( pow( JERsf(vjet.eta()) , 2 ) - 1. ) * unc;
                                double JER_smear = r3.Gaus(0,JER_smear_width);
                                globalJESRescale = 1. + JER_smear/vjet.pt();
                                if(globalJESRescale < 0.1 ) globalJESRescale = 0.1;
                                if(globalJESRescale > 10. ) globalJESRescale = 10.;
                                //cout<<JER_smear_width<<" "<<JER_smear<<" "<<globalJESRescale<<endl;
                            }
                            v_goodBtagJets_p4.push_back(pfjets_p4()[v_goodBtagJets[i]]*pfjets_corL1FastL2L3()[v_goodBtagJets[i]]*globalJESRescale);
                        }
                    }


                    //make a vector of corrected Non-Btag jets


                    for (unsigned int i = 0; i < v_goodNonBtagJets.size(); i++)
                    {

                        if (usecaloJets)
                            v_goodNonBtagJets_p4.push_back(jets_p4()[v_goodNonBtagJets[i]]*jets_cor()[v_goodNonBtagJets[i]]*globalJESRescale);

                        else if (usejptJets)
                            v_goodNonBtagJets_p4.push_back(jpts_p4()[v_goodNonBtagJets[i]]*jpts_corL1FastL2L3()[v_goodNonBtagJets[i]]*globalJESRescale);

                        else if (usepfJets)
                        {
                            LorentzVector vjet   = pfjets_p4()[v_goodNonBtagJets[i]] * pfjets_corL1FastL2L3()[v_goodNonBtagJets[i]];
                            pfUncertainty->setJetEta(vjet.eta());
                            pfUncertainty->setJetPt(vjet.pt());   // here you must use the CORRECTED jet pt
                            double unc = pfUncertainty->getUncertainty(true);
                            if ( scaleJESMETUp && isData )globalJESRescale = 1 + unc;
                            if ( scaleJESMETDown && isData )globalJESRescale = 1 - unc;
                            if ( scaleJER && !isData ) {
                                TRandom3 r3(  evt_event()*1000 + v_goodNonBtagJets[i]  );
                                double JER_smear_width = sqrt( pow( JERsf(vjet.eta()) , 2 ) - 1. ) * unc;
                                double JER_smear = r3.Gaus(0,JER_smear_width);
                                globalJESRescale = 1. + JER_smear/vjet.pt();
                                if(globalJESRescale < 0.1 ) globalJESRescale = 0.1;
                                if(globalJESRescale > 10. ) globalJESRescale = 10.;
                                //cout<<JER_smear_width<<" "<<JER_smear<<" "<<globalJESRescale<<endl;
                            }
                            v_goodNonBtagJets_p4.push_back(pfjets_p4()[v_goodNonBtagJets[i]]*pfjets_corL1FastL2L3()[v_goodNonBtagJets[i]]*globalJESRescale);
                        }
                    }



                    theSumBtagJetPt = 0.;
                    for (unsigned int ijet = 0; ijet < v_goodBtagJets_p4.size(); ijet++)
                    {
                        theSumBtagJetPt += v_goodBtagJets_p4.at(ijet).Pt();
                    }

                    //for MT2
                    vector<LorentzVector> v_goodJetsNoEtaCut_p4;
                    for (unsigned int i = 0; i < v_goodJetsNoEtaCut.size(); i++)
                    {
                        if (usecaloJets)
                            v_goodJetsNoEtaCut_p4.push_back(jets_p4()[v_goodJetsNoEtaCut[i]]*jets_cor()[v_goodJetsNoEtaCut[i]]*globalJESRescale);
                        else if (usejptJets)
                        {
                            v_goodJetsNoEtaCut_p4.push_back(jpts_p4()[v_goodJetsNoEtaCut[i]]*jpts_corL1FastL2L3()[v_goodJetsNoEtaCut[i]]*globalJESRescale);
                        }
                        else if (usepfJets)
                        {
                            v_goodJetsNoEtaCut_p4.push_back(pfjets_p4()[v_goodJetsNoEtaCut[i]]*pfjets_corL1FastL2L3()[v_goodJetsNoEtaCut[i]]*globalJESRescale);
                        }
                    }

                    //apply HT CUT
                    if (applyHTCut)
                    {
                        if (theSumJetPt < 100) continue;
                    }


                    //sort the jets
                    std::sort(v_goodJets_p4.begin(), v_goodJets_p4.end(), sortByPt);
                    std::sort(v_goodBtagJets_p4.begin(), v_goodBtagJets_p4.end(), sortByPt);
                    std::sort(v_goodNonBtagJets_p4.begin(), v_goodNonBtagJets_p4.end(), sortByPt);


                    //apply b tagging weights
                    double weightb = 1.;
                    double weightbnottagged = 1.;
                    double weightmistag = 1.;
                    double weightnotmistagged = 1.;
                    double weightcmistag = 1.;
                    double weightcnotmistagged = 1.;
                    if (!isData)
                    {
                        vector<LorentzVector>  v_b_p4;
                        vector<LorentzVector>  v_c_p4;
                        for (unsigned int i = 0; i < genps_p4().size(); i++)
                        {
                            if ( fabs(genps_id()[i]) == 5) v_b_p4.push_back(genps_p4()[i]);
                            if ( fabs(genps_id()[i]) == 4) v_c_p4.push_back(genps_p4()[i]);
                        }
                        int nGenb = v_b_p4.size();
                        int nGenc = v_c_p4.size();

                        int i_matching_b[nBtagJets];
                        double dR_matching_b[nBtagJets];
                        for (unsigned int j = 0; j < nBtagJets; j++)
                        {
                            i_matching_b[j] = -1;
                            dR_matching_b[j] = 9999.;
                            i_matching_b[j] = match4vector(v_goodBtagJets_p4.at(j) , v_b_p4 , 0.4 );
                            if (i_matching_b[j] > -1) dR_matching_b[j] = ROOT::Math::VectorUtil::DeltaR( v_goodBtagJets_p4.at(j), v_b_p4.at(i_matching_b[j] ) );
                        }

                        vector<LorentzVector>  v_b_p4_unmatched;
                        int n_closest_match[nGenb];
                        for (unsigned int i = 0; i < nGenb; i++)
                        {
                            n_closest_match[i] = 0;
                            for (unsigned int j = 0; j < nBtagJets; j++)
                            {
                                if (i == i_matching_b[j]) n_closest_match[i]++;
                            }
                            if (n_closest_match[i] == 0) v_b_p4_unmatched.push_back(v_b_p4.at(i));
                        }
                        /*
                        //cout check
                        for(unsigned int i = 0; i < nGenb; i++) if(n_closest_match[i]>1) {
                            for(unsigned int j = 0; j< nGenb; j++) cout<<"n_closest_match "<<i<<" "<<j<<" "<<n_closest_match[j]<<endl;
                            for(unsigned int j = 0; j < nBtagJets; j++) cout<<"i_matching_b >1 closest "<<j<<" "<<i_matching_b[j]<<endl;
                        }
                        */
                        //if more than 1 b-tag jets have the same closest b, associate the closest matching pair
                        for (unsigned int j = 0; j < nBtagJets; j++)
                        {
                            for (unsigned int i = 0; i < nBtagJets; i++)
                            {
                                if (i != j && i_matching_b[i] != -1 && i_matching_b[i] == i_matching_b[j])
                                {
                                    //cout<<"more than 1 b-tag jets have the same closest b. nBtagJets = "<< nBtagJets<<" i_matching_b = "<<i_matching_b[i]<<endl;
                                    if (dR_matching_b[i] < dR_matching_b[j])
                                    {
                                        i_matching_b[j] =  match4vector(v_goodBtagJets_p4.at(j) , v_b_p4_unmatched , 0.4 );
                                        if (i_matching_b[j] > -1)
                                        {
                                            dR_matching_b[j] = ROOT::Math::VectorUtil::DeltaR( v_goodBtagJets_p4.at(j), v_b_p4_unmatched.at(i_matching_b[j] ) );
                                            i_matching_b[j] += 100; //so it can't have the same index as a b from v_b_p4
                                        }
                                    }
                                    else
                                    {
                                        i_matching_b[i] = match4vector(v_goodBtagJets_p4.at(i) , v_b_p4_unmatched , 0.4 );
                                        if (i_matching_b[i] > -1)
                                        {
                                            dR_matching_b[i] = ROOT::Math::VectorUtil::DeltaR( v_goodBtagJets_p4.at(i), v_b_p4_unmatched.at(i_matching_b[i] ) );
                                            i_matching_b[i] += 100; //so it can't have the same index as a b from v_b_p4
                                        }
                                    }

                                }
                            }
                        }

                        /*
                        //cout check again
                        for(unsigned int i = 0; i < nGenb; i++) if(n_closest_match[i]>1) {
                            for(unsigned int j = 0; j < nBtagJets; j++) cout<<"new i_matching_b >1 closest "<<j<<" "<<i_matching_b[j]<<endl;
                        }
                        */
                        //if more than 1 b-tag jets still have the same closest b, associate the closest matching pair
                        for (unsigned int j = 0; j < nBtagJets; j++)
                        {
                            for (unsigned int i = 0; i < nBtagJets; i++)
                            {
                                if (i != j && i_matching_b[i] != -1 && i_matching_b[i] == i_matching_b[j])
                                {
                                    cout << "more than 1 b-tag jets still have the same closest b. Only possible if nBtagJets>2. nBtagJets = " << nBtagJets << " i_matching_b = " << i_matching_b[i] << endl;
                                    if (dR_matching_b[i] < dR_matching_b[j]) i_matching_b[j] = -1;
                                    else i_matching_b[i] = -1;
                                }
                            }
                        }

                        //fill vector of b-tagged jets that aren't matched to a b. These are mistagged cs or light jets.
                        int nUnmatchedBtagJets = 0;
                        int j_unmatched[nBtagJets];
                        vector<LorentzVector>  v_goodBtagJets_unmatched_p4;
                        for (unsigned int j = 0; j < nBtagJets; j++)
                        {
                            j_unmatched[j] = -1;
                            if (i_matching_b[j] == -1)
                            {
                                v_goodBtagJets_unmatched_p4.push_back(v_goodBtagJets_p4.at(j));
                                j_unmatched[j] = nUnmatchedBtagJets;
                                nUnmatchedBtagJets++;
                            }
                        }





                        //now repeat all the matching, but with v_c_p4 instead of v_b_p4
                        int i_matching_c[nUnmatchedBtagJets];
                        if (nUnmatchedBtagJets > 0)
                        {
                            //if(nGenc>0) cout<<"trying to match to cs"<<endl;
                            double dR_matching_c[nUnmatchedBtagJets];
                            for (unsigned int j = 0; j < nUnmatchedBtagJets; j++)
                            {
                                i_matching_c[j] = -1;
                                dR_matching_c[j] = 9999.;
                                i_matching_c[j] = match4vector(v_goodBtagJets_unmatched_p4.at(j) , v_c_p4 , 0.4 );
                                if (i_matching_c[j] > -1) dR_matching_c[j] = ROOT::Math::VectorUtil::DeltaR( v_goodBtagJets_unmatched_p4.at(j), v_c_p4.at(i_matching_c[j] ) );
                            }

                            vector<LorentzVector>  v_c_p4_unmatched;
                            int n_closest_match_c[nGenc];
                            for (unsigned int i = 0; i < nGenc; i++)
                            {
                                n_closest_match_c[i] = 0;
                                for (unsigned int j = 0; j < nUnmatchedBtagJets; j++)
                                {
                                    if (i == i_matching_c[j]) n_closest_match_c[i]++;
                                }
                                if (n_closest_match_c[i] == 0) v_c_p4_unmatched.push_back(v_c_p4.at(i));
                            }
                            /*
                            //cout check
                            for(unsigned int i = 0; i < nGenc; i++) if(n_closest_match_c[i]>1) {
                                for(unsigned int j = 0; j< nGenc; j++) cout<<"n_closest_match_c "<<i<<" "<<j<<" "<<n_closest_match_c[j]<<endl;
                                for(unsigned int j = 0; j < nUnmatchedBtagJets; j++) cout<<"i_matching_c >1 closest "<<j<<" "<<i_matching_c[j]<<endl;
                            }
                            */
                            //if more than 1 b-tag jets have the same closest c, associate the closest matching pair
                            for (unsigned int j = 0; j < nUnmatchedBtagJets; j++)
                            {
                                for (unsigned int i = 0; i < nUnmatchedBtagJets; i++)
                                {
                                    if (i != j && i_matching_c[i] != -1 && i_matching_c[i] == i_matching_c[j])
                                    {
                                        //cout<<"more than 1 b-tag jets have the same closest c. nUnmatchedBtagJets = "<< nUnmatchedBtagJets<<" i_matching_c = "<<i_matching_c[i]<<endl;
                                        if (dR_matching_c[i] < dR_matching_c[j])
                                        {
                                            i_matching_c[j] =  match4vector(v_goodBtagJets_unmatched_p4.at(j) , v_c_p4_unmatched , 0.4 );
                                            if (i_matching_c[j] > -1)
                                            {
                                                dR_matching_c[j] = ROOT::Math::VectorUtil::DeltaR( v_goodBtagJets_unmatched_p4.at(j), v_c_p4_unmatched.at(i_matching_c[j] ) );
                                                i_matching_c[j] += 100; //so it can't have the same index as a c from v_c_p4
                                            }
                                        }
                                        else
                                        {
                                            i_matching_c[i] = match4vector(v_goodBtagJets_unmatched_p4.at(i) , v_c_p4_unmatched , 0.4 );
                                            if (i_matching_c[i] > -1)
                                            {
                                                dR_matching_c[i] = ROOT::Math::VectorUtil::DeltaR( v_goodBtagJets_unmatched_p4.at(i), v_c_p4_unmatched.at(i_matching_c[i] ) );
                                                i_matching_c[i] += 100; //so it can't have the same index as a c from v_c_p4
                                            }
                                        }

                                    }
                                }
                            }

                            /*
                            //cout check again
                            for(unsigned int i = 0; i < nGenc; i++) if(n_closest_match_c[i]>1) {
                                for(unsigned int j = 0; j < nUnmatchedBtagJets; j++) cout<<"new i_matching_c >1 closest "<<j<<" "<<i_matching_c[j]<<endl;
                            }
                            */
                            //if more than 1 b-tag jets still have the same closest c, associate the closest matching pair
                            for (unsigned int j = 0; j < nUnmatchedBtagJets; j++)
                            {
                                for (unsigned int i = 0; i < nUnmatchedBtagJets; i++)
                                {
                                    if (i != j && i_matching_c[i] != -1 && i_matching_c[i] == i_matching_c[j])
                                    {
                                        cout << "more than 1 b-tag jets still have the same closest c. Only possible if nUnmatchedBtagJets>2. nUnmatchedBtagJets = " << nUnmatchedBtagJets << " i_matching_c = " << i_matching_c[i] << endl;
                                        if (dR_matching_c[i] < dR_matching_c[j]) i_matching_c[j] = -1;
                                        else i_matching_c[i] = -1;
                                    }
                                }
                            }
                        }











                        double btageffdata = getBTagEff("CSV", 0.679, false);
                        double ctageffdata = getCTagEff("CSV", 0.679, false);
                        double  nonb_pt, nonb_eta, bjet_pt, bjet_eta;

                        //b tagged jet weighting, using number of matching real bs, cs and mistags
                        int nBtagJets_real = 0;
                        int nMistags = 0;
                        int ncMistags = 0;
                        double btag_sf = 1.0;
                        double ctag_sf = 1.0;

                        if (BTagAlgTCHE)
                        {
                            btag_sf = getBTagSF("TCHE", 3.3);
                            // ctag_sf equal to btag_sf but error is different
                            // assumption true for CSV, implementing it here for code safety, needs to be checked if TCHE is used
                            ctag_sf = btag_sf;
                            // systematics, add or subtract one sigma error
                            if ( scaleBTAGSFup ) {
                              btag_sf += getBTagSF_Err("TCHE");
                            } else if ( scaleBTAGSFdown ) {
                              btag_sf -= getBTagSF_Err("TCHE");
                            }
                            // systematics, add or subsctract one sigma error
                            // at the same time when btag sf is scaled
                            if ( scaleBTAGSFup ) {
                              ctag_sf += getCTagSF_Err("TCHE");
                            } else if ( scaleBTAGSFdown ) {
                              ctag_sf -= getCTagSF_Err("TCHE");
                            }
                          }
                        else
                        {
                            btag_sf = getBTagSF("CSV", 0.679);
                            // ctag_sf equal to btag_sf but error is different
                            ctag_sf = btag_sf;
                            // systematics, add or subtract one sigma error
                            if ( scaleBTAGSFup ) {
                              btag_sf += getBTagSF_Err("CSV");
                            } else if ( scaleBTAGSFdown ) {
                              btag_sf -= getBTagSF_Err("CSV");
                            }
                            // systematics, add or subsctract one sigma error
                            // at the same time when btag sf is scaled
                            if ( scaleBTAGSFup ) {
                              ctag_sf += getCTagSF_Err("CSV");
                            } else if ( scaleBTAGSFdown ) {
                              ctag_sf -= getCTagSF_Err("CSV");
                            }
                        }

                        for (unsigned int j = 0; j < nBtagJets; j++)
                        {
                            bjet_pt =  v_goodBtagJets_p4.at(j).Pt();
                            bjet_eta =  fabs(v_goodBtagJets_p4.at(j).Eta());
                            //to keep within the range of the function
                            if (bjet_pt > 499.0)
                            {
                                bjet_pt =  499.0;
                            }
                            if (bjet_eta > 2.3)
                            {
                                bjet_eta =  2.3;
                            }
                            if (i_matching_b[j] > -1)
                            {
                                nBtagJets_real++;
                                weightb *= btag_sf;
                                if ( j_unmatched[j] != -1) cout << "something went wrong filling v_goodBtagJets_unmatched_p4" << endl;
                            }
                            else if (i_matching_c[j_unmatched[j]] > -1)
                            {
                                //cout<<"nUnmatched, j, j_unmatched[j], i_matching_c "<<nUnmatchedBtagJets<<" "<<j<<" "<<j_unmatched[j]<<" "<<i_matching_c[j_unmatched[j]]<<endl;
                                ncMistags++;
                                weightcmistag *= ctag_sf;
                                if ( j_unmatched[j] == -1) cout << "something went wrong filling v_goodBtagJets_unmatched_p4" << endl;
                            }
                            else
                            {
                                nMistags++;
                                weightmistag *= getMisTagSF( bjet_pt, bjet_eta  , btag_algo_name.Data());
                                if ( j_unmatched[j] == -1) cout << "something went wrong filling v_goodBtagJets_unmatched_p4" << endl;
                            }
                        }

                        //if(nGenc>0) cout<<"nGenb, nBtagJets, nBtagJets_real, nGenc, ncMistags, nMistags "<<nGenb<<" "<<nBtagJets<<" "<<nBtagJets_real<<" "<<nGenc<<" "<<ncMistags<<" "<<nMistags<<endl;
                        if (nBtagJets_real > nGenb) cout << "***problem with b matching in btag weighting***" << endl;
                        if (ncMistags > nGenc) cout << "***problem with c matching in btag weighting***" << endl;

                        //weighting for the non b-jets
                        int nNonBtagJets = v_goodNonBtagJets_p4.size();
                        int i_nonb_matching_b[nNonBtagJets];
                        double dR_nonb_matching_b[nNonBtagJets];
                        int i_nonb_matching_c[nNonBtagJets];
                        double dR_nonb_matching_c[nNonBtagJets];

                        for (unsigned int j = 0; j < nNonBtagJets; j++)
                        {
                            //match non b-jets to bs and cs.
                            i_nonb_matching_b[j] = -1;
                            dR_nonb_matching_b[j] = 9999.;
                            i_nonb_matching_c[j] = -1;
                            dR_nonb_matching_c[j] = 9999.;
                            i_nonb_matching_b[j] = match4vector( v_goodNonBtagJets_p4.at(j) , v_b_p4 , 0.4 );
                            i_nonb_matching_c[j] = match4vector( v_goodNonBtagJets_p4.at(j) , v_c_p4 , 0.4 );
                            if (i_nonb_matching_b[j] > -1) dR_nonb_matching_b[j] = ROOT::Math::VectorUtil::DeltaR( v_goodNonBtagJets_p4.at(j), v_b_p4.at(i_nonb_matching_b[j] ) );
                            if (i_nonb_matching_c[j] > -1) dR_nonb_matching_c[j] = ROOT::Math::VectorUtil::DeltaR( v_goodNonBtagJets_p4.at(j), v_c_p4.at(i_nonb_matching_c[j] ) );
                        }

                        //if more than 1 non b-tag jets have the same closest b or c, associate the closest matching pair
                        for (unsigned int j = 0; j < nNonBtagJets; j++)
                        {
                            for (unsigned int i = 0; i < nNonBtagJets; i++)
                            {
                                if (i != j && i_nonb_matching_b[i] != -1 && i_nonb_matching_b[i] == i_nonb_matching_b[j])
                                {
                                    //cout<<"more than 1 non b-tag jets have the same closest b.  nNonBtagJets = "<< nNonBtagJets<<" i_nonb_matching_b = "<<i_nonb_matching_b[i]<<endl;
                                    if (dR_nonb_matching_b[i] < dR_nonb_matching_b[j]) i_nonb_matching_b[j] = -1;
                                    else i_nonb_matching_b[i] = -1;
                                }
                                if (i != j && i_nonb_matching_c[i] != -1 && i_nonb_matching_c[i] == i_nonb_matching_c[j])
                                {
                                    //cout<<"more than 1 non b-tag jets have the same closest c.  nNonBtagJets = "<< nNonBtagJets<<" i_nonb_matching_c = "<<i_nonb_matching_c[i]<<endl;
                                    if (dR_nonb_matching_c[i] < dR_nonb_matching_c[j]) i_nonb_matching_c[j] = -1;
                                    else i_nonb_matching_c[i] = -1;
                                }
                            }
                        }

                        for (unsigned int i = 0; i < v_goodNonBtagJets_p4.size(); i++)
                        {
                            //weight the untagged jets based on whether they match a b, a c, or neither.
                            nonb_pt =  v_goodNonBtagJets_p4.at(i).Pt();
                            nonb_eta =  fabs(v_goodNonBtagJets_p4.at(i).Eta());
                            //to keep within the range of the function
                            if (nonb_pt > 499.0)
                            {
                                nonb_pt =  499.0;
                            }
                            if (nonb_eta > 2.3)
                            {
                                nonb_eta =  2.3;
                            }
                            double bjetfrdata =  getMisTagRate( nonb_pt, nonb_eta  , btag_algo_name.Data());

                            //if the matched real b or c is already matched to a btag jet, the matching to the nonb is overridden
                            bool alreadymatched = false;
                            bool alreadymatchedc = false;
                            for (unsigned int j = 0; j < nBtagJets; j++)
                            {
                                if (i_nonb_matching_b[i] != -1 && i_nonb_matching_b[i] == i_matching_b[j]) alreadymatched = true;
                                if (j_unmatched[j] > -1) if (i_nonb_matching_c[i] != -1 && i_nonb_matching_c[i] == i_matching_c[j_unmatched[j]]) alreadymatchedc = true;
                            }
                            //if(alreadymatchedc) cout<<"alreadymatchedc"<<endl;
                            if ( !(i_nonb_matching_b[i] == -1 || alreadymatched) ) weightbnottagged *= (1. - btageffdata) / (1. - btageffdata / btag_sf );
                            if ( (i_nonb_matching_b[i] == -1 || alreadymatched) && !(i_nonb_matching_c[i] == -1 || alreadymatchedc) ) weightcnotmistagged *= (1. - ctageffdata) / (1. - ctageffdata / ctag_sf );
                            if ( (i_nonb_matching_b[i] == -1 || alreadymatched) && (i_nonb_matching_c[i] == -1 || alreadymatchedc) ) weightnotmistagged *= (1. - bjetfrdata) / (1. - bjetfrdata / getMisTagSF( nonb_pt, nonb_eta  , btag_algo_name.Data()));
                        }
                        //if(nGenc>0) cout<<"weightb, weightbnottagged, weightcmistag, weightcnotmistagged, weightmistag, weightnotmistagged: "<<weightb<<" , "<<weightbnottagged <<" , "<<weightcmistag<<" , "<<weightcnotmistagged <<" , "<<weightmistag<<" , "<<weightnotmistagged <<endl;
                        weight = weight * weightb * weightbnottagged * weightcmistag * weightcnotmistagged * weightmistag * weightnotmistagged;
                        if (doBFR && nBtagJets == 2 && nMistags == 0) continue;
                    }




                    float BFRweight = 1 ;
                    if (doBFR)
                    {
                        BFRweight = getBFRWeight(hypIdx, v_goodNonBtagJets_p4, v_goodBtagJets_p4, isData);

                        // if getFRWeight returns less then 1 it means the current hyp
                        // does not fulfill the FO selections.
                        if (BFRweight < -1.)
                            continue;
                        cout << "weight =  " << BFRweight << endl;
                        weight =  weight * BFRweight;
                    }



                    //make the list of jets which will be comined with leptons
                    //vector<LorentzVector>  v_goodJets_cand_p4;
                    vector<LorentzVector> v_goodJets_p4_tmp;

                    v_goodJets_p4_tmp = v_goodJets_p4;

                    if ( nBtagJets == 0 && veto2Jets)
                    {
                        if (sortJetCandidatesbyPt)
                        {
                            v_goodJets_cand_p4.push_back(v_goodJets_p4[0]);
                            v_goodJets_cand_p4.push_back(v_goodJets_p4[1]);
                        }
                        else if (sortJetCandidatesbyDR)
                            v_goodJets_cand_p4 = v_goodJets_p4;
                    }
                    else if (nBtagJets == 1)
                    {
                        int overlap_jet_idx = -1;
                        overlap_jet_idx = match4vector(v_goodBtagJets_p4[0], v_goodJets_p4);
                        v_goodJets_p4_tmp.erase(v_goodJets_p4_tmp.begin() + overlap_jet_idx);
                        if (sortJetCandidatesbyPt)
                        {
                            v_goodJets_cand_p4.push_back(v_goodBtagJets_p4[0]);
                            v_goodJets_cand_p4.push_back(v_goodJets_p4_tmp[0]);
                        }
                        else if (sortJetCandidatesbyDR)
                        {
                            v_goodJets_cand_p4.push_back(v_goodBtagJets_p4[0]);
                            for (unsigned int i = 0; i < v_goodJets_p4_tmp.size(); i++)
                            {
                                v_goodJets_cand_p4.push_back(v_goodJets_p4_tmp[i]);
                            }
                        }
                    }
                    else if (nBtagJets >= 2)
                    {
                        if (sortJetCandidatesbyPt)
                        {
                            v_goodJets_cand_p4.push_back(v_goodBtagJets_p4[0]);
                            v_goodJets_cand_p4.push_back(v_goodBtagJets_p4[1]);
                        }
                        else if (sortJetCandidatesbyDR)
                            v_goodJets_cand_p4 = v_goodBtagJets_p4;
                    }

                    npreSelectedEvents++;

                    i_ltbjet = -1;
                    i_llbjet = -1;
                    //  vector<LorentzVector>  v_goodJets_cand_p4_tmp;
                    //  v_goodJets_cand_p4_tmp = v_goodJets_cand_p4;
                    double dr_ltb, dr_llb;
                    //double mass_ltb, mass_llb;
                    double mass_ltbs, mass_llbs;
                    thefirstJet_pt = -999;
                    thesecondJet_pt = -999;
                    double m_1, m_2, m_3, m_4;
                    double mass_lb_min, mass_lb_min_otherpair;

                    if (v_goodJets_cand_p4.size() > 1)
                    {


                        i_ltbjet = match4vector(lt_p4, v_goodJets_cand_p4);
                        i_llbjet = match4vector(ll_p4, v_goodJets_cand_p4);
                        if (matchLeptonJetbyMaxDR)
                        {
                            i_ltbjet = antimatch4vector(lt_p4, v_goodJets_cand_p4);
                            i_llbjet = antimatch4vector(ll_p4, v_goodJets_cand_p4);
                        }

                        if (i_ltbjet == i_llbjet) // the closest jets are the same.
                        {
                            dr_ltb = ROOT::Math::VectorUtil::DeltaR(lt_p4, v_goodJets_cand_p4.at(i_ltbjet));
                            dr_llb = ROOT::Math::VectorUtil::DeltaR(ll_p4, v_goodJets_cand_p4.at(i_llbjet));
                            if (dr_ltb > dr_llb)i_ltbjet = antimatch4vector(lt_p4, v_goodJets_cand_p4);
                            else i_llbjet =  antimatch4vector(ll_p4, v_goodJets_cand_p4);
                            if (matchLeptonJetbyMaxDR)
                            {
                                if (dr_ltb < dr_llb)i_ltbjet = match4vector(lt_p4, v_goodJets_cand_p4);
                                else i_llbjet =  match4vector(ll_p4, v_goodJets_cand_p4);
                            }
                        }
                        mass_ltb = (lt_p4 + v_goodJets_cand_p4.at(i_ltbjet)).M();
                        mass_llb = (ll_p4 + v_goodJets_cand_p4.at(i_llbjet)).M();

                        // if(applyLeptonJetInvMassCut170 && !(mass_ltb > 170 && mass_llb> 170)) continue;

                        thefirstJet_pt = v_goodJets_cand_p4.at(i_ltbjet).Pt();
                        thesecondJet_pt = v_goodJets_cand_p4.at(i_llbjet).Pt();
                        massllbb = (lt_p4 + v_goodJets_cand_p4.at(i_ltbjet) + ll_p4 + v_goodJets_cand_p4.at(i_llbjet)).M();
                    }

                    //float m_top;
                    //vector <float> AMWTweight, AMWTweight_nojetsmear;
                    //vector <TLorentzVector> top1_p4, top2_p4, top1_nojetsmear_p4, top2_nojetsmear_p4;
                    //TLorentzVector cms, cms_nojetsmear, lepPlus,lepMinus, jet1,jet2;

                    //first solve with no jet smearing for comparison (to check bias caused by jet smearing)
                    m_top_nojetsmear = getTopMassEstimate(d_llsol, hypIdx, v_goodJets_cand_p4, p_met.first, p_met.second, 1, top1_nojetsmear_p4, top2_nojetsmear_p4, AMWTweight_nojetsmear, AMWTmass);
                    if ( m_top_nojetsmear > 0 && (fabs(m_top_nojetsmear - top1_nojetsmear_p4[0].M()) > 0.5 || fabs(m_top_nojetsmear - top2_nojetsmear_p4[0].M()) > 0.5) ) cout << "*** mass solution mismatch (no smearing) *** " << m_top_nojetsmear << " " << top1_nojetsmear_p4[0].M() << " " << top2_nojetsmear_p4[0].M() << endl;
                    tt_mass_nojetsmear = -999.0;
                    if (m_top_nojetsmear > 0) tt_mass_nojetsmear = (top1_nojetsmear_p4[0] + top2_nojetsmear_p4[0]).M();

                    //cout<<AMWTweight_nojetsmear.size()<<endl;

                    //now repeat using jet smearing
                    m_top = getTopMassEstimate(d_llsol, hypIdx, v_goodJets_cand_p4, p_met.first, p_met.second, 100, top1_p4, top2_p4, AMWTweight, AMWTmass);

                    Nsolns = AMWTweight.size();
                    imaxAMWTweight = -999;
                    maxAMWTweight = -999;
                    sumAMWTweight = 0.;
                    for (int ia = 0; ia < Nsolns; ++ia)
                    {
                        sumAMWTweight += AMWTweight[ia];
                        if (AMWTweight[ia] > maxAMWTweight)
                        {
                            maxAMWTweight = AMWTweight[ia];
                            imaxAMWTweight = ia;
                        }
                    }
                    //if(Nsolns<100 && Nsolns>0) cout<<"maxAMWTweight: "<<maxAMWTweight<<" i: "<<imaxAMWTweight<<" size: "<<Nsolns<<endl;
                    //cout << "got to line: " << __LINE__ <<endl;

                    if (Nsolns < 1) Nsolns = 1;
                    if (useOnlyMaxWsoln) Nsolns = 1;
                    weight = weight / double(Nsolns);
                    aveAMWTweight = sumAMWTweight / double(Nsolns);


                }//!applynocuts

                if (applyNoCuts) Nsolns = 1;


                
                if ( applyTopPtWeighting && (prefix == "ttdil" || prefix == "ttotr") )
                {

                    TLorentzVector topplus_genp_p4(0, 0, 0, 0), topminus_genp_p4(0, 0, 0, 0);
                    bool hastop = false;
                    bool hastbar = false;

                    for (unsigned int i = 0; i < genps_p4().size(); i++)
                    {
                        if (genps_status()[i] == 3)
                        {

                            if (genps_id()[i] == 6 )
                            {
                                topplus_genp_p4.SetXYZT( genps_p4()[i].x(),
                                                         genps_p4()[i].y(),
                                                         genps_p4()[i].z(),
                                                         genps_p4()[i].t()
                                                       );
                                hastop = true;
                            }
                            else if (genps_id()[i] == -6 )
                            {
                                topminus_genp_p4.SetXYZT( genps_p4()[i].x(),
                                                          genps_p4()[i].y(),
                                                          genps_p4()[i].z(),
                                                          genps_p4()[i].t()
                                                        );
                                hastbar = true;
                            }

                        }
                    }

                    float pT_topplus_gen = topplus_genp_p4.Pt();
                    float pT_topminus_gen = topminus_genp_p4.Pt();
                    if (hastop && hastbar) {
                      if ( useReweightingUncorrelated ) {
                        weight = weight * TopPtWeight(pT_topplus_gen) * TopPtWeight(pT_topminus_gen);
                      } else if ( useReweightingLeadingObject ) {
                        weight = weight * TopPtWeight(pT_topplus_gen);
                      } else {
                        weight = weight * sqrt( TopPtWeight(pT_topplus_gen) * TopPtWeight(pT_topminus_gen) );                       }
                    } 
                }


                if ( applyLeptonPtWeighting && (prefix == "ttdil" || prefix == "ttotr") ){
                  if ( hyp_lt_id()[hypIdx] < 0 )
                  {
                      lepPlus.SetXYZT(
                          hyp_lt_p4()[hypIdx].x(),
                          hyp_lt_p4()[hypIdx].y(),
                          hyp_lt_p4()[hypIdx].z(),
                          hyp_lt_p4()[hypIdx].t()
                      );

                      lepMinus.SetXYZT(
                          hyp_ll_p4()[hypIdx].x(),
                          hyp_ll_p4()[hypIdx].y(),
                          hyp_ll_p4()[hypIdx].z(),
                          hyp_ll_p4()[hypIdx].t()
                      );

                  }
                  else
                  {
                      lepPlus.SetXYZT(
                          hyp_ll_p4()[hypIdx].x(),
                          hyp_ll_p4()[hypIdx].y(),
                          hyp_ll_p4()[hypIdx].z(),
                          hyp_ll_p4()[hypIdx].t()
                      );

                      lepMinus.SetXYZT(
                          hyp_lt_p4()[hypIdx].x(),
                          hyp_lt_p4()[hypIdx].y(),
                          hyp_lt_p4()[hypIdx].z(),
                          hyp_lt_p4()[hypIdx].t()
                      );

                  }
                  float lepPlus_Pt = lepPlus.Pt();
                  float lepMinus_Pt = lepMinus.Pt();
                  if ( useReweightingUncorrelated ) {
                    weight = weight * LeptonPtWeight(lepPlus_Pt) * LeptonPtWeight(lepMinus_Pt);                            
                  } else if ( useReweightingLeadingObject ) {
                    weight = weight * LeptonPtWeight(lepPlus_Pt);
                  } else {
                    weight = weight * sqrt( LeptonPtWeight(lepPlus_Pt) * LeptonPtWeight(lepMinus_Pt) );                            
                  }
                }
                if ( applyJetPtWeighting && (prefix == "ttdil" || prefix == "ttotr") ){
                  if ( useReweightingUncorrelated ) {
                    weight = weight * JetPtWeight(thefirstJet_pt) * JetPtWeight(thesecondJet_pt);
                  } else if ( useReweightingLeadingObject ) {
                    weight = weight * JetPtWeight(thefirstJet_pt);
                  } else {
                    weight = weight * sqrt( JetPtWeight(thefirstJet_pt) * JetPtWeight(thesecondJet_pt) );
                  }
                }


                for (int i_smear = 0; i_smear < Nsolns; ++i_smear)
                {

                    if (!applyNoCuts)
                    {
                        if (useOnlyMaxWsoln) i_smear = ( imaxAMWTweight > 0 ? imaxAMWTweight : 0 );

                        if ( m_top > 0 && (fabs(m_top - top1_p4[i_smear].M()) > 0.5 || fabs(m_top - top2_p4[i_smear].M()) > 0.5) ) cout << "*** mass solution mismatch *** " << m_top << " " << top1_p4[i_smear].M() << " " << top2_p4[i_smear].M() << endl;
                        tt_mass = -999.0;
                        if (m_top > 0) tt_mass = (top1_p4[i_smear] + top2_p4[i_smear]).M();

                        //cout<<m_top<<" "<<m_top_nojetsmear<<" "<<tt_mass<<" "<<tt_mass_nojetsmear<<endl;

                        //float ttRapidity = top1_p4[i_smear].Eta()+top2_p4[i_smear].Eta();

                        ttRapidity = -999.0;
                        if (m_top > 0) ttRapidity = top1_p4[i_smear].Rapidity() + top2_p4[i_smear].Rapidity();
                        ttRapidity2 = -999.0;
                        if (m_top > 0) ttRapidity2 = (top1_p4[i_smear] + top2_p4[i_smear]).Rapidity();
                        //if(m_top < 0) continue;
                        if ((applyLeptonJetInvMassCut450 || applyTopSystEta ) && m_top < 0) continue;
                        if (applyLeptonJetInvMassCut450 && (tt_mass < 450 )) continue;
                        if (applyTopSystEta &&  (fabs(ttRapidity) < 2.0) ) continue;

                        top_rapiditydiff_cms = -999.0;
                        if (m_top > 0) top_rapiditydiff_cms = (top1_p4[i_smear].Rapidity() - top2_p4[i_smear].Rapidity()) * (top1_p4[i_smear].Rapidity() + top2_p4[i_smear].Rapidity());

                        top_pseudorapiditydiff_cms = -999.0;
                        if (m_top > 0) top_pseudorapiditydiff_cms = abs(top1_p4[i_smear].Eta()) - abs(top2_p4[i_smear].Eta());

                        top_rapiditydiff_Marco = -999.0;
                        if (m_top > 0) top_rapiditydiff_Marco = abs(top1_p4[i_smear].Rapidity()) - abs(top2_p4[i_smear].Rapidity());

                        float top_rapiditydiff2_cms = -999.0;
                        if (m_top > 0) top_rapiditydiff2_cms = (top1_p4[i_smear].Rapidity() - top2_p4[i_smear].Rapidity());

                        float top_pseudorapiditydiff2_cms = -999.0;
                        if (m_top > 0) top_pseudorapiditydiff2_cms = (top1_p4[i_smear].Eta()) - (top2_p4[i_smear].Eta());



                        //ttbar solution in lab frame for bias check
                        float top1_pt_nojetsmear = -999.0;
                        float top2_pt_nojetsmear = -999.0;
                        //float top1_eta_nojetsmear = -999.0;
                        //float top1_phi_nojetsmear = -999.0;
                        //float top2_eta_nojetsmear = -999.0;
                        //float top2_phi_nojetsmear = -999.0;
                        if (m_top_nojetsmear > 0)
                        {
                            top1_pt_nojetsmear =  top1_nojetsmear_p4[0].Pt();
                            top2_pt_nojetsmear =  top2_nojetsmear_p4[0].Pt();
                            //top1_eta_nojetsmear =  top1_nojetsmear_p4[0].Eta();
                            //top1_phi_nojetsmear =  top1_nojetsmear_p4[0].Phi();
                            //top2_eta_nojetsmear =  top2_nojetsmear_p4[0].Eta();
                            //top2_phi_nojetsmear =  top2_nojetsmear_p4[0].Phi();
                        }


                        if (m_top_nojetsmear > 0) cms_nojetsmear = top1_nojetsmear_p4[0] + top2_nojetsmear_p4[0];
                        tt_pT_nojetsmear = -999.0;
                        if (m_top_nojetsmear > 0) tt_pT_nojetsmear = cms_nojetsmear.Pt();
                        if (m_top_nojetsmear > 0) top1_nojetsmear_p4[0].Boost(-cms_nojetsmear.BoostVector());
                        if (m_top_nojetsmear > 0) top2_nojetsmear_p4[0].Boost(-cms_nojetsmear.BoostVector());


                        //ttbar solution in CM frame for bias check
                        float top1_p_CM_nojetsmear = -999.0;
                        float top2_p_CM_nojetsmear = -999.0;
                        if (m_top_nojetsmear > 0)
                        {
                            top1_p_CM_nojetsmear =  top1_nojetsmear_p4[0].P();
                            top2_p_CM_nojetsmear =  top2_nojetsmear_p4[0].P();
                        }



                        //ttbar solution in lab frame for bias check
                        float top1_pt = -999.0;
                        float top2_pt = -999.0;
                        //float top1_eta = -999.0;
                        //float top1_phi = -999.0;
                        //float top2_eta = -999.0;
                        //float top2_phi = -999.0;
                        if (m_top > 0)
                        {
                            top1_pt =  top1_p4[i_smear].Pt();
                            top2_pt =  top2_p4[i_smear].Pt();
                            //top1_eta =  top1_p4[i_smear].Eta();
                            //top1_phi =  top1_p4[i_smear].Phi();
                            //top2_eta =  top2_p4[i_smear].Eta();
                            //top2_phi =  top2_p4[i_smear].Phi();
                        }


                        if (m_top > 0) cms = top1_p4[i_smear] + top2_p4[i_smear];
                        tt_pT = -999.0;
                        if (m_top > 0) tt_pT = cms.Pt();
                        if (m_top > 0) top1_p4[i_smear].Boost(-cms.BoostVector());
                        if (m_top > 0) top2_p4[i_smear].Boost(-cms.BoostVector());

                        //ttbar solution in CM frame for bias check
                        float top1_p_CM = -999.0;
                        float top2_p_CM = -999.0;
                        if (m_top > 0)
                        {
                            top1_p_CM =  top1_p4[i_smear].P();
                            top2_p_CM =  top2_p4[i_smear].P();
                        }


                        top_costheta_cms = -999.0;
                        if (m_top > 0) top_costheta_cms = top1_p4[i_smear].Vect().Dot(cms.Vect()) / (top1_p4[i_smear].Vect().Mag() * cms.Vect().Mag());

                        if ( hyp_lt_id()[hypIdx] < 0 )
                        {
                            lepPlus.SetXYZT(
                                hyp_lt_p4()[hypIdx].x(),
                                hyp_lt_p4()[hypIdx].y(),
                                hyp_lt_p4()[hypIdx].z(),
                                hyp_lt_p4()[hypIdx].t()
                            );

                            // lepton energy scale systematic variation
                            if (fabs(hyp_lt_id()[hypIdx]) == 11) lepPlus*=leptonEnergyScaleFactor;

                            lepMinus.SetXYZT(
                                hyp_ll_p4()[hypIdx].x(),
                                hyp_ll_p4()[hypIdx].y(),
                                hyp_ll_p4()[hypIdx].z(),
                                hyp_ll_p4()[hypIdx].t()
                            );

                            // lepton energy scale systematic variation
                            if (fabs(hyp_ll_id()[hypIdx]) == 11) lepMinus*=leptonEnergyScaleFactor;

                        }
                        else
                        {
                            lepPlus.SetXYZT(
                                hyp_ll_p4()[hypIdx].x(),
                                hyp_ll_p4()[hypIdx].y(),
                                hyp_ll_p4()[hypIdx].z(),
                                hyp_ll_p4()[hypIdx].t()
                            );

                            // lepton energy scale systematic variation
                            if (fabs(hyp_ll_id()[hypIdx]) == 11) lepPlus*=leptonEnergyScaleFactor;

                            lepMinus.SetXYZT(
                                hyp_lt_p4()[hypIdx].x(),
                                hyp_lt_p4()[hypIdx].y(),
                                hyp_lt_p4()[hypIdx].z(),
                                hyp_lt_p4()[hypIdx].t()
                            );

                            // lepton energy scale systematic variation
                            if (fabs(hyp_lt_id()[hypIdx]) == 11) lepMinus*=leptonEnergyScaleFactor;

                        }


                        jet1.SetXYZT(
                            v_goodJets_cand_p4[0].x(),
                            v_goodJets_cand_p4[0].y(),
                            v_goodJets_cand_p4[0].z(),
                            v_goodJets_cand_p4[0].t()
                        );

                        jet2.SetXYZT(
                            v_goodJets_cand_p4[1].x(),
                            v_goodJets_cand_p4[1].y(),
                            v_goodJets_cand_p4[1].z(),
                            v_goodJets_cand_p4[1].t()
                        );                      

                        lep_charge_asymmetry = -999.0;
                        lep_charge_asymmetry = abs(lepPlus.Eta()) - abs(lepMinus.Eta());

                        lep_azimuthal_asymmetry = -999.0;
                        lep_azimuthal_asymmetry = cos(lepPlus.DeltaPhi(lepMinus));

                        lep_azimuthal_asymmetry_2 = -999.0;
                        lep_azimuthal_asymmetry_2 = acos(lep_azimuthal_asymmetry);

                        lep_pseudorap_diff = -999.0;
                        lep_pseudorap_diff = (lepPlus.Eta()) - (lepMinus.Eta());

                        float lep_cosalpha =  lepPlus.Vect().Dot( lepMinus.Vect() ) / (lepPlus.Vect().Mag() * lepMinus.Vect().Mag());
                        float lepPlus_phi = lepPlus.Phi();
                        float lepMinus_phi = lepMinus.Phi();
                        float lepPlus_Eta = lepPlus.Eta();
                        float lepMinus_Eta = lepMinus.Eta();
                        float lepPlus_Pt = lepPlus.Pt();
                        float lepMinus_Pt = lepMinus.Pt();

                        float jet_azimuthal_asymmetry = cos(jet1.DeltaPhi(jet2));
                        float jet_pseudorap_diff = jet1.Eta() - jet2.Eta();
                        float jet_cosalpha =  jet1.Vect().Dot( jet2.Vect() ) / (jet1.Vect().Mag() * jet2.Vect().Mag());
                        float jet1_phi = jet1.Phi();
                        float jet2_phi = jet2.Phi();

                        lepPlus.Boost(-cms.BoostVector());
                        lepMinus.Boost(-cms.BoostVector());

                        jet1.Boost(-cms.BoostVector());
                        jet2.Boost(-cms.BoostVector());

                        float lep_cosalpha_cms = lepPlus.Vect().Dot( lepMinus.Vect() ) / (lepPlus.Vect().Mag() * lepMinus.Vect().Mag());
                        float jet_cosalpha_cms =  jet1.Vect().Dot( jet2.Vect() ) / (jet1.Vect().Mag() * jet2.Vect().Mag());


                        if (m_top > 0) lepPlus.Boost(-top1_p4[i_smear].BoostVector());
                        if (m_top > 0) lepMinus.Boost(-top2_p4[i_smear].BoostVector());

                        lepPlus_costheta_cms = -999.0;
                        if (m_top > 0) lepPlus_costheta_cms = lepPlus.Vect().Dot(top1_p4[i_smear].Vect()) / (lepPlus.Vect().Mag() * top1_p4[i_smear].Vect().Mag());
                        lepMinus_costheta_cms = -999.0;
                        if (m_top > 0) lepMinus_costheta_cms = lepMinus.Vect().Dot(top2_p4[i_smear].Vect()) / (lepMinus.Vect().Mag() * top2_p4[i_smear].Vect().Mag());

                        top_spin_correlation = -999.0;
                        if (m_top > 0) top_spin_correlation = lepPlus_costheta_cms * lepMinus_costheta_cms;
                        //if we have gotten here, then all cuts have been passed
                        
                        nSelectedEvents = nSelectedEvents + 1.0 * weight;

                        //  weight_ = weight ;
                        //nSelectedEvents++;
                        ///some checks for events in the signal region

                        //  if(prefix=="data" && applyLeptonJetInvMassCut450) {

                        //    for(unsigned int i = 0; i < v_goodJets.size(); i++) {
                        //      int idx = v_goodJets[i];
                        //      LorentzVector goodjetp4 = pfjets_p4()[idx]*pfjets_corL1FastL2L3()[idx]*globalJESRescale;
                        //      //  if(passbTagging(v_goodJets[i],  "pfJets", "simpleSecondaryVertexHighEffBJetTag") ) {
                        //      cout <<"b tag Pt "<< goodjetp4.Pt()<<endl;
                        //      cout <<"b tag eta "<<goodjetp4.Eta()<<endl;
                        //      cout <<"b tag phi "<< goodjetp4.Phi()<<endl;
                        //      cout <<"TCHEL           "<<pfjets_trackCountingHighEffBJetTag()[v_goodJets[i]]<<endl;
                        //      cout <<"SSVHEM          "<<pfjets_simpleSecondaryVertexHighEffBJetTag()[v_goodJets[i]]<<endl;


                        //    }

                        //   if (abs(hyp_ll_id()[hypIdx]) == 13   ) cout << "ll muon err" <<cms2.mus_ptErr().at(hyp_ll_index()[hypIdx])<<endl;
                        //           if (abs(hyp_lt_id()[hypIdx]) == 13  ) cout << "lt muon err" <<cms2.mus_ptErr().at(hyp_lt_index()[hypIdx])<<endl;

                        //           //OSV1
                        //           if (abs(hyp_ll_id()[hypIdx]) == 11 ) cout << "ll isSpike" <<isSpikeElectron(hyp_ll_index()[hypIdx])<<endl;
                        //           if (abs(hyp_lt_id()[hypIdx]) == 11 ) cout << "lt isSpike" <<isSpikeElectron(hyp_lt_index()[hypIdx])<<endl;
                        //}
                        //
                        //  cout <<"type =    " <<type <<"  id_lt =   " <<id_lt <<"  id_ll =   " <<id_ll <<"  lt is from W = " <<leptonIsFromW(idx_lt, id_lt, false) << "  ll is from W = " << leptonIsFromW(idx_ll, id_ll, false)<<endl;

                        if (prefix == "ttdil")
                        {

                            // cout <<"run number =  " <<cms2.evt_run()<< ";  event number =  " << cms2.evt_event()<<endl;
                            // dumpDocLines();
                        }

                        int nvtx = 0;

                        for (size_t v = 0; v < cms2.vtxs_position().size(); ++v)
                        {
                            if (isGoodVertex(v)) ++nvtx;
                        }

                        //int ndavtx = 0;
                        ndavtx = 0;

                        for (size_t v = 0; v < cms2.davtxs_position().size(); ++v)
                        {
                            if (isGoodDAVertex(v)) ++ndavtx;
                        }

                        //cout << "ndavtx = " << ndavtx <<endl;
                        //unsigned int jetBin = v_goodJets.size();
                        jetBin = v_goodJets.size();
                        unsigned int nJets  = v_goodJets.size();
                        if (jetBin > 2)
                            jetBin = 2;


                        //which channel?
                        //int myType = 99;
                        myType = 99;
                        if (hyp_type()[hypIdx] == 3)                                myType = 0; // ee
                        if (hyp_type()[hypIdx] == 0)                                myType = 1; // mm
                        if (hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2)     myType = 2; // em
                        if (myType == 99)
                        {
                            cout << "Skipping unknown dilepton type = " << hyp_type()[hypIdx] << endl;
                            continue;
                        }

                        //hnJet[myType]                     ->Fill(nJets,               weight);
                        //hnJet[3]                          ->Fill(nJets,               weight);
                        //hnBtagJet[myType]                 ->Fill(nBtagJets,            weight);
                        //hnBtagJet[3]                      ->Fill(nBtagJets,            weight);

                        //hnVtx[myType]                     ->Fill(ndavtx,               weight);
                        //hnVtx[3]                          ->Fill(ndavtx,               weight);


                        fillHistos( hnJet, nJets,  weight, myType, jetBin, Nsolns);
                        fillHistos( hnBtagJet, nBtagJets,  weight, myType, jetBin, Nsolns);
                        fillHistos( hnVtx, ndavtx,  weight, myType, jetBin, Nsolns);

                        int Nsolns_for_hist = AMWTweight.size();
                        if (Nsolns_for_hist == 0) Nsolns_for_hist = -999;
                        fillHistos( hNsolns, Nsolns_for_hist,  weight, myType, jetBin, Nsolns);
                        fillHistos( hmaxAMWTweight, ( maxAMWTweight > 0 ? maxAMWTweight : -999 ),  weight, myType, jetBin, Nsolns);
                        fillHistos( hsumAMWTweight, ( sumAMWTweight > 0 ? sumAMWTweight : -999 ),  weight, myType, jetBin, Nsolns);
                        fillHistos( haveAMWTweight, ( aveAMWTweight > 0 ? aveAMWTweight : -999 ),  weight, myType, jetBin, Nsolns);
                        fillHistos( hAMWTweightnojetsmear, (m_top_nojetsmear > 0 ? AMWTweight_nojetsmear[0] : -999 ),  weight, myType, jetBin, Nsolns);

                        fillHistos( httMass, tt_mass ,  weight, myType, jetBin, Nsolns);
                        fillHistos( httMass_nojetsmear, tt_mass_nojetsmear ,  weight, myType, jetBin, Nsolns);
                        fillHistos( httpT, tt_pT ,  weight, myType, jetBin, Nsolns);
                        fillHistos( httpT_nojetsmear, tt_pT_nojetsmear ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepChargeAsym, lep_charge_asymmetry ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepAzimAsym, lep_azimuthal_asymmetry ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepAzimAsym2, lep_azimuthal_asymmetry_2 ,  weight, myType, jetBin, Nsolns);
                        if (m_top > 0)
                        {
                            fillHistos( htopSpinCorr, top_spin_correlation  ,  weight, myType, jetBin, Nsolns);
                            fillHistos( htopCosTheta, top_costheta_cms   ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hpseudorapiditydiff, top_pseudorapiditydiff_cms ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hrapiditydiff, top_rapiditydiff_cms ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hrapiditydiffMarco, top_rapiditydiff_Marco ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hlepCosTheta, lepPlus_costheta_cms  ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hlepCosTheta, lepMinus_costheta_cms  ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hlepPlusCosTheta,  lepPlus_costheta_cms, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepMinusCosTheta,  lepMinus_costheta_cms, weight, myType, jetBin, Nsolns);
                            fillHistos( httRapidity, ttRapidity ,  weight, myType, jetBin, Nsolns);
                            fillHistos( httRapidity2, ttRapidity2 ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hlepAngleBetweenCMS,  lep_cosalpha_cms, weight, myType, jetBin, Nsolns);
                            fillHistos( hjetAngleBetweenCMS,  jet_cosalpha_cms, weight, myType, jetBin, Nsolns);
                            fillHistos( hpseudorapiditydiff2,  top_pseudorapiditydiff2_cms, weight, myType, jetBin, Nsolns);
                            fillHistos( hrapiditydiff2,  top_rapiditydiff2_cms, weight, myType, jetBin, Nsolns);
                        }
                        fillHistos( htheleadinglepPt, lt_p4.Pt()  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hthesecondlepPt, ll_p4.Pt()  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepEta, lt_p4.Eta()  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepEta, ll_p4.Eta()  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hMET, p_met.first  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( htopMass, m_top ,  weight, myType, jetBin, Nsolns);
                        fillHistos( htopMass_nojetsmear, m_top_nojetsmear ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepRapDiff,  lep_pseudorap_diff, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepAngleBetween,  lep_cosalpha, weight, myType, jetBin, Nsolns);
                        fillHistos( hjetAzimAsym,  jet_azimuthal_asymmetry, weight, myType, jetBin, Nsolns);
                        fillHistos( hjetRapDiff,  jet_pseudorap_diff, weight, myType, jetBin, Nsolns);
                        fillHistos( hjetAngleBetween,  jet_cosalpha, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepPhi,  lepPlus_phi, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepPhi,  lepMinus_phi, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepPlusPhi,  lepPlus_phi, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepMinusPhi,  lepMinus_phi, weight, myType, jetBin, Nsolns);
                        fillHistos( hjetPhi,  jet1_phi, weight, myType, jetBin, Nsolns);
                        fillHistos( hjetPhi,  jet2_phi, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepPlusEta,  lepPlus_Eta, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepMinusEta,  lepMinus_Eta, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepPlusPt,  lepPlus_Pt, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepMinusPt,  lepMinus_Pt, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepPt,  lepPlus_Pt, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepPt,  lepMinus_Pt, weight, myType, jetBin, Nsolns);


                        if (v_goodJets_cand_p4.size() > 1)
                        {

                            // DYEst histos
                            // fill the mass histograms
                            double DYEst_mass = hyp_p4()[hypIdx].mass();
                            fillHistos( hdilMassNoMetDYEst , DYEst_mass ,  weight, myType, jetBin, Nsolns);
                            if ( p_met.first >= 30.  )
                                fillHistos( hdilMassWithMetDYEst  , DYEst_mass ,  weight, myType, jetBin, Nsolns);
                            // fill the met histograms for "in" and "out" regions
                            if (inZmassWindow(DYEst_mass) )
                            {
                                fillHistos( hmetInDYEst  , p_met.first ,  weight, myType, jetBin), Nsolns;
                            }
                            else
                            {
                                fillHistos( hmetOutDYEst  , p_met.first ,  weight, myType, jetBin, Nsolns);
                            }


                            fillHistos( hmassltb,  mass_ltb ,  weight, myType, jetBin), Nsolns;
                            fillHistos( hmassllb,  mass_llb ,  weight, myType, jetBin, Nsolns);
                            if (mass_llb > 170) fillHistos( hmassltb1Dmasscut,  mass_ltb ,  weight, myType, jetBin, Nsolns);
                            if (mass_ltb > 170) fillHistos( hmassllb1Dmasscut,  mass_llb ,  weight, myType, jetBin, Nsolns);
                            //randomly pick up mass_lb
                            int evtnumber = evt_event();
                            int evt_rad = rand() % 2;
                            //cout <<"random number " <<evt_rad <<endl;

                            if (evt_rad == 1)
                            {
                                fillHistos( hmasslb_2d,  mass_ltb , mass_llb,    weight, myType, jetBin, Nsolns);
                            }
                            else  //if(evt_rad ==0){
                            {
                                fillHistos( hmasslb_2d,  mass_llb , mass_ltb,    weight, myType, jetBin, Nsolns);
                            }


                            fillHistos( htheSumJetPt,      thefirstJet_pt + thesecondJet_pt ,  weight, myType, jetBin, Nsolns);
                            fillHistos( htheSumLepPt,    lt_p4.Pt() + ll_p4.Pt()   ,  weight, myType, jetBin, Nsolns);
                            fillHistos( htheSumBtagJetPt,  theSumBtagJetPt ,    weight, myType, jetBin, Nsolns);
                            fillHistos( hthefirstJetPt,    thefirstJet_pt  ,    weight, myType, jetBin, Nsolns);
                            fillHistos( hthesecondJetPt,   thesecondJet_pt  ,   weight, myType, jetBin, Nsolns);
                            fillHistos( hjetPt,  thefirstJet_pt ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hjetPt,  thesecondJet_pt ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hjetEta, v_goodJets_cand_p4.at(i_ltbjet).Eta()  ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hjetEta, v_goodJets_cand_p4.at(i_llbjet).Eta()  ,  weight, myType, jetBin, Nsolns);


                            fillHistos( httmasssm_2d,  tt_mass_nojetsmear , tt_mass,    weight, myType, jetBin, Nsolns);
                            fillHistos( htopmasssm_2d,  m_top_nojetsmear  , m_top,    weight, myType, jetBin, Nsolns);

                            fillHistos( htop1pCMsm_2d,   top1_p_CM_nojetsmear  , top1_p_CM,     weight, myType, jetBin, Nsolns);
                            fillHistos( htop2pCMsm_2d,   top2_p_CM_nojetsmear  , top2_p_CM,     weight, myType, jetBin, Nsolns);

                            fillHistos( htop1pTsm_2d,   top1_pt_nojetsmear  , top1_pt,     weight, myType, jetBin, Nsolns);
                            //fillHistos( htop1etasm_2d,  top1_eta_nojetsmear , top1_eta,    weight, myType, jetBin, Nsolns);
                            //fillHistos( htop1phism_2d,  top1_phi_nojetsmear , top1_phi,    weight, myType, jetBin, Nsolns);

                            fillHistos( htop2pTsm_2d,   top2_pt_nojetsmear  , top2_pt,     weight, myType, jetBin, Nsolns);
                            //fillHistos( htop2etasm_2d,  top2_eta_nojetsmear , top2_eta,    weight, myType, jetBin, Nsolns);
                            //fillHistos( htop2phism_2d,  top2_phi_nojetsmear , top2_phi,    weight, myType, jetBin, Nsolns);

                            //none weighted histograms
                            //fillHistos( habcd_2d,  mass_ltb  , mass_llb,    1./double(Nsolns), myType, jetBin, Nsolns);

                        }
                    }//!applyNoCuts

                    float dr_ltjet_gen = -999.0;
                    float dr_lljet_gen = -999.0;
                    float tt_mass_gen;
                    float tt_pT_gen;
                    float ttRapidity_gen;
		    float ttRapidity2_gen;
                    float top_costheta_cms_gen;
                    float lep_charge_asymmetry_gen;
                    float lep_azimuthal_asymmetry_gen;
                    float lep_azimuthal_asymmetry_gen2;
                    float top_spin_correlation_gen;
                    float lepPlus_costheta_cms_gen;
                    float lepMinus_costheta_cms_gen;
                    float massllbb_gen;
                    float llbbRapidityQuark_gen;
                    float llbbRapidityGluon_gen;
                    float top_rapiditydiff_cms_gen;
                    float top_rapiditydiff_Marco_gen;
                    float top_pseudorapiditydiff_cms_gen;
                    // generator level plots
                    //if(!isData && (prefix == "ttdil"|| prefix == "wprime400"|| prefix == "wprime600" || prefix == "axigluonR")){
                    if (!isData)
                    {



                        TLorentzVector topplus_genp_p4(0, 0, 0, 0), topminus_genp_p4(0, 0, 0, 0), cms_gen(0, 0, 0, 0), lepPlus_gen(0, 0, 0, 0), lepMinus_gen(0, 0, 0, 0), bPlus_gen(0, 0, 0, 0), bMinus_gen(0, 0, 0, 0);

                        bool from_gluon = true;
                        for (unsigned int i = 0; i < genps_p4().size(); i++)
                        {
                            if (genps_status()[i] == 3)
                            {
                                if ((genps_id_mother()[i] == 6 || genps_id_mother()[i] == 24 ))
                                {
                                    if ( (genps_id()[i] == -11 || genps_id()[i] == -13 ||  genps_id()[i] == -15) )
                                    {
                                        lepPlus_gen.SetXYZT(genps_p4()[i].x(),
                                                            genps_p4()[i].y(),
                                                            genps_p4()[i].z(),
                                                            genps_p4()[i].t()
                                                           );
                                    }
                                    else if ( genps_id()[i] == 5)
                                    {
                                        bPlus_gen.SetXYZT(genps_p4()[i].x(),
                                                          genps_p4()[i].y(),
                                                          genps_p4()[i].z(),
                                                          genps_p4()[i].t()
                                                         );
                                    }
                                }
                                else if ( (genps_id_mother()[i] == -6 || genps_id_mother()[i] == -24 ))
                                {
                                    if ( (genps_id()[i] == 11 || genps_id()[i] == 13 ||  genps_id()[i] == 15) )
                                    {

                                        lepMinus_gen.SetXYZT( genps_p4()[i].x(),
                                                              genps_p4()[i].y(),
                                                              genps_p4()[i].z(),
                                                              genps_p4()[i].t()
                                                            );
                                    }
                                    else if ( genps_id()[i] == -5)
                                    {
                                        bMinus_gen.SetXYZT(genps_p4()[i].x(),
                                                           genps_p4()[i].y(),
                                                           genps_p4()[i].z(),
                                                           genps_p4()[i].t()
                                                          );
                                    }
                                }

                                if (genps_id()[i] == 6 )
                                {
                                    topplus_genp_p4.SetXYZT( genps_p4()[i].x(),
                                                             genps_p4()[i].y(),
                                                             genps_p4()[i].z(),
                                                             genps_p4()[i].t()
                                                           );
                                    if (abs(genps_id_mother()[i]) == 21)
                                    {
                                        from_gluon = true;
                                    }
                                    else
                                    {
                                        from_gluon = false;
                                    }
                                }
                                else if (genps_id()[i] == -6 )
                                {
                                    topminus_genp_p4.SetXYZT( genps_p4()[i].x(),
                                                              genps_p4()[i].y(),
                                                              genps_p4()[i].z(),
                                                              genps_p4()[i].t()
                                                            );
                                    if (abs(genps_id_mother()[i]) == 21)
                                    {
                                        from_gluon = true;
                                    }
                                    else
                                    {
                                        from_gluon = false;
                                    }

                                }

                            }
                        }
                        massllbb_gen = (lepPlus_gen + bPlus_gen + lepMinus_gen + bMinus_gen).M();

                        float m_topminus_gen = topminus_genp_p4.M();
                        float m_topplus_gen = topplus_genp_p4.M();

                        tt_mass_gen = (topplus_genp_p4 + topminus_genp_p4).M();
                        ttRapidity_gen = topplus_genp_p4.Rapidity() + topminus_genp_p4.Rapidity();
			ttRapidity2_gen = (topplus_genp_p4 + topminus_genp_p4).Rapidity();
                        //ttRapidity_gen = topplus_genp_p4.Eta() + topminus_genp_p4.Eta();

                        top_rapiditydiff_cms_gen = (topplus_genp_p4.Rapidity() - topminus_genp_p4.Rapidity()) * (topplus_genp_p4.Rapidity() + topminus_genp_p4.Rapidity());
                        top_pseudorapiditydiff_cms_gen = abs(topplus_genp_p4.Eta()) - abs(topminus_genp_p4.Eta());
                        top_rapiditydiff_Marco_gen = abs(topplus_genp_p4.Rapidity()) - abs(topminus_genp_p4.Rapidity());


                        cms_gen = topplus_genp_p4 + topminus_genp_p4;
                        tt_pT_gen = cms_gen.Pt();
                        topplus_genp_p4.Boost(-cms_gen.BoostVector());
                        topminus_genp_p4.Boost(-cms_gen.BoostVector());
                        top_costheta_cms_gen = topplus_genp_p4.Vect().Dot(cms_gen.Vect()) / (topplus_genp_p4.Vect().Mag() * cms_gen.Vect().Mag());


                        lep_charge_asymmetry_gen = abs(lepPlus_gen.Eta()) - abs(lepMinus_gen.Eta());
                        lep_azimuthal_asymmetry_gen = cos(lepPlus_gen.DeltaPhi(lepMinus_gen));
                        lep_azimuthal_asymmetry_gen2 = acos(lep_azimuthal_asymmetry_gen);

                        lepPlus_gen.Boost(-cms_gen.BoostVector());
                        lepPlus_gen.Boost(-topplus_genp_p4.BoostVector());
                        lepMinus_gen.Boost(-cms_gen.BoostVector());
                        lepMinus_gen.Boost(-topminus_genp_p4.BoostVector());

                        lepPlus_costheta_cms_gen = lepPlus_gen.Vect().Dot(topplus_genp_p4.Vect()) / (lepPlus_gen.Vect().Mag() * topplus_genp_p4.Vect().Mag());
                        lepMinus_costheta_cms_gen = lepMinus_gen.Vect().Dot(topminus_genp_p4.Vect()) / (lepMinus_gen.Vect().Mag() * topminus_genp_p4.Vect().Mag());

                        top_spin_correlation_gen = lepPlus_costheta_cms_gen * lepMinus_costheta_cms_gen;

                        if (!applyNoCuts)
                        {
                            float tt_mass_pull = (tt_mass - tt_mass_gen) / tt_mass_gen;
                            float tt_mass_pull_nojetsmear = (tt_mass_nojetsmear - tt_mass_gen) / tt_mass_gen;
                            float tt_mass_diff = (tt_mass - tt_mass_gen);
                            float tt_mass_diff_nojetsmear = (tt_mass_nojetsmear - tt_mass_gen);

                            float top_mass_diff_plus = (m_top - m_topplus_gen);
                            float top_mass_diff_minus = (m_top - m_topminus_gen);
                            float top_mass_nojetsmear_diff_plus = (m_top_nojetsmear - m_topplus_gen);
                            float top_mass_nojetsmear_diff_minus = (m_top_nojetsmear - m_topminus_gen);

                            float top_p_diff_plus  = -999;
                            float top_p_diff_minus = -999;
                            if (m_top > 0) top_p_diff_plus = (top1_p4[i_smear].P() - topplus_genp_p4.P());
                            if (m_top > 0) top_p_diff_minus = (top2_p4[i_smear].P() - topminus_genp_p4.P());

                            fillHistos( httMass_2d, tt_mass_gen , tt_mass,  weight, myType, jetBin, Nsolns);
                            fillHistos( httMass_nojetsmear_2d, tt_mass_gen , tt_mass_nojetsmear,  weight, myType, jetBin, Nsolns);
                            fillHistos( httpT_2d, tt_pT_gen , tt_pT,  weight, myType, jetBin, Nsolns);
                            fillHistos( httpT_nojetsmear_2d, tt_pT_gen , tt_pT_nojetsmear,  weight, myType, jetBin, Nsolns);
                            if (m_top > 0) fillHistos( htopP_2d, topplus_genp_p4.P() , top1_p4[i_smear].P(),  weight, myType, jetBin, Nsolns);
                            if (m_top_nojetsmear > 0) fillHistos( htopP_nojetsmear_2d, topplus_genp_p4.P() , top1_nojetsmear_p4[0].P(),  weight, myType, jetBin, Nsolns);

                            if (m_top > 0)
                            {
                                fillHistos( httMass_pull, tt_mass_pull,  weight, myType, jetBin, Nsolns);
                                fillHistos( httMass_nojetsmear_pull, tt_mass_pull_nojetsmear,  weight, myType, jetBin, Nsolns);
                                fillHistos( httMass_diff, tt_mass_diff,  weight, myType, jetBin, Nsolns);
                                fillHistos( httMass_nojetsmear_diff, tt_mass_diff_nojetsmear,  weight, myType, jetBin, Nsolns);
                                fillHistos( htopMass_diff_plus, top_mass_diff_plus,  weight, myType, jetBin, Nsolns);
                                fillHistos( htopMass_diff_minus, top_mass_diff_minus,  weight, myType, jetBin, Nsolns);
                                fillHistos( htopMass_nojetsmear_diff_plus, top_mass_nojetsmear_diff_plus,  weight, myType, jetBin, Nsolns);
                                fillHistos( htopMass_nojetsmear_diff_minus, top_mass_nojetsmear_diff_minus,  weight, myType, jetBin, Nsolns);
                                fillHistos( htopPCM_diff_plus, top_p_diff_plus,  weight, myType, jetBin, Nsolns);
                                fillHistos( htopPCM_diff_minus, top_p_diff_minus,  weight, myType, jetBin, Nsolns);
                            }



                        }
                        fillHistos( htopMass_plus_gen, m_topplus_gen ,  weight, myType, jetBin, Nsolns);
                        fillHistos( htopMass_minus_gen, m_topminus_gen ,  weight, myType, jetBin, Nsolns);
                        fillHistos( httMass_gen, tt_mass_gen ,  weight, myType, jetBin, Nsolns);
                        fillHistos( httpT_gen, tt_pT_gen ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepChargeAsym_gen, lep_charge_asymmetry_gen ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepAzimAsym_gen, lep_azimuthal_asymmetry_gen ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepAzimAsym2_gen, acos(lep_azimuthal_asymmetry_gen) ,  weight, myType, jetBin, Nsolns);
                        if (m_top > 0 || applyNoCuts)
                        {
                            fillHistos( htopSpinCorr_gen, top_spin_correlation_gen  ,  weight, myType, jetBin, Nsolns);
                            fillHistos( htopCosTheta_gen, top_costheta_cms_gen   ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hlepCosTheta_gen, lepPlus_costheta_cms_gen  ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hlepCosTheta_gen, lepMinus_costheta_cms_gen  ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hlepPlusCosTheta_gen, lepPlus_costheta_cms_gen  ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hlepMinusCosTheta_gen, lepMinus_costheta_cms_gen  ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hpseudorapiditydiff_gen, top_pseudorapiditydiff_cms_gen ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hrapiditydiff_gen, top_rapiditydiff_cms_gen ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hrapiditydiffMarco_gen, top_rapiditydiff_Marco_gen ,  weight, myType, jetBin, Nsolns);
                        }


                        if (m_top > 0 || applyNoCuts)
                        {
                            //2D unfolding requires ttbar solution even for the purely leptonic variables
                            fillHistos( hlepChargeAsym_gen2d, lep_charge_asymmetry_gen ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepAzimAsym_gen2d, lep_azimuthal_asymmetry_gen ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepAzimAsym2_gen2d, acos(lep_azimuthal_asymmetry_gen) ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( htopSpinCorr_gen2d, top_spin_correlation_gen  ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( htopCosTheta_gen2d, top_costheta_cms_gen   ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepCosTheta_gen2d, lepPlus_costheta_cms_gen  ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepCosTheta_gen2d, lepMinus_costheta_cms_gen  ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepPlusCosTheta_gen2d, lepPlus_costheta_cms_gen  ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepMinusCosTheta_gen2d, lepMinus_costheta_cms_gen  ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hpseudorapiditydiff_gen2d, top_pseudorapiditydiff_cms_gen ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hrapiditydiff_gen2d, top_rapiditydiff_cms_gen ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hrapiditydiffMarco_gen2d, top_rapiditydiff_Marco_gen ,  tt_mass_gen, weight, myType, jetBin, Nsolns);

			    fillHistos( hlepChargeAsym_ttpT_gen2d, lep_charge_asymmetry_gen ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepAzimAsym_ttpT_gen2d, lep_azimuthal_asymmetry_gen ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepAzimAsym2_ttpT_gen2d, acos(lep_azimuthal_asymmetry_gen) ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( htopSpinCorr_ttpT_gen2d, top_spin_correlation_gen  ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( htopCosTheta_ttpT_gen2d, top_costheta_cms_gen   ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepCosTheta_ttpT_gen2d, lepPlus_costheta_cms_gen  ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepCosTheta_ttpT_gen2d, lepMinus_costheta_cms_gen  ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepPlusCosTheta_ttpT_gen2d, lepPlus_costheta_cms_gen  ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hlepMinusCosTheta_ttpT_gen2d, lepMinus_costheta_cms_gen  ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hpseudorapiditydiff_ttpT_gen2d, top_pseudorapiditydiff_cms_gen ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hrapiditydiff_ttpT_gen2d, top_rapiditydiff_cms_gen ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                            fillHistos( hrapiditydiffMarco_ttpT_gen2d, top_rapiditydiff_Marco_gen ,  tt_pT_gen, weight, myType, jetBin, Nsolns);

			    fillHistos( hlepChargeAsym_ttRapidity2_gen2d, lep_charge_asymmetry_gen ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                            fillHistos( hlepAzimAsym_ttRapidity2_gen2d, lep_azimuthal_asymmetry_gen ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                            fillHistos( hlepAzimAsym2_ttRapidity2_gen2d, acos(lep_azimuthal_asymmetry_gen) ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                            fillHistos( htopSpinCorr_ttRapidity2_gen2d, top_spin_correlation_gen  ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                            fillHistos( htopCosTheta_ttRapidity2_gen2d, top_costheta_cms_gen   ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                            fillHistos( hlepCosTheta_ttRapidity2_gen2d, lepPlus_costheta_cms_gen  ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                            fillHistos( hlepCosTheta_ttRapidity2_gen2d, lepMinus_costheta_cms_gen  ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                            fillHistos( hlepPlusCosTheta_ttRapidity2_gen2d, lepPlus_costheta_cms_gen  ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                            fillHistos( hlepMinusCosTheta_ttRapidity2_gen2d, lepMinus_costheta_cms_gen  ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                            fillHistos( hpseudorapiditydiff_ttRapidity2_gen2d, top_pseudorapiditydiff_cms_gen ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                            fillHistos( hrapiditydiff_ttRapidity2_gen2d, top_rapiditydiff_cms_gen ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                            fillHistos( hrapiditydiffMarco_ttRapidity2_gen2d, top_rapiditydiff_Marco_gen ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        }


                        if (!applyNoCuts)
                        {
                            fillHistos( hlepChargeAsym_2d, lep_charge_asymmetry_gen , lep_charge_asymmetry,  weight, myType, jetBin, Nsolns);
                            fillHistos( hlepAzimAsym_2d, lep_azimuthal_asymmetry_gen , lep_azimuthal_asymmetry,  weight, myType, jetBin, Nsolns);
                            fillHistos( htopSpinCorr_2d, top_spin_correlation_gen, top_spin_correlation ,  weight, myType, jetBin, Nsolns);
                            fillHistos( htopCosTheta_2d, top_costheta_cms_gen, top_costheta_cms   ,  weight, myType, jetBin, Nsolns);
                            fillHistos( hlepCosTheta_2d, lepPlus_costheta_cms_gen  , lepPlus_costheta_cms, weight, myType, jetBin, Nsolns);
                        }

                        if (!applyNoCuts)
                        {
                            fillHistos( hlepChargeAsymGenDiff, -lep_charge_asymmetry_gen + lep_charge_asymmetry,  weight, myType, jetBin, Nsolns);
                            fillHistos( hlepAzimAsymGenDiff, -lep_azimuthal_asymmetry_gen + lep_azimuthal_asymmetry,  weight, myType, jetBin, Nsolns);
                            fillHistos( hlepAzimAsym2GenDiff, -acos(lep_azimuthal_asymmetry_gen) + acos(lep_azimuthal_asymmetry),  weight, myType, jetBin, Nsolns);
                            if ( m_top > 0 )
                            {
                                fillHistos( htopSpinCorrGenDiff, -top_spin_correlation_gen + top_spin_correlation ,  weight, myType, jetBin, Nsolns);
                                fillHistos( htopCosThetaGenDiff, -top_costheta_cms_gen + top_costheta_cms   ,  weight, myType, jetBin, Nsolns);
                                fillHistos( hlepCosThetaGenDiff, -lepPlus_costheta_cms_gen + lepPlus_costheta_cms, weight, myType, jetBin, Nsolns);
                                fillHistos( hlepCosThetaGenDiff, -lepMinus_costheta_cms_gen + lepMinus_costheta_cms, weight, myType, jetBin, Nsolns);
                                fillHistos( hlepPlusCosThetaGenDiff, -lepPlus_costheta_cms_gen + lepPlus_costheta_cms, weight, myType, jetBin, Nsolns);
                                fillHistos( hlepMinusCosThetaGenDiff, -lepMinus_costheta_cms_gen + lepMinus_costheta_cms, weight, myType, jetBin, Nsolns);
                                fillHistos( hpseudorapiditydiffGenDiff, -top_pseudorapiditydiff_cms_gen + top_pseudorapiditydiff_cms , weight, myType, jetBin, Nsolns);
                                fillHistos( hrapiditydiffGenDiff, -top_rapiditydiff_cms_gen + top_rapiditydiff_cms , weight, myType, jetBin, Nsolns);
                                fillHistos( hrapiditydiffMarcoGenDiff, -top_rapiditydiff_Marco_gen + top_rapiditydiff_Marco , weight, myType, jetBin, Nsolns);
                            }
                        }

                        if (from_gluon == true)
                        {
                            fillHistos( httMassGluongenp, tt_mass_gen  ,  weight, myType, jetBin, Nsolns);
                            fillHistos( httRapidityGluongenp, ttRapidity_gen  ,  weight, myType, jetBin, Nsolns);
                            llbbRapidityGluon_gen =  (lepPlus_gen + bPlus_gen).Rapidity() + (lepMinus_gen + bMinus_gen).Rapidity();
                            fillHistos( hllbbRapidityGluongenp, llbbRapidityGluon_gen  ,  weight, myType, jetBin, Nsolns);
                        }
                        else
                        {
                            fillHistos( httMassQuarkgenp, tt_mass_gen  ,  weight, myType, jetBin, Nsolns);
                            fillHistos( httRapidityQuarkgenp, ttRapidity_gen  ,  weight, myType, jetBin, Nsolns);
                            llbbRapidityQuark_gen =  (lepPlus_gen + bPlus_gen).Rapidity() + (lepMinus_gen + bMinus_gen).Rapidity();
                            fillHistos( hllbbRapidityQuarkgenp,  llbbRapidityQuark_gen ,  weight, myType, jetBin, Nsolns);
                        }

                    }//only for mc

                    // fill ntuples
                    if (createBabyNtuples && !applyNoCuts)
                    {
                        InitBabyNtuple();
                        run_ = cms2.evt_run();
                        ls_ = cms2.evt_lumiBlock();
                        evt_ = cms2.evt_event();
                        weight_ = weight;
                        Nsolns_ = Nsolns;
                        massltb_ = mass_ltb;
                        massllb_ = mass_llb;
                        dr_ltjet_gen_ = dr_ltjet_gen ;
                        dr_lljet_gen_ = dr_lljet_gen ;
                        ndavtx_ = ndavtx;
                        tt_mass_ = tt_mass;
                        t_mass_  = m_top;
                        ttRapidity_ = ttRapidity;
			ttRapidity2_ = ttRapidity2;
			ttRapidity2_gen_ = ttRapidity2_gen;
                        ttPt_ = tt_pT;
			ttPt_gen_ = tt_pT_gen;
			lep_charge_asymmetry_ = lep_charge_asymmetry;
                        lep_pseudorap_diff_ =  lep_pseudorap_diff;
                        lep_azimuthal_asymmetry_ = lep_azimuthal_asymmetry;
                        lep_azimuthal_asymmetry2_ = lep_azimuthal_asymmetry_2;
                        top_spin_correlation_ = top_spin_correlation;
                        top_costheta_cms_     = top_costheta_cms;
                        top_rapiditydiff_cms_ = top_rapiditydiff_cms;
                        top_rapiditydiff_Marco_ = top_rapiditydiff_Marco;
                        top_pseudorapiditydiff_cms_ = top_pseudorapiditydiff_cms;
                        top_rapiditydiff_cms_gen_ = top_rapiditydiff_cms_gen;
                        top_rapiditydiff_Marco_gen_ = top_rapiditydiff_Marco_gen;
                        top_pseudorapiditydiff_cms_gen_ = top_pseudorapiditydiff_cms_gen;
                        lepPlus_costheta_cms_ = lepPlus_costheta_cms;
                        lepMinus_costheta_cms_ = lepMinus_costheta_cms;
                        tt_mass_gen_ = tt_mass_gen;
                        ttRapidity_gen_ = ttRapidity_gen;
                        lep_charge_asymmetry_gen_ = lep_charge_asymmetry_gen;
                        lep_azimuthal_asymmetry_gen_ = lep_azimuthal_asymmetry_gen;
                        lep_azimuthal_asymmetry_gen2_ = lep_azimuthal_asymmetry_gen2;
                        top_spin_correlation_gen_ = top_spin_correlation_gen;
                        top_costheta_cms_gen_     = top_costheta_cms_gen;
                        lepPlus_costheta_cms_gen_ = lepPlus_costheta_cms_gen;
                        lepMinus_costheta_cms_gen_ = lepMinus_costheta_cms_gen;
                        massllbb_ = massllbb;
                        massllbb_gen_ = massllbb_gen;
                        //llbbRapidityQuark_gen_=llbbRapidityQuark_gen;
                        //  llbbRapidityGluon_gen_=llbbRapidityGluon_gen;
                        babyTree_->Fill();

                    }

                } //jet smearing loop

            }//good hypothesis loop

        } // closes loop over events

        if (applyNoCuts) cout << "number of events (no cuts) before and after vertex weighting              =  " << nEvents_noCuts_novtxweight << "   " << nEvents_noCuts << endl;
        if (applyNoCuts) cout << "number of dilepton events (no cuts) before and after vertex weighting              =  " << nEvents_noCuts_novtxweight_dil << "   " << nEvents_noCuts_dil << endl;

        //float nEvents_primary = cms2.evt_nEvts();
        //cout << "acceptance                       =  " << (1.0*nSelectedEvents)/(nEvents_primary*kFactor * evt_scale1fb() * lumi) <<endl;


    }  // closes loop over files

    if (createBabyNtuples && !applyNoCuts)CloseBabyNtuple();
    return;

} // closes myLooper function


//------------------------------------------
// Initialize baby ntuple variables
//------------------------------------------
void topAFB_looper::InitBabyNtuple ()
{
    run_ = -1;
    ls_  = -1;
    evt_ = -1;
    t_mass_ = -999;
    weight_ = 1;
    Nsolns_ = -1;
    massltb_ = -999;
    massllb_ = -999;
    dr_ltjet_gen_ = -999.0;
    dr_lljet_gen_  = -999.0;
    ndavtx_ = -999;
    tt_mass_ = -999.0;
    ttRapidity_ = -999.0;
    ttRapidity2_ = -999.0;
    ttRapidity2_gen_ = -999.0;
    ttPt_ = -999.0;
    ttPt_gen_ = -999.0;
    lep_charge_asymmetry_ = -999.0;
    lep_pseudorap_diff_ = -999.0;
    lep_azimuthal_asymmetry_ = -999.0;
    lep_azimuthal_asymmetry2_ = -999.0;
    top_spin_correlation_ = -999.0;
    top_costheta_cms_     = -999.0;
    top_rapiditydiff_cms_           = -999.0;
    top_rapiditydiff_Marco_           = -999.0;
    top_pseudorapiditydiff_cms_     = -999.0;
    top_rapiditydiff_Marco_gen_           = -999.0;
    top_rapiditydiff_cms_gen_           = -999.0;
    top_pseudorapiditydiff_cms_gen_     = -999.0;
    lepPlus_costheta_cms_ = -999.0;
    lepMinus_costheta_cms_ = -999.0;
    tt_mass_gen_ = -999.0;
    ttRapidity_gen_ = -999.0;
    lep_charge_asymmetry_gen_ = -999.0;
    lep_azimuthal_asymmetry_gen_ = -999.0;
    lep_azimuthal_asymmetry_gen2_ = -999.0;
    top_spin_correlation_gen_ = -999.0;
    top_costheta_cms_gen_     = -999.0;
    lepPlus_costheta_cms_gen_ = -999.0;
    lepMinus_costheta_cms_gen_ = -999.0;
    massllbb_gen_ = -999.0;
    //  llbbRapidityQuark_gen_=-999.0;
    // llbbRapidityGluon_gen_=-999.0;
    massllbb_ = -999.0;

}
//-------------------------------------
// Book the baby ntuple
//-------------------------------------
void topAFB_looper::MakeBabyNtuple(const char *babyFilename)
{
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();
    babyFile_ = new TFile(Form("%s", babyFilename), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");

    babyTree_->Branch("run",                   &run_,                 "run/I"                  );
    babyTree_->Branch("ls",                    &ls_,                  "ls/I"                   );
    babyTree_->Branch("evt",                   &evt_,                 "evt/I"                  );
    babyTree_->Branch("t_mass",                &t_mass_,              "t_mass/F"               );
    babyTree_->Branch("weight",                &weight_,              "weight/D"               );
    babyTree_->Branch("Nsolns",                &Nsolns_,              "Nsolns/I"               );
    babyTree_->Branch("massltb",               &massltb_,             "massltb/F"              );
    babyTree_->Branch("massllb",               &massllb_,             "massllb/F"              );
    babyTree_->Branch("dr_ltjet_gen",          &dr_ltjet_gen_,        "dr_ltjet_gen/F"         );
    babyTree_->Branch("dr_lljet_gen",          &dr_lljet_gen_,        "dr_lljet_gen/F"         );
    babyTree_->Branch("ndavtx",                &ndavtx_,              "ndavtx/I"               );
    babyTree_->Branch("tt_mass",               &tt_mass_,             "tt_mass/F"              );
    babyTree_->Branch("ttRapidity",            &ttRapidity_,          "ttRapidity/F"           );
    babyTree_->Branch("ttRapidity2",            &ttRapidity2_,          "ttRapidity2/F"           );
    babyTree_->Branch("ttRapidity2_gen",            &ttRapidity2_gen_,          "ttRapidity2_gen/F"           );
    babyTree_->Branch("ttPt",                   &ttPt_,          "ttPt/F"           );
    babyTree_->Branch("ttPt_gen",                   &ttPt_gen_,          "ttPt_gen/F"           );
    babyTree_->Branch("lep_charge_asymmetry",  &lep_charge_asymmetry_, "lep_charge_asymmetry/F" );
    babyTree_->Branch("lep_pseudorap_diff",    &lep_pseudorap_diff_,  "lep_pseudorap_diff/F"   );
    babyTree_->Branch("lep_azimuthal_asymmetry", &lep_azimuthal_asymmetry_,  "lep_azimuthal_asymmetry/F"   );
    babyTree_->Branch("lep_azimuthal_asymmetry2", &lep_azimuthal_asymmetry2_,  "lep_azimuthal_asymmetry2/F"   );
    babyTree_->Branch("top_spin_correlation",  &top_spin_correlation_, "top_spin_correlation/F" );
    babyTree_->Branch("top_costheta_cms",      &top_costheta_cms_,    "top_costheta_cms/F"     );
    babyTree_->Branch("lep_costheta_cms",      &lepPlus_costheta_cms_, "lep_costheta_cms/F"     );
    babyTree_->Branch("lepMinus_costheta_cms",      &lepMinus_costheta_cms_, "lepMinus_costheta_cms/F"     );
    babyTree_->Branch("top_rapidtiydiff_cms",      &top_rapiditydiff_cms_, "top_rapiditydiff_cms/F"     );
    babyTree_->Branch("top_rapidtiydiff_Marco",      &top_rapiditydiff_Marco_, "top_rapiditydiff_Marco/F"     );
    babyTree_->Branch("top_pseudorapidtiydiff_cms",      &top_pseudorapiditydiff_cms_, "top_pseudorapiditydiff_cms/F"     );
    babyTree_->Branch("top_rapidtiydiff_cms_gen",      &top_rapiditydiff_cms_gen_, "top_rapiditydiff_cms_gen/F"     );
    babyTree_->Branch("top_rapidtiydiff_Marco_gen",      &top_rapiditydiff_Marco_gen_, "top_rapiditydiff_Marco_gen/F"     );
    babyTree_->Branch("top_pseudorapidtiydiff_cms_gen",      &top_pseudorapiditydiff_cms_gen_, "top_pseudorapiditydiff_cms_gen/F"     );
    babyTree_->Branch("tt_mass_gen",           &tt_mass_gen_,          "tt_mass_gen/F"              );
    babyTree_->Branch("ttRapidity_gen",            &ttRapidity_gen_,          "ttRapidity_gen/F"            );
    babyTree_->Branch("lep_charge_asymmetry_gen",  &lep_charge_asymmetry_gen_, "lep_charge_asymmetry_gen_/F" );
    babyTree_->Branch("lep_azimuthal_asymmetry_gen",    &lep_azimuthal_asymmetry_gen_,  "lep_azimuthal_asymmetry_gen_/F"   );
    babyTree_->Branch("lep_azimuthal_asymmetry_gen2",    &lep_azimuthal_asymmetry_gen2_,  "lep_azimuthal_asymmetry_gen2_/F"   );
    babyTree_->Branch("top_spin_correlation_gen",  &top_spin_correlation_gen_, "top_spin_correlation_gen/F"  );
    babyTree_->Branch("top_costheta_cms_gen",      &top_costheta_cms_gen_,    "top_costheta_cms_gen/F"      );
    babyTree_->Branch("lep_costheta_cms_gen",      &lepPlus_costheta_cms_gen_, "lep_costheta_cms_gen/F"      );
    babyTree_->Branch("lepMinus_costheta_cms_gen",      &lepMinus_costheta_cms_gen_, "lepMinus_costheta_cms_gen/F"      );
    babyTree_->Branch("massllbb",               &massllbb_,             "massllbb/F"              );
    babyTree_->Branch("massllbb_gen",               &massllbb_gen_,             "massllbb_gen/F"              );
    //  babyTree_->Branch("llbbRapidityQuark_gen",      &llbbRapidityQuark_gen_,    "llbbRapidityQuark_gen/F"     );
    // babyTree_->Branch("llbbRapidityGluon_gen",      &llbbRapidityGluon_gen_,    "llbbRapidityGluon_gen/F"     );


}
//----------------------------------
// Fill the baby
//----------------------------------
void topAFB_looper::FillBabyNtuple()
{
    babyTree_->Fill();
}
//--------------------------------
// Close the baby
//--------------------------------
void topAFB_looper::CloseBabyNtuple()
{
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();
}


void topAFB_looper::fillUnderOverFlow(TH1D *h1, float value, double weight, int Nsolns)
{
    double min = h1->GetXaxis()->GetXmin();
    double max = h1->GetXaxis()->GetXmax();

    if (value >= max) value = h1->GetBinCenter(h1->GetNbinsX());
    if (value <= min) value = h1->GetBinCenter(1);

    int bin_number = h1->FindBin(value);
    double orig_content = h1->GetBinContent(bin_number);
    double orig_error = h1->GetBinError(bin_number);

    //h1->Fill(value, weight);
    h1->SetBinContent( bin_number, orig_content + weight );
    h1->SetBinError( bin_number, sqrt( orig_error * orig_error + weight * weight * double(Nsolns) ) );
}

//--------------------------------------------------------------------

void topAFB_looper::fillUnderOverFlow(TH2D *h2, float xvalue, float yvalue, double weight, int Nsolns)
{
    double maxx = h2->GetXaxis()->GetXmax();
    double minx = h2->GetXaxis()->GetXmin();
    double maxy = h2->GetYaxis()->GetXmax();
    double miny = h2->GetYaxis()->GetXmin();

    if (xvalue >= maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
    if (xvalue <= minx) xvalue = h2->GetXaxis()->GetBinCenter(1);
    if (yvalue >= maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());
    if (yvalue <= miny) yvalue = h2->GetYaxis()->GetBinCenter(1);

    int bin_number = h2->FindBin(xvalue, yvalue);
    double orig_content = h2->GetBinContent(bin_number);
    double orig_error = h2->GetBinError(bin_number);

    //h2->Fill(xvalue, yvalue, weight);
    h2->SetBinContent( bin_number, orig_content + weight );
    h2->SetBinError( bin_number, sqrt( orig_error * orig_error + weight * weight * double(Nsolns) ) );
}


//--------------------------------------------------------------------

void topAFB_looper::fillOverFlow(TH1D *h1, float value, float weight)
{
    float max = h1->GetXaxis()->GetXmax();
    if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
    h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void topAFB_looper::fillOverFlow(TH2D *h2, float xvalue, float yvalue, float weight)
{
    float maxx = h2->GetXaxis()->GetXmax();
    float maxy = h2->GetYaxis()->GetXmax();

    if (xvalue > maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
    if (yvalue > maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());

    h2->Fill(xvalue, yvalue, weight);
}

//--------------------------------------------------------------------

void topAFB_looper::fillHistos(TH1D *h1[4][4], float value, double weight, int myType, int nJetsIdx, int Nsolns)
{
    //fillUnderOverFlow(h1[myType][nJetsIdx], value, weight, Nsolns);
    fillUnderOverFlow(h1[myType][3],        value, weight, Nsolns);
    //fillUnderOverFlow(h1[3][nJetsIdx],      value, weight, Nsolns);
    fillUnderOverFlow(h1[3][3],             value, weight, Nsolns);
}

//--------------------------------------------------------------------

void topAFB_looper::fillHistos(TH2D *h2[4][4], float xvalue, float yvalue, double weight, int myType, int nJetsIdx, int Nsolns)
{
    //fillUnderOverFlow(h2[myType][nJetsIdx], xvalue, yvalue, weight, Nsolns);
    fillUnderOverFlow(h2[myType][3],        xvalue, yvalue, weight, Nsolns);
    //fillUnderOverFlow(h2[3][nJetsIdx],      xvalue, yvalue, weight, Nsolns);
    fillUnderOverFlow(h2[3][3],             xvalue, yvalue, weight, Nsolns);
}

//--------------------------------------------------------------------

void topAFB_looper::fillHistos(TProfile *h2[4][4], float xvalue, float yvalue, int myType, int nJetsIdx)
{
    //h2[myType][nJetsIdx] -> Fill(xvalue, yvalue);
    h2[myType][3]        -> Fill(xvalue, yvalue);
    //h2[3][nJetsIdx]      -> Fill(xvalue, yvalue);
    h2[3][3]             -> Fill(xvalue, yvalue);
}

int topAFB_looper::antimatch4vector(const LorentzVector &lvec,
                                    const vector<LorentzVector> &vec)
{

    if ( vec.size() == 0 ) return -1;
    //cout << "size of vec = " << vec.size() << endl;
    double dR = 0.0;
    double x;
    int iret = -1;
    for ( unsigned int i = 0; i < vec.size(); ++i)
    {
        x = ROOT::Math::VectorUtil::DeltaR(lvec, vec[i]);
        if (x > dR )
        {
            dR = x;
            iret = i;
        }
    }
    return iret;
}
