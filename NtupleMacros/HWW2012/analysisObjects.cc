#include "analysisObjects.h"
#include "analysisSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/electronSelections.h"

#include "Math/VectorUtil.h"

//
// hyps
//

HypTypeInNtuples hypType(unsigned int i_hyp){
    HypTypeInNtuples type = HypTypeInNtuples(cms2.hyp_type().at(i_hyp));
    return type;
}

//
// met
//

double metValue()
{    
    return cms2.evt_pfmet(); 
}

double metPhiValue()
{ 
    return cms2.evt_pfmetPhi(); 
}

double sumetValue()
{    
    return cms2.evt_pfsumet(); 
}

//
// jets
//

WWJetType jetType(){
    return pfJet;
}

    std::vector<JetPair>
getJets(WWJetType type, int i_hyp, double etThreshold, double etaMax, bool applyJEC, FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3, bool sortJets, bool btag)
{

    return getJets(type, cms2.hyp_lt_p4()[i_hyp], cms2.hyp_ll_p4()[i_hyp], etThreshold, etaMax, applyJEC, jet_corrector_pfL1FastJetL2L3, sortJets, btag);

}

    std::vector<JetPair>
getJets(WWJetType type, LorentzVector &lt, LorentzVector &ll, double etThreshold, double etaMax, bool applyJEC, FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3, bool sortJets, bool btag)
{
    std::vector<JetPair> jets;
    const double vetoCone = 0.3;

    switch ( type ){
        case jptJet:
            for ( unsigned int i=0; i < cms2.jpts_p4().size(); ++i) {
                double jec = 1.0;
                if ( cms2.jpts_p4()[i].pt() * jec < etThreshold ) continue;
                if ( btag && !defaultBTag(type,i, jec) ) continue;
                if ( TMath::Abs(cms2.jpts_p4()[i].eta()) > etaMax ) continue;
                if ( (lt.Pt() > 0 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(lt, cms2.jpts_p4()[i])) < vetoCone) ||
                        (ll.Pt() > 0 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ll, cms2.jpts_p4()[i])) < vetoCone) ) continue;
                jets.push_back(JetPair(cms2.jpts_p4()[i] * jec,i));
            }
            break;
        case pfJet:
            for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
                double jec = 1.0;
                // cout << cms2.evt_event() << " \traw pt: " << cms2.pfjets_p4().at(i).pt() << endl;
                if(applyJEC){
		  			jet_corrector_pfL1FastJetL2L3->setRho(cms2.evt_ww_rho());
                    jet_corrector_pfL1FastJetL2L3->setJetA(cms2.pfjets_area().at(i));
                    jet_corrector_pfL1FastJetL2L3->setJetPt(cms2.pfjets_p4()[i].pt());
                    jet_corrector_pfL1FastJetL2L3->setJetEta(cms2.pfjets_p4()[i].eta());
                    double corr = jet_corrector_pfL1FastJetL2L3->getCorrection();
                    jec *= corr;
                    // cout << " \tL1FastJetL2L3 corr: " << corr << " \tpt: " << cms2.pfjets_p4().at(i).pt()*corr << endl;

                    // jec *= jetCorrection(cms2.pfjets_p4()[i], jet_corrector_pf);
                    // cout << " \tL2L3 corr: " << jetCorrection(cms2.pfjets_p4()[i], jet_corrector_pf) << endl;
                }
                //       if(applyFastJetCorrection){
                //  jec *= (1-cms2.evt_rho()*cms2.pfjets_area().at(i)/cms2.pfjets_p4().at(i).pt()); //*cms2.pfjets_cor().at(i);// It's full L1Fast*L2*L3
                //  cout << " \tL1 corr: " << (1-cms2.evt_rho()*cms2.pfjets_area().at(i)/cms2.pfjets_p4().at(i).pt()) << " \t" << 
                //    (1-cms2.evt_rho()*cms2.pfjets_area().at(i)/cms2.pfjets_p4().at(i).pt())*cms2.pfjets_p4().at(i).pt() << 
                //    " \t" << cms2.evt_rho() << " \t" << cms2.pfjets_area().at(i) << endl;
                //       }
                if ( cms2.pfjets_p4()[i].pt() * jec < etThreshold ) continue;
                if ( btag && !defaultBTag(type,i, jec) ) continue;
                if ( TMath::Abs(cms2.pfjets_p4()[i].eta()) > etaMax ) continue;
                //if (btag && cms2.pfjets_p4()[i].pt() * jec < 30. && fabs(jetDz(i,0))>2.) continue;	//newcuts
                if ( (lt.Pt() > 0 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(lt, cms2.pfjets_p4()[i])) < vetoCone) ||
                        (ll.Pt() > 0 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ll, cms2.pfjets_p4()[i])) < vetoCone) ) continue;

				if ( !passMVAJetId( cms2.pfjets_p4()[i].pt() * jec, cms2.pfjets_p4()[i].eta(), cms2.pfjets_mvavalue()[i], 2) ) continue;

                // cout << " \tpassed all cuts" << endl;
                jets.push_back(JetPair(cms2.pfjets_p4()[i] * jec,i));
            }
            break;
        case GenJet:
            for ( unsigned int i=0; i < cms2.genjets_p4().size(); ++i) {
                if ( cms2.genjets_p4()[i].pt() < etThreshold ) continue;
                if ( btag && !defaultBTag(type,i) ) continue;
                if ( TMath::Abs(cms2.genjets_p4()[i].eta()) > etaMax ) continue;
                if ( (lt.Pt() > 0 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(lt, cms2.genjets_p4()[i])) < vetoCone) ||
                        (ll.Pt() > 0 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ll, cms2.genjets_p4()[i])) < vetoCone) ) continue;

                jets.push_back(JetPair(cms2.genjets_p4()[i],i));
            }
            break;
            //      case CaloJet:
            //        for ( unsigned int i=0; i < cms2.jets_pat_jet_p4().size(); ++i) {
            //   if ( cms2.jets_pat_jet_p4()[i].pt() < etThreshold ) continue; // note that this is already corrected
            //   if ( btag && !defaultBTag(type,i) ) continue;
            //   if ( TMath::Abs(cms2.jets_pat_jet_p4()[i].eta()) > etaMax ) continue;
            //   if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(lt, cms2.jets_pat_jet_p4()[i])) < vetoCone ||
            //        TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ll, cms2.jets_pat_jet_p4()[i])) < vetoCone ) continue;
            //   jets.push_back(cms2.jets_pat_jet_p4()[i]);
            //        }
            //        break;
        case TrkJet:
            for ( unsigned int i=0; i < cms2.trkjets_p4().size(); ++i) {
                double jec = 1.0;
                if ( cms2.trkjets_p4()[i].pt() < etThreshold ) continue;
                if ( btag && !defaultBTag(type,i, jec) ) continue;
                if ( TMath::Abs(cms2.trkjets_p4()[i].eta()) > etaMax ) continue;
                if ( (lt.Pt() > 0 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(lt, cms2.genjets_p4()[i])) < vetoCone) ||
                        (ll.Pt() > 0 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ll, cms2.genjets_p4()[i])) < vetoCone) ) continue;

                jets.push_back(JetPair(cms2.trkjets_p4()[i] * jec,i));
            }
            break;
        default:
            std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
    }
    if ( sortJets ) std::sort(jets.begin(), jets.end(), comparePt);
    return jets;
}

std::vector<JetPair> getDefaultJets(unsigned int i_hyp, bool applyJEC, FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3, bool btagged){
    return getJets(jetType(), i_hyp, 30, 4.7, applyJEC, jet_corrector_pfL1FastJetL2L3, false, btagged); // V1 // new cut
}

unsigned int numberOfJets(unsigned int i_hyp, bool applyJEC, FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3){
    return getDefaultJets(i_hyp, applyJEC, jet_corrector_pfL1FastJetL2L3, false).size();
}

Bool_t comparePt(JetPair lv1, JetPair lv2) {
    return lv1.first.pt() > lv2.first.pt();
}

//
// Btagging
//

double BTag(LorentzVector jetP4){
    int refJet = -1;
    for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
        if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(jetP4,cms2.pfjets_p4()[i])) > 0.3 ) continue;
        refJet = i;
    }
    if (refJet == -1){
        // std::cout << "Warning: failed to find a matching jet for b-tagging." << std::endl; 
        return 0.0;
    }
    return cms2.pfjets_trackCountingHighEffBJetTag().at(refJet);
}

double BTag(WWJetType type, unsigned int iJet){
    switch ( type ) {
        case jptJet:
            return BTag(cms2.jpts_p4().at(iJet));
            break;
        case CaloJet:
            return cms2.jets_trackCountingHighEffBJetTag()[iJet];
            break;
        case pfJet:
            return cms2.pfjets_trackCountingHighEffBJetTag()[iJet];
            break;
        default:
            std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
            assert(0);
    }
    return 0;
}

double BTag(WWJetType type, unsigned int iJet, float corjetpt){
    switch ( type ) {
        case jptJet:
            return BTag(cms2.jpts_p4().at(iJet));
            break;
        case CaloJet:
            return cms2.jets_trackCountingHighEffBJetTag()[iJet];
            break;
        case pfJet:
            return  cms2.pfjets_trackCountingHighEffBJetTag()[iJet]; 
            break;
        default:
            std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
            assert(0);
    }
    return 0;
}

bool defaultBTag(WWJetType type, unsigned int iJet, float jec) {
    
	switch ( type ) {
        case pfJet:
			if ( cms2.pfjets_trackCountingHighEffBJetTag()[iJet] > 2.1) return true;
            break;
        default:
            std::cout << "ERROR: Please use PFJets : " << type << " FixIt!" << std::endl;
            assert(0);
    }
    return 0;
}

