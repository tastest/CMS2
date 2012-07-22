#include "analysisSelections.h"

#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "Math/VectorUtil.h"
#include "TROOT.h"
#include "DYMVA.h"

using namespace std;

#ifndef __CINT__
#include "CORE/CMS2.h"
#include "CORE/electronSelections.h"
#include "CORE/muonSelections.h"
#include "CORE/jetSelections.h"
#include "CORE/metSelections.h"
#include "../Tools/MuonEffectiveArea.h"
#endif

//
// Electron ID
//
/*
bool goodElectronTMVA(ElectronIDMVA *mva, int useMVAeleId, unsigned int i) {
  //cout << "electronIdMVA.MVAValue=" << electronIdMVA->MVAValue(i, 0) << endl;
  
  //Find MVA Bin
  float pt = cms2.els_p4().at(i).pt();
  float etaSC = cms2.els_etaSC().at(i);
  if (useMVAeleId==2) {
    //preselection
    if (fabs(etaSC)<1.479) {
      if (cms2.els_sigmaIEtaIEta().at(i)>0.01 || 
	  fabs(cms2.els_dEtaIn().at(i))>0.007 ||
	  fabs(cms2.els_dPhiIn().at(i))>0.15 ||
	  cms2.els_hOverE().at(i)>0.12 ||
	  cms2.els_tkIso().at(i)/pt>0.2 ||
	  TMath::Max(cms2.els_ecalIso().at(i) - 1.0, 0.0)/pt>0.20 ||
	  //cms2.els_ecalIso().at(i)/pt>0.20 ||//fixme
	  cms2.els_hcalIso().at(i)/pt>0.20 ) return 0;
    } else {
      if (cms2.els_sigmaIEtaIEta().at(i)>0.03 || 
	  fabs(cms2.els_dEtaIn().at(i))>0.009 ||
	  fabs(cms2.els_dPhiIn().at(i))>0.10 ||
	  cms2.els_hOverE().at(i)>0.10 ||
	  cms2.els_tkIso().at(i)/pt>0.2 ||
	  //TMath::Max(cms2.els_ecalIso().at(i) - 1.0, 0.0)/pt>0.20 ||
	  cms2.els_ecalIso().at(i)/pt>0.20 ||
	  cms2.els_hcalIso().at(i)/pt>0.20 ) return 0;
    }
    //selection
    if (pt<20){
      if (fabs(etaSC)<1. && mva->MVAValue(i, 0)>0.139) return 1;
      else if (fabs(etaSC)>=1. && fabs(etaSC)<1.479 && mva->MVAValue(i, 0)>0.525) return 1;
      //else if (fabs(etaSC)>=1.479 && fabs(etaSC)<2.5 && mva->MVAValue(i, 0)>0.543) return 1;
      else if (fabs(etaSC)>=1.479 && mva->MVAValue(i, 0)>0.543) return 1;
    } else {
      if (fabs(etaSC)<1. && mva->MVAValue(i, 0)>0.947) return 1;
      else if (fabs(etaSC)>=1. && fabs(etaSC)<1.479 && mva->MVAValue(i, 0)>0.950) return 1;
      //else if (fabs(etaSC)>=1.479 && fabs(etaSC)<2.5 && mva->MVAValue(i, 0)>0.884) return 1;
      else if (fabs(etaSC)>=1.479 && mva->MVAValue(i, 0)>0.884) return 1;
    }
    return 0;

  } else if (useMVAeleId==3) {    
    int subdet = 0;
    if (fabs(etaSC) < 1.0) subdet = 0;
    else if (fabs(etaSC) < 1.479) subdet = 1;
    else subdet = 2;
    int ptBin = 0;
    if (pt > 20.0) ptBin = 1;
    int MVABin = -1;
    if (subdet == 0 && ptBin == 0) MVABin = 0;
    if (subdet == 1 && ptBin == 0) MVABin = 1;
    if (subdet == 2 && ptBin == 0) MVABin = 2;
    if (subdet == 0 && ptBin == 1) MVABin = 3;
    if (subdet == 1 && ptBin == 1) MVABin = 4;
    if (subdet == 2 && ptBin == 1) MVABin = 5;  
    double MVACut = -999.;
    double mvaValue=mva->MVAValue(i, 0);
    //WP with same Eff as LP2011 Cut based
    if (MVABin == 0) MVACut = 0.4202;
    if (MVABin == 1) MVACut = 0.6206;
    if (MVABin == 2) MVACut = 0.619; 
    if (MVABin == 3) MVACut = 0.959;
    if (MVABin == 4) MVACut = 0.9586;
    if (MVABin == 5) MVACut = 0.9278;
    //Explicitly Apply V4 Denominator Cuts
    bool pass = true;
    if (fabs(cms2.els_p4().at(i).eta()) >= 2.5) pass = false;
    //Barrel 
    if (fabs(etaSC)<1.479) {
      if (! ( (0==0)
	      && cms2.els_sigmaIEtaIEta().at(i)<0.01  
	      && fabs(cms2.els_dEtaIn().at(i))<0.007 
	      && fabs(cms2.els_dPhiIn().at(i))<0.15 
	      && cms2.els_hOverE().at(i)<0.12 
	      && cms2.els_tkIso().at(i)/pt<0.2 
	      && TMath::Max(cms2.els_ecalIso().at(i) - 1.0, 0.0)/pt<0.20 
	      && cms2.els_hcalIso().at(i)/pt<0.20 
	      //&& ele->nExpHitsInner <= 0
	      //&& passConversionVeto(ele->isConv)
	      //&& fabs(ele->d0) < 0.02
	      //&& fabs(ele->dz) < 0.1
	      && mvaValue > MVACut
	      )
	  ) {
	pass = false;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
	       && cms2.els_sigmaIEtaIEta().at(i)<0.03  
	       && fabs(cms2.els_dEtaIn().at(i))<0.009 
	       && fabs(cms2.els_dPhiIn().at(i))<0.10 
	       && cms2.els_hOverE().at(i)<0.10 
	       && cms2.els_tkIso().at(i)/pt<0.2 
	       && TMath::Max(cms2.els_ecalIso().at(i) - 1.0, 0.0)/pt<0.20 
	       && cms2.els_hcalIso().at(i)/pt<0.20 
	       //&& ele->nExpHitsInner <= 0
	       //&& passConversionVeto(ele->isConv)
	       //&& fabs(ele->d0) < 0.02
	       //&& fabs(ele->dz) < 0.1
	       && mvaValue > MVACut            
	       )
	  ) {
	pass = false;
      }
    } 
    
    return pass;
  } else {
    return false;
  }
}
*/

bool goodElectronWithoutIsolation(unsigned int i,  bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator){
    return ww_elBase(i) && ww_elId(i, useLHeleId, useMVAeleId, egammaMvaEleEstimator) && ww_eld0PV(i) && ww_eldZPV(i);
}

bool goodElectronIsolated(unsigned int i,  bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator, bool lockToCoreSelectors){
    bool ptcut = cms2.els_p4().at(i).pt() >= 10.0;
    bool core = ptcut && pass_electronSelection( i, electronSelection_smurfV6);
    bool internal = ww_elBase(i) && ww_elId(i, useLHeleId, useMVAeleId, egammaMvaEleEstimator) && ww_eld0PV(i) && ww_eldZPV(i) && ww_elIso(i);
    assert(!lockToCoreSelectors || core==internal);
    return internal;
}

bool ElectronFOIdV4(unsigned int i) {

	float pt = cms2.els_p4().at(i).pt();
	float etaSC = cms2.els_etaSC().at(i);

	if (fabs(etaSC)<1.479) {
		if (cms2.els_sigmaIEtaIEta().at(i)>0.01		||
			fabs(cms2.els_dEtaIn().at(i))>0.007 	||
			fabs(cms2.els_dPhiIn().at(i))>0.15 		||
			cms2.els_hOverE().at(i)>0.12 			||
			cms2.els_tkIso().at(i)/pt>0.2 			||
			(cms2.els_ecalIso().at(i) - 1.0)/pt>0.2 ||
			cms2.els_hcalIso().at(i)/pt>0.2 ) return false;
	} else {
		if (cms2.els_sigmaIEtaIEta().at(i)>0.03		|| 
			fabs(cms2.els_dEtaIn().at(i))>0.009 	||
			fabs(cms2.els_dPhiIn().at(i))>0.10 		|| 
			cms2.els_hOverE().at(i)>0.10 			||
			cms2.els_tkIso().at(i)/pt>0.2 			||
			cms2.els_ecalIso().at(i)/pt>0.2 		||
			cms2.els_hcalIso().at(i)/pt>0.2 ) return false;
	}

    // MIT conversion
	if ( isFromConversionMIT(i) ) return false;
	// conversion rejection - hit based
	if ( cms2.els_exp_innerlayers().at(i) > 0 ) return false;
	
	return true;
} 

bool ElectronFOV4(unsigned int i){
    return ww_elBase(i) && ElectronFOIdV4(i) && ww_eld0PV(i) && ww_eldZPV(i);
}

bool fakableElectron(unsigned int i, EleFOTypes type){
    if ( cms2.els_p4().at(i).pt() < 10.0 ) return false;
    switch (type){
        case EleFOV1: return pass_electronSelection( i, electronSelectionFO_el_smurf_v1);
        case EleFOV2: return pass_electronSelection( i, electronSelectionFO_el_smurf_v2);
        case EleFOV3: return pass_electronSelection( i, electronSelectionFO_el_smurf_v3);
        //case EleFOV4: return pass_electronSelection( i, electronSelectionFO_el_smurf_v4);
        case EleFOV4: return ElectronFOV4(i);
    }
    return false;
}

//
// Muon ID
//

bool goodMuonTMVA(MuonIDMVA* mva, unsigned int i) {
  //Find MVA Bin
  int subdet = 0;
  if (fabs(cms2.mus_p4().at(i).eta()) < 1.479) subdet = 0;
  else subdet = 1;
  int ptBin = 0;
  if (cms2.mus_p4().at(i).pt() > 14.5) ptBin = 1;
  if (cms2.mus_p4().at(i).pt() > 20.0) ptBin = 2;

  int MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;
  if (subdet == 0 && ptBin == 2) MVABin = 4;
  if (subdet == 1 && ptBin == 2) MVABin = 5;

  double MVACut = -999.;
  //same signal eff as cut-based (using V10 - Detector Based Iso)
  if (MVABin == 0) MVACut = -0.5618;
  if (MVABin == 1) MVACut = -0.3002;
  if (MVABin == 2) MVACut = -0.4642;
  if (MVABin == 3) MVACut = -0.2478;
  if (MVABin == 4) MVACut = 0.1706;
  if (MVABin == 5) MVACut = 0.8146;

  double mvaValue=mva->MVAValue(i, 0);

  //Isolation
  double iso03 = 0;
  iso03 = muonIsoValuePF(i,0,0.3);

  //Explicitly Apply M2 Denominator Cuts
  bool pass = true;
  if (cms2.mus_p4().at(i).pt() < 10) pass = false;
  if (fabs(cms2.mus_p4().at(i).eta()) >= 2.4) pass = false;

  if (! ( (0==0)
         &&
          (
           (((cms2.mus_type().at(i)) & (1<<1)) == (1<<1) 
            && cms2.mus_gfit_chi2().at(i)/cms2.mus_gfit_ndof().at(i) < 10.0
            && (cms2.mus_gfit_validSTAHits().at(i) > 0)
            && (cms2.mus_nmatches().at(i) > 1 )
            )
           || 
           ( ((cms2.mus_type().at(i)) & (1<<2)) == (1<<2)   
             && cms2.mus_pid_TMLastStationTight().at(i) == 1
             )
           )
          && ((cms2.mus_type().at(i)) & (1<<2)) == (1<<2)
          && cms2.mus_validHits().at(i) > 10
          && (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(i)) > 0)          
          //&& fabs(mu->d0) < 0.2
          //&& fabs(mu->dz) < 0.1
          && iso03 < 0.4
          && ( cms2.mus_ptErr().at(i)/cms2.mus_p4().at(i).pt() < 0.1)
          && cms2.mus_trkKink().at(i) < 20.
          && mvaValue > MVACut
          )
      ) {
    pass = false;
  }
  return pass;
}

bool goodMuonWithoutIsolation(unsigned int i, bool useMVAmuId, MuonIDMVA *mva, 
								MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle){
  return ww_muBase(i) && ww_mud0PV(i) && ww_mudZPV(i) && ww_muId(i, useMVAmuId, mva) && passMuonRingsMVAFO(i, muonMVAEstimator, IdentifiedMu, IdentifiedEle);
}

bool goodMuonIsolated(	unsigned int i, bool lockToCoreSelectors, bool useMVAmuId, MuonIDMVA *mva, 
						MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle ){
    bool ptcut = cms2.mus_p4().at(i).pt() >= 10.0;
    bool core = ptcut && muonId(i, NominalSmurfV6);
    bool internal = ww_muBase(i) && ww_mud0PV(i) && ww_mudZPV(i) && ww_muId(i, useMVAmuId, mva) && ww_muIso(i, muonMVAEstimator, IdentifiedMu,  IdentifiedEle); 
    assert(!lockToCoreSelectors || core==internal);
    return internal;
}


//
// Electron Id
//

bool ww_elBase(unsigned int index){
    if (cms2.els_p4().at(index).pt() < 10.0) return false;
    if (fabs(cms2.els_p4().at(index).eta()) > 2.5) return false;
    return true;
}
bool ww_elId(unsigned int index, bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator) {

    if (useLHeleId) {
        if (cms2.els_p4().at(index).pt()>20 && (passLikelihoodId_v2(index,cms2.els_lh().at(index),0) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false; 
        if (cms2.els_p4().at(index).pt()<20 && (passLikelihoodId_v2(index,cms2.els_lh().at(index),0) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false;
    }
    if (useMVAeleId>0){
      if (!goodElectronTMVA(egammaMvaEleEstimator, useMVAeleId, index)) return false;
    } else {
        if (!pass_electronSelection(index, electronSelection_smurfV3_id, false, false) ) return false;
    }

    // MIT conversion
    if ( isFromConversionMIT(index) ) return false;
    // conversion rejection - hit based
    if ( cms2.els_exp_innerlayers().at(index) > 0 ) return false;

    return true;
}

bool ww_eld0(unsigned int index){
    return fabs(cms2.els_d0corr()[index]) < 0.02;
}


double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
    return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
}

bool ww_eld0PV(unsigned int index){
    int vtxIndex = primaryVertex();
    if (vtxIndex<0) return false;
    double dxyPV = cms2.els_d0()[index]-
        cms2.vtxs_position()[vtxIndex].x()*sin(cms2.els_trk_p4()[index].phi())+
        cms2.vtxs_position()[vtxIndex].y()*cos(cms2.els_trk_p4()[index].phi());
    return fabs(dxyPV) < 0.02;
}

bool ww_eldZPV(unsigned int index){
    int vtxIndex = primaryVertex();
    if (vtxIndex<0) return false;
    // double dzPV = cms2.els_z0corr()[index]-cms2.vtxs_position()[iMax].z();
    double dzpv = dzPV(cms2.els_vertex_p4()[index], cms2.els_trk_p4()[index], cms2.vtxs_position()[vtxIndex]);
    return fabs(dzpv)<0.1;
}

double ww_elIsoVal(unsigned int index){
	return electronIsoValuePF2012_FastJetEffArea_HWW( index );
}

bool ww_elIso(unsigned int index){
	float pfiso = ww_elIsoVal( index ); 
	return pfiso<0.15;
}


//
// Muon Id
//

bool ww_muBase(unsigned int index){
    if (cms2.mus_p4().at(index).pt() < 10.0) return false;
    if (fabs(cms2.mus_p4().at(index).eta()) > 2.4) return false;
    if (cms2.mus_type().at(index) == 8) return false; // not STA
    return true;
}
bool ww_mud0(unsigned int index){
    return fabs(cms2.mus_d0corr()[index]) < 0.02;
}
double ww_mud0ValuePV(unsigned int index){
    int vtxIndex = primaryVertex();
    if (vtxIndex<0) return 9999;
    double dxyPV = cms2.mus_d0()[index]-
        cms2.vtxs_position()[vtxIndex].x()*sin(cms2.mus_trk_p4()[index].phi())+
        cms2.vtxs_position()[vtxIndex].y()*cos(cms2.mus_trk_p4()[index].phi());
    return fabs(dxyPV);
}

bool ww_mud0PV(unsigned int index){
    if ( cms2.mus_p4().at(index).pt() < 20. ) return ww_mud0ValuePV(index) < 0.01;
    return ww_mud0ValuePV(index) < 0.02;
}

bool ww_mudZPV(unsigned int index, float cut){
    int vtxIndex = primaryVertex();
    if (vtxIndex<0) return false;
    // double dzpv = cms2.mus_z0corr()[index]-cms2.vtxs_position()[iMax].z();
    double dzpv = dzPV(cms2.mus_vertex_p4()[index], cms2.mus_trk_p4()[index], cms2.vtxs_position()[vtxIndex]);
    return fabs(dzpv)<cut;
}

bool ww_muId(unsigned int index, bool useMVAmuId, MuonIDMVA *mva){ 
    if (useMVAmuId){
      if (!goodMuonTMVA(mva,index)) return false;
      return true;
    }
    
	if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
    if (cms2.trks_nlayers().at(cms2.mus_trkidx().at(index)) < 6) return false; // # of tracker hits 
    if (cms2.mus_ptErr().at(index)/cms2.mus_trk_p4().at(index).pt()>0.1) return false; // Does pt come from track?
    if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(index))==0) return false;
    if (cms2.mus_trkKink().at(index) > 20.) return false; //kink finder
    if (!cms2.mus_pid_PFMuon().at(index)) return false; // should be a pfmuon
    // global muon
    bool goodMuonGlobalMuon = false;
    if (((cms2.mus_type().at(index)) & (1<<1)) == (1<<1)){
        goodMuonGlobalMuon = true;
        if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) goodMuonGlobalMuon = false; //glb fit chisq
        if (cms2.mus_gfit_validSTAHits().at(index)==0 ) goodMuonGlobalMuon = false;
        if (cms2.mus_nmatches().at(index)<2) goodMuonGlobalMuon = false;
    }
    return goodMuonGlobalMuon || 
        cms2.mus_pid_TMLastStationTight().at(index) == 1; // TM id
}

double ww_muIsoVal(unsigned int index){
    double sum =  cms2.mus_iso03_sumPt().at(index) +
        cms2.mus_iso03_emEt().at(index)  +
        cms2.mus_iso03_hadEt().at(index);
    double pt  = cms2.mus_p4().at(index).pt();
    return sum/pt;
}

bool ww_muIso(unsigned int index){
    if (cms2.mus_p4().at(index).pt()>20) {
        if (TMath::Abs(cms2.mus_p4()[index].eta())<1.479) 
            return muonIsoValuePF(index,0,0.3) < 0.13;
        else 
            return muonIsoValuePF(index,0,0.3) < 0.09;
    } else {
        if (TMath::Abs(cms2.mus_p4()[index].eta())<1.479) 
            return muonIsoValuePF(index,0,0.3) < 0.06;
        else 
            return muonIsoValuePF(index,0,0.3) < 0.05;
    }
}

bool ww_muIso(unsigned int index, MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle){
	return passMuonRingsMVA(index, muonMVAEstimator, IdentifiedMu, IdentifiedEle); 
}


/////
unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated,
        const std::vector<JetPair>& vetojets)
{
    unsigned int nMuons = 0;
    for (int imu=0; imu < int(cms2.mus_charge().size()); ++imu) {
        // quality cuts
        if (  ((cms2.mus_goodmask()[imu]) & (1<<19)) == 0 ) continue; // TMLastStationAngTight
        if ( cms2.mus_p4()[imu].pt() < 3 ) continue;
        if ( ww_mud0ValuePV(imu) > 0.2) continue;
        if ( ! ww_mudZPV(imu,0.2) ) continue; //newcuts, was 0.1
        //if ( cms2.mus_validHits()[imu] < 11) continue;
        if (cms2.trks_nlayers().at(cms2.mus_trkidx().at(imu)) < 6) return false; // # of tracker hits 
        if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == imu ) continue;
        if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == imu ) continue;
        if ( nonisolated && ww_muIsoVal(imu)<0.1 && cms2.mus_p4()[imu].pt()>20 ) continue;
        bool skip = false;
        for ( std::vector<JetPair>::const_iterator ijet = vetojets.begin();
                ijet != vetojets.end(); ++ijet )
            if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ijet->first,cms2.mus_p4()[imu])) < 0.3 ) skip=true;
        if ( skip ) continue;
        ++nMuons;
    }
    return nMuons;
}

std::vector<LeptonPair> getExtraLeptons(int i_hyp, double minPt,  bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator, bool useMVAmuId, MuonIDMVA *mumva,
										MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle ){

	std::vector<LeptonPair> leptons;
    for (int i=0; i < int(cms2.mus_charge().size()); ++i) {
        if ( cms2.mus_p4()[i].pt() < minPt ) continue;
        if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == i ) continue;
        if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == i ) continue;
        //if ( ! (ww_mud0PV(i) && ww_muId(i, useMVAmuId, mumva) && ww_muIso(i)&&
        if ( ! (ww_mud0PV(i) && ww_muId(i, useMVAmuId, mumva) && ww_muIso(i, muonMVAEstimator, IdentifiedMu,  IdentifiedEle) &&
                    fabs(cms2.mus_p4().at(i).eta()) <2.4) ) continue;
        leptons.push_back(LeptonPair(true,i));
    }
    for (int i=0; i < int(cms2.els_charge().size()); ++i) {
        if ( cms2.els_p4()[i].pt() < minPt ) continue;
        if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.els_p4().at(i)) <0.1) ) continue;
        if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.els_p4().at(i)) <0.1) ) continue;
        if ( !(ww_elId(i, useLHeleId, useMVAeleId, egammaMvaEleEstimator) && ww_eld0PV(i) && ww_elIso(i) &&
                    fabs(cms2.els_p4().at(i).eta()) < 2.5) ) continue;
        leptons.push_back(LeptonPair(false,i));
    }
    return leptons;
}

unsigned int numberOfExtraLeptons(	int i_hyp, double minPt, bool useLHeleId, int useMVAeleId, EGammaMvaEleEstimator* egammaMvaEleEstimator, bool useMVAmuId, MuonIDMVA *mumva,
									MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle){
	return getExtraLeptons(i_hyp, minPt, useLHeleId, useMVAeleId, egammaMvaEleEstimator, useMVAmuId, mumva,muonMVAEstimator,IdentifiedMu,IdentifiedEle).size();
}


//
// Triger
//

bool passedTriggerRegExp(TPMERegexp trigName)
{

    for (unsigned int tidx = 0; tidx < cms2.hlt_trigNames().size(); tidx++) {
        if (trigName.Match(cms2.hlt_trigNames().at(tidx)) != 0) {
            return true;
        }
    }

    return false;

}

bool passedTrigger(TString trigName, unsigned int minRun, unsigned int maxRun) 
{
  if ( cms2.evt_run() < minRun || cms2.evt_run() > maxRun ) return false;
  if ( find(cms2.hlt_trigNames().begin(), cms2.hlt_trigNames().end(), trigName)
       == cms2.hlt_trigNames().end() ) return false;
  return cms2.passHLTTrigger(trigName);
}

bool passedTriggerRequirements() {
    
	// 
	//  2012 ICHEP triggers
	// 

    // double electron 
    if(  	passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15")   ||
     		passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16") 	|| 
     		passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17")  
	  ) return true;

    // double muon 
	if(  	passedTrigger("HLT_Mu17_Mu8_v16") 	||
			passedTrigger("HLT_Mu17_Mu8_v17") 	||
			passedTrigger("HLT_Mu17_Mu8_v18") 	||
			passedTrigger("HLT_Mu17_TkMu8_v9") 	||
			passedTrigger("HLT_Mu17_TkMu8_v10")	||
			passedTrigger("HLT_Mu17_TkMu8_v11")
	  ) return true;

	// emu
	if(  	passedTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4") 	||
	  		passedTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5") 	||
	  		passedTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6") 	||
	  		passedTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7") 	||
	  		passedTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4") 	||
	  		passedTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5") 	||
	  		passedTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6") 	||
	  		passedTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7") 	
	  ) return true;

	// single ele
	if(  	passedTrigger("HLT_Ele27_WP80_v8") 	||
	  		passedTrigger("HLT_Ele27_WP80_v9") 	||
	  		passedTrigger("HLT_Ele27_WP80_v10") 	
	  ) return true;

	// single ele
	if(  	passedTrigger("HLT_IsoMu24_eta2p1_v11") 	||
	 		passedTrigger("HLT_IsoMu24_eta2p1_v12") 	||
	 		passedTrigger("HLT_IsoMu24_eta2p1_v13") 	
	  ) return true;
    
    return false;
}

bool passedTriggerRequirementsWithRuns() {
  // return true; // no trigger requirements
  // return cms2.filter_ele10mu10IsoId_passed();
  if ( passedTrigger("HLT_Mu17_Ele8_CaloIdL_v1",1,175972) || 
       passedTrigger("HLT_Mu17_Ele8_CaloIdL_v2",1,175972) ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdL_v3",1,175972) ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdL_v4",1,175972) ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdL_v5",1,175972) ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdL_v6",1,175972) ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdL_v8",1,175972) ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v1",175973,999999) ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v3",175973,999999) ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4",175973,999999) ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7",175973,999999) ||
       passedTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8",175973,999999) ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdL_v1",1,167913) ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdL_v2",1,167913) ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdL_v3",1,167913) ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdL_v4",1,167913) ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdL_v5",1,167913) || 
       passedTrigger("HLT_Mu8_Ele17_CaloIdL_v6",1,167913) ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v1",167914,999999) ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3",167914,999999) ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4",167914,999999) ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7",167914,999999) ||
       passedTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8",167914,999999) 
      ) return true;
  if ( passedTrigger("HLT_DoubleMu7_v1",1,164237) ||
       passedTrigger("HLT_DoubleMu7_v2",1,164237) ||
       passedTrigger("HLT_Mu13_Mu8_v2",165085,178419) ||
       passedTrigger("HLT_Mu13_Mu8_v3",165085,178419) ||
       passedTrigger("HLT_Mu13_Mu8_v4",165085,178419) ||
       passedTrigger("HLT_Mu13_Mu8_v6",165085,178419) ||
       passedTrigger("HLT_Mu13_Mu8_v7",165085,178419) ||
       passedTrigger("HLT_Mu17_Mu8_v2",178420,999999) ||
       passedTrigger("HLT_Mu17_Mu8_v3",178420,999999) ||
       passedTrigger("HLT_Mu17_Mu8_v4",178420,999999) ||
       passedTrigger("HLT_Mu17_Mu8_v6",178420,999999) ||
       passedTrigger("HLT_Mu17_Mu8_v7",178420,999999) ||
       passedTrigger("HLT_Mu17_Mu8_v10",178420,999999) ||
       passedTrigger("HLT_Mu17_Mu8_v11",178420,999999) ||
       passedTrigger("HLT_Mu17_TkMu8_v3") ||
       passedTrigger("HLT_Mu17_TkMu8_v4") ) return true;
  if ( passedTrigger("HLT_Mu15_v2",1,163261) ||
       passedTrigger("HLT_Mu24_v1",163262,164237) ||
       passedTrigger("HLT_Mu24_v2",163262,164237) ||
       passedTrigger("HLT_Mu30_v1",165085,173235) ||
       passedTrigger("HLT_Mu30_v2",165085,173235) ||
       passedTrigger("HLT_Mu30_v3",165085,173235) ||
       passedTrigger("HLT_Mu30_v4",165085,173235) ||
       passedTrigger("HLT_Mu30_v5",165085,173235) ||
       passedTrigger("HLT_Mu30_v7",165085,173235) ||
       passedTrigger("HLT_Mu40_v5",173236,175972) ||
       passedTrigger("HLT_Mu40_eta2p1_v1",175973,999999) ||
       passedTrigger("HLT_Mu40_eta2p1_v4",175973,999999) ||
       passedTrigger("HLT_Mu40_eta2p1_v5",175973,999999) ||
       passedTrigger("HLT_IsoMu17_v5",163262,167043) ||
       passedTrigger("HLT_IsoMu17_v6",163262,167043) ||
       passedTrigger("HLT_IsoMu17_v8",163262,167043) ||
       passedTrigger("HLT_IsoMu17_v9",163262,167043) ||
       passedTrigger("HLT_IsoMu17_v10",163262,167043) ||
       passedTrigger("HLT_IsoMu17_v11",163262,167043) ||
       passedTrigger("HLT_IsoMu17_eta2p1_v1",167044,167913) ||
       passedTrigger("HLT_IsoMu20_v8",170053,175910) ||
       passedTrigger("HLT_IsoMu24_v1",175911,175921) ||
       passedTrigger("HLT_IsoMu24_v2",175911,175921) ||
       passedTrigger("HLT_IsoMu24_v4",175911,175921) ||
       passedTrigger("HLT_IsoMu24_v5",175911,175921) ||
       passedTrigger("HLT_IsoMu24_v6",175911,175921) ||
       passedTrigger("HLT_IsoMu24_v7",175911,175921) ||
       passedTrigger("HLT_IsoMu24_v8",175911,175921) ||
       passedTrigger("HLT_IsoMu24_eta2p1_v3",175922,176544) ||
       passedTrigger("HLT_IsoMu24_eta2p1_v6",175922,176544) ||
       passedTrigger("HLT_IsoMu24_eta2p1_v7",175922,176544) ||
       passedTrigger("HLT_IsoMu30_eta2p1_v3",176545,999999) ||
       passedTrigger("HLT_IsoMu30_eta2p1_v6",176545,999999) ||
       passedTrigger("HLT_IsoMu30_eta2p1_v7",176545,999999) ) return true;
  if ( passedTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",1,170052) ||
       passedTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",1,170052) ||
       passedTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",1,170052) ||
       passedTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4",1,170052) ||
       passedTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5",1,170052) ||
       passedTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6",1,170052) ||
       passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5",170053,999999) ||
       passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",170053,999999) ||
       passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",170053,999999) ||
       passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",170053,999999) ||
       passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",170053,999999) ||
       passedTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",170053,999999) ) return true;
  if ( passedTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1",1,164237) ||
       passedTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",1,164237) ||
       passedTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",1,164237) ||
       passedTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1",165085,166967) ||
       passedTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",165085,166967) ||
       passedTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",165085,166967) ||
       passedTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4",165085,166967) ||
       passedTrigger("HLT_Ele52_CaloIdVT_TrkIdT_v1",166968,170901) ||
       passedTrigger("HLT_Ele52_CaloIdVT_TrkIdT_v2",166968,170901) ||
       passedTrigger("HLT_Ele52_CaloIdVT_TrkIdT_v3",166968,170901) ||
       passedTrigger("HLT_Ele65_CaloIdVT_TrkIdT_v3",170902,178419) ||
       passedTrigger("HLT_Ele65_CaloIdVT_TrkIdT_v4",170902,178419) ||
       passedTrigger("HLT_Ele80_CaloIdVT_TrkIdT_v2",178420,999999) ||
       passedTrigger("HLT_Ele80_CaloIdVT_TrkIdT_v3",178420,999999) ) return true;
  return false;

}

//
// MET
//

double minmet(unsigned int i_hyp) {
    metStruct trkMET = trackerMET(i_hyp,0.1); 
    double pMet = std::min(projectedMet(i_hyp, metValue(), metPhiValue()),
            projectedMet(i_hyp, trkMET.met, trkMET.metphi));
    return pMet;
}

bool passedMetRequirements(unsigned int i_hyp, FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3) {
   
    bool applyJEC = true;
    int njets = numberOfJets(i_hyp, applyJEC, jet_corrector_pfL1FastJetL2L3); 
	
	return passedMetRequirements(i_hyp, njets, getJets( jetType(), i_hyp, 0., 4.7, applyJEC, jet_corrector_pfL1FastJetL2L3, true, false));
}

bool passedMetRequirements(unsigned int i_hyp, unsigned int njets, std::vector<JetPair> jets ){
	float dymva=-999.;
    HypothesisType type = getHypothesisTypeNew(i_hyp);
    metStruct trkMET = trackerMET(i_hyp,0.1); //,&jets);
    double pMet = std::min(projectedMet(i_hyp, metValue(), metPhiValue()),
            projectedMet(i_hyp, trkMET.met, trkMET.metphi));
    if ( njets<2 && pMet < 20 ) return false;
    if (type == EE || type == MM) {
		if(njets==0) {
			dymva =  DYMVA(i_hyp, njets, jets );
			if( dymva < 0.6 ) return false;
		}
		else if(njets==1) {
			dymva =  DYMVA(i_hyp, njets, jets ); 
			if( dymva < -0.01 ) return false;
		}
		else {
        	double threshold = 40 + nGoodVertex()/2.0;
        	if ( metValue() < threshold ) return false;
		}
    }
	return true;
}

double nearestDeltaPhi(double Phi, int i_hyp)
{
    double tightDPhi = fabs(cms2.hyp_lt_p4()[i_hyp].Phi() - Phi);
    tightDPhi = std::min(2*TMath::Pi() - tightDPhi, tightDPhi);
    double looseDPhi = fabs(cms2.hyp_ll_p4()[i_hyp].Phi() - Phi);
    looseDPhi = std::min(2*TMath::Pi() - looseDPhi, looseDPhi);
    return TMath::Min(tightDPhi, looseDPhi);
}

double projectedMet(unsigned int i_hyp, double met, double phi)
{
    double DeltaPhi = nearestDeltaPhi(phi,i_hyp);
    if (DeltaPhi < TMath::Pi()/2) return met*TMath::Sin(DeltaPhi);
    return met;
}

//
// Vertex selections
//

bool isGoodVertex(size_t ivtx) {
    if (cms2.vtxs_isFake()[ivtx]) return false;
    if (cms2.vtxs_ndof()[ivtx] < 4.) return false;
    if (cms2.vtxs_position()[ivtx].Rho() > 2.0) return false;
    if (fabs(cms2.vtxs_position()[ivtx].Z()) > 24.0) return false;
    return true;
}

unsigned int nGoodVertex() {
    unsigned int nVtx = 0;
    for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
        // if (cms2.vtxs_isFake()[i]) continue;
        if (!isGoodVertex(i)) continue;
        nVtx++;
    }
    return nVtx;
}

int primaryVertex(){
    //  double sumPtMax = -1;
    //   int iMax = -1;
    //   for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
    //     // if (cms2.vtxs_isFake()[i]) continue;
    //     if (!isGoodVertex(i)) continue;
    //     if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
    //       iMax = i;
    //       sumPtMax = cms2.vtxs_sumpt().at(i);
    //     }
    //   }
    //   if (iMax<0) return false;
    return 0;
}

//
// other cuts
//

bool inZmassWindow(float mass, double delta){
    // return ( mass > 76. && mass < 106. );
    return fabs(mass - 91.1876) < delta;
}

//
// event top tagging
//

bool toptag(WWJetType type, int i_hyp, double minPt, FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3,
        std::vector<JetPair> ignoreJets)
{
    const double vetoCone    = 0.3;

    switch ( type ){
        case pfJet:
            for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
                if ( cms2.pfjets_p4()[i].pt() < minPt ) continue;
                bool ignoreJet = false;
                for ( std::vector<JetPair>::const_iterator ijet = ignoreJets.begin();
                        ijet != ignoreJets.end(); ++ijet )
                    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ijet->first,cms2.pfjets_p4()[i])) < vetoCone ) ignoreJet=true;
                if ( ignoreJet ) continue;
                if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ||
						TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ) continue;

				double jec = 1.0;
				jet_corrector_pfL1FastJetL2L3->setRho(cms2.evt_ww_rho()); 
				jet_corrector_pfL1FastJetL2L3->setJetA(cms2.pfjets_area().at(i));
				jet_corrector_pfL1FastJetL2L3->setJetPt(cms2.pfjets_p4()[i].pt());
				jet_corrector_pfL1FastJetL2L3->setJetEta(cms2.pfjets_p4()[i].eta());
				double corr = jet_corrector_pfL1FastJetL2L3->getCorrection();
				jec *= corr;

				if ( !passMVAJetId( cms2.pfjets_p4()[i].pt() * jec, cms2.pfjets_p4()[i].eta(), cms2.pfjets_mvavalue()[i], 2) ) continue;
               
				if ( !defaultBTag(type,i, jec) ) continue;
				
                return true;
            }
            break;
        case CaloJet:
            for ( unsigned int i=0; i < cms2.jets_p4().size(); ++i) {
                if ( cms2.jets_p4()[i].pt() < minPt ) continue;
                bool ignoreJet = false;
                for ( std::vector<JetPair>::const_iterator ijet = ignoreJets.begin();
                        ijet != ignoreJets.end(); ++ijet )
                    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(ijet->first,cms2.jets_p4()[i])) < vetoCone ) ignoreJet=true;
                if ( ignoreJet ) continue;
                if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jets_p4()[i])) < vetoCone ||
                        TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jets_p4()[i])) < vetoCone ) continue;
                //     if ( defaultBTag(type,i) && ignoreJets.size()==1 ){
                //       cout << "b-tagged jet pt: " << cms2.jets_p4()[i].pt() << " \teta: " << cms2.jets_p4()[i].eta() <<
                //         " \tphi: " << cms2.jets_p4()[i].phi() << endl;
                //     }
                if ( defaultBTag(type,i) ) return true;
            }
            break;
        default:
            std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
    }
    return false;
}

// tightness : 2=loose 1=medium 0=tight
bool passMVAJetId(double corjetpt, double jeteta, double mvavalue, unsigned int tightness)         
{
	if(tightness<0 || tightness>2)
	{
		cout << "ERROR : tightness should be 0, 1, or 2. " << endl;
		return false;
	}

	double fMVACut[3][4][4];
/*
	// original cuts for 52X (used for MVA Met)
	// Do not use these cuts for Jet Id in 52X
	//Tight Id
	fMVACut[0][0][0] =  0.5; fMVACut[0][0][1] = 0.6; fMVACut[0][0][2] = 0.6; fMVACut[0][0][3] = 0.9;
	fMVACut[0][1][0] = -0.2; fMVACut[0][1][1] = 0.2; fMVACut[0][1][2] = 0.2; fMVACut[0][1][3] = 0.6;
	fMVACut[0][2][0] =  0.3; fMVACut[0][2][1] = 0.4; fMVACut[0][2][2] = 0.7; fMVACut[0][2][3] = 0.8;
	fMVACut[0][3][0] =  0.5; fMVACut[0][3][1] = 0.4; fMVACut[0][3][2] = 0.8; fMVACut[0][3][3] = 0.9;
	//Medium id
	fMVACut[1][0][0] =  0.2; fMVACut[1][0][1] = 0.4; fMVACut[1][0][2] = 0.2; fMVACut[1][0][3] = 0.6;
	fMVACut[1][1][0] = -0.3; fMVACut[1][1][1] = 0. ; fMVACut[1][1][2] = 0. ; fMVACut[1][1][3] = 0.5;
	fMVACut[1][2][0] =  0.2; fMVACut[1][2][1] = 0.2; fMVACut[1][2][2] = 0.5; fMVACut[1][2][3] = 0.7;
	fMVACut[1][3][0] =  0.3; fMVACut[1][3][1] = 0.2; fMVACut[1][3][2] = 0.7; fMVACut[1][3][3] = 0.8;
	//Loose Id 
	fMVACut[2][0][0] = -0.2; fMVACut[2][0][1] =  0. ; fMVACut[2][0][2] =  0.2; fMVACut[2][0][3] =  0.5;
	fMVACut[2][1][0] =  0.2; fMVACut[2][1][1] = -0.6; fMVACut[2][1][2] = -0.6; fMVACut[2][1][3] = -0.4;
	fMVACut[2][2][0] =  0.2; fMVACut[2][2][1] = -0.6; fMVACut[2][2][2] = -0.6; fMVACut[2][2][3] = -0.4;
	fMVACut[2][3][0] =  0.2; fMVACut[2][3][1] = -0.8; fMVACut[2][3][2] = -0.8; fMVACut[2][3][3] = -0.4;
*/	

	// These cuts are for 42X but used for 52X jet Id
	//Tight Id
	fMVACut[0][0][0] =  0.5; fMVACut[0][0][1] = 0.6; fMVACut[0][0][2] = 0.6; fMVACut[0][0][3] = 0.9;
	fMVACut[0][1][0] = -0.2; fMVACut[0][1][1] = 0.2; fMVACut[0][1][2] = 0.2; fMVACut[0][1][3] = 0.6;
	fMVACut[0][2][0] =  0.3; fMVACut[0][2][1] = 0.4; fMVACut[0][2][2] = 0.7; fMVACut[0][2][3] = 0.8;
	fMVACut[0][3][0] =  0.5; fMVACut[0][3][1] = 0.4; fMVACut[0][3][2] = 0.8; fMVACut[0][3][3] = 0.9;
	//Medium id
	fMVACut[1][0][0] =  0.2; fMVACut[1][0][1] = 0.4; fMVACut[1][0][2] = 0.2; fMVACut[1][0][3] = 0.6;
	fMVACut[1][1][0] = -0.3; fMVACut[1][1][1] = 0. ; fMVACut[1][1][2] = 0. ; fMVACut[1][1][3] = 0.5;
	fMVACut[1][2][0] =  0.2; fMVACut[1][2][1] = 0.2; fMVACut[1][2][2] = 0.5; fMVACut[1][2][3] = 0.7;
	fMVACut[1][3][0] =  0.3; fMVACut[1][3][1] = 0.2; fMVACut[1][3][2] = 0.7; fMVACut[1][3][3] = 0.8;
	//Loose Id 
	fMVACut[2][0][0] =  0. ; fMVACut[2][0][1] =  0. ; fMVACut[2][0][2] =  0. ; fMVACut[2][0][3] = 0.2;
	fMVACut[2][1][0] = -0.4; fMVACut[2][1][1] = -0.4; fMVACut[2][1][2] = -0.4; fMVACut[2][1][3] = 0.4;
	fMVACut[2][2][0] =  0. ; fMVACut[2][2][1] =  0. ; fMVACut[2][2][2] =  0.2; fMVACut[2][2][3] = 0.6;
	fMVACut[2][3][0] =  0. ; fMVACut[2][3][1] =  0. ; fMVACut[2][3][2] =  0.6; fMVACut[2][3][3] = 0.2;


	// pT categorization
	int ptId = 0;
	if( corjetpt > 10 && corjetpt < 20 ) ptId = 1;
	if( corjetpt > 20 && corjetpt < 30 ) ptId = 2;
	if( corjetpt > 30                  ) ptId = 3;

	// eta categorization
	int etaId = 0;
	if( fabs(jeteta) > 2.5  && fabs(jeteta) < 2.75 ) etaId = 1;
	if( fabs(jeteta) > 2.75 && fabs(jeteta) < 3.0  ) etaId = 2;
	if( fabs(jeteta) > 3.0  && fabs(jeteta) < 5.0  ) etaId = 3;

	// return  
	if( mvavalue > fMVACut[tightness][ptId][etaId] ) return true;
	return false;
}

bool goodElectronTMVA(EGammaMvaEleEstimator* egammaMvaEleEstimator, int useMVAeleId, unsigned int i) 
{  

	float pt = cms2.els_p4().at(i).pt();
  	float etaSC = cms2.els_etaSC().at(i);

	//preselection
	if (fabs(etaSC)<1.479) {
		if (cms2.els_sigmaIEtaIEta().at(i)>0.01 || 
				fabs(cms2.els_dEtaIn().at(i))>0.007 ||
				fabs(cms2.els_dPhiIn().at(i))>0.15 ||
				cms2.els_hOverE().at(i)>0.12 ||
				cms2.els_tkIso().at(i)/pt>0.2 ||
				TMath::Max(cms2.els_ecalIso().at(i) - 1.0, 0.0)/pt>0.20 ||
				cms2.els_hcalIso().at(i)/pt>0.20 ) return 0;
	} else {
		if (cms2.els_sigmaIEtaIEta().at(i)>0.03 || 
				fabs(cms2.els_dEtaIn().at(i))>0.009 ||
				fabs(cms2.els_dPhiIn().at(i))>0.10 ||
				cms2.els_hOverE().at(i)>0.10 ||
				cms2.els_tkIso().at(i)/pt>0.2 ||
				cms2.els_ecalIso().at(i)/pt>0.20 ||
				cms2.els_hcalIso().at(i)/pt>0.20 ) return 0;
	}

	// MIT conversion
	if ( isFromConversionMIT(i) ) return false;
	// conversion rejection - hit based
	if ( cms2.els_exp_innerlayers().at(i) > 0 ) return false;

	double mvavalue =  egammaMvaEleEstimator->mvaValue(i,false);

	if( pt > 20 ) {
		if( fabs(etaSC)>=1.479 && mvavalue>0.92)  return true;
		if( fabs(etaSC)>=0.8 && fabs(etaSC)<1.479 && mvavalue>0.85)  return true;
		if( fabs(etaSC)<0.8 && mvavalue>0.94)  return true;
		return false;
	}
	else {
		if( fabs(etaSC)>=1.479 && mvavalue>0.62)  return true;
		if( fabs(etaSC)>=0.8 && fabs(etaSC)<1.479 && mvavalue>0.1)  return true;
		if( fabs(etaSC)<0.8 && 
			mvavalue>0.0)  return true;
		return false;
	}

	cout << "Something is wrong. You should not see this! " << endl; 
}

bool passMuonRingsMVA(unsigned int mu, MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle)
{
	double mvavalue = muonMVAEstimator->mvaValueIso( mu, cms2.evt_ww_rho(), MuonEffectiveArea::kMuEAFall11MC,
	                                               	 IdentifiedEle, IdentifiedMu, false );

	double pt 	= cms2.mus_trk_p4()[mu].pt();
	double eta 	= cms2.mus_trk_p4()[mu].eta();

	if( pt>20 ) {
		if( fabs(eta)>=1.479 && fabs(eta)<2.4 && mvavalue>0.86 )  return true;
		if( fabs(eta)<1.479 && mvavalue>0.82 )  return true;
		return false;
	}
	else {
		if( fabs(eta)>=1.479 && fabs(eta)<2.4 && mvavalue>0.82 )  return true;
		if( fabs(eta)<1.479 && mvavalue>0.86 )  return true;
		return false;
	}	
	
	cout << "Something is wrong. You should not see this! " << endl; 
}

bool passMuonRingsMVAFO(unsigned int mu, MuonMVAEstimator* muonMVAEstimator, std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle)
{
	double mvavalue = muonMVAEstimator->mvaValueIso( mu, cms2.evt_ww_rho(), MuonEffectiveArea::kMuEAFall11MC,
	                                               	 IdentifiedEle, IdentifiedMu, false );

	if( mvavalue>-0.6 )  return true;
	return false;
}

bool MuonFOV2(	unsigned int i, MuonMVAEstimator* muonMVAEstimator, 
				std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle){

	if (((cms2.mus_type().at(i)) & (1<<2)) == 0)    return false; // tracker muon
    if (cms2.trks_nlayers().at(cms2.mus_trkidx().at(i)) < 6) return false; // # of tracker hits 
    if (cms2.mus_ptErr().at(i)/cms2.mus_trk_p4().at(i).pt()>0.1) return false; // Does pt come from track?
    if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(i))==0) return false;
    if (cms2.mus_trkKink().at(i) > 20.) return false; //kink finder
    if (!cms2.mus_pid_PFMuon().at(i)) return false; // should be a pfmuon
    // global muon
    bool goodMuonGlobalMuon = false;
    if (((cms2.mus_type().at(i)) & (1<<1)) == (1<<1)) {
        goodMuonGlobalMuon = true;
        if (cms2.mus_gfit_chi2().at(i)/cms2.mus_gfit_ndof().at(i) >= 10) goodMuonGlobalMuon = false; //glb fit chisq
        if (cms2.mus_gfit_validSTAHits().at(i)==0 ) goodMuonGlobalMuon = false;
        if (cms2.mus_nmatches().at(i)<2) goodMuonGlobalMuon = false;
    }

    return 	(goodMuonGlobalMuon || cms2.mus_pid_TMLastStationTight().at(i) == 1) 	&& // ---> Id
			ww_muBase(i) 															&& 
			ww_mud0ValuePV(i)<0.2 													&& 
			ww_mudZPV(i) 															&& 
  			passMuonRingsMVAFO(i, muonMVAEstimator, IdentifiedMu, IdentifiedEle);
}


bool fakableMuon(unsigned int i, MuFOTypes type,MuonMVAEstimator* muonMVAEstimator,
                std::vector<Int_t> IdentifiedMu, std::vector<Int_t> IdentifiedEle){
    if ( cms2.mus_p4().at(i).pt() < 10.0 ) return false;
    switch (type){
        case MuFOV1: return muonId(i, muonSelectionFO_mu_smurf_10);
        //case MuFOV2: return muonId(i, muonSelectionFO_mu_smurf_04);
        case MuFOV2: return MuonFOV2(i, muonMVAEstimator, IdentifiedMu, IdentifiedEle);
    }
    return false;
}

