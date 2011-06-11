//
// electron selections
//

// CMS2 includes
#include "electronSelections.h"
#include "eventSelections.h"
//#include "MITConversionUtilities.cc"
#include "MITConversionUtilities.h"
#include "trackSelections.h"
#include "utilities.h"
#include "muonSelections.h"

bool pass_electronSelectionCompareMask(const cuts_t cuts_passed, const cuts_t selectionType)
{
    if ((cuts_passed & selectionType) == selectionType) return true;
    return false;
}

bool pass_electronSelection(const unsigned int index, const cuts_t selectionType, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, bool useGsfTrack, int vertex_index)
{
  checkElectronSelections();
  cuts_t cuts_passed = electronSelection(index, applyAlignmentCorrection, removedEtaCutInEndcap, vertex_index);
    if ((cuts_passed & selectionType) == selectionType) return true;
    return false;
}

cuts_t electronSelection(const unsigned int index, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, bool useGsfTrack, int vertex_index) 
{
    // keep track of which cuts passed
    cuts_t cuts_passed = 0;

    // isolation

    // relative isolation non truncated
    if( electronIsolation_rel_v1(index, true ) < 0.10) cuts_passed |= (1ll<<ELEISO_RELNT010);       // Relative Isolation
    if( electronIsolation_rel_v1(index, true ) < 0.15) cuts_passed |= (1ll<<ELEISO_RELNT015);       //
    if( electronIsolation_rel_v1(index, true ) < 0.40) cuts_passed |= (1ll<<ELEISO_RELNT040);       //
    if( electronIsolation_rel_v1(index, false) < 0.20) cuts_passed |= (1ll<<ELEISO_TRK_RELNT020);   // Tracker Relative Isolation
    if( electronIsolation_ECAL_rel_v1(index)   < 0.20) cuts_passed |= (1ll<<ELEISO_ECAL_RELNT020);  // ECAL    Relative Isolation
    if( electronIsolation_HCAL_rel_v1(index)   < 0.20) cuts_passed |= (1ll<<ELEISO_HCAL_RELNT020);  // HCAL    Relative Isolation
    if (electronIsolation_ECAL_rel_v1(index, false) < 0.20) cuts_passed |= (1ll<<ELEISO_ECAL_RELNT020_NPS); // ECAL Relative Isolation, no ped sub in EB
    if( electronIsolation_ECAL_rel(index)      < 0.20) cuts_passed |= (1ll<<ELEISO_ECAL_REL020);    // ECAL    Relative Isolation (truncated)
    if( electronIsolation_HCAL_rel(index)      < 0.20) cuts_passed |= (1ll<<ELEISO_HCAL_REL020);    // HCAL    Relative Isolation (truncated)

    //relative isolation truncated
    if (electronIsolation_rel_FastJet(index, true) < 0.05) cuts_passed |= (1ll<<ELEISO_FASTJET_REL005); // ADDED
    if (electronIsolation_rel_FastJet(index, true) < 0.10) cuts_passed |= (1ll<<ELEISO_FASTJET_REL010); // ADDED
    if (electronIsolation_rel_FastJet(index, true) < 0.15) cuts_passed |= (1ll<<ELEISO_FASTJET_REL015); // ADDED

    //relative isolation truncated
    if (electronIsolation_rel(index, true) < 0.10) cuts_passed |= (1ll<<ELEISO_REL010);
    if (electronIsolation_rel(index, true) < 0.15) cuts_passed |= (1ll<<ELEISO_REL015);
    if (electronIsolation_rel(index, true) < 0.40) cuts_passed |= (1ll<<ELEISO_REL040);
    if (electronIsolation_rel(index, true) < 1.00) cuts_passed |= (1ll<<ELEISO_REL100);
    if (electronIsolation_rel_ww(index, true) < 0.10) cuts_passed |= (1ll<<ELEISO_REL010_WW);
    if (electronIsolation_rel_ww(index, true) < 0.40) cuts_passed |= (1ll<<ELEISO_REL040_WW);
    if (electronIsolation_rel_ww(index, true) < 1.00) cuts_passed |= (1ll<<ELEISO_REL100_WW);

    //sliding isolation
//     if( cms2.els_p4()[index].pt() < 15.0 && electronIsolation_rel_v1(index, true) < 0.05) cuts_passed |= (1ll<<ELEISO_SMURFV1);
//     else if( cms2.els_p4()[index].pt() < 20.0 && electronIsolation_rel_v1(index, true) < 0.07) cuts_passed |= (1ll<<ELEISO_SMURFV1);
//     else if( electronIsolation_rel_v1(index, true) < 0.10) cuts_passed |= (1ll<<ELEISO_SMURFV1);

    //pf iso
    float pfiso = electronIsoValuePF(index,0);
    if (fabs(cms2.els_p4()[index].eta()) < 1.479){
      if (pfiso<0.15) cuts_passed |= (1ll<<ELEISO_SMURFV4);
      if (pfiso<0.13) cuts_passed |= (1ll<<ELEISO_SMURFV5);
    } else if (pfiso<0.09) {
      cuts_passed |= (1ll<<ELEISO_SMURFV4);
      cuts_passed |= (1ll<<ELEISO_SMURFV5);
    }
    // ip

    if (fabs(cms2.els_d0corr()[index]) < 0.02) cuts_passed |= (1ll<<ELEIP_200);
    if (fabs(cms2.els_d0corr()[index]) < 0.04) cuts_passed |= (1ll<<ELEIP_400);
    if (fabs(electron_d0PV(index)) < 0.02) cuts_passed |= (1ll<<ELEIP_PV_200);
    if (fabs(electron_d0PV_wwV1(index)) < 0.02 && fabs(electron_dzPV_wwV1(index)) < 1.0 ) cuts_passed |= (1ll<<ELEIP_PV_wwV1);
    if (fabs(electron_d0PV_smurfV3(index)) < 0.02 && fabs(electron_dzPV_smurfV3(index)) < 0.2 ) cuts_passed |= (1ll<<ELEIP_PV_SMURFV3);
    if (fabs(electron_dzPV_smurfV3(index)) < 0.1 ) cuts_passed |= (1ll<<ELEIP_PV_DZ_1MM);
    if (fabs(electron_d0PV_smurfV3(index)) < 0.04 && fabs(electron_dzPV_smurfV3(index)) < 1.0 ) cuts_passed |= (1ll<<ELEIP_PV_OSV2);
    if (fabs(electron_d0PV_smurfV3(index)) < 0.20 && fabs(electron_dzPV_smurfV3(index)) < 1.0 ) cuts_passed |= (1ll<<ELEIP_PV_OSV2_FO);

    if (vertex_index < 0) {
        int vtxidx = firstGoodDAvertex();
        if (vtxidx >= 0) {
            if (useGsfTrack) {
                if (fabs(gsftrks_d0_pv(cms2.els_gsftrkidx()[index], vtxidx, true).first) < 0.02)
                    cuts_passed |= (1ll<<ELEIP_SS200);
            }
            else if (cms2.els_trkidx()[index] >= 0) {            
                if (fabs(trks_d0_pv(cms2.els_trkidx()[index], vtxidx, true).first) < 0.02)
                    cuts_passed |= (1ll<<ELEIP_SS200);  
            }
        }
        else if (fabs(cms2.els_d0corr()[index]) < 0.02)
            cuts_passed |= (1ll<<ELEIP_SS200);  
    }
    else {
        if (useGsfTrack) {
            if (fabs(gsftrks_d0_pv(cms2.els_gsftrkidx()[index], vertex_index, true).first) < 0.02)
                cuts_passed |= (1ll<<ELEIP_SS200);
        }
        else if (cms2.els_trkidx()[index] < 0) {
            if (fabs(cms2.els_d0corr()[index]) < 0.02)
                cuts_passed |= (1ll<<ELEIP_SS200);
        }
        else {
            if (fabs(trks_d0_pv(cms2.els_trkidx()[index], vertex_index, true).first) < 0.02)
                cuts_passed |= (1ll<<ELEIP_SS200);
        }
    }

    // id

    // SMURF ID
    if (electronId_smurf_v1(index)) cuts_passed |= (1ll<<ELEID_SMURFV1_EXTRA);
    if (electronId_smurf_v2(index)) cuts_passed |= (1ll<<ELEID_SMURFV2_EXTRA);
    if (electronId_smurf_v3(index)) cuts_passed |= (1ll<<ELEID_SMURFV3_EXTRA);
    if (electronId_smurf_v1ss(index)) cuts_passed |= (1ll<<ELEID_SMURFV1SS_EXTRA);
    if (electronId_smurf_v2ss(index)) cuts_passed |= (1ll<<ELEID_SMURFV2SS_EXTRA);

    //
    // VBTF ID
    //
    electronIdComponent_t answer_vbtf = 0;
    // VBTF95 (optimised in 35X)
    answer_vbtf = electronId_VBTF(index, VBTF_35X_95, applyAlignmentCorrection, removedEtaCutInEndcap);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_35X_95);
    // VBTF90 (optimised in 35X)
    answer_vbtf = electronId_VBTF(index, VBTF_35X_90, applyAlignmentCorrection, removedEtaCutInEndcap);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_35X_90);
    // VBTF80 (optimised in 35X)
    answer_vbtf = electronId_VBTF(index, VBTF_35X_80, applyAlignmentCorrection, removedEtaCutInEndcap);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_35X_80);
    // VBTF70 (optimised in 35X)
    answer_vbtf = electronId_VBTF(index, VBTF_35Xr2_70, applyAlignmentCorrection, removedEtaCutInEndcap);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_35X_70);
    // VBTF85 no H/E in endcap
    answer_vbtf = electronId_VBTF(index, VBTF_85_NOHOEEND, applyAlignmentCorrection, removedEtaCutInEndcap);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_85_NOHOEEND);
    // VBTF85
    answer_vbtf = electronId_VBTF(index, VBTF_85 , applyAlignmentCorrection, removedEtaCutInEndcap);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_85);
    // VBTF80 no H/E in endcap
    answer_vbtf = electronId_VBTF(index, VBTF_80_NOHOEEND, applyAlignmentCorrection, removedEtaCutInEndcap);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_80_NOHOEEND);
    // VBTF70 no H/E in endcap
    answer_vbtf = electronId_VBTF(index, VBTF_70_NOHOEEND, applyAlignmentCorrection, removedEtaCutInEndcap);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_70_NOHOEEND);
    // VBTF90 with H/E and dPhiIn tuned to match HLT
    answer_vbtf = electronId_VBTF(index, VBTF_90_HLT, applyAlignmentCorrection, removedEtaCutInEndcap);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_90_HLT);
    // VBTF90 with H/E and dPhiIn tuned to match HLT (CaloIdT+TrkIdVL)
    answer_vbtf = electronId_VBTF(index, VBTF_90_HLT_CALOIDT_TRKIDVL, applyAlignmentCorrection, removedEtaCutInEndcap);
    if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_VBTF_90_HLT_CALOIDT_TRKIDVL);
    //
    // CIC ID  
    //
    // MEDIUM (V03 optimisation)
    electronIdComponent_t answer_cic = electronId_CIC(index, 3, CIC_MEDIUM, applyAlignmentCorrection, removedEtaCutInEndcap);
    if ((answer_cic & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) cuts_passed |= (1ll<<ELEID_CIC_V03_MEDIUM);

    //
    // conversion rejection cuts
    //

    if (!isFromConversionPartnerTrack(index)) cuts_passed |= (1ll<<ELENOTCONV_DISTDCOT002);
    if (!isFromConversionHitPattern(index)) cuts_passed |= (1ll<<ELENOTCONV_HITPATTERN);
    if (cms2.els_exp_innerlayers().at(index) == 0) cuts_passed |= (1ll<<ELENOTCONV_HITPATTERN_0MHITS);
    if(!isFromConversionMIT(index)) cuts_passed |= (1ll<<ELENOTCONV_MIT);

    //
    // fiduciality/other cuts
    //

    if ((cms2.els_type()[index] & (1ll<<ISECALDRIVEN))) cuts_passed |= (1ll<<ELESEED_ECAL);
    if (fabs(cms2.els_p4()[index].eta()) < 2.5) cuts_passed |= (1ll<<ELEETA_250);
    if (fabs(cms2.els_p4()[index].eta()) < 2.4) cuts_passed |= (1ll<<ELEETA_240);
    if (electronId_noMuon(index)) cuts_passed |= (1ll<<ELENOMUON_010);
    if (electronId_noMuon_SS(index)) cuts_passed |= (1ll<<ELENOMUON_010_SS);

    //
    // Pt
    //
    if( cms2.els_p4()[index].pt() > 10.0 ) cuts_passed |= (1ll<<ELEPT_010);
    if( cms2.els_p4()[index].pt() > 15.0 ) cuts_passed |= (1ll<<ELEPT_015);
    if( cms2.els_p4()[index].pt() > 20.0 ) cuts_passed |= (1ll<<ELEPT_020);

    //
    // Super Cluster Et
    //
    if( ( cms2.els_eSC()[index] / cosh(cms2.els_etaSC()[index]) ) > 10.0 ) cuts_passed |= (1ll<<ELESCET_010);
    if( ( cms2.els_eSC()[index] / cosh(cms2.els_etaSC()[index]) ) > 15.0 ) cuts_passed |= (1ll<<ELESCET_015);

    //
    // Veto electron in transition region
    //
    if( fabs(cms2.els_etaSC()[index]) < 1.4442 || fabs(cms2.els_etaSC()[index]) > 1.566 )  cuts_passed |= (1ll<<ELE_NOT_TRANSITION);

    //
    // chargeflip
    //

    if (!isChargeFlip(index)) cuts_passed |= (1ll<<ELECHARGE_NOTFLIP); 
    if (!isChargeFlip3agree(index)) cuts_passed |= (1ll<<ELECHARGE_NOTFLIP3AGREE);

    //
    // spike rejection
    //
    if (!isSpikeElectron(index)) cuts_passed |= (1ll<<ELENOSPIKE_SWISS005);

    //
    // return which selections passed
    //

    return cuts_passed;

}

//
// if fbrem is low then cut on e/p_in
//
bool electronId_extra(const unsigned int index)
{
    if (cms2.els_fbrem()[index] < 0.2) {
        if (cms2.els_eOverPIn()[index] < 0.7 || cms2.els_eOverPIn()[index] > 1.5) return false;
    }

    return true;
}

//
// remove if close to a muon
//
bool electronId_noMuon(const unsigned int index)
{
    if ( cms2.els_closestMuon().at(index) != -1) return false;
    return true;
}

//
// remove if close to a good muon
//
bool electronId_noMuon_SS(const unsigned int index)
{
    unsigned int idx = cms2.els_closestMuon().at(index);

    if (idx < 0) return true;
    if (muonId(idx, NominalSSv4)) return false;
    return true;
}

//
// candidate electron id function
//
bool electronId_cand(const unsigned int index, const cand_tightness tightness, bool applyAlignementCorrection, bool removedEtaCutInEndcap)
{

    std::vector<double> relisoThresholds;
    std::vector<double> dEtaInThresholds;
    std::vector<double> dPhiInThresholds;
    std::vector<double> hoeThresholds;
    std::vector<double> latThresholds;

    eidGetCand(tightness, dEtaInThresholds, dPhiInThresholds, hoeThresholds, latThresholds);

    //
    // get corrected dEtaIn and dPhiIn
    //

    float dEtaIn = cms2.els_dEtaIn()[index];
    float dPhiIn = cms2.els_dPhiIn()[index];
    if (applyAlignementCorrection) electronCorrection_pos(index, dEtaIn, dPhiIn);

    //
    // apply cuts
    //
    if (fabs(cms2.els_etaSC()[index]) < 1.479) {
        if (fabs(dEtaIn) > dEtaInThresholds[0]) 	return false;
        if (fabs(dPhiIn) > dPhiInThresholds[0]) 	return false;
        if (cms2.els_hOverE()[index] > hoeThresholds[0]) 		    return false;
        if ((cms2.els_e2x5Max()[index]/cms2.els_e5x5()[index]) < latThresholds[0]) return false;
    }
    if (fabs(cms2.els_etaSC()[index]) > 1.479) {
        if (!removedEtaCutInEndcap && (fabs(dEtaIn) > dEtaInThresholds[1])) 	return false;
        if (fabs(dPhiIn) > dPhiInThresholds[1]) 	return false;
        if (cms2.els_hOverE()[index] > hoeThresholds[1]) 		    return false;
        if (cms2.els_sigmaIEtaIEta()[index] > latThresholds[1])  return false;	
    }

    return true;

}

bool electronId_smurf_v1(const unsigned int index)
{

  if (cms2.els_p4()[index].pt() > 20.0) return true;

  if (cms2.els_fbrem()[index] > 0.15) return true;
  
  if (fabs(cms2.els_etaSC()[index]) < 1.) {
    if (cms2.els_eOverPIn()[index] > 0.95 && fabs(cms2.els_dPhiIn()[index]*cms2.els_charge()[index]) < 0.006) return true; 
  }

  return false;
}

bool electronId_smurf_v2(const unsigned int index)	 
 {	 
 	 
   if (cms2.els_p4()[index].pt() > 20.0) return true;	 
 	 
   if (cms2.els_fbrem()[index] > 0.15) return true;	 
 	 
   if (fabs(cms2.els_etaSC()[index]) < 1.) {	 
     if (cms2.els_eOverPIn()[index] > 0.95) return true;	 
   }	 
 	 
   return false;	 
 }	 
 	 
 bool electronId_smurf_v3(const unsigned int index)	 
 {	 
   if (fabs(cms2.els_etaSC()[index]) > 1.479 && cms2.els_hOverE()[index] > 0.1) return false;
 	 
   if (cms2.els_p4()[index].pt() > 20.0) return true;	 
 	 
   electronIdComponent_t answer_vbtf = 0;	 
   answer_vbtf = electronId_VBTF(index, VBTF_70_NOHOEEND, false, false);	 
   if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) {	 
 	 
     if (cms2.els_fbrem()[index] > 0.15) return true;	 
 	 
     if (fabs(cms2.els_etaSC()[index]) < 1.) {	 
       if (cms2.els_eOverPIn()[index] > 0.95) return true;	 
     }	 
 	 
   }	 
 	 
   return false;	 
 }

 bool electronId_smurf_v1ss(const unsigned int index)
 {

   if (cms2.els_p4()[index].pt() > 20.0) return true;

   electronIdComponent_t answer_vbtf = 0;
   answer_vbtf = electronId_VBTF(index, VBTF_80_NOHOEEND, false, false);
   if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) {

     if (cms2.els_fbrem()[index] > 0.15) return true;

     if (fabs(cms2.els_etaSC()[index]) < 1.) {
       if (cms2.els_eOverPIn()[index] > 0.95) return true;
     }

   }

   return false;
 }

 bool electronId_smurf_v2ss(const unsigned int index)
 {

   if (cms2.els_p4()[index].pt() > 20.0) return true;

   electronIdComponent_t answer_vbtf = 0;
   answer_vbtf = electronId_VBTF(index, VBTF_80_NOHOEEND, false, false);
   if ((answer_vbtf & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID)) {

       if (fabs(cms2.els_etaSC()[index]) > 1.479)
           if (cms2.els_hOverE()[index] < 0.15) return true;

     if (cms2.els_fbrem()[index] > 0.15) return true;

     if (fabs(cms2.els_etaSC()[index]) < 1.) {
       if (cms2.els_eOverPIn()[index] > 0.95) return true;
     }

   }

   return false;
 }



//
// class based id that is new/experimental
//

electronIdComponent_t electronId_CIC(const unsigned int index, const unsigned int version, const cic_tightness tightness, bool applyAlignementCorrection, bool removedEtaCutInEndcap)
{

    // check that a valid version number was supplied
    if (version != 2 && version != 3 && version != 4 && version != 6) {
        std::cout << "[electronId_CIC] Error! Version must be 2, 3, 4 or 6 - fail" << std::endl;
        return 0;
    }

    //
    // get the variables that are going to be cut on
    //

    double scTheta = (2*atan(exp(-1*cms2.els_etaSC()[index])));
    double scEt = cms2.els_eSC()[index]*sin(scTheta); 
    double scEta = cms2.els_etaSC()[index];
    double eSeedOverPin = cms2.els_eSeedOverPIn()[index];
    double fBrem = cms2.els_fbrem()[index];
    double hOverE = cms2.els_hOverE()[index];
    double sigmaee = cms2.els_sigmaIEtaIEta()[index];
    int mishits = cms2.els_exp_innerlayers()[index];
    double dist = (cms2.els_conv_dist()[index] == -9999.? 9999:cms2.els_conv_dist()[index]);
    double dcot = (cms2.els_conv_dcot()[index] == -9999.? 9999:cms2.els_conv_dcot()[index]);
    float dcotdistcomb = ((0.04 - std::max(dist, dcot)) > 0?(0.04 - std::max(dist, dcot)):0);
    double tkIso = cms2.els_tkIso()[index];
    double ecalIso = cms2.els_ecalIso04()[index];
    double hcalIso = cms2.els_hcalIso04()[index];

    double ip = fabs(cms2.els_d0()[index]);
    if (cms2.vtxs_sumpt().size() > 0) {
        const float vx = cms2.vtxs_position()[0].x();
        const float vy = cms2.vtxs_position()[0].y();
        const float px = cms2.els_trk_p4()[index].px();
        const float py = cms2.els_trk_p4()[index].py();
        ip = fabs(-1*cms2.els_d0()[index] + (vx*py - vy*px) / cms2.els_trk_p4()[index].pt());
    }

    //
    // get corrected dEtaIn and dPhiIn
    //

    float deltaEtaIn = cms2.els_dEtaIn()[index];
    float deltaPhiIn = cms2.els_dPhiIn()[index];
    if (applyAlignementCorrection) electronCorrection_pos(index, deltaEtaIn, deltaPhiIn);

    // find the catagory for this electron
    unsigned int cat = eidClassify(version, index);

    // determine if in EB or EE
    int eb;
    if (cms2.els_fiduciality()[index] & (1<<ISEB)) eb = 0;
    else eb = 1;

    //
    // Version V02
    //

    if (version == 2) {

        //
        // set the parameters for the chosen tightness
        //

        std::vector<double> cutdeta;
        std::vector<double> cutdphi;
        std::vector<double> cuteopin;
        std::vector<double> cutet;
        std::vector<double> cuthoe;
        std::vector<double> cutip;
        std::vector<double> cutisoecal;
        std::vector<double> cutisohcal;
        std::vector<double> cutisotk;
        std::vector<double> cutmishits;
        std::vector<double> cutsee;

        eidGetCIC_V02(tightness, cutdeta, cutdphi, cuteopin, cutet, cuthoe, 
                cutip, cutisoecal, cutisohcal, cutisotk, cutmishits, cutsee);


        unsigned int result = 0;
        int bin = 0;
        if (scEt < 20.)
            bin = 2;
        else if (scEt > 30.)
            bin = 0;
        else
            bin = 1;

        if (fBrem > 0)
            eSeedOverPin = eSeedOverPin + fBrem;

        if (bin != 2) {     
            tkIso = tkIso*pow(40./scEt, 2); 
            ecalIso = ecalIso*pow(40./scEt, 2); 
            hcalIso = hcalIso*pow(40./scEt, 2); 
        }

        if ((tkIso < cutisotk[cat+3*eb+bin*6]) &&
                (ecalIso < cutisoecal[cat+3*eb+bin*6]) &&
                (hcalIso < cutisohcal[cat+3*eb+bin*6]))
            result |= (1<<ELEID_ISO);

        if (fBrem < -2)
            return result;

        bool passdEtaCut = fabs(deltaEtaIn) < cutdeta[cat+3*eb+bin*6];


        if(removedEtaCutInEndcap && eb == 1) passdEtaCut = true;
        if (hOverE < cuthoe[cat+3*eb+bin*6] &&
                sigmaee < cutsee[cat+3*eb+bin*6] &&
                fabs(deltaPhiIn) < cutdphi[cat+3*eb+bin*6] &&
                passdEtaCut &&
                eSeedOverPin > cuteopin[cat+3*eb+bin*6])
            result |= (1<<ELEID_ID);

        if (ip < cutip[cat+3*eb+bin*6])
            result |= (1<<ELEID_IP);

        if (mishits < cutmishits[cat+3*eb+bin*6])
            result |= (1<<ELEID_CONV);

        return result;
    }

    //
    // version V03, V04 or V05
    //

    if (version == 3 || version == 4 || version == 5) {

        //
        // set the parameters for the chosen tightness
        //
        std::vector<double> cutdcotdist;
        std::vector<double> cutdetain;
        std::vector<double> cutdphiin;
        std::vector<double> cuteseedopcor;
        std::vector<double> cutet;
        std::vector<double> cutfmishits;
        std::vector<double> cuthoe;
        std::vector<double> cutip_gsf;
        std::vector<double> cutiso_sum;
        std::vector<double> cutiso_sumoet;
        std::vector<double> cutsee;

        // V03 uses Et binning
        bool wantBinning = true;

        if (version == 3) {
            eidGetCIC_V03(tightness, cutdcotdist, cutdetain, cutdphiin, cuteseedopcor, cutet,
                    cutfmishits, cuthoe, cutip_gsf, cutiso_sum, cutiso_sumoet, cutsee);
        }

        if (version == 4) {
            eidGetCIC_V04(tightness, cutdcotdist, cutdetain, cutdphiin, cuteseedopcor, cutet,
                    cutfmishits, cuthoe, cutip_gsf, cutiso_sum, cutiso_sumoet, cutsee);
            // V04 does not use Et binning
            wantBinning = false;
        }

        // this is certainly true for V03
        // but not sure the meaning of V04,05
        unsigned int result = 0;
        int bin = 0;
        if (wantBinning) {
            if (scEt < 20.)
                bin = 2;
            else if (scEt > 30.)
                bin = 0;
            else
                bin = 1;
        }

        if (fBrem > 0)
            eSeedOverPin = eSeedOverPin + fBrem;

        float iso_sum = tkIso + ecalIso + hcalIso;
        float iso_sum_corrected = iso_sum*pow(40./scEt, 2);

        if ((iso_sum < cutiso_sum[cat+bin*9]) &&
                (iso_sum_corrected < cutiso_sumoet[cat+bin*9]))
            result |= (1<<ELEID_ISO);

        if (fBrem > -2) {
            if ((hOverE < cuthoe[cat+bin*9]) and
                    (sigmaee < cutsee[cat+bin*9]) and
                    (fabs(deltaPhiIn) < cutdphiin[cat+bin*9]) and
                    (fabs(deltaEtaIn) < cutdetain[cat+bin*9]) and
                    (eSeedOverPin > cuteseedopcor[cat+bin*9]) and
                    (scEt > cutet[cat+bin*9]))
                result |= (1<<ELEID_ID);
        }

        if (ip < cutip_gsf[cat+bin*9])
            result |= (1<<ELEID_IP);

        if ((mishits < cutfmishits[cat+bin*9]) and 
                (dcotdistcomb < cutdcotdist[cat+bin*9]))
            result |= (1<<ELEID_CONV);

        return result;
    }

    if (version == 6) {

        std::vector<double> cutIsoSum;
        std::vector<double> cutIsoSumCorr;
        std::vector<double> cuthoe;
        std::vector<double> cutsee;
        std::vector<double> cutdphi;
        std::vector<double> cutdeta;
        std::vector<double> cuteopin;
        std::vector<double> cutmishits;
        std::vector<double> cutdcotdist;
        std::vector<double> cutip;
        std::vector<double> cutIsoSumCorrl;
        std::vector<double> cuthoel;
        std::vector<double> cutseel;
        std::vector<double> cutdphil;
        std::vector<double> cutdetal;
        std::vector<double> cutipl;

        eidGetCIC_V06(tightness, cutIsoSum, cutIsoSumCorr, cuthoe, cutsee, cutdphi, 
                cutdeta, cuteopin, cutmishits, cutdcotdist, cutip, cutIsoSumCorrl, 
                cuthoel, cutseel, cutdphil, cutdetal, cutipl);

        int result = 0;
        const int ncuts = 10;
        std::vector<bool> cut_results(ncuts, false);

        float iso_sum = tkIso + ecalIso + hcalIso;
        if(fabs(scEta)>1.5) 
            iso_sum += (fabs(scEta)-1.5)*1.09;

        float iso_sumoet = iso_sum*(40./scEt);
        float eseedopincor = eSeedOverPin + fBrem;
        if(fBrem < 0) eseedopincor = eSeedOverPin;

        for (int cut=0; cut<ncuts; cut++) {

            switch (cut) {
                case 0:
                    cut_results[cut] = eidComputeCut(fabs(deltaEtaIn), scEt, cutdetal[cat], cutdeta[cat]);
                    break;
                case 1:
                    cut_results[cut] = eidComputeCut(fabs(deltaPhiIn), scEt, cutdphil[cat], cutdphi[cat]);
                    break;
                case 2:
                    cut_results[cut] = (eseedopincor > cuteopin[cat]);
                    break;
                case 3:
                    cut_results[cut] = eidComputeCut(hOverE, scEt, cuthoel[cat], cuthoe[cat]);
                    break;
                case 4:
                    cut_results[cut] = eidComputeCut(sigmaee, scEt, cutseel[cat], cutsee[cat]);
                    break;
                case 5:
                    cut_results[cut] = eidComputeCut(iso_sumoet, scEt, cutIsoSumCorrl[cat], cutIsoSumCorr[cat]);
                    break;
                case 6:
                    cut_results[cut] = (iso_sum < cutIsoSum[cat]);
                    break;
                case 7:
                    cut_results[cut] = eidComputeCut(fabs(ip), scEt, cutipl[cat], cutip[cat]);
                    break;
                case 8:
                    cut_results[cut] = (mishits < cutmishits[cat]);
                    break;
                case 9:
                    cut_results[cut] = (dcotdistcomb < cutdcotdist[cat]);
                    break;
            }

        }

        // ID part
        if (cut_results[0] & cut_results[1] & cut_results[2] & cut_results[3] & cut_results[4]) result |= (1<<ELEID_ID);   

        // ISO part
        if (cut_results[5] & cut_results[6]) result |= (1<<ELEID_ISO);

        // IP part
        if (cut_results[7]) result |= (1<<ELEID_IP);

        // Conversion part
        if (cut_results[8] and cut_results[9]) result |= (1<<ELEID_CONV);

        return result;

    }

    std::cout << "[electronId_CIC] Error! got to the end and didn't apply any cuts - fail" << std::endl;
    return 0;

}

unsigned int eidClassify(const unsigned int version, const unsigned int index) {

    // variables used for classification
    double eta = fabs(cms2.els_etaSC()[index]);
    double eOverP = cms2.els_eOverPIn()[index];
    double fBrem = cms2.els_fbrem()[index];
    bool isEB = (cms2.els_fiduciality()[index] & (1<<ISEB));
    bool isEE = (cms2.els_fiduciality()[index] & (1<<ISEE));
    bool ecalDriven = (cms2.els_type()[index] & (1<<ISECALDRIVEN));
    bool trackerDriven = (cms2.els_type()[index] & (1<<ISTRACKERDRIVEN));

    int cat = -1;

    //
    // version V00 or V01
    //
    if (version == 0 || version == 1) {
        if((isEB && fBrem<0.06) || (isEE && fBrem<0.1)) 
            cat=1;
        else if (eOverP < 1.2 && eOverP > 0.8) 
            cat=0;
        else 
            cat=2;
        return cat;
    }

    //
    // version V02
    //
    if (version == 2) {
        if (isEB) {       // BARREL
            if(fBrem < 0.12)
                cat=1;
            else if (eOverP < 1.2 && eOverP > 0.9) 
                cat=0;
            else 
                cat=2;
        } else {                     // ENDCAP
            if(fBrem < 0.2)
                cat=1;
            else if (eOverP < 1.22 && eOverP > 0.82) 
                cat=0;
            else 
                cat=2;
        }
        return cat;
    }

    //
    // version V03, V04 or V05
    // took this from V00-03-07-02 of 
    // ElectronIdentification/src/CutBasedElectronID.cc
    // NOTE - DLE - now adding version 6
    // this means newCategories is switched on
    //


    if (version == 3 || version == 4 || version == 5 || version == 6) {

        // this is certainly true for V03 and V04
        // not sure about V05
        bool newCategories = false;
        if (version == 6) newCategories = true;

        if (isEB) {
            if ((fBrem >= 0.12) && (eOverP > 0.9) && (eOverP < 1.2))
                cat = 0;
            else if (((eta >  .445   && eta <  .45  ) ||
                        (eta >  .79    && eta <  .81  ) ||
                        (eta > 1.137   && eta < 1.157 ) ||
                        (eta > 1.47285 && eta < 1.4744)) && newCategories)
                cat = 6;
            else if (trackerDriven && !ecalDriven && newCategories)
                cat = 8;
            else if (fBrem < 0.12)
                cat = 1;
            else
                cat = 2;
        } else {
            if ((fBrem >= 0.2) && (eOverP > 0.82) && (eOverP < 1.22))
                cat = 3;
            else if (eta > 1.5 && eta <  1.58 && newCategories)
                cat = 7;
            else if (trackerDriven && !ecalDriven && newCategories)
                cat = 8;
            else if (fBrem < 0.2)
                cat = 4;
            else
                cat = 5;
        }

        return cat;
    }


    return cat;

}

bool eidComputeCut(double x, double et, double cut_min, double cut_max, bool gtn)
{

  float et_min = 10;
  float et_max = 40;
  bool accept = false;
  float cut = cut_max; //  the cut at et=40 GeV

  if(et < et_max) {
    cut = cut_min + (1/et_min - 1/et)*(cut_max - cut_min)/(1/et_min - 1/et_max);
  } 
  
  if(et < et_min) {
    cut = cut_min;
  } 

  if(gtn) {   // useful for e/p cut which is gt
    accept = (x >= cut);
  } 
  else {
    accept = (x <= cut);
  }

  //std::cout << x << " " << cut_min << " " << cut << " " << cut_max << " " << et << " " << accept << std::endl;
  return accept;
}

//
// VBTF stuff
//

electronIdComponent_t electronId_VBTF(const unsigned int index, const vbtf_tightness tightness, bool applyAlignementCorrection, bool removedEtaCutInEndcap)
{

    unsigned int answer = 0;

    std::vector<double> relisoThresholds;
    std::vector<double> dEtaInThresholds;
    std::vector<double> dPhiInThresholds;
    std::vector<double> hoeThresholds;
    std::vector<double> sigmaIEtaIEtaThresholds;

    eidGetVBTF(tightness, dEtaInThresholds, dPhiInThresholds, hoeThresholds, 
            sigmaIEtaIEtaThresholds, relisoThresholds);

    //
    // get corrected dEtaIn and dPhiIn
    //

    float dEtaIn = cms2.els_dEtaIn()[index];
    float dPhiIn = cms2.els_dPhiIn()[index];
    if (applyAlignementCorrection) electronCorrection_pos(index, dEtaIn, dPhiIn);

    // barrel
    if (fabs(cms2.els_etaSC()[index]) < 1.479) {

        if (electronIsolation_rel(index, true) < relisoThresholds[0])
            answer |= (1<<ELEID_ISO);

        if (fabs(dEtaIn) < dEtaInThresholds[0] &&
                fabs(dPhiIn) < dPhiInThresholds[0] &&
                cms2.els_hOverE()[index] < hoeThresholds[0] &&
                cms2.els_sigmaIEtaIEta()[index] < sigmaIEtaIEtaThresholds[0])
            answer |= (1<<ELEID_ID);
    }

    // endcap
    if (fabs(cms2.els_etaSC()[index]) > 1.479) {
        if (electronIsolation_rel(index, true) < relisoThresholds[1])
            answer |= (1<<ELEID_ISO);
        bool passdEtaCut = fabs(dEtaIn) < dEtaInThresholds[1];
        if(removedEtaCutInEndcap) passdEtaCut = true;
        if ( passdEtaCut &&
                fabs(dPhiIn) < dPhiInThresholds[1] &&
                cms2.els_hOverE()[index] < hoeThresholds[1] &&
                cms2.els_sigmaIEtaIEta()[index] < sigmaIEtaIEtaThresholds[1])
            answer |= (1<<ELEID_ID);
    }

    return answer;

}

//
// electron isolation definitions 
//

//relative truncated
float electronIsolation_rel(const unsigned int index, bool use_calo_iso)
{
    float sum = cms2.els_tkIso().at(index);
    if (use_calo_iso) {
        if (fabs(cms2.els_etaSC().at(index)) > 1.479) sum += cms2.els_ecalIso().at(index);
        if (fabs(cms2.els_etaSC().at(index)) <= 1.479) sum += max(0., (cms2.els_ecalIso().at(index) -1.));
        sum += cms2.els_hcalIso().at(index);
    }
    double pt = cms2.els_p4().at(index).pt();
    return sum/max(pt, 20.);
}

//relative truncated, fast-jet corrected
float electronIsolation_rel_FastJet(const unsigned int index, bool use_calo_iso)
{
    float sum = cms2.els_tkIso().at(index);

    float offset = cms2.evt_rho() * TMath::Pi() * pow( 0.3 , 2 );
    
    if (use_calo_iso) {
        float caloiso = 0.;

        if (fabs(cms2.els_etaSC().at(index)) >  1.479) caloiso += cms2.els_ecalIso().at(index);
        if (fabs(cms2.els_etaSC().at(index)) <= 1.479) caloiso += max(0., (cms2.els_ecalIso().at(index) -1.));
        caloiso += cms2.els_hcalIso().at(index);
        
        caloiso -= offset;
        if( caloiso > 0 ) sum += caloiso;
    }
    double pt = cms2.els_p4().at(index).pt();
    return sum/max(pt, 20.);
}


/*
//relative non-truncated
float electronIsolation_rel_v1Original(const unsigned int index, bool use_calo_iso)
{
    float sum = cms2.els_tkIso().at(index);
    if (use_calo_iso) {
        if (fabs(cms2.els_etaSC().at(index)) > 1.479) sum += cms2.els_ecalIso().at(index);
        if (fabs(cms2.els_etaSC().at(index)) <= 1.479) sum += max(0., (cms2.els_ecalIso().at(index) -1.));
        sum += cms2.els_hcalIso().at(index);
    }
    double pt = cms2.els_p4().at(index).pt();
    return sum/pt;
}
*/

// Relative Isolation, Non-Truncated
float electronIsolation_rel_v1(const unsigned int index, bool use_calo_iso){
    float pt               = cms2.els_p4().at(index).pt();          // Electron Pt
    float TRCK_sum_over_pt = cms2.els_tkIso().at(index) / pt;       // Tracker Relative Isolation, Non-Truncated
    float ECAL_sum_over_pt = electronIsolation_ECAL_rel_v1(index);  // ECAL    Relative Isolation, Non-Truncated
    float HCAL_sum_over_pt = electronIsolation_HCAL_rel_v1(index);  // HCAL    Relative Isolation, Non-Truncated

    float sum_over_pt      = TRCK_sum_over_pt;                      // Combined Subdetector Relative Isolation, Non-Truncated
    if(use_calo_iso){
      sum_over_pt += ECAL_sum_over_pt;
      sum_over_pt += HCAL_sum_over_pt;
    }
    return sum_over_pt;
}

// Relative Isolation, Non-Truncated, FastJet-corrected
float electronIsolation_rel_v1_FastJet(const unsigned int index, bool use_calo_iso){
    float pt               = cms2.els_p4().at(index).pt();          // Electron Pt
    float TRCK_sum_over_pt = cms2.els_tkIso().at(index) / pt;       // Tracker Relative Isolation, Non-Truncated
    float ECAL_sum_over_pt = electronIsolation_ECAL_rel_v1(index);  // ECAL    Relative Isolation, Non-Truncated
    float HCAL_sum_over_pt = electronIsolation_HCAL_rel_v1(index);  // HCAL    Relative Isolation, Non-Truncated
    float offset           = el_fastjet_rel_offset(index);          // fastjet offset = pi X dR^2 X rho
    float sum_over_pt      = TRCK_sum_over_pt;                      // Combined Subdetector Relative Isolation, Non-Truncated
    if(use_calo_iso){
      float calo_iso = ECAL_sum_over_pt + HCAL_sum_over_pt - offset;
      if( calo_iso > 0 ) sum_over_pt += calo_iso;
    }
    return sum_over_pt;
}

// ECAL Relative Isolation, Non-Truncated
float electronIsolation_ECAL_rel_v1(const unsigned int index, bool useEBps){
  float pt               = cms2.els_p4().at(index).pt();                                                                  // Electron Pt
  float ecal_sum_over_pt = 0.0;                                                                                           // ECAL Relative Isolation, NT
  if( fabs(cms2.els_etaSC().at(index)) > 1.479  ) ecal_sum_over_pt += cms2.els_ecalIso().at(index);                       // EE: Ecal Endcap  
  if( fabs(cms2.els_etaSC().at(index)) <= 1.479 ) {
      if (useEBps)
          ecal_sum_over_pt += max( 0.0, ( cms2.els_ecalIso().at(index) - 1.0 ) ); // EB: Ecal Barrel
      else
          ecal_sum_over_pt += cms2.els_ecalIso().at(index); // EB: Ecal Barrel
  }
  ecal_sum_over_pt /= pt;
  return ecal_sum_over_pt;
}

// HCAL Relative Isolation, Non-Truncated
float electronIsolation_HCAL_rel_v1(const unsigned int index){
  float pt               = cms2.els_p4().at(index).pt();      // Electron Pt
  float hcal_sum_over_pt = cms2.els_hcalIso().at(index) / pt; // HCAL Relative Isolation, NT
  return hcal_sum_over_pt;
}

float el_fastjet_rel_offset(const unsigned int index){
  double pt     = cms2.els_p4().at(index).pt();
  double offset = TMath::Pi() * pow( 0.3 , 2 ) * cms2.evt_rho();
  return offset / pt;
}

// ECAL Relative Isolation, Truncated
float electronIsolation_ECAL_rel(const unsigned int index){
  float pt               = cms2.els_p4().at(index).pt();                                                                  // Electron Pt
  float ecal_sum_over_pt = 0.0;                                                                                           // ECAL Relative Isolation
  if( fabs(cms2.els_etaSC().at(index)) > 1.479  ) ecal_sum_over_pt += cms2.els_ecalIso().at(index);                       // EE: Ecal Endcap
  if( fabs(cms2.els_etaSC().at(index)) <= 1.479 ) ecal_sum_over_pt += max( 0.0, ( cms2.els_ecalIso().at(index) - 1.0 ) ); // EB: Ecal Barrel
  ecal_sum_over_pt /= max(pt,(float)20.0);
  return ecal_sum_over_pt;
}

// HCAL Relative Isolation, Truncated
float electronIsolation_HCAL_rel(const unsigned int index){
  float pt               = max( cms2.els_p4().at(index).pt() , (float) 20.0 ) ;      // Electron Pt
  float hcal_sum_over_pt = cms2.els_hcalIso().at(index) / pt;                        // HCAL Relative Isolation, NT
  return hcal_sum_over_pt;
}



//
// electron isolation definitions for WW analysis
//
float electronIsolation_rel_ww(const unsigned int index, bool use_calo_iso)
{
    float sum = cms2.els_tkIso().at(index);
    if(use_calo_iso)
        sum += max(0., (cms2.els_ecalIso().at(index) -1.));
    sum += cms2.els_hcalIso().at(index);
    double pt = cms2.els_p4().at(index).pt();
    return sum/max(pt, 20.);
}


float electronIsoValuePF(const unsigned int iel, unsigned int idavtx, float coner, float minptn, float dzcut, 
			 float footprintdr, float gammastripveto, float elestripveto){

  int elgsftkid = cms2.els_gsftrkidx().at(iel);
  int eltkid = cms2.els_trkidx().at(iel);
  //take dz from gsf, and if it does not exist (should always exist) take it from ctf track
  float eldz = elgsftkid>=0 ? gsftrks_dz_dapv( elgsftkid,idavtx ).first : trks_dz_dapv(eltkid,idavtx).first;
  float eleta = cms2.els_p4().at(iel).eta();

  float pfciso = 0.;
  float pfniso = 0.;
  float pffootprint = 0.;
  float pfjurveto = 0.;
  float pfjurvetoq = 0.;
  for (unsigned int ipf=0; ipf<cms2.pfcands_p4().size(); ++ipf){

    float dR = dRbetweenVectors(cms2.pfcands_p4().at(ipf),cms2.els_p4().at(iel));
    if (dR>coner) continue;

    float pfpt = cms2.pfcands_p4().at(ipf).pt();    
    float pfeta = cms2.pfcands_p4().at(ipf).eta();    
    float deta = fabs(pfeta - eleta);
    int pfid = abs(cms2.pfcands_particleId().at(ipf));
    if (cms2.pfcands_charge().at(ipf)==0) {
      //neutrals
      if (pfpt>minptn) {
	pfniso+=pfpt;
	if (dR<footprintdr && pfid==130) pffootprint+=pfpt;
	if (deta<gammastripveto && pfid==22)  pfjurveto+=pfpt;
      }
    } else {
      //charged  
      //avoid double counting of electron itself
      //if either the gsf or the ctf track are shared with the candidate, skip it
      int pftkid = cms2.pfcands_trkidx().at(ipf);
      if (eltkid>=0 && pftkid>=0 && eltkid==pftkid) continue;
      if (pfid==11 && cms2.pfcands_pfelsidx().at(ipf)>=0 && cms2.pfels_elsidx().at(cms2.pfcands_pfelsidx().at(ipf))>=0) {
	int pfgsfid = cms2.els_gsftrkidx().at(cms2.pfels_elsidx().at(cms2.pfcands_pfelsidx().at(ipf))); 
	if (elgsftkid>=0 && pfgsfid>=0 && elgsftkid==pfgsfid) continue;
      }
      //check electrons with gsf track
      if (pfid==11 && cms2.pfcands_pfelsidx().at(ipf)>=0 && cms2.pfels_elsidx().at(cms2.pfcands_pfelsidx().at(ipf))>=0) {
	int gsfid = cms2.els_gsftrkidx().at(cms2.pfels_elsidx().at(cms2.pfcands_pfelsidx().at(ipf))); 
	if (gsfid>=0) { 
	  if(fabs(gsftrks_dz_dapv( gsfid,idavtx ).first - eldz )<dzcut) {//dz cut
	    pfciso+=pfpt;
	    if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
	  }   
	  continue;//and avoid double counting
	}
      }      
      //then check anything that has a ctf track
      if (pftkid>=0) {//charged (with a ctf track)
	if(fabs( trks_dz_dapv(cms2.pfcands_trkidx().at(ipf),idavtx).first - eldz )<dzcut) {//dz cut
	  pfciso+=pfpt;
	  if (deta<elestripveto && pfid==11) pfjurvetoq+=pfpt;
	}
      }
    } 
  }

  return (pfciso+pfniso-pffootprint-pfjurveto-pfjurvetoq)/cms2.els_p4().at(iel).pt();

}



//
//conversion rejection
//
bool isFromConversionHitPattern(const unsigned int index)
{
    if(cms2.els_exp_innerlayers().at(index) > 1) return true;
    return false;
}

bool isFromConversionPartnerTrack(const unsigned int index) {

    if(fabs(cms2.els_conv_dist().at(index)) < 0.02 &&
            fabs(cms2.els_conv_dcot().at(index)) < 0.02)
        return true;

    return false;

}

int getChargeUsingMajorityLogic(int elIdx, float minFracSharedHits) {

    if(cms2.els_sccharge()[elIdx]*cms2.els_trk_charge()[elIdx] > 0 || (cms2.els_trkidx()[elIdx] < 0 || cms2.els_trkshFrac()[elIdx] < minFracSharedHits))
        return cms2.els_sccharge()[elIdx];
    else 
        return  cms2.trks_charge().at(cms2.els_trkidx().at(elIdx));

}

bool isChargeFlip(int elIndex){
    //true if electron is likely to be a charge flip
    if ((cms2.els_trkidx().at(elIndex) >= 0) && (cms2.els_trk_charge().at(elIndex) != cms2.trks_charge().at(cms2.els_trkidx().at(elIndex))) ) return true;
    if ((cms2.els_trkidx().at(elIndex) < 0)  && (cms2.els_trk_charge().at(elIndex) != cms2.els_sccharge().at(elIndex))) return true;

    return false;
}

bool isChargeFlip3agree(int elIndex){

    if (cms2.els_trkidx().at(elIndex) >= 0) 
    {
        // false if 3 charge measurements agree
        if(
                (cms2.els_trk_charge().at(elIndex)          // gsf
                 == cms2.trks_charge().at(cms2.els_trkidx().at(elIndex))) // ctf 
                &&
                (cms2.trks_charge().at(cms2.els_trkidx().at(elIndex)) // ctf 
                 == cms2.els_sccharge().at(elIndex)) )      // sc
            return false;  
    }

    return true;
}


//
// spike rejection for electrons
//

bool isSpikeElectron(const unsigned int index) {

    const int scidx = cms2.els_scindex()[index];
    bool isSpike = false;
    if (scidx != -1) {
        //subtract twice max since max is in both 1x3 and 3x1, and we want neither
        const float r4 = (cms2.scs_e1x3()[scidx] + cms2.scs_e3x1()[scidx] - 2*cms2.scs_eMax()[scidx])/cms2.scs_eMax()[scidx];
        if (r4 < 0.05) isSpike = true;
    }

    return isSpike;

}

//
// position correction for electrons
//

void electronCorrection_pos(const unsigned int index, float &dEtaIn, float &dPhiIn)
{

    //
    // uncorrected dEtaIn and dPhiIn
    //

    dEtaIn = cms2.els_dEtaIn()[index];
    dPhiIn = cms2.els_dPhiIn()[index];

    //
    // if configered not to apply correction
    // or in barrel or no valid super cluster
    // return uncorrected values
    //

    if (!(cms2.els_fiduciality()[index] & 1<<ISEE)) return;
    if (cms2.els_scindex()[index] == -1) return;

    //
    // set up correction parameters for EE+ and EE-
    // RecoEgamma/EgammaTools/python/correctedElectronsProducer_cfi.py?revision=1.2
    //

    //                                      X',     Y',     Z'
    float scPositionCorrectionEEP[3] = {   0.52,   -0.81,  0.81};
    float scPositionCorrectionEEM[3] = {    -0.02,  -0.81,  -0.94};

    LorentzVector initial_pos = cms2.scs_pos_p4()[cms2.els_scindex()[index]];
    LorentzVector corrected_pos;

    //
    // work out corrected position
    //

    if (cms2.els_etaSC()[index] < 0) {
        corrected_pos = LorentzVector(  initial_pos.x() + scPositionCorrectionEEM[0],
                initial_pos.y() + scPositionCorrectionEEM[1],
                initial_pos.z() + scPositionCorrectionEEM[2], 0.0);
    }
    if (cms2.els_etaSC()[index] > 0) {
        corrected_pos = LorentzVector(  initial_pos.x() + scPositionCorrectionEEP[0],
                initial_pos.y() + scPositionCorrectionEEP[1],
                initial_pos.z() + scPositionCorrectionEEP[2], 0.0);
    }

    //
    // work out correction to dEtaIn and dPhiIn
    //

    float deta_sc = corrected_pos.Eta() - initial_pos.Eta();
    float dphi_sc = acos(cos(corrected_pos.Phi() - initial_pos.Phi()));
    dEtaIn = deta_sc + cms2.els_dEtaIn()[index];
    dPhiIn = acos(cos(dphi_sc + cms2.els_dPhiIn()[index]));

}


double electron_d0PV(unsigned int index){
    if ( cms2.vtxs_sumpt().empty() ) return false;
    unsigned int iMax = 0;
    double sumPtMax = cms2.vtxs_sumpt().at(0);
    for ( unsigned int i = iMax+1; i < cms2.vtxs_sumpt().size(); ++i )
        if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = cms2.vtxs_sumpt().at(i);
        }
    double dxyPV = cms2.els_d0()[index]-
        cms2.vtxs_position()[iMax].x()*sin(cms2.els_trk_p4()[index].phi())+
        cms2.vtxs_position()[iMax].y()*cos(cms2.els_trk_p4()[index].phi());
    return dxyPV;
}


double electron_d0PV_smurfV3(unsigned int index){
  int vtxIndex = 0;
  double dxyPV = cms2.els_d0()[index]-
    cms2.davtxs_position()[vtxIndex].x()*sin(cms2.els_trk_p4()[index].phi())+
    cms2.davtxs_position()[vtxIndex].y()*cos(cms2.els_trk_p4()[index].phi());
  return dxyPV;
}

double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
  return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
}
double electron_dzPV_smurfV3(unsigned int index){
  int vtxIndex = 0;
  double dzpv = dzPV(cms2.els_vertex_p4()[index], cms2.els_trk_p4()[index], cms2.davtxs_position()[vtxIndex]);
  return dzpv;
}


double electron_dzPV_wwV1(unsigned int index){ 
    if ( cms2.vtxs_sumpt().empty() ) return 9999.;
    double sumPtMax = -1;
    int iMax = -1;
    for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
        if (!isGoodVertex(i)) continue;
        if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = cms2.vtxs_sumpt().at(i);
        }
    }
    if (iMax<0) return 9999.;

    const LorentzVector& vtx = cms2.els_vertex_p4()[index];
    const LorentzVector& p4 = cms2.els_trk_p4()[index];
    const LorentzVector& pv = cms2.vtxs_position()[iMax]; 
    return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt(); 
    /* directly from NtupleMacros/WW/doAnalysis.cc
       double dzpv = dzPV(cms2.els_vertex_p4()[index], cms2.els_trk_p4()[index], cms2.vtxs_position()[iMax]);
       double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
       return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
       }*/
}


double electron_d0PV_wwV1(unsigned int index){ 
    if ( cms2.vtxs_sumpt().empty() ) return 9999.;
    double sumPtMax = -1;
    int iMax = -1;
    for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
        if (!isGoodVertex(i)) continue;
        if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
            iMax = i;
            sumPtMax = cms2.vtxs_sumpt().at(i);
        }
    }
    if (iMax<0) return 9999.;
    double dxyPV = cms2.els_d0()[index]-
        cms2.vtxs_position()[iMax].x()*sin(cms2.els_trk_p4()[index].phi())+
        cms2.vtxs_position()[iMax].y()*cos(cms2.els_trk_p4()[index].phi());
    return dxyPV;
}

double electron_d0PV_mindz(unsigned int index){ 
    if ( cms2.vtxs_sumpt().empty() ) return 9999.;
    double minDz = 999.;
    int iMin = -1;
    for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
        if (!isGoodVertex(i)) continue;
      
      const LorentzVector& vtx = cms2.gsftrks_vertex_p4().at(cms2.els_gsftrkidx()[index]);
      const LorentzVector& p4 = cms2.gsftrks_p4().at(cms2.els_gsftrkidx()[index]);
      const LorentzVector& pv = cms2.vtxs_position()[i]; 
      float dz = fabs((vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt()); 
      
      if ( dz < minDz ){
	iMin = i;
	minDz = dz;
      }
    }
    if (iMin<0) return 9999.;
    double dxyPV = cms2.els_d0()[index]-
        cms2.vtxs_position()[iMin].x()*sin(cms2.els_trk_p4()[index].phi())+
        cms2.vtxs_position()[iMin].y()*cos(cms2.els_trk_p4()[index].phi());
    return dxyPV;
}

bool isFromConversionMIT(const unsigned int index){
  return isMITConversion(index, 0,   1e-6,   2.0,   true,  false);
}

electronIdComponent_t passLikelihoodId(unsigned int index, float lhValue, int workingPoint) {
  unsigned int answer = 0;
  float etaSC = cms2.els_etaSC().at(index);
  unsigned int nbrem = cms2.els_nSeed().at(index);
  if (workingPoint==95) 
    if (fabs(etaSC)<1.479&&nbrem==0&&lhValue>-4.274 || fabs(etaSC)<1.479&&nbrem>=1&&lhValue>-3.773 || fabs(etaSC)>1.479&&nbrem==0&&lhValue>-5.092 || fabs(etaSC)>1.479&&nbrem>=1&&lhValue>-2.796)
      answer |= (1<<ELEID_ID);
  else if (workingPoint==90)
    if (fabs(etaSC)<1.479&&nbrem==0&&lhValue>-1.497 || fabs(etaSC)<1.479&&nbrem>=1&&lhValue>-1.521 || fabs(etaSC)>1.479&&nbrem==0&&lhValue>-2.571 || fabs(etaSC)>1.479&&nbrem>=1&&lhValue>-0.657)
      answer |= (1<<ELEID_ID);
  else if (workingPoint==85)
    if (fabs(etaSC)<1.479&&nbrem==0&&lhValue>+0.163 || fabs(etaSC)<1.479&&nbrem>=1&&lhValue>+0.065 || fabs(etaSC)>1.479&&nbrem==0&&lhValue>-0.683 || fabs(etaSC)>1.479&&nbrem>=1&&lhValue>+1.564)
      answer |= (1<<ELEID_ID);
  else if (workingPoint==80)
    if (fabs(etaSC)<1.479&&nbrem==0&&lhValue>+1.193 || fabs(etaSC)<1.479&&nbrem>=1&&lhValue>+1.345 || fabs(etaSC)>1.479&&nbrem==0&&lhValue>+0.810 || fabs(etaSC)>1.479&&nbrem>=1&&lhValue>+3.021)
      answer |= (1<<ELEID_ID);
  else if (workingPoint==70)
    if (fabs(etaSC)<1.479&&nbrem==0&&lhValue>+1.781 || fabs(etaSC)<1.479&&nbrem>=1&&lhValue>+2.397 || fabs(etaSC)>1.479&&nbrem==0&&lhValue>+2.361 || fabs(etaSC)>1.479&&nbrem>=1&&lhValue>+4.052)
      answer |= (1<<ELEID_ID);
  else {
    cout << "Error! Likelihood WP not supported: " << workingPoint << ". Please choose 70, 80, 85, 90, 95" << endl;
  }
  return answer;
}
