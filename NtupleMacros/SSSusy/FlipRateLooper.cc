#include <math.h>
#include <sstream>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"
#include "Tools/chargeflip.h"
#include "Tools/fliprate_egun.h"

FlipRateLooper::FlipRateLooper (Sample s, cuts_t c, const char *fname) 
     : Looper(s, c, fname)
{
     // zero out the candidate counters (don't comment this out)
     memset(cands_passing_syst_hi, 0, sizeof(cands_passing_syst_hi));
     memset(cands_passing_syst_lo, 0, sizeof(cands_passing_syst_lo));
}

void FlipRateLooper::BookHistos 	()
{
//      fake_syst = new TH2F("fake_syst", "fake syst uncertainty;#eta;pt", 
// 			  fakeRate().GetNbinsX(), fakeRate().GetXaxis()->GetXbins()->GetArray(),
// 			  fakeRate().GetNbinsY(), fakeRate().GetYaxis()->GetXbins()->GetArray());

     Looper::BookHistos();


}

cuts_t FlipRateLooper::DilepSelect (int i_hyp)
{
     // this doesn't do anything special
     return Looper::DilepSelect(i_hyp);
}

void FlipRateLooper::FillDilepHistos (int i_hyp)
{
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     cuts_t cuts_passed = DilepSelect(i_hyp);
     
     if ((cuts_passed & cuts_) != cuts_)
	  return;

     // the logic is as follows:
     //

     bool usePGunFlips = true;
     double weight = Weight(i_hyp);
     
     double fr_ll = 0;
     switch (abs(cms2.hyp_ll_id()[i_hyp])) {
     case 11:
       if(usePGunFlips)
         fr_ll = getSingleEleFlipRate(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_ll_p4()[i_hyp].eta());
       else
         fr_ll = getZtoEEMCFlipRate(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_ll_p4()[i_hyp].eta());
       break;
     case 13:
       break;
     default:
       assert(0);
     }

     double fr_lt = 0;
     switch (abs(cms2.hyp_lt_id()[i_hyp])) {
     case 11:
       if(usePGunFlips)
         fr_lt = getSingleEleFlipRate(cms2.hyp_lt_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].eta());
       else
         fr_lt = getZtoEEMCFlipRate(cms2.hyp_lt_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].eta());
       break;
     case 13:
       break;
     default:
       assert(0);
     }
     // a la Claudio talk
     weight *= fr_lt * (1 - fr_ll) + fr_ll * (1 - fr_lt);
     // a la Frank (taking SS in OS into account)
     //     weight *= fr_lt / (1 - fr_lt) + fr_ll / (1 - fr_ll);

     cands_passing_[myType] += weight;
     cands_passing_w2_[myType] += weight * weight;
     cands_count_[myType]++;
     cands_passing_[DILEPTON_ALL] += weight;
     cands_passing_w2_[DILEPTON_ALL] += weight * weight;
     cands_count_[DILEPTON_ALL]++;

     //     hnJet->Fill(cuts_passed, myType, cms2.hyp_njets()[i_hyp], weight); //hyp jets
     hnJet->Fill(cuts_passed, myType, caloJets.size() , weight); // change to plotting caloJets size
     hnHyp->Fill(cuts_passed, myType, cms2.hyp_p4().size(), weight);
     hdilMass->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].mass(), weight);
     hmet->Fill(cuts_passed, myType, cms2.evt_tcmet(), weight);

     //     hnJet3D_[myType]->Fill(cms2.hyp_njets()[i_hyp],eta,pt, weight * err);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
       helPt->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
       helEta->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
       heleRelIso->Fill(cuts_passed, myType, inv_el_relsusy_iso(cms2.hyp_lt_index()[i_hyp], true), weight);
       helPdgId->Fill(cuts_passed, myType, abs(cms2.els_mc_id()[ cms2.hyp_lt_index()[i_hyp] ]), weight);
       helMoPdgId->Fill(cuts_passed, myType, abs(cms2.els_mc_motherid()[ cms2.hyp_lt_index()[i_hyp] ]), weight);
//        const double err = elFakeProb(cms2.hyp_lt_index()[i_hyp], 1) -
//          elFakeProb(cms2.hyp_lt_index()[i_hyp], 0);
//        helPt3D_[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(),eta,pt, weight * err);
//        helEta3D_[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(),eta,pt, weight * err);
//      //     hmet3D_[myType]->Fill(cms2.evt_tcmet(),eta,pt, weight * err);
       if(
          (abs(cms2.els_mc_id()[ cms2.hyp_lt_index()[i_hyp] ])==11 && abs(cms2.els_mc_motherid()[ cms2.hyp_lt_index()[i_hyp] ]) == 22) ||
          (abs(cms2.els_mc_id()[ cms2.hyp_lt_index()[i_hyp] ])==22)                                                                    ||
          (abs(cms2.els_mc_id()[ cms2.hyp_lt_index()[i_hyp] ]) > 100 && (abs(cms2.els_mc_id()[ cms2.hyp_lt_index()[i_hyp] ]) < 200))
          ) {
         helPdgIdCat->Fill(cuts_passed, myType, 1, weight);
       }
       else if((abs(cms2.els_mc_id()[ cms2.hyp_lt_index()[i_hyp] ]) > 200 && (abs(cms2.els_mc_id()[ cms2.hyp_lt_index()[i_hyp] ]) < 400))){
         helPdgIdCat->Fill(cuts_passed, myType, 2, weight);
       }
       else if((abs(cms2.els_mc_id()[ cms2.hyp_lt_index()[i_hyp] ]) == 11 && abs(cms2.els_mc_motherid()[ cms2.hyp_lt_index()[i_hyp] ]) >=400 )){
         helPdgIdCat->Fill(cuts_passed, myType, 3, weight);
       }
       else {
         helPdgIdCat->Fill(cuts_passed, myType, 4, weight);
       }
        // 	  heleRelIsoTrk->Fill(cuts_passed, myType, reliso_lt(i_hyp, false), weight);
     } else {
       //        hmuPt->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
       //        hmuEta->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
       helPt->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
       helEta->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
       heleRelIso->Fill(cuts_passed, myType, inv_el_relsusy_iso(cms2.hyp_ll_index()[i_hyp], true), weight);
       helPdgId->Fill(cuts_passed, myType, abs(cms2.els_mc_id()[ cms2.hyp_ll_index()[i_hyp] ]), weight);
       helMoPdgId->Fill(cuts_passed, myType, abs(cms2.els_mc_motherid()[ cms2.hyp_ll_index()[i_hyp] ]), weight);
//        const double err = elFakeProb(cms2.hyp_ll_index()[i_hyp], 1) -
//          elFakeProb(cms2.hyp_ll_index()[i_hyp], 0);
//        helPt3D_[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(),eta,pt, weight * err);
//        helEta3D_[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(),eta,pt, weight * err);
//      //     hmet3D_[myType]->Fill(cms2.evt_tcmet(),eta,pt, weight * err);
       if(
          (abs(cms2.els_mc_id()[ cms2.hyp_ll_index()[i_hyp] ])==11 && abs(cms2.els_mc_motherid()[ cms2.hyp_ll_index()[i_hyp] ]) == 22) ||
          (abs(cms2.els_mc_id()[ cms2.hyp_ll_index()[i_hyp] ])==22)                                                                    ||
          (abs(cms2.els_mc_id()[ cms2.hyp_ll_index()[i_hyp] ]) > 100 && (abs(cms2.els_mc_id()[ cms2.hyp_ll_index()[i_hyp] ]) < 200))
          ) {
         helPdgIdCat->Fill(cuts_passed, myType, 1, weight);
       }
       else if((abs(cms2.els_mc_id()[ cms2.hyp_ll_index()[i_hyp] ]) > 200 && (abs(cms2.els_mc_id()[ cms2.hyp_ll_index()[i_hyp] ]) < 400))){
         helPdgIdCat->Fill(cuts_passed, myType, 2, weight);
       }
       else if((abs(cms2.els_mc_id()[ cms2.hyp_ll_index()[i_hyp] ]) == 11 && abs(cms2.els_mc_motherid()[ cms2.hyp_ll_index()[i_hyp] ]) >=400 )){
         helPdgIdCat->Fill(cuts_passed, myType, 3, weight);
       }
       else {
         helPdgIdCat->Fill(cuts_passed, myType, 4, weight);
       }
       // 	  heleRelIsoTrk->Fill(cuts_passed, myType, reliso_ll(i_hyp, false), weight);
     } else {
       //        hmuPt->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
       //        hmuEta->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
     }

}

