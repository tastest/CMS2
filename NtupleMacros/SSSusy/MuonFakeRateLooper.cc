#include <math.h>
#include <sstream>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Tools/fakerates.h"
#include "Looper.h"

MuonFakeRateLooper::MuonFakeRateLooper (Sample s, cuts_t c, const char *fname) 
     : Looper(s, c, fname)
{
     // zero out the candidate counters (don't comment this out)
     memset(cands_passing_syst_hi, 0, sizeof(cands_passing_syst_hi));
     memset(cands_passing_syst_lo, 0, sizeof(cands_passing_syst_lo));
}

void MuonFakeRateLooper::BookHistos 	()
{
     fake_syst = new TH2F("fake_syst", "fake syst uncertainty;#eta;pt", 
			  fakeRate().GetNbinsX(), fakeRate().GetXaxis()->GetXbins()->GetArray(),
			  fakeRate().GetNbinsY(), fakeRate().GetYaxis()->GetXbins()->GetArray());

     Looper::BookHistos();

//      for (unsigned int bucket = 0;
//           bucket < 4;
//           ++bucket ) {
//        hnJet3D_[bucket] = book3DVarHist(Form("%s_%s_%s",sample_.name.c_str(),"nJet3D",dilepton_hypo_names[bucket]),
//                                         Form("%s_%s_%s",sample_.name.c_str(),"nJet3D",dilepton_hypo_names[bucket]),
//                                         hnJet->GetNBinsX(), hnJet->GetXmin(), hnJet->GetXmax(),
//                                         fakeRate().GetNbinsX(), fakeRate().GetXaxis()->GetXbins()->GetArray(),
//                                         fakeRate().GetNbinsY(), fakeRate().GetYaxis()->GetYbins()->GetArray(),
//                                         "n_{jet}","#eta", "p_{T} [GeV]", sample_.histo_color);

//        hmuPt3D_[bucket] = book3DVarHist(Form("%s_%s_%s",sample_.name.c_str(),"muPt3D",dilepton_hypo_names[bucket]),
//                                         Form("%s_%s_%s",sample_.name.c_str(),"muPt3D",dilepton_hypo_names[bucket]),
//                                         hmuPt->GetNBinsX(), hmuPt->GetXmin(), hmuPt->GetXmax(),
//                                         fakeRate().GetNbinsX(), fakeRate().GetXaxis()->GetXbins()->GetArray(),
//                                         fakeRate().GetNbinsY(), fakeRate().GetYaxis()->GetYbins()->GetArray(),
// //                                         ptNBins,ptBins,
// //                                         fakeXNBins,fakeXBins,
// //                                         fakeYNBins,fakeYBins,
//                                         "p_{T}^{e} [GeV]","#eta", "p_{T} [GeV]",sample_.histo_color);
//        hmuEta3D_[bucket] = book3DVarHist(Form("%s_%s_%s",sample_.name.c_str(),"muEta3D",dilepton_hypo_names[bucket]),
//                                          Form("%s_%s_%s",sample_.name.c_str(),"muEta3D",dilepton_hypo_names[bucket]),
//                                         hmuEta->GetNBinsX(), hmuEta->GetXmin(), hmuEta->GetXmax(),
//                                         fakeRate().GetNbinsX(), fakeRate().GetXaxis()->GetXbins()->GetArray(),
//                                         fakeRate().GetNbinsY(), fakeRate().GetYaxis()->GetYbins()->GetArray(),
// //                                          etaNBins,etaBins,
// //                                          fakeXNBins,fakeXBins,
// //                                          fakeYNBins,fakeYBins,
//                                          "#eta^{e} [GeV]","#eta", "p_{T} [GeV]",sample_.histo_color);
//        //        hmet3D_[bucket] = book3DVarHist(Form("%s_%s_%s",sample_.name.c_str(),"met3D",dilepton_hypo_names[bucket]),
//        //                                        Form("%s_%s_%s",sample_.name.c_str(),"met3D",dilepton_hypo_names[bucket]),
//        //                                        metNBins,metBins,
//        //                                        fakeXNBins,fakeXBins,
//        //                                        fakeYNBins,fakeYBins,
//        //                                        "MET [GeV]","#eta", "p_{T} [GeV]",sample_.histo_color);
//      }

}

cuts_t MuonFakeRateLooper::DilepSelect (int i_hyp)
{
     // this doesn't do anything special
     return Looper::DilepSelect(i_hyp);
}

void MuonFakeRateLooper::FillDilepHistos (int i_hyp)
{
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     cuts_t cuts_passed = DilepSelect(i_hyp);
     
     if ((cuts_passed & cuts_) != cuts_)
	  return;

     // the logic is as follows:
     //
     // first, we assume the lt to be real (by testing whether it
     // passes the reco+truth tags) and calculate the probability that the
     // ll is fake.
     //
     // then we switch it around and assume the ll to be real (by
     // testing whether it passes the reco+truth tags) and calculate the probability that the
     // lt is fake.
     //
     // this double counts double fakes (but our W+jets shouldn't contain a
     // lot of those anyway)
     //
     // it also neglects the horrible things that will happen if the
     // trigger cuts are tighter than the denominator cuts 

     double weight = Weight(i_hyp);
     
     double fr_ll = 0;
     // muon truth&ID/iso tags are already in global cuts
     // 	  const double eta = cms2.hyp_ll_p4()[i_hyp].eta();
     // 	  const double pt = cms2.hyp_ll_p4()[i_hyp].pt();
     switch (abs(cms2.hyp_ll_id()[i_hyp])) {
     case 11:
       break;
     case 13:
       if (isFakeableMuon(cms2.hyp_ll_index()[i_hyp]) &&
           not isNumeratorMuon(cms2.hyp_ll_index()[i_hyp]))
         fr_ll = muFakeProb(cms2.hyp_ll_index()[i_hyp], 0);
       break;
       // 	       if (isFakeableMuon(cms2.hyp_ll_index()[i_hyp]) &&
       // 		   not isNumeratorMuon(cms2.hyp_ll_index()[i_hyp]))
       // 		    fr_ll = muFakeProb(cms2.hyp_ll_index()[i_hyp], 0);
       break;
     default:
       assert(0);
     }

     double fr_lt = 0;
     // 	  const double eta = cms2.hyp_lt_p4()[i_hyp].eta();
     // 	  const double pt = cms2.hyp_lt_p4()[i_hyp].pt();
     switch (abs(cms2.hyp_lt_id()[i_hyp])) {
     case 11:
       break;
     case 13:
       if (isFakeableMuon(cms2.hyp_lt_index()[i_hyp]) &&
           not isNumeratorMuon(cms2.hyp_lt_index()[i_hyp]))
         fr_lt = muFakeProb(cms2.hyp_lt_index()[i_hyp], 0);
       // 	       if (isFakeableMuon(cms2.hyp_lt_index()[i_hyp]) &&
       // 		   not isNumeratorMuon(cms2.hyp_lt_index()[i_hyp]))
       // 		    fr_lt = muFakeProb(cms2.hyp_lt_index()[i_hyp], 0);
       break;
     default:
       assert(0);
     }

     weight *= fr_lt / (1 - fr_lt) + fr_ll / (1 - fr_ll);

     // just to make sure: only one of fr_lt and fr_ll should ever be non-zero
     // printf("fr_lt = %g\tfr_ll = %g\n", fr_lt, fr_ll);
     assert(not ((fr_lt != 0) && (fr_ll != 0)));

     cands_passing_[myType] += weight;
     cands_passing_w2_[myType] += weight * weight;
     cands_count_[myType]++;
     cands_passing_[DILEPTON_ALL] += weight;
     cands_passing_w2_[DILEPTON_ALL] += weight * weight;
     cands_count_[DILEPTON_ALL]++;


     hnJet->Fill(cuts_passed, myType, cms2.hyp_njets()[i_hyp], weight);
     //     hnJet3D_[myType]->Fill(cms2.hyp_njets()[i_hyp],eta,pt, weight * err);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
       hmuPt->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
       hmuEta->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
       hmuRelIso->Fill(cuts_passed, myType, inv_mu_relsusy_iso(cms2.hyp_lt_index()[i_hyp]), weight);
       hmuPdgId->Fill(cuts_passed, myType, abs(cms2.mus_mc_id()[ cms2.hyp_lt_index()[i_hyp] ]), weight);
       hmuMoPdgId->Fill(cuts_passed, myType, abs(cms2.mus_mc_motherid()[ cms2.hyp_lt_index()[i_hyp] ]), weight);
//        const double err = muFakeProb(cms2.hyp_lt_index()[i_hyp], 1) -
//          muFakeProb(cms2.hyp_lt_index()[i_hyp], 0);
//        hmuPt3D_[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(),eta,pt, weight * err);
//        hmuEta3D_[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(),eta,pt, weight * err);
//      //     hmet3D_[myType]->Fill(cms2.evt_tcmet(),eta,pt, weight * err);
       if( // punchthrough / sailthrough
          (abs(cms2.mus_mc_id()[ cms2.hyp_lt_index()[i_hyp]  ]) != 13 )
          ) {
         hmuPdgIdCat->Fill(cuts_passed, myType, 1, weight);
       }
       else if( (abs(cms2.mus_mc_id()[ cms2.hyp_lt_index()[i_hyp]  ]) == 13) && abs(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]]) < 400 ){
         hmuPdgIdCat->Fill(cuts_passed, myType, 2, weight);
       }
       else if(abs(cms2.mus_mc_id()[ cms2.hyp_lt_index()[i_hyp]  ]) == 13 && (idIsBeauty( abs(cms2.mus_mc_motherid()[ cms2.hyp_lt_index()[i_hyp] ])) || idIsCharm( abs(cms2.mus_mc_motherid()[ cms2.hyp_lt_index()[i_hyp] ]))) ){
         hmuPdgIdCat->Fill(cuts_passed, myType, 3, weight);
       }
       else {
         hmuPdgIdCat->Fill(cuts_passed, myType, 4, weight);
       }
        // 	  hmueRelIsoTrk->Fill(cuts_passed, myType, reliso_lt(i_hyp, false), weight);
     } else {
       //        hmuPt->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
       //        hmuEta->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
       hmuPt->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
       hmuEta->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
       hmuRelIso->Fill(cuts_passed, myType, inv_mu_relsusy_iso(cms2.hyp_ll_index()[i_hyp]), weight);
       hmuPdgId->Fill(cuts_passed, myType, abs(cms2.mus_mc_id()[ cms2.hyp_ll_index()[i_hyp] ]), weight);
       hmuMoPdgId->Fill(cuts_passed, myType, abs(cms2.mus_mc_motherid()[ cms2.hyp_ll_index()[i_hyp] ]), weight);
//        const double err = muFakeProb(cms2.hyp_ll_index()[i_hyp], 1) -
//          muFakeProb(cms2.hyp_ll_index()[i_hyp], 0);
//        hmuPt3D_[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(),eta,pt, weight * err);
//        hmuEta3D_[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(),eta,pt, weight * err);
//      //     hmet3D_[myType]->Fill(cms2.evt_tcmet(),eta,pt, weight * err);
       if( // punchthrough / sailthrough
          (abs(cms2.mus_mc_id()[ cms2.hyp_ll_index()[i_hyp]  ]) != 13 )
          ) {
         hmuPdgIdCat->Fill(cuts_passed, myType, 1, weight);
       }
       else if( (abs(cms2.mus_mc_id()[ cms2.hyp_ll_index()[i_hyp]  ]) == 13) && abs(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]]) < 400 ){
         hmuPdgIdCat->Fill(cuts_passed, myType, 2, weight);
       }
       else if(abs(cms2.mus_mc_id()[ cms2.hyp_ll_index()[i_hyp]  ]) == 13 && (idIsBeauty( abs(cms2.mus_mc_motherid()[ cms2.hyp_ll_index()[i_hyp] ])) || idIsCharm( abs(cms2.mus_mc_motherid()[ cms2.hyp_ll_index()[i_hyp] ]))) ){
         hmuPdgIdCat->Fill(cuts_passed, myType, 3, weight);
       }
       else {
         hmuPdgIdCat->Fill(cuts_passed, myType, 4, weight);
       }
       // 	  hmuRelIsoTrk->Fill(cuts_passed, myType, reliso_ll(i_hyp, false), weight);
     } else {
       //        hmuPt->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
       //        hmuEta->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
     }

}

// bool MuonFakeRateLooper::fillErrorInPrediction(TH1F* prediction,
//                                            TH3F* predictionError,
//                                            bool addStatisticalError) {
//   //
//   // calculate error for prediction from predictionError
//   //

//   for (int predictionBin = 0;
//        predictionBin <= prediction->GetNbinsX()+1;
//        ++predictionBin ) {
//     float err2 = 0;
//     for ( int fakeXBin = 0;
//           fakeXBin <= fakeRate().GetNbinsX()+1;
//           ++fakeXBin ) {
//       for ( int fakeYBin = 0;
//             fakeYBin <= fakeRate().GetNbinsY()+1;
//             ++fakeYBin ) {
//         err2 += predictionError->GetBinContent(predictionBin,fakeXBin,fakeYBin) *
//           predictionError->GetBinContent(predictionBin,fakeXBin,fakeYBin);
//       }
//     }
//     float err = 0;
//     if ( addStatisticalError ) {
//       err = sqrt( prediction->GetBinError(predictionBin) *
//                   prediction->GetBinError(predictionBin) +
//                   err2 );
//     } else {
//       err = sqrt(err2);
//     }

//     prediction->SetBinError(predictionBin,err);
//   }
//   return true;
// }


// void MuonFakeRateLooper::End ()
// {
//   //------------------------------------------------------------
//   //Example status message at the end of a looper; edit for your
//   //application
//   //------------------------------------------------------------

//   // treat errors correctly
//   for (unsigned int bucket = 0;
//        bucket < 4;
//        ++bucket ) {
//     fillErrorInPrediction(hnJet_[bucket],hnJet3D_[bucket]);
//     fillErrorInPrediction(hmuPt_[bucket],hmuPt3D_[bucket]);
//     fillErrorInPrediction(hmuEta_[bucket],hmuEta3D_[bucket]);
//     fillErrorInPrediction(hmet_[bucket],hmet3D_[bucket]);
//   }

//   // finish with std End
//   Looper::End();
// }
