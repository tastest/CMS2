#include <math.h>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Tools/fakerates.h"
#include "Looper.h"

Looper::Looper (Sample s, cuts_t c, const char *fname) 
  : LooperBase(s, c, fname)
{
  // zero out the candidate counters (don't comment this out)
  memset(cands_passing_	, 0, sizeof(cands_passing_ ));
  memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_ ));
  memset(cands_count_		, 0, sizeof(cands_count_ ));
}

void Looper::BookHistos ()
{

  // njet plots
  hnJet_ = new NMinus1Hist(sample_, "nJet;n_{jets};Dilepton Cand." , 6 , -0.5, 5 , cuts_, 0 );
  hnCaloJet_ = new NMinus1Hist(sample_, "nCaloJet;n_{jets};Dilepton Cand." , 6 , -0.5, 5 , cuts_, 0 );
  hnTrackJet_ = new NMinus1Hist(sample_, "nTrackJet;n_{jets};Dilepton Cand." , 6 , -0.5, 5 , cuts_, 0 );

  // general lepton plots

  std::vector<cuts_t> lt_cuts;
  lt_cuts.push_back((CUT_BIT(CUT_LT_PT)));
  lt_cuts.push_back((CUT_BIT(CUT_LT_GOOD))); 
  lt_cuts.push_back((CUT_BIT(CUT_LT_CALOISO)));

  std::vector<std::string> lt_cut_names;
  lt_cut_names.push_back("-Cut-On-PT");
  lt_cut_names.push_back("-Cut-On-GOOD");
  lt_cut_names.push_back("-Cut-On-CALOISO");

  hltPt_ = new NMinus1Hist(sample_, "ltPt;p_{T} [GeV];Dilepton Cand." , 150 , 0, 150 , cuts_, lt_cuts, lt_cut_names);
  hltEta_ = new NMinus1Hist(sample_, "ltEta;#eta;Dilepton Cand." , 12 , -3, 3 , cuts_, lt_cuts, lt_cut_names);
  hltCaloIso_ = new NMinus1Hist(sample_, "ltCaloIso;calo iso;Dilepton Cand." , 50 , 0., 1., cuts_, lt_cuts, lt_cut_names);

  std::vector<cuts_t> ll_cuts;
  ll_cuts.push_back((CUT_BIT(CUT_LL_PT)));
  ll_cuts.push_back((CUT_BIT(CUT_LL_GOOD))); 
  ll_cuts.push_back((CUT_BIT(CUT_LL_CALOISO)));

  std::vector<std::string> ll_cut_names;
  ll_cut_names.push_back("-Cut-On-PT");
  ll_cut_names.push_back("-Cut-On-GOOD");
  ll_cut_names.push_back("-Cut-On-CALOISO");

  hllPt_ = new NMinus1Hist(sample_, "llPt;p_{T} [GeV];Dilepton Cand." , 75 , 0, 150 , cuts_, ll_cuts, ll_cut_names);
  hllEta_ = new NMinus1Hist(sample_, "llEta;#eta;Dilepton Cand." , 12 , -3, 3 , cuts_, ll_cuts, ll_cut_names);
  hllCaloIso_ = new NMinus1Hist(sample_, "llCaloIso;calo iso;Dilepton Cand." , 50 , 0., 1., cuts_, ll_cuts, ll_cut_names);

  std::vector<cuts_t> lepton_cuts;
  lepton_cuts.push_back((CUT_BIT(CUT_LT_PT)));
  lepton_cuts.push_back((CUT_BIT(CUT_LT_GOOD))); 
  lepton_cuts.push_back((CUT_BIT(CUT_LT_CALOISO)));
  lepton_cuts.push_back((CUT_BIT(CUT_LL_PT)));
  lepton_cuts.push_back((CUT_BIT(CUT_LL_GOOD))); 
  lepton_cuts.push_back((CUT_BIT(CUT_LL_CALOISO)));

  std::vector<std::string> lepton_cut_names;
  lepton_cut_names.push_back("-Cut-On-lt-PT");
  lepton_cut_names.push_back("-Cut-On-lt-GOOD");
  lepton_cut_names.push_back("-Cut-On-lt-CALOISO");
  lepton_cut_names.push_back("-Cut-On-ll-PT");
  lepton_cut_names.push_back("-Cut-On-ll-GOOD");
  lepton_cut_names.push_back("-Cut-On-ll-CALOISO");

  helPt_ = new NMinus1Hist(sample_, "elPt;p_{T} [GeV];Dilepton Cand." , 75 , 0, 150 , cuts_, lepton_cuts, lepton_cut_names);
  hmuPt_ = new NMinus1Hist(sample_, "muPt;p_{T} [GeV];Dilepton Cand." , 75 , 0, 150 , cuts_, lepton_cuts, lepton_cut_names);
  helEta_ = new NMinus1Hist(sample_, "elEta;#eta;Dilepton Cand." , 12 , -3, 3 , cuts_, lepton_cuts, lepton_cut_names);
  hmuEta_ = new NMinus1Hist(sample_, "muEta;#eta;Dilepton Cand." , 12 , -3, 3 , cuts_, lepton_cuts, lepton_cut_names);
  helCaloIso_ = new NMinus1Hist(sample_, "elCaloIso;calo iso;Dilepton Cand." , 50 , 0., 1. , cuts_, lepton_cuts, lepton_cut_names);
  hmuCaloIso_ = new NMinus1Hist(sample_, "muCaloIso;calo iso;Dilepton Cand." , 50 , 0., 1. , cuts_, lepton_cuts, lepton_cut_names);

  // track plots
  hnTrack_ = new NMinus1Hist(sample_, "nTrack;n_{track};Dilepton Cand." , 200 , 0., 200. , cuts_, (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS)) );

  // dilepton plots
  hdphiLep_ = new NMinus1Hist(sample_, "dphiLep;#Delta#Phi;Dilepton Cand." , 50 , 0, M_PI , cuts_, 0 );
  hdilMass_ = new NMinus1Hist(sample_, "dilMass;inv. mass [GeV];Dilepton Cand." , 100 , 0, 300 , cuts_, 0 );
  hdilPt_ = new NMinus1Hist(sample_, "dilPt;p_{T} [GeV];Dilepton Cand." , 100 , 0, 300 , cuts_, 0 );

  // MET plots
  std::vector<cuts_t> met_cuts;
  met_cuts.push_back((CUT_BIT(CUT_PASS2_METCORR)));
  met_cuts.push_back((CUT_BIT(CUT_PASS4_METCORR)));
  met_cuts.push_back((CUT_BIT(CUT_PASS2_METCORR)) | (CUT_BIT(CUT_PASS4_METCORR)));
  met_cuts.push_back((CUT_BIT(CUT_PASS_METCORR_10)));
  met_cuts.push_back((CUT_BIT(CUT_PASS_METCORR_1)));

  std::vector<std::string> met_cut_names;
  met_cut_names.push_back("-Cut-On-PASS2-MET");
  met_cut_names.push_back("-Cut-On-PASS4-MET");
  met_cut_names.push_back("-Cut-On-PASS2-AND-PASS4-MET");
  met_cut_names.push_back("-Cut-On-MET-10");
  met_cut_names.push_back("-Cut-On-MET-1");

  hmet_ = new NMinus1Hist(sample_, "met;MET [GeV];Dilepton Cand." , 100 , 0, 200 , cuts_, met_cuts, met_cut_names );
  hmetSpec_ = new NMinus1Hist(sample_, "metSpec;MET [GeV];Dilepton Cand." , 100 , 0, 200 , cuts_, met_cuts, met_cut_names );
  hmetTrkCorr_ = new NMinus1Hist(sample_, "metTrkCorr;MET [GeV];Dilepton Cand." , 100 , 0, 200 , cuts_, met_cuts, met_cut_names );

  // SumET plots
  
  std::vector<cuts_t> sumet_cuts;
  sumet_cuts.push_back((CUT_BIT(CUT_PASS_SUMET_10)));
  sumet_cuts.push_back((CUT_BIT(CUT_PASS_SUMET_1)));

  std::vector<std::string> sumet_cut_names;
  sumet_cut_names.push_back("-Cut-On-SUMET-10");
  sumet_cut_names.push_back("-Cut-On-SUMET-1");

  hsumet_ = new NMinus1Hist(sample_, "sum-et;#Sum E_{T} [GeV];Dilepton Cand." , 300 , 0, 1500 , cuts_, sumet_cuts, sumet_cut_names );

  hmeff_ = new NMinus1Hist(sample_, "nMeff;M_{eff} [GeV];Dilepton Cand." , 100 , 0., 4000 , cuts_, 0 );
  hmeffcorr_ = new NMinus1Hist(sample_, "nMeffCorr;M_{eff, tcCorr} [GeV];Dilepton Cand." , 100 , 0., 4000 , cuts_, 0 );

  // Thrust plots
  hthrust_ = new NMinus1Hist(sample_, "thrust;Thrust;Dilepton Cand." , 150 , 0.4, 1.1 , cuts_, 0 );
  hjptthrust_ = new NMinus1Hist(sample_, "jptthrust;Jptthrust;Dilepton Cand." , 150 , 0.4, 1.1 , cuts_, 0 );

  pdg = new TDatabasePDG();
}


bool Looper::FilterEvent()
{ 

  //
  // duplicate filter, based on trk information and dilepton hyp
  //
  // comment in following lines
  // 

  if (cms2.trks_d0().size() == 0)
    return true;
  DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.trks_d0()[0], 
			      cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
  if (is_duplicate(id)) {
    duplicates_total_n_++;
    duplicates_total_weight_ += cms2.evt_scale1fb();
    //     cout << "Filtered duplicate run: " << cms2.evt_run() << " event: " << cms2.evt_event() << endl;
    return true;
  }

  return false; 
}



cuts_t Looper::EventSelect ()
{
  //------------------------------------------------------------
  // In an event-based analysis, you would make your cuts here
  //------------------------------------------------------------

  cuts_t ret = 0;
  return ret;
}

cuts_t Looper::DilepSelect (int i_hyp)
{
  //------------------------------------------------------------
  // Example dilepton cuts; edit for your application
  //------------------------------------------------------------

  // cuts are failed until proven otherwise
  cuts_t ret = 0;
  const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);

  // temporarily calculate ONE global uncorrected jetSumEt:
  // prepare uncorrected jet collection
  vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > calo_jets_p4(cms2.jets_p4());
  sumEt_ = 0.;

  for (unsigned int icalojet=0; icalojet<calo_jets_p4.size(); ++icalojet) {
    // remove electron jets:
    if ((TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],calo_jets_p4[icalojet]) < 0.4)||
	(TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],calo_jets_p4[icalojet]) < 0.4)
	) continue;
    TLorentzVector p(calo_jets_p4[icalojet].Px()*cms2.jets_pat_noCorrF()[icalojet],
		     calo_jets_p4[icalojet].Py()*cms2.jets_pat_noCorrF()[icalojet],
		     calo_jets_p4[icalojet].Pz()*cms2.jets_pat_noCorrF()[icalojet],
		     calo_jets_p4[icalojet].E()*cms2.jets_pat_noCorrF()[icalojet]);//p is now uncorrected jet energy
    //       if (p.Perp() < jetet) continue;
    //       if (TMath::Abs(p.Eta()) > jeteta) continue;
    //        if (p.Perp() < 28 ) cout << p.Perp() << endl;
    if (TMath::Abs(p.Eta()) > 3.0 ) continue;
    if (p.Perp() < 15.) continue;
    //                 calo_jets->push_back(p);
    sumEt_ += p.Perp();
  }

  thrust = -999;
  jptthrust = -999;

  if(42 != 42) { // disabled 090702 for speed
  // calculate thrust (IBL variant)
  // ********** Claculate the thrust for the given event ****************************************** - begin >>>
  // *** initialize thrust value:
  //    thrust_mu1_ptrel =  50. * ( gRandom->Rndm(eventnr) );
  //    thrust_mu2_ptrel =  50. * ( gRandom->Rndm(eventnr) );

  // *** initialize thrust axis:
  thrust_axis.SetXYZ( 1., 1., 1.);
  thrust_axis_phi   = thrust_axis.Phi();
  thrust_axis_theta = thrust_axis.Theta();
    
  // *** find axis and calculate thrust:
  find_thrust_axis( &thrust_axis_phi, &thrust_axis_theta,  &thrust );
    
  // *** give found thrust theta and phi to thrust axis:
  thrust_axis.SetXYZ(TMath::Sin(thrust_axis_theta)*TMath::Cos(thrust_axis_phi) , TMath::Sin(thrust_axis_phi)*TMath::Sin(thrust_axis_theta) ,  TMath::Cos(thrust_axis_theta) );

  //  cout<<"Final thrust is: "<<thrust<<" wth eta / phi: "<<thrust_axis.Eta()<<" / "<<thrust_axis.Phi()<<endl;
  //  cout<<"Final thrust is: "<<thrust<<" wth theta / phi: "<<thrust_axis_theta<<" / "<<thrust_axis_phi<<endl;

  //      thrust_axis.SetMag(1.);
  //      thrust_axis.SetPhi(   thrust_axis_phi   );
  //      thrust_axis.SetTheta( thrust_axis_theta );
  // ********** Calculate the thrust for the given event ****************************************** - <<< end
  }

  if(42 != 42) { // disabled 090702 for speed
  // calculate jptthrust (IBL variant)
  // ********** Claculate the jptthrust for the given event ****************************************** - begin >>>
  // *** initialize jptthrust value:

  //    jptthrust_mu1_ptrel =  50. * ( gRandom->Rndm(eventnr) );
  //    jptthrust_mu2_ptrel =  50. * ( gRandom->Rndm(eventnr) );

  // *** initialize jptthrust axis:
  jptthrust_axis.SetXYZ( 1., 1., 1.);
  jptthrust_axis_phi   = jptthrust_axis.Phi();
  jptthrust_axis_theta = jptthrust_axis.Theta();
    
  // *** find axis and calculate jptthrust:
  find_thrust_axis( &jptthrust_axis_phi, &jptthrust_axis_theta,  &jptthrust, "JPT" );
    
  // *** give found jptthrust theta and phi to jptthrust axis:
  jptthrust_axis.SetXYZ(TMath::Sin(jptthrust_axis_theta)*TMath::Cos(jptthrust_axis_phi) , TMath::Sin(jptthrust_axis_phi)*TMath::Sin(jptthrust_axis_theta) ,  TMath::Cos(jptthrust_axis_theta) );

  //  cout<<"Final jptthrust is: "<<jptthrust<<" wth eta / phi: "<<jptthrust_axis.Eta()<<" / "<<jptthrust_axis.Phi()<<endl;
  //  cout<<"Final jptthrust is: "<<jptthrust<<" wth theta / phi: "<<jptthrust_axis_theta<<" / "<<jptthrust_axis_phi<<endl;

  //      jptthrust_axis.SetMag(1.);
  //      jptthrust_axis.SetPhi(   jptthrust_axis_phi   );
  //      jptthrust_axis.SetTheta( jptthrust_axis_theta );
  // ********** Calculate the jptthrust for the given event ****************************************** - <<< end
  }

  // enough tracks?
  if (cms2.trks_trk_p4().size() > 2)
    ret |= CUT_BIT(CUT_MORE_THAN_TWO_TRACKS);

  if(42 != 42) {
    // all these cuts included in hyp_jets collection
    // already. Just have to require that it
    // contains >0 entries.
//     // leading jet pt cut
//     if (cms2.jets_p4()[0].pt() > 30.)
//       ret |= CUT_BIT(CUT_LJET_ET);
    
//     // leading jet eta cut
//     if (TMath::Abs(cms2.jets_p4()[0].eta()) < 2.4)
//       ret |= CUT_BIT(CUT_LJET_ETA);
    
//     // leading jet emfrac cut
//     if (TMath::Abs(cms2.jets_emFrac[0]) > 0.1)
//       ret |= CUT_BIT(CUT_LJET_EMF);
    
//     // jet-lepton cleaning is also done on ntuple level
  }

  // jet cuts
  if (cms2.hyp_jets_p4()[i_hyp].size() > 0) {
    //    if (cms2.hyp_jets_p4()[i_hyp][0].pt() > 100.) // TEMP cut on 100 GeV jet
    if (cms2.hyp_jets_p4()[i_hyp][0].pt() > 30.) // TEMP cut on 30 GeV leading jet
      ret |= CUT_BIT(CUT_LJETS);
  }

  // pt cuts
  if (cms2.hyp_lt_p4()[i_hyp].pt() > 10.0) 
    ret |= (CUT_BIT(CUT_LT_PT));
  if (cms2.hyp_ll_p4()[i_hyp].pt() > 10.0) 
    ret |= (CUT_BIT(CUT_LL_PT));

  // sign cuts
  if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 ) 
    ret |= (CUT_BIT(CUT_OPP_SIGN));
  if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) 
    ret |= (CUT_BIT(CUT_SAME_SIGN));

//   // track corrected MET
//   const TVector3 trkCorr = correctMETforTracks();
//   if (pass4Met(i_hyp, trkCorr))
//     ret |= (CUT_BIT(CUT_PASS4_METCORR));
//   if (pass2Met(i_hyp, trkCorr))
//     ret |= (CUT_BIT(CUT_PASS2_METCORR));

//   // MET
//   if (pass4Met(i_hyp, TVector3()))
//     ret |= (CUT_BIT(CUT_PASS4_MET));
//   if (pass2Met(i_hyp, TVector3()))
//     ret |= (CUT_BIT(CUT_PASS2_MET));

  // track corrected MET
  const TVector3 trkCorr;
  // cheat a bit = make pass2 and pass4 plain tcMET cuts
  if (metSimple(80, trkCorr))
    ret |= (CUT_BIT(CUT_PASS4_METCORR));
  if (metSimple(80, trkCorr))
    ret |= (CUT_BIT(CUT_PASS2_METCORR));

  // MET
  if (metSimple(80, trkCorr))
    ret |= (CUT_BIT(CUT_PASS4_MET));
  if (metSimple(80, trkCorr))
    ret |= (CUT_BIT(CUT_PASS2_MET));


  // SameSigns grid cuts
  if (met10(i_hyp, TVector3())) 
    ret |= (CUT_BIT(CUT_PASS_MET_10));
  if (met1(i_hyp, TVector3()))
    ret |= (CUT_BIT(CUT_PASS_MET_1));
  if (met10(i_hyp, trkCorr)) 
    ret |= (CUT_BIT(CUT_PASS_METCORR_10));
  if (met1(i_hyp, trkCorr))
    ret |= (CUT_BIT(CUT_PASS_METCORR_1));

  if (sumEt10( sumEt_ )) 
    ret |= (CUT_BIT(CUT_PASS_SUMET_10));
  if (sumEt1( sumEt_ ))
    ret |= (CUT_BIT(CUT_PASS_SUMET_1));

  // muon quality
//   if(SameSigns_STD) {
//    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ) 
//      ret |= CUT_BIT(CUT_LT_GOOD);
//    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ) 
//      ret |= CUT_BIT(CUT_LL_GOOD);
//   // muID from Susy stuff
//   }
//  else {
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && GoodSusyLeptonID(13,cms2.hyp_lt_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_LT_GOOD);
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && GoodSusyLeptonID(13,cms2.hyp_ll_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_LL_GOOD);
  //  }
  // electron quality
  //   if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) )
  //     ret |= CUT_BIT(CUT_LT_GOOD);
  //   if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) )
  //     ret |= CUT_BIT(CUT_LL_GOOD);
  // eleID from Susy stuff
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && GoodSusyLeptonID(11,cms2.hyp_lt_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_LT_GOOD);
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && GoodSusyLeptonID(11,cms2.hyp_ll_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_LL_GOOD);
  
  // calo iso
  //   if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && passMuonIsolation(cms2.hyp_lt_index()[i_hyp]) ) {
  //     ret |= (CUT_BIT(CUT_LT_CALOISO));
  //   }
  //   if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && passMuonIsolation(cms2.hyp_ll_index()[i_hyp]) ) {
  //     ret |= (CUT_BIT(CUT_LL_CALOISO));
  //   }
  //   if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && passElectronIsolation(cms2.hyp_lt_index()[i_hyp], true)) {
  //     ret |= CUT_BIT(CUT_LT_CALOISO);
  //   }
  //   if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && passElectronIsolation(cms2.hyp_ll_index()[i_hyp], true)) {
  //     ret |= CUT_BIT(CUT_LL_CALOISO);
  //   } 
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && GoodSusyLeptonIsolation(13, cms2.hyp_lt_index()[i_hyp]) ) {
    ret |= (CUT_BIT(CUT_LT_CALOISO));
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && GoodSusyLeptonIsolation(13, cms2.hyp_ll_index()[i_hyp]) ) {
    ret |= (CUT_BIT(CUT_LL_CALOISO));
  }
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && GoodSusyLeptonIsolation(11, cms2.hyp_lt_index()[i_hyp]) ) {
    ret |= CUT_BIT(CUT_LT_CALOISO);
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && GoodSusyLeptonIsolation(11, cms2.hyp_ll_index()[i_hyp]) ) {
    ret |= CUT_BIT(CUT_LL_CALOISO);
  } 
  //  
//   // muon b tag, with 20 GeV upper cut on the muon pt
//   if (passMuonBVeto(i_hyp, true))
//     ret |= (CUT_BIT(CUT_PASS_MUON_B_VETO));
//   else ret |= (CUT_BIT(CUT_MUON_TAGGED));
//   // muon b tag, with no upper cut on the muon pt
//   if (passMuonBVeto(i_hyp, false))
//     ret |= (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT));
//   else ret |= (CUT_BIT(CUT_MUON_TAGGED_WITHOUT_PTCUT));

  // Z veto
  if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2)
    ret |= (CUT_BIT(CUT_PASS_ZVETO));
  else if (not inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))
    ret |= (CUT_BIT(CUT_PASS_ZVETO));
  else ret |= (CUT_BIT(CUT_IN_Z_WINDOW));

  // Z veto using additional leptons in the event
  if (not additionalZveto())
    ret |= (CUT_BIT(CUT_PASS_ADDZVETO));

  // track z veto
  // use one of CUT_Z_TRACK_VETO_HYP, CUT_Z_TRACK_VETO_TRACKS or CUT_Z_TRACK_VETO_HYP_OR_TRACKS
  
  // int modus = passTrackZVeto(i_hyp);
  int modus = 999;

  // modus 1: combine trk-isolated tracks with hyp ll or lt, set CUT_Z_TRACK_VETO_HYP
  if ( ! ( modus == 1 || modus == 3 ) ) ret |= CUT_BIT(CUT_Z_TRACK_VETO_HYP);

  // modus 2: combine trk-isolated tracks with another track, set CUT_Z_TRACK_VETO_TRACKS
  if ( ! ( modus == 2 || modus == 3 ) ) ret |= CUT_BIT(CUT_Z_TRACK_VETO_TRACKS);

  // modus 3: combine trk-isolated tracks with another track or hyp ll or lt , set CUT_Z_TRACK_VETO_HYP_OR_TRACKS
  if ( ! ( modus == 1 || modus == 2 || modus ==3 ) ) ret |= CUT_BIT(CUT_Z_TRACK_VETO_HYP_OR_TRACKS);

  // *****************************************************************
  // special handling for the fake rate cuts for now, because they
  // only work for emu
  // *****************************************************************
  if (myType != DILEPTON_EMU)
    return ret;
  // in addition, for the muons, check that they pass tight+iso
  // (since the fake rate is electron only right now)
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
    if ((ret & (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_CALOISO))) != 
	(CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_CALOISO)))
      return ret;
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
    if ((ret & (CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_CALOISO))) != 
	(CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_CALOISO)))
      return ret;
  }
  // now set the fake flags for the electron
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
    if (isFakeable(cms2.hyp_lt_index()[i_hyp]))
      ret |= CUT_BIT(CUT_ELFAKE_FAKEABLE_OBJECT);
    if (isNumeratorElectron(cms2.hyp_lt_index()[i_hyp]))
      ret |= CUT_BIT(CUT_ELFAKE_NUMERATOR);
    else ret |= CUT_BIT(CUT_ELFAKE_NOT_NUMERATOR);
  } else {
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
      if (isFakeable(cms2.hyp_ll_index()[i_hyp]))
	ret |= CUT_BIT(CUT_ELFAKE_FAKEABLE_OBJECT);
      if (isNumeratorElectron(cms2.hyp_ll_index()[i_hyp]))
	ret |= CUT_BIT(CUT_ELFAKE_NUMERATOR);
      else ret |= CUT_BIT(CUT_ELFAKE_NOT_NUMERATOR);
    }
  }

  return ret;
}

cuts_t Looper::TrilepSelect (int i_hyp)
{
  //------------------------------------------------------------
  // In a trilepton analysis, you would make your cuts here
  //------------------------------------------------------------

  cuts_t ret = 0;
  return ret;
}

cuts_t Looper::QuadlepSelect (int i_hyp)
{
  //------------------------------------------------------------
  // In a quadlepton analysis, you would make your cuts here
  //------------------------------------------------------------

  cuts_t ret = 0;
  return ret;
}

void Looper::FillEventHistos ()
{
  //------------------------------------------------------------
  // In an event-based analysis, you would fill your histos here
  //------------------------------------------------------------
}

void Looper::FillDilepHistos (int i_hyp)
{
  //------------------------------------------------------------
  // Example dilepton histo filling; edit for your application
  //------------------------------------------------------------

  // everybody histogram needs to know what hypothesis he is 
  const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
  // and what the event weight is 
  const double weight = Weight(i_hyp);

  // these are the cuts that the candidate passes:
  cuts_t cuts_passed = DilepSelect(i_hyp);

  // this is how to test that the candidate passes the cuts (which
  // we specified in the constructor when we made the looper)
  // (*note: the parentheses are important*):
  if ((cuts_passed & cuts_) == cuts_) {
    //     cout << "run: " << cms2.evt_run() << " event: " << cms2.evt_event() << endl;
    // if the candidate passed, we count it
    cands_passing_[myType] += weight;
    cands_passing_w2_[myType] += weight * weight;
    cands_count_[myType]++;
    cands_passing_[DILEPTON_ALL] += weight;
    cands_passing_w2_[DILEPTON_ALL] += weight * weight;
    cands_count_[DILEPTON_ALL]++;
    
    // do some printout for SameSign events:
    // penner
    //  TCut corCharge("(els_charge == -1 && els_mc_id == 11) || (els_charge == 1 && els_mc_id == -11)");

//     bool isBeauty = false;
//     bool isCharm  = false;

    // look at electron pairs
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
      if(
         (
          ( cms2.els_charge()[cms2.hyp_lt_index()[i_hyp]] == -1 && cms2.els_mc_id()[cms2.hyp_lt_index()[i_hyp]] == 11) &&
          ( cms2.els_charge()[cms2.hyp_ll_index()[i_hyp]] == -1 && cms2.els_mc_id()[cms2.hyp_ll_index()[i_hyp]] == 11) 
          ) ||
         (
          ( cms2.els_charge()[cms2.hyp_lt_index()[i_hyp]] == 1 && cms2.els_mc_id()[cms2.hyp_lt_index()[i_hyp]] == -11) &&
          ( cms2.els_charge()[cms2.hyp_ll_index()[i_hyp]] == 1 && cms2.els_mc_id()[cms2.hyp_ll_index()[i_hyp]] == -11) 
          )
         ) {
        std::cout<<"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"<<std::endl;
        std::cout<<"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"<<std::endl;
        std::cout<<"We have a genuine SS ee pair IDs are "<<" DS/lum/r/e: "<< cms2.evt_dataset()<<" / "<< cms2.evt_lumiBlock() <<" / "<<cms2.evt_run() <<" / "<< cms2.evt_event()<<" scale1fb "<<cms2.evt_scale1fb() <<std::endl<<
          " lt: "<<cms2.els_mc_id()[cms2.hyp_lt_index()[i_hyp]]<<" Mother "<<pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->GetName()<<" e pt "<<cms2.hyp_lt_p4()[i_hyp].pt()<<" e phi "<<cms2.hyp_lt_p4()[i_hyp].phi()<<" e eta "<<cms2.hyp_lt_p4()[i_hyp].eta()<<std::endl<<
          " ll: "<<cms2.els_mc_id()[cms2.hyp_ll_index()[i_hyp]]<<" Mother "<<pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->GetName()<<" e pt "<<cms2.hyp_ll_p4()[i_hyp].pt()<<" e phi "<<cms2.hyp_ll_p4()[i_hyp].phi()<<" e eta "<<cms2.hyp_ll_p4()[i_hyp].eta()<<std::endl;
        dumpDocLines();
        cout<<" MId lt "<<cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]]<<" C? "<<idIsCharm(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])<<" B? "<<idIsBeauty(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])<<endl;
        cout<<" MId ll "<<cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]]<<" C? "<<idIsCharm(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])<<" B? "<<idIsBeauty(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])<<endl;
        //         cout<<"lt Beauty? "<<pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Beauty()<<" MId "<<cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]]<<endl;
        //          cout<<"lt Charm?  "<<pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Charm()<<endl;
        //          cout<<"ll Beauty? "<<pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Beauty()<<endl;
        //          cout<<"ll Charm?  "<<pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Charm()<<endl;
        std::cout<<"eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"<<std::endl;
        std::cout<<"eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"<<std::endl;

        if( pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Beauty() ) {
          cout<<"Guggma: lt is Beauty"<<endl;
        }
        if( pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Beauty() ) {
          cout<<"Guggma: ll is Beauty"<<endl;
        }
        if( pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Charm() ) {
          cout<<"Guggma: lt is Charm"<<endl;
        }
        if( pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Charm() ) {
          cout<<"Guggma: ll is Charm"<<endl;
        }

      }
      else if ( 
               (
                ( cms2.els_charge()[cms2.hyp_lt_index()[i_hyp]] == -1 && cms2.els_mc_id()[cms2.hyp_lt_index()[i_hyp]] ==  11) &&
                ( cms2.els_charge()[cms2.hyp_ll_index()[i_hyp]] ==  1 && cms2.els_mc_id()[cms2.hyp_ll_index()[i_hyp]] == -11) 
                ) ||
               (
                ( cms2.els_charge()[cms2.hyp_lt_index()[i_hyp]] ==  1 && cms2.els_mc_id()[cms2.hyp_lt_index()[i_hyp]] == -11) &&
                ( cms2.els_charge()[cms2.hyp_ll_index()[i_hyp]] == -1 && cms2.els_mc_id()[cms2.hyp_ll_index()[i_hyp]] ==  11) 
                )
               ) {
        //      std::cout<<"We have a genuine OS ee pair IDs are lt: "<<cms2.els_mc_id()[cms2.hyp_lt_index()[i_hyp]]<<" ll: "<<cms2.els_mc_id()[cms2.hyp_ll_index()[i_hyp]]<<std::endl;
      }
      else {
        //      std::cout<<"Not a genuine e pair IDs are lt: "<<cms2.els_mc_id()[cms2.hyp_lt_index()[i_hyp]]<<" ll: "<<cms2.els_mc_id()[cms2.hyp_ll_index()[i_hyp]]<<std::endl;
        std::cout<<"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"<<std::endl;
        std::cout<<"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"<<std::endl;
        std::cout<<"We have a false-charge SS ee pair IDs are "<<" DS/lum/r/e: "<< cms2.evt_dataset()<<" / "<< cms2.evt_lumiBlock() <<" / "<<cms2.evt_run() <<" / "<< cms2.evt_event()<<" scale1fb "<<cms2.evt_scale1fb() <<std::endl<<
          " lt: "<<cms2.els_mc_id()[cms2.hyp_lt_index()[i_hyp]]<<" Mother "<<pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->GetName()<<" e pt "<<cms2.hyp_lt_p4()[i_hyp].pt()<<" e phi "<<cms2.hyp_lt_p4()[i_hyp].phi()<<" e eta "<<cms2.hyp_lt_p4()[i_hyp].eta()<<std::endl<<
          " ll: "<<cms2.els_mc_id()[cms2.hyp_ll_index()[i_hyp]]<<" Mother "<<pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->GetName()<<" e pt "<<cms2.hyp_ll_p4()[i_hyp].pt()<<" e phi "<<cms2.hyp_ll_p4()[i_hyp].phi()<<" e eta "<<cms2.hyp_ll_p4()[i_hyp].eta()<<std::endl;
        dumpDocLines();
        //         cout<<" MId lt "<<cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]]<<endl;
        //         cout<<" MId ll "<<cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]]<<endl;
        cout<<" MId lt "<<cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]]<<" C? "<<idIsCharm(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])<<" B? "<<idIsBeauty(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])<<endl;
        cout<<" MId ll "<<cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]]<<" C? "<<idIsCharm(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])<<" B? "<<idIsBeauty(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])<<endl;
//          cout<<"lt Beauty? "<<pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Beauty()<<endl;
//          cout<<"lt Charm?  "<<pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Charm()<<endl;
//          cout<<"ll Beauty? "<<pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Beauty()<<endl;
//          cout<<"ll Charm?  "<<pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Charm()<<endl;
        std::cout<<"eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"<<std::endl;
        std::cout<<"eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"<<std::endl;
        if( pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Beauty() ) {
          cout<<"Guggma: lt is Beauty"<<endl;
        }
        if( pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Beauty() ) {
          cout<<"Guggma: ll is Beauty"<<endl;
        }
        if( pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Charm() ) {
          cout<<"Guggma: lt is Charm"<<endl;
        }
        if( pdg->GetParticle(cms2.els_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Charm() ) {
          cout<<"Guggma: ll is Charm"<<endl;
        }
      }
    }
    
    // look at muon pairs
    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
      if(
         (
          ( cms2.mus_charge()[cms2.hyp_lt_index()[i_hyp]] == -1 && cms2.mus_mc_id()[cms2.hyp_lt_index()[i_hyp]] == 13) &&
          ( cms2.mus_charge()[cms2.hyp_ll_index()[i_hyp]] == -1 && cms2.mus_mc_id()[cms2.hyp_ll_index()[i_hyp]] == 13) 
          ) ||
         (
          ( cms2.mus_charge()[cms2.hyp_lt_index()[i_hyp]] == 1 && cms2.mus_mc_id()[cms2.hyp_lt_index()[i_hyp]] == -13) &&
          ( cms2.mus_charge()[cms2.hyp_ll_index()[i_hyp]] == 1 && cms2.mus_mc_id()[cms2.hyp_ll_index()[i_hyp]] == -13) 
          )
         ) {
         std::cout<<"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"<<std::endl;
         std::cout<<"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"<<std::endl;
         std::cout<<"We have a genuine SS mu mu pair IDs are: "<<" DS/lum/r/e: "<< cms2.evt_dataset()<<" / "<< cms2.evt_lumiBlock() <<" / "<<cms2.evt_run() <<" / "<< cms2.evt_event()<<" scale1fb "<<cms2.evt_scale1fb() <<std::endl<<
           " lt: "<<cms2.mus_mc_id()[cms2.hyp_lt_index()[i_hyp]]<<" Mother "<<pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->GetName()<<" mu pt "<<cms2.hyp_lt_p4()[i_hyp].pt()<<" mu phi "<<cms2.hyp_lt_p4()[i_hyp].phi()<<" mu eta "<<cms2.hyp_lt_p4()[i_hyp].eta()<<" type "<<cms2.mus_type()[cms2.hyp_lt_index()[i_hyp]]<<std::endl<<
           " ll: "<<cms2.mus_mc_id()[cms2.hyp_ll_index()[i_hyp]]<<" Mother "<<pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->GetName()<<" mu pt "<<cms2.hyp_ll_p4()[i_hyp].pt()<<" mu phi "<<cms2.hyp_ll_p4()[i_hyp].phi()<<" mu eta "<<cms2.hyp_ll_p4()[i_hyp].eta()<<" type "<<cms2.mus_type()[cms2.hyp_ll_index()[i_hyp]]<<std::endl;
         dumpDocLines();
//         cout<<" MId lt "<<cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]]<<endl;
//        cout<<" MId ll "<<cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]]<<endl;
        cout<<" MId lt "<<cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]]<<" C? "<<idIsCharm(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])<<" B? "<<idIsBeauty(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])<<endl;
        cout<<" MId ll "<<cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]]<<" C? "<<idIsCharm(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])<<" B? "<<idIsBeauty(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])<<endl;
        //          cout<<"lt Beauty? "<<pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Beauty()<<endl;
        //          cout<<"lt Charm?  "<<pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Charm()<<endl;
        //          cout<<"ll Beauty? "<<pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Beauty()<<endl;
        //          cout<<"ll Charm?  "<<pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Charm()<<endl;
         std::cout<<"eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"<<std::endl;
         std::cout<<"eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"<<std::endl;
        if( pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Beauty() ) {
          cout<<"Guggma: lt is Beauty"<<endl;
        }
        if( pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Beauty() ) {
          cout<<"Guggma: ll is Beauty"<<endl;
        }
        if( pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Charm() ) {
          cout<<"Guggma: lt is Charm"<<endl;
        }
        if( pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Charm() ) {
          cout<<"Guggma: ll is Charm"<<endl;
        }
      }
      else if ( 
               (
                ( cms2.mus_charge()[cms2.hyp_lt_index()[i_hyp]] == -1 && cms2.mus_mc_id()[cms2.hyp_lt_index()[i_hyp]] ==  13) &&
                ( cms2.mus_charge()[cms2.hyp_ll_index()[i_hyp]] ==  1 && cms2.mus_mc_id()[cms2.hyp_ll_index()[i_hyp]] == -13) 
                ) ||
               (
                ( cms2.mus_charge()[cms2.hyp_lt_index()[i_hyp]] ==  1 && cms2.mus_mc_id()[cms2.hyp_lt_index()[i_hyp]] == -13) &&
                ( cms2.mus_charge()[cms2.hyp_ll_index()[i_hyp]] == -1 && cms2.mus_mc_id()[cms2.hyp_ll_index()[i_hyp]] ==  13) 
                )
               ) {
        //      std::cout<<"We have a genuine OS mu mu pair"<<std::endl;
      }
      else {
        //      std::cout<<"Not a genuine mu mu pair. IDs are lt: "<<cms2.els_mc_id()[cms2.hyp_lt_index()[i_hyp]]<<" ll: "<<cms2.els_mc_id()[cms2.hyp_ll_index()[i_hyp]]<<std::endl;
         std::cout<<"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"<<std::endl;
         std::cout<<"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"<<std::endl;
         std::cout<<"We have a false-charge SS mu mu pair IDs are: "<<" DS/lum/r/e: "<< cms2.evt_dataset()<<" / "<< cms2.evt_lumiBlock() <<" / "<<cms2.evt_run() <<" / "<< cms2.evt_event()<<" scale1fb "<<cms2.evt_scale1fb() <<std::endl<<
           " lt: "<<cms2.mus_mc_id()[cms2.hyp_lt_index()[i_hyp]]<<" Mother "<<pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->GetName()<<" mu pt "<<cms2.hyp_lt_p4()[i_hyp].pt()<<" mu phi "<<cms2.hyp_lt_p4()[i_hyp].phi()<<" mu eta "<<cms2.hyp_lt_p4()[i_hyp].eta()<<" type "<<cms2.mus_type()[cms2.hyp_lt_index()[i_hyp]]<<std::endl<<
           " ll: "<<cms2.mus_mc_id()[cms2.hyp_ll_index()[i_hyp]]<<" Mother "<<pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->GetName()<<" mu pt "<<cms2.hyp_ll_p4()[i_hyp].pt()<<" mu phi "<<cms2.hyp_ll_p4()[i_hyp].phi()<<" mu eta "<<cms2.hyp_ll_p4()[i_hyp].eta()<<" type "<<cms2.mus_type()[cms2.hyp_ll_index()[i_hyp]]<<std::endl;
         dumpDocLines();
//         cout<<" MId lt "<<cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]]<<endl;
//         cout<<" MId ll "<<cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]]<<endl;
        cout<<" MId lt "<<cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]]<<" C? "<<idIsCharm(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])<<" B? "<<idIsBeauty(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])<<endl;
        cout<<" MId ll "<<cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]]<<" C? "<<idIsCharm(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])<<" B? "<<idIsBeauty(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])<<endl;
        //          cout<<"lt Beauty? "<<pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Beauty()<<endl;
        //          cout<<"lt Charm?  "<<pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Charm()<<endl;
        //          cout<<"ll Beauty? "<<pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Beauty()<<endl;
        //          cout<<"ll Charm?  "<<pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Charm()<<endl;
        std::cout<<"eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"<<std::endl;
         std::cout<<"eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"<<std::endl;
        if( pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Beauty() != 0) {
          cout<<"Guggma: lt is Beauty"<<endl;
        }
        if( pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Beauty() != 0 ) {
          cout<<"Guggma: ll is Beauty"<<endl;
        }
        if( pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_lt_index()[i_hyp]])->Charm()  != 0)  {
          cout<<"Guggma: lt is Charm"<<endl;
        }
        if( pdg->GetParticle(cms2.mus_mc_motherid()[cms2.hyp_ll_index()[i_hyp]])->Charm()  != 0 ) {
          cout<<"Guggma: ll is Charm"<<endl;
        }
      }
    }
    
  }

  // jet plots
  hnJet_->Fill(cuts_passed, myType, cms2.hyp_njets()[i_hyp] + nTrkJets(i_hyp), weight);
  hnCaloJet_->Fill(cuts_passed, myType, cms2.hyp_njets()[i_hyp], weight);
  hnTrackJet_->Fill(cuts_passed, myType, nTrkJets(i_hyp), weight);

  // lepton plots
  hltPt_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
  hllPt_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
  hltEta_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
  hllEta_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);

  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
    helPt_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
    helEta_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
    hltCaloIso_->Fill(cuts_passed, myType, el_rel_iso(cms2.hyp_lt_index()[i_hyp],true), weight);
    helCaloIso_->Fill(cuts_passed, myType, el_rel_iso(cms2.hyp_lt_index()[i_hyp],true), weight);
  } else {
    hmuPt_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
    hmuEta_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
    hltCaloIso_->Fill(cuts_passed, myType, mu_rel_iso(cms2.hyp_lt_index()[i_hyp]), weight);
    hmuCaloIso_->Fill(cuts_passed, myType, mu_rel_iso(cms2.hyp_lt_index()[i_hyp]), weight);
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
    helPt_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
    helEta_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
    hllCaloIso_->Fill(cuts_passed, myType, el_rel_iso(cms2.hyp_ll_index()[i_hyp],true), weight);
    helCaloIso_->Fill(cuts_passed, myType, el_rel_iso(cms2.hyp_ll_index()[i_hyp],true), weight);
  } else {
    hmuPt_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
    hmuEta_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
    hllCaloIso_->Fill(cuts_passed, myType, mu_rel_iso(cms2.hyp_ll_index()[i_hyp]), weight);
    hmuCaloIso_->Fill(cuts_passed, myType, mu_rel_iso(cms2.hyp_ll_index()[i_hyp]), weight);
  }

  // track plots
  hnTrack_->Fill(cuts_passed, myType, cms2.trks_trk_p4().size(), weight);

  // dilepton mass
  hdilMass_->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].mass(), weight);
    
  // delta phi btw leptons
  double dphi = TMath::Abs(cms2.hyp_lt_p4()[i_hyp].phi() - cms2.hyp_ll_p4()[i_hyp].phi());
  if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
  hdphiLep_->Fill(cuts_passed, myType, dphi, weight);

  // dilepton pt
  hdilPt_->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].pt(), weight);
    
  // MET plots
  hmet_->Fill(cuts_passed, myType, cms2.evt_tcmet(), weight);      
  hmetSpec_->Fill(cuts_passed, myType, 
		  MetSpecial(cms2.evt_tcmet(), cms2.evt_tcmetPhi(), i_hyp),
		  weight);
  // track correction to the met
  //  const TVector3 trkCorr = correctMETforTracks();
  //  TVector3 hyp_met;
  //  hyp_met.SetPtEtaPhi(cms2.evt_tcmet(), 0, cms2.hyp_metPhi()[i_hyp]);
  //  hyp_met += trkCorr;
  hmetTrkCorr_->Fill(cuts_passed, myType, cms2.evt_tcmet(), weight);

  // sumet plots
  hsumet_->Fill(cuts_passed, myType, sumEt_, weight);      

  // thrust plots
  hthrust_->Fill(cuts_passed, myType, thrust, weight);      
  hjptthrust_->Fill(cuts_passed, myType, jptthrust, weight);      

  // effective mass
  float meff = 0.;
  for ( unsigned int jet = 0;
	jet < cms2.hyp_jets_p4()[i_hyp].size();
	++jet ) {
    meff += cms2.hyp_jets_p4()[i_hyp][jet].pt();
  }
  meff += cms2.hyp_lt_p4()[i_hyp].pt();
  meff += cms2.hyp_ll_p4()[i_hyp].pt();

  hmeff_->Fill(cuts_passed, myType, meff+cms2.evt_tcmet(),weight);
  hmeffcorr_->Fill(cuts_passed, myType, meff+cms2.evt_tcmet(),weight);

}

void Looper::FillTrilepHistos (int i_hyp)
{
  //------------------------------------------------------------
  // In a trilepton analysis, you would fill your histos here
  //------------------------------------------------------------
}

void Looper::FillQuadlepHistos (int i_hyp)
{
  //------------------------------------------------------------
  // In a quadlepton analysis, you would fill your histos here
  //------------------------------------------------------------
}

void Looper::End ()
{
  //------------------------------------------------------------
  //Example status message at the end of a looper; edit for your
  //application
  //------------------------------------------------------------

  int ret = fprintf(logfile_, 
		    "Sample %10s: Total candidate count (ee em mm all): %8u %8u %8u %8u."
		    " Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n", 
		    sample_.name.c_str(),
		    CandsCount(DILEPTON_EE), CandsCount(DILEPTON_EMU), CandsCount(DILEPTON_MUMU), CandsCount(DILEPTON_ALL), 
		    CandsPassing(DILEPTON_EE) , RMS(DILEPTON_EE), 
		    CandsPassing(DILEPTON_EMU) , RMS(DILEPTON_EMU), 
		    CandsPassing(DILEPTON_MUMU), RMS(DILEPTON_MUMU), 
		    CandsPassing(DILEPTON_ALL) , RMS(DILEPTON_ALL));
  if (ret < 0)
    perror("writing to log file");
}
//__________________________________________________________________________________

Double_t Looper::calc_thrust(Double_t * thrust_axis_phi, Double_t * thrust_axis_theta, TString EFO ) { // start calc_thrust

  Double_t thrust;
  Double_t p_dot_thr_sum, p_sum;
  Double_t bufferphi, buffertheta;

  bufferphi   = *thrust_axis_phi;
  buffertheta = *thrust_axis_theta;

  TVector3 zufo_pi(1.,1.,1.);
  TVector3 ThAxLoc(1.,1.,1.);

  ThAxLoc.SetXYZ( TMath::Sin(buffertheta)*TMath::Cos(bufferphi) , TMath::Sin(bufferphi)*TMath::Sin(buffertheta) ,  TMath::Cos(buffertheta) );
  //  ThAxLoc.SetTheta(buffertheta);
  //  ThAxLoc.SetMag(1.);
  thrust        = -999;
  p_dot_thr_sum = 0.;
  p_sum         = 0.;

  // *** calc thrust form zufos here (at the moment with all zufos. Intend to prepare a list of 'chosen' zufos for 3d version):
  // start with tracks, try JPT also, keep kicking Kostas butt for doing efos ;)
  //  for(Int_t zufonr =0; zufonr<nzufos; zufonr++) {
  if( EFO.Contains("track") ) {
    for( unsigned int trkCount = 0; trkCount < cms2.trks_trk_p4().size(); ++trkCount ) {
      
      
      zufo_pi.SetXYZ(cms2.trks_trk_p4()[trkCount].x(), cms2.trks_trk_p4()[trkCount].y(), cms2.trks_trk_p4()[trkCount].z()); // activcated this line 040311 to allow theta iteration.
      //    zufo_pi.SetXYZ(zufo[zufonr][0], zufo[zufonr][1], 0.0); // removed this line 040311 to allow theta iteration.
      
      // *** sum up all p_i * n (with theta>10. deg):
      if( zufo_pi.Theta() > 10.*(TMath::Pi()/180.) || 42 == 42 ){ // throw out all zufoses with theta<10. deg - deactivated10 deg cut  090625 (using all EFOs)
        p_dot_thr_sum += TMath::Abs(zufo_pi.Dot( ThAxLoc ));
        p_sum += TMath::Sqrt( zufo_pi.Dot(zufo_pi)  );
      }
    }
  }
  else if( EFO.Contains("JPT") ) {
    //    for( unsigned int trkCount = 0; trkCount < cms2.trks_trk_p4().size(); ++trkCount ) {
      for ( unsigned int ijptjet = 0; ijptjet < cms2.jpts_p4().size(); ++ijptjet) {
      
      zufo_pi.SetXYZ(cms2.jpts_p4()[ijptjet].x(), cms2.jpts_p4()[ijptjet].y(), cms2.jpts_p4()[ijptjet].z()); // activcated this line 040311 to allow theta iteration.
      //    zufo_pi.SetXYZ(zufo[zufonr][0], zufo[zufonr][1], 0.0); // removed this line 040311 to allow theta iteration.
      
      // *** sum up all p_i * n (with theta>10. deg):
      if( zufo_pi.Theta() > 10.*(TMath::Pi()/180.) || 42 == 42 ){ // throw out all zufoses with theta<10. deg - deactivated10 deg cut  090625 (using all EFOs)
        p_dot_thr_sum += TMath::Abs(zufo_pi.Dot( ThAxLoc ));
        p_sum += TMath::Sqrt( zufo_pi.Dot(zufo_pi)  );
      }
    }
  }
  else { // default to track based calculation
    cout<<"Thrust defaults - is calcualted based on tracks."<<endl;
    for( unsigned int trkCount = 0; trkCount < cms2.trks_trk_p4().size(); ++trkCount ) {
      
      
      zufo_pi.SetXYZ(cms2.trks_trk_p4()[trkCount].x(), cms2.trks_trk_p4()[trkCount].y(), cms2.trks_trk_p4()[trkCount].z()); // activcated this line 040311 to allow theta iteration.
      //    zufo_pi.SetXYZ(zufo[zufonr][0], zufo[zufonr][1], 0.0); // removed this line 040311 to allow theta iteration.
      
      // *** sum up all p_i * n (with theta>10. deg):
      if( zufo_pi.Theta() > 10.*(TMath::Pi()/180.) || 42 == 42 ){ // throw out all zufoses with theta<10. deg - deactivated10 deg cut  090625 (using all EFOs)
        p_dot_thr_sum += TMath::Abs(zufo_pi.Dot( ThAxLoc ));
        p_sum += TMath::Sqrt( zufo_pi.Dot(zufo_pi)  );
      }
    }
  }
    //  cout<<endl<<" after zufo loop: "<<" p_sum: "<<p_sum<<" p_dot_thr_sum: "<<p_dot_thr_sum<<endl;
  thrust = p_dot_thr_sum / p_sum;
  return thrust;
}
//__________________________________________________________________________________

Double_t Looper::find_thrust_axis(Double_t * thrust_axis_phi, Double_t * thrust_axis_theta, Double_t * thrust, TString EFO) { // start find_thrust_axis

  Int_t phisteps       = 20;
  //  Int_t phisteps          = 60; // changed from 20 to 60 040225
  //  Int_t phisteps       = 3240; // changed from 60 to 3240 to test theta/phi scan option. -> takes too long!!
  //  Int_t thetasteps        = 60;
  Int_t thetasteps        = 20;
  Double_t phi            = -999.;
  Double_t phi_max        = -999.;
  Double_t theta          = TMath::Pi()/2.;
  Double_t theta_max      = -999.;
  Double_t pi             = TMath::Pi();
  Double_t thrust_iter    = -999.;
  Double_t thrust_max     = -999.;
  Double_t thrust_max_phi = -999.;

  // calc the thrust for phi, not changing theta:  
  for(Int_t phistep = 0; phistep < phisteps; phistep++) {

    phi = (Double_t)phistep*(1.*pi/((Double_t)phisteps-1.));

    *thrust_axis_phi   = phi;
    *thrust_axis_theta = theta;

    // *** calc the thrust for given phi
    thrust_iter = calc_thrust( thrust_axis_phi, thrust_axis_theta, EFO );

    // *** compare calculated thrust with previous values, store maximal thrust - phi based here:
    if(thrust_iter > thrust_max ) {
      thrust_max = thrust_iter;
      thrust_max_phi = thrust_max;
      phi_max    = phi;            // or elegantly use: thrust_axis_local.Phi();
      theta_max  = theta;          // or elegantly use: thrust_axis_local.Theta();
    }
  }
  //  cout<<"find_thrust_axis found - phi: "<<phi_max*180./TMath::Pi()<<" theta: "<<theta_max*180./TMath::Pi()<<" thrust both / phi only:" <<thrust_max<<" / "<<thrust_max_phi<<endl;

  // refine the thrust for theta, keeping the found optimal phi:  
  for(Int_t thetastep = 0; thetastep < thetasteps; thetastep++) {

    theta = (Double_t)thetastep*(1.*pi/((Double_t)thetasteps-1.));

    *thrust_axis_phi   = phi_max;
    *thrust_axis_theta = theta;

    // *** calc the thrust for given theta
    thrust_iter = calc_thrust( thrust_axis_phi, thrust_axis_theta );
    
    //        cout<<"thrust new is "<<thrust_iter<<" thrust_max: "<<thrust_max<<" at theta"<<theta<<" in theta loop"<<endl;
        //    getchar();


    // *** compare calculated thrust with previous values, store maximal thrust - theta based here:
    if(thrust_iter > thrust_max ) {
      thrust_max = thrust_iter;
      theta_max  = theta;          // or elegantly use: thrust_axis_local.Theta();
      //      cout<<"theta_max is "<<theta_max<<" in theta_max step loop"<<endl;
    }
  }

  // *** set output values:
  *thrust            = thrust_max;
  *thrust_axis_phi   = phi_max;
  *thrust_axis_theta = theta_max;

  //  thrust_axis_local.SetTheta(theta_max);  // keeping rho and phi
  //  thrust_axis_local.SetPhi(phi_max);      // keeping rho and theta
  //  thrust_axis_local.SetMag(1.);           // keeping theta and phi
  //  cout<<"find_thrust_axis found - phi: "<<phi_max*180./TMath::Pi()<<" theta: "<<theta_max*180./TMath::Pi()<<" thrust both / phi only:" <<thrust_max<<" / "<<thrust_max_phi<<endl;
  //  getchar();

  if(42)  return *thrust; // dummy - prep for error messages
  else  return -1.;

}

bool Looper::idIsCharm(int id) {
  id = abs(id);
  if (
      id == 411     ||
      id == 421     ||
      id == 10411   ||
      id == 10421   ||
      id == 413     ||
      id == 423     ||
      id == 10413   ||
      id == 10423   ||
      id == 20413   ||
      id == 20423   ||
      id == 415     ||
      id == 425     ||
      id == 431     ||
      id == 10431   ||
      id == 433     ||
      id == 10433   ||
      id == 20433   ||
      id == 435     ||
      id == 441     ||
      id == 10441   ||
      id == 100441  ||
      id == 443     ||
      id == 10443   ||
      id == 20443   ||
      id == 100443  ||
      id == 30443   ||
      id == 9000443 ||
      id == 9010443 ||
      id == 9020443 ||
      id == 445     ||
      id == 9000445 ||
      id == 4122    ||
      id == 4222    ||
      id == 4212    ||
      id == 4112    ||
      id == 4224    ||
      id == 4214    ||
      id == 4114    ||
      id == 4232    ||
      id == 4132    ||
      id == 4322    ||
      id == 4312    ||
      id == 4324    ||
      id == 4314    ||
      id == 4332    ||
      id == 4334    ||
      id == 4412    ||
      id == 4422    ||
      id == 4414    ||
      id == 4424    ||
      id == 4432    ||
      id == 4434    ||
      id == 4444
      ) {
    return true;
  }
  else return false;
}

bool Looper::idIsBeauty(int id) {
  id = abs(id);
  if (
      id == 5       ||
      id == 511     ||
      id == 521     ||
      id == 10511   ||
      id == 10521   ||
      id == 513     ||
      id == 523     ||
      id == 10513   ||
      id == 10523   ||
      id == 20513   ||
      id == 20523   ||
      id == 515     ||
      id == 525     ||
      id == 531     ||
      id == 10531   ||
      id == 533     ||
      id == 10533   ||
      id == 20533   ||
      id == 535     ||
      id == 541     ||
      id == 10541   ||
      id == 543     ||
      id == 10543   ||
      id == 20543   ||
      id == 545     ||
      id == 551     ||
      id == 10551   ||
      id == 100551  ||
      id == 110551  ||
      id == 200551  ||
      id == 210551  ||
      id == 553     ||
      id == 10553   ||
      id == 20553   ||
      id == 30553   ||
      id == 100553  ||
      id == 110553  ||
      id == 120553  ||
      id == 130553  ||
      id == 200553  ||
      id == 210553  ||
      id == 220553  ||
      id == 300553  ||
      id == 9000553 ||
      id == 9010553 ||
      id == 555     ||
      id == 10555   ||
      id == 20555   ||
      id == 100555  ||
      id == 110555  ||
      id == 120555  ||
      id == 200555  ||
      id == 557     ||
      id == 100557  ||
      id == 5122    || 
      id == 5112    ||
      id == 5212    ||
      id == 5222    ||
      id == 5114    ||
      id == 5214    ||
      id == 5224    ||
      id == 5132    ||
      id == 5232    ||
      id == 5312    ||
      id == 5322    ||
      id == 5314    ||
      id == 5324    ||
      id == 5332    ||
      id == 5334    ||
      id == 5142    ||
      id == 5242    ||
      id == 5412    ||
      id == 5422    ||
      id == 5414    ||
      id == 5424    ||
      id == 5342    ||
      id == 5432    ||
      id == 5434    ||
      id == 5442    ||
      id == 5444    ||
      id == 5512    ||
      id == 5522    ||
      id == 5514    ||
      id == 5524    ||
      id == 5532    ||
      id == 5534    ||
      id == 5542    ||
      id == 5544    ||
      id == 5554 
      ) {
    return true;
  }
  else return false;
}
