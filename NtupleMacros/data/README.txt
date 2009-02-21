This directory contains different versions of electron fake rate histograms.

v2_2:
------
Using Denominator and Numerator definitions below. 
Only QCD ptHat bins 80-300 were used for these fakeRates.

bool isDenominatorElectron(int index) {
  //
  // returns true if input fulfills certain cuts
  //
  // cut definition
  float et_cut        = 0.;
  float pt_cut        = 15.;
  float eta_cut       = 2.5;
  float iso_ratio_cut = 0.92;
  float eOverP_cut    = 999999.99;
  float hOverE_cut    = 0.2;

  float iso_ratio = 0.0;

  if( (els_p4->at(index).Pt()+els_tkIso->at(index)) > 0.0 ) iso_ratio = els_p4->at(index).Pt()/(els_p4->at(index).Pt()+els_tkIso->at(index));
  else iso_ratio = 0.0; // reject events with 0 momentum - do we have thses at all?

  bool result = true;

  if ( els_ESc->at(index)      < et_cut )                 result = false;
  if ( els_p4->at(index).Pt()  < pt_cut )                 result = false;
  if ( std::abs(els_p4->at(index).Eta()) > eta_cut )      result = false;
  if ( iso_ratio               < iso_ratio_cut )          result = false;
  if ( els_eOverPIn->at(index) > eOverP_cut )             result = false;
  if ( els_hOverE->at(index)   > hOverE_cut )             result = false;

  return result;

}

bool isNumeratorElectron(int index, int type=0) { // 0=loose, 1=tight, for pass4: 1=loose, 2=tight
  //
  // returns true if input fulfills certain cuts
  //
  // cut definition
  float et_cut        = 0.;
  float pt_cut        = 15;
  float eta_cut       = 2.5;
  float iso_ratio_cut = 0.92; // use alternatively to iso cut above
  float eOverP_cut    = 999999.99;
  float hOverE_cut    = 0.2;
  float d0_cut        = 0.025;

  float iso_ratio = els_p4->at(index).Pt()/(els_p4->at(index).Pt()+els_tkIso->at(index));
  
  bool result = true;
  
  if ( els_ESc->at(index)      < et_cut )                 result = false;
  if ( els_p4->at(index).Pt()  < pt_cut )                 result = false;
  if ( std::abs(els_p4->at(index).Eta()) > eta_cut )      result = false;
  if ( iso_ratio               < iso_ratio_cut )          result = false;
  if ( els_eOverPIn->at(index) > eOverP_cut )             result = false;
  if ( els_hOverE->at(index)   > hOverE_cut )             result = false;
  if ( std::abs(els_d0->at(index))  > d0_cut )            result = false;
  
  // _pass4 has 3 types - 0=robust, 1=loose, 2=tight
  // - need to adjust in all places here
  bool IdCuts = electron_selection_pass4(index, type); // eleID from code from Avi
  if (!IdCuts) result = false;
  
  return result;
  
}


v2_2_allpt:
------
same as v2_2
using ALL xsec weighted QCD pt-hat samples to produce this guy. This is 
daring, as the low pthat bins are severely statistics limited (the weights are HUGE).


v2_3_allpt_ljet_gt_30: (not in cvs)
------
same as v2_2_allpt but additionally requiring a leading jet with 
jetpt>30 GeV (lowest HLTrigger cut) for both the numerator and denominator


v4_0:
------
same as v2_3 using only ptHat bins 0 to 600 BUT:
Using the ele ID from ntuples
Using uncorrected Jets for the leading Jet Veto
Using ptHat bins 0 to 600
Using code from cvs :)
see http://omega.physics.ucsb.edu/twiki/bin/viewfile/CMS/20081014AgendaMinutes?rev=2;filename=EFakes_ogu_ibl_081014.pdf

v5_0: VERY preliminary
------
same as v2_3 using only ptHat bins 0 to 600 BUT:
Using the ele ID from ntuples
Using uncorrected Jets for the leading Jet Veto
Using ptHat bins 0 to 600
Using CaloIso+trackIso!
result: overestimation

v5_1:
------
same as v5_0 BUT:
Using CaloIso+trackIso on Numerator
Using trackIso         on Denominator
result: underestimation

v5_2:
------
same as v5_0 BUT:
Using CaloIso+trackIso                on Numerator
Using CaloIso(cut at 0.8)+trackIso    on Denominator
result: overestimation (pred 91, obs 81)

v5_3:
------
same as v5_0 BUT:
Using CaloIso+trackIso                on Numerator
Using CaloIso(cut at 0.85)+trackIso   on Denominator
result: overestimation (pred 107, obs 81)

v5_4:
------
same as v5_0 BUT:
Using CaloIso+trackIso                on Numerator
Using CaloIso(cut at 0.75)+trackIso   on Denominator
result: good estimation (pred 85, obs 81)

v5_5:
-----
same as v5_4 BUT:
use muon veto near electron           on Denominator
use absolute value of eta in numerator and denominator to increase statistics
result: better estimation (pred 79, obs 81)

v6_1:
-----
same as v5_5 BUT:
use d0corr
use V01-02-01 new ntuples with bug in d0corr
use EMenriched samples

v6_2:
-----
same selection as v6_1 BUT:
use new ntuples V01-02-06 with fixed d0corr

v6_3:
-----
same as v6_2 BUT:
add BCtoE samples to EMenriched samples

v6_4:
-----
same as v6_2 BUT:
use QCDpt30 sample

v6_5:
-----
same as v6_4 BUT:
remove trigger simulation (don't require at least one jet with uncorrected pt > 30 GeV)
only to be used with fake rates without trigger bias removal (take all
electrons, do not discard electron near trigger jet)

v6_6:
-----
same as v6_3 BUT:
new final WW selection cuts (cms2.els_tightId22XMinMatteo and tcMET) from 2/12/09

v6_7:
-----
same as v6_5 BUT:
new final WW selection cuts (cms2.els_tightId22XMinMatteo and tcMET) from 2/12/09

v7_0:
-----
Twiki: http://omega.physics.ucsb.edu/twiki/bin/view/CMS/FakeRateIncQCDpt30FinalWWcutsStatisticsOptimized

        Determination:
        - selections as for v6_7
        - uses incQCDpt30 sample
        - does NOT have leading jet > 30 GeV requirement (i.e. Trigger simulation disabled)
        - does HAVE fewer bins (3 in eta (-2.5, -1.1, 1.1, 2.5), 3 in pT (0, 20, 60, 150))
        Application:
        - changed fakerates.cc standard for v7 so the used histogram
          does NOT have FO close to leading jets removed (i.e. Trigger bias removal disabled)

v7_1:
-----
Twiki: http://omega.physics.ucsb.edu/twiki/bin/view/CMS/FakeRateIncQCDpt30FinalWWcutsStatisticsOptimizedII

        Determination:
        - selections as for v6_7
        - uses incQCDpt30 sample
        - does NOT have leading jet > 30 GeV requirement (i.e. Trigger simulation disabled)
ONLY CHANGE to 7_0:
        - does HAVE fewer bins (3 in eta (-2.5, -1.479, 1.479, 2.5), 3 in pT (0, 20, 60, 150))
        Application:
        - not yet changed default to 7_1 -> decide in next TaS meeting.
