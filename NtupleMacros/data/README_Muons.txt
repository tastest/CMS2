v1_0:
Twiki: http://omega.physics.ucsb.edu/twiki/bin/view/CMS/FakeRateAllQCDpt30FinalWWcutsV1_0
Talk: http://omega.physics.ucsb.edu/twiki/pub/CMS/20090224AgendaMinutes/090224_electron_fakes.pdf
initial attempt at Muon fake rate.
*****************
Numerator:
*****************
bool isFakeNumeratorMuon_v1 (int index, int type)
{
  //
  // 1=loose, 2=tight
  //
  // returns true if input fulfills certain cuts
  //

  // cut definition
  float pt_cut                  = 20;
  float eta_cut                 = 2.5;

  bool result = true;

  // need:
  // - global muon
  // - pt > 20 GeV
  // - d0corr < 0.1 cm
  // - chi2/ndf < 10
  // - isoSumPt < 3 -> needs revision - would go for same def as for electrons?

  if ( cms2.mus_p4()[index].Pt()  < pt_cut )              result = false;
  if ( TMath::Abs(cms2.mus_p4()[index].Eta()) > eta_cut ) result = false;
  if ( !passMuonIsolation(index) )                        result = false;
  if ( type == 1 ) {
    // loose
    // Currently this is the same as tight
    if ( !goodMuonWithoutIsolation(index) )               result = false;
  } else if ( type == 2 ) {
    // tight
    if ( !goodMuonWithoutIsolation(index) )               result = false;
  } else {
    cout << "WARNING: wrong muon type detected, please select loose (1) or tight (2)" << endl;
  }

  return result;

}
*****************
Denominator:
*****************
bool isFakeDenominatorMuon_v1 (int index)
{
  //
  // returns true if input fulfills certain cuts
  //

  // cut definition
  float pt_cut                  = 20.;
  float eta_cut                 = 2.5;

  bool result = true;

  // need:

  // - global muon
  if ( (cms2.mus_type().at(index)&0x2)==0 )                                result = false;

  // - pt > 20 GeV
  if ( cms2.mus_p4()[index].Pt()  < pt_cut )                               result = false;

  // - d0corr < 0.1 cm -> loosened to 0.2!!
  if (TMath::Abs(cms2.mus_d0corr().at(index))   > 0.2)                     result = false;

  // - chi2/ndf < 20 (?)
  if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) > 20.) result = false;

  // - isoSumPt < 15 -> needs revision - would go for same def as for electrons?
  // - changed iso to mu_rel_iso > 0.75
    if ( !passMuonIsolationLoose(index) )                                  result = false;

  //  if (cms2.mus_validHits().at(index) < 11)                             result = false;

  if ( TMath::Abs(cms2.mus_p4()[index].Eta()) > eta_cut )                  result = false;

  return result;

}
