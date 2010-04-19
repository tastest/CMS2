// -*- C++ -*-
#ifndef FAKERATES_H
#define FAKERATES_H

// list of available FR versions:
enum fakeRateVersion {
  el_ttbar_v1,   
  el_ttbar_v2,   
  el_ttbar_v3,   

  el_v1_cand01, // pure cand01 selection
  el_v1_cand02, // cand02 + extra
  el_v1_cand02flip, // cand02 + extra + chargeflipVeto
  el_v2_cand01,
  el_v2_cand02,
  el_v2_cand02flip,
  el_v3_cand01,
  el_v3_cand02,
  el_v3_cand02flip,
  mu_v1,
  mu_ttbar
};

bool isFakeableMuon (int i_mu, enum fakeRateVersion version);
double muFakeProb (int i_mu, enum fakeRateVersion version);
double muFakeProbErr (int i_mu, enum fakeRateVersion version);

bool isFakeableElectron (int i_el, enum fakeRateVersion version);
double elFakeProb (int i_el, enum fakeRateVersion version);
double elFakeProbErr (int i_el, enum fakeRateVersion version);

#endif
