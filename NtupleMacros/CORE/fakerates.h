// -*- C++ -*-
#ifndef FAKERATES_H
#define FAKERATES_H

// list of available FR versions:
enum fakeRateVersion {
  el_v1_cand01,
  el_v1_cand02,
  el_v2_cand01,
  el_v2_cand02,
  el_v3_cand01,
  el_v3_cand02,
  mu_v1,
};

bool isFakeableMuon (int i_mu, enum fakeRateVersion version);
double muFakeProb (int i_mu, enum fakeRateVersion version);
double muFakeProbErr (int i_mu, enum fakeRateVersion version);

bool isFakeableElectron (int i_el, enum fakeRateVersion version);
double elFakeProb (int i_el, enum fakeRateVersion version);
double elFakeProbErr (int i_el, enum fakeRateVersion version);

int elFakeMCCategory(int i_el);
int muFakeMCCategory(int i_mu);

#endif
