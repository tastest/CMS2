// -*- C++ -*-

// $Id: fakerates.h,v 1.1 2010/02/06 00:00:37 spadhi Exp $

// for electron fake rates
bool isFakeable (int i_el);
double elFakeProb (int i_el, int add_error_times = 0);
bool isNumeratorElectron (int index, int type = 0);

int elFakeMCCategory(int i_el);

// for muon fake rates
bool isFakeableMuon (int i_mu);
double muFakeProb (int i_mu, int add_error_times = 0);
bool isNumeratorMuon (int index, int type = 0);
double FakeProb_v1 (int i_el, int add_error_times, int id);

int muFakeMCCategory(int i_mu);


class TH2F &fakeRate ();
class TH2F &fakeRateError ();
class TH2F &fakeRateMuon ();
class TH2F &fakeRateErrorMuon ();
