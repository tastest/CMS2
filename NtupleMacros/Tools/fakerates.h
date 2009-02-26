// -*- C++ -*-

// for electron fake rates
bool isFakeable (int i_el);
double elFakeProb (int i_el, int add_error_times = 0);
bool isNumeratorElectron (int index, int type = 0);

// for muon fake rates
bool isFakeableMuon (int i_mu);
double muFakeProb (int i_mu, int add_error_times = 0);
bool isNumeratorMuon (int index, int type = 0);

class TH2F &fakeRate ();
class TH2F &fakeRateError ();
class TH2F &fakeRateMuon ();
class TH2F &fakeRateErrorMuon ();
