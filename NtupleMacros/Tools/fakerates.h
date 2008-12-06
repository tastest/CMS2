// -*- C++ -*-

bool isFakeable (int i_el);
double elFakeProb (int i_el, int add_error_times = 0);
bool isNumeratorElectron (int index, int type = 0);

class TH2F &fakeRate ();
class TH2F &fakeRateError ();
