// Lepton selections
bool goodElectronWithoutIsolation(int index);
bool goodMuonWithoutIsolation(int index);
bool passMuonIsolation(int index);
bool goodMuonIsolated(int index);
bool passElectronIsolation(int index);
bool goodElectronIsolated(int index);
//loose lepton (20GeV, Itk > 0.5, Ical > 0.5; |eta|<2.4
bool looseElectron20WithIsolation05(int index);
bool looseMuon20WithIsolation05(int index);

// Event level stuff
bool inZmassWindow(float mass);
bool passMet();
bool additionalZveto();
// Auxiliary stuff
int getDrellYanType();
bool isDYee();
bool isDYmm();
bool isDYtt();
void dumpDocLines();
void leptonCount(int& nele, int& nmuon, int& ntau);

bool  haveExtraMuon();
bool  haveExtraMuon5();
