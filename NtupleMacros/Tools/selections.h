// Lepton selections
bool goodElectronWithoutIsolation(int index);
bool goodMuonWithoutIsolation(int index);
bool passMuonIsolation(int index);
bool goodMuonIsolated(int index);
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



