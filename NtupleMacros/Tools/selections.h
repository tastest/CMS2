#ifndef SELECTION_H
#define SELECTION_H

bool inZmassWindow (float mass);
bool goodElectronWithoutIsolation(int index);
bool goodMuonWithoutIsolation(int index);
bool passElectronIsolation(int index);
bool passMuonIsolation(int index);
bool goodMuonIsolated(int index);
bool goodElectronIsolated(int index);
bool pass2Met(int index);
double nearestDeltaPhi(double Phi, int index);
double MetSpecial(double MET, double MetPhi, int index);
bool pass4Met(int index);
void correctMETmuons_crossedE(double& met, double& metPhi, 
			      double muon_pt, double muon_phi,
			      double muon_track_theta, double muon_track_phi,
			      double mu_crossedem_dep, double mu_crossedhad_dep, double mu_crossedho_dep );
#endif
