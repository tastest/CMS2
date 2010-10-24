#ifndef MUON_SELECTIONS_H
#define MUON_SELECTIONS_H
// #include "CMS2.h"
//
// muon id
// see http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/MuonID2010
//

enum SelectionType { 
	 FakableObject = 0, 
	 Nominal = 1, 
	 Tight = 2,
	 NominalTTbar = 3,
	 muonSelectionFO_mu_v1 = 4,
	 muonSelectionFO_mu_ttbar = 5,
	 NominalTTbarV2 = 6,
	 muonSelectionFO_mu_ttbar_iso10 = 7,
	 NominalWWV0 = 8,
	 muonSelectionFO_mu_ww = 9,
	 muonSelectionFO_mu_ww_iso10 = 10,
	 NominalSS = 11,
	 NominalSSd0PV = 12,
	 NominalSSnod0 = 13,
	 muonSelectionFO_mu_ss = 14,
	 muonSelectionFO_mu_ss_iso10 = 15,
}; 

//
// desined for high Pt muons (above 10 GeV). Nominal selection
// is tight enough for most if not all analysis. Tight 
// version is a backup.
//
bool muonId(unsigned int index, SelectionType type = Nominal);
bool muonIdNotIsolated(unsigned int index, SelectionType type = Nominal);

//
// relative isolation
//
double muonIsoValue(unsigned int index);

//
// Trigger
// for now we don't match muons to triggered objects
//
bool passedMuonTriggerRequirements();

//
// cosmics rejection
//
bool isCosmics(unsigned int index);

//
// Muon d0 corrected by PV
//

double mud0PV(unsigned int index);

#endif

