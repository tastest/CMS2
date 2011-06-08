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
    NominalWWV1 = 16,
    muonSelectionFO_mu_wwV1 = 17,
    muonSelectionFO_mu_wwV1_iso10 = 18,
    OSGeneric_v1 = 19,
    OSZ_v1 = 20,
    NominalSSv2 = 21,
    muonSelectionFO_mu_ssV2 = 22,
    muonSelectionFO_mu_ssV2_iso10 = 23,
    muonSelectionFO_mu_wwV1_iso10_d0 = 24,
    NominalTTbar_pass6 = 25,
    muonSelectionFO_mu_ttbar_pass6 = 26,
    OSGeneric_v2 = 27,
    NominalSmurfV3 = 28,
    muonSelectionFO_mu_smurf_04 = 29,
    muonSelectionFO_mu_smurf_10 = 30,
    NominalSSv3 = 31,
    muonSelectionFO_ssV3 = 32,
    OSZ_v2 = 33,
    NominalSmurfV4 = 34,
    OSGeneric_v2_FO = 35,
    NominalSmurfV5 = 36,
    OSGeneric_v3 = 37,
    OSGeneric_v3_FO = 38,
    NominalSSv4 = 39,
    muonSelectionFO_ssV4 = 40,
    NominalSmurfV6 = 41
}; 

//
// desined for high Pt muons (above 10 GeV). Nominal selection
// is tight enough for most if not all analysis. Tight 
// version is a backup.
//
bool muonId(unsigned int index, SelectionType type = Nominal, int vertex_index = -1);
bool muonIdNotIsolated(unsigned int index, SelectionType type = Nominal, int vertex_index = -1);

//
// relative isolation
//
double muonIsoValueOriginal( unsigned int index );
double muonIsoValue     ( unsigned int, bool = true );
double muonIsoValue_TRK ( unsigned int, bool = true );
double muonIsoValue_ECAL( unsigned int, bool = true );
double muonIsoValue_HCAL( unsigned int, bool = true );
double muonIsoValue_FastJet  ( unsigned int , bool = true );
double mu_fastjet_rel_offset ( unsigned int , bool = true );

double muonIsoValuePF   ( unsigned int imu, unsigned int idavtx, float coner=0.4, float minptn=1.0, float dzcut=0.1);

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
double mud0PV_wwV1(unsigned int index);
double mudzPV_wwV1(unsigned int index);
double mud0PV_smurfV3(unsigned int index);
double mudzPV_smurfV3(unsigned int index);

//
// checks if muon is also pfmuon, and pfmuon pt = reco muon pt
//
bool isPFMuon( int index , bool requireSamePt = true , float dpt_max = 1.0 );

#endif

