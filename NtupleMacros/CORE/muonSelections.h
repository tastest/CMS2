#ifndef MUON_SELECTIONS_H
#define MUON_SELECTIONS_H

///////////////
// Selectors //
///////////////

enum SelectionType { 
    NominalWWV0 = 8,
    muonSelectionFO_mu_ww = 9,
    muonSelectionFO_mu_ww_iso10 = 10,
    NominalWWV1 = 16,
    muonSelectionFO_mu_wwV1 = 17,
    muonSelectionFO_mu_wwV1_iso10 = 18,
    muonSelectionFO_mu_wwV1_iso10_d0 = 24,
    NominalSmurfV3 = 28,
    muonSelectionFO_mu_smurf_04 = 29,
    muonSelectionFO_mu_smurf_10 = 30,
    NominalSSv3 = 31,
    muonSelectionFO_ssV3 = 32,
    OSZ_v2 = 33,
    NominalSmurfV4 = 34,
    NominalSmurfV5 = 36,
    OSGeneric_v3 = 37,
    OSGeneric_v3_FO = 38,
    NominalSSv4 = 39,
    muonSelectionFO_ssV4 = 40,
    NominalSmurfV6 = 41
}; 

////////////////////
// Identification //
////////////////////

bool muonId           ( unsigned int index, SelectionType type, int vertex_index = -1);
bool muonIdNotIsolated( unsigned int index, SelectionType type, int vertex_index = -1);

///////////////
// Isolation //
///////////////
double muonIsoValue          ( unsigned int , bool = true );
double muonIsoValue_TRK      ( unsigned int , bool = true );
double muonIsoValue_ECAL     ( unsigned int , bool = true );
double muonIsoValue_HCAL     ( unsigned int , bool = true );
double muonIsoValue_FastJet  ( unsigned int , bool = true );
double mu_fastjet_rel_offset ( unsigned int , bool = true );

double muonIsoValuePF        ( unsigned int imu, unsigned int idavtx, float coner=0.4, float minptn=1.0, float dzcut=0.1);

///////////////////////
// Cosmics Rejection //
///////////////////////

bool isCosmics(unsigned int index);

/////////////////////////////
// Muon d0 corrected by PV //
/////////////////////////////

double mud0PV         (unsigned int index);
double mud0PV_wwV1    (unsigned int index);
double mudzPV_wwV1    (unsigned int index);
double mud0PV_smurfV3 (unsigned int index);
double mudzPV_smurfV3 (unsigned int index);

/////////////////////////////////////////////////////////////////
// checks if muon is also pfmuon, and pfmuon pt = reco muon pt //
/////////////////////////////////////////////////////////////////

bool isPFMuon( int index , bool requireSamePt = true , float dpt_max = 1.0 );

#endif

