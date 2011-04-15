#ifndef EVENTSELECTIONS_H
#define EVENTSELECTIONS_H
#include <sys/types.h>
//----------------------------------------------------------------
// A ridicolusly simple function, but since the Z veto is used 
// in two places, might as well centralize it to keep consistency
//----------------------------------------------------------------
bool inZmassWindow (float mass);

//----------------------------------------------------------------
// standard event cleaning
// for data
//----------------------------------------------------------------
bool cleaning_standard(bool isData);

//----------------------------------------------------------------
// 20 October 2010
// standard event cleaning used for SS analysis
//----------------------------------------------------------------
bool cleaning_standardOctober2010();

//----------------------------------------------------------------
// 5 August 2010
// standard event cleaning
// for low pt dilepton / fake rate data studies
//----------------------------------------------------------------
bool cleaning_standardAugust2010(bool isdata);

//----------------------------------------------------------------
// standard event cleaning
// for low pt dilepton / fake rate data studies
//----------------------------------------------------------------
bool cleaning_standardNoBSC(bool isData);

//----------------------------------------------------------------
// require BPTX
//----------------------------------------------------------------
bool cleaning_BPTX(bool isData);

//----------------------------------------------------------------
// require bit 40 or 41 passed
//----------------------------------------------------------------
bool cleaning_BSC();

//----------------------------------------------------------------
// require bits 36-39 DIDN't pass
// to reject beam halo
//----------------------------------------------------------------
bool cleaning_beamHalo();

// ----------------------------------------------------------------
// 5 August 2010
// at least 1 good vertex
// z position increased from 15 to 24 cm
// ----------------------------------------------------------------
bool cleaning_goodVertexAugust2010();

//----------------------------------------------------------------
// if >= 10 tracks, require at least 25% high purity
//----------------------------------------------------------------
bool cleaning_goodTracks();

//----------------------------------------------------------------
// checks whether a vertex is good or not
//----------------------------------------------------------------
bool isGoodVertex(size_t ivtx);
bool isGoodDAVertex(size_t ivtx);

//----------------------------------------------------------------
// checks whether the leptons of a given
// hypothesis come from the same good vertex
// by checking if both leptons are within dz
// of 1cm of the same PV
//----------------------------------------------------------------
bool hypsFromSameVtx(size_t hypIdx);

//----------------------------------------------------------------
// checks whether the leptons of a given
// hypothesis come from the same good vertex
// by checking if both leptons are within dz
// of 0.2 cm of the same PV and if that PV is
// the closest vertex to each lepton
//----------------------------------------------------------------
bool hypsFromSameVtx2011(size_t hypIdx, float dz = 0.2, bool requireClosest = true);


#endif

