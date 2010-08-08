#ifndef EVENTSELECTIONS_H
#define EVENTSELECTIONS_H

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
// 5 August 2010
// standard event cleaning
// for low pt dilepton / fake rate data studies
//----------------------------------------------------------------
bool cleaning_standardAugust2010(bool isData);

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

//----------------------------------------------------------------
// at least 1 good vertex
// ----------------------------------------------------------------
bool cleaning_goodVertex();

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

#endif

