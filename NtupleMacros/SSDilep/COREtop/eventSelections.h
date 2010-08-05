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
// fake rate event cleaning
// for data
//----------------------------------------------------------------
bool cleaning_fakerate(bool isData);

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

//----------------------------------------------------------------
// if >= 10 tracks, require at least 25% high purity
//----------------------------------------------------------------
bool cleaning_goodTracks();

#endif

