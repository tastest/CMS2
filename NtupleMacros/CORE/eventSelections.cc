
#include "eventSelections.h"
#include "trackSelections.h"
#include "Math/LorentzVector.h"

#include "CMS2.h"

//----------------------------------------------------------------
// A ridicolusly simple function, but since the Z veto is used 
// in two places, might as well centralize it to keep consistency
//----------------------------------------------------------------
bool inZmassWindow (float mass) {
    if (mass > 76. && mass < 106.) return true;
    return false;
}

//
// standard event cleaning
// for data
//
bool cleaning_standard(bool isData)
{
    if (!cleaning_BPTX(isData)) return false;
    if (!cleaning_BSC())        return false;
    if (!cleaning_beamHalo())   return false;
    if (!cleaning_goodVertex()) return false;
    if (!cleaning_goodTracks()) return false;
    return true;
}

//
// require bit 40 or 41 passed
//
bool cleaning_BSC() 
{
    if (!(cms2.l1_techbits2() & (1<<8) || cms2.l1_techbits2() & (1<<9))) return false;
    return true;
}

//
// require bits 36-39 DIDN't pass
// to reject beam halo
//
bool cleaning_beamHalo()
{
    if (cms2.l1_techbits2() & (1<<7) || cms2.l1_techbits2() & (1<<6) ||
        cms2.l1_techbits2() & (1<<5) || cms2.l1_techbits2() & (1<<4)) return false;
    return true;
}

//
// require BPTX
//
bool cleaning_BPTX(bool isData)
{
    if (isData && !(cms2.l1_techbits1() & (1<<0))) return false;
    return true;
}

//
// at least 1 good vertex
// 
bool cleaning_goodVertex()
{             
    int nGoodVertex = 0;
    for (size_t v = 0; v < cms2.vtxs_position().size(); ++v) 
    {
        if (cms2.vtxs_isFake()[v]) continue;
        if (cms2.vtxs_ndof()[v] < 4) continue;
        if (fabs(cms2.vtxs_position()[v].Z()) > 15.0) continue;
        if (cms2.vtxs_position()[v].Rho() > 2.0) continue;
        nGoodVertex ++;
    }
    if (nGoodVertex == 0) return false;
    return true;
}

//
// if >= 10 tracks, require at least 25% high purity
//
bool cleaning_goodTracks()
{
    if (cms2.trks_ndof().size() >= 10) {
        int nHighPurityTracks = 0;
        for (size_t t = 0; t < cms2.trks_ndof().size(); ++t)
        {
            if (isTrackQuality(t, (1<<highPurity))) nHighPurityTracks ++;
        }
        if (cms2.trks_ndof().size()/float(nHighPurityTracks) < 0.25) return false;
    }
    return true;
}

