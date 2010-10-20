
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
    if (!cleaning_goodVertexAugust2010()) return false;
    if (!cleaning_goodTracks()) return false;
    return true;
}

//
// standard event cleaning
// for low pt dilepton / fake rate data studies
//
bool cleaning_standardNoBSC(bool isData)
{
    if (!cleaning_BPTX(isData)) return false;
    if (!cleaning_beamHalo())   return false;
    if (!cleaning_goodVertexAugust2010()) return false;
    if (!cleaning_goodTracks()) return false;
    return true;
}

//
// 20 October 2010
// standard event cleaning
// for 2010 SS analysis
//
bool cleaning_standardOctober2010()
{
    if (!cleaning_goodVertexAugust2010()) return false;
    if (!cleaning_goodTracks()) return false;
    
    return true;
}

//
// 5 August 2010
// standard event cleaning
// for low pt dilepton / fake rate data studies
//
bool cleaning_standardAugust2010(bool isdata)
{
    if (!cleaning_goodVertexAugust2010()) return false;
    if (!cleaning_goodTracks()) return false;
    if(isdata) {
      if (cms2.evt_hbheFilter()==0) 
	return false;
    } else {
      if(cms2.hcalnoise_minE2Over10TS()<0.7) return false;
      //if(cms2.hcalnoise_maxE2Over10TS()>maxRatio_) return false; //don't have this in the MC ntuples :(
      if(cms2.hcalnoise_maxHPDHits()>=17) return false;
      if(cms2.hcalnoise_maxRBXHits()>=999) return false;
      //if(cms2.hcalnoise_maxHPDNoOtherHits()>=10) return false; //don't have this in the MC ntuples :(
      if(cms2.hcalnoise_maxZeros()>=10) return false;
      if(cms2.hcalnoise_min25GeVHitTime()<-9999.0) return false;
      if(cms2.hcalnoise_max25GeVHitTime()>9999.0) return false;
      if(cms2.hcalnoise_minRBXEMF()<-999.0) return false;
    }
    
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
// 5 August 2010
// at least 1 good vertex
// 
bool cleaning_goodVertexAugust2010()
{             
    int nGoodVertex = 0;
    for (size_t v = 0; v < cms2.vtxs_position().size(); ++v) 
    {

      if(isGoodVertex(v))
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
        if (float(nHighPurityTracks)/cms2.trks_ndof().size() < 0.25) return false;
    }
    return true;
}

//
// function to select a good vertex
// 
bool isGoodVertex(size_t ivtx) {

  if (cms2.vtxs_isFake()[ivtx]) return false;
  if (cms2.vtxs_ndof()[ivtx] < 4.) return false;
  if (cms2.vtxs_position()[ivtx].Rho() > 2.0) return false;
  if (fabs(cms2.vtxs_position()[ivtx].Z()) > 24.0) return false;
  return true;

}



//
// function to check whether or not both the hypotheses 
// are from the same vertex
//

bool hypsFromSameVtx(size_t hypIdx) {

  Float_t lt_vz = 9999.;
  Float_t ll_vz = 9999.;

  if(abs(cms2.hyp_lt_id()[hypIdx]) == 11)
    lt_vz = cms2.els_vertex_p4()[cms2.hyp_lt_index()[hypIdx]].Z();
  if(abs(cms2.hyp_lt_id()[hypIdx]) == 13)
    lt_vz = cms2.mus_vertex_p4()[cms2.hyp_lt_index()[hypIdx]].Z();

  if(abs(cms2.hyp_ll_id()[hypIdx]) == 11)
    ll_vz = cms2.els_vertex_p4()[cms2.hyp_ll_index()[hypIdx]].Z();
  if(abs(cms2.hyp_ll_id()[hypIdx]) == 13)
    ll_vz = cms2.mus_vertex_p4()[cms2.hyp_ll_index()[hypIdx]].Z();


  for (size_t v = 0; v < cms2.vtxs_position().size(); ++v)  {
    if(!isGoodVertex(v))
      continue;
    if(fabs(lt_vz - cms2.vtxs_position()[v].Z()) > 1)
      continue;
    if(fabs(ll_vz - cms2.vtxs_position()[v].Z()) > 1)
      continue;

    //if we've gotten here, then the vertex is good
    //and both leptons belong to it
    return true;
  }

  return false;
}

