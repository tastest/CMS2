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
// 26 April 2011
// standard event cleaning
// for 2011 OS analysis
//
bool cleaning_standardApril2011()
{
    if (!cleaning_goodDAVertexApril2011()) return false;
    if (!cleaning_goodTracks())            return false;
    
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
// 26 April 2011
// at least 1 good DA vertex
// 
bool cleaning_goodDAVertexApril2011()
{             
    int nGoodDAVertex = 0;
    for (size_t v = 0; v < cms2.davtxs_position().size(); ++v) 
    {

      if(isGoodDAVertex(v))
        nGoodDAVertex ++;
    }
    if (nGoodDAVertex == 0) return false;
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
// function to select a good vertex
// 
bool isGoodDAVertex(size_t ivtx) {

  if (cms2.davtxs_isFake()[ivtx]) return false;
  if (cms2.davtxs_ndof()[ivtx] < 4.) return false;
  if (cms2.davtxs_position()[ivtx].Rho() > 2.0) return false;
  if (fabs(cms2.davtxs_position()[ivtx].Z()) > 24.0) return false;
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
    if(fabs(lt_vz - cms2.vtxs_position()[v].Z()) > 1.)
      continue;
    if(fabs(ll_vz - cms2.vtxs_position()[v].Z()) > 1.)
      continue;

    //if we've gotten here, then the vertex is good
    //and both leptons belong to it
    return true;
  }

  return false;
}

//----------------------------------------------------------------
// checks whether the leptons of a given
// hypothesis come from the same good vertex
// by checking if both leptons are within dz
// of 0.2 cm of the same PV and if that PV is
// the closest vertex to each lepton
//----------------------------------------------------------------
int hypsFromSameVtx2011(size_t hypIdx, float dz, bool useDAvtxs, bool requireClosest)
{
    int lt_trkidx = -1;
    int ll_trkidx = -1;
    bool lt_isGsf = false;
    bool ll_isGsf = false;

    if (abs(cms2.hyp_lt_id()[hypIdx]) == 13)
        lt_trkidx = cms2.mus_trkidx()[cms2.hyp_lt_index()[hypIdx]];
    if (abs(cms2.hyp_lt_id()[hypIdx]) == 11) {
        lt_trkidx = cms2.els_trkidx()[cms2.hyp_lt_index()[hypIdx]] > -1 ? cms2.els_trkidx()[cms2.hyp_lt_index()[hypIdx]] : cms2.els_gsftrkidx()[cms2.hyp_lt_index()[hypIdx]];
        lt_isGsf  = cms2.els_trkidx()[cms2.hyp_lt_index()[hypIdx]] > -1 ? false : true;
    }
    if (abs(cms2.hyp_ll_id()[hypIdx]) == 13)
        ll_trkidx = cms2.mus_trkidx()[cms2.hyp_ll_index()[hypIdx]];
    if (abs(cms2.hyp_ll_id()[hypIdx]) == 11) {
        ll_trkidx = cms2.els_trkidx()[cms2.hyp_ll_index()[hypIdx]] > -1 ? cms2.els_trkidx()[cms2.hyp_ll_index()[hypIdx]] : cms2.els_gsftrkidx()[cms2.hyp_ll_index()[hypIdx]];
        ll_isGsf  = cms2.els_trkidx()[cms2.hyp_ll_index()[hypIdx]] > -1 ? false : true;
    }

    if (lt_trkidx < 0 || ll_trkidx < 0)
        return -1;

    // figure out which vertex collection to use
    std::vector<LorentzVector> vtxP4s = useDAvtxs ? cms2.davtxs_position() : cms2.vtxs_position();

    if (!requireClosest) {
        for (size_t v = 0; v < vtxP4s.size(); ++v)  {

            bool vertexIsGood = useDAvtxs ? isGoodDAVertex(v) : isGoodVertex(v);
            if (!vertexIsGood)
                continue;
            if (lt_isGsf) {
                if (gsftrks_dz_pv(lt_trkidx, v, useDAvtxs).first > dz)
                    continue;
            }
            else {
                if (trks_dz_pv(lt_trkidx, v, useDAvtxs).first > dz)
                    continue;
            }
            if (ll_isGsf) {
                if (gsftrks_dz_pv(ll_trkidx, v, useDAvtxs).first > dz)
                    continue;
            }
            else {
                if (trks_dz_pv(ll_trkidx, v, useDAvtxs).first > dz)
                    continue;
            }
            //if we've gotten here, then the vertex is good
            //and both leptons belong to it
            return v;
        }

        return -1;
    }

    float lt_dz = 999.;
    float ll_dz = 999.;
    int lt_vidx = -999;
    int ll_vidx = -999;

    for (unsigned int vtxi = 0; vtxi < vtxP4s.size(); vtxi++) {
        
        bool vertexIsGood = useDAvtxs ? isGoodDAVertex(vtxi) : isGoodVertex(vtxi);
        if (!vertexIsGood)
            continue;

        // first take care of lt
        if (lt_isGsf) {
            if (fabs(gsftrks_dz_pv(lt_trkidx, vtxi, useDAvtxs).first) < lt_dz) {
                lt_dz = fabs(gsftrks_dz_pv(lt_trkidx, vtxi, useDAvtxs).first);
                lt_vidx = vtxi;
            }            
        }
        else {
            if (fabs(trks_dz_pv(lt_trkidx, vtxi, useDAvtxs).first) < lt_dz) {
                lt_dz = fabs(trks_dz_pv(lt_trkidx, vtxi, useDAvtxs).first);
                lt_vidx = vtxi;
            }
        }

        // now same thing for ll
        if (ll_isGsf) {
            if (fabs(gsftrks_dz_pv(ll_trkidx, vtxi, useDAvtxs).first) < ll_dz) {
                ll_dz = fabs(gsftrks_dz_pv(ll_trkidx, vtxi, useDAvtxs).first);
                ll_vidx = vtxi;
            }            
        }
        else {
            if (fabs(trks_dz_pv(ll_trkidx, vtxi, useDAvtxs).first) < ll_dz) {
                ll_dz = fabs(trks_dz_pv(ll_trkidx, vtxi, useDAvtxs).first);
                ll_vidx = vtxi;
            }
        }
    } // end loop over vertices

    if (lt_vidx < 0 || ll_vidx < 0)
        return -1;

    if (lt_vidx != ll_vidx)
        return -1;

    if (fabs(lt_dz) > dz || fabs(ll_dz) > dz)
        return -1;

    return lt_vidx;
}
