#include <math.h>
#include "TFile.h"
#include "TH2.h"
#include "TVector3.h"
#include "CMS2.h"

static TFile *metcorr_file = TFile::Open("data/metcorr.root", "read");
static TH2D *rfhist = dynamic_cast<TH2D*>(metcorr_file->Get("rf_pt_mbin"));

TVector3 correctMETforTracks() 
{
  // initialize temporary met variables
  double met_x = 0;
  double met_y = 0;

  // set kinematic and quality cuts
  const int    nhits_cut  = 7;
  const double hoe_cut    = 0.1;
  const double d0_cut     = 0.05;
  const double nchisq_cut = 5.;
  const double eta_cut    = 2.4;
  const double low_pt_cut = 2.;
  const double hi_pt_cut  = 100.;

  for( unsigned int trkCount = 0; trkCount < cms2.trks_trk_p4().size(); ++trkCount ) {
    // skip track if matched to a global muon
    if( cms2.trk_musidx()[trkCount] != -999 && cms2.trk_musdr()[trkCount] < 0.1) {
      if( cms2.mus_trkidx()[cms2.trk_musidx()[trkCount]]  == trkCount && cms2.mus_trkdr()[cms2.trk_musidx()[trkCount]] < 0.1 ) continue;
    }

    // skip track if matched to an "electron"
    if( cms2.trk_elsidx()[trkCount] != -999 && cms2.trk_elsdr()[trkCount] < 0.1) {
      if( cms2.els_hOverE()[ cms2.trk_elsidx()[trkCount] ] < hoe_cut ) {
	if( cms2.els_trkidx()[cms2.trk_elsidx()[trkCount]]  == trkCount && cms2.els_trkdr()[cms2.trk_elsidx()[trkCount]] < 0.1 ) continue;
      }
    }

    // skip tracks at large eta or with large pt
    if( fabs( cms2.trks_trk_p4()[trkCount].eta() ) > eta_cut || cms2.trks_trk_p4()[trkCount].pt() > hi_pt_cut ) continue;

    // skip tracks that do no pass quality cuts
    if( cms2.trks_validHits()[trkCount] < nhits_cut || ( cms2.trks_chi2()[trkCount] / cms2.trks_ndof()[trkCount] ) > nchisq_cut || fabs( cms2.trks_d0()[trkCount] ) > d0_cut ) continue;

    // correct tracks w/ pt < 2 setting RF = 0
    if( cms2.trks_trk_p4()[trkCount].pt() < low_pt_cut ) {
      met_x -= cms2.trks_trk_p4()[trkCount].pt() * cos( cms2.trks_trk_p4()[trkCount].phi() );
      met_y -= cms2.trks_trk_p4()[trkCount].pt() * sin( cms2.trks_trk_p4()[trkCount].phi() );
      continue;
    }

    // skip any remaining tracks that don't have outerEta, outerPhi information
    if( cms2.trks_outerEta()[trkCount] == -999 || cms2.trks_outerPhi()[trkCount] == -999 ) continue;

    // if we've made it this far, get the response from the histogram
    int bin   = rfhist->FindBin( cms2.trks_trk_p4()[trkCount].eta(), cms2.trks_trk_p4()[trkCount].pt() );
    double rf = rfhist->GetBinContent( bin );

    // now, correct MET for track using RF
    met_x +=  ( rf * cms2.trks_trk_p4()[trkCount].P() * ( 1 / cosh( cms2.trks_outerEta()[trkCount] ) ) * cos( cms2.trks_outerPhi()[trkCount] ) - cms2.trks_trk_p4()[trkCount].pt() * cos( cms2.trks_trk_p4()[trkCount].phi() ) );
    met_y +=  ( rf * cms2.trks_trk_p4()[trkCount].P() * ( 1 / cosh( cms2.trks_outerEta()[trkCount] ) ) * sin( cms2.trks_outerPhi()[trkCount] ) - cms2.trks_trk_p4()[trkCount].pt() * sin( cms2.trks_trk_p4()[trkCount].phi() ) );
  }

  // fill MET vector
  TVector3 metvec( met_x, met_y, 0);
  
  // return corrected MET
  return metvec;
}
