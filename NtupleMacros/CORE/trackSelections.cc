#include <math.h>
#include "CMS2.h"
#include "trackSelections.h"

// return a pair of d0, d0err with respect to a primary vertex
std::pair<double, double> trks_d0_pv (int itrk, int ipv)
{
     const double bx = cms2.vtxs_position().at(ipv).x();
     const double by = cms2.vtxs_position().at(ipv).y();
     // assume the layout of the covariance matrix is (Vxx, Vxy, Vxz)
     //						      (Vyx, Vyy, ...)
     const double vxx = cms2.vtxs_covMatrix().at(ipv).at(0);
     const double vxy = cms2.vtxs_covMatrix().at(ipv).at(1);
     const double vyy = cms2.vtxs_covMatrix().at(ipv).at(4);
     const double phi = cms2.trks_trk_p4().at(itrk).phi();
     const double d0vtx = cms2.trks_d0().at(itrk) - bx * sin(phi) + by * cos(phi);
     const double d0err = cms2.trks_d0Err().at(itrk);
     const double phierr = cms2.trks_phiErr().at(itrk);
     const double d0phicov = cms2.trks_d0phiCov().at(itrk);
     // we will let the optimizer take care of subexpression
     // elimination for this one...
     const double d0err2vtx = d0err * d0err 
	  - 2 * (bx * cos(phi) + by * sin(phi)) * d0phicov
	  + (bx * cos(phi) + by * sin(phi)) * (bx * cos(phi) + by * sin(phi)) * phierr * phierr
	  + sin(phi) * sin(phi) * vxx + cos(phi) * cos(phi) * vyy
	  - 2 * sin(phi) * cos(phi) * vxy;
     if (d0err2vtx >= 0) 
	  return std::pair<double, double>(d0vtx, sqrt(d0err2vtx));
     std::cerr << "Oh no!  sigma^2(d0corr) < 0!" << std::endl;
     return std::pair<double, double>(d0vtx, -sqrt(-d0err2vtx));
}

std::pair<double, double> trks_dz_pv (int itrk, int ipv)
{
  LorentzVector tkp = cms2.trks_trk_p4().at(itrk);
  LorentzVector tkv = cms2.trks_vertex_p4().at(itrk);
  LorentzVector pv = cms2.vtxs_position().at(ipv);
  double pt = cms2.trks_trk_p4().at(itrk).pt();
  double phi = cms2.trks_trk_p4().at(itrk).phi();
  double theta = cms2.trks_trk_p4().at(itrk).theta();

  double ddzdpvx = cos(phi)*1./tan(theta);
  double ddzdpvy = sin(phi)*1./tan(theta);
  double ddzdphi = -1*pv.x()*sin(phi)*1./tan(theta) + pv.y()*cos(phi)*1./tan(theta);
  double ddzdtheta = -1*1/sin(theta)*1/sin(theta) * (pv.x()*cos(phi) + pv.y()*sin(phi));
  ddzdpvx   *= ddzdpvx;
  ddzdpvy   *= ddzdpvy;
  ddzdphi   *= ddzdphi;
  ddzdtheta *= ddzdtheta;
  double z0Err  = cms2.trks_z0Err().at(itrk);
  double phiErr = cms2.trks_phiErr().at(itrk);
  double thetaErr = cms2.trks_etaErr().at(itrk)*sin(theta);
  double pvxErr = cms2.vtxs_xError().at(ipv);
  double pvyErr = cms2.vtxs_yError().at(ipv);
  double pvzErr = cms2.vtxs_zError().at(ipv);
  z0Err    *= z0Err;
  phiErr   *= phiErr;
  thetaErr *= thetaErr;
  pvxErr *= pvxErr;
  pvyErr *= pvyErr;
  pvzErr *= pvzErr;

  double value = cms2.trks_z0().at(itrk) - pv.z() + (pv.x()*cos(phi) + pv.y()*sin(phi) )*1./tan(theta);
  //note that the error does not account for correlations since we do not store the track covariance matrix
  double error = sqrt(z0Err + pvzErr + ddzdpvx*pvxErr + ddzdpvy*pvyErr + ddzdphi*phiErr + ddzdtheta*thetaErr);

  return std::pair<double, double>(value, error);
}

std::pair<double, double> gsftrks_dz_pv (int itrk, int ipv)
{
  LorentzVector tkp = cms2.gsftrks_p4().at(itrk);
  LorentzVector tkv = cms2.gsftrks_vertex_p4().at(itrk);
  LorentzVector pv = cms2.vtxs_position().at(ipv);
  double pt = cms2.gsftrks_p4().at(itrk).pt();
  double phi = cms2.gsftrks_p4().at(itrk).phi();
  double theta = cms2.gsftrks_p4().at(itrk).theta();

  double ddzdpvx = cos(phi)*1./tan(theta);
  double ddzdpvy = sin(phi)*1./tan(theta);
  double ddzdphi = -1*pv.x()*sin(phi)*1./tan(theta) + pv.y()*cos(phi)*1./tan(theta);
  double ddzdtheta = -1*1/sin(theta)*1/sin(theta) * (pv.x()*cos(phi) + pv.y()*sin(phi));
  ddzdpvx   *= ddzdpvx;
  ddzdpvy   *= ddzdpvy;
  ddzdphi   *= ddzdphi;
  ddzdtheta *= ddzdtheta;
  double z0Err  = cms2.gsftrks_z0Err().at(itrk);
  double phiErr = cms2.gsftrks_phiErr().at(itrk);
  double thetaErr = cms2.gsftrks_etaErr().at(itrk)*sin(theta);
  double pvxErr = cms2.vtxs_xError().at(ipv);
  double pvyErr = cms2.vtxs_yError().at(ipv);
  double pvzErr = cms2.vtxs_zError().at(ipv);
  z0Err    *= z0Err;
  phiErr   *= phiErr;
  thetaErr *= thetaErr;
  pvxErr *= pvxErr;
  pvyErr *= pvyErr;
  pvzErr *= pvzErr;

  double value = cms2.gsftrks_z0().at(itrk) - pv.z() + (pv.x()*cos(phi) + pv.y()*sin(phi) )*1./tan(theta);
  //note that the error does not account for correlations since we do not store the track covariance matrix
  double error = sqrt(z0Err + pvzErr + ddzdpvx*pvxErr + ddzdpvy*pvyErr + ddzdphi*phiErr + ddzdtheta*thetaErr);

  return std::pair<double, double>(value, error);
}


//----------------------------------------------------------------
// Simple function that tells you whether or not a track passed 
// a particular quality flag.
//----------------------------------------------------------------
bool isTrackQuality( int index, int cuts ) {
     return ( ( cms2.trks_qualityMask().at(index) & cuts ) == cuts );
}

///////////////////////////////////////////////////////////////////
//////////////////////////////DA VTX///////////////////////////////
///////////////////////////////////////////////////////////////////

// return a pair of d0, d0err with respect to a primary vertex
std::pair<double, double> trks_d0_dapv (int itrk, int idapv)
{
     const double bx = cms2.davtxs_position().at(idapv).x();
     const double by = cms2.davtxs_position().at(idapv).y();
     // assume the layout of the covariance matrix is (Vxx, Vxy, Vxz)
     //						      (Vyx, Vyy, ...)
     const double vxx = cms2.davtxs_covMatrix().at(idapv).at(0);
     const double vxy = cms2.davtxs_covMatrix().at(idapv).at(1);
     const double vyy = cms2.davtxs_covMatrix().at(idapv).at(4);
     const double phi = cms2.trks_trk_p4().at(itrk).phi();
     const double d0vtx = cms2.trks_d0().at(itrk) - bx * sin(phi) + by * cos(phi);
     const double d0err = cms2.trks_d0Err().at(itrk);
     const double phierr = cms2.trks_phiErr().at(itrk);
     const double d0phicov = cms2.trks_d0phiCov().at(itrk);
     // we will let the optimizer take care of subexpression
     // elimination for this one...
     const double d0err2vtx = d0err * d0err 
	  - 2 * (bx * cos(phi) + by * sin(phi)) * d0phicov
	  + (bx * cos(phi) + by * sin(phi)) * (bx * cos(phi) + by * sin(phi)) * phierr * phierr
	  + sin(phi) * sin(phi) * vxx + cos(phi) * cos(phi) * vyy
	  - 2 * sin(phi) * cos(phi) * vxy;
     if (d0err2vtx >= 0) 
	  return std::pair<double, double>(d0vtx, sqrt(d0err2vtx));
     std::cerr << "Oh no!  sigma^2(d0corr) < 0!" << std::endl;
     return std::pair<double, double>(d0vtx, -sqrt(-d0err2vtx));
}

std::pair<double, double> trks_dz_dapv (int itrk, int idapv)
{
  LorentzVector tkp = cms2.trks_trk_p4().at(itrk);
  LorentzVector tkv = cms2.trks_vertex_p4().at(itrk);
  LorentzVector dapv = cms2.davtxs_position().at(idapv);
  double pt = cms2.trks_trk_p4().at(itrk).pt();
  double phi = cms2.trks_trk_p4().at(itrk).phi();
  double theta = cms2.trks_trk_p4().at(itrk).theta();

  double ddzddapvx = cos(phi)*1./tan(theta);
  double ddzddapvy = sin(phi)*1./tan(theta);
  double ddzdphi = -1*dapv.x()*sin(phi)*1./tan(theta) + dapv.y()*cos(phi)*1./tan(theta);
  double ddzdtheta = -1*1/sin(theta)*1/sin(theta) * (dapv.x()*cos(phi) + dapv.y()*sin(phi));
  ddzddapvx   *= ddzddapvx;
  ddzddapvy   *= ddzddapvy;
  ddzdphi   *= ddzdphi;
  ddzdtheta *= ddzdtheta;
  double z0Err  = cms2.trks_z0Err().at(itrk);
  double phiErr = cms2.trks_phiErr().at(itrk);
  double thetaErr = cms2.trks_etaErr().at(itrk)*sin(theta);
  double dapvxErr = cms2.davtxs_xError().at(idapv);
  double dapvyErr = cms2.davtxs_yError().at(idapv);
  double dapvzErr = cms2.davtxs_zError().at(idapv);
  z0Err    *= z0Err;
  phiErr   *= phiErr;
  thetaErr *= thetaErr;
  dapvxErr *= dapvxErr;
  dapvyErr *= dapvyErr;
  dapvzErr *= dapvzErr;

  double value = cms2.trks_z0().at(itrk) - dapv.z() + (dapv.x()*cos(phi) + dapv.y()*sin(phi) )*1./tan(theta);
  //note that the error does not account for correlations since we do not store the track covariance matrix
  double error = sqrt(z0Err + dapvzErr + ddzddapvx*dapvxErr + ddzddapvy*dapvyErr + ddzdphi*phiErr + ddzdtheta*thetaErr);

  return std::pair<double, double>(value, error);
}

std::pair<double, double> gsftrks_dz_dapv (int itrk, int idapv)
{
  LorentzVector tkp = cms2.gsftrks_p4().at(itrk);
  LorentzVector tkv = cms2.gsftrks_vertex_p4().at(itrk);
  LorentzVector dapv = cms2.davtxs_position().at(idapv);
  double pt = cms2.gsftrks_p4().at(itrk).pt();
  double phi = cms2.gsftrks_p4().at(itrk).phi();
  double theta = cms2.gsftrks_p4().at(itrk).theta();

  double ddzddapvx = cos(phi)*1./tan(theta);
  double ddzddapvy = sin(phi)*1./tan(theta);
  double ddzdphi = -1*dapv.x()*sin(phi)*1./tan(theta) + dapv.y()*cos(phi)*1./tan(theta);
  double ddzdtheta = -1*1/sin(theta)*1/sin(theta) * (dapv.x()*cos(phi) + dapv.y()*sin(phi));
  ddzddapvx   *= ddzddapvx;
  ddzddapvy   *= ddzddapvy;
  ddzdphi   *= ddzdphi;
  ddzdtheta *= ddzdtheta;
  double z0Err  = cms2.gsftrks_z0Err().at(itrk);
  double phiErr = cms2.gsftrks_phiErr().at(itrk);
  double thetaErr = cms2.gsftrks_etaErr().at(itrk)*sin(theta);
  double dapvxErr = cms2.davtxs_xError().at(idapv);
  double dapvyErr = cms2.davtxs_yError().at(idapv);
  double dapvzErr = cms2.davtxs_zError().at(idapv);
  z0Err    *= z0Err;
  phiErr   *= phiErr;
  thetaErr *= thetaErr;
  dapvxErr *= dapvxErr;
  dapvyErr *= dapvyErr;
  dapvzErr *= dapvzErr;

  double value = cms2.gsftrks_z0().at(itrk) - dapv.z() + (dapv.x()*cos(phi) + dapv.y()*sin(phi) )*1./tan(theta);
  //note that the error does not account for correlations since we do not store the track covariance matrix
  double error = sqrt(z0Err + dapvzErr + ddzddapvx*dapvxErr + ddzddapvy*dapvyErr + ddzdphi*phiErr + ddzdtheta*thetaErr);

  return std::pair<double, double>(value, error);
}

