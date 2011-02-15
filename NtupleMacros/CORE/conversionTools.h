#include <iostream>
#include <vector>

#include "Math/Point3D.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > XYZPoint;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class ConversionInfo
{
public:

  ConversionInfo()
  {
    dist_ = -9999.;
    dcot_ = -9999.;
    radiusOfConversion_ = -9999.;
    pointOfConversion_ = XYZPoint(-9999.,-9999.,-9999);
    ctfPartnerIndex_ = -9999;
    gsfPartnerIndex_ = -9999;
    deltaMissingHits_ = -9999;
    flag_ = -9999;
  }
  ~ConversionInfo() {}
  ConversionInfo
  ( double dist, double dcot,
    double radiusOfConversion, XYZPoint pointOfConversion,      
    int ctfPartnerIndex,
    int gsfPartnerIndex,
    int deltaMissingHits,
    int flag )
  {
    dist_ = dist;
    dcot_ = dcot;
    radiusOfConversion_ = radiusOfConversion;
    pointOfConversion_ = pointOfConversion;
    ctfPartnerIndex_ = ctfPartnerIndex;
    gsfPartnerIndex_  = gsfPartnerIndex;
    deltaMissingHits_ = deltaMissingHits;
    flag_ = flag;
  }
  double dist() const {return dist_ ; }
  double dcot() const {return dcot_ ; }
  double radiusOfConversion() const { return radiusOfConversion_ ; }
  int ctfPartnerIndex() const { return ctfPartnerIndex_; } 
  int gsfPartnerIndex() const { return gsfPartnerIndex_; }
  XYZPoint pointOfConversion() const { return pointOfConversion_ ; }
     

  //if the partner track is found in the  CTF track collection,
  //we return a ref to the CTF partner track
  int deltaMissingHits() const { return deltaMissingHits_ ; }
  
  //if(flag == 0) //Partner track found in the CTF collection using the electron's CTF track
  //if(flag == 1) //Partner track found in the CTF collection using the electron's GSF track
  //if(flag == 2) //Partner track found in the GSF collection using the electron's CTF track
  //if(flag == 3) //Partner track found in the GSF collection using the electron's GSF track
  
  int flag() const { return flag_ ; }
private:

  double dist_;
  double dcot_;
  double radiusOfConversion_;
  XYZPoint pointOfConversion_;
  int ctfPartnerIndex_;
  int gsfPartnerIndex_;
  int deltaMissingHits_;
  int flag_;

} ;





ConversionInfo getConversionInfo(const LorentzVector el_tk_p4, 
				 const int el_charge,
				 const float el_d0,
				 const float el_dz,
				 const LorentzVector cand_p4,
				 const int cand_charge,
				 const float cand_d0,
				 const double bFieldAtOrigin);


ConversionInfo arbitrateConversionPartnersbyR(const std::vector<ConversionInfo>& v_convCandidates);


ConversionInfo findBestConversionMatch(const std::vector<ConversionInfo>& v_convCandidates);

std::vector<ConversionInfo> getConversionInfos(const int gsfElectronIdx,				
					       const double bFieldAtOrigin,
					       const double minFracSharedHits);
