#ifndef PH_METUTILITIES_H
#define PH_METUTILITIES_H

#include <utility>
#include <vector>
#include "Math/LorentzVector.h"
#include "Math/Vector3D.h"
#include "TMath.h"

class MetUtilities {
 public:
//  typedef math::XYZVector         Vector;
//  typedef math::XYZTLorentzVector LorentzVector;
  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > Vector;
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
  struct JetInfo {
    LorentzVector p4;
    double        mva;
    double        neutFrac;  
  };

  MetUtilities();
  virtual ~MetUtilities();

  bool              passMVA  (std::pair<LorentzVector,double> iJet);
  LorentzVector    *leadPt   (std::vector<JetInfo> &iJets,bool iFirst);
  int               NJets    (std::vector<JetInfo> &iJets,double iPt);
  double            deltaR   (LorentzVector &iVec1,LorentzVector &iVec2);
  void              cleanJets(std::vector<LorentzVector> &iVis,std::vector<JetInfo> &iJets);
  std::pair<LorentzVector,double> TKMet   (std::vector<std::pair<LorentzVector,double> > &iCands,double iDZ,int iLowDz);
  std::pair<LorentzVector,double> JetMet  (std::vector<JetInfo> &iJets ,bool iPassMVA);
  std::pair<LorentzVector,double> NoPUMet (std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<JetInfo> &iJets,double iDZ);
  std::pair<LorentzVector,double> PUMet   (std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<JetInfo> &iJets,double iDZ);
  std::pair<LorentzVector,double> PUCMet  (std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<JetInfo> &iJets,double iDZ);

  std::pair<LorentzVector,double> PFRecoil  (double iSumEt,LorentzVector iVis,std::vector<std::pair<LorentzVector,double> > &iCands,double iDZ);
  std::pair<LorentzVector,double> TKRecoil  (double iSumEt,LorentzVector iVis,std::vector<std::pair<LorentzVector,double> > &iCands,double iDZ);
  std::pair<LorentzVector,double> NoPURecoil(double iSumEt,LorentzVector iVis,
					     std::vector<std::pair<LorentzVector,double> > &iCands,std::vector<JetInfo> &iJets,double iDZ);
  std::pair<LorentzVector,double> PUCRecoil(double iSumEt,LorentzVector iVis,
					    std::vector<std::pair<LorentzVector,double> > &iCands ,std::vector<JetInfo> &iJets,double iDZ);
 protected:
  // PU jet identifier 
  double fMVACut[3][4][4];  //Jet Id MVA
};

#endif
