#ifndef BTAGRESHAPHING_H
#define BTAGRESHAPHING_H
#include <utility>
#include <math.h>

#include <TH1F.h>
#include <TFile.h>
#include <TGraph.h>
#include <vector>
#include <string>
#include <iostream>
#include <TSpline.h>
#include "Math/Interpolator.h"
#define MAXPOINTS 200

#include "btag_payload_b.h"
#include "btag_payload_light.h"


class BTagShape 
{
public: 
  
  BTagShape();
  BTagShape( TFile *file, const std::string name, 
			 const std::vector < std::pair <float, float> > &cutsAndSF, 
			 float boundX, float boundY );

  ~BTagShape();

  float eval( float x );

private:
  ROOT::Math::Interpolator * m_i;

};

class EtaPtBin
{
public:
  EtaPtBin();
  EtaPtBin( float emin, float emax, float ptmin, float ptmax );

  ~EtaPtBin();

  bool contains( float eta, float pt );
  float centerEta();
  float centerPt();

  float etaMin;
  float etaMax;
  float ptMin;
  float ptMax;

};

class BinnedBTagShape
{
public:
  BinnedBTagShape();
  BinnedBTagShape( std::vector<EtaPtBin> &bins, 
				   std::vector< std::vector<std::pair<float, float> > > &cutsAndSF, 
				   TFile * f, const std::string name,
				   float boundX, float boundY );
  
  ~BinnedBTagShape();

  float eval(float eta,float pt,float x); 

  std::vector<BTagShape> m_shapes;
  std::vector<EtaPtBin> m_bins; 

};

class BTagShapeInterface
{
public:
  BTagShapeInterface();
  BTagShapeInterface( const std::string file, 
					  float scaleBC, float scaleL, 
					  bool use4points = false, 
					  float boundX = 1.001, float boundY = 1.001,
					  unsigned int maxbins = 9999);

  ~BTagShapeInterface();

  float reshape( float eta, float pt, float csv, int flav);
 
  TFile * m_file; 
  BinnedBTagShape * m_b;
  BinnedBTagShape * m_c;
  BinnedBTagShape * m_l;
 
};

#endif
