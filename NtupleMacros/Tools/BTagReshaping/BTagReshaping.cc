// root -l BTagReshaping.C+
// reshape(1.2, 30, 0.3, 5)
#include <utility>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TGraph.h>
#include <TSpline.h>

#include "TFile.h"
#include "TH1F.h"
#include "Math/Interpolator.h"

#include "BTagReshaping.h"
#include "btag_payload_b.h"
#include "btag_payload_light.h"

// BTagShapeInterface sh("csvdiscr.root",0,0);

BTagShapeInterface& sh()
{
  static BTagShapeInterface sh_local("csvdiscr.root",0,0);
  return sh_local;
}


float reshape(float eta, float pt, float csv, int flav)
{
  return  sh().reshape(eta, pt, csv, flav);

}

//Libraries for BTagShape class
BTagShape::BTagShape()
{
  m_i = NULL;
}

BTagShape::BTagShape( TFile *file , const std::string name, const std::vector< std::pair<float, float> > &cutsAndSF, float boundX, float boundY )
{
  TH1F * m_h = (TH1F *) file->Get(name.c_str());

  //compute equivalents
  std::vector < std::pair <float,float> > eq;
  int lastbin =2001;
  float integral =  m_h->Integral( -1, lastbin );

  for( size_t i = 0; i < cutsAndSF.size(); i++ )
    {
	  float oldCut           = cutsAndSF.at(i).first;
	  float sf               = cutsAndSF.at(i).second;
      float originalIntegral = m_h->Integral(      m_h->FindBin( oldCut ), lastbin );
      float originalLowEdge  = m_h->GetBinLowEdge( m_h->FindBin( oldCut ));
      //std::cout << std::endl<<    " Scale Factor : " << sf << std::endl;
	  //      float target=originalIntegral/sf;
      float target = originalIntegral*sf;
      //std::cout << " Target " << target << " orig " << originalIntegral << std::endl;
      for(int j = lastbin; j > -1; j-- )
		{
		  if( m_h->Integral( j, lastbin ) >= target )
			{
			  //equivalents.push_back(std::pair<float,float>(originalLowEdge,h->GetBinLowEdge(j))); 
			  eq.push_back( std::pair <float,float> ( m_h->GetBinLowEdge(j), originalLowEdge ) ); 
			  //std::cout << "Found at " << j << " was " << m_h->FindBin(oldCut) <<  std::endl;
			  //std::cout << m_h->GetBinLowEdge(j) << " was " << originalLowEdge << " cut: " << oldCut <<  std::endl;
			  break;
			}
		}
 
    }


  //Interpolator

  std::vector <double> x;
  x.push_back(0.0);

  std::vector <double> y;
  y.push_back(0.0);

  for( size_t i = 0 ; i < eq.size(); i++ )
	{
	  x.push_back( eq.at(eq.size()-i-1).first  );
	  y.push_back( eq.at(eq.size()-i-1).second );
	}
  x.push_back(boundX);
  y.push_back(boundY);

  m_i = new ROOT::Math::Interpolator( x, y, ROOT::Math::Interpolation::kLINEAR );
}

float BTagShape::eval(float x) 
{ 
  return m_i->Eval(x); 
} 

BTagShape::~BTagShape()
{
  // delete m_i;
}


//Libraries for EtaPtBin Class
EtaPtBin::EtaPtBin()
{
  etaMin = -1.0;
  etaMax = -1.0;
  ptMin = -1.0;
  ptMax = -1.0; 
}

EtaPtBin::EtaPtBin( float emin, float emax, float ptmin, float ptmax ) 
  :etaMin(emin),
   etaMax(emax),
   ptMin(ptmin),
   ptMax(ptmax) 
{
}

EtaPtBin::~EtaPtBin()
{}

bool EtaPtBin::contains( float eta, float pt ) 
{
  return eta < etaMax && eta >= etaMin && pt < ptMax && pt >= ptMin; 
} 

float EtaPtBin::centerEta() 
{ 
  return (etaMax+etaMin)/2.;
}
 
float EtaPtBin::centerPt() 
{ 
  return (ptMax+ptMin)/2.;
}


//Libraries for BinnedBTagShape class
BinnedBTagShape::BinnedBTagShape()
{}

BinnedBTagShape::BinnedBTagShape( std::vector <EtaPtBin> &bins, 
								  std::vector < std::vector < std::pair <float, float> > > &cutsAndSF, 
								  TFile * f, const std::string name,
								  float boundX, float boundY )
  :m_bins(bins)

{
    
  for( size_t i = 0; i < bins.size(); i++ )
	{
	  m_shapes.push_back( BTagShape( f, name, cutsAndSF[i], boundX, boundY ) );
	}

}

float BinnedBTagShape::eval( float eta, float pt, float x ) 
{
  for( size_t i = 0; i < m_bins.size(); i++)
	{
	  if( m_bins.at(i).contains( fabs(eta), pt ) ) return m_shapes.at(i).eval(x);
	}
  //    std::cout << "Cannot reshape eta pt discr "  << eta << " " << pt << " " << x << std::endl; 
  return x;
}

BinnedBTagShape::~BinnedBTagShape()
{}

//Libraries for BTagShapeInterface class
BTagShapeInterface::BTagShapeInterface()
{
  m_file = NULL; 

  m_b = NULL;
  m_c = NULL;
  m_l = NULL;
}

BTagShapeInterface::BTagShapeInterface( const std::string file, 
										float scaleBC, float scaleL, 
										bool use4points, 
										float boundX, float boundY,
										unsigned int maxbins )
  :m_file(new TFile(file.c_str(), "READ"))
   
{
  
  // beff::seteta();
  // beff::setpt();
  // beff::setbins();
  // beff::setCSVL_SFb_error();
  // beff::setCSVM_SFb_error();
  // beff::setCSVT_SFb_error();
    
  beff::bins = 16;

  beff::ptmin.clear();
  beff::ptmin.push_back(20.0);
  beff::ptmin.push_back(30.0);
  beff::ptmin.push_back(40.0);
  beff::ptmin.push_back(50.0);
  beff::ptmin.push_back(60.0);
  beff::ptmin.push_back(70.0);
  beff::ptmin.push_back(80.0);
  beff::ptmin.push_back(100.0);
  beff::ptmin.push_back(120.0);
  beff::ptmin.push_back(160.0);
  beff::ptmin.push_back(210.0);
  beff::ptmin.push_back(260.0);
  beff::ptmin.push_back(320.0);
  beff::ptmin.push_back(400.0);
  beff::ptmin.push_back(500.0);
  beff::ptmin.push_back(600.0);

  beff::ptmax.clear();
  beff::ptmax.push_back(30.0);
  beff::ptmax.push_back(40.0);
  beff::ptmax.push_back(50.0);
  beff::ptmax.push_back(60.0);
  beff::ptmax.push_back(70.0);
  beff::ptmax.push_back(80.0);
  beff::ptmax.push_back(100.0);
  beff::ptmax.push_back(120.0);
  beff::ptmax.push_back(160.0);
  beff::ptmax.push_back(210.0);
  beff::ptmax.push_back(260.0);
  beff::ptmax.push_back(320.0);
  beff::ptmax.push_back(400.0);
  beff::ptmax.push_back(500.0);
  beff::ptmax.push_back(600.0);
  beff::ptmax.push_back(800.0);

  beff::etamin.clear();
  beff::etamin.push_back(0.0);
  beff::etamin.push_back(0.5);
  beff::etamin.push_back(0.8);
  beff::etamin.push_back(1.0);
  beff::etamin.push_back(1.5);

  beff::etamax.clear();
  beff::etamax.push_back(0.5);
  beff::etamax.push_back(0.8);
  beff::etamax.push_back(1.0);
  beff::etamax.push_back(1.5);
  beff::etamax.push_back(2.5);

  beff::CSVL_SFb_error.clear();
  beff::CSVL_SFb_error.push_back(0.0484285);
  beff::CSVL_SFb_error.push_back(0.0126178);
  beff::CSVL_SFb_error.push_back(0.0120027);
  beff::CSVL_SFb_error.push_back(0.0141137);
  beff::CSVL_SFb_error.push_back(0.0145441);
  beff::CSVL_SFb_error.push_back(0.0131145);
  beff::CSVL_SFb_error.push_back(0.0168479);
  beff::CSVL_SFb_error.push_back(0.0160836);
  beff::CSVL_SFb_error.push_back(0.0126209);
  beff::CSVL_SFb_error.push_back(0.0136017);
  beff::CSVL_SFb_error.push_back(0.019182);
  beff::CSVL_SFb_error.push_back(0.0198805);
  beff::CSVL_SFb_error.push_back(0.0386531);
  beff::CSVL_SFb_error.push_back(0.0392831);
  beff::CSVL_SFb_error.push_back(0.0481008);
  beff::CSVL_SFb_error.push_back(0.0474291); 

  beff::CSVM_SFb_error.clear();
  beff::CSVM_SFb_error.push_back(0.0554504);
  beff::CSVM_SFb_error.push_back(0.0209663);
  beff::CSVM_SFb_error.push_back(0.0207019);
  beff::CSVM_SFb_error.push_back(0.0230073);
  beff::CSVM_SFb_error.push_back(0.0208719);
  beff::CSVM_SFb_error.push_back(0.0200453);
  beff::CSVM_SFb_error.push_back(0.0264232);
  beff::CSVM_SFb_error.push_back(0.0240102);
  beff::CSVM_SFb_error.push_back(0.0229375);
  beff::CSVM_SFb_error.push_back(0.0184615);
  beff::CSVM_SFb_error.push_back(0.0216242);
  beff::CSVM_SFb_error.push_back(0.0248119);
  beff::CSVM_SFb_error.push_back(0.0465748);
  beff::CSVM_SFb_error.push_back(0.0474666);
  beff::CSVM_SFb_error.push_back(0.0718173);
  beff::CSVM_SFb_error.push_back(0.0717567);

  beff::CSVT_SFb_error.clear();
  beff::CSVT_SFb_error.push_back(0.0567059);
  beff::CSVT_SFb_error.push_back(0.0266907);
  beff::CSVT_SFb_error.push_back(0.0263491);
  beff::CSVT_SFb_error.push_back(0.0342831);
  beff::CSVT_SFb_error.push_back(0.0303327);
  beff::CSVT_SFb_error.push_back(0.024608);
  beff::CSVT_SFb_error.push_back(0.0333786);
  beff::CSVT_SFb_error.push_back(0.0317642);
  beff::CSVT_SFb_error.push_back(0.031102);
  beff::CSVT_SFb_error.push_back(0.0295603);
  beff::CSVT_SFb_error.push_back(0.0474663);
  beff::CSVT_SFb_error.push_back(0.0503182);
  beff::CSVT_SFb_error.push_back(0.0580424);
  beff::CSVT_SFb_error.push_back(0.0575776);
  beff::CSVT_SFb_error.push_back(0.0769779);
  beff::CSVT_SFb_error.push_back(0.0898199);

  std::vector <EtaPtBin> binsBC;
  std::vector < std::vector < std::pair <float, float> > > cutsAndSFB;
  std::vector < std::vector < std::pair <float, float> > > cutsAndSFC;
  float charmFactor = 2. -1.0 ; //additional uncertainty for charm
  if( maxbins > beff::bins )  maxbins = beff::bins;


  for( size_t i = 0; i < maxbins; i++ )
	{
	  EtaPtBin bin( -2.5, 2.5, beff::ptmin.at(i), beff::ptmax.at(i) );
	  binsBC.push_back(bin);

	  std::vector < std::pair <float, float> > cutsAndSFbinB;
	  std::vector < std::pair <float, float> > cutsAndSFbinC;

	  float sft = 0.94;
	  if(use4points)
		{
		  sft+=scaleBC * beff::CSVT_SFb_error[i]; // add error
		  cutsAndSFbinB.push_back( std::pair <float, float> ( 0.98, sft ) );
		  sft+=scaleBC * beff::CSVT_SFb_error[i]*charmFactor; // charm additional error
		  cutsAndSFbinC.push_back( std::pair <float, float> ( 0.98, sft ) );
		}

	  sft = beff::CSVT_SFb(bin.centerPt());
	  sft+=scaleBC * beff::CSVT_SFb_error[i]; // add error
	  cutsAndSFbinB.push_back( std::pair <float, float> ( 0.898, sft ) );
	  sft+=scaleBC * beff::CSVT_SFb_error[i]*charmFactor; // charm additional error
	  cutsAndSFbinC.push_back( std::pair <float, float> ( 0.898, sft ) );  

	  float sfm = beff::CSVM_SFb(bin.centerPt());
	  sfm+=scaleBC * beff::CSVM_SFb_error[i]; // add error
	  cutsAndSFbinB.push_back( std::pair <float, float> ( 0.679, sfm ) );
	  sfm+=scaleBC * beff::CSVM_SFb_error[i]*charmFactor; // charm additional error
	  cutsAndSFbinC.push_back( std::pair <float, float> ( 0.679, sfm ) );

	  float sfl = beff::CSVL_SFb(bin.centerPt());
	  sfl+=scaleBC * beff::CSVL_SFb_error[i]; // add error
	  cutsAndSFbinB.push_back(std::pair<float, float>(0.244,sfl));
	  sfl+=scaleBC * beff::CSVL_SFb_error[i]*charmFactor; // charm additional error
	  cutsAndSFbinC.push_back(std::pair<float, float>(0.244,sfl));

	  //std::cout << "SFs "  << i << " " << sfl << " " << sfm << " " << sft << std::endl;
	  cutsAndSFB.push_back(cutsAndSFbinB);
	  cutsAndSFC.push_back(cutsAndSFbinC);
	}  

  //underflow:
  {
	std::vector<std::pair< float, float > > cutsAndSFbinC;
	std::vector<std::pair< float, float > > cutsAndSFbinB;

	binsBC.push_back( EtaPtBin( -2.5, 2.5, -9e99, beff::ptmin[0] ) );
	float sft = beff::CSVT_SFb( beff::ptmin[0] );
	sft += scaleBC * 0.12;// add error
	cutsAndSFbinB.push_back( std::pair <float, float> ( 0.898, sft ) );
	sft += scaleBC * 0.12*charmFactor; // charm additional error
	cutsAndSFbinC.push_back( std::pair <float, float> ( 0.898, sft ) );
   
	float sfm = beff::CSVM_SFb( beff::ptmin[0] );
	sfm += scaleBC * 0.12; // add error
	cutsAndSFbinB.push_back( std::pair <float, float> ( 0.679, sfm ) );
	sfm += scaleBC * 0.12*charmFactor; // charm additional error
	cutsAndSFbinC.push_back( std::pair <float, float> ( 0.679, sfm ) );
      
	float sfl = beff::CSVL_SFb( beff::ptmin[0] );
	sfl += scaleBC * 0.12; // add error
	cutsAndSFbinB.push_back( std::pair <float, float> ( 0.244, sfl ) );
	sfl += scaleBC * 0.12*charmFactor; // charm additional error
	cutsAndSFbinC.push_back( std::pair <float, float> ( 0.244, sfl ) );

	//std::cout << "Firstbin SFs " << sfl << " " << sfm << " " << sft << std::endl;
	cutsAndSFB.push_back(cutsAndSFbinB);
	cutsAndSFC.push_back(cutsAndSFbinC);
  }

  //overflow:
  {
	std::vector<std::pair< float, float > > cutsAndSFbinC;
	std::vector<std::pair< float, float > > cutsAndSFbinB;

	binsBC.push_back( EtaPtBin( -2.5, 2.5, beff::ptmax[maxbins-1], 9e99 ) );
	float sft = beff::CSVT_SFb( beff::ptmax[maxbins-1] );
	sft += scaleBC * beff::CSVT_SFb_error[maxbins-1]*2;// add error
	cutsAndSFbinB.push_back( std::pair <float, float> ( 0.898, sft ) );
	sft += scaleBC * beff::CSVT_SFb_error[maxbins-1]*charmFactor; // charm additional error
	cutsAndSFbinC.push_back( std::pair <float, float> ( 0.898, sft ) );
      
	float sfm = beff::CSVM_SFb( beff::ptmax[maxbins-1] );
	sfm += scaleBC * beff::CSVM_SFb_error[maxbins-1]*2; // add error
	cutsAndSFbinB.push_back( std::pair <float, float> ( 0.679, sfm ) );
	sfm += scaleBC * beff::CSVM_SFb_error[maxbins-1]*charmFactor; // charm additional error
	cutsAndSFbinC.push_back( std::pair <float, float> ( 0.679, sfm ) );
      
	float sfl = beff::CSVL_SFb( beff::ptmax[maxbins-1] );
	sfl += scaleBC * beff::CSVL_SFb_error[maxbins-1]*2; // add error
	cutsAndSFbinB.push_back( std::pair <float, float> ( 0.244, sfl ) );
	sfl += scaleBC * beff::CSVL_SFb_error[maxbins-1]*charmFactor; // charm additional error
	cutsAndSFbinC.push_back( std::pair <float, float> ( 0.244, sfl ) );

	//std::cout << "Lastbin SFs " << sfl << " " << sfm << " " << sft << std::endl;
	cutsAndSFB.push_back(cutsAndSFbinB);
	cutsAndSFC.push_back(cutsAndSFbinC);
  }
   
  m_b = new BinnedBTagShape( binsBC, cutsAndSFB, m_file, "hb", boundX, boundY );
  m_c = new BinnedBTagShape( binsBC, cutsAndSFC, m_file, "hc", boundX, boundY );

  std::vector <EtaPtBin> binsL;
  std::vector< std::vector < std::pair <float, float> > > cutsAndSFL;

  for( size_t j = 0; j < beff::etamin.size(); j++ )
	{
   
	  for( size_t i = 0; i < beff::bins; i++ )
		{
		  EtaPtBin bin( beff::etamin.at(j), beff::etamax.at(j), beff::ptmin.at(i), beff::ptmax.at(i) );
		  binsL.push_back(bin);
		  std::vector < std::pair < float, float > > cutsAndSFbinL;
		  float sft = mistag_CSVT( bin.centerEta(), bin.centerPt(), scaleL );
		  cutsAndSFbinL.push_back( std::pair <float, float> ( 0.898, sft ) );

		  float sfm = mistag_CSVM( bin.centerEta(), bin.centerPt(), scaleL );
		  cutsAndSFbinL.push_back( std::pair <float, float> ( 0.679, sfm ) );

		  float sfl = mistag_CSVL( bin.centerEta(), bin.centerPt(), scaleL );
		  cutsAndSFbinL.push_back( std::pair <float, float> ( 0.244, sfl ) );
   
		  //std::cout << "SFs light " << j << " " << i << " " << sfl << " " << sfm << " " << sft << std::endl;
		  cutsAndSFL.push_back(cutsAndSFbinL);
		}

	  //overflow:
	  {
		std::vector<std::pair< float, float > > cutsAndSFbinL;
		binsL.push_back( EtaPtBin( beff::etamin.at(j), beff::etamax.at(j), beff::ptmax.at(beff::bins-1), 9e99 ) );

		float sft = mistag_CSVT( (beff::etamin.at(j) + beff::etamax.at(j))/2.0, beff::ptmax.at(beff::bins-1), scaleL*2 );
		cutsAndSFbinL.push_back( std::pair <float, float>( 0.898, sft ) );

		float sfm = mistag_CSVM( (beff::etamin.at(j) + beff::etamax.at(j))/2.0, beff::ptmax.at(beff::bins-1), scaleL*2 );
		cutsAndSFbinL.push_back( std::pair <float, float>( 0.679, sfm ) );

		float sfl = mistag_CSVL( (beff::etamin.at(j) + beff::etamax.at(j))/2.0, beff::ptmax.at(beff::bins-1), scaleL*2 );
		cutsAndSFbinL.push_back( std::pair <float, float>( 0.244, sfl ) );

		//std::cout << "SFs light " << sfl << " " << sfm << " " << sft << std::endl;
		cutsAndSFL.push_back(cutsAndSFbinL);
	  }

	}
 
  m_l = new BinnedBTagShape( binsL, cutsAndSFL, m_file, "hl", boundX, boundY );

}

  
float BTagShapeInterface::reshape( float eta, float pt, float csv, int flav)
{
  if(csv < 0) return csv;
  if(csv > 1) return csv; 
  if(flav == 0) return csv;  
  if(fabs(flav) == 5) return  m_b->eval(eta,pt,csv); 
  if(fabs(flav) == 4) return  m_c->eval(eta,pt,csv);
  if(fabs(flav) != 4 && fabs(flav) != 5) return m_l->eval(eta,pt,csv);
  return -10000;    
}


BTagShapeInterface::~BTagShapeInterface()
{
  // delete m_file; 
  // delete m_b;
  // delete m_c;
  // delete m_l;
}

void BTagReshaping()
{
 
}

