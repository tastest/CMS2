//---------------------------------------------//
// Custom class to read pid, |P|, and costheta //
// derived from Ztautau MC and stored in text  //
// file for application to electrons & muons   //
// in data ( data driven Ztautau prediction )  //
//---------------------------------------------//

// all includes in header
#include "tauify.h"  

using namespace std;

// private Get Methods - access by index
unsigned int Tauify::TauSize(void){
  return tau_data.size();
}
int Tauify::ParticleId(int index){
  return tau_data[index].first;
}
double Tauify::MomentumCM(int index){
  return tau_data[index].second.first;
}
double Tauify::CosThetaCM(int index){
  return tau_data[index].second.second;
}

// public Get Methods
LorentzVector Tauify::TauP4(void){
  return p4_tau_lepton_lab;
}
int Tauify::ParticleId(void){
  return id;
}
float Tauify::MomentumCM(void){
  return p_cm;
}
float Tauify::CosThetaCM(void){
  return costheta_cm;

}
float Tauify::TauMET(void){
  return tau_met;
}
float Tauify::TauIso(void){
  return tau_iso;
}
float Tauify::TauIP(void){
  return tau_ip;
}

// Set Methods
void Tauify::SetLepton( LorentzVector p4_arg, float met_arg, float iso_arg, float ip_arg){

  // set the 4 vector of the lepton in the LAB
  p4_lepton_lab = p4_arg;

  // how to handle lepton id? random for now...
  // assume all tau decays are equally likely, pick a random P and cos(theta) from the file
  TRandom3 *rand1 = new TRandom3(0);
  int index       = rand1->Integer( TauSize() );
  id              = ParticleId(index);
  p_cm            = MomentumCM(index);
  costheta_cm     = CosThetaCM(index);
  double theta    = acos(costheta_cm); 

  // 3 vectors
  Vector p3_lepton_lab;
  Vector para;
  Vector perp;

  // 3 momentum of lepton in the LAB
  p3_lepton_lab.SetCoordinates( (double)p4_lepton_lab.Px(), (double)p4_lepton_lab.Py(), (double)p4_lepton_lab.Pz() );
   
  // 3 momentum paralell to lepton
  para = p3_lepton_lab;
  para /= sqrt( p3_lepton_lab.Mag2() );
  para *= p_cm*costheta_cm;

  // the 3 momentum perpendicular to the lepton is free to lie anywhere in the normal plane
  // taking the cross product of the paralell vector and any other vector will give a vector that lies in that plane
  // easiest to choose one of the axis unit vectors
  // but, make sure the unit vector isn't equal or close to the normal vector
  // to do this, use the unit vector corresponding to the smallest normal component

  Vector p3_cross;
  double smallest = min( para.X(), min( para.Y(), para.Z() ) );
  if( smallest == para.X() ){
    p3_cross.SetCoordinates( 1.0, 0.0, 0.0 );
  } 
  else if( smallest == para.Y() ){
    p3_cross.SetCoordinates( 0.0, 1.0, 0.0 );
  } 
  else if( smallest == para.Z() ){
    p3_cross.SetCoordinates( 0.0, 0.0, 1.0 );
  } 
  else {
    cout << "ERROR! Couldn't find smallest component... exiting." << endl; 
    exit(1);
  } 

  // 3 momentum perpendicular to lepton
  perp = para.Cross( p3_cross );
  perp /= sqrt( perp.Mag2() );
  perp *= p_cm*sin(theta);
  
  // random rotation
  TRandom3 *rand2 = new TRandom3(0);
  double phi      = rand2->Uniform( 0, 2*TMath::Pi() );
  ROOT::Math::AxisAngle rot( para, phi );
  perp = rot( perp );

  //
  Vector p3_tau_lepton_cm = para + perp;




  /* sanity checks */
  double epsilon = 1.0/1000000000000000.0;

  // check perp and para are orthogonal
  if( fabs(perp.Dot(para)) > epsilon ){ 
    cout << "ERROR! perp & para not orthogonal para.perp = " << perp.Dot(para) << endl;
    exit(1);
  }
  
  // check para magnitude
  double Mag_para = sqrt( para.Mag2() );
  if( fabs(Mag_para - fabs(p_cm*costheta_cm)) > epsilon ){ 
    cout << "ERROR! paralell magnitude is wrong: " << Mag_para << " != " << fabs(p_cm*costheta_cm) << endl;
    exit(1);
  }

  // check perp magnitude
  double Mag_perp = sqrt( perp.Mag2() );
  if( fabs(Mag_perp - fabs(p_cm*sin(theta))) > epsilon ){ 
    cout << "ERROR! perpendicular magnitude is wrong: " << Mag_perp << " != " << fabs(p_cm*sin(theta)) << endl;
    exit(1);
  }

  // check p^2 = ( para + perp )^2
  double test = sqrt( perp.Dot(perp) + para.Dot(para) + 2*para.Dot(perp) );
  if( fabs(test - p_cm) > epsilon ){ 
    cout << "ERROR! p^2 != ( para + perp)^2: " << p_cm << " != " << test << endl;
    exit(1);
  }



  /* boost back to the lab */

  // get the lepton boost from the LAB to the CM
  ROOT::Math::Boost boost_cm( p4_lepton_lab.BoostToCM().x(), p4_lepton_lab.BoostToCM().y(), p4_lepton_lab.BoostToCM().z() );

  // get the lepton boost from the CM to the LAB
  ROOT::Math::Boost boost_lab = boost_cm.Inverse();

  // do something smarter about the tau mass later
  p4_tau_lepton_cm.SetPx( p3_tau_lepton_cm.X() );
  p4_tau_lepton_cm.SetPy( p3_tau_lepton_cm.Y() );
  p4_tau_lepton_cm.SetPz( p3_tau_lepton_cm.Z() );
  p4_tau_lepton_cm.SetE( 1.7 );

  // boost from CM to LAB
  p4_tau_lepton_lab = boost_lab*p4_tau_lepton_cm;

  //
  tau_met = met_arg;
  tau_iso = iso_arg;
  tau_ip  = ip_arg;

  tau_met = -999;
  tau_iso = -999;
  tau_ip  = -999;

  return;
}

// constructor
Tauify::Tauify( const char* infile, bool verbose /* default false in header */ ){
 
  /* "verbose = true" useful to check that input file is read correctly */

  //---------------------------------------------------------------------------------//
  // parse (3) space delimited columns from text file into scalar types (int, float) //
  // functionality is probably generally useful and could be put somwhere common...  //
  // also probably already available as some kind of library                         //
  //---------------------------------------------------------------------------------//
   
  // parse input text file
  ifstream input_stream(infile);
  if( !input_stream.is_open() ){  // check file exists
    cout << endl << "Error! Could not open file: " << infile << " ... exiting." << endl << endl;
    exit(1);
  }
  else {  // read file
    int index = 0;
    while (! input_stream.eof() ){
      // read line
      string line;
      if( getline (input_stream, line) ){
        // temp variables to read text into int, float 
        // initialize large and negative ( -999 may be a valid pid )
        int   particle_id = numeric_limits<int>::min();
        float momentum    = numeric_limits<float>::min();
        float costheta    = numeric_limits<float>::min();
        // split line delmited by N white spaces
        int iColumn = 0;
        //cout << line << endl;
        // initial tokenization
        char *pch = strtok( (char*)line.c_str(), " ");
        while( pch != NULL ){
          // string to int, float conversion
          iColumn++;
          if( iColumn == 1 ){
            particle_id = atoi( pch );
          }
          if( iColumn == 2 ){
            momentum = atof( pch );
          }
          if( iColumn == 3 ){
            costheta = atof( pch );
          }
          if( iColumn > 3 ){ 
            cout << endl << "Error! Only 3 columns expected... exiting." << endl << endl;
            exit(1);
          }
          //cout << "\t:" << pch << ":" << endl;
          // tokenize
          pch = strtok( NULL, " ");
        }

        //  
        if(verbose) cout << particle_id << " " << momentum << " " << costheta << endl;

        // pair of floats
        pair<float, float> pair_ff;
        pair_ff.first = momentum;
        pair_ff.second = costheta;

        // pair of an int and a pair of floats
        pair<int, pair<float, float> > pair_iff;
        pair_iff.first = particle_id;
        pair_iff.second = pair_ff;
       
        tau_data[index] = pair_iff;

        index++;

      }
    }
    input_stream.close();
  }

  return;
}
