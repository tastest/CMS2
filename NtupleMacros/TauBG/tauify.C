//---------------------------------------------//
// Custom class to read pid, |P|, and costheta //
// derived from Ztautau MC and stored in text  //
// file for application to electrons & muons   //
// in data ( data driven Ztautau prediction )  //
//---------------------------------------------//

// header
#include "tauify.h"  

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <limits>

using namespace std;

// Get Methods
LorentzVector TauP4(void){
  LorentzVector p4_return;
  p4_return.SetPx(-999.0);
  p4_return.SetPy(-999.0);
  p4_return.SetPz(-999.0);
  p4_return.SetE(-999.0);
  return p4_return;
}
float TauMET(void){
  return -999.0;
}
float TauIso(void){
  return -999.0;
}
float TauIP(void){
  return -999.0;
}

// Set Methods
void Tauify::SetLepton( LorentzVector p4_lepton, float met, float iso, float d0){
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

      }
    }
    input_stream.close();
  }


  return;
}
