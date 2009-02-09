// -*- C++ -*-
#ifndef SAMPLE_H
#define SAMPLE_H

#include <string>

enum Process { WW, WZ, ZZ, Wjets, DYee, DYmm, DYtt, ttbar, tW, LM1, LM2, LM4,
	       InclusiveMu5Pt50, InclusiveMuPt15, 
	       QCDBCtoEPt20to30, QCDBCtoEPt30to80, QCDBCtoEPt80to170, 
	       QCDEMenrichedPt20to30, QCDEMenrichedPt30to80, QCDEMenrichedPt80to170};

class TChain;
// struct that contains all the necessary information about a sample
class Sample {
public:
     TChain 		*chain;
     enum Process 	process;
     int		histo_color;
     double		kFactor;
     std::string	name;
     bool               sm;
};

// helper functions that provide samples from their default locations
// (takes the guesswork out of data access...)
Sample fWW	();
Sample fWZ	();
// Sample fWZ_incl	();
Sample fZZ	();
Sample fWjets	();
Sample fDYee 	();
Sample fDYmm 	();
Sample fDYtt 	();
Sample fttbar	();
Sample fttbar_taula	();
Sample ftW	();
Sample fSingleTop_tChannel	();
Sample fSingleTop_sChannel	();
Sample fLM1     ();
Sample fLM2     ();
Sample fLM4     ();

// QCD samples
Sample fInclusiveMu5Pt50	();
Sample fInclusiveMuPt15	        ();
Sample fQCDBCtoEPt20to30	();
Sample fQCDBCtoEPt30to80	();
Sample fQCDBCtoEPt80to170	();
Sample fQCDEMenrichedPt20to30	();
Sample fQCDEMenrichedPt30to80	();
Sample fQCDEMenrichedPt80to170  ();


// filter events by process
bool filterByProcess (enum Process sample);

#endif
