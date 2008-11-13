// -*- C++ -*-
#ifndef SAMPLE_H
#define SAMPLE_H

#include <string>

enum Process { WW, WZ, ZZ, Wjets, DYee, DYmm, DYtt, ttbar, tW, LM1, LM2, LM4 };

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
Sample fZZ	();
Sample fWjets	();
Sample fDYee 	();
Sample fDYmm 	();
Sample fDYtt 	();
Sample fttbar	();
Sample ftW	();
Sample fLM1     ();
Sample fLM2     ();
Sample fLM4     ();

// filter events by process
bool filterByProcess (enum Process sample);

#endif
