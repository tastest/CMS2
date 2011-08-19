// -*- C++ -*-
#ifndef HISTPARM_H
#define HISTPARM_H

#include <string>

// all you need to know to make a histo, really
struct HistParm {
     std::string 	name;
     int 		nbins;
     double 		min;
     double 		max;
};

#endif
