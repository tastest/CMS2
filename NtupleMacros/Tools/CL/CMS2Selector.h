// -*- C++ -*-

// $Id: CMS2Selector.h,v 1.1 2010/04/16 19:01:22 jmuelmen Exp $

#ifndef CMS2SELECTOR_H
#define CMS2SELECTOR_H

#include "TSelectorDraw.h"

class CMS2Selector : public TSelectorDraw {
public:
     CMS2Selector ();
     virtual void Begin (TTree*);
     virtual void Init (TTree*);
     virtual Bool_t Notify ();
     virtual void ProcessFill (Long64_t);
     virtual void ProcessFillMultiple (Long64_t);
     virtual void ProcessFillObject (Long64_t);
};

#endif
