// -*- C++ -*-

// $Id: CMS2TreePlayer.h,v 1.1 2010/04/16 19:01:25 jmuelmen Exp $

#ifndef CMS2TREEPLAYER_H
#define CMS2TREEPLAYER_H

#include "TTreePlayer.h"

class CMS2TreePlayer : public TTreePlayer {
private:
     // no copy ctor or assignment operator
     CMS2TreePlayer (const CMS2TreePlayer &) { }
     CMS2TreePlayer &operator = (const CMS2TreePlayer &);

public:
     CMS2TreePlayer ();

public:
     ClassDef(CMS2TreePlayer, 1)
};

#endif
