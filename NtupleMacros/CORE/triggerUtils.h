#ifndef triggerUtils_h
#define triggerUtils_h

#include "CMS2.h"

bool goodEGTrigger5July2010 (bool);
int nHLTObjects(const char*);
LorentzVector p4HLTObject(const char*, int) ;
void PrintTriggers();
#endif
