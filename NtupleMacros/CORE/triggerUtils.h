#ifndef triggerUtils_h
#define triggerUtils_h
bool goodEGTrigger5July2010 (bool is_mc);
int nHLTObjects(const char*, bool quiet = false);
LorentzVector p4HLTObject(const char*, int) ;
void PrintTriggers();
#endif
