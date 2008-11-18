#ifndef UTILITIES_H
#define UTILITIES_H

#include "TH1F.h"
#include "TVector3.h"
#include <algorithm>
#include <set>
#include "../CORE/selections.h"
#include "Math/VectorUtil.h"

unsigned int encodeTriLeptonCand(unsigned int bucket,unsigned int first, unsigned int second, unsigned int third);
unsigned int decodeBucket(unsigned int cand);
unsigned int decodeFirst(unsigned int cand);
unsigned int decodeSecond(unsigned int cand);
unsigned int decodeThird(unsigned int cand);
TH1F* book1DHist(const char* name, const char* title, unsigned int nbins, float low, float high, const char* xtitle, const char* ytitle);
TH1F* book1DVarHist(const char* name, const char* title, unsigned int nbins, float* bins, const char* xtitle, const char* ytitle);
float mee(int i, int j);
float mmm(int i, int j);
bool goodLeptonIsolated(int bucket, int first, int second, int third);
float ptLowestPtLepton(int bucket, int first, int second, int third);
bool passTriggerLeptonMinPtCut(int bucket, int first, int second, int third, float triggerLeptonMinPtCut);
TString printCand(int bucket, int first, int second, int third);

struct DorkyEventIdentifier {
     // this is a workaround for not having unique event id's in MC 
     unsigned long int run, event;
     float trks_d0;
     float hyp_lt_pt, hyp_lt_eta, hyp_lt_phi;
     bool operator < (const DorkyEventIdentifier &) const;
     bool operator == (const DorkyEventIdentifier &) const;
};
extern std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id);

void saveHist(const char* filename, const char* pat="*");

extern class TDirectory *histo_directory;
#endif
