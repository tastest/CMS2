// -*- C++ -*-

#ifndef TOOLS_H
#define TOOLS_H

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TVector3.h"
#include <algorithm>
#include <vector>
#include <set>
#include "Math/VectorUtil.h"

class CMS2;
using std::vector;

unsigned int encodeTriLeptonCand(unsigned int bucket,unsigned int first, unsigned int second, unsigned int third);
unsigned int decodeBucket(unsigned int cand);
unsigned int decodeFirst(unsigned int cand);
unsigned int decodeSecond(unsigned int cand);
unsigned int decodeThird(unsigned int cand);
TH1F* book1DHist(const char* name, const char* title, unsigned int nbins, float low, float high, const char* xtitle, const char* ytitle, int color = 1);
TH1F* book1DVarHist(const char* name, const char* title, unsigned int nbins, float* bins, const char* xtitle, const char* ytitle, int color = 1);
TH1F* book1DVarHist(const char* name, const char* title, vector<float> &bins, const char* xtitle, const char* ytitle, int color = 1);
TH2F* book2DHist(const char* name, const char* title, unsigned int nxbins, float xlow, float xhigh, unsigned int nybins, float ylow, float yhigh, const char* xtitle, const char* ytitle, const char* ztitle, int color = 1);
TH2F* book2DVarHist(const char* name, const char* title, unsigned int nxbins, float* xbins, unsigned int nybins, float* ybins, const char* xtitle, const char* ytitle, const char* ztitle, int color = 1);
TH2F* book2DVarHist(const char* name, const char* title, vector<float> &xbins, vector<float> &ybins, const char* xtitle, const char* ytitle, const char* ztitle, int color = 1);
TH3F* book3DHist(const char* name, const char* title, unsigned int nxbins, float xlow, float xhigh, unsigned int nybins, float ylow, float yhigh, unsigned int nzbins, float zlow, float zhigh, const char* xtitle, const char* ytitle, const char* ztitle, int color = 1);
TH3F* book3DVarHist(const char* name, const char* title, unsigned int nxbins, float* xbins, unsigned int nybins, float* ybins, unsigned int nzbins, float* zbins, const char* xtitle, const char* ytitle, const char* ztitle, int color = 1);
TH3F* book3DVarHist(const char* name, const char* title, vector<float> &xbins, vector<float> &zbins, vector<float> &ybins, const char* xtitle, const char* ytitle, const char* ztitle, int color = 1);
float mee(int i, int j);
float mmm(int i, int j);
bool goodLeptonIsolated(int bucket, int first, int second, int third);
float ptLowestPtLepton(int bucket, int first, int second, int third);
bool passTriggerLeptonMinPtCut(int bucket, int first, int second, int third, float triggerLeptonMinPtCut);
TString printCand(int bucket, int first, int second, int third);

namespace duplicate_removal {
     struct DorkyEventIdentifier {
	  // this is a workaround for not having unique event id's in MC 
	  DorkyEventIdentifier (CMS2 &cms2);
	  DorkyEventIdentifier (unsigned long int r, unsigned long int e, unsigned long int l);
	  unsigned long int run, event, lumi_section;
	  float trks_d0;
	  float trks_pt, trks_eta, trks_phi;
	  bool operator < (const DorkyEventIdentifier &) const;
	  bool operator == (const DorkyEventIdentifier &) const;
     };
     extern std::set<DorkyEventIdentifier> already_seen;
     bool is_duplicate (const DorkyEventIdentifier &id);
     
     struct Status3Identifier {
	  Status3Identifier (CMS2 &cms2);
	  unsigned long int run, event, lumi_section;
	  std::multiset<float> pts;
	  bool operator < (const Status3Identifier &) const;
	  bool operator > (const Status3Identifier &) const;
	  bool operator == (const Status3Identifier &) const;
     };
     extern std::set<Status3Identifier> already_seen_stat3;
     bool is_duplicate (const Status3Identifier &id);
     bool operator == (const Status3Identifier &, const Status3Identifier &);
}

void saveHist(const char* filename, const char* pat="*");

double correctd0Phi(int trk_idx);

bool goodrun (unsigned int run, unsigned int lumi_block);

extern class TDirectory *histo_directory;
#endif
