#ifndef WW_monitor_h
#define WW_monitor_h
#include <vector>
#include <string>
#include "wwtypes.h"

class CMS2;
struct MonitorEventId { 
  unsigned long int run, event, lumi; 
  // --------------------------------------------------------------- //
  MonitorEventId(CMS2&);
  MonitorEventId();
  bool operator < (const MonitorEventId& id) const{
    if (run != id.run) return run < id.run;
    if (lumi != id.lumi) return lumi < id.lumi;
    return event < id.event;
  }
  bool operator == (const MonitorEventId& id) const{
    return (run==id.run) && (lumi==id.lumi) && (event==id.event); 
  }
  bool operator != (const MonitorEventId& id) const{
    return ! operator == (id);
  }
};

struct Entry {
  unsigned int nhyp[4];
  unsigned int nevt[4];
  double nhyp_weighted[4];
  double nevt_weighted[4];
  bool seen[4];
  MonitorEventId lastEvent;
  std::string name;
  // -------------------------------------------------------------- //
  Entry();
};

struct hypo_monitor{
  void count(CMS2&, HypothesisType type, const char* name, double weight=1.0);
  void print() const;
  void makeHistograms(const char* prefix) const;
  hypo_monitor():nEvtProcessed(0){}
  // -------------------------------------- //
  std::vector<Entry> counters;
  unsigned int nEvtProcessed;
};
#endif
