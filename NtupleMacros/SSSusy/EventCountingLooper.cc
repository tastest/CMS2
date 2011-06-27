#include "Looper.h"

EventCountingLooper::EventCountingLooper (Sample s, cuts_t c, const char *fname) 
     : Looper(s, c, fname)
{
     // zero out the event counters (don't comment this out)
     memset(events_passing_, 0, sizeof(events_passing_));
     memset(events_passing_w2_, 0, sizeof(events_passing_w2_));
}

void EventCountingLooper::End ()
{
     // this isn't elegant, but it's fast to write:
     // replace the candidate counts with event counts 
     memcpy(cands_passing_, events_passing_, sizeof(events_passing_));
     memcpy(cands_passing_w2_, events_passing_w2_, sizeof(events_passing_w2_));
     // then print the default tables and stuff (which just prints out
     // the cand statistics)
     Looper::End();
}

void EventCountingLooper::BeforeDilepHistos ()
{
     memcpy(cands_passing_prev_, cands_passing_, sizeof(cands_passing_));
}

void EventCountingLooper::AfterDilepHistos ()
{
     const double weight = Weight(0);
     for (int i = 0; i < 4; ++i) {
	  if (cands_passing_[i] != cands_passing_prev_[i]) {
	       events_passing_[i] += weight;
	       events_passing_w2_[i] += weight * weight;
	  }
     }
}
