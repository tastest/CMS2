#include <iostream>
#include <vector>
#include <set>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEventList.h"
#include "EventListCreator.h"

TEventList* EventListCreator::getAllUniqueEvents( TTree* tree, 
				      const char* name /*= "AllUniqueEvents"*/, 
				      const char* title /*= "All events without duplicates"*/ )
{
  Init(tree);
  TEventList* list = new TEventList(name,title);
  std::set<ULong64_t> known_events;
  unsigned int nEvents = tree->GetEntries();
  for( unsigned int event = 0; event < nEvents; ++event) {
    GetEntry(event);
    ULong64_t key = int((evt_met()-int(evt_met()))*10000) + 
      ULong64_t(evt_run()%100000)*10000 + 
      ULong64_t(evt_event()%100000)*1000000000;
    if ( known_events.insert(key).second ) list->Enter( event );
  }
  return list;
}

