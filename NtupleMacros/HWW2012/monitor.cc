#include "monitor.h"
#include "../CORE/CMS2.h"
#include <fstream>

MonitorEventId::MonitorEventId(CMS2& cms2){
  run = cms2.evt_run();
  event = cms2.evt_event();
  lumi = cms2.evt_lumiBlock();
}

MonitorEventId::MonitorEventId(){
  run = 0;
  event = 0;
  lumi = 0;
}

Entry::Entry()
{
  for (unsigned int i=0; i<5; ++i){
    nhyp[i] = 0;
    nevt[i] = 0;
    seen[i] = false;
  }
}

void hypo_monitor::count(CMS2& cms2, HypothesisType type, const char* name, double weight)
{
  std::vector<Entry>::iterator itr = counters.begin();
  while (itr != counters.end() && itr->name != name) itr++;
  Entry* entry(0);
  if ( itr == counters.end() ){
    counters.push_back(Entry());
    entry = &counters.back();
    entry->name = name;
  } else {
    entry = &*itr;
  }
  MonitorEventId id(cms2);
  entry->nhyp[type]++;
  entry->nhyp[ALL]++;
  if (id != entry->lastEvent){
    if (keepEventList) entry->events.push_back(id);
    for (unsigned int i=0; i<5; ++i) entry->seen[i] = false;
    entry->nevt[type]++;
    entry->nevt[ALL]++;
    entry->nevt_weighted[type]+=weight;
    entry->nevt_weighted[ALL]+=weight;
    entry->seen[type] = true;
    entry->lastEvent = id;
//   } else {
//     if ( !entry->seen[type] ){
//       entry->nevt[type]++;
//       entry->nevt[ALL]++;
//       entry->nevt_weighted[type]+=weight;
//       entry->nevt_weighted[ALL]+=weight;
//       entry->seen[type] = true;
//     }
  }
}

void hypo_monitor::print() const
{
  std::cout << "Total number of processed events: \t" << nEvtProcessed << std::endl;
  for ( unsigned int i=0; i<counters.size(); ++i ){
    std::cout << Form("%-40s \thyps: %u/%u/%u/%u/%u \tnevts: %u/%u/%u/%u/%u", counters[i].name.c_str(),
		      counters[i].nhyp[MM],counters[i].nhyp[EE],counters[i].nhyp[EM],counters[i].nhyp[ME],counters[i].nhyp[ALL],
		      counters[i].nevt[MM],counters[i].nevt[EE],counters[i].nevt[EM],counters[i].nevt[ME],counters[i].nevt[ALL]) 	      
	      << std::endl;
    std::ofstream cut_file(Form("cut-%d.txt",i));
    cut_file << Form("%-40s \tnevts: %u/%u/%u/%u/%u", counters[i].name.c_str(),
		     counters[i].nevt[MM],counters[i].nevt[EE],counters[i].nevt[EM],counters[i].nevt[ME],counters[i].nevt[ALL]) << "\n";
    for ( std::vector<MonitorEventId>::const_iterator id=counters[i].events.begin();
	  id!=counters[i].events.end(); ++id ){
      cut_file << id->run << "\t" << id->lumi << "\t" << id->event <<"\n";
    }
    cut_file.close();
  }
  for ( unsigned int i=0; i<counters.size(); ++i ){
    std::cout << Form("%-40s \thyps: %f/%f/%f/%f/%f \tnevts: %f/%f/%f/%f/%f", counters[i].name.c_str(),
		      counters[i].nhyp[MM]/((float)nEvtProcessed),counters[i].nhyp[EE]/((float)nEvtProcessed),counters[i].nhyp[EM]/((float)nEvtProcessed),counters[i].nhyp[ME]/((float)nEvtProcessed),counters[i].nhyp[ALL]/((float)nEvtProcessed),
		      counters[i].nevt[MM]/((float)nEvtProcessed),counters[i].nevt[EE]/((float)nEvtProcessed),counters[i].nevt[EM]/((float)nEvtProcessed),counters[i].nevt[ME]/((float)nEvtProcessed),counters[i].nevt[ALL]/((float)nEvtProcessed)) 
	      << std::endl;
  }
}

void hypo_monitor::makeHistograms(const char* prefix) const
{
  TH1F* hist[4];
  TH1F* histw[4];
  for (unsigned int i=0; i<4; i++){
    hist[i]  = new TH1F(Form("%s_hcuts_%s", prefix, HypothesisTypeName(i)), 
			"Number of events vs cuts", counters.size(), 0, counters.size() );	
    histw[i] = new TH1F(Form("%s_hcuts_weighted_%s", prefix, HypothesisTypeName(i)), 
		       "Number of weighted events vs cuts", counters.size(), 0, counters.size() );	
    for ( unsigned int j=0; j<counters.size(); ++j ){
      hist[i]->SetBinContent(j+1,counters[j].nevt[i]);
      hist[i]->GetXaxis()->SetBinLabel(j+1,counters[j].name.c_str());
      histw[i]->SetBinContent(j+1,counters[j].nevt_weighted[i]);
      histw[i]->GetXaxis()->SetBinLabel(j+1,counters[j].name.c_str());
    }
  }
}
