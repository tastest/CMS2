#include "TDirectory.h"
#include "TrileptonHist.h"
#include "HistParm.h"
#include "Sample.h"
#include "tools.h"

using std::string;

TrileptonHist::TrileptonHist (const Sample &sample, const string &var_name, 
			    int bins, double min, double max)
{
     build(sample, var_name, bins, min, max);
}

TrileptonHist::TrileptonHist (const Sample &sample, const HistParm &parm)
{
     build(sample, parm.name, parm.nbins, parm.min, parm.max);
}

void TrileptonHist::build (const Sample &sample, const string &var_name, 
			    int bins, double min, double max)
{
//      TDirectory *old_gDirectory = gDirectory;
     gDirectory = histo_directory;
     for (int i = 0; i < 21; ++i) {
          // remove trailing axis label descriptions
       string local_var_name = var_name;
	  if ( var_name.find_first_of(";") < var_name.size() ) {
	    local_var_name = var_name.substr(0,var_name.find_first_of(";"));
	  }
	  string name = sample.name + "_" + local_var_name + "_";
	  name += trilepton_hypo_names[i];
	  // need to make this thing first, without the directory
	  // finding out about it
// 	  bool dir_stat = TH1::AddDirectoryStatus();
	  TH1::AddDirectory(true);
	  histos[i] = new H_t(name.c_str(), var_name.c_str(), bins, min, max);
	  histos[i]->Sumw2();
	  histos[i]->SetFillColor(sample.histo_color);
	  histos[i]->SetLineColor(sample.histo_color);
// 	  TH1::AddDirectory(dir_stat);
     }
//      gDirectory = old_gDirectory;
     memset(entries, 0, sizeof(entries));
     memset(integral, 0, sizeof(integral));
}

TrileptonHist::~TrileptonHist ()
{
//      for (int i = 0; i < 4; delete histos[i++]);
}

int TrileptonHist::Fill (enum TrileptonHypType i, double x, double w)
{
     histos[i]->Fill(x, w);
     entries[i]++;
     integral[i] += w;
     int ret = histos[TRILEPTON_ALL]->Fill(x, w);
     entries[TRILEPTON_ALL]++;
     integral[TRILEPTON_ALL] += w;
     return ret;
}
