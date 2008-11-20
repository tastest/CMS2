#include "TDirectory.h"
#include "DileptonHist.h"
#include "HistParm.h"
#include "Sample.h"
#include "tools.h"

using std::string;

DileptonHist::DileptonHist (const Sample &sample, const string &var_name, 
			    int bins, double min, double max)
{
     build(sample, var_name, bins, min, max);
}

DileptonHist::DileptonHist (const Sample &sample, const HistParm &parm)
{
     build(sample, parm.name, parm.nbins, parm.min, parm.max);
}

void DileptonHist::build (const Sample &sample, const string &var_name, 
			    int bins, double min, double max)
{
     TDirectory *old_gDirectory = gDirectory;
     gDirectory = histo_directory;
     for (int i = 0; i < 4; ++i) {
	  string name = sample.name + "_" + var_name + "_";
	  name += dilepton_hypo_names[i];
	  // need to make this thing first, without the directory
	  // finding out about it
	  bool dir_stat = TH1::AddDirectoryStatus();
	  TH1::AddDirectory(true);
	  histos[i] = new H_t(name.c_str(), var_name.c_str(), bins, min, max);
	  histos[i]->Sumw2();
	  histos[i]->SetFillColor(sample.histo_color);
	  histos[i]->SetLineColor(sample.histo_color);
	  TH1::AddDirectory(dir_stat);
     }
     gDirectory = old_gDirectory;
     memset(entries, 0, sizeof(entries));
     memset(integral, 0, sizeof(integral));
}

DileptonHist::~DileptonHist ()
{
//      for (int i = 0; i < 4; delete histos[i++]);
}

int DileptonHist::Fill (enum DileptonHypType i, double x, double w)
{
     histos[i]->Fill(x, w);
     entries[i]++;
     integral[i] += w;
     int ret = histos[DILEPTON_ALL]->Fill(x, w);
     entries[DILEPTON_ALL]++;
     integral[DILEPTON_ALL] += w;
     return ret;
}
