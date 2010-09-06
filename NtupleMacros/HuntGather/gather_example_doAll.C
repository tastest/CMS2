#include "cuts.h"
#include "gather.h"
#include "goodrun.h"

#include "TChain.h"
#include "TH1F.h"

#include <iostream>

void gather_example_doAll()
{
    std::cout << "Using json.txt for goodruns\n";
    set_goodrun_file_json("json.txt");

    // calculate integrated luminosity in /fb
    // note that I input /nb and so get out /nb
    // then scale to /fb which is what I need
    // 1e-3*GetIntLumi(/pb) is the same thing
    // it matters not so long as it corresponds
    // to the lumi of the json file
    float f_intlumifb = 1e-6*GetIntLumi(2790);
    std::cout << "Integrated luminosity: " << f_intlumifb << "/fb\n";

    // data babies
    TChain *cdata_emu = new TChain("tree");
    cdata_emu->Add("/tas05/disk00/jribnik/hunt/emu_baby/*.root");
    TChain *cdata_dilep = new TChain("tree");
    cdata_dilep->Add("/tas05/disk00/jribnik/hunt/dilep_baby/*.root");
    TChain *cdata_trilep = new TChain("tree");
    cdata_trilep->Add("/tas05/disk00/jribnik/hunt/trilep_baby/*.root");

    // mc babies
    TChain *czjets_emu = new TChain("tree");
    czjets_emu->Add("/tas05/disk00/jribnik/huntmc/ZJets-madgraph_Spring10-START3X_V26_S09-v1/emu_baby/*.root");
    TChain *czjets_dilep = new TChain("tree");
    czjets_dilep->Add("/tas05/disk00/jribnik/huntmc/ZJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root");
    TChain *czjets_trilep = new TChain("tree");
    czjets_trilep->Add("/tas05/disk00/jribnik/huntmc/ZJets-madgraph_Spring10-START3X_V26_S09-v1/trilep_baby/*.root");

    // make plots
    TH1F* hdata = Plot("data", cdata_dilep, "mass",inclusivez_dilep,"",100,60,120);
    TH1F* hmc   = Plot("zjets",czjets_dilep,"mass",inclusivez_dilep,"",f_intlumifb,1.27,100,60,120);

    float ymax = hdata->GetMaximum() > hmc->GetMaximum() ? hdata->GetMaximum() : hmc->GetMaximum();
    hdata->SetMarkerStyle(20);
    hmc->SetFillColor(kRed);
    hmc->SetLineColor(kRed);

    hmc->SetMaximum(ymax*1.1);
    hmc->Draw();
    hdata->Draw("pesames");
}
