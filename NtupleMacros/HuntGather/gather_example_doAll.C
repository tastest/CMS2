#include "cuts.h"
#include "gather.h"
#include "goodrun.h"

#include "TChain.h"
#include "TH1F.h"

#include <iostream>

void gather_example_doAll()
{
    // this sets the json file obviously; it's
    // easiest to just  soft link to  the real
    // file as json.txt
    std::cout << "Using json.txt for goodruns\n";
    set_goodrun_file_json("json.txt");

    // calculate  integrated luminosity  in /fb
    // note that I input /nb and so get out /nb
    // then scale to /fb which  is what  I need
    // 1e-3*GetIntLumi(/pb)  is the  same thing
    // doesn't matter so long as it corresponds
    // to the lumi of the json file and that it
    // is correctly scaled to /fb afterward
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
    // here is a breakdown of the Plot method
    // TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel, float intlumifb, float kfactor,
    //                   unsigned int nbins, float xlo, float xhi);
    // the arguments are:
    //     prefix for histogram names
    //     chain of files
    //     the baby-ntuple-field you want to draw
    //     a selection, generally one used by the hunt defined in cuts.h
    //     a pre-selection, this is applied in addition to the selection, it is meant for
    //         things like stitching together MC, i.e. "pthat<30"; I imagine it can  also
    //         be used to add OS,SS,OF,SF cuts, i.e. "abs(eormu1)!=abs(eormu2) w/o having
    //         to create a new TCut in cuts.h
    //     the integrated luminosity in /fb; scale1fb is of course multiplied by this no.
    //     whatever additional scale factor you need
    //     then the standard nbins and histogram range
    TH1F* hmc   = Plot("zjets",czjets_dilep,"mass",inclusivez_dilep,"",f_intlumifb,1.27,100,60,120,false,false);

    // TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel,
    //                   unsigned int nbins, float xlo, float xhi);
    // same as above but as this is data there is no need to specify intlumifb or kfactor
    TH1F* hdata = Plot("data", cdata_dilep, "mass",inclusivez_dilep,"",100,60,120,false,true);

    // there is another Plot implementation for data as well:
    // TH1F* Plot(const char *pfx, TChain *chain, const char *field, TCut sel, TCut presel, float kfactor,
    //                   unsigned int nbins, float xlo, float xhi);
    // which is the same as above but allows specification of a kfactor, which I can imagine using

    float ymax = hdata->GetMaximum() > hmc->GetMaximum() ? hdata->GetMaximum() : hmc->GetMaximum();
    hdata->SetMarkerStyle(20);
    hmc->SetFillColor(kRed);
    hmc->SetLineColor(kRed);

    hmc->SetMaximum(ymax*1.1);
    hmc->Draw();
    hdata->Draw("pesames");
}
