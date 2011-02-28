#include "cuts.h"

void makeGatherPlotsElectrons(const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making Electron Monitor...\n";

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, base_dilep);

    //
    // Define stuff
    //

    Int_t run_start = 146000;
    Int_t run_end = 150000;

    // ok, so it's actually 4 runs per bin
    // but configure as needed
    Int_t n_run_bins = (run_end - run_start)/4;

    // event selections
    TCut cut_trg_zee("trg_zee", inclusivez_dilep+ee_dilep);
    TCut cut_trg_zmm("trg_zmm", inclusivez_dilep+mm_dilep);

    // the single object triggers
    // note that trg_single_x is a mask
    // bit zero means LT passed and bit one means LL passed
    TCut cut_trg_single_e("trig_single_e", "(trg_single_e & (1<<0)) || (trg_single_e & (1<<1))");
    TCut cut_trg_single_mu("trig_single_mu", "(trg_single_mu & (1<<0)) || (trg_single_mu & (1<<1))");

    //
    // Make the plots
    //

    // this one is complicated, but it basically says that either the first LT or LL must pass tight criteria
    // where that means single trigger + offline ID
    TCut cut_tp_tag("cut_tp_tag", "(abs(eormu1) == 11 && e1_vbtf90full && (trg_single_e & (1<<0)) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic && (trg_single_e & (1<<0)))) || (abs(eormu2) == 11 && e2_vbtf90full && (trg_single_e & (1<<1)) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic && (trg_single_mu & (1<<1))))");
    TCut cut_tp_base("cut_tp_base", base_dilep1+sf_dilep+inclusivez_dilep1+cut_trg_single_e+cut_tp_tag);

    TCut var_probe_sigmaIEtaIEta("probe_sigmaIEtaIEta", "TMath::Max((e1_sigieie * (trg_single_e & (1<<1))/2.), (e2_sigieie * (trg_single_e & (1<<0))/1.))");

    DrawAll(var_probe_sigmaIEtaIEta, "electron_sigmaietaieta", cut_tp_base+cut_trg_single_e, luminosity, 100, 0, 0.05, 0, babyVector);
    DrawAll("mass", "electron_mass", cut_tp_base+cut_trg_single_e, luminosity, 100, 0, 200, 0, babyVector);

}

