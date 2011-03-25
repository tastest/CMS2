#include "cuts.h"

void makeGatherPlotsMuons(const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making Muon Monitor...\n";

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    // tag criteria
    // generally, do not change
    TCut tp_tag1("tp_tag1", "pt1>20 && abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic && ((trg_single_mu & (1<<0)) || !isdata)");
    TCut tp_tag2("tp_tag2", "pt2>20 && abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic && ((trg_single_mu & (1<<1)) || !isdata)");

    //
    // probe criteria
    //

    // for id studies
    TCut tp_probe1("tp_probe1", "iso1<0.15 && pt1>20 && abs(eormu1) == 13");
    TCut tp_probe2("tp_probe2", "iso2<0.15 && pt2>20 && abs(eormu2) == 13");

    // for trigger studies
    TCut tp_trg_probe1("tp_trg_probe1", "pt1>10 && abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic");
    TCut tp_trg_probe2("tp_trg_probe2", "pt2>10 && abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic");

    TCut tp_base = inclusivez_dilep1 && (tp_tag1 || tp_tag2);
    PreselectBabies(babyVector, tp_base);

    //
    // Define stuff
    //

    // example for TTbarV2 ID efficiency
    TCut tp_id_1("tp_id_1", "mu1_muonid && ! mu1_cosmic");
    TCut tp_id_2("tp_id_2", "mu2_muonid && ! mu2_cosmic");

    // double trigger
    TCut cut_trg_DoubleMu7("trg_DoubleMu7", "(trg_double_mu & (1<<0)) == (1<<0)");
    
    // single trigger
    TCut tp_trg_single1("tp_trg_single1", "(trg_single_mu & (1<<0))");
    TCut tp_trg_single2("tp_trg_single2", "(trg_single_mu & (1<<1))");


    // no run range restriction
    TCut tp_event("tp_event", "");
//    TagAndProbe("tp_mu_trg_eta", tp_event, "eta1", "eta2", tp_tag1, tp_tag2, tp_trg_probe1, tp_trg_probe2, cut_trg_DoubleMu7, cut_trg_DoubleMu7, luminosity, 15, -3.0, 3.0, false, babyVector);
    TagAndProbe("tp_mu_trg_single_eta", tp_event, "eta1", "eta2", tp_tag1, tp_tag2, tp_trg_probe1+"pt1>20", tp_trg_probe2+"pt2>20", tp_trg_single1, tp_trg_single2, luminosity, 30, -3.0, 3.0, false, babyVector);
    TagAndProbe("tp_mu_trg_single_nvtx", tp_event, "nvtx", "nvtx", tp_tag1, tp_tag2, tp_trg_probe1+"pt1>20", tp_trg_probe2+"pt2>20", tp_trg_single1, tp_trg_single2, luminosity, 20, -0.5, 19.5, false, babyVector);



    // run range restriction 
    //tp_event = TCut("tp_event", "!isdata || run > 149000");
    //TagAndProbe("tp_mu149000_id_pt", tp_event, "pt1", "pt2", tp_tag1, tp_tag2, tp_probe1, tp_probe2, tp_id_1, tp_id_2, luminosity, 25, 0.0, 100.0, false, babyVector);


}

