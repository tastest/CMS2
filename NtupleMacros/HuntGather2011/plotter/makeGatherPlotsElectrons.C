#include "cuts.h"

void makeGatherPlotsElectrons(const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making Electron Monitor...\n";

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    // tag criteria
    // generally, do not change
    TCut tp_tag1("tp_tag1", "pt1>20 && abs(eormu1) == 11 && e1_vbtf90full && (trg_single_e & (1<<0) || !isdata)");
    TCut tp_tag2("tp_tag2", "pt2>20 && abs(eormu2) == 11 && e2_vbtf90full && (trg_single_e & (1<<1) || !isdata)");

    // probe criteria
    // generally, do not change
    TCut tp_probe1("tp_probe1", "iso1<0.15 && pt1>20 && abs(eormu1) == 11");
    TCut tp_probe2("tp_probe2", "iso2<0.15 && pt2>20 && abs(eormu2) == 11");

    TCut tp_base = inclusivez_dilep1 && (tp_tag1 || tp_tag2 || TCut("notdata", "!isdata"));
    PreselectBabies(babyVector, tp_base);

    //
    // Define stuff
    //

    // cuts
    TCut tp_id_1("tp_id_1", "e1_vbtf90");
    TCut tp_id_2("tp_id_2", "e2_vbtf90");
    TagAndProbe("tp_ele_id_pt", "pt1", "pt2", tp_tag1, tp_tag2, tp_probe1, tp_probe2, tp_id_1, tp_id_2, luminosity, 25, 0.0, 100.0, false, babyVector);

}

