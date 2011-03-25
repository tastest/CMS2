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
    TCut tp_tag1("tp_tag1", "e1_scet>30 && abs(eormu1) == 11 && e1_vbtf80 && iso1 < 0.15 && ((trg_single_e & (1<<0)) || !isdata)");
    TCut tp_tag2("tp_tag2", "e2_scet>30 && abs(eormu2) == 11 && e2_vbtf80 && iso2 < 0.15 && ((trg_single_e & (1<<1)) || !isdata)");

    //
    // probe criteria
    //

    // id
    TCut tp_probe1("tp_probe1", "iso1<0.15 && e1_scet>20 && abs(eormu1) == 11");
    TCut tp_probe2("tp_probe2", "iso2<0.15 && e2_scet>20 && abs(eormu2) == 11");
    
    // trigger
    TCut tp_trg_probe1("tp_trg_probe1", "abs(eormu1) == 11 && e1_vbtf80 + iso1 < 0.15");
    TCut tp_trg_probe2("tp_trg_probe2", "abs(eormu2) == 11 && e2_vbtf80 + iso2 < 0.15");

    TCut tp_base = inclusivez_dilep1 && (tp_tag1 || tp_tag2);
    PreselectBabies(babyVector, tp_base);
    const BabySample *data = babyVector[0]; // fix me

    //
    // Define stuff
    //

    // example for VBTF90 efficiency
    TCut tp_id_1("tp_id_1", "e1_vbtf80");
    TCut tp_id_2("tp_id_2", "e2_vbtf80");

    // trigger
    TCut tp_trg_single1("tp_trg_single1", "(trg_single_e & (1<<0))");
    TCut tp_trg_single2("tp_trg_single2", "(trg_single_e & (1<<1))");
    TCut tp_trg_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL("tp_trg_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL", "(trg_double_e & (1<<0))");


    // no run range restriction
    //TCut tp_event("tp_event", "");
    //TagAndProbe("tp_ele_id_pt", tp_event, "pt1", "pt2", tp_tag1, tp_tag2, tp_probe1, tp_probe2, tp_id_1, tp_id_2, luminosity, 25, 0.0, 100.0, false, babyVector);
    //TagAndProbe("tp_ele_id_eta", tp_event, "eta1", "eta2", tp_tag1, tp_tag2, tp_probe1, tp_probe2, tp_id_1, tp_id_2, luminosity, 30, -3.0, 3.0, false, babyVector);
    //TagAndProbe("tp_ele_trg_single_nvtx", tp_event, "nvtx", "nvtx", tp_tag1, tp_tag2, tp_trg_probe1, tp_trg_probe2, tp_trg_single1, tp_trg_single2, luminosity, 20, -0.5, 19.5, false, babyVector);


    //
    // trigger before and after spike killing
    // before is < 161310
    // after is >= 161310 
    //
    tp_nokill = TCut("tp_nokill", "!isdata || run < 161310");
    tp_kill = TCut("tp_kill", "!isdata || run >= 161310");

    // single trigger

    TCut abs_eta1("abs_eta1", "abs(eta1)");
    TCut abs_eta2("abs_eta2", "abs(eta2)");

    std::cout << "photons" << std::endl;
    //TagAndProbe("tp_HLT_Photon26_IsoVL_Photon18_nvtx", "run>0", "nvtx", "nvtx", tp_tag1, tp_tag2, tp_trg_probe1, tp_trg_probe2, 
    //        "trg_double_e1 & (1<<4)", "trg_double_e2 & (1<<4)", luminosity, 20, -0.5, 19.5, false, babyVector);

    TriggerMonitor("tm_HLT_Photon26_IsoVL_Photon18_nvtx", "nvtx", "mass>76&&mass<106"+((tp_tag1+tp_trg_probe2+"e2_scet>20")||(tp_tag2+tp_trg_probe1+"e1_scet>20")), (tp_tag2+"(trg_double_e1 & (1<<4))")||(tp_tag1+"(trg_double_e2 & (1<<4))") , luminosity, 20, -0.5, 19.5, false, data);



/*
    std::cout << "SINGLE - NO KILL" << std::endl;
    TagAndProbe("tp_single_kill_eta", tp_nokill, abs_eta1, abs_eta2, tp_tag1, tp_tag2, tp_trg_probe1+"e1_scet>30", tp_trg_probe2+"e2_scet>30",
            tp_trg_single1, tp_trg_single2, luminosity, 2, 0.0, 3.0, false, babyVector);

    std::cout << "SINGLE - KILL" << std::endl;
    TagAndProbe("tp_single_kill_eta", tp_kill, abs_eta1, abs_eta2, tp_tag1, tp_tag2, tp_trg_probe1+"e1_scet>30", tp_trg_probe2+"e2_scet>30",
           tp_trg_single1, tp_trg_single2, luminosity, 2, 0.0, 3.0, false, babyVector);
*/

    // double trigger

//    std::cout << "DOUBLE - NO KILL" << std::endl;
//    TagAndProbe("tp_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_nokill_ptgt10_eta", tp_nokill, abs_eta1, abs_eta2, tp_tag1, tp_tag2, tp_trg_probe1+"pt1>10", tp_trg_probe2+"pt2>10", 
//            tp_trg_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, tp_trg_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, luminosity, 2, 0.0, 3.0, false, babyVector);

//    std::cout << "DOUBLE - KILL" << std::endl;
//    TagAndProbe("tp_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_kill_ptgt10_eta", tp_kill, abs_eta1, abs_eta2, tp_tag1, tp_tag2, tp_trg_probe1+"e1_scet>10", tp_trg_probe2+"e2_scet>10",
//            tp_trg_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, tp_trg_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, luminosity, 2, 0.0, 3.0, false, babyVector);

//    TagAndProbe("tp_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_kill_pt", tp_kill, "e1_scet", "e2_scet", tp_tag1, tp_tag2, tp_trg_probe1+"e1_scet>10", tp_trg_probe2+"e2_scet>10",
//            tp_trg_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, tp_trg_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, luminosity, 25, 0.0, 100.0, false, babyVector);
//    TagAndProbe("tp_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_nokill_pt", tp_nokill, "e1_scet", "e2_scet", tp_tag1, tp_tag2, tp_trg_probe1+"e1_scet>10", tp_trg_probe2+"e2_scet>10",
//            tp_trg_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, tp_trg_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, luminosity, 25, 0.0, 100.0, false, babyVector);

}

