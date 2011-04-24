#include "cuts.h"

void makeGatherPlotsZMet(TString prefix, const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making Z+Met plots...\n";

    //
    // Define specific cuts for ZMet plots
    //

    TCut osz_pt2010("(pt1 > 20 && pt2 > 20)");
    TCut osz_dilep2("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
    TCut osz_dilep3("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
    TCut osz_dilep4("ntrks > 2");
    TCut osz_dilep5("eormu1*eormu2<0");
    TCut osz_dilep_all("base_dilep",osz_pt2010+osz_dilep2+osz_dilep3+osz_dilep4+osz_dilep5);
    TCut osz_zmass("(hyp_type == 0 || hyp_type == 3) && mass > 81. && mass < 101.");

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, osz_dilep_all);

    //
    // Make the plots
    //

    //non-integral plots
    DrawAll("mass",     prefix+"_osz_mass_ee_2jets",             TCut("ee {NJets>=2}", osz_dilep_all + "njets>=2 && hyp_type==3")             , luminosity, 40, 0., 200, 0, babyVector);
    DrawAll("mass",     prefix+"_osz_mass_ee_2jets_pfmet60",     TCut("ee {NJets>=2, PFMet>60}", osz_dilep_all + "njets>=2 && hyp_type==3 && pfmet>60") , luminosity, 40, 0., 200, 0, babyVector);
    DrawAll("mass",     prefix+"_osz_mass_mm_2jets",             TCut("mm {NJets>=2}", osz_dilep_all + "njets>=2 && hyp_type==0")             , luminosity, 40, 0., 200, 0, babyVector);
    DrawAll("mass",     prefix+"_osz_mass_mm_2jets_pfmet60",     TCut("mm {NJets>=2, PFMet>60}", osz_dilep_all + "njets>=2 && hyp_type==0 && pfmet>60") , luminosity, 40, 0., 200, 0, babyVector);
    DrawAll("pfmet",    prefix+"_osz_pfmet_ee_2jets",            TCut("ee {NJets>=2}", osz_dilep_all + osz_zmass + "njets>=2 && hyp_type==3") , luminosity, 40, 0., 200, 0, babyVector);
    DrawAll("pfmet",    prefix+"_osz_pfmet_mm_2jets",            TCut("mm {NJets>=2}", osz_dilep_all + osz_zmass + "njets>=2 && hyp_type==0") , luminosity, 40, 0., 200, 0, babyVector);

    //same as above, but integral plots
    DrawAll("mass",     prefix+"_osz_mass_ee_2jets_int",             TCut("ee {NJets>=2}", osz_dilep_all + "njets>=2 && hyp_type==3")             , luminosity, 40, 0., 200, 1, babyVector);
    DrawAll("mass",     prefix+"_osz_mass_ee_2jets_pfmet60_int",     TCut("ee {NJets>=2, PFMet>60}", osz_dilep_all + "njets>=2 && hyp_type==3 && pfmet>60") , luminosity, 40, 0., 200, 1, babyVector);
    DrawAll("mass",     prefix+"_osz_mass_mm_2jets_int",             TCut("mm {NJets>=2}", osz_dilep_all + "njets>=2 && hyp_type==0")             , luminosity, 40, 0., 200, 1, babyVector);
    DrawAll("mass",     prefix+"_osz_mass_mm_2jets_pfmet60_int",     TCut("mm {NJets>=2, PFMet>60}", osz_dilep_all + "njets>=2 && hyp_type==0 && pfmet>60") , luminosity, 40, 0., 200, 1, babyVector);
    DrawAll("pfmet",    prefix+"_osz_pfmet_ee_2jets_int",            TCut("ee {NJets>=2}", osz_dilep_all + osz_zmass + "njets>=2 && hyp_type==3") , luminosity, 40, 0., 200, 1, babyVector);
    DrawAll("pfmet",    prefix+"_osz_pfmet_mm_2jets_int",            TCut("mm {NJets>=2}", osz_dilep_all + osz_zmass + "njets>=2 && hyp_type==0") , luminosity, 40, 0., 200, 1, babyVector);

}

