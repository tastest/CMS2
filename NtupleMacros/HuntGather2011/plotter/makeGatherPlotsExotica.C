#include "cuts.h"

void makeGatherPlotsExotica(TString prefix, const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making Exotica plots...\n";

    //
    // Define specific cuts for Exotica plots
    //

    TCut cut_z_dijets         ("z_dijets", inclusivez_dilep+"njets>=2");
    TCut cut_z_highptdijets   ("z_highptdijets", inclusivez_dilep+"njets>=2&&jet1pt>150.&&jet2pt>150.");

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, inclusivez_dilep);

    //
    // Make the plots
    //

    //
    // Misc
    //

    DrawAll("jetmass",prefix+"_exotica_z_highptdijets",cut_z_highptdijets,luminosity,40,0.,2000.,0, babyVector);
    DrawAll("jetmass",prefix+"_exotica_z_dijets",cut_z_dijets,luminosity,40,0.,2000.,0, babyVector);

    //
    // High pT Z (AN2011/115)
    //

    //
    // Region 1 (njets >=0)
    //

    TCut cut_r1ee("Region 1 ee", ee_dilep+inclusivez_dilep+"njets>=0");
    TCut cut_r1mm("Region 1 mm", mm_dilep+inclusivez_dilep+"njets>=0");
    DrawAll("dilpt",prefix+"_exotica_R1_zee_pt", cut_r1ee, luminosity, 32, 40., 200., 0, babyVector);
    DrawAll("dilpt",prefix+"_exotica_R1_zmm_pt", cut_r1mm, luminosity, 32, 40., 200., 0, babyVector);

    //
    // Region 2 (njets == 1 and |eta_Z| < 1.0)
    //

    // z pT distributions 
    TCut cut_r2ee("Region 2 ee", ee_dilep+inclusivez_dilep+"njets==1 && abs(dileta) < 1.0");
    TCut cut_r2mm("Region 2 mm", mm_dilep+inclusivez_dilep+"njets==1 && abs(dileta) < 1.0");
    DrawAll("dilpt",prefix+"_exotica_R2_zee_pt", cut_r2ee,luminosity,32,40.,200.,0, babyVector);
    DrawAll("dilpt",prefix+"_exotica_R2_zmm_pt", cut_r2mm,luminosity,32,40.,200.,0, babyVector);

    // mll-j mass distributions for 60 < pTZ < 140
    TCut cut_r2ee_1("Region 2 ee 60<pTZ<140", ee_dilep+inclusivez_dilep+"njets==1 && abs(dileta) < 1.0 && dilpt>60&&dilpt<140");
    TCut cut_r2mm_1("Region 2 mm 60<pTZ<140", mm_dilep+inclusivez_dilep+"njets==1 && abs(dileta) < 1.0 && dilpt>60&&dilpt<140");
    DrawAll("mllj",prefix+"_exotica_R2_zee_mllj", cut_r2ee_1, luminosity, 40, 0., 500., 0, babyVector);
    DrawAll("mllj",prefix+"_exotica_R2_zmm_mllj", cut_r2mm_1, luminosity, 40, 0., 500., 0, babyVector);


}

