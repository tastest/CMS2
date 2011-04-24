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

    PreselectBabies(babyVector, base_dilep);

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
    //TCut cut_r1ll("Region 1 ll", inclusivez_dilep+"njets>=0");
    DrawAll("dilpt",prefix+"_exotica_R1_zee_pt", cut_r1ee, luminosity, 32, 40., 200., 0, babyVector);
    DrawAll("dilpt",prefix+"_exotica_R1_zmm_pt", cut_r1mm, luminosity, 32, 40., 200., 0, babyVector);
    //DrawAll("dilpt",prefix+"_exotica_R1_zll_pt", cut_r1ll, luminosity, 32, 40., 200., 0, babyVector);
 
    // mll-j mass distributions for 60 < pTZ < 140
    TCut cut_r1ee_1("Region 1 ee 60<pTZ<140", ee_dilep+inclusivez_dilep+"njets>=0 && dilpt>60&&dilpt<140");
    TCut cut_r1mm_1("Region 1 mm 60<pTZ<140", mm_dilep+inclusivez_dilep+"njets>=0 && dilpt>60&&dilpt<140");
    //TCut cut_r1ll_1("Region 1 ll 60<pTZ<140", inclusivez_dilep+"njets>=0 && dilpt>60&&dilpt<140");
    DrawAll("mllj",prefix+"_exotica_R1_zee_mllj", cut_r1ee_1, luminosity, 40, 0., 500., 0, babyVector);
    DrawAll("mllj",prefix+"_exotica_R1_zmm_mllj", cut_r1mm_1, luminosity, 40, 0., 500., 0, babyVector);
    //DrawAll("mllj",prefix+"_exotica_R1_zll_mllj", cut_r1ll_1, luminosity, 40, 0., 500., 0, babyVector);
 
    // 
    // Region 2 (njets == 1 and |eta_Z| < 1.0)
    //

    // z pT distributions 
    TCut cut_r2ee("Region 2 ee", ee_dilep+inclusivez_dilep+"njets==1 && abs(dileta) < 1.0");
    TCut cut_r2mm("Region 2 mm", mm_dilep+inclusivez_dilep+"njets==1 && abs(dileta) < 1.0");
    //TCut cut_r2ll("Region 2 ll", inclusivez_dilep+"njets==1 && abs(dileta) < 1.0");
    DrawAll("dilpt",prefix+"_exotica_R2_zee_pt", cut_r2ee,luminosity,32,40.,200.,0, babyVector);
    DrawAll("dilpt",prefix+"_exotica_R2_zmm_pt", cut_r2mm,luminosity,32,40.,200.,0, babyVector);
    //DrawAll("dilpt",prefix+"_exotica_R2_zll_pt", cut_r2ll,luminosity,32,40.,200.,0, babyVector);

    // mll-j mass distributions for 60 < pTZ < 140
    TCut cut_r2ee_1("Region 2 ee 60<pTZ<140", ee_dilep+inclusivez_dilep+"njets==1 && abs(dileta) < 1.0 && dilpt>60&&dilpt<140");
    TCut cut_r2mm_1("Region 2 mm 60<pTZ<140", mm_dilep+inclusivez_dilep+"njets==1 && abs(dileta) < 1.0 && dilpt>60&&dilpt<140");
    //TCut cut_r2ll_1("Region 2 ll 60<pTZ<140", inclusivez_dilep+"njets==1 && abs(dileta) < 1.0 && dilpt>60&&dilpt<140");
    DrawAll("mllj",prefix+"_exotica_R2_zee_mllj", cut_r2ee_1, luminosity, 40, 0., 500., 0, babyVector);
    DrawAll("mllj",prefix+"_exotica_R2_zmm_mllj", cut_r2mm_1, luminosity, 40, 0., 500., 0, babyVector);
    //DrawAll("mllj",prefix+"_exotica_R2_zll_mllj", cut_r2ll_1, luminosity, 40, 0., 500., 0, babyVector);

    // Slavas special Z massage
    TCut cut_r3("Region 3 Slavas cuts", "eormu1*eormu2==-169&& mu1_muonidfull && mu2_muonidfull && abs(mass-91)<10 && njets==1 && abs(dileta)<1 && pt1>15 && pt2>15 && acos(cos(phi2-phi1))>1 &&acos(cos(phi2-phi1))<2.25 && pfmet<20 && mu1_emVetoDep<4&&mu2_emVetoDep<4 &&mu1_hadVetoDep<6&&mu2_hadVetoDep<6 && mu1_saHits>10&&mu2_saHits>10 && ecalIso1<2 && ecalIso2<2 && hcalIso1<2 && hcalIso2<2 && trkIso1<3 && trkIso2<3");

    DrawAll("dilpt",prefix+"_exotica_R3_zmm", cut_r3, luminosity, 32, 40., 200., 0, babyVector);

/*
    TCanvas *c1 = DrawAll("dilpt",prefix+"_exotica_R3_zmm", cut_r3, luminosity, 32, 40., 200., 0, babyVector);
    TH1F *h1_bg = (TH1F*)c1->FindObject("dyeemm_Region 3 Slavas cuts_dilpt_0")->Clone("h1_bg");
    TH1F *h1_data = (TH1F*)c1->FindObject("data_Region 3 Slavas cuts_dilpt_0")->Clone("h1_data");

    Float_t mcInt = h1_bg->Integral(0, 6);
    Float_t dataInt = h1_data->Integral(0, 6);
    TH1F *h1_bg_scaled = (TH1F*)h1_bg->Clone("h1_bg_scaled");
    h1_bg_scaled->Scale(dataInt/mcInt);

    TCanvas *c2 = new TCanvas(prefix+"_exotica_R3_zmm_scaled");
    c2->Divide(1, 2);
    c2->cd(1);
    h1_bg->Draw("HIST");
    h1_data->Draw("SAME E1");
    h1_bg->GetYaxis()->SetRangeUser(0, 25);
    c2->cd(2);
    h1_bg_scaled->Draw("HIST");
    h1_data->Draw("SAME E1");
    h1_bg_scaled->GetYaxis()->SetRangeUser(0, 25);
    Float_t mcIntScaled_signal = h1_bg->Integral((85-40)/5, (105-40)/5);
    Float_t dataIntScaled_signal = h1_data->Integral((85-40)/5, (105-40)/5);
    std::cout << "85-105: " << mcIntScaled_signal << ", " << dataIntScaled_signal << std::endl;
*/

/*
    TCut zprime_1("zprime_1", "pt1 > 30 && (abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
    TCut zprime_2("zprime_2", "pt2 > 30 && (abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
    TCut zprime = zprime_1 + zprime_2;

    DrawAll("mass",prefix+"_exotica_zprime", zprime, luminosity, 50, 0.0, 1000., 0, babyVector);
    DrawAll("mass",prefix+"_exotica_zprime_highmass", zprime, luminosity, 40, 200.0, 1000., 0, babyVector);
    DrawAll("mass",prefix+"_exotica_zprime_highmass_ee", zprime+ee_dilep, luminosity, 40, 200.0, 1000., 0, babyVector);
    DrawAll("mass",prefix+"_exotica_zprime_highmass_mm", zprime+mm_dilep, luminosity, 40, 200.0, 1000., 0, babyVector);
*/

}

