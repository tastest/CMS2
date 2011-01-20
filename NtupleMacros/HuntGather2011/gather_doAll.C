#include "cuts.h"

void gather_doAll() {
    gROOT->ProcessLine(".L tdrstyle.C");
    setTDRStyle();

    //gROOT ->SetStyle("Plain");
    //gStyle->SetHistMinimumZero();
    //gStyle->SetOptStat(0);

    gROOT->ProcessLine(".L ../Tools/goodrun.cc+");
    gROOT->ProcessLine(".L gather.C+");

    //
    // configuration
    //

    //const char *goodrunlist = "runlists/Cert_TopOct6_Merged_135059-146729_allPVT_extra_146804-147116.txt";
    //float goodrunlumi = 5860;
    
    // const char *goodrunlist = "runlists/Cert_TopOct8_Merged_135059-147116_allPVT.txt";
    // float goodrunlumi = 7090;
    
    //   const char *goodrunlist = "runlists/Cert_TopOct15_Merged_135821-147454_allPVT.txt";
    // float goodrunlumi = 11060;
    //const char *goodrunlist = "runlists/Cert_TopOct22_Merged_135821-148058_allPVT.txt";
    //    float goodrunlumi = 15210;
    const char *goodrunlist = "runlists/json_135821_148864_22.02pb.txt";
    float goodrunlumi = 22020;

    // this sets the json file obviously; it's
    // easiest to just  soft link to  the real
    // file as json.txt
    std::cout << "Using " << goodrunlist << " for goodruns\n";
    set_goodrun_file_json(goodrunlist);
    unsigned int lastgoodrun = max_run();
    unsigned int lastgoodlumi = max_run_max_lumi();

    // calculate  integrated luminosity  in /fb
    // note that I input /nb and so get out /nb
    // then scale to /fb which  is what  I need
    // 1e-3*GetIntLumi(/pb)  is the  same thing
    // doesn't matter so long as it corresponds
    // to the lumi of the json file and that it
    // is correctly scaled to /fb afterward
    float f_intlumifb = 1e-6*GetIntLumi(goodrunlumi);

    std::cout << "Integrated luminosity: " << 1e3*f_intlumifb << "/pb\n";

    //
    // validation plots
    //

    std::cout << "Making validation plots...\n";

    // lumi that is in the good run list
    float f_goodruns_intlumifb = 1e-6*goodrunlumi;
    // lumi that is out of the good run list
    float f_newruns_intlumifb = f_intlumifb - f_goodruns_intlumifb;

    TCut validation_ee ("validation_ee", base_dilep+ee_dilep);
    TCut validation_mm ("validation_mm", base_dilep+mm_dilep);

    DrawAll("mass", "validation_mass_goodruns_ee", validation_ee, Form("!isdata||(run < %i || (run == %i && ls <= %i))", lastgoodrun, lastgoodrun, lastgoodlumi), f_goodruns_intlumifb, 50,0., 200., 0);
    DrawAll("mass", "validation_mass_newruns_ee", validation_ee, Form("!isdata||(run > %i || (run == %i && ls > %i))", lastgoodrun, lastgoodrun, lastgoodlumi), f_newruns_intlumifb, 50,0., 200., 0);
    DrawAll("mass", "validation_mass_goodruns_mm", validation_mm, Form("!isdata||(run < %i || (run == %i && ls <= %i))", lastgoodrun, lastgoodrun, lastgoodlumi), f_goodruns_intlumifb, 50,0., 200., 0);
    DrawAll("mass", "validation_mass_newruns_mm", validation_mm, Form("!isdata||(run > %i || (run == %i && ls > %i))", lastgoodrun, lastgoodrun, lastgoodlumi), f_newruns_intlumifb, 50,0., 200., 0);

    //
    // OS PLOTS
    //
    std::cout << "Making OS plots...\n";

    TCut osanal_dilep   ("osanal_dilep"   ,base_dilep+os_dilep+"tcmet>50.&&sumjetpt>150.");
    TCut osanal_of_dilep("osanal_of_dilep",osanal_dilep+of_dilep);
    TCut osanal_sf_dilep("osanal_sf_dilep",osanal_dilep+sf_dilep);

    DrawAll("mass","os_of_mass",osanal_of_dilep,"",f_intlumifb,40,0.,500.,0);
    DrawAll("mass","os_sf_mass",osanal_sf_dilep,"",f_intlumifb,40,0.,500.,0);
    DrawAll("sumjetpt","os_sumjetpt_int",osanal_dilep,"",f_intlumifb,40,0.,800.,1);
    DrawAll("tcmet","os_tcmet_int",osanal_dilep,"",f_intlumifb,40,0.,300.,1);
    
    // 4 Lepton events search
    // TCut exotic_4lepton   ("exotic_4lepton"   ,exotic_dilep);
    // DrawAll("mass","exotic_dilep_mass",exotic_4lepton,"",f_intlumifb,40,0.,500.,0);
    //
    // SS PLOTS
    //

    std::cout << "Making SS plots...\n";

    TCut ssanal_lepid1 = "(abs(eormu1)==11&&e1_ctfCharge==e1_scCharge&&e1_ctfCharge==e1_gsfCharge&&e1_vbtf70) || abs(eormu1)==13";
    TCut ssanal_lepid2 = "(abs(eormu2)==11&&e2_ctfCharge==e2_scCharge&&e2_ctfCharge==e2_gsfCharge&&e2_vbtf70) || abs(eormu2)==13";
    TCut ssanal_lepid  = ssanal_lepid1+ssanal_lepid2;
    TCut ssanal_dilep   ("ssanal_dilep"   ,base_dilep+ss_dilep+ssanal_lepid);
    TCut ssanal_ee_dilep("ssanal_ee_dilep",ssanal_dilep+ee_dilep);
    TCut ssanal_mm_dilep("ssanal_mm_dilep",ssanal_dilep+mm_dilep);

    DrawAll("mass","ss_ee_mass",ssanal_ee_dilep,"",f_intlumifb,40,0.,500.,0);
    DrawAll("mass","ss_mm_mass",ssanal_mm_dilep,"",f_intlumifb,40,0.,500.,0);
    DrawAll("sumjetpt","ss_sumjetpt_int",ssanal_dilep,"",f_intlumifb,40,0.,800.,1);
    DrawAll("tcmet","ss_tcmet_int",ssanal_dilep,"",f_intlumifb,40,0.,300.,1);

    //
    // Z+MET
    //

    std::cout << "Making Z+MET plots...\n";
    TCut zmet_os_0j_dilep("zmet_os_0j_dilep",inclusivez_dilep+os_dilep+"njetsClean==0");
    TCut zmet_os_1j_dilep("zmet_os_1j_dilep",inclusivez_dilep+os_dilep+"njetsClean==1");
    TCut zmet_os_2j_dilep("zmet_os_2j_dilep",inclusivez_dilep+os_dilep+"njetsClean>=2");
    TCut zmet_sumjetptgt100_os_ee_dilep("zmet_sumjetptgt100_os_ee_dilep",inclusivez_dilep+os_dilep+ee_dilep+"sumjetpt>100.");
    TCut zmet_sumjetptgt100_os_mm_dilep("zmet_sumjetptgt100_os_ee_dilep",inclusivez_dilep+os_dilep+mm_dilep+"sumjetpt>100.");
    DrawAll("tcmet","zmet_os_0j_tcmet_int",zmet_os_0j_dilep,"",f_intlumifb,40,0.,300.,1);
    DrawAll("tcmet","zmet_os_1j_tcmet_int",zmet_os_1j_dilep,"",f_intlumifb,40,0.,300.,1);
    DrawAll("tcmet","zmet_os_2j_tcmet_int",zmet_os_2j_dilep,"",f_intlumifb,40,0.,300.,1);
    DrawAll("pfmet","zmet_os_0j_pfmet_int",zmet_os_0j_dilep,"",f_intlumifb,40,0.,300.,1);
    DrawAll("pfmet","zmet_os_1j_pfmet_int",zmet_os_1j_dilep,"",f_intlumifb,40,0.,300.,1);
    DrawAll("pfmet","zmet_os_2j_pfmet_int",zmet_os_2j_dilep,"",f_intlumifb,40,0.,300.,1);
    DrawAll("tcmet","zmet_sumjetptgt100_os_ee_tcmet_int",zmet_sumjetptgt100_os_ee_dilep,"",f_intlumifb,40,0.,300.,1);
    DrawAll("pfmet","zmet_sumjetptgt100_os_ee_pfmet_int",zmet_sumjetptgt100_os_ee_dilep,"",f_intlumifb,40,0.,300.,1);
    DrawAll("tcmet","zmet_sumjetptgt100_os_mm_tcmet_int",zmet_sumjetptgt100_os_mm_dilep,"",f_intlumifb,40,0.,300.,1);
    DrawAll("pfmet","zmet_sumjetptgt100_os_mm_pfmet_int",zmet_sumjetptgt100_os_mm_dilep,"",f_intlumifb,40,0.,300.,1);

    //
    // Effective Mass
    //

    std::cout << "Making other plots...\n";

    TCut base_tcmeffgt400_dilep  ("base_tcmeffgt400_dilep",base_dilep+"tcmeff>400.");
    TCut base_tcmetgt50_dilep    ("base_tcmetgt50_dilep",base_dilep+"tcmet>50.");
    TCut base_sumjetptgt200_dilep("base_sumjetptgt200_dilep",base_dilep+"sumjetpt>200.");
    TCut base_dilptgt100_dilep   ("base_dilptgt100_dilep",base_dilep+"dilpt>100.");
    DrawAll("njetsClean","meff_njetsclean",base_dilep,"",f_intlumifb,10,-0.5,9.5,0);
    DrawAll("njetsClean","meff_tcmeffgt400_njetsclean",base_tcmeffgt400_dilep,"",f_intlumifb,10,-0.5,9.5,0);
    DrawAll("tcmeff","meff_tcmeff_int",base_dilep,"",f_intlumifb,40,0.,1000.,1);
    DrawAll("tcmeff","meff_tcmetgt50_tcmeff_int",base_tcmetgt50_dilep,"",f_intlumifb,40,0.,1000.,1);
    DrawAll("tcmeff","meff_sumjetptgt200_tcmeff_int",base_sumjetptgt200_dilep,"",f_intlumifb,40,0.,1000.,1);
    DrawAll("tcmeff","meff_dilptgt100_tcmeff_int",base_dilptgt100_dilep,"",f_intlumifb,40,0.,1000.,1);

    //
    // Exotica 
    //

    std::cout << "Making additional plots...\n";

    TCut cut_z_dijets         ("z_dijets", inclusivez_dilep+"njetsClean>=2");
    TCut cut_z_highptdijets   ("z_highptdijets", inclusivez_dilep+"njetsClean>=2&&jet1pt>150.&&jet2pt>150.");
    DrawAll("jetmass","exotica_z_highptdijets",cut_z_highptdijets,"",f_intlumifb,40,0.,2000.,0);
    DrawAll("jetmass","exotica_z_dijets",cut_z_dijets,"",f_intlumifb,40,0.,2000.,0);


    TSeqCollection *list = gROOT->GetListOfCanvases();
    TIterator *iter = list->MakeIterator();
    TCanvas *c1 = 0;

    while ((c1 = (TCanvas*)iter->Next())) {
        c1->Print(Form("%s.png", c1->GetName()));
        c1->SetLogy(1);
        c1->Print(Form("%s_log.png", c1->GetName()));
    }
}
