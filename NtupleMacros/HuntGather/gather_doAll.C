#include "cuts.h"

void gather_doAll() {
    gROOT->ProcessLine(".L tdrstyle.C");
    setTDRStyle();

    //gROOT ->SetStyle("Plain");
    //gStyle->SetHistMinimumZero();
    //gStyle->SetOptStat(0);

    gROOT->ProcessLine(".L goodrun.cc+");
    gROOT->ProcessLine(".L BabyDorkIdentifier.C+");
    gROOT->ProcessLine(".L gather.C+");

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
    std::cout << "Integrated luminosity: " << 1e3*f_intlumifb << "/pb\n";

    //
    // OS PLOTS
    //

    std::cout << "Making OS plots...\n";

    TCut osanal_dilep   ("osanal_dilep"   ,base_dilep+os_dilep+"tcmet>50.&&(jet1pt+jet2pt+jet3pt)>150.");
    TCut osanal_of_dilep("osanal_of_dilep",osanal_dilep+of_dilep);
    TCut osanal_sf_dilep("osanal_sf_dilep",osanal_dilep+sf_dilep);

    DrawAll("mass",osanal_of_dilep,f_intlumifb,40,0.,500.,0);
    DrawAll("mass",osanal_sf_dilep,f_intlumifb,40,0.,500.,0);
    DrawAll("jet1pt+jet2pt+jet3pt",osanal_dilep,f_intlumifb,40,0.,800.,1);
    DrawAll("tcmet",osanal_dilep,f_intlumifb,40,0.,300.,1);

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

    DrawAll("mass",ssanal_ee_dilep,f_intlumifb,40,0.,500.,0);
    DrawAll("mass",ssanal_mm_dilep,f_intlumifb,40,0.,500.,0);
    DrawAll("jet1pt+jet2pt+jet3pt",ssanal_dilep,f_intlumifb,40,0.,800.,1);
    DrawAll("tcmet",ssanal_dilep,f_intlumifb,40,0.,300.,1);

    //
    // Z+MET
    //

    std::cout << "Making Z+MET plots...\n";

    TCut zmet_os_0j_dilep("zmet_os_0j_dilep",inclusivez_dilep+os_dilep+"njetsClean==0");
    TCut zmet_os_1j_dilep("zmet_os_1j_dilep",inclusivez_dilep+os_dilep+"njetsClean==1");
    TCut zmet_os_2j_dilep("zmet_os_2j_dilep",inclusivez_dilep+os_dilep+"njetsClean>=2");

    DrawAll("tcmet",zmet_os_0j_dilep,f_intlumifb,40,0.,300.,1);
    DrawAll("tcmet",zmet_os_1j_dilep,f_intlumifb,40,0.,300.,1);
    DrawAll("tcmet",zmet_os_2j_dilep,f_intlumifb,40,0.,300.,1);
    DrawAll("pfmet",zmet_os_0j_dilep,f_intlumifb,40,0.,300.,1);
    DrawAll("pfmet",zmet_os_1j_dilep,f_intlumifb,40,0.,300.,1);
    DrawAll("pfmet",zmet_os_2j_dilep,f_intlumifb,40,0.,300.,1);

    //
    // Effective Mass
    //

    std::cout << "Making other plots...\n";

    TCut base_tcmeffgt400_dilep  ("base_tcmeffgt400_dilep",base_dilep+"tcmeff>400.");
    TCut base_tcmetgt50_dilep    ("base_tcmetgt50_dilep",base_dilep+"tcmet>50.");
    TCut base_sumjetptgt200_dilep("base_sumjetptgt200_dilep",base_dilep+"(jet1pt+jet2pt+jet3pt)>50.");
    TCut base_dilptgt100_dilep   ("base_dilptgt100_dilep",base_dilep+"dilpt>100.");

    DrawAll("njetsClean",base_dilep,f_intlumifb,10,-0.5,9.5,0);
    DrawAll("njetsClean",base_tcmeffgt400_dilep,f_intlumifb,10,-0.5,9.5,0);
    DrawAll("tcmeff",base_dilep,f_intlumifb,40,0.,1000.,1);
    DrawAll("tcmeff",base_tcmetgt50_dilep,f_intlumifb,40,0.,1000.,1);
    DrawAll("tcmeff",base_sumjetptgt200_dilep,f_intlumifb,40,0.,1000.,1);
    DrawAll("tcmeff",base_dilptgt100_dilep,f_intlumifb,40,0.,1000.,1);

    TSeqCollection *list = gROOT->GetListOfCanvases();
    TIterator *iter = list->MakeIterator();
    TCanvas *c1 = 0;

    while ((c1 = (TCanvas*)iter->Next())) {
        c1->Print(Form("%s.png", c1->GetName()));
        c1->SetLogy(1);
        c1->Print(Form("%s_log.png", c1->GetName()));
    }
}
