#ifndef CUTS_H
#define CUTS_H

#include "TCut.h"

TCut wplusjetsv1_emu1("pt > 20.");
TCut wplusjetsv1_emu2("(abs(eormu) == 11 && e_vbtf90full) || (abs(eormu) == 13 && mu_muonidfull && ! mu_cosmic)");
TCut wplusjetsv1_emu3("pfmet > 20. && tcmet > 20.");
TCut wplusjetsv1_emu4("mt > 50.");
TCut wplusjetsv1_emu5("njets >= 2");
TCut wplusjetsv1_emu6("ntrks > 2");
TCut wplusjetsv1_emu("wplusjetsv1_emu",wplusjetsv1_emu1+wplusjetsv1_emu2+wplusjetsv1_emu3+wplusjetsv1_emu4+wplusjetsv1_emu5+wplusjetsv1_emu6);

TCut wprime100_emu1("pt > 20.");
TCut wprime100_emu2("(abs(eormu) == 11 && e_vbtf90full) || (abs(eormu) == 13 && mu_muonidfull && ! mu_cosmic)");
TCut wprime100_emu3("pfmet > 100. && tcmet > 100.");
TCut wprime100_emu4("mt > 100.");
TCut wprime100_emu5("ntrks > 2");
TCut wprime100_emu("wprime100_emu",wprime100_emu1+wprime100_emu2+wprime100_emu3+wprime100_emu4+wprime100_emu5);

TCut semileptonictop_emu1("pt > 20.");
TCut semileptonictop_emu2("(abs(eormu) == 11 && e_vbtf90full) || (abs(eormu) == 13 && mu_muonidfull && ! mu_cosmic)");
TCut semileptonictop_emu3("njets >= 3");
TCut semileptonictop_emu4("ngoodlep < 2");
TCut semileptonictop_emu5("ntrks > 2");
TCut semileptonictop_emu("semileptonictop_emu",semileptonictop_emu1+semileptonictop_emu2+semileptonictop_emu3+semileptonictop_emu4+semileptonictop_emu5);

TCut semileptonictopv2_emu1("pt > 20.");
TCut semileptonictopv2_emu2("(abs(eormu) == 11 && e_vbtf90full) || (abs(eormu) == 13 && mu_muonidfull && ! mu_cosmic)");
TCut semileptonictopv2_emu3("(njets >= 4 && neffbtags>=1) || (njets >=3 && neffbtags>=2)");
TCut semileptonictopv2_emu4("ngoodlep < 2");
TCut semileptonictopv2_emu5("ntrks > 2");
TCut semileptonictopv2_emu("semileptonictopv2_emu",semileptonictopv2_emu1+semileptonictopv2_emu2+semileptonictopv2_emu3+semileptonictopv2_emu4+semileptonictopv2_emu5);

TCut semileptonictopv3_emu1("pt > 20.");
TCut semileptonictopv3_emu2("(abs(eormu) == 11 && e_vbtf90full) || (abs(eormu) == 13 && mu_muonidfull && ! mu_cosmic)");
TCut semileptonictopv3_emu3("njets >= 4 && npurbtags > 1");
TCut semileptonictopv3_emu4("ngoodlep < 2");
TCut semileptonictopv3_emu5("ntrks > 2");
TCut semileptonictopv3_emu("semileptonictopv3_emu",semileptonictopv3_emu1+semileptonictopv3_emu2+semileptonictopv3_emu3+semileptonictopv3_emu4+semileptonictopv3_emu5);

TCut highptmuon_emu1("pt > 150.");
TCut highptmuon_emu2("abs(eormu) == 13 && mu_muonid && !mu_cosmic");
TCut highptmuon_emu3("ntrks > 2");
TCut highptmuon_emu("highptmuon_emu",highptmuon_emu1+highptmuon_emu2+highptmuon_emu3);

TCut highptmuonv2_emu1("pt > 200.");
TCut highptmuonv2_emu2("abs(eormu) == 13 && mu_muonid && !mu_cosmic");
TCut highptmuonv2_emu3("ntrks > 2");
TCut highptmuonv2_emu("highptmuonv2_emu",highptmuonv2_emu1+highptmuonv2_emu2+highptmuonv2_emu3);

TCut multimuon_emu1("ngoodmus > 2");
TCut multimuon_emu2("ntrks > 2");
TCut multimuon_emu("multimuon_emu",multimuon_emu1+multimuon_emu2);

//
// Dilepton+X
//

TCut base_dilep1("pt1 > 20. && pt2 > 20.");
TCut base_dilep2("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut base_dilep3("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut base_dilep4("ntrks > 2");
TCut base_dilep("base_dilep",base_dilep1+base_dilep2+base_dilep3+base_dilep4);

TCut ss_dilep("eormu1*eormu2>0");
TCut os_dilep("eormu1*eormu2<0");
TCut sf_dilep("abs(eormu1)==abs(eormu2)");
TCut of_dilep("abs(eormu1)!=abs(eormu2)");

TCut ee_dilep("abs(eormu1)==11&&abs(eormu2)==11");
TCut mm_dilep("abs(eormu1)==13&&abs(eormu2)==13");
TCut em_dilep("(abs(eormu1)==11&&abs(eormu2)==13)||(abs(eormu1)==13&&abs(eormu2)==11)");

//TCut exotic_dilep1("(pt1 > 20. && pt2 > 10.) || (pt1 > 10. && pt2 > 20.)");

TCut inclusivez_dilep1("mass > 76. && mass < 106.");
TCut inclusivez_dilep("inclusivez_dilep",base_dilep+sf_dilep+inclusivez_dilep1);

TCut inclusivenonz_dilep1("(((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 < 0) && abs(mass-91.) > 15.) || ((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 > 0) || (hyp_type == 1 || hyp_type == 2)");
TCut inclusivenonz_dilep("inclusivenonz_dilep",base_dilep+inclusivenonz_dilep1);

TCut zplusjetsv1_dilep1("abs(eormu1) == abs(eormu2)");
TCut zplusjetsv1_dilep2("pt1 > 20. && pt2 > 20.");
TCut zplusjetsv1_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut zplusjetsv1_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut zplusjetsv1_dilep5("mass > 76. && mass < 106.");
TCut zplusjetsv1_dilep6("njets >= 1");
TCut zplusjetsv1_dilep7("ntrks > 2");
TCut zplusjetsv1_dilep("zplusjetsv1_dilep",zplusjetsv1_dilep1+zplusjetsv1_dilep2+zplusjetsv1_dilep3+zplusjetsv1_dilep4+zplusjetsv1_dilep5+zplusjetsv1_dilep6+zplusjetsv1_dilep7);

TCut zplusjetsv2_dilep1("abs(eormu1) == abs(eormu2)");
TCut zplusjetsv2_dilep2("pt1 > 20. && pt2 > 20.");
TCut zplusjetsv2_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut zplusjetsv2_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut zplusjetsv2_dilep5("mass > 76. && mass < 106.");
TCut zplusjetsv2_dilep6("njets >= 3 || jet1pt > 100 || (njets>=2 && ((jet1pt+jet2pt)>100))");
TCut zplusjetsv2_dilep7("ntrks > 2");
TCut zplusjetsv2_dilep("zplusjetsv2_dilep",zplusjetsv2_dilep1+zplusjetsv2_dilep2+zplusjetsv2_dilep3+zplusjetsv2_dilep4+zplusjetsv2_dilep5+zplusjetsv2_dilep6+zplusjetsv2_dilep7);

TCut zplusjetsv3_dilep1("abs(eormu1) == abs(eormu2)");
TCut zplusjetsv3_dilep2("pt1 > 20. && pt2 > 20.");
TCut zplusjetsv3_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut zplusjetsv3_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut zplusjetsv3_dilep5("mass > 76. && mass < 106.");
TCut zplusjetsv3_dilep6("(jet1pt + jet2pt + jet3pt) > 200.");
TCut zplusjetsv3_dilep7("ntrks > 2");
TCut zplusjetsv3_dilep("zplusjetsv3_dilep",zplusjetsv3_dilep1+zplusjetsv3_dilep2+zplusjetsv3_dilep3+zplusjetsv3_dilep4+zplusjetsv3_dilep5+zplusjetsv3_dilep6+zplusjetsv3_dilep7);

TCut zplusbjetsv1_dilep1("abs(eormu1) == abs(eormu2)");
TCut zplusbjetsv1_dilep2("pt1 > 20. && pt2 > 20.");
TCut zplusbjetsv1_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut zplusbjetsv1_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut zplusbjetsv1_dilep5("mass > 76. && mass < 106.");
TCut zplusbjetsv1_dilep6("njets>=1 && neffbtags>0");
TCut zplusbjetsv1_dilep7("ntrks > 2");
TCut zplusbjetsv1_dilep("zplusbjetsv1_dilep",zplusbjetsv1_dilep1+zplusbjetsv1_dilep2+zplusbjetsv1_dilep3+zplusbjetsv1_dilep4+zplusbjetsv1_dilep5+zplusbjetsv1_dilep6+zplusbjetsv1_dilep7);

TCut zplusbjetsv2_dilep1("abs(eormu1) == abs(eormu2)");
TCut zplusbjetsv2_dilep2("pt1 > 20. && pt2 > 20.");
TCut zplusbjetsv2_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut zplusbjetsv2_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut zplusbjetsv2_dilep5("mass > 76. && mass < 106.");
TCut zplusbjetsv2_dilep6("njets>=1 && (npurbtags>0 || neffbtags>1)");
TCut zplusbjetsv2_dilep7("ntrks > 2");
TCut zplusbjetsv2_dilep("zplusbjetsv2_dilep",zplusbjetsv2_dilep1+zplusbjetsv2_dilep2+zplusbjetsv2_dilep3+zplusbjetsv2_dilep4+zplusbjetsv2_dilep5+zplusbjetsv2_dilep6+zplusbjetsv2_dilep7);

TCut zmet30_dilep1("abs(eormu1) == abs(eormu2)");
TCut zmet30_dilep2("pt1 > 20. && pt2 > 20.");
TCut zmet30_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut zmet30_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut zmet30_dilep5("mass > 76. && mass < 106.");
TCut zmet30_dilep6("pfmet > 30. || tcmet > 30.");
TCut zmet30_dilep7("ntrks > 2");
TCut zmet30_dilep("zmet30_dilep",zmet30_dilep1+zmet30_dilep2+zmet30_dilep3+zmet30_dilep4+zmet30_dilep5+zmet30_dilep6+zmet30_dilep7);

TCut zptz50_dilep1("abs(eormu1) == abs(eormu2)");
TCut zptz50_dilep2("pt1 > 20. && pt2 > 20.");
TCut zptz50_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut zptz50_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut zptz50_dilep5("mass > 76. && mass < 106.");
TCut zptz50_dilep6("dilpt > 50.");
TCut zptz50_dilep7("ntrks > 2");
TCut zptz50_dilep("zptz50_dilep",zptz50_dilep1+zptz50_dilep2+zptz50_dilep3+zptz50_dilep4+zptz50_dilep5+zptz50_dilep6+zptz50_dilep7);

TCut zptz100_dilep1("abs(eormu1) == abs(eormu2)");
TCut zptz100_dilep2("pt1 > 20. && pt2 > 20.");
TCut zptz100_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut zptz100_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut zptz100_dilep5("mass > 76. && mass < 106.");
TCut zptz100_dilep6("dilpt > 100.");
TCut zptz100_dilep7("ntrks > 2");
TCut zptz100_dilep("zptz100_dilep",zptz100_dilep1+zptz100_dilep2+zptz100_dilep3+zptz100_dilep4+zptz100_dilep5+zptz100_dilep6+zptz100_dilep7);

TCut inclusivenonzof_dilep1("pt1 > 15. && pt2 > 15.");
TCut inclusivenonzof_dilep2("abs(eormu1) != abs(eormu2)");
TCut inclusivenonzof_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut inclusivenonzof_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut inclusivenonzof_dilep5("(((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 < 0) && abs(mass-91.) > 15.) || ((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 > 0) || (hyp_type == 1 || hyp_type == 2)");
TCut inclusivenonzof_dilep6("ntrks > 2");
TCut inclusivenonzof_dilep("inclusivenonzof_dilep",inclusivenonzof_dilep1+inclusivenonzof_dilep2+inclusivenonzof_dilep3+inclusivenonzof_dilep4+inclusivenonzof_dilep5+inclusivenonzof_dilep6);

TCut inclusivenonzhighmass_dilep1("pt1 > 15. && pt2 > 15.");
TCut inclusivenonzhighmass_dilep2("abs(eormu1) == abs(eormu2)");
TCut inclusivenonzhighmass_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut inclusivenonzhighmass_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut inclusivenonzhighmass_dilep5("(((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 < 0) && abs(mass-91.) > 15.) || ((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 > 0) || (hyp_type == 1 || hyp_type == 2)");
TCut inclusivenonzhighmass_dilep6("mass > 120.");
TCut inclusivenonzhighmass_dilep7("ntrks > 2"); 
TCut inclusivenonzhighmass_dilep("inclusivenonzhighmass_dilep",inclusivenonzhighmass_dilep1+inclusivenonzhighmass_dilep2+inclusivenonzhighmass_dilep3+inclusivenonzhighmass_dilep4+inclusivenonzhighmass_dilep5+inclusivenonzhighmass_dilep6+inclusivenonzhighmass_dilep7);

TCut inclusivenonzhighmassv2_dilep1("pt1 > 20. && pt2 > 20.");
TCut inclusivenonzhighmassv2_dilep2("abs(eormu1) == abs(eormu2)");
TCut inclusivenonzhighmassv2_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut inclusivenonzhighmassv2_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut inclusivenonzhighmassv2_dilep5("(((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 < 0) && abs(mass-91.) > 15.) || ((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 > 0) || (hyp_type == 1 || hyp_type == 2)");
TCut inclusivenonzhighmassv2_dilep6("mass > 150.");
TCut inclusivenonzhighmassv2_dilep7("ntrks > 2"); 
TCut inclusivenonzhighmassv2_dilep("inclusivenonzhighmassv2_dilep",inclusivenonzhighmassv2_dilep1+inclusivenonzhighmassv2_dilep2+inclusivenonzhighmassv2_dilep3+inclusivenonzhighmassv2_dilep4+inclusivenonzhighmassv2_dilep5+inclusivenonzhighmassv2_dilep6+inclusivenonzhighmassv2_dilep7);

TCut inclusivenonzjets_dilep1("pt1 > 15. && pt2 > 15.");
TCut inclusivenonzjets_dilep2("abs(eormu1) == abs(eormu2)");
TCut inclusivenonzjets_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut inclusivenonzjets_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut inclusivenonzjets_dilep5("(((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 < 0) && abs(mass-91.) > 15.) || ((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 > 0) || (hyp_type == 1 || hyp_type == 2)");
TCut inclusivenonzjets_dilep6("njets>=2");
TCut inclusivenonzjets_dilep7("ntrks > 2"); 
TCut inclusivenonzjets_dilep("inclusivenonzjets_dilep",inclusivenonzjets_dilep1+inclusivenonzjets_dilep2+inclusivenonzjets_dilep3+inclusivenonzjets_dilep4+inclusivenonzjets_dilep5+inclusivenonzjets_dilep6+inclusivenonzjets_dilep7);

TCut inclusivenonzjetsv2_dilep1("pt1 > 20. && pt2 > 20.");
TCut inclusivenonzjetsv2_dilep2("abs(eormu1) == abs(eormu2)");
TCut inclusivenonzjetsv2_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut inclusivenonzjetsv2_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut inclusivenonzjetsv2_dilep5("(((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 < 0) && abs(mass-91.) > 15.) || ((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 > 0) || (hyp_type == 1 || hyp_type == 2)");
TCut inclusivenonzjetsv2_dilep6("njets>=2");
TCut inclusivenonzjetsv2_dilep7("pfmeff > 300. || tcmeff > 300.");
TCut inclusivenonzjetsv2_dilep8("ntrks > 2"); 
TCut inclusivenonzjetsv2_dilep("inclusivenonzjetsv2_dilep",inclusivenonzjetsv2_dilep1+inclusivenonzjetsv2_dilep2+inclusivenonzjetsv2_dilep3+inclusivenonzjetsv2_dilep4+inclusivenonzjetsv2_dilep5+inclusivenonzjetsv2_dilep6+inclusivenonzjetsv2_dilep7+inclusivenonzjetsv2_dilep8);

TCut inclusivenonzmet20_dilep1("pt1 > 15. && pt2 > 15.");
TCut inclusivenonzmet20_dilep2("abs(eormu1) == abs(eormu2)");
TCut inclusivenonzmet20_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut inclusivenonzmet20_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut inclusivenonzmet20_dilep5("(((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 < 0) && abs(mass-91.) > 15.) || ((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 > 0) || (hyp_type == 1 || hyp_type == 2)");
TCut inclusivenonzmet20_dilep6("tcmet>20 || pfmet>20");
TCut inclusivenonzmet20_dilep7("ntrks > 2"); 
TCut inclusivenonzmet20_dilep("inclusivenonzmet20_dilep",inclusivenonzmet20_dilep1+inclusivenonzmet20_dilep2+inclusivenonzmet20_dilep3+inclusivenonzmet20_dilep4+inclusivenonzmet20_dilep5+inclusivenonzmet20_dilep6+inclusivenonzmet20_dilep7);

TCut inclusivenonzmet20nojets_dilep1("pt1 > 15. && pt2 > 15.");
TCut inclusivenonzmet20nojets_dilep2("abs(eormu1) == abs(eormu2)");
TCut inclusivenonzmet20nojets_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut inclusivenonzmet20nojets_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut inclusivenonzmet20nojets_dilep5("(((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 < 0) && abs(mass-91.) > 15.) || ((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 > 0) || (hyp_type == 1 || hyp_type == 2)");
TCut inclusivenonzmet20nojets_dilep6("tcmet>20 || pfmet>20");
TCut inclusivenonzmet20nojets_dilep7("njets == 0");
TCut inclusivenonzmet20nojets_dilep8("ntrks > 2"); 
TCut inclusivenonzmet20nojets_dilep("inclusivenonzmet20nojets_dilep",inclusivenonzmet20nojets_dilep1+inclusivenonzmet20nojets_dilep2+inclusivenonzmet20nojets_dilep3+inclusivenonzmet20nojets_dilep4+inclusivenonzmet20nojets_dilep5+inclusivenonzmet20nojets_dilep6+inclusivenonzmet20nojets_dilep7+inclusivenonzmet20nojets_dilep8);

TCut inclusivenonzmet40jets_dilep1("pt1 > 15. && pt2 > 15.");
TCut inclusivenonzmet40jets_dilep2("abs(eormu1) == abs(eormu2)");
TCut inclusivenonzmet40jets_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut inclusivenonzmet40jets_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut inclusivenonzmet40jets_dilep5("(((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 < 0) && abs(mass-91.) > 15.) || ((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 > 0) || (hyp_type == 1 || hyp_type == 2)");
TCut inclusivenonzmet40jets_dilep6("tcmet>40 || pfmet>40");
TCut inclusivenonzmet40jets_dilep7("njets > 1");
TCut inclusivenonzmet40jets_dilep8("ntrks > 2");	
TCut inclusivenonzmet40jets_dilep("inclusivenonzmet40jets_dilep",inclusivenonzmet40jets_dilep1+inclusivenonzmet40jets_dilep2+inclusivenonzmet40jets_dilep3+inclusivenonzmet40jets_dilep4+inclusivenonzmet40jets_dilep5+inclusivenonzmet40jets_dilep6+inclusivenonzmet40jets_dilep7+inclusivenonzmet40jets_dilep8);

TCut inclusivenonzptz50_dilep1("pt1 > 15. && pt2 > 15.");
TCut inclusivenonzptz50_dilep2("abs(eormu1) == abs(eormu2)");
TCut inclusivenonzptz50_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut inclusivenonzptz50_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut inclusivenonzptz50_dilep5("(((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 < 0) && abs(mass-91.) > 15.) || ((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 > 0) || (hyp_type == 1 || hyp_type == 2)");
TCut inclusivenonzptz50_dilep6("dilpt > 50.");
TCut inclusivenonzptz50_dilep7("ntrks > 2"); 
TCut inclusivenonzptz50_dilep("inclusivenonzptz50_dilep",inclusivenonzptz50_dilep1+inclusivenonzptz50_dilep2+inclusivenonzptz50_dilep3+inclusivenonzptz50_dilep4+inclusivenonzptz50_dilep5+inclusivenonzptz50_dilep6+inclusivenonzptz50_dilep7);

TCut inclusivenonzptz100_dilep1("pt1 > 15. && pt2 > 15.");
TCut inclusivenonzptz100_dilep2("abs(eormu1) == abs(eormu2)");
TCut inclusivenonzptz100_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut inclusivenonzptz100_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut inclusivenonzptz100_dilep5("(((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 < 0) && abs(mass-91.) > 15.) || ((hyp_type == 0 || hyp_type == 3) && eormu1*eormu2 > 0) || (hyp_type == 1 || hyp_type == 2)");
TCut inclusivenonzptz100_dilep6("dilpt > 100.");
TCut inclusivenonzptz100_dilep7("ntrks > 2"); 
TCut inclusivenonzptz100_dilep("inclusivenonzptz100_dilep",inclusivenonzptz100_dilep1+inclusivenonzptz100_dilep2+inclusivenonzptz100_dilep3+inclusivenonzptz100_dilep4+inclusivenonzptz100_dilep5+inclusivenonzptz100_dilep6+inclusivenonzptz100_dilep7);

TCut dileptonictop_dilep1("pt1 > 15. && pt2 > 15.");
TCut dileptonictop_dilep2("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut dileptonictop_dilep3("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut dileptonictop_dilep4("((hyp_type == 0 || hyp_type == 3) && abs(mass-91.) > 15.) || (hyp_type == 1 || hyp_type == 2)");
TCut dileptonictop_dilep5("((hyp_type == 0 || hyp_type == 3) && (pfmet > 20. || tcmet > 20.)) || (hyp_type == 1 || hyp_type == 2)");
TCut dileptonictop_dilep6("njets >= 1");
TCut dileptonictop_dilep7("dileptonictop_dilep7","ntrks > 2");
TCut dileptonictop_dilep("dileptonictop_dilep",dileptonictop_dilep1+dileptonictop_dilep2+dileptonictop_dilep3+dileptonictop_dilep4+dileptonictop_dilep5+dileptonictop_dilep6+dileptonictop_dilep7);

TCut dileptonictopv2_dilep1("pt1 > 20. && pt2 > 20.");
TCut dileptonictopv2_dilep2("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut dileptonictopv2_dilep3("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut dileptonictopv2_dilep4("((hyp_type == 0 || hyp_type == 3) && abs(mass-91.) > 15.) || (hyp_type == 1 || hyp_type == 2)");
TCut dileptonictopv2_dilep5("((hyp_type == 0 || hyp_type == 3) && (pfmet > 20. || tcmet > 20.)) || (hyp_type == 1 || hyp_type == 2)");
TCut dileptonictopv2_dilep6("npurbtags > 1");
TCut dileptonictopv2_dilep7("pfmet > 30. && tcmet > 30.");
TCut dileptonictopv2_dilep8("ntrks > 2");
TCut dileptonictopv2_dilep("dileptonictopv2_dilep",dileptonictopv2_dilep1+dileptonictopv2_dilep2+dileptonictopv2_dilep3+dileptonictopv2_dilep4+dileptonictopv2_dilep5+dileptonictopv2_dilep6+dileptonictopv2_dilep7+dileptonictopv2_dilep8);

TCut dileptonictopv3_dilep1("pt1 > 20. && pt2 > 20.");
TCut dileptonictopv3_dilep2("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut dileptonictopv3_dilep3("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut dileptonictopv3_dilep4("((hyp_type == 0 || hyp_type == 3) && abs(mass-91.) > 15.) || (hyp_type == 1 || hyp_type == 2)");
TCut dileptonictopv3_dilep5("((hyp_type == 0 || hyp_type == 3) && (pfmet > 20. || tcmet > 20.)) || (hyp_type == 1 || hyp_type == 2)");
TCut dileptonictopv3_dilep6("npurbtags > 1");
TCut dileptonictopv3_dilep7("njets >= 2");
TCut dileptonictopv3_dilep8("pfmet > 30. && tcmet > 30.");
TCut dileptonictopv3_dilep9("ntrks > 2");
TCut dileptonictopv3_dilep("dileptonictopv3_dilep",dileptonictopv3_dilep1+dileptonictopv3_dilep2+dileptonictopv3_dilep3+dileptonictopv3_dilep4+dileptonictopv3_dilep5+dileptonictopv3_dilep6+dileptonictopv3_dilep7+dileptonictopv3_dilep8+dileptonictopv3_dilep9);

TCut dileptonictopv4_dilep1("pt1 > 20. && pt2 > 20.");
TCut dileptonictopv4_dilep2("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut dileptonictopv4_dilep3("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut dileptonictopv4_dilep4("((hyp_type == 0 || hyp_type == 3) && abs(mass-91.) > 15.) || (hyp_type == 1 || hyp_type == 2)");
TCut dileptonictopv4_dilep5("((hyp_type == 0 || hyp_type == 3) && (pfmet > 30. || tcmet > 30.)) || ((hyp_type == 1 || hyp_type == 2))");
TCut dileptonictopv4_dilep6("njets > 1");
TCut dileptonictopv4_dilep7("ntrks > 2");
TCut dileptonictopv4_dilep("dileptonictopv4_dilep",dileptonictopv4_dilep1+dileptonictopv4_dilep2+dileptonictopv4_dilep3+dileptonictopv4_dilep4+dileptonictopv4_dilep5+dileptonictopv4_dilep6+dileptonictopv4_dilep7);
TCut dileptonictopv4_dilep_mnjets("dileptonictopv4_dilep",dileptonictopv4_dilep1+dileptonictopv4_dilep2+dileptonictopv4_dilep3+dileptonictopv4_dilep4+dileptonictopv4_dilep5+dileptonictopv4_dilep7);

TCut samesigninclusive_dilep1("eormu1*eormu2 > 0");
TCut samesigninclusive_dilep2("pt1 > 15. && pt2 > 15.");
TCut samesigninclusive_dilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut samesigninclusive_dilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut samesigninclusive_dilep5("ntrks > 2");
TCut samesigninclusive_dilep("samesigninclusive_dilep",samesigninclusive_dilep1+samesigninclusive_dilep2+samesigninclusive_dilep3+samesigninclusive_dilep4+samesigninclusive_dilep5);

TCut zplusworzprime_dilep1("abs(eormu1) == abs(eormu2)");
TCut zplusworzprime_dilep2("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut zplusworzprime_dilep3("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut zplusworzprime_dilep4("mass > 76. && mass < 106.");
TCut zplusworzprime_dilep5("jetmass > 200.");
TCut zplusworzprime_dilep6("ntrks > 2"); 
TCut zplusworzprime_dilep("zplusworzprime_dilep",zplusworzprime_dilep1+zplusworzprime_dilep2+zplusworzprime_dilep3+zplusworzprime_dilep4+zplusworzprime_dilep5+zplusworzprime_dilep6);

TCut meff_dilep1("pt1 > 15. && pt2 > 15.");
TCut meff_dilep2("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut meff_dilep3("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut meff_dilep4("(tcmeff > 500. || pfmeff > 500.)");
TCut meff_dilep5("ntrks > 2");
TCut meff_dilep("meff_dilep",meff_dilep1+meff_dilep2+meff_dilep3+meff_dilep4+meff_dilep5);

//
// Trilepton+X
//

TCut inclusive_trilep1("max(max(pt1, pt2), pt3) > 10.");
TCut inclusive_trilep2("min(min(pt1, pt2), pt3) > 5.");
TCut inclusive_trilep3("(abs(eormu1) == 11 && e1_vbtf90full) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
TCut inclusive_trilep4("(abs(eormu2) == 11 && e2_vbtf90full) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
TCut inclusive_trilep5("(abs(eormu3) == 11 && e3_vbtf90full) || (abs(eormu3) == 13 && mu3_muonidfull && ! mu3_cosmic)");
TCut inclusive_trilep6("ntrks > 2");
TCut inclusive_trilep("inclusive_trilep",inclusive_trilep1+inclusive_trilep2+inclusive_trilep3+inclusive_trilep4+inclusive_trilep5+inclusive_trilep6);

#endif
