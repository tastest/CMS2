#include "globals.h"
//#include "cbpdfs.h"
//#include "cruijffpdfs.h"
#include "histpdfs.h"

#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"

#include "TCanvas.h"
#include "TH2F.h"
#include "TPaveText.h"

#include <iostream>

void jakes() {
    // in globals.C
    setVars();
    setDataSets();
    // in histpdfs.C
    setPDFs();

    float lambdazs[20] = {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
        -0.05,-0.1,-0.15,-0.2,-0.25,-0.3,-0.35,-0.4,-0.45,-0.5};

    //
    // signal+background pdfs
    //

    //RooRealVar* nsig = new RooRealVar("nsig","nsig",0.,10000.);
    //RooRealVar* nbkg = new RooRealVar("nbkg","nbkg",0.,10000.);
    //RooAddPdf*  sume = new RooAddPdf("sume","sume",RooArgList(*pdf_ww_pt1,*pdf_tt_pt1),RooArgList(*nsig,*nbkg));

    //RooRealVar* nsig_lz[20];
    //RooRealVar* nbkg_lz[20];
    //RooAddPdf*  sume_lz[20];
    //for(unsigned int i = 0; i < 20; ++i) {
    //    nsig_lz[i] = new RooRealVar(Form("nsig_lz%i",i),Form("nsig_lz%i",i),0.,10000.);
    //    nbkg_lz[i] = new RooRealVar(Form("nbkg_lz%i",i),Form("nbkg_lz%i",i),0.,10000.);
    //    sume_lz[i] = new RooAddPdf(Form("sume_lz%i",i),Form("sume_lz%i",i),RooArgList(*pdf_lz_pt1[i],*pdf_tt_pt1),RooArgList(*nsig_lz[i],*nbkg_lz[i]));
    //}

    //
    // ww vs. ttbar
    //

    TH2F* hSigmavNevts_deta  = new TH2F("hSigmavNevts_deta","hSigmavNevts_deta",11,-5.,105.,100,0.,100.);
    TH2F* hSigmavNevts_dilpt = new TH2F("hSigmavNevts_dilpt","hSigmavNevts_dilpt",11,-5.,105.,100,0.,100.);
    TH2F* hSigmavNevts_dphi  = new TH2F("hSigmavNevts_dphi","hSigmavNevts_dphi",11,-5.,105.,100,0.,100.);
    TH2F* hSigmavNevts_mass  = new TH2F("hSigmavNevts_mass","hSigmavNevts_mass",11,-5.,105.,100,0.,100.);
    TH2F* hSigmavNevts_met   = new TH2F("hSigmavNevts_met","hSigmavNevts_met",11,-5.,105.,100,0.,100.);
    TH2F* hSigmavNevts_pt1   = new TH2F("hSigmavNevts_pt1","hSigmavNevts_pt1",11,-5.,105.,100,0.,100.);
    TH2F* hSigmavNevts_pt2   = new TH2F("hSigmavNevts_pt2","hSigmavNevts_pt2",11,-5.,105.,100,0.,100.);
    TH2F* hSigmavNevts_pt1vdeta  = new TH2F("hSigmavNevts_pt1vdeta","hSigmavNevts_pt1vdeta",11,-5.,105.,100,0.,100.);
    TH2F* hSigmavNevts_pt1vdilpt = new TH2F("hSigmavNevts_pt1vdilpt","hSigmavNevts_pt1vdilpt",11,-5.,105.,100,0.,100.);
    TH2F* hSigmavNevts_pt1vdphi  = new TH2F("hSigmavNevts_pt1vdphi","hSigmavNevts_pt1vdphi",11,-5.,105.,100,0.,100.);
    TH2F* hSigmavNevts_pt1vmass  = new TH2F("hSigmavNevts_pt1vmass","hSigmavNevts_pt1vmass",11,-5.,105.,100,0.,100.);
    TH2F* hSigmavNevts_pt1vmet   = new TH2F("hSigmavNevts_pt1vmet","hSigmavNevts_pt1vmet",11,-5.,105.,100,0.,100.);
    TH2F* hSigmavNevts_pt1vpt2   = new TH2F("hSigmavNevts_pt1vpt2","hSigmavNevts_pt1vpt2",11,-5.,105.,100,0.,100.);
    RooDataSet* ww_red = 0;
    RooAbsReal *nll_ww = 0, *nll_tt = 0;
    for(int ii = 10; ii <= 100; ii+=10) {
        std::cout << "\nii=" << ii << " : ";
        std::cout.flush();

        for(int jj = 0; jj < 1000; ++jj) {
            if (!((jj+1)%100)) {
                std::cout << jj+1 << " ";
                std::cout.flush();
            }

            double nll_ww_val, nll_tt_val;

            ww_red = pdf_ww_deta->generate(*var_deta,ii);
            nll_ww = pdf_ww_deta->createNLL(*ww_red);
            nll_tt = pdf_tt_deta->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_deta->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;

            ww_red = pdf_ww_dilpt->generate(*var_dilpt,ii);
            nll_ww = pdf_ww_dilpt->createNLL(*ww_red);
            nll_tt = pdf_tt_dilpt->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_dilpt->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;

            ww_red = pdf_ww_dphi->generate(*var_dphi,ii);
            nll_ww = pdf_ww_dphi->createNLL(*ww_red);
            nll_tt = pdf_tt_dphi->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_dphi->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;

            ww_red = pdf_ww_mass->generate(*var_mass,ii);
            nll_ww = pdf_ww_mass->createNLL(*ww_red);
            nll_tt = pdf_tt_mass->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_mass->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;

            ww_red = pdf_ww_met->generate(*var_met,ii);
            nll_ww = pdf_ww_met->createNLL(*ww_red);
            nll_tt = pdf_tt_met->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_met->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;

            ww_red = pdf_ww_pt1->generate(*var_pt1,ii);
            nll_ww = pdf_ww_pt1->createNLL(*ww_red);
            nll_tt = pdf_tt_pt1->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_pt1->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;

            ww_red = pdf_ww_pt2->generate(*var_pt2,ii);
            nll_ww = pdf_ww_pt2->createNLL(*ww_red);
            nll_tt = pdf_tt_pt2->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_pt2->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;

            /*
            ww_red = pdf_ww_pt1vdeta->generate(RooArgSet(*var_pt1,*var_deta),ii);
            nll_ww = pdf_ww_pt1vdeta->createNLL(*ww_red);
            nll_tt = pdf_tt_pt1vdeta->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_pt1vdeta->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;

            ww_red = pdf_ww_pt1vdilpt->generate(RooArgSet(*var_pt1,*var_dilpt),ii);
            nll_ww = pdf_ww_pt1vdilpt->createNLL(*ww_red);
            nll_tt = pdf_tt_pt1vdilpt->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_pt1vdilpt->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;

            ww_red = pdf_ww_pt1vdphi->generate(RooArgSet(*var_pt1,*var_dphi),ii);
            nll_ww = pdf_ww_pt1vdphi->createNLL(*ww_red);
            nll_tt = pdf_tt_pt1vdphi->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_pt1vdphi->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;

            ww_red = pdf_ww_pt1vmass->generate(RooArgSet(*var_pt1,*var_mass),ii);
            nll_ww = pdf_ww_pt1vmass->createNLL(*ww_red);
            nll_tt = pdf_tt_pt1vmass->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_pt1vmass->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;

            ww_red = pdf_ww_pt1vmet->generate(RooArgSet(*var_pt1,*var_met),ii);
            nll_ww = pdf_ww_pt1vmet->createNLL(*ww_red);
            nll_tt = pdf_tt_pt1vmet->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_pt1vmet->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;

            ww_red = pdf_ww_pt1vpt2->generate(RooArgSet(*var_pt1,*var_pt2),ii);
            nll_ww = pdf_ww_pt1vpt2->createNLL(*ww_red);
            nll_tt = pdf_tt_pt1vpt2->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            nll_tt->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_tt_val = nll_tt->getVal();
            hSigmavNevts_pt1vpt2->Fill(ii,sqrt(2.*fabs(nll_ww_val-nll_tt_val)));

            delete nll_ww;
            delete nll_tt;
            delete ww_red;
            ww_red = 0;
            nll_ww = 0;
            nll_tt = 0;
            */
        }
    }

    //
    // ww + ttbar at a ratio of 4:1
    //

    TH2F* hDiffNSigvNevts = new TH2F("hDiffNSigvNevts","hDiffNSigvNevts",11,-5.,105.,100,-10.,10.);
    TH2F* hDiffNBkgvNevts = new TH2F("hDiffNBkgvNevts","hDiffNBkgvNevts",11,-5.,105.,100,-10.,10.);
    RooDataSet* ww_red = 0;
    RooDataSet* tt_red = 0;
    for(int ii = 10; ii <= 100; ii+=10) {
        std::cout << "\nii=" << ii << " : ";
        std::cout.flush();

        for(int jj = 0; jj < 1000; ++jj) {
            if (!((jj+1)%100)) {
                std::cout << jj+1 << " ";
                std::cout.flush();
            }

            ww_red = pdf_ww_pt1->generate(*var_pt1,ii);

            float tmp = (float)ii/4.;
            tt_red = pdf_tt_pt1->generate(*var_pt1,(int)tmp);

            ww_red->append(*tt_red);
            sume->fitTo(*ww_red);

            Double_t fitnsig  = nsig->getVal();
            Double_t fitnbkg  = nbkg->getVal();
            Double_t fitnsige = nsig->getError() ? nsig->getError() : 1E-6;
            Double_t fitnbkge = nbkg->getError() ? nbkg->getError() : 1E-6;

            float nsigpull = ((float)ii-fitnsig)/fitnsige;
            float nbkgpull = (tmp-fitnbkg)/fitnbkge;

            hDiffNSigvNevts->Fill(ii,nsigpull);
            hDiffNBkgvNevts->Fill(ii,nbkgpull);

            delete ww_red;
            delete tt_red;
            ww_red = 0;
            tt_red = 0;
        }
    }

    //
    // 10 -> ~30/pb, 30 -> ~100/pb, 100 -> ~300/pb
    //

    TH2F* hSigma10vLambdaz_deta   = new TH2F("hSigma10vLambdaz_deta","hSigma10vLambdaz_deta",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_deta   = new TH2F("hSigma30vLambdaz_deta","hSigma30vLambdaz_deta",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_deta  = new TH2F("hSigma100vLambdaz_deta","hSigma100vLambdaz_deta",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma10vLambdaz_dilpt  = new TH2F("hSigma10vLambdaz_dilpt","hSigma10vLambdaz_dilpt",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_dilpt  = new TH2F("hSigma30vLambdaz_dilpt","hSigma30vLambdaz_dilpt",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_dilpt = new TH2F("hSigma100vLambdaz_dilpt","hSigma100vLambdaz_dilpt",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma10vLambdaz_dphi   = new TH2F("hSigma10vLambdaz_dphi","hSigma10vLambdaz_dphi",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_dphi   = new TH2F("hSigma30vLambdaz_dphi","hSigma30vLambdaz_dphi",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_dphi  = new TH2F("hSigma100vLambdaz_dphi","hSigma100vLambdaz_dphi",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma10vLambdaz_mass   = new TH2F("hSigma10vLambdaz_mass","hSigma10vLambdaz_mass",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_mass   = new TH2F("hSigma30vLambdaz_mass","hSigma30vLambdaz_mass",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_mass  = new TH2F("hSigma100vLambdaz_mass","hSigma100vLambdaz_mass",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma10vLambdaz_met    = new TH2F("hSigma10vLambdaz_met","hSigma10vLambdaz_met",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_met    = new TH2F("hSigma30vLambdaz_met","hSigma30vLambdaz_met",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_met   = new TH2F("hSigma100vLambdaz_met","hSigma100vLambdaz_met",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma10vLambdaz_pt1    = new TH2F("hSigma10vLambdaz_pt1","hSigma10vLambdaz_pt1",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_pt1    = new TH2F("hSigma30vLambdaz_pt1","hSigma30vLambdaz_pt1",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_pt1   = new TH2F("hSigma100vLambdaz_pt1","hSigma100vLambdaz_pt1",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma10vLambdaz_pt2    = new TH2F("hSigma10vLambdaz_pt2","hSigma10vLambdaz_pt2",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_pt2    = new TH2F("hSigma30vLambdaz_pt2","hSigma30vLambdaz_pt2",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_pt2   = new TH2F("hSigma100vLambdaz_pt2","hSigma100vLambdaz_pt2",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma10vLambdaz_pt1vdeta   = new TH2F("hSigma10vLambdaz_pt1vdeta","hSigma10vLambdaz_pt1vdeta",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_pt1vdeta   = new TH2F("hSigma30vLambdaz_pt1vdeta","hSigma30vLambdaz_pt1vdeta",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_pt1vdeta  = new TH2F("hSigma100vLambdaz_pt1vdeta","hSigma100vLambdaz_pt1vdeta",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma10vLambdaz_pt1vdilpt  = new TH2F("hSigma10vLambdaz_pt1vdilpt","hSigma10vLambdaz_pt1vdilpt",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_pt1vdilpt  = new TH2F("hSigma30vLambdaz_pt1vdilpt","hSigma30vLambdaz_pt1vdilpt",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_pt1vdilpt = new TH2F("hSigma100vLambdaz_pt1vdilpt","hSigma100vLambdaz_pt1vdilpt",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma10vLambdaz_pt1vdphi   = new TH2F("hSigma10vLambdaz_pt1vdphi","hSigma10vLambdaz_pt1vdphi",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_pt1vdphi   = new TH2F("hSigma30vLambdaz_pt1vdphi","hSigma30vLambdaz_pt1vdphi",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_pt1vdphi  = new TH2F("hSigma100vLambdaz_pt1vdphi","hSigma100vLambdaz_pt1vdphi",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma10vLambdaz_pt1vmass   = new TH2F("hSigma10vLambdaz_pt1vmass","hSigma10vLambdaz_pt1vmass",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_pt1vmass   = new TH2F("hSigma30vLambdaz_pt1vmass","hSigma30vLambdaz_pt1vmass",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_pt1vmass  = new TH2F("hSigma100vLambdaz_pt1vmass","hSigma100vLambdaz_pt1vmass",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma10vLambdaz_pt1vmet    = new TH2F("hSigma10vLambdaz_pt1vmet","hSigma10vLambdaz_pt1vmet",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_pt1vmet    = new TH2F("hSigma30vLambdaz_pt1vmet","hSigma30vLambdaz_pt1vmet",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_pt1vmet   = new TH2F("hSigma100vLambdaz_pt1vmet","hSigma100vLambdaz_pt1vmet",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma10vLambdaz_pt1vpt2    = new TH2F("hSigma10vLambdaz_pt1vpt2","hSigma10vLambdaz_pt1vpt2",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz_pt1vpt2    = new TH2F("hSigma30vLambdaz_pt1vpt2","hSigma30vLambdaz_pt1vpt2",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz_pt1vpt2   = new TH2F("hSigma100vLambdaz_pt1vpt2","hSigma100vLambdaz_pt1vpt2",21,-0.505,0.505,100,0.,100.);
    RooDataSet* ww10      = 0;
    RooDataSet* ww30      = 0;
    RooDataSet* ww100     = 0;
    RooAbsReal *nll_ww10  = 0, *nll_lz10  = 0;
    RooAbsReal *nll_ww30  = 0, *nll_lz30  = 0;
    RooAbsReal *nll_ww100 = 0, *nll_lz100 = 0;
    for(int ii = 0; ii < 20; ++ii) { // lambdazs
        std::cout << "\nii=" << ii << " : ";
        std::cout.flush();

        for(int jj = 0; jj < 1000; ++jj) {
            if (!((jj+1)%100)) {
                std::cout << jj+1 << " ";
                std::cout.flush();
            }

            double nll_ww_val, nll_lz_val;

            ww10 = pdf_ww_deta->generate(*var_deta,10);
            nll_ww10 = pdf_ww_deta->createNLL(*ww10);
            nll_lz10 = pdf_lz_deta[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_deta->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_deta->generate(*var_deta,30);
            nll_ww30 = pdf_ww_deta->createNLL(*ww30);
            nll_lz30 = pdf_lz_deta[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_deta->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_deta->generate(*var_deta,100);
            nll_ww100 = pdf_ww_deta->createNLL(*ww100);
            nll_lz100 = pdf_lz_deta[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_deta->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;

            ww10 = pdf_ww_dilpt->generate(*var_dilpt,10);
            nll_ww10 = pdf_ww_dilpt->createNLL(*ww10);
            nll_lz10 = pdf_lz_dilpt[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_dilpt->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_dilpt->generate(*var_dilpt,30);
            nll_ww30 = pdf_ww_dilpt->createNLL(*ww30);
            nll_lz30 = pdf_lz_dilpt[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_dilpt->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_dilpt->generate(*var_dilpt,100);
            nll_ww100 = pdf_ww_dilpt->createNLL(*ww100);
            nll_lz100 = pdf_lz_dilpt[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_dilpt->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;

            ww10 = pdf_ww_dphi->generate(*var_dphi,10);
            nll_ww10 = pdf_ww_dphi->createNLL(*ww10);
            nll_lz10 = pdf_lz_dphi[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_dphi->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_dphi->generate(*var_dphi,30);
            nll_ww30 = pdf_ww_dphi->createNLL(*ww30);
            nll_lz30 = pdf_lz_dphi[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_dphi->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_dphi->generate(*var_dphi,100);
            nll_ww100 = pdf_ww_dphi->createNLL(*ww100);
            nll_lz100 = pdf_lz_dphi[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_dphi->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;

            ww10 = pdf_ww_mass->generate(*var_mass,10);
            nll_ww10 = pdf_ww_mass->createNLL(*ww10);
            nll_lz10 = pdf_lz_mass[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_mass->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_mass->generate(*var_mass,30);
            nll_ww30 = pdf_ww_mass->createNLL(*ww30);
            nll_lz30 = pdf_lz_mass[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_mass->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_mass->generate(*var_mass,100);
            nll_ww100 = pdf_ww_mass->createNLL(*ww100);
            nll_lz100 = pdf_lz_mass[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_mass->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;

            ww10 = pdf_ww_met->generate(*var_met,10);
            nll_ww10 = pdf_ww_met->createNLL(*ww10);
            nll_lz10 = pdf_lz_met[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_met->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_met->generate(*var_met,30);
            nll_ww30 = pdf_ww_met->createNLL(*ww30);
            nll_lz30 = pdf_lz_met[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_met->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_met->generate(*var_met,100);
            nll_ww100 = pdf_ww_met->createNLL(*ww100);
            nll_lz100 = pdf_lz_met[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_met->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;

            ww10 = pdf_ww_pt1->generate(*var_pt1,10);
            nll_ww10 = pdf_ww_pt1->createNLL(*ww10);
            nll_lz10 = pdf_lz_pt1[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_pt1->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_pt1->generate(*var_pt1,30);
            nll_ww30 = pdf_ww_pt1->createNLL(*ww30);
            nll_lz30 = pdf_lz_pt1[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_pt1->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_pt1->generate(*var_pt1,100);
            nll_ww100 = pdf_ww_pt1->createNLL(*ww100);
            nll_lz100 = pdf_lz_pt1[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_pt1->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;

            ww10 = pdf_ww_pt2->generate(*var_pt2,10);
            nll_ww10 = pdf_ww_pt2->createNLL(*ww10);
            nll_lz10 = pdf_lz_pt2[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_pt2->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_pt2->generate(*var_pt2,30);
            nll_ww30 = pdf_ww_pt2->createNLL(*ww30);
            nll_lz30 = pdf_lz_pt2[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_pt2->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_pt2->generate(*var_pt2,100);
            nll_ww100 = pdf_ww_pt2->createNLL(*ww100);
            nll_lz100 = pdf_lz_pt2[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_pt2->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;

            /*
            ww10 = pdf_ww_pt1vdeta->generate(RooArgSet(*var_pt1,*var_deta),10);
            nll_ww10 = pdf_ww_pt1vdeta->createNLL(*ww10);
            nll_lz10 = pdf_lz_pt1vdeta[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_pt1vdeta->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_pt1vdeta->generate(RooArgSet(*var_pt1,*var_deta),30);
            nll_ww30 = pdf_ww_pt1vdeta->createNLL(*ww30);
            nll_lz30 = pdf_lz_pt1vdeta[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_pt1vdeta->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_pt1vdeta->generate(RooArgSet(*var_pt1,*var_deta),100);
            nll_ww100 = pdf_ww_pt1vdeta->createNLL(*ww100);
            nll_lz100 = pdf_lz_pt1vdeta[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_pt1vdeta->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;

            ww10 = pdf_ww_pt1vdilpt->generate(RooArgSet(*var_pt1,*var_dilpt),10);
            nll_ww10 = pdf_ww_pt1vdilpt->createNLL(*ww10);
            nll_lz10 = pdf_lz_pt1vdilpt[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_pt1vdilpt->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_pt1vdilpt->generate(RooArgSet(*var_pt1,*var_dilpt),30);
            nll_ww30 = pdf_ww_pt1vdilpt->createNLL(*ww30);
            nll_lz30 = pdf_lz_pt1vdilpt[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_pt1vdilpt->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_pt1vdilpt->generate(RooArgSet(*var_pt1,*var_dilpt),100);
            nll_ww100 = pdf_ww_pt1vdilpt->createNLL(*ww100);
            nll_lz100 = pdf_lz_pt1vdilpt[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_pt1vdilpt->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;

            ww10 = pdf_ww_pt1vdphi->generate(RooArgSet(*var_pt1,*var_dphi),10);
            nll_ww10 = pdf_ww_pt1vdphi->createNLL(*ww10);
            nll_lz10 = pdf_lz_pt1vdphi[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_pt1vdphi->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_pt1vdphi->generate(RooArgSet(*var_pt1,*var_dphi),30);
            nll_ww30 = pdf_ww_pt1vdphi->createNLL(*ww30);
            nll_lz30 = pdf_lz_pt1vdphi[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_pt1vdphi->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_pt1vdphi->generate(RooArgSet(*var_pt1,*var_dphi),100);
            nll_ww100 = pdf_ww_pt1vdphi->createNLL(*ww100);
            nll_lz100 = pdf_lz_pt1vdphi[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_pt1vdphi->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;

            ww10 = pdf_ww_pt1vmass->generate(RooArgSet(*var_pt1,*var_mass),10);
            nll_ww10 = pdf_ww_pt1vmass->createNLL(*ww10);
            nll_lz10 = pdf_lz_pt1vmass[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_pt1vmass->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_pt1vmass->generate(RooArgSet(*var_pt1,*var_mass),30);
            nll_ww30 = pdf_ww_pt1vmass->createNLL(*ww30);
            nll_lz30 = pdf_lz_pt1vmass[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_pt1vmass->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_pt1vmass->generate(RooArgSet(*var_pt1,*var_mass),100);
            nll_ww100 = pdf_ww_pt1vmass->createNLL(*ww100);
            nll_lz100 = pdf_lz_pt1vmass[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_pt1vmass->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;

            ww10 = pdf_ww_pt1vmet->generate(RooArgSet(*var_pt1,*var_met),10);
            nll_ww10 = pdf_ww_pt1vmet->createNLL(*ww10);
            nll_lz10 = pdf_lz_pt1vmet[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_pt1vmet->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_pt1vmet->generate(RooArgSet(*var_pt1,*var_met),30);
            nll_ww30 = pdf_ww_pt1vmet->createNLL(*ww30);
            nll_lz30 = pdf_lz_pt1vmet[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_pt1vmet->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_pt1vmet->generate(RooArgSet(*var_pt1,*var_met),100);
            nll_ww100 = pdf_ww_pt1vmet->createNLL(*ww100);
            nll_lz100 = pdf_lz_pt1vmet[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_pt1vmet->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;

            ww10 = pdf_ww_pt1vpt2->generate(RooArgSet(*var_pt1,*var_pt2),10);
            nll_ww10 = pdf_ww_pt1vpt2->createNLL(*ww10);
            nll_lz10 = pdf_lz_pt1vpt2[ii]->createNLL(*ww10);
            nll_ww10->addServer(*var_dummy);
            nll_lz10->addServer(*var_dummy);
            nll_ww_val = nll_ww10->getVal();
            nll_lz_val = nll_lz10->getVal();
            hSigma10vLambdaz_pt1vpt2->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_pt1vpt2->generate(RooArgSet(*var_pt1,*var_pt2),30);
            nll_ww30 = pdf_ww_pt1vpt2->createNLL(*ww30);
            nll_lz30 = pdf_lz_pt1vpt2[ii]->createNLL(*ww30);
            nll_ww30->addServer(*var_dummy);
            nll_lz30->addServer(*var_dummy);
            nll_ww_val = nll_ww30->getVal();
            nll_lz_val = nll_lz30->getVal();
            hSigma30vLambdaz_pt1vpt2->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_pt1vpt2->generate(RooArgSet(*var_pt1,*var_pt2),100);
            nll_ww100 = pdf_ww_pt1vpt2->createNLL(*ww100);
            nll_lz100 = pdf_lz_pt1vpt2[ii]->createNLL(*ww100);
            nll_ww100->addServer(*var_dummy);
            nll_lz100->addServer(*var_dummy);
            nll_ww_val = nll_ww100->getVal();
            nll_lz_val = nll_lz100->getVal();
            hSigma100vLambdaz_pt1vpt2->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete nll_ww10;
            delete nll_lz10;
            delete nll_ww30;
            delete nll_lz30;
            delete nll_ww100;
            delete nll_lz100;
            delete ww10;
            delete ww30;
            delete ww100;
            nll_ww10  = 0;
            nll_lz10  = 0;
            nll_ww30  = 0;
            nll_lz30  = 0;
            nll_ww100 = 0;
            nll_lz100 = 0;
            ww10      = 0;
            ww30      = 0;
            ww100     = 0;
            */
        }
    }

    //TCanvas *cww = new TCanvas("cww","cww");
    //RooPlot *plotww = var_pt1->frame(RooFit::Bins(40),RooFit::Title("ww leading lepton pT"));
    //ds_ww_pt1 ->plotOn(plotww);
    //ds_ww_pt1 ->statOn(plotww);
    //pdf_ww_pt1->plotOn(plotww);
    //plotww->Draw();
    //int   nfloatpars = 0;
    //int   ndf = 40-nfloatpars;
    //float chi2ndf = plotww->chiSquare(nfloatpars);
    //float chi2 = chi2ndf*ndf;
    //TPaveText *pav = new TPaveText(0.65,0.65,0.99,0.75,"NDC");
    //pav->AddText(Form("chi2/ndf = %.2f/%i",chi2,ndf));
    //plotww->addObject(pav);
    //plotww->Draw();
    //cww->Draw();
}
