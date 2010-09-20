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

void jakes() {
    // in globals.C
    setDataSets();
    // in histpdfs.C
    setPDFs();

    float lambdazs[20] = {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
        -0.05,-0.1,-0.15,-0.2,-0.25,-0.3,-0.35,-0.4,-0.45,-0.5};

    //
    // signal+background pdfs
    //
    RooRealVar* nsig = new RooRealVar("nsig","nsig",0.,10000.);
    RooRealVar* nbkg = new RooRealVar("nbkg","nbkg",0.,10000.);
    RooAddPdf*  sume = new RooAddPdf("sume","sume",RooArgList(*pdf_ww_pt,*pdf_tt_pt),RooArgList(*nsig,*nbkg));

    //RooRealVar* nsig_lz[20];
    //RooRealVar* nbkg_lz[20];
    //RooAddPdf*  sume_lz[20];
    //for(unsigned int i = 0; i < 20; ++i) {
    //    nsig_lz[i] = new RooRealVar(Form("nsig_lz%i",i),Form("nsig_lz%i",i),0.,10000.);
    //    nbkg_lz[i] = new RooRealVar(Form("nbkg_lz%i",i),Form("nbkg_lz%i",i),0.,10000.);
    //    sume_lz[i] = new RooAddPdf(Form("sume_lz%i",i),Form("sume_lz%i",i),RooArgList(*pdf_lz_pt[i],*pdf_tt_pt),RooArgList(*nsig_lz[i],*nbkg_lz[i]));
    //}

    /*
    //
    // ww vs. ttbar
    //
    TH2F* hSigmavNevts = new TH2F("hSigmavNevts","hSigmavNevts",11,-5.,105.,100,0.,10.);
    RooDataSet* ww_red = 0;
    for(int ii = 10; ii <= 100; ii+=10) {
        for(int jj = 0; jj < 1000; ++jj)
        {
            ww_red = pdf_ww_pt->generate(*var_pt,ii);

            RooAbsReal* nll_ww = pdf_ww_pt->createNLL(*ww_red);
            nll_ww->addServer(*var_dummy);
            RooAbsReal* nll_tt = pdf_tt_pt->createNLL(*ww_red);
            nll_tt->addServer(*var_dummy);

            double nll_ww_val = nll_ww->getVal();
            double nll_tt_val = nll_tt->getVal();

            hSigmavNevts->Fill(ii,sqrt(2*fabs(nll_ww_val-nll_tt_val)));

            delete ww_red;
            ww_red = 0;
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
        for(int jj = 0; jj < 1000; ++jj)
        {
            ww_red = pdf_ww_pt->generate(*var_pt,ii);

            float tmp = (float)ii/4.;
            tt_red = pdf_tt_pt->generate(*var_pt,(int)tmp);

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
    TH2F* hSigma10vLambdaz  = new TH2F("hSigma10vLambdaz","hSigma10vLambdaz",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma30vLambdaz  = new TH2F("hSigma30vLambdaz","hSigma30vLambdaz",21,-0.505,0.505,100,0.,100.);
    TH2F* hSigma100vLambdaz = new TH2F("hSigma100vLambdaz","hSigma100vLambdaz",21,-0.505,0.505,100,0.,100.);
    RooDataSet* ww10  = 0;
    RooDataSet* ww30  = 0;
    RooDataSet* ww100 = 0;
    for(int ii = 0; ii < 20; ++ii) { // lambdazs
        for(int jj = 0; jj < 1000; ++jj)
        {
            RooAbsReal *nll_ww, *nll_lz;
            double nll_ww_val, nll_lz_val;

            ww10 = pdf_ww_pt->generate(*var_pt,10);
            nll_ww = pdf_ww_pt->createNLL(*ww10);
            nll_lz = pdf_lz_pt[ii]->createNLL(*ww10);
            nll_ww->addServer(*var_dummy);
            nll_lz->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_lz_val = nll_lz->getVal();
            hSigma10vLambdaz->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww30 = pdf_ww_pt->generate(*var_pt,30);
            nll_ww = pdf_ww_pt->createNLL(*ww30);
            nll_lz = pdf_lz_pt[ii]->createNLL(*ww30);
            nll_ww->addServer(*var_dummy);
            nll_lz->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_lz_val = nll_lz->getVal();
            hSigma30vLambdaz->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            ww100 = pdf_ww_pt->generate(*var_pt,100);
            nll_ww = pdf_ww_pt->createNLL(*ww100);
            nll_lz = pdf_lz_pt[ii]->createNLL(*ww100);
            nll_ww->addServer(*var_dummy);
            nll_lz->addServer(*var_dummy);
            nll_ww_val = nll_ww->getVal();
            nll_lz_val = nll_lz->getVal();
            hSigma100vLambdaz->Fill(lambdazs[ii],sqrt(2.*fabs(nll_ww_val-nll_lz_val)));

            delete ww10;
            delete ww30;
            delete ww100;
            ww10  = 0;
            ww30  = 0;
            ww100 = 0;
        }
    }
    */

    //TCanvas *cww = new TCanvas("cww","cww");
    //RooPlot *plotww = var_pt->frame(RooFit::Bins(40),RooFit::Title("ww leading lepton pT"));
    //ds_ww_pt ->plotOn(plotww);
    //ds_ww_pt ->statOn(plotww);
    //pdf_ww_pt->plotOn(plotww);
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
