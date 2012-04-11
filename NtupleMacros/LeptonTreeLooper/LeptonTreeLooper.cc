
#include "LeptonTreeLooper.h"
#include "../../../Smurf/Core/LeptonTree.h"
#include "../Tools/goodrun.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TChain.h"

#include <algorithm>
#include <utility>

LeptonTreeLooper::LeptonTreeLooper()
{
    runlistIsSet_ = false;
}

LeptonTreeLooper::~LeptonTreeLooper()
{
}

void LeptonTreeLooper::setGoodRunList(const char *runlist)
{
    set_goodrun_file(runlist);
    runlistIsSet_ = true;
}

void LeptonTreeLooper::loop(TChain *chain, TString name)
{

    printf("[LeptonTreeLooper::loop] %s\n", name.Data());

    //
    // check for valid chain
    //

    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    if (listOfFiles->GetEntries() == 0) {
        std::cout << "[LeptonTreeLooper::loop] no files in chain" << std::endl;
        return;
    }

    //
    // load nvtx weight
    //

    //TFile f_weight("nvtxweight.root", "READ");
    //gROOT->cd();
    //TH1F *h1_nvtxweight = (TH1F*)f_weight.Get("weight")->Clone("weight");
    //f_weight.Close();

    //
    // set up N-1 cuts
    //

enum CutType {
    DETAIN          = (1<<0),
    DPHIIN          = (1<<1),
    SIGMAIETAIETA   = (1<<2),
    HOE             = (1<<3),
    OOEMOOP         = (1<<4),
    D0VTX           = (1<<5),
    DZVTX           = (1<<6),
    ISO             = (1<<7),
    VTXFIT          = (1<<8),
    MHITS           = (1<<9),
    EOPFBREM        = (1<<10)
};

    unsigned int all = DETAIN | DPHIIN | SIGMAIETAIETA | HOE | OOEMOOP | D0VTX | DZVTX | ISO | VTXFIT | MHITS;
    unsigned int nm1_detain     = all & ~DETAIN;
    unsigned int nm1_dphiin     = all & ~DPHIIN;
    unsigned int nm1_sieie      = all & ~SIGMAIETAIETA;
    unsigned int nm1_hoe        = all & ~HOE;
    unsigned int nm1_ooemoop    = all & ~OOEMOOP;
    unsigned int nm1_d0vtx      = all & ~D0VTX;
    unsigned int nm1_dzvtx      = all & ~DZVTX;
    unsigned int nm1_iso        = all & ~ISO;
    unsigned int nm1_vfitprob   = all & ~VTXFIT;
    unsigned int nm1_mhit       = all & ~MHITS;

    //
    // set up histograms
    //

    gROOT->cd();

    // tag and probe
    TH1F *h1_tp_nvtx = new TH1F(Form("%s_h1_tp_nvtx", name.Data()), "nvtx", 40, -0.5, 39.5);
    TH1F *h1_tp_tagAndProbeMass = new TH1F(Form("%s_h1_tp_tagAndProbeMass", name.Data()), "tagAndProbeMass", 80, 20.0, 160.0);
    h1_tp_nvtx->Sumw2();
    h1_tp_tagAndProbeMass->Sumw2();

    TH1F *h1_tp_detain = new TH1F(Form("%s_h1_tp_detain", name.Data()), "detain", 50, 0.0, 0.01);
    TH1F *h1_tp_dphiin = new TH1F(Form("%s_h1_tp_dphiin", name.Data()), "dphiin", 50, 0.0, 0.1);
    TH1F *h1_tp_sieie = new TH1F(Form("%s_h1_tp_sieie", name.Data()), "sieie", 50, 0.0, 0.05);
    TH1F *h1_tp_hoe = new TH1F(Form("%s_h1_tp_hoe", name.Data()), "hoe", 50, 0.0, 0.5);
    TH1F *h1_tp_ooemoop = new TH1F(Form("%s_h1_tp_ooemoop", name.Data()), "ooemoop", 50, 0.0, 0.05);
    TH1F *h1_tp_d0vtx = new TH1F(Form("%s_h1_tp_d0vtx", name.Data()), "d0vtx", 50, 0.0, 0.1);
    TH1F *h1_tp_dzvtx = new TH1F(Form("%s_h1_tp_dzvtx", name.Data()), "dzvtx", 50, 0.0, 0.2);
    TH1F *h1_tp_vfitprob = new TH1F(Form("%s_h1_tp_vfitprob", name.Data()), "vfitprob", 2, -0.5, 1.5);
    TH1F *h1_tp_mhit = new TH1F(Form("%s_h1_tp_mhit", name.Data()), "mhit", 5, -0.5, 4.5);
    h1_tp_detain->Sumw2();
    h1_tp_dphiin->Sumw2();
    h1_tp_sieie->Sumw2();
    h1_tp_hoe->Sumw2();
    h1_tp_ooemoop->Sumw2();
    h1_tp_d0vtx->Sumw2();
    h1_tp_dzvtx->Sumw2();
    h1_tp_vfitprob->Sumw2();
    h1_tp_mhit->Sumw2();

    // fake rate
    TH1F *h1_fr_nvtx = new TH1F(Form("%s_h1_fr_nvtx", name.Data()), "nvtx", 40, -0.5, 39.5);
    h1_fr_nvtx->Sumw2();

    TH1F *h1_fr_detain = new TH1F(Form("%s_h1_fr_detain", name.Data()), "detain", 50, 0.0, 0.01);
    TH1F *h1_fr_dphiin = new TH1F(Form("%s_h1_fr_dphiin", name.Data()), "dphiin", 50, 0.0, 0.1);
    TH1F *h1_fr_sieie = new TH1F(Form("%s_h1_fr_sieie", name.Data()), "sieie", 50, 0.0, 0.05);
    TH1F *h1_fr_hoe = new TH1F(Form("%s_h1_fr_hoe", name.Data()), "hoe", 50, 0.0, 0.5);
    TH1F *h1_fr_ooemoop = new TH1F(Form("%s_h1_fr_ooemoop", name.Data()), "ooemoop", 50, 0.0, 0.05);
    TH1F *h1_fr_d0vtx = new TH1F(Form("%s_h1_fr_d0vtx", name.Data()), "d0vtx", 50, 0.0, 0.1);
    TH1F *h1_fr_dzvtx = new TH1F(Form("%s_h1_fr_dzvtx", name.Data()), "dzvtx", 50, 0.0, 0.2);
    TH1F *h1_fr_vfitprob = new TH1F(Form("%s_h1_fr_vfitprob", name.Data()), "vfitprob", 2, -0.5, 1.5);
    TH1F *h1_fr_mhit = new TH1F(Form("%s_h1_fr_mhit", name.Data()), "mhit", 5, -0.5, 4.5);
    h1_fr_detain->Sumw2();
    h1_fr_dphiin->Sumw2();
    h1_fr_sieie->Sumw2();
    h1_fr_hoe->Sumw2();
    h1_fr_ooemoop->Sumw2();
    h1_fr_d0vtx->Sumw2();
    h1_fr_dzvtx->Sumw2();
    h1_fr_vfitprob->Sumw2();
    h1_fr_mhit->Sumw2();

    //
    // file loop
    //

    unsigned int nEventsChain=0;
    unsigned int nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;
    int i_permille_old = 0;

    while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

        //
        // load the lepton tree
        //

        LeptonTree *tree = new LeptonTree();
        tree->LoadTree(currentFile->GetTitle());
        tree->InitTree();

        //
        // event loop
        //

        ULong64_t nEvents = tree->tree_->GetEntries();
        for(ULong64_t event = 0; event < nEvents; ++event) {
            tree->tree_->GetEntry(event);

            //
            // incrimenet counters
            //

            ++nEventsTotal;
            int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
            if (i_permille != i_permille_old) {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1)) {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                            "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
                    fflush(stdout);
                }
                i_permille_old = i_permille;
            }

            //
            // good run list
            //

            if (runlistIsSet_) {
                if (!goodrun_json(tree->run_, tree->lumi_)) continue;
            }

            //
            // nvtx weight
            //

            float weight = 1.0;
            //int nvtxbin = h1_nvtxweight->FindBin(TMath::Min(int(tree->nvtx_), 19));
            //if (name == "dyee") weight = h1_nvtxweight->GetBinContent(nvtxbin);

            //
            // dataset specific cuts
            //

            if (name == "data_2012A") {
                if (tree->run_ < 190659) continue;
            }

            //
            // fill tag and probe histograms
            //

            if ((tree->eventSelection_ & LeptonTree::ZeeTagAndProbe) && tree->probe_.Pt() > 10.0)
            {
                h1_tp_nvtx->Fill(tree->nvtx_, weight);
                h1_tp_tagAndProbeMass->Fill(tree->tagAndProbeMass_, weight);

                if (abs(tree->tagAndProbeMass_ - 91) < 15.0) {
                    h1_tp_nvtx->Fill(tree->nvtx_, weight);
                    if ((tree->looseId_ & nm1_detain) == nm1_detain)     h1_tp_detain->Fill(    fabs(tree->detain_));
                    if ((tree->looseId_ & nm1_dphiin) == nm1_dphiin)     h1_tp_dphiin->Fill(    fabs(tree->dphiin_));
                    if ((tree->looseId_ & nm1_sieie) == nm1_sieie)       h1_tp_sieie->Fill(     tree->sieie_);
                    if ((tree->looseId_ & nm1_hoe) == nm1_hoe)           h1_tp_hoe->Fill(       tree->hoe_);
                    if ((tree->looseId_ & nm1_ooemoop) == nm1_ooemoop)   h1_tp_ooemoop->Fill(   fabs(tree->ooemoop_));
                    if ((tree->looseId_ & nm1_d0vtx) == nm1_d0vtx)       h1_tp_d0vtx->Fill(     fabs(tree->d0vtx_));
                    if ((tree->looseId_ & nm1_dzvtx) == nm1_dzvtx)       h1_tp_dzvtx->Fill(     fabs(tree->dzvtx_));
                    if ((tree->looseId_ & nm1_vfitprob) == nm1_vfitprob) h1_tp_vfitprob->Fill(  tree->vfitprob_);
                    if ((tree->looseId_ & nm1_mhit) == nm1_mhit)         h1_tp_mhit->Fill(      tree->mhit_);
                }

            }

            //
            // fill FR histograms
            //

            if (tree->tagAndProbeMass_ < 0.0) {
                if (tree->probe_.Pt() < 30.0) {
                    h1_fr_nvtx->Fill(tree->nvtx_, weight);
                    if ((tree->looseId_ & nm1_detain) == nm1_detain)     h1_fr_detain->Fill(    fabs(tree->detain_));
                    if ((tree->looseId_ & nm1_dphiin) == nm1_dphiin)     h1_fr_dphiin->Fill(    fabs(tree->dphiin_));
                    if ((tree->looseId_ & nm1_sieie) == nm1_sieie)       h1_fr_sieie->Fill(     tree->sieie_);
                    if ((tree->looseId_ & nm1_hoe) == nm1_hoe)           h1_fr_hoe->Fill(       tree->hoe_);
                    if ((tree->looseId_ & nm1_ooemoop) == nm1_ooemoop)   h1_fr_ooemoop->Fill(   fabs(tree->ooemoop_));
                    if ((tree->looseId_ & nm1_d0vtx) == nm1_d0vtx)       h1_fr_d0vtx->Fill(     fabs(tree->d0vtx_));
                    if ((tree->looseId_ & nm1_dzvtx) == nm1_dzvtx)       h1_fr_dzvtx->Fill(     fabs(tree->dzvtx_));
                    if ((tree->looseId_ & nm1_vfitprob) == nm1_vfitprob) h1_fr_vfitprob->Fill(  tree->vfitprob_);
                    if ((tree->looseId_ & nm1_mhit) == nm1_mhit)         h1_fr_mhit->Fill(      tree->mhit_);
                }
            }

        } // end event loop

        delete tree;

    } // end file loop

    //
    // finish
    //

    gROOT->cd();

}

