
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

void LeptonTreeLooper::unsetGoodRunList()
{
    runlistIsSet_ = false;
}

void LeptonTreeLooper::loop(TChain *chain, TString name, unsigned int plotBin)
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

  TFile *pufile_el    = 0;
  TFile *pufile_mu    = 0;
  TH1D  *puWeights_el = 0;
  TH1D  *puWeights_mu = 0;
  //const TString pufname2012_el("/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_el.root");
  //const TString pufname2012_mu("/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_mu.root");
  //const TString pufname2012_el("/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_3000ipb.root");
  //const TString pufname2012_mu("/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_3000ipb.root");
    const TString pufname2012_el("/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_2400ipb.root");
    const TString pufname2012_mu("/smurf/data/Run2012_Summer12_SmurfV9_52X/auxiliar/puWeights_Summer12_2400ipb.root");

  pufile_el = new TFile(pufname2012_el);
  pufile_mu = new TFile(pufname2012_mu);
  puWeights_el = (TH1D*)pufile_el->Get("puWeights");
  puWeights_mu = (TH1D*)pufile_mu->Get("puWeights");

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

    unsigned int all_nohoe_noiso = DETAIN | DPHIIN | SIGMAIETAIETA | OOEMOOP | D0VTX | DZVTX | VTXFIT | MHITS;

    //
    // set up histograms
    //

    gROOT->cd();

    // general
    TH1F *h1_run = new TH1F(Form("%s_h1_run", name.Data()), "run", 200, 190600, 190600 + 200);

    // tag and probe
    TH1F *h1_mu_tp_nvtx = new TH1F(Form("%s_h1_mu_tp_nvtx", name.Data()), "nvtx; Number of PV", 40, -0.5, 39.5);
    TH1F *h1_el_tp_nvtx = new TH1F(Form("%s_h1_el_tp_nvtx", name.Data()), "nvtx; Number of PV", 40, -0.5, 39.5);
    h1_mu_tp_nvtx->Sumw2();
    h1_el_tp_nvtx->Sumw2();

    // electron BDT output studies
    TH1F *h1_el_bdt = new TH1F(Form("%s_h1_el_bdt", name.Data()), "BDT; BDT", 50, -1.0, 1.0);
    TH1F *h1_el_eopin = new TH1F(Form("%s_h1_el_eopin", name.Data()), "E/p; E/p", 40, 0.0, 2);
    TH1F *h1_el_hoe = new TH1F(Form("%s_h1_el_hoe", name.Data()), "H/E; H/E", 40, 0.0, 0.2);
    TH1F *h1_el_detain = new TH1F(Form("%s_h1_el_detain", name.Data()), "dEtaIn; dEtaIn", 40, -0.02, 0.02);
    TH1F *h1_el_dphiin = new TH1F(Form("%s_h1_el_dphiin", name.Data()), "dPhiIn; dPhiIn", 40, -0.2, 0.2);
    TH1F *h1_el_sieie = new TH1F(Form("%s_h1_el_sieie", name.Data()), "#sigma_{i#etai#eta}; #sigma_{i#etai#eta}", 50, 0.0, 0.05);

    h1_el_bdt->Sumw2();
    h1_el_eopin->Sumw2();
    h1_el_hoe->Sumw2();
    h1_el_detain->Sumw2();
    h1_el_dphiin->Sumw2();
    h1_el_sieie->Sumw2();


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

        unsigned int HLT_Mu8_probe_ = 1;
        unsigned int HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_ = 1;
        unsigned int HLT_IsoMu24_eta2p1_tag_ = 1;
        unsigned int HLT_Ele27_WP80_tag_ = 1;
        if (name != "dy_ee" && name != "dy_mm") {
            tree->tree_->SetBranchAddress("HLT_Mu8_probe", &HLT_Mu8_probe_);
            tree->tree_->SetBranchAddress("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe", &HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_);
            tree->tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_tag", &HLT_IsoMu24_eta2p1_tag_);
            tree->tree_->SetBranchAddress("HLT_Ele27_WP80_tag", &HLT_Ele27_WP80_tag_);
        }


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
            // plot bin
            //

            if (!testPlotBin(plotBin, tree)) continue;


            //
            // good run list
            //

            if (runlistIsSet_) {
                if (!goodrun_json(tree->run_, tree->lumi_)) continue;
            }

            h1_run->Fill(tree->run_, 1.0);

            //
            // nvtx weight
            //

            float weight = 1.0;
            if (name == "dy_ee") weight *= puWeights_el->GetBinContent(tree->npu_+1);
            if (name == "dy_mm") weight *= puWeights_mu->GetBinContent(tree->npu_+1);

            //
            // dataset specific cuts
            //

            if ((tree->eventSelection_ & LeptonTree::ZmmTagAndProbe) && HLT_IsoMu24_eta2p1_tag_ > 0
                    && tree->qTag_ * tree->qProbe_ < 0)
            {
                h1_mu_tp_nvtx->Fill(tree->nvtx_, weight);

            }


            if ((tree->eventSelection_ & LeptonTree::ZeeTagAndProbe) 
                    && tree->probe_.Pt() > 10.0 
                    && HLT_Ele27_WP80_tag_ > 0
                    && tree->qTag_ * tree->qProbe_ < 0)
            {

                h1_el_tp_nvtx->Fill(tree->nvtx_, weight);


                // pass isolation
                // pass id

                if (passElectronFO2012(tree) && passElectronIso2012(tree)
                            && fabs(tree->tagAndProbeMass_ - 91) < 10.0) {
                    h1_el_bdt->Fill(tree->egammaPOG2012MVA_, weight);

                    if (tree->egammaPOG2012MVA_ > 0.6) {
                        h1_el_eopin->Fill(tree->eopin_, weight);
                        h1_el_hoe->Fill(tree->hoe_, weight);
                        h1_el_detain->Fill(tree->detain_, weight);
                        h1_el_dphiin->Fill(tree->dphiin_, weight);
                        h1_el_sieie->Fill(tree->sieie_, weight);
                    }
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

bool LeptonTreeLooper::passElectronIso2012(const LeptonTree *leptonTree)
{
          float EffectiveArea = 0.0;
          if (fabs(leptonTree->sceta_) >= 0.0 && fabs(leptonTree->sceta_) < 1.0 ) EffectiveArea = 0.176;
          if (fabs(leptonTree->sceta_) >= 1.0 && fabs(leptonTree->sceta_) < 1.479 ) EffectiveArea = 0.206;
          if (fabs(leptonTree->sceta_) >= 1.479 && fabs(leptonTree->sceta_) < 2.0 ) EffectiveArea = 0.094;
          if (fabs(leptonTree->sceta_) >= 2.2 && fabs(leptonTree->sceta_) < 2.2 ) EffectiveArea = 0.172;
          if (fabs(leptonTree->sceta_) >= 2.3 && fabs(leptonTree->sceta_) < 2.3 ) EffectiveArea = 0.244;
          if (fabs(leptonTree->sceta_) >= 2.4 && fabs(leptonTree->sceta_) < 2.4 ) EffectiveArea = 0.333;
          if (fabs(leptonTree->sceta_) >= 2.4 ) EffectiveArea = 0.348;

        if ((leptonTree->pfchiso04_ + TMath::Max(float(0.0), leptonTree->pfemiso04_ + leptonTree->pfnhiso04_
            - EffectiveArea * TMath::Max(float(0.0), leptonTree->rhoIsoAll_)))/leptonTree->probe_.Pt() > 0.15) return false;
    return true;
}

bool LeptonTreeLooper::passElectronFO2012(const LeptonTree *leptonTree)
{

        float pt = leptonTree->probe_.Pt();
        float d0 = leptonTree->d0vtx_;
        float dz = leptonTree->dzvtx_;
        if (leptonTree->trkiso_/pt               > 0.2)      return false;
        if (leptonTree->hcaliso_/pt              > 0.2)      return false;
        if (fabs(d0) > 0.02)                                return false;
        if (fabs(dz) > 0.1)                                 return false;

        unsigned int mhits = leptonTree->mhit_;
        bool conv = leptonTree->vfitprob_;
        if (mhits > 0)      return false;
        if (conv)           return false;

        if (fabs(leptonTree->sceta_) < 1.479) {
            if (leptonTree->sieie_               > 0.01)  return false;
            if (fabs(leptonTree->detain_)        > 0.007) return false;
            if (fabs(leptonTree->dphiin_)        > 0.15)  return false;
            if (leptonTree->hoe_                 > 0.12)  return false;
            if ((leptonTree->ecaliso_ - 1.0)/pt  > 0.2)   return false;
        } else {
            if (leptonTree->sieie_               > 0.03)  return false;
            if (fabs(leptonTree->detain_)        > 0.009) return false;
            if (fabs(leptonTree->dphiin_)        > 0.10)  return false;
            if (leptonTree->hoe_                 > 0.10)  return false;
            if ((leptonTree->ecaliso_)/pt        > 0.2)   return false;
        }

    return true;

}

bool LeptonTreeLooper::testPlotBin(const unsigned int plotBin, const LeptonTree *tree)
{

    float pt = tree->probe_.Pt();

    float eta = fabs(tree->sceta_);
    if ((tree->eventSelection_ & LeptonTree::QCDFakeMu) || (tree->eventSelection_ & LeptonTree::ZmmTagAndProbe))
        eta = fabs(tree->probe_.Eta());

    if ((plotBin & plotbin::EB) == plotbin::EB && !(eta < 1.479))                     return false;
    if ((plotBin & plotbin::CENTRAL) == plotbin::CENTRAL && !(eta < 0.80))            return false;
    if ((plotBin & plotbin::TRANSITION) == plotbin::TRANSITION && !(eta >= 0.80 && eta < 1.479))            return false;
    if ((plotBin & plotbin::EE) == plotbin::EE && !(eta >= 1.479))                    return false;
    if ((plotBin & plotbin::EELO) == plotbin::EELO && !(eta >= 1.479 && eta < 2.0))                    return false;
    if ((plotBin & plotbin::EEHI) == plotbin::EEHI && !(eta >= 2.0))                    return false;

    if ((plotBin & plotbin::PT1015) == plotbin::PT1015 && !(pt >= 10.0 && pt < 15.0)) return false;
    if ((plotBin & plotbin::PT1020) == plotbin::PT1020 && !(pt >= 10.0 && pt < 20.0)) return false;
    if ((plotBin & plotbin::PT20UP) == plotbin::PT20UP && !(pt >= 20.0))              return false;
    if ((plotBin & plotbin::PT10UP) == plotbin::PT10UP && !(pt >= 10.0))              return false;

    return true;
}

