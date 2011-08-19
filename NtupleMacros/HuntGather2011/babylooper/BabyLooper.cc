
#include "BabyLooper.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

#include "../../Tools/DileptonHypType.h"
#include "../../Tools/goodrun.cc"
#include "../plotter/BabyDorkIdentifier.cc"

//
// set good run list
//

enum DileptonHypType hyp_typeToHypType (int hyp_type)
{
     switch (hyp_type) {
     case 0:
      return DILEPTON_MUMU;
     case 1: case 2:
      return DILEPTON_EMU;
     case 3:
      return DILEPTON_EE;
     default:
      assert(hyp_type < 4);
     }
     return DILEPTON_ALL;
}

void BabyLooper::setGoodRunList(std::string fname) 
{
    set_goodrun_file(fname.c_str());
}

void BabyLooper::Fill(TH1F** hist, const unsigned int hyp, const float &val, const float &weight)
{
    hist[hyp]->Fill(val, weight);
    hist[DILEPTON_ALL]->Fill(val, weight);
}

void BabyLooper::Fill2D(TH2F** hist, const unsigned int hyp, const float &valx, const float &valy, const float &weight)
{
    hist[hyp]->Fill(valx, valy, weight);
    hist[DILEPTON_ALL]->Fill(valx, valy, weight);
}

void BabyLooper::FormatHist(TH1F** hist, std::string sampleName, std::string name, int n, float min, float max)
{
    // loop on EB, EE
    for (unsigned int i = 0; i < 4; ++i)
    {
        std::string str = dilepton_hypo_names[i];
        std::string title = name + "_" + str;
        hist[i] = new TH1F(Form("%s_%s_%s", sampleName.c_str(), name.c_str(), str.c_str()),
                title.c_str(), n, min, max);
        hist[i]->GetXaxis()->SetTitle(name.c_str());
        hist[i]->Sumw2();
    }
}

void BabyLooper::FormatHist2D(TH2F** hist, std::string sampleName, std::string name,
        int nx, float minx, float maxx, int ny, float miny, float maxy)
{
    // loop on EB, EE 
    for (unsigned int i = 0; i < 4; ++i)
    {
        std::string str = dilepton_hypo_names[i];
        std::string title = name + "_" + str;
        hist[i] = new TH2F(Form("%s_%s_%s", sampleName.c_str(), name.c_str(), str.c_str()),
                title.c_str(), nx, minx, maxx, ny, miny, maxy);
        hist[i]->GetXaxis()->SetTitle(name.c_str());
        hist[i]->Sumw2();
    }
}

enum JetBin {
    JET_0J,
    JET_1J,
    JET_GT2J,
    JET_ANY,
};
static const char jetbin_names[][128] = { "0j", "1j", "2j", "allj"};

void BabyLooper::Loop(std::string sampleName, TChain *chain, unsigned int runMin, unsigned int runMax)
{

    //
    // set up histograms
    //

    TH1F *h1_hyp_njets_etaincl_[4];
    TH1F *h1_hyp_njets_eta1_[4];
    FormatHist(h1_hyp_njets_etaincl_, sampleName, "h1_hyp_njets_etaincl", 5, -0.5, 4.5);
    FormatHist(h1_hyp_njets_eta1_, sampleName, "h1_hyp_njets_eta1", 5, -0.5, 4.5);

    TH1F *h1_hyp_njets_sig_etaincl_[4];
    TH1F *h1_hyp_njets_sig_eta1_[4];
    FormatHist(h1_hyp_njets_sig_etaincl_, sampleName, "h1_hyp_njets_sig_etaincl", 5, -0.5, 4.5);
    FormatHist(h1_hyp_njets_sig_eta1_, sampleName, "h1_hyp_njets_sig_eta1", 5, -0.5, 4.5);

    TH1F *h1_hyp_pt_etaincl_[4][4];
    TH1F *h1_hyp_pt_eta1_[4][4];

    TH1F *h1_hyp_mllj_etaincl_[4][4];
    TH1F *h1_hyp_mllj_eta1_[4][4];

    TH1F *h1_hyp_dphi_sig_etaincl_[4][4];
    TH1F *h1_hyp_dphi_sig_eta1_[4][4];


    for (unsigned int j = 0; j < 4; ++j) {
        std::string jetbin = jetbin_names[j];

        FormatHist(h1_hyp_pt_etaincl_[j], sampleName, "h1_hyp_pt_etaincl_" + jetbin, 50, 0.0, 400.0);
        FormatHist(h1_hyp_pt_eta1_[j], sampleName, "h1_hyp_pt_eta1_" + jetbin, 50, 0.0, 400.0);

        FormatHist(h1_hyp_mllj_etaincl_[j], sampleName, "h1_hyp_mllj_etaincl_" + jetbin, 50, 0.0, 500.0);
        FormatHist(h1_hyp_mllj_eta1_[j], sampleName, "h1_hyp_mllj_eta1_" + jetbin, 50, 0.0, 500.0);

        // delta phi between the two leptons
        // - one plot with |eta_Z| < 1.0
        // - other plot with no restriction of |eta_Z|
        FormatHist(h1_hyp_dphi_sig_etaincl_[j], sampleName, "h1_hyp_dphi_sig_etaincl_" + jetbin, 40, -4.0, 4.0);
        FormatHist(h1_hyp_dphi_sig_eta1_[j], sampleName, "h1_hyp_dphi_sig_eta1_" + jetbin, 40, -4.0, 4.0);

    }

    //
    // reset the duplicate list
    //

    reset_babydorkidentifier();

    //
    // do the looping
    //

    
    std::cout << sampleName << ":\t (" << runMin << ", " << runMax << ")" << std::endl;
    TObjArray *listOfFiles = chain->GetListOfFiles();
    unsigned int nEventsChain=0;
    nEventsChain = chain->GetEntries();
    unsigned int nEventsTotal = 0;

    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;
    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        TFile f(currentFile->GetTitle());
        TTree *tree = (TTree*)f.Get("tree");
        Init(tree);

        tree->GetEntry(0);

        //
        // Event Loop
        //

        unsigned int nEvents = tree->GetEntries();
        for( unsigned int event = 0; event < nEvents; ++event) {
            tree->GetEntry(event);
            ++nEventsTotal;

            //
            // define basic selection
            //

            if (isdata) {
                if (!goodrun(run, ls)) continue;
                if (run < runMin) continue;
                if (run > runMax) continue;
            }

            // event with more than two tracks
            if (ntrks <= 2) continue;

            // 20, 20 with good leptons
            if (pt1 < 20.0 || pt2 < 20.0) continue;
            if (abs(eormu1) == 11 && !e1_vbtf90full) continue;
            if (abs(eormu2) == 11 && !e2_vbtf90full) continue;
            if (abs(eormu1) == 13 && !(mu1_muonidfull && ! mu1_cosmic)) continue;
            if (abs(eormu2) == 13 && !(mu2_muonidfull && ! mu2_cosmic)) continue;

            // z
            if (mass < 76.0 || mass > 106.0) continue;
            if (abs(eormu1) != abs(eormu2)) continue;

            //
            // do duplicate check
            //

            if (isdata)
                if (is_duplicate(run,evt,ls,pt1,pt2)) continue;

            //
            // fill plots
            //

            JetBin jetBin = JET_0J;
            if (njets == 1) jetBin = JET_1J;
            if (njets >= 2) jetBin = JET_GT2J;

            DileptonHypType hypType = hyp_typeToHypType (hyp_type);
            float weight = scale1fb;
            if (isdata) weight = 1.0;


            // |eta_Z| inclusive region

            Fill(h1_hyp_njets_etaincl_, hypType, njets, weight);

            Fill(h1_hyp_pt_etaincl_[JET_ANY], hypType, dilpt, weight);
            Fill(h1_hyp_pt_etaincl_[jetBin], hypType, dilpt, weight);

            Fill(h1_hyp_mllj_etaincl_[JET_ANY],  hypType, mllj, weight);
            Fill(h1_hyp_mllj_etaincl_[jetBin],  hypType, mllj, weight);

            if (dilpt > 80 && dilpt < 100) {
                Fill(h1_hyp_njets_sig_etaincl_, hypType, njets, weight);
                Fill(h1_hyp_dphi_sig_etaincl_[JET_ANY], hypType, phi1 - phi2, weight);
                Fill(h1_hyp_dphi_sig_etaincl_[jetBin], hypType,phi1 - phi2, weight);
            }

            // |eta_Z| < 1.0 region
            if (fabs(dileta) < 1.0) {

                Fill(h1_hyp_njets_eta1_, hypType, njets, weight);

                Fill(h1_hyp_pt_eta1_[JET_ANY], hypType, dilpt, weight);
                Fill(h1_hyp_pt_eta1_[jetBin], hypType, dilpt, weight);

                Fill(h1_hyp_mllj_eta1_[JET_ANY],  hypType, mllj, weight);
                Fill(h1_hyp_mllj_eta1_[jetBin],  hypType, mllj, weight);

                if (dilpt > 80 && dilpt < 100) {
                    Fill(h1_hyp_njets_sig_eta1_, hypType, njets, weight);
                    Fill(h1_hyp_dphi_sig_eta1_[JET_ANY], hypType, phi1 - phi2, weight);
                    Fill(h1_hyp_dphi_sig_eta1_[jetBin], hypType, phi1 - phi2, weight);
                }


            }
    

        } // end loop on events

    } //end loop on files

}

BabyLooper::BabyLooper()
{
    std::cout << "[BabyLooper::BabyLooper]" << std::endl;
}

BabyLooper::~BabyLooper()
{
    std::cout << "[BabyLooper::~BabyLooper]" << std::endl;
}

void BabyLooper::Init(TTree *tree)
{   
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("rndm", &rndm, &b_rndm);
    fChain->SetBranchAddress("dataset", dataset, &b_dataset);
    fChain->SetBranchAddress("run", &run, &b_run);
    fChain->SetBranchAddress("ls", &ls, &b_ls);
    fChain->SetBranchAddress("evt", &evt, &b_evt);
    fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
    fChain->SetBranchAddress("isdata", &isdata, &b_isdata);
    fChain->SetBranchAddress("evt_clean082010", &evt_clean082010, &b_evt_clean082010);
    fChain->SetBranchAddress("evt_clean102010", &evt_clean102010, &b_evt_clean102010);
    fChain->SetBranchAddress("evt_clean042011", &evt_clean042011, &b_evt_clean042011);
    fChain->SetBranchAddress("scale1fb", &scale1fb, &b_scale1fb);
    fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
    fChain->SetBranchAddress("hyp_type", &hyp_type, &b_hyp_type);
    fChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
    fChain->SetBranchAddress("tcmet", &tcmet, &b_tcmet);
    fChain->SetBranchAddress("calotcmet", &calotcmet, &b_calotcmet);
    fChain->SetBranchAddress("proj_pfmet", &proj_pfmet, &b_proj_pfmet);
    fChain->SetBranchAddress("proj_tcmet", &proj_tcmet, &b_proj_tcmet);
    fChain->SetBranchAddress("ntrks", &ntrks, &b_ntrks);
    fChain->SetBranchAddress("njets", &njets, &b_njets);
    fChain->SetBranchAddress("njets25", &njets25, &b_njets25);
    fChain->SetBranchAddress("njetsSS", &njetsSS, &b_njetsSS);
    fChain->SetBranchAddress("jet1pt", &jet1pt, &b_jet1pt);
    fChain->SetBranchAddress("jet2pt", &jet2pt, &b_jet2pt);
    fChain->SetBranchAddress("jet3pt", &jet3pt, &b_jet3pt);
    fChain->SetBranchAddress("sumjetpt", &sumjetpt, &b_sumjetpt);
    fChain->SetBranchAddress("sumjetpt25", &sumjetpt25, &b_sumjetpt25);
    fChain->SetBranchAddress("sumjetptSS", &sumjetptSS, &b_sumjetptSS);
    fChain->SetBranchAddress("jet1eta", &jet1eta, &b_jet1eta);
    fChain->SetBranchAddress("jet2eta", &jet2eta, &b_jet2eta);
    fChain->SetBranchAddress("jet3eta", &jet3eta, &b_jet3eta);
    fChain->SetBranchAddress("jet1phi", &jet1phi, &b_jet1phi);
    fChain->SetBranchAddress("jet2phi", &jet2phi, &b_jet2phi);
    fChain->SetBranchAddress("jet3phi", &jet3phi, &b_jet3phi);
    fChain->SetBranchAddress("jetmass", &jetmass, &b_jetmass);
    fChain->SetBranchAddress("pfmth", &pfmth, &b_pfmth);
    fChain->SetBranchAddress("tcmth", &tcmth, &b_tcmth);
    fChain->SetBranchAddress("jet1isBtag", &jet1isBtag, &b_jet1isBtag);
    fChain->SetBranchAddress("jet2isBtag", &jet2isBtag, &b_jet2isBtag);
    fChain->SetBranchAddress("jet3isBtag", &jet3isBtag, &b_jet3isBtag);
    fChain->SetBranchAddress("dphipfmetjet", &dphipfmetjet, &b_dphipfmetjet);
    fChain->SetBranchAddress("dphitcmetjet", &dphitcmetjet, &b_dphitcmetjet);
    fChain->SetBranchAddress("deltaphi", &deltaphi, &b_deltaphi);
    fChain->SetBranchAddress("ntchelbtags", &ntchelbtags, &b_ntchelbtags);
    fChain->SetBranchAddress("nssvhembtags", &nssvhembtags, &b_nssvhembtags);
    fChain->SetBranchAddress("nssvhetbtags", &nssvhetbtags, &b_nssvhetbtags);
    fChain->SetBranchAddress("nssvhptbtags", &nssvhptbtags, &b_nssvhptbtags);
    fChain->SetBranchAddress("pfmeff", &pfmeff, &b_pfmeff);
    fChain->SetBranchAddress("tcmeff", &tcmeff, &b_tcmeff);
    fChain->SetBranchAddress("intLumiPerLS", &intLumiPerLS, &b_intLumiPerLS);
    fChain->SetBranchAddress("ngoodlep", &ngoodlep, &b_ngoodlep);
    fChain->SetBranchAddress("ngoodmus", &ngoodmus, &b_ngoodmus);
    fChain->SetBranchAddress("ngoodels", &ngoodels, &b_ngoodels);
    fChain->SetBranchAddress("ngoodlepSS", &ngoodlepSS, &b_ngoodlepSS);
    fChain->SetBranchAddress("ngoodmusSS", &ngoodmusSS, &b_ngoodmusSS);
    fChain->SetBranchAddress("ngoodelsSS", &ngoodelsSS, &b_ngoodelsSS);
    fChain->SetBranchAddress("dilpt", &dilpt, &b_dilpt);
    fChain->SetBranchAddress("dileta", &dileta, &b_dileta);
    fChain->SetBranchAddress("dilphi", &dilphi, &b_dilphi);
    fChain->SetBranchAddress("mass", &mass, &b_mass);
    fChain->SetBranchAddress("mlljj", &mlljj, &b_mlljj);
    fChain->SetBranchAddress("mllj", &mllj, &b_mllj);
    fChain->SetBranchAddress("eormu1", &eormu1, &b_eormu1);
    fChain->SetBranchAddress("type1", &type1, &b_type1);
    fChain->SetBranchAddress("ngenels", &ngenels, &b_ngenels);
    fChain->SetBranchAddress("ngenmus", &ngenmus, &b_ngenmus);
    fChain->SetBranchAddress("ngentaus", &ngentaus, &b_ngentaus);
    fChain->SetBranchAddress("pt1", &pt1, &b_pt1);
    fChain->SetBranchAddress("eta1", &eta1, &b_eta1);
    fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
    fChain->SetBranchAddress("iso1", &iso1, &b_iso1);
    fChain->SetBranchAddress("ntiso1", &ntiso1, &b_ntiso1);
    fChain->SetBranchAddress("d0corr1", &d0corr1, &b_d0corr1);
    fChain->SetBranchAddress("d0vtx1", &d0vtx1, &b_d0vtx1);
    fChain->SetBranchAddress("dphipfmet1", &dphipfmet1, &b_dphipfmet1);
    fChain->SetBranchAddress("dphitcmet1", &dphitcmet1, &b_dphitcmet1);
    fChain->SetBranchAddress("drjet1", &drjet1, &b_drjet1);
    fChain->SetBranchAddress("mcid1", &mcid1, &b_mcid1);
    fChain->SetBranchAddress("mcmotherid1", &mcmotherid1, &b_mcmotherid1);
    fChain->SetBranchAddress("eormu2", &eormu2, &b_eormu2);
    fChain->SetBranchAddress("type2", &type2, &b_type2);
    fChain->SetBranchAddress("pt2", &pt2, &b_pt2);
    fChain->SetBranchAddress("eta2", &eta2, &b_eta2);
    fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
    fChain->SetBranchAddress("iso2", &iso2, &b_iso2);
    fChain->SetBranchAddress("ntiso2", &ntiso2, &b_ntiso2);
    fChain->SetBranchAddress("d0corr2", &d0corr2, &b_d0corr2);
    fChain->SetBranchAddress("d0vtx2", &d0vtx2, &b_d0vtx2);
    fChain->SetBranchAddress("dphipfmet2", &dphipfmet2, &b_dphipfmet2);
    fChain->SetBranchAddress("dphitcmet2", &dphitcmet2, &b_dphitcmet2);
    fChain->SetBranchAddress("drjet2", &drjet2, &b_drjet2);
    fChain->SetBranchAddress("mcid2", &mcid2, &b_mcid2);
    fChain->SetBranchAddress("mcmotherid2", &mcmotherid2, &b_mcmotherid2);
    fChain->SetBranchAddress("mt2", &mt2, &b_mt2);
    fChain->SetBranchAddress("mt2j", &mt2j, &b_mt2j);
    fChain->SetBranchAddress("extraZveto", &extraZveto, &b_extraZveto);
    fChain->SetBranchAddress("trkIso1", &trkIso1, &b_trkIso1);
    fChain->SetBranchAddress("ecalIso1", &ecalIso1, &b_ecalIso1);
    fChain->SetBranchAddress("hcalIso1", &hcalIso1, &b_hcalIso1);
    fChain->SetBranchAddress("trkIso2", &trkIso2, &b_trkIso2);
    fChain->SetBranchAddress("ecalIso2", &ecalIso2, &b_ecalIso2);
    fChain->SetBranchAddress("hcalIso2", &hcalIso2, &b_hcalIso2);
    fChain->SetBranchAddress("ecalIso1ps", &ecalIso1ps, &b_ecalIso1ps);
    fChain->SetBranchAddress("ecalIso2ps", &ecalIso2ps, &b_ecalIso2ps);
    fChain->SetBranchAddress("rho", &rho, &b_rho);
    fChain->SetBranchAddress("lepsFromSameVtx", &lepsFromSameVtx, &b_lepsFromSameVtx);
    fChain->SetBranchAddress("lep1isFromW", &lep1isFromW, &b_lep1isFromW);
    fChain->SetBranchAddress("lep2isFromW", &lep2isFromW, &b_lep2isFromW);
    fChain->SetBranchAddress("mu1_numSSv3", &mu1_numSSv3, &b_mu1_numSSv3);
    fChain->SetBranchAddress("mu1_foSSv3", &mu1_foSSv3, &b_mu1_foSSv3);
    fChain->SetBranchAddress("mu1_muonidfull", &mu1_muonidfull, &b_mu1_muonidfull);
    fChain->SetBranchAddress("mu1_muonid", &mu1_muonid, &b_mu1_muonid);
    fChain->SetBranchAddress("mu1_muonidfullV1", &mu1_muonidfullV1, &b_mu1_muonidfullV1);
    fChain->SetBranchAddress("mu1_muonidV1", &mu1_muonidV1, &b_mu1_muonidV1);
    fChain->SetBranchAddress("mu1_goodmask", &mu1_goodmask, &b_mu1_goodmask);
    fChain->SetBranchAddress("mu1_gfitchi2", &mu1_gfitchi2, &b_mu1_gfitchi2);
    fChain->SetBranchAddress("mu1_cosmic", &mu1_cosmic, &b_mu1_cosmic);
    fChain->SetBranchAddress("mu1_siHits", &mu1_siHits, &b_mu1_siHits);
    fChain->SetBranchAddress("mu1_saHits", &mu1_saHits, &b_mu1_saHits);
    fChain->SetBranchAddress("mu1_emVetoDep", &mu1_emVetoDep, &b_mu1_emVetoDep);
    fChain->SetBranchAddress("mu1_hadVetoDep", &mu1_hadVetoDep, &b_mu1_hadVetoDep);
    fChain->SetBranchAddress("mu1_isPFmuon", &mu1_isPFmuon, &b_mu1_isPFmuon);
    fChain->SetBranchAddress("mu2_numSSv3", &mu2_numSSv3, &b_mu2_numSSv3);
    fChain->SetBranchAddress("mu2_foSSv3", &mu2_foSSv3, &b_mu2_foSSv3);
    fChain->SetBranchAddress("mu2_muonidfull", &mu2_muonidfull, &b_mu2_muonidfull);
    fChain->SetBranchAddress("mu2_muonid", &mu2_muonid, &b_mu2_muonid);
    fChain->SetBranchAddress("mu2_muonidfullV1", &mu2_muonidfullV1, &b_mu2_muonidfullV1);
    fChain->SetBranchAddress("mu2_muonidV1", &mu2_muonidV1, &b_mu2_muonidV1);
    fChain->SetBranchAddress("mu2_goodmask", &mu2_goodmask, &b_mu2_goodmask);
    fChain->SetBranchAddress("mu2_gfitchi2", &mu2_gfitchi2, &b_mu2_gfitchi2);
    fChain->SetBranchAddress("mu2_cosmic", &mu2_cosmic, &b_mu2_cosmic);
    fChain->SetBranchAddress("mu2_siHits", &mu2_siHits, &b_mu2_siHits);
    fChain->SetBranchAddress("mu2_saHits", &mu2_saHits, &b_mu2_saHits);
    fChain->SetBranchAddress("mu2_emVetoDep", &mu2_emVetoDep, &b_mu2_emVetoDep);
    fChain->SetBranchAddress("mu2_hadVetoDep", &mu2_hadVetoDep, &b_mu2_hadVetoDep);
    fChain->SetBranchAddress("mu2_isPFmuon", &mu2_isPFmuon, &b_mu2_isPFmuon);
    fChain->SetBranchAddress("e1_numSSv3", &e1_numSSv3, &b_e1_numSSv3);
    fChain->SetBranchAddress("e1_foSSv3", &e1_foSSv3, &b_e1_foSSv3);
    fChain->SetBranchAddress("e1_vbtf90full", &e1_vbtf90full, &b_e1_vbtf90full);
    fChain->SetBranchAddress("e1_vbtf90", &e1_vbtf90, &b_e1_vbtf90);
    fChain->SetBranchAddress("e1_vbtf85", &e1_vbtf85, &b_e1_vbtf85);
    fChain->SetBranchAddress("e1_vbtf80", &e1_vbtf80, &b_e1_vbtf80);
    fChain->SetBranchAddress("e1_vbtf70", &e1_vbtf70, &b_e1_vbtf70);
    fChain->SetBranchAddress("e1_smurfV3", &e1_smurfV3, &b_e1_smurfV3);
    fChain->SetBranchAddress("e1_scet", &e1_scet, &b_e1_scet);
    fChain->SetBranchAddress("e1_eopin", &e1_eopin, &b_e1_eopin);
    fChain->SetBranchAddress("e1_hoe", &e1_hoe, &b_e1_hoe);
    fChain->SetBranchAddress("e1_dphiin", &e1_dphiin, &b_e1_dphiin);
    fChain->SetBranchAddress("e1_detain", &e1_detain, &b_e1_detain);
    fChain->SetBranchAddress("e1_e25Me55", &e1_e25Me55, &b_e1_e25Me55);
    fChain->SetBranchAddress("e1_sigieie", &e1_sigieie, &b_e1_sigieie);
    fChain->SetBranchAddress("e1_eMe55", &e1_eMe55, &b_e1_eMe55);
    fChain->SetBranchAddress("e1_nmHits", &e1_nmHits, &b_e1_nmHits);
    fChain->SetBranchAddress("e1_dcot", &e1_dcot, &b_e1_dcot);
    fChain->SetBranchAddress("e1_dist", &e1_dist, &b_e1_dist);
    fChain->SetBranchAddress("e1_drmu", &e1_drmu, &b_e1_drmu);
    fChain->SetBranchAddress("e1_isspike", &e1_isspike, &b_e1_isspike);
    fChain->SetBranchAddress("e1_ctfCharge", &e1_ctfCharge, &b_e1_ctfCharge);
    fChain->SetBranchAddress("e1_gsfCharge", &e1_gsfCharge, &b_e1_gsfCharge);
    fChain->SetBranchAddress("e1_scCharge", &e1_scCharge, &b_e1_scCharge);
    fChain->SetBranchAddress("e1_fbrem", &e1_fbrem, &b_e1_fbrem);
    fChain->SetBranchAddress("e1_mitConv", &e1_mitConv, &b_e1_mitConv);
    fChain->SetBranchAddress("e2_numSSv3", &e2_numSSv3, &b_e2_numSSv3);
    fChain->SetBranchAddress("e2_foSSv3", &e2_foSSv3, &b_e2_foSSv3);
    fChain->SetBranchAddress("e2_vbtf90full", &e2_vbtf90full, &b_e2_vbtf90full);
    fChain->SetBranchAddress("e2_vbtf90", &e2_vbtf90, &b_e2_vbtf90);
    fChain->SetBranchAddress("e2_vbtf85", &e2_vbtf85, &b_e2_vbtf85);
    fChain->SetBranchAddress("e2_vbtf80", &e2_vbtf80, &b_e2_vbtf80);
    fChain->SetBranchAddress("e2_vbtf70", &e2_vbtf70, &b_e2_vbtf70);
    fChain->SetBranchAddress("e2_smurfV3", &e2_smurfV3, &b_e2_smurfV3);
    fChain->SetBranchAddress("e2_scet", &e2_scet, &b_e2_scet);
    fChain->SetBranchAddress("e2_eopin", &e2_eopin, &b_e2_eopin);
    fChain->SetBranchAddress("e2_hoe", &e2_hoe, &b_e2_hoe);
    fChain->SetBranchAddress("e2_dphiin", &e2_dphiin, &b_e2_dphiin);
    fChain->SetBranchAddress("e2_detain", &e2_detain, &b_e2_detain);
    fChain->SetBranchAddress("e2_e25Me55", &e2_e25Me55, &b_e2_e25Me55);
    fChain->SetBranchAddress("e2_sigieie", &e2_sigieie, &b_e2_sigieie);
    fChain->SetBranchAddress("e2_eMe55", &e2_eMe55, &b_e2_eMe55);
    fChain->SetBranchAddress("e2_nmHits", &e2_nmHits, &b_e2_nmHits);
    fChain->SetBranchAddress("e2_dcot", &e2_dcot, &b_e2_dcot);
    fChain->SetBranchAddress("e2_dist", &e2_dist, &b_e2_dist);
    fChain->SetBranchAddress("e2_drmu", &e2_drmu, &b_e2_drmu);
    fChain->SetBranchAddress("e2_isspike", &e2_isspike, &b_e2_isspike);
    fChain->SetBranchAddress("e2_ctfCharge", &e2_ctfCharge, &b_e2_ctfCharge);
    fChain->SetBranchAddress("e2_gsfCharge", &e2_gsfCharge, &b_e2_gsfCharge);
    fChain->SetBranchAddress("e2_scCharge", &e2_scCharge, &b_e2_scCharge);
    fChain->SetBranchAddress("e2_fbrem", &e2_fbrem, &b_e2_fbrem);
    fChain->SetBranchAddress("e2_mitConv", &e2_mitConv, &b_e2_mitConv);
    fChain->SetBranchAddress("trg_single_mu1", &trg_single_mu1, &b_trg_single_mu1);
    fChain->SetBranchAddress("trg_single_mu2", &trg_single_mu2, &b_trg_single_mu2);
    fChain->SetBranchAddress("trg_single_e1", &trg_single_e1, &b_trg_single_e1);
    fChain->SetBranchAddress("trg_single_e2", &trg_single_e2, &b_trg_single_e2);
    fChain->SetBranchAddress("trg_double_mu1", &trg_double_mu1, &b_trg_double_mu1);
    fChain->SetBranchAddress("trg_double_mu2", &trg_double_mu2, &b_trg_double_mu2);
    fChain->SetBranchAddress("trg_double_e1", &trg_double_e1, &b_trg_double_e1);
    fChain->SetBranchAddress("trg_double_e2", &trg_double_e2, &b_trg_double_e2);
    fChain->SetBranchAddress("trg_cross_emu", &trg_cross_emu, &b_trg_cross_emu);
    fChain->SetBranchAddress("trg_had_double_e1", &trg_had_double_e1, &b_trg_had_double_e1);
    fChain->SetBranchAddress("trg_had_double_e2", &trg_had_double_e2, &b_trg_had_double_e2);
    fChain->SetBranchAddress("trg_had_double_mu1", &trg_had_double_mu1, &b_trg_had_double_mu1);
    fChain->SetBranchAddress("trg_had_double_mu2", &trg_had_double_mu2, &b_trg_had_double_mu2);
    fChain->SetBranchAddress("trg_had_cross_emu", &trg_had_cross_emu, &b_trg_had_cross_emu);
}

