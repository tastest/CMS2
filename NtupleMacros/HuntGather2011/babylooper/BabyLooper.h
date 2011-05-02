
#ifndef BabyLooper_h
#define BabyLooper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
class TH1F;
class TH2F;

#include <iostream>

class BabyLooper {
    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent; //!current Tree number in a TChain

        // Declaration of leaf types
        Float_t         rndm;
        Char_t          dataset[200];
        Int_t           run;
        Int_t           ls;
        Int_t           evt;
        Int_t           nvtx;
        Int_t           isdata;
        Int_t           evt_clean082010;
        Int_t           evt_clean102010;
        Int_t           evt_clean042011;
        Float_t         scale1fb;
        Float_t         pthat;
        Int_t           hyp_type;
        Float_t         pfmet;
        Float_t         tcmet;
        Float_t         calotcmet;
        Float_t         proj_pfmet;
        Float_t         proj_tcmet;
        Int_t           ntrks;
        Int_t           njets;
        Int_t           njets25;
        Int_t           njetsSS;
        Float_t         jet1pt;
        Float_t         jet2pt;
        Float_t         jet3pt;
        Float_t         sumjetpt;
        Float_t         sumjetpt25;
        Float_t         sumjetptSS;
        Float_t         jet1eta;
        Float_t         jet2eta;
        Float_t         jet3eta;
        Float_t         jet1phi;
        Float_t         jet2phi;
        Float_t         jet3phi;
        Float_t         jetmass;
        Float_t         pfmth;
        Float_t         tcmth;
        Int_t           jet1isBtag;
        Int_t           jet2isBtag;
        Int_t           jet3isBtag;
        Float_t         dphipfmetjet;
        Float_t         dphitcmetjet;
        Float_t         deltaphi;
        Int_t           ntchelbtags;
        Int_t           nssvhembtags;
        Int_t           nssvhetbtags;
        Int_t           nssvhptbtags;
        Float_t         pfmeff;
        Float_t         tcmeff;
        Float_t         intLumiPerLS;
        Int_t           ngoodlep;
        Int_t           ngoodmus;
        Int_t           ngoodels;
        Int_t           ngoodlepSS;
        Int_t           ngoodmusSS;
        Int_t           ngoodelsSS;
        Float_t         dilpt;
        Float_t         dileta;
        Float_t         dilphi;
        Float_t         mass;
        Float_t         mlljj;
        Float_t         mllj;
        Int_t           eormu1;
        Int_t           type1;
        Int_t           ngenels;
        Int_t           ngenmus;
        Int_t           ngentaus;
        Float_t         pt1;
        Float_t         eta1;
        Float_t         phi1;
        Float_t         iso1;
        Float_t         ntiso1;
        Float_t         d0corr1;
        Float_t         d0vtx1;
        Float_t         dphipfmet1;
        Float_t         dphitcmet1;
        Float_t         drjet1;
        Int_t           mcid1;
        Int_t           mcmotherid1;
        Int_t           eormu2;
        Int_t           type2;
        Float_t         pt2;
        Float_t         eta2;
        Float_t         phi2;
        Float_t         iso2;
        Float_t         ntiso2;
        Float_t         d0corr2;
        Float_t         d0vtx2;
        Float_t         dphipfmet2;
        Float_t         dphitcmet2;
        Float_t         drjet2;
        Int_t           mcid2;
        Int_t           mcmotherid2;
        Float_t         mt2;
        Float_t         mt2j;
        Bool_t          extraZveto;
        Float_t         trkIso1;
        Float_t         ecalIso1;
        Float_t         hcalIso1;
        Float_t         trkIso2;
        Float_t         ecalIso2;
        Float_t         hcalIso2;
        Float_t         ecalIso1ps;
        Float_t         ecalIso2ps;
        Float_t         rho;
        Bool_t          lepsFromSameVtx;
        Int_t           lep1isFromW;
        Int_t           lep2isFromW;
        Bool_t          mu1_numSSv3;
        Bool_t          mu1_foSSv3;
        Bool_t          mu1_muonidfull;
        Bool_t          mu1_muonid;
        Bool_t          mu1_muonidfullV1;
        Bool_t          mu1_muonidV1;
        Int_t           mu1_goodmask;
        Float_t         mu1_gfitchi2;
        Bool_t          mu1_cosmic;
        Int_t           mu1_siHits;
        Int_t           mu1_saHits;
        Float_t         mu1_emVetoDep;
        Float_t         mu1_hadVetoDep;
        Bool_t          mu1_isPFmuon;
        Bool_t          mu2_numSSv3;
        Bool_t          mu2_foSSv3;
        Bool_t          mu2_muonidfull;
        Bool_t          mu2_muonid;
        Bool_t          mu2_muonidfullV1;
        Bool_t          mu2_muonidV1;
        Int_t           mu2_goodmask;
        Float_t         mu2_gfitchi2;
        Bool_t          mu2_cosmic;
        Int_t           mu2_siHits;
        Int_t           mu2_saHits;
        Float_t         mu2_emVetoDep;
        Float_t         mu2_hadVetoDep;
        Bool_t          mu2_isPFmuon;
        Bool_t          e1_numSSv3;
        Bool_t          e1_foSSv3;
        Bool_t          e1_vbtf90full;
        Bool_t          e1_vbtf90;
        Bool_t          e1_vbtf85;
        Bool_t          e1_vbtf80;
        Bool_t          e1_vbtf70;
        Bool_t          e1_smurfV3;
        Float_t         e1_scet;
        Float_t         e1_eopin;
        Float_t         e1_hoe;
        Float_t         e1_dphiin;
        Float_t         e1_detain;
        Float_t         e1_e25Me55;
        Float_t         e1_sigieie;
        Float_t         e1_eMe55;
        Int_t           e1_nmHits;
        Float_t         e1_dcot;
        Float_t         e1_dist;
        Float_t         e1_drmu;
        Bool_t          e1_isspike;
        Int_t           e1_ctfCharge;
        Int_t           e1_gsfCharge;
        Int_t           e1_scCharge;
        Float_t         e1_fbrem;
        Bool_t          e1_mitConv;
        Bool_t          e2_numSSv3;
        Bool_t          e2_foSSv3;
        Bool_t          e2_vbtf90full;
        Bool_t          e2_vbtf90;
        Bool_t          e2_vbtf85;
        Bool_t          e2_vbtf80;
        Bool_t          e2_vbtf70;
        Bool_t          e2_smurfV3;
        Float_t         e2_scet;
        Float_t         e2_eopin;
        Float_t         e2_hoe;
        Float_t         e2_dphiin;
        Float_t         e2_detain;
        Float_t         e2_e25Me55;
        Float_t         e2_sigieie;
        Float_t         e2_eMe55;
        Int_t           e2_nmHits;
        Float_t         e2_dcot;
        Float_t         e2_dist;
        Float_t         e2_drmu;
        Bool_t          e2_isspike;
        Int_t           e2_ctfCharge;
        Int_t           e2_gsfCharge;
        Int_t           e2_scCharge;
        Float_t         e2_fbrem;
        Bool_t          e2_mitConv;
        Int_t           trg_single_mu1;
        Int_t           trg_single_mu2;
        Int_t           trg_single_e1;
        Int_t           trg_single_e2;
        Int_t           trg_double_mu1;
        Int_t           trg_double_mu2;
        Int_t           trg_double_e1;
        Int_t           trg_double_e2;
        Int_t           trg_cross_emu;
        Int_t           trg_had_double_e1;
        Int_t           trg_had_double_e2;
        Int_t           trg_had_double_mu1;
        Int_t           trg_had_double_mu2;
        Int_t           trg_had_cross_emu;

        // List of branches
        TBranch        *b_rndm;   //!
        TBranch        *b_dataset;   //!
        TBranch        *b_run;   //!
        TBranch        *b_ls;   //!
        TBranch        *b_evt;   //!
        TBranch        *b_nvtx;   //!
        TBranch        *b_isdata;   //!
        TBranch        *b_evt_clean082010;   //!
        TBranch        *b_evt_clean102010;   //!
        TBranch        *b_evt_clean042011;   //!
        TBranch        *b_scale1fb;   //!
        TBranch        *b_pthat;   //!
        TBranch        *b_hyp_type;   //!
        TBranch        *b_pfmet;   //!
        TBranch        *b_tcmet;   //!
        TBranch        *b_calotcmet;   //!
        TBranch        *b_proj_pfmet;   //!
        TBranch        *b_proj_tcmet;   //!
        TBranch        *b_ntrks;   //!
        TBranch        *b_njets;   //!
        TBranch        *b_njets25;   //!
        TBranch        *b_njetsSS;   //!
        TBranch        *b_jet1pt;   //!
        TBranch        *b_jet2pt;   //!
        TBranch        *b_jet3pt;   //!
        TBranch        *b_sumjetpt;   //!
        TBranch        *b_sumjetpt25;   //!
        TBranch        *b_sumjetptSS;   //!
        TBranch        *b_jet1eta;   //!
        TBranch        *b_jet2eta;   //!
        TBranch        *b_jet3eta;   //!
        TBranch        *b_jet1phi;   //!
        TBranch        *b_jet2phi;   //!
        TBranch        *b_jet3phi;   //!
        TBranch        *b_jetmass;   //!
        TBranch        *b_pfmth;   //!
        TBranch        *b_tcmth;   //!
        TBranch        *b_jet1isBtag;   //!
        TBranch        *b_jet2isBtag;   //!
        TBranch        *b_jet3isBtag;   //!
        TBranch        *b_dphipfmetjet;   //!
        TBranch        *b_dphitcmetjet;   //!
        TBranch        *b_deltaphi;   //!
        TBranch        *b_ntchelbtags;   //!
        TBranch        *b_nssvhembtags;   //!
        TBranch        *b_nssvhetbtags;   //!
        TBranch        *b_nssvhptbtags;   //!
        TBranch        *b_pfmeff;   //!
        TBranch        *b_tcmeff;   //!
        TBranch        *b_intLumiPerLS;   //!
        TBranch        *b_ngoodlep;   //!
        TBranch        *b_ngoodmus;   //!
        TBranch        *b_ngoodels;   //!
        TBranch        *b_ngoodlepSS;   //!
        TBranch        *b_ngoodmusSS;   //!
        TBranch        *b_ngoodelsSS;   //!
        TBranch        *b_dilpt;   //!
        TBranch        *b_dileta;   //!
        TBranch        *b_dilphi;   //!
        TBranch        *b_mass;   //!
        TBranch        *b_mlljj;   //!
        TBranch        *b_mllj;   //!
        TBranch        *b_eormu1;   //!
        TBranch        *b_type1;   //!
        TBranch        *b_ngenels;   //!
        TBranch        *b_ngenmus;   //!
        TBranch        *b_ngentaus;   //!
        TBranch        *b_pt1;   //!
        TBranch        *b_eta1;   //!
        TBranch        *b_phi1;   //!
        TBranch        *b_iso1;   //!
        TBranch        *b_ntiso1;   //!
        TBranch        *b_d0corr1;   //!
        TBranch        *b_d0vtx1;   //!
        TBranch        *b_dphipfmet1;   //!
        TBranch        *b_dphitcmet1;   //!
        TBranch        *b_drjet1;   //!
        TBranch        *b_mcid1;   //!
        TBranch        *b_mcmotherid1;   //!
        TBranch        *b_eormu2;   //!
        TBranch        *b_type2;   //!
        TBranch        *b_pt2;   //!
        TBranch        *b_eta2;   //!
        TBranch        *b_phi2;   //!
        TBranch        *b_iso2;   //!
        TBranch        *b_ntiso2;   //!
        TBranch        *b_d0corr2;   //!
        TBranch        *b_d0vtx2;   //!
        TBranch        *b_dphipfmet2;   //!
        TBranch        *b_dphitcmet2;   //!
        TBranch        *b_drjet2;   //!
        TBranch        *b_mcid2;   //!
        TBranch        *b_mcmotherid2;   //!
        TBranch        *b_mt2;   //!
        TBranch        *b_mt2j;   //!
        TBranch        *b_extraZveto;   //!
        TBranch        *b_trkIso1;   //!
        TBranch        *b_ecalIso1;   //!
        TBranch        *b_hcalIso1;   //!
        TBranch        *b_trkIso2;   //!
        TBranch        *b_ecalIso2;   //!
        TBranch        *b_hcalIso2;   //!
        TBranch        *b_ecalIso1ps;   //!
        TBranch        *b_ecalIso2ps;   //!
        TBranch        *b_rho;   //!
        TBranch        *b_lepsFromSameVtx;   //!
        TBranch        *b_lep1isFromW;   //!
        TBranch        *b_lep2isFromW;   //!
        TBranch        *b_mu1_numSSv3;   //!
        TBranch        *b_mu1_foSSv3;   //!
        TBranch        *b_mu1_muonidfull;   //!
        TBranch        *b_mu1_muonid;   //!
        TBranch        *b_mu1_muonidfullV1;   //!
        TBranch        *b_mu1_muonidV1;   //!
        TBranch        *b_mu1_goodmask;   //!
        TBranch        *b_mu1_gfitchi2;   //!
        TBranch        *b_mu1_cosmic;   //!
        TBranch        *b_mu1_siHits;   //!
        TBranch        *b_mu1_saHits;   //!
        TBranch        *b_mu1_emVetoDep;   //!
        TBranch        *b_mu1_hadVetoDep;   //!
        TBranch        *b_mu1_isPFmuon;   //!
        TBranch        *b_mu2_numSSv3;   //!
        TBranch        *b_mu2_foSSv3;   //!
        TBranch        *b_mu2_muonidfull;   //!
        TBranch        *b_mu2_muonid;   //!
        TBranch        *b_mu2_muonidfullV1;   //!
        TBranch        *b_mu2_muonidV1;   //!
        TBranch        *b_mu2_goodmask;   //!
        TBranch        *b_mu2_gfitchi2;   //!
        TBranch        *b_mu2_cosmic;   //!
        TBranch        *b_mu2_siHits;   //!
        TBranch        *b_mu2_saHits;   //!
        TBranch        *b_mu2_emVetoDep;   //!
        TBranch        *b_mu2_hadVetoDep;   //!
        TBranch        *b_mu2_isPFmuon;   //!
        TBranch        *b_e1_numSSv3;   //!
        TBranch        *b_e1_foSSv3;   //!
        TBranch        *b_e1_vbtf90full;   //!
        TBranch        *b_e1_vbtf90;   //!
        TBranch        *b_e1_vbtf85;   //!
        TBranch        *b_e1_vbtf80;   //!
        TBranch        *b_e1_vbtf70;   //!
        TBranch        *b_e1_smurfV3;   //!
        TBranch        *b_e1_scet;   //!
        TBranch        *b_e1_eopin;   //!
        TBranch        *b_e1_hoe;   //!
        TBranch        *b_e1_dphiin;   //!
        TBranch        *b_e1_detain;   //!
        TBranch        *b_e1_e25Me55;   //!
        TBranch        *b_e1_sigieie;   //!
        TBranch        *b_e1_eMe55;   //!
        TBranch        *b_e1_nmHits;   //!
        TBranch        *b_e1_dcot;   //!
        TBranch        *b_e1_dist;   //!
        TBranch        *b_e1_drmu;   //!
        TBranch        *b_e1_isspike;   //!
        TBranch        *b_e1_ctfCharge;   //!
        TBranch        *b_e1_gsfCharge;   //!
        TBranch        *b_e1_scCharge;   //!
        TBranch        *b_e1_fbrem;   //!
        TBranch        *b_e1_mitConv;   //!
        TBranch        *b_e2_numSSv3;   //!
        TBranch        *b_e2_foSSv3;   //!
        TBranch        *b_e2_vbtf90full;   //!
        TBranch        *b_e2_vbtf90;   //!
        TBranch        *b_e2_vbtf85;   //!
        TBranch        *b_e2_vbtf80;   //!
        TBranch        *b_e2_vbtf70;   //!
        TBranch        *b_e2_smurfV3;   //!
        TBranch        *b_e2_scet;   //!
        TBranch        *b_e2_eopin;   //!
        TBranch        *b_e2_hoe;   //!
        TBranch        *b_e2_dphiin;   //!
        TBranch        *b_e2_detain;   //!
        TBranch        *b_e2_e25Me55;   //!
        TBranch        *b_e2_sigieie;   //!
        TBranch        *b_e2_eMe55;   //!
        TBranch        *b_e2_nmHits;   //!
        TBranch        *b_e2_dcot;   //!
        TBranch        *b_e2_dist;   //!
        TBranch        *b_e2_drmu;   //!
        TBranch        *b_e2_isspike;   //!
        TBranch        *b_e2_ctfCharge;   //!
        TBranch        *b_e2_gsfCharge;   //!
        TBranch        *b_e2_scCharge;   //!
        TBranch        *b_e2_fbrem;   //!
        TBranch        *b_e2_mitConv;   //!
        TBranch        *b_trg_single_mu1;   //!
        TBranch        *b_trg_single_mu2;   //!
        TBranch        *b_trg_single_e1;   //!
        TBranch        *b_trg_single_e2;   //!
        TBranch        *b_trg_double_mu1;   //!
        TBranch        *b_trg_double_mu2;   //!
        TBranch        *b_trg_double_e1;   //!
        TBranch        *b_trg_double_e2;   //!
        TBranch        *b_trg_cross_emu;   //!
        TBranch        *b_trg_had_double_e1;   //!
        TBranch        *b_trg_had_double_e2;   //!
        TBranch        *b_trg_had_double_mu1;   //!
        TBranch        *b_trg_had_double_mu2;   //!
        TBranch        *b_trg_had_cross_emu;   //!

        BabyLooper();
        virtual ~BabyLooper();
        virtual void     Init(TTree *tree);
        virtual void     Loop(std::string sampleName, TChain *chain, unsigned int runMin=0, unsigned int runMax=9999999);                                      

        //enum DileptonHypType hyp_typeToHypType (int hyp_type);
        void setGoodRunList(std::string fname);

        void Fill(TH1F** hist, const unsigned int hyp, const float &val, const float &weight);
        void Fill2D(TH2F** hist, const unsigned int hyp, const float &valx, const float &valy, const float &weight);
        void FormatHist(TH1F** hist, std::string sampleName, std::string name, int n, float min, float max);
        void FormatHist2D(TH2F** hist, std::string sampleName, std::string name,
                int nx, float minx, float maxx, int ny, float miny, float maxy);

};

#endif

