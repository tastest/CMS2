#ifndef ossusy_looper_gmsb_h
#define ossusy_looper_gmsb_h

#include <vector>
#include <map>
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"
typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;
typedef map<unsigned int, unsigned int> m_uiui;

class  TChain;
class  TH1F;
class  TH2F;
class  TRandom3;
class  TTree;
struct metStruct;

class ossusy_looper_gmsb
{
    public: 
        ossusy_looper_gmsb();
        ~ossusy_looper_gmsb() {}

        enum JetTypeEnum { e_JPT = 0, e_calo };
        // e_JPT     :   jpt jets
        // e_calo    :   l1 and l2 corrected calo jets
        enum MetTypeEnum { e_tcmet = 0, e_muon, e_muonjes };
        // e_tcmet   :   track corrected met
        // e_muon    :   calo met with muon corrections
        // e_muonjes :   calo met with muon and jet energy scale corrections
        enum ZVetoEnum   { e_standard = 0, e_allzveto, e_nozveto, e_selectz };
        // e_standard:   apply Z-veto to same-flavor pairs
        // e_allzveto:   apply Z-veto regardless of lepton flavor
        // e_nozveto :   no Z-veto
        // e_selectz :   select Z by requiring SF OS pair in Z mass window

        int  ScanChain(TChain *chain, char *prefix = "", float kFactor = 1., int prescale = 1., 
                       JetTypeEnum jetType = e_JPT, 
                       MetTypeEnum metType = e_tcmet,
                       ZVetoEnum zveto = e_standard,
                       bool doFakeApp = false,
                       bool calculateTCMET = false);
        void BookHistos (char *prefix);
      
        // Set globals
        void set_susybaseline (bool b) { g_susybaseline = b; }
        void set_createTree   (bool b) { g_createTree   = b; }
        void set_useBitMask   (bool b) { g_useBitMask   = b; }

        // Baby ntuple methods
        void makeTree (char *prefix);
        void closeTree ();



    private:

        // Globals
        bool g_susybaseline;
        bool g_createTree;
        bool g_useBitMask;
        TRandom3 *random3_;

        TFile* outFile;
        TTree* outTree;

        TH1F* hnleptons_all_pass;
        TH1F* hnleptons_ss_pass;
        TH1F* hnleptons_os_pass;
        TH1F* hsign;
        TH1F* hnleptons;
        TH1F* hnelectrons;
        TH1F* hnmuons;
        TH1F* hmet;
        TH1F* hnjets;
        TH1F* hsumjetpt;
        TH1F* hpt;
        TH1F* hpt1;
        TH1F* hpt2;
        TH1F* hpt3;
        TH1F* hpt4;
        TH1F* hdr_all;
        TH1F* hdilmass_all;
        TH1F* hdr_ptmax;
        TH1F* hdilmass_ptmax;
   
};

#endif
