
#ifndef MYSCANCHAIN_H
#define MYSCANCHAIN_H

// C++ includes
#include <iostream>
#include <vector>

#include "TROOT.h"
#include "../../CORE/CMS2.h"

class TChain;
class TDirectory;
class TProfile;
class TFile;
class TTree;
class LikelihoodUtil;

class MyScanChain {

    public:

        MyScanChain() {};
        ~MyScanChain() {};

        int ScanChain(bool isData, std::string sampleName, TChain *chain, float kFactor = 1.0, int nEvents = -1, std::string skimFilePrefix="");
        void setGoodRunList(std::string fname);

    private:

        // baby ntuple
        // init baby branches with dummy values
        void InitBaby();
        TFile *babyFile_;
        TTree *babyTree_;

        Int_t run_;
        Int_t ls_;
        Int_t evt_;
        Float_t weight_;
        Int_t type_;

        Float_t reco_ptlt_ ;
        Float_t reco_etalt_;
        Float_t reco_philt_;
        Float_t reco_isolt_;
        Int_t reco_typelt_;
        Int_t reco_mctypelt_;
        Int_t reco_algolt_;

        Float_t reco_ptll_ ;
        Float_t reco_etall_;
        Float_t reco_phill_;
        Float_t reco_isoll_;
        Int_t reco_typell_;
        Int_t reco_mctypell_;
        Int_t reco_algoll_;

        Float_t reco_mdil_;
        Float_t reco_vdilpt_;
        Float_t reco_vdilphi_;

        // pf based jets and met
        Float_t reco_tcmet_;
        Float_t reco_pfmet_;
        Int_t reco_npfjets_;

        // what passed
        Float_t reco_lhlt_;
        Float_t reco_lhll_;
        Float_t reco_pfmvalt_;
        Float_t reco_pfmvall_;

        //
        // likelihood util
        //

        LikelihoodUtil *likelihoodUtil_;

        //
        // analysis functions
        //

        void AnalyseDilepton(const float &weight);

        // is it running on real data
        bool isData_;

};

#endif

