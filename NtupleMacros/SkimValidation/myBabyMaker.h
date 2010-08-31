#ifndef myBabyMaker_h
#define myBabyMaker_h

#include "TFile.h"
#include "TTree.h"

class myBabyMaker
{
public:
  myBabyMaker();
  ~myBabyMaker() {
	delete babyFile_;
	delete babyTree_;
  };
  void MakeBabyNtuple(const char *);
  void InitBabyNtuple();
  void FillBabyNtuple();
  void CloseBabyNtuple();
  void ScanChain(TChain *, const char *, const char* GoodRunFile, const bool doels, const bool domus);
  vector<int> doElectrons();
  vector<int> doMuons();
  void doJets();
  void doMet();
  vector<LorentzVector> CleanJets(const vector<LorentzVector> vect_p4_jets,
								  const vector<int> elidxs,
								  const vector<int> muidxs,
								  const bool docor = false, //only calo
								  const float jet_pt_threshold = 30.0,
								  const float jet_eta_threshold = 2.40,
								  const float jet_lepton_dR_veto_cone = 0.4
								  );

private:

  //TRIGGERS
  const string mutrig_;
  const string eltrig_;

  //cuts
  //const float minelpt_ = 10.; //this is in the standard selection already
  //const float minmupt_ = 10.; //check this

  // BABY NTUPLE VARIABLES
  TFile *babyFile_;
  TTree *babyTree_;

  //trigger/dataset decisions--no need for both params and data members
  //bool doels;
  //bool domus;

  //ints
  unsigned int run_;
  unsigned int event_;
  unsigned int lumi_;

  Int_t   nJets_;
  Int_t   nels_;
  Int_t   nmus_;
  vector<int>   elscharge_; //0 or 1 for now (0 if disagree, 1 if agree), but if it becomes a bitmask, it'll be int
  vector<int>   muscharge_; //0 or 1 for now (0 if disagree, 1 if agree), but if it becomes a bitmask, it'll be int

  // floats--els
  vector<float> elsd0corr_;
  vector<float> elsreliso_;
  vector<float> elstrkiso_;
  vector<float> elsecliso_;
  vector<float> elshcliso_;

  // floats--mus
  vector<float> musd0corr_;
  vector<float> musreliso_;
  vector<float> mustrkiso_;
  vector<float> musecliso_;
  vector<float> mushcliso_;
  vector<bool>  musId_;
  vector<bool>  elsId_;

  // Jets
  vector<LorentzVector> jets_;
  vector<LorentzVector> pfjets_;
  vector<LorentzVector> trkjets_;
  //vector<LorentzVector> hypjets_;

  // floats--met (no vectors)
  Float_t clmet_;
  Float_t tcmet_;
  Float_t pfmet_;
  Float_t clmetphi_;
  Float_t tcmetphi_;
  Float_t pfmetphi_;

  // Lorentz Vectors--jets
  vector<LorentzVector> musp4_;
  vector<LorentzVector> elsp4_;  

};

#endif
