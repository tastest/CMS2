//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 30 04:58:46 2010 by ROOT version 5.22/00d
// from TTree tree/A Baby Ntuple
// found on file: validate_mus_before.root
//////////////////////////////////////////////////////////

#ifndef xsecLoop_h
#define xsecLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <iomanip>
#include <iostream>
#include <stdio.h>

using std::vector;

class xsecLoop {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nel;
   vector<bool>    *elsid;
   vector<float>   *elsd0corr;
   vector<int>     *elscharge;
   Int_t           elsp4_;
   Float_t         elsp4_fCoordinates_fX[123];   //[_]
   Float_t         elsp4_fCoordinates_fY[123];   //[_]
   Float_t         elsp4_fCoordinates_fZ[123];   //[_]
   Float_t         elsp4_fCoordinates_fT[123];   //[_]
   vector<float>   *elsiso;
   vector<float>   *elstrkiso;
   vector<float>   *elsecliso;
   vector<float>   *elshcliso;
   Int_t           nmu;
   vector<bool>    *musid;
   vector<float>   *musd0corr;
   vector<int>     *muscharge;
   Int_t           musp4_;
   Float_t         musp4_fCoordinates_fX[123];   //[_]
   Float_t         musp4_fCoordinates_fY[123];   //[_]
   Float_t         musp4_fCoordinates_fZ[123];   //[_]
   Float_t         musp4_fCoordinates_fT[123];   //[_]
   vector<float>   *musiso;
   vector<float>   *mustrkiso;
   vector<float>   *musecliso;
   vector<float>   *mushcliso;
   Float_t         clmet;
   Float_t         tcmet;
   Float_t         pfmet;
   Float_t         clmetphi;
   Float_t         tcmetphi;
   Float_t         pfmetphi;
   UInt_t          run;
   UInt_t          event;
   UInt_t          lumi;
   Int_t           pfjets_;
   Float_t         pfjets_fCoordinates_fX[123];   //[_]
   Float_t         pfjets_fCoordinates_fY[123];   //[_]
   Float_t         pfjets_fCoordinates_fZ[123];   //[_]
   Float_t         pfjets_fCoordinates_fT[123];   //[_]
   Int_t           trkjets_;
   Float_t         trkjets_fCoordinates_fX[123];   //[_]
   Float_t         trkjets_fCoordinates_fY[123];   //[_]
   Float_t         trkjets_fCoordinates_fZ[123];   //[_]
   Float_t         trkjets_fCoordinates_fT[123];   //[_]
   Int_t           jets_;
   Float_t         jets_fCoordinates_fX[123];   //[_]
   Float_t         jets_fCoordinates_fY[123];   //[_]
   Float_t         jets_fCoordinates_fZ[123];   //[_]
   Float_t         jets_fCoordinates_fT[123];   //[_]

   // List of branches
   TBranch        *b_nel;   //!
   TBranch        *b_elsid;   //!
   TBranch        *b_elsd0corr;   //!
   TBranch        *b_elscharge;   //!
   TBranch        *b_elsp4;   //!
   TBranch        *b_elsp4_fCoordinates_fX;   //!
   TBranch        *b_elsp4_fCoordinates_fY;   //!
   TBranch        *b_elsp4_fCoordinates_fZ;   //!
   TBranch        *b_elsp4_fCoordinates_fT;   //!
   TBranch        *b_elsiso;   //!
   TBranch        *b_elstrkiso;   //!
   TBranch        *b_elsecliso;   //!
   TBranch        *b_elshcliso;   //!
   TBranch        *b_nmu;   //!
   TBranch        *b_musid;   //!
   TBranch        *b_musd0corr;   //!
   TBranch        *b_muscharge;   //!
   TBranch        *b_musp4;   //!
   TBranch        *b_musp4_fCoordinates_fX;   //!
   TBranch        *b_musp4_fCoordinates_fY;   //!
   TBranch        *b_musp4_fCoordinates_fZ;   //!
   TBranch        *b_musp4_fCoordinates_fT;   //!
   TBranch        *b_musiso;   //!
   TBranch        *b_mustrkiso;   //!
   TBranch        *b_musecliso;   //!
   TBranch        *b_mushcliso;   //!
   TBranch        *b_clmet;   //!
   TBranch        *b_tcmet;   //!
   TBranch        *b_pfmet;   //!
   TBranch        *b_clmetphi;   //!
   TBranch        *b_tcmetphi;   //!
   TBranch        *b_pfmetphi;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_pfjets;   //!
   TBranch        *b_pfjets_fCoordinates_fX;   //!
   TBranch        *b_pfjets_fCoordinates_fY;   //!
   TBranch        *b_pfjets_fCoordinates_fZ;   //!
   TBranch        *b_pfjets_fCoordinates_fT;   //!
   TBranch        *b_trkjets;   //!
   TBranch        *b_trkjets_fCoordinates_fX;   //!
   TBranch        *b_trkjets_fCoordinates_fY;   //!
   TBranch        *b_trkjets_fCoordinates_fZ;   //!
   TBranch        *b_trkjets_fCoordinates_fT;   //!
   TBranch        *b_jets;   //!
   TBranch        *b_jets_fCoordinates_fX;   //!
   TBranch        *b_jets_fCoordinates_fY;   //!
   TBranch        *b_jets_fCoordinates_fZ;   //!
   TBranch        *b_jets_fCoordinates_fT;   //!

   xsecLoop(TTree *tree=0);
   virtual ~xsecLoop();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::vector<int>& numLepVsRun, Bool_t doMu);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   //ibl/jmu
   virtual void     readFile(TString, std::vector<double>&);
};

#endif

#ifdef xsecLoop_cxx
xsecLoop::xsecLoop(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("validate_mus_before.root");
      if (!f) {
         f = new TFile("validate_mus_before.root");
      }
      tree = (TTree*)gDirectory->Get("tree");

   }
   Init(tree);
}

xsecLoop::~xsecLoop()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t xsecLoop::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t xsecLoop::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void xsecLoop::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   elsid = 0;
   elsd0corr = 0;
   elscharge = 0;
   elsiso = 0;
   elstrkiso = 0;
   elsecliso = 0;
   elshcliso = 0;
   musid = 0;
   musd0corr = 0;
   muscharge = 0;
   musiso = 0;
   mustrkiso = 0;
   musecliso = 0;
   mushcliso = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nel", &nel, &b_nel);
   fChain->SetBranchAddress("elsid", &elsid, &b_elsid);
   fChain->SetBranchAddress("elsd0corr", &elsd0corr, &b_elsd0corr);
   fChain->SetBranchAddress("elscharge", &elscharge, &b_elscharge);
   fChain->SetBranchAddress("elsp4", &elsp4_, &b_elsp4);
   fChain->SetBranchAddress("elsp4.fCoordinates.fX", &elsp4_fCoordinates_fX, &b_elsp4_fCoordinates_fX);
   fChain->SetBranchAddress("elsp4.fCoordinates.fY", &elsp4_fCoordinates_fY, &b_elsp4_fCoordinates_fY);
   fChain->SetBranchAddress("elsp4.fCoordinates.fZ", &elsp4_fCoordinates_fZ, &b_elsp4_fCoordinates_fZ);
   fChain->SetBranchAddress("elsp4.fCoordinates.fT", &elsp4_fCoordinates_fT, &b_elsp4_fCoordinates_fT);
   fChain->SetBranchAddress("elsiso", &elsiso, &b_elsiso);
   fChain->SetBranchAddress("elstrkiso", &elstrkiso, &b_elstrkiso);
   fChain->SetBranchAddress("elsecliso", &elsecliso, &b_elsecliso);
   fChain->SetBranchAddress("elshcliso", &elshcliso, &b_elshcliso);
   fChain->SetBranchAddress("nmu", &nmu, &b_nmu);
   fChain->SetBranchAddress("musid", &musid, &b_musid);
   fChain->SetBranchAddress("musd0corr", &musd0corr, &b_musd0corr);
   fChain->SetBranchAddress("muscharge", &muscharge, &b_muscharge);
   fChain->SetBranchAddress("musp4", &musp4_, &b_musp4);
   fChain->SetBranchAddress("musp4.fCoordinates.fX", musp4_fCoordinates_fX, &b_musp4_fCoordinates_fX);
   fChain->SetBranchAddress("musp4.fCoordinates.fY", musp4_fCoordinates_fY, &b_musp4_fCoordinates_fY);
   fChain->SetBranchAddress("musp4.fCoordinates.fZ", musp4_fCoordinates_fZ, &b_musp4_fCoordinates_fZ);
   fChain->SetBranchAddress("musp4.fCoordinates.fT", musp4_fCoordinates_fT, &b_musp4_fCoordinates_fT);
   fChain->SetBranchAddress("musiso", &musiso, &b_musiso);
   fChain->SetBranchAddress("mustrkiso", &mustrkiso, &b_mustrkiso);
   fChain->SetBranchAddress("musecliso", &musecliso, &b_musecliso);
   fChain->SetBranchAddress("mushcliso", &mushcliso, &b_mushcliso);
   fChain->SetBranchAddress("clmet", &clmet, &b_clmet);
   fChain->SetBranchAddress("tcmet", &tcmet, &b_tcmet);
   fChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
   fChain->SetBranchAddress("clmetphi", &clmetphi, &b_clmetphi);
   fChain->SetBranchAddress("tcmetphi", &tcmetphi, &b_tcmetphi);
   fChain->SetBranchAddress("pfmetphi", &pfmetphi, &b_pfmetphi);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("pfjets", &pfjets_, &b_pfjets);
   fChain->SetBranchAddress("pfjets.fCoordinates.fX", pfjets_fCoordinates_fX, &b_pfjets_fCoordinates_fX);
   fChain->SetBranchAddress("pfjets.fCoordinates.fY", pfjets_fCoordinates_fY, &b_pfjets_fCoordinates_fY);
   fChain->SetBranchAddress("pfjets.fCoordinates.fZ", pfjets_fCoordinates_fZ, &b_pfjets_fCoordinates_fZ);
   fChain->SetBranchAddress("pfjets.fCoordinates.fT", pfjets_fCoordinates_fT, &b_pfjets_fCoordinates_fT);
   fChain->SetBranchAddress("trkjets", &trkjets_, &b_trkjets);
   fChain->SetBranchAddress("trkjets.fCoordinates.fX", trkjets_fCoordinates_fX, &b_trkjets_fCoordinates_fX);
   fChain->SetBranchAddress("trkjets.fCoordinates.fY", trkjets_fCoordinates_fY, &b_trkjets_fCoordinates_fY);
   fChain->SetBranchAddress("trkjets.fCoordinates.fZ", trkjets_fCoordinates_fZ, &b_trkjets_fCoordinates_fZ);
   fChain->SetBranchAddress("trkjets.fCoordinates.fT", trkjets_fCoordinates_fT, &b_trkjets_fCoordinates_fT);
   fChain->SetBranchAddress("jets", &jets_, &b_jets);
   fChain->SetBranchAddress("jets.fCoordinates.fX", jets_fCoordinates_fX, &b_jets_fCoordinates_fX);
   fChain->SetBranchAddress("jets.fCoordinates.fY", jets_fCoordinates_fY, &b_jets_fCoordinates_fY);
   fChain->SetBranchAddress("jets.fCoordinates.fZ", jets_fCoordinates_fZ, &b_jets_fCoordinates_fZ);
   fChain->SetBranchAddress("jets.fCoordinates.fT", jets_fCoordinates_fT, &b_jets_fCoordinates_fT);
   Notify();
}

Bool_t xsecLoop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void xsecLoop::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t xsecLoop::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef xsecLoop_cxx


#ifndef __CINT__
template<typename T, typename T2>
std::vector<float> rollingThunder (const std::vector<T> &v, const std::vector<T2> &vLumi, int n_comb)
{
  std::vector<float> ret(v.size(), 0);
  for (unsigned int i = 0; i < v.size(); ++i) {
	  // skip over empty entries
	  if (vLumi[i] == 0) {
      ret[i] = 0;
      continue;
	  }
	  T tot = 0;
	  int n_back = 0;
	  // go backward n_comb steps, skipping over empty entries, but not going past 0
	  for (int j = i - 1; j >= 0 && n_back < n_comb; --j) {
      if (vLumi[j] == 0)
		    continue;
      tot += v[j];
      n_back++;
	  }
	  int n_fwd = 0;
	  // go forward n_comb steps, skipping over empty entries, but not going past v.size()
	  for (unsigned int j = i + 1; j < v.size() && n_fwd < n_comb; ++j) {
      if (vLumi[j] == 0)
		    continue;
      tot += v[j];
      n_fwd++;
	  }
	  // and count the present entry as well
	  tot += v[i];
	  ret[i] = tot / float(n_back + n_fwd + 1);
  }
  return ret;
}


#endif

