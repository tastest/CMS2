#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TProfile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Math/LorentzVector.h"

#include "ossusy_looper.h"
#include "getMt2.C"

// #include "CORE/CMS2.cc"
// //#include "CMS2.cc" //include SParm branches
// #include "CORE/electronSelections.cc"
// #include "CORE/electronSelectionsParameters.cc"
// #include "CORE/SimpleFakeRate.cc"
// #include "CORE/mcSelections.cc"
// #include "CORE/metSelections.cc"
// #include "CORE/MT2/MT2.cc"
// #include "CORE/muonSelections.cc"
// #include "CORE/trackSelections.cc"
// #include "CORE/utilities.cc"
// #include "Tools/goodrun.cc"
// #include "CORE/ttbarSelections.cc"
// #include "CORE/jetSelections.cc"
// #include "CORE/triggerUtils.cc"

#include "CORE/CMS2.h"
#include "CORE/metSelections.h"
#include "CORE/trackSelections.h"
#include "CORE/eventSelections.h"
#include "CORE/electronSelections.h"
#include "CORE/electronSelectionsParameters.h"
#include "CORE/SimpleFakeRate.h"
#include "CORE/mcSelections.h"
#include "CORE/muonSelections.h"
#include "CORE/MT2/MT2.h"
#include "Tools/goodrun.cc"
#include "CORE/utilities.cc"
#include "histtools.h"
#include "CORE/ttbarSelections.cc"
#include "CORE/jetSelections.cc"
#include "CORE/triggerUtils.cc"


using namespace std;
using namespace tas;

//mSUGRA scan parameters-----------------------------

const int   nm0points    = 81;
const float m0min        = 0.;
const float m0max        = 4050.;
const int   nm12points   = 26;
const float m12min       = 100.;
const float m12max       = 620.;

TH1F* hsusydilPt[nm0points][nm12points]; 
TH1F* hsusytcmet[nm0points][nm12points]; 
TH2F* hsusy_met_sumjetpt[nm0points][nm12points]; 

//---------------------------------------------------

void fillUnderOverFlow(TH1F *h1, float value, float weight = 1.);
void fillUnderOverFlow(TH2F *h2, float xvalue, float yvalue, float weight = 1.);
//void fillUnderOverFlow(TProfile *h2, float xvalue, float yvalue);
void fillOverFlow(TH1F *h1, float value, float weight = 1.);
void fillOverFlow(TH2F *h2, float xvalue, float yvalue, float weight = 1.);
void fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx);
void fillHistos(TH2F *h2[4][4],float xvalue, float yvalue, float weight, int myType, int nJetsIdx);
void fillHistos(TProfile *h2[4][4],float xvalue, float yvalue,  int myType, int nJetsIdx);
float returnSigma(float sumJetPt, ossusy_looper::MetTypeEnum metType);
float returnBias(float sumJetPt, ossusy_looper::MetTypeEnum metType);

//--------------------------------------------------------------------

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

//--------------------------------------------------------------------

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

//--------------------------------------------------------------------

int getIndexFromM0(float m0){
  
  float binsize = (m0max - m0min) / (float) nm0points;
  int index     = (int)(m0 - m0min) / binsize;
  return index;
    
}

//--------------------------------------------------------------------

int getIndexFromM12(float m12){
  
  float binsize = (m12max - m12min) / (float) nm12points;
  int index     = (int)(m12 - m12min) / binsize;
  return index;
    
}

//--------------------------------------------------------------------

void ossusy_looper::makeTree(char *prefix)
{
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  //Super compressed ntuple here
  outFile   = new TFile(Form("ntp/%s_smallTree.root",prefix), "RECREATE");
  outFile->cd();
  outTree = new TTree("t","Tree");

  //Set branch addresses
  //variables must be declared in ossusy_looper.h
  outTree->Branch("weight",      &weight_,     "weight/F");
  outTree->Branch("proc",        &proc_,       "proc/I");
  outTree->Branch("leptype",     &leptype_,    "leptype/I");
  outTree->Branch("dilmass",     &dilmass_,    "dilmass/F");
  outTree->Branch("tcmet",       &tcmet_,      "tcmet/F");
  outTree->Branch("tcmetphi",    &tcmetphi_,   "tcmetphi/F");
  outTree->Branch("tcsumet",     &tcsumet_,    "tcsumet/F");
  outTree->Branch("mt2",         &mt2_,        "mt2/F");  
  outTree->Branch("mt2j",        &mt2j_,       "mt2j/F");  
  outTree->Branch("sumjetpt",    &sumjetpt_,   "sumjetpt/F");
  outTree->Branch("dileta",      &dileta_,     "dileta/F");
  outTree->Branch("dilpt",       &dilpt_,      "dilpt/F");
  outTree->Branch("dildphi",     &dildphi_,    "dildphi/F");
  outTree->Branch("njets",       &njets_,      "njets/I");
  outTree->Branch("vecjetpt",    &vecjetpt_,   "vecjetpt/F");
  outTree->Branch("pass",        &pass_,       "pass/I");
  outTree->Branch("passz",       &passz_,      "passz/I");
  outTree->Branch("m0",          &m0_,         "m0/F");
  outTree->Branch("m12",         &m12_,        "m12/F");
  outTree->Branch("ptl1",        &ptl1_,       "ptl1/F");
  outTree->Branch("ptl2",        &ptl2_,       "ptl2/F");
  outTree->Branch("ptj1",        &ptj1_,       "ptj1/F");
  outTree->Branch("ptj2",        &ptj2_,       "ptj2/F");
  outTree->Branch("meff",        &meff_,       "meff/F");
  outTree->Branch("mt",          &mt_,         "mt/F");


}

//--------------------------------------------------------------------

int getProcessType(char *prefix)
{
  int proc = -1;

  if(strcmp(prefix,"Zjets")  == 0) proc = 1;
  if(strcmp(prefix,"ttdil")  == 0) proc = 2;
  if(strcmp(prefix,"ttotr")  == 0) proc = 3;
  if(strcmp(prefix,"ww")     == 0) proc = 4;
  if(strcmp(prefix,"wz")     == 0) proc = 5;
  if(strcmp(prefix,"zz")     == 0) proc = 6;
  if(strcmp(prefix,"wjets")  == 0) proc = 7;
  if(strcmp(prefix,"tW")     == 0) proc = 8;
  if(strcmp(prefix,"LM0")    == 0) proc = 10;
  if(strcmp(prefix,"LM1")    == 0) proc = 11;
  if(strcmp(prefix,"LM2")    == 0) proc = 12;
  if(strcmp(prefix,"LM3")    == 0) proc = 13;
  if(strcmp(prefix,"LM4")    == 0) proc = 14;
  if(strcmp(prefix,"LM5")    == 0) proc = 15;
  if(strcmp(prefix,"LM6")    == 0) proc = 16;
  if(strcmp(prefix,"LM7")    == 0) proc = 17;
  if(strcmp(prefix,"LM8")    == 0) proc = 18;
  if(strcmp(prefix,"LM9")    == 0) proc = 19;
  if(strcmp(prefix,"LM10")   == 0) proc = 20;
  if(strcmp(prefix,"LM11")   == 0) proc = 21;
  if(strcmp(prefix,"LM12")   == 0) proc = 22;

  return proc;
}

//--------------------------------------------------------------------

void ossusy_looper::closeTree()
{
  outFile->cd();
  outTree->Write();
  outFile->Close();
  delete outFile;
}

//--------------------------------------------------------------------

bool nkcut(const unsigned value,const int n,
	   const int x0 =-1, const int x1 =-1, const int x2 =-1, const int x3 =-1,
	   const int x4 =-1, const int x5 =-1, const int x6 =-1, const int x7 =-1,
	   const int x8 =-1, const int x9 =-1, const int x10=-1, const int x11=-1,
	   const int x12=-1, const int x13=-1, const int x14=-1, const int x15=-1,
	   const int x16=-1, const int x17=-1, const int x18=-1, const int x19=-1,
	   const int x20=-1, const int x21=-1, const int x22=-1, const int x23=-1,
	   const int x24=-1, const int x25=-1, const int x26=-1, const int x27=-1,
	   const int x28=-1, const int x29=-1, const int x30=-1)
{
  //if (value<0) return false;
  for (int i=0;i<n;++i) {
    if (i==x0 ||i==x1 ||i==x2 ||i==x3 ||i==x4 ||
        i==x5 ||i==x6 ||i==x7 ||i==x8 ||i==x9 ||
        i==x10||i==x11||i==x12||i==x13||i==x14||
        i==x15||i==x16||i==x17||i==x18||i==x19||
        i==x20||i==x21||i==x22||i==x23||i==x24||
        i==x25||i==x26||i==x27||i==x28||i==x29||i==x30) continue;
    if (value&(1<<i)) {} else { return false;}
  }
  return true;
}

//--------------------------------------------------------------------

ossusy_looper::ossusy_looper()
{
  g_susybaseline = false;
  g_createTree   = false;
  g_useBitMask   = false;
  random3_ = new TRandom3(1);
}

//--------------------------------------------------------------------

int ossusy_looper::ScanChain(TChain* chain, char *prefix, float kFactor, int prescale,
                             JetTypeEnum jetType, MetTypeEnum metType, ZVetoEnum zveto, bool doFakeApp, bool calculateTCMET)
{

  set_goodrun_file( "Cert_TopAug30_Merged_135059-144114_recover_noESDCS_goodruns.txt" );

  SimpleFakeRate *fr_el=0;
  SimpleFakeRate *fr_mu=0;

  if(doFakeApp) {
    std::cout<<"**************************"<<std::endl;
    std::cout<<"Running FR application job"<<std::endl;
    std::cout<<"**************************"<<std::endl;
    fr_el = new SimpleFakeRate("FakeRates31May.root","eFRv215u"); 
    fr_mu = new SimpleFakeRate("FakeRates31May.root","muFR15u"); 
  }

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  BookHistos(prefix);

  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;
  // map isn't needed for this purpose, vector is sufficient
  // better would be to use a struct with run, lb, event
  map<int,int> m_events;

  // loop over files
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TChainElement* currentFile = 0;

  int nGoodEl = 0;
  int nGoodMu = 0;
  int nSkip_els_conv_dist = 0;

  /*
  // btag discriminant working points
  // { loose, medium, tight }
  float trackCountingHighEffBJetTag_wp[3]    = { 2.03,  4.38,  14.2  };
  float trackCountingHighPurBJetTag_wp[3]    = { 1.47,  2.36,  5.36  };
  float jetProbabilityBJetTag_wp[3]          = { 0.241, 0.49,  0.795 };
  float jetBProbabilityBJetTag_wp[3]         = { 1.1,   1.37,  1.39  };
  float simpleSecondaryVertexBJetTag_wp[3]   = { 1.25,  2.05,  4.07  };
  float combinedSecondaryVertexBJetTag_wp[3] = { 0.387, 0.838, 0.94  };
  */

  if(g_createTree) makeTree(prefix);

  while((currentFile = (TChainElement*)fileIter.Next())) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);

    unsigned int nEntries = tree->GetEntries();

    for(unsigned int z = 0; z < nEntries; ++z) {
      ++nEventsTotal;
      //if (!((nEventsTotal+1)%10000))
      //  cout << "Processing event " << nEventsTotal+1 << " of " << nEventsChain << " in " << prefix << endl;

      // progress feedback to user
      if (nEventsTotal % 1000 == 0){
        
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)){
                
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }

      cms2.GetEntry(z);

      bool isData = false;
      if( strcmp( prefix , "data" ) == 0 ){
        isData = true;
      }

      // skip duplicates
      if( isData ) {
 	DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
 	if (is_duplicate(id) )
 	  continue;
      }
   
      //skip events with bad els_conv_dist 
      bool skipEvent = false;
      for( unsigned int iEl = 0 ; iEl < els_conv_dist().size() ; ++iEl ){
        if( els_conv_dist().at(iEl) != els_conv_dist().at(iEl) ){
          skipEvent = true;
        }
      }
    
      if( skipEvent ){
        nSkip_els_conv_dist++;
        continue;
      }

      //goodrun list + event cleaning
      if( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;
      if( !cleaning_standardAugust2010( isData) )                    continue;

      for(unsigned int hypIdx = 0; hypIdx < hyp_p4().size(); ++hypIdx) {

        //store dilepton type in myType
        int myType = 99;
        if (hyp_type()[hypIdx] == 3) myType = 0;                          // ee
        if (hyp_type()[hypIdx] == 0) myType = 1;                          // mm
        if (hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) myType=2; // em
        if (myType == 99) {
          cout << "Skipping unknown dilepton type = " << hyp_type()[hypIdx] << endl;
          continue;
        }
      
        //trigger selection (synced with ttdil)
        bool passMu = passHLTTrigger("HLT_Mu9");
        bool runningOnMC = isData ? false : true;
        bool passEl = passEGTrigger(runningOnMC);
      
        int type = hyp_type()[hypIdx];
        if(type == 0 && !passMu)                                    continue;
        if(type == 3 && !passEl)                                    continue;
        if((type == 1 || type == 2) && !passMu && !passEl)          continue;
      
        //check that hyp leptons come from same vertex
        if(!hypsFromSameVtx(hypIdx))   continue;

        bool isNumeratorElectron_ll = false;
        bool isNumeratorElectron_lt = false;
        bool isNumeratorMuon_ll     = false;
        bool isNumeratorMuon_lt     = false;

        bool goodFakeApplicationCand = false;

        if(doFakeApp) { // only flag numerator leptons to use in FR application

          //UPDATE ME TO MATCH CURRENT ELECTRON/MUON IS!!!!!!!!!!!!!!!!!!!!!!

          //muon ID
          //           if (abs(hyp_ll_id()[hypIdx]) == 13  && (fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_ll_index()[hypIdx])))   isNumeratorMuon_ll = true;
          //           if (abs(hyp_lt_id()[hypIdx]) == 13  && (fabs(hyp_lt_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_lt_index()[hypIdx])))   isNumeratorMuon_lt = true;
          if (abs(hyp_ll_id()[hypIdx]) == 13  && (fabs(hyp_ll_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_ll_index()[hypIdx],NominalTTbar)))   isNumeratorMuon_ll = true;
          if (abs(hyp_lt_id()[hypIdx]) == 13  && (fabs(hyp_lt_p4()[hypIdx].eta()) < 2.4 && muonId(hyp_lt_index()[hypIdx],NominalTTbar)))   isNumeratorMuon_lt = true;

          //ttbar
          //           if (abs(hyp_ll_id()[hypIdx]) == 11  && ( pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_cand01 )))  isNumeratorElectron_ll = true;
          //           if (abs(hyp_lt_id()[hypIdx]) == 11  && ( pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_cand01 )))  isNumeratorElectron_lt = true;  
          if (abs(hyp_ll_id()[hypIdx]) == 11  && ( pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbar )))  isNumeratorElectron_ll = true;
          if (abs(hyp_lt_id()[hypIdx]) == 11  && ( pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbar )))  isNumeratorElectron_lt = true;  

          // selection statements for Denominator selection (missing in zeroth implementation)
          // more longer statements for better readability
          // Check that one of the leptons is FO&!Den and the other is a Num: 
          // section for fake electrons:
          if (abs(hyp_ll_id()[hypIdx]) == 11 && abs(hyp_lt_id()[hypIdx]) == 11) {
            if( 
               //                ( isFakeableElectron(cms2.hyp_ll_index()[hypIdx],el_v2_cand01)  && !isNumeratorElectron_ll && isNumeratorElectron_lt) ||
               //                ( isFakeableElectron(cms2.hyp_lt_index()[hypIdx],el_v2_cand01)  && !isNumeratorElectron_lt && isNumeratorElectron_ll)    
               ( pass_electronSelection (cms2.hyp_ll_index()[hypIdx], electronSelectionFO_el_ttbarV1_v2) && !isNumeratorElectron_ll && isNumeratorElectron_lt) ||
               ( pass_electronSelection (cms2.hyp_lt_index()[hypIdx], electronSelectionFO_el_ttbarV1_v2) && !isNumeratorElectron_lt && isNumeratorElectron_ll)    
               ) goodFakeApplicationCand = true;
            else continue;
          }
          if (abs(hyp_ll_id()[hypIdx]) == 11 && abs(hyp_lt_id()[hypIdx]) == 13) {
            if( 
               //                ( isFakeableElectron(cms2.hyp_ll_index()[hypIdx],el_v2_cand01)  && !isNumeratorElectron_ll && isNumeratorMuon_lt ) ||
               //                ( isFakeableMuon(cms2.hyp_lt_index()[hypIdx], mu_v1) && !isNumeratorMuon_lt && isNumeratorElectron_ll )
               ( pass_electronSelection (cms2.hyp_ll_index()[hypIdx], electronSelectionFO_el_ttbarV1_v2) && !isNumeratorElectron_ll && isNumeratorMuon_lt ) ||
 	       ( muonId (cms2.hyp_lt_index()[hypIdx], muonSelectionFO_mu_ttbar)                          && !isNumeratorMuon_lt && isNumeratorElectron_ll )
               ) goodFakeApplicationCand = true;
            else continue;
          }
          // section for fake muons
          if (abs(hyp_ll_id()[hypIdx]) == 13 && abs(hyp_lt_id()[hypIdx]) == 11) {
            if(
               //                ( isFakeableMuon(cms2.hyp_ll_index()[hypIdx], mu_v1) && !isNumeratorMuon_ll && isNumeratorElectron_lt ) ||
               //                ( isFakeableElectron(cms2.hyp_lt_index()[hypIdx],el_v2_cand01)  && !isNumeratorElectron_lt && isNumeratorMuon_ll )
               ( muonId (cms2.hyp_ll_index()[hypIdx], muonSelectionFO_mu_ttbar) && !isNumeratorMuon_ll   && isNumeratorElectron_lt ) ||
               ( pass_electronSelection (cms2.hyp_lt_index()[hypIdx], electronSelectionFO_el_ttbarV1_v2) && !isNumeratorElectron_lt && isNumeratorMuon_ll )
               ) goodFakeApplicationCand = true;
            else continue;
          }
          if (abs(hyp_ll_id()[hypIdx]) == 13 && abs(hyp_lt_id()[hypIdx]) == 13) {
            if(
               //                ( isFakeableMuon(cms2.hyp_ll_index()[hypIdx], mu_v1) && !isNumeratorMuon_ll && isNumeratorMuon_lt ) ||
               //                ( isFakeableMuon(cms2.hyp_lt_index()[hypIdx], mu_v1) && !isNumeratorMuon_lt && isNumeratorMuon_ll )
               ( muonId (cms2.hyp_ll_index()[hypIdx], muonSelectionFO_mu_ttbar) && !isNumeratorMuon_ll && isNumeratorMuon_lt ) ||
               ( muonId (cms2.hyp_lt_index()[hypIdx], muonSelectionFO_mu_ttbar) && !isNumeratorMuon_lt && isNumeratorMuon_ll )
               ) goodFakeApplicationCand = true;
            else continue;
          }
          //    if( goodFakeApplicationCand ) cout<< "Found a proper FR application case" << endl;

        }
        else { // default cuts
        
          //ttbarV2 muon ID
          if (abs(hyp_ll_id()[hypIdx]) == 13  && (! muonId(hyp_ll_index()[hypIdx],NominalTTbarV2)))   continue;
          if (abs(hyp_lt_id()[hypIdx]) == 13  && (! muonId(hyp_lt_index()[hypIdx],NominalTTbarV2)))   continue;
        
          //ttbarV2 electron ID
          if (abs(hyp_ll_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_ll_index()[hypIdx] , electronSelection_ttbarV2 , isData , true ))) continue;
          if (abs(hyp_lt_id()[hypIdx]) == 11  && (! pass_electronSelection( hyp_lt_index()[hypIdx] , electronSelection_ttbarV2 , isData , true ))) continue;

        }

        //OS, pt > (20,10) GeV
        if( hyp_lt_id()[hypIdx] * hyp_ll_id()[hypIdx] > 0 )                             continue;
        if( TMath::Max( hyp_ll_p4()[hypIdx].pt() , hyp_lt_p4()[hypIdx].pt() ) < 20. )   continue;
        if( TMath::Min( hyp_ll_p4()[hypIdx].pt() , hyp_lt_p4()[hypIdx].pt() ) < 10. )   continue;
      
        //for tt, check if 2 leptons are from W's
        if(strcmp(prefix,"ttdil") == 0 && ttbarconstituents(hypIdx) != 1 ) continue;
        if(strcmp(prefix,"ttotr") == 0 && ttbarconstituents(hypIdx) == 1 ) continue;

        //cout<<"pass lepton/trigger selection"<<endl;
        // check if it's a correct genp-event (deprecated)
        //std::string prefixStr(prefix);
        //if (prefixStr == "ttdil"    && genpCountPDGId(11,13,15) != 2) continue;
        //if (prefixStr == "ttotr"    && genpCountPDGId(11,13,15) == 2) continue;
        //if (prefixStr == "DYee"     && genpCountPDGId(11)       != 2) continue;
        //if (prefixStr == "DYmm"     && genpCountPDGId(13)       != 2) continue;
        //if (prefixStr == "DYtautau" && genpCountPDGId(15)       != 2) continue;


        if( hyp_p4()[hypIdx].mass() > 76 && hyp_p4()[hypIdx].mass() < 106 && hyp_ll_p4()[hypIdx].pt() > 20 && hyp_lt_p4()[hypIdx].pt() > 20 ){
          if( myType == 0 ) nGoodEl+=weight_;
          if( myType == 1 ) nGoodMu+=weight_;
        }
      
        int id_lt = hyp_lt_id()[hypIdx];
        int id_ll = hyp_ll_id()[hypIdx];

        // met and jet cuts down below...

        /*
          if (dilTruthMatch) {
          //this better be in the selections.cc
          bool isTrueLepton_ll = false;
          bool isTrueLepton_lt = false;
          isTrueLepton_ll = ( (abs(hyp_ll_id()[hypIdx]) == abs(hyp_ll_mc_id()[hypIdx]) &&
          abs(hyp_ll_mc_motherid()[hypIdx]) < 50 //I wish I could match to W or Z explicitely, not in MGraph
          )
          || (hyp_ll_mc_id()[hypIdx]==22 && 
          TMath::Abs(ROOT::Math::VectorUtil::DeltaR(hyp_ll_p4()[hypIdx],hyp_ll_mc_p4()[hypIdx])) <0.05
          && abs(hyp_ll_id()[hypIdx]) == abs(hyp_ll_mc_motherid()[hypIdx])
          )
          );
          isTrueLepton_lt = ( (abs(hyp_lt_id()[hypIdx]) == abs(hyp_lt_mc_id()[hypIdx]) &&
          abs(hyp_lt_mc_motherid()[hypIdx]) < 50 //I wish I could match to W or Z explicitely, not in MGraph
          )
          || (hyp_lt_mc_id()[hypIdx]==22 && 
          TMath::Abs(ROOT::Math::VectorUtil::DeltaR(hyp_lt_p4()[hypIdx],hyp_lt_mc_p4()[hypIdx])) <0.05
          && abs(hyp_lt_id()[hypIdx]) == abs(hyp_lt_mc_motherid()[hypIdx])
          )
          );
          if (!isTrueLepton_lt && !isTrueLepton_ll) continue;
          }
        */

        //
        // JETS
        //
        // All jets are stored corrected
        //
        // jpts
        // ===
        // NOT electron cleaned
        // abs(eta) < 999999
        // pt > 5
        // correction is applied
        //
        // hyp_jets 
        // ========
        // corrected JPT jets
        // hyp lepton cleaned
        // corrected pt > 30
        //
        // jets
        // ==============
        // siscone calojets
        // L2/L3 corrections applied
        //
        VofP4 vjets_noetacut_p4;
        VofP4 vjets_p4;
        VofP4 vjpts_p4;
        VofP4 vhyp_jets_p4;
        VofP4 vpfjets_p4;
        LorentzVector bah; 
        LorentzVector bahtot;
        vector<float> vjets_cor;

        // Count jets with corrected Pt > 30 && abs(eta) < 2.4
        // Use e/mu cleaning for all jets

        for (unsigned int ijet = 0; ijet < jets_p4().size(); ijet++) {
          LorentzVector vjet = jets_p4().at(ijet);
          LorentzVector vlt  = hyp_lt_p4()[hypIdx];
          LorentzVector vll  = hyp_ll_p4()[hypIdx];
          if (dRbetweenVectors(vjet, vll) < 0.4) continue;
          if (dRbetweenVectors(vjet, vlt) < 0.4) continue;

          bah = jets_p4().at(ijet);
          //USE CORRECTED JET PT
          if (bah.pt() * jets_cor().at(ijet) > 30. && fabs(bah.eta()) < 2.4) {
            vjets_p4.push_back(bah);
            bahtot += bah;
            vjets_cor.push_back(jets_cor().at(ijet));
          }
          if (bah.pt() * jets_cor().at(ijet) > 30.){
            vjets_noetacut_p4.push_back(bah);
          }
        }

        float vecjetpt = bahtot.pt();

        for (unsigned int ijet = 0; ijet < jpts_p4().size(); ijet++) {
          LorentzVector vjet = jpts_p4().at(ijet);
          LorentzVector vlt  = hyp_lt_p4()[hypIdx];
          LorentzVector vll  = hyp_ll_p4()[hypIdx];
          if (dRbetweenVectors(vjet, vll) < 0.4) continue;
          if (dRbetweenVectors(vjet, vlt) < 0.4) continue;

          bah = jpts_p4().at(ijet);
          if (bah.pt() > 30. && fabs(bah.eta()) < 2.4) {
            vjpts_p4.push_back(bah);
          }
        }

        for (unsigned int ijet = 0; ijet < hyp_jets_p4()[hypIdx].size(); ijet++) {
          LorentzVector vjet = hyp_jets_p4()[hypIdx].at(ijet);
          LorentzVector vlt  = hyp_lt_p4()[hypIdx];
          LorentzVector vll  = hyp_ll_p4()[hypIdx];
          if (dRbetweenVectors(vjet, vll) < 0.4) continue;
          if (dRbetweenVectors(vjet, vlt) < 0.4) continue;

          bah = hyp_jets_p4()[hypIdx].at(ijet);
          if (bah.pt() > 30. && fabs(bah.eta()) < 2.4) {
            vhyp_jets_p4.push_back(bah);
          }
        }

        for (unsigned int ijet = 0 ; ijet < pfjets_p4().size() ; ijet++) {
          
          LorentzVector vjet = pfjets_cor().at(ijet) * pfjets_p4().at(ijet);
          LorentzVector vlt  = hyp_lt_p4()[hypIdx];
          LorentzVector vll  = hyp_ll_p4()[hypIdx];
          
          if( dRbetweenVectors(vjet, vll) < 0.4 )  continue;
          if( dRbetweenVectors(vjet, vlt) < 0.4 )  continue;
          if( fabs( vjet.eta() ) > 2.5 )           continue;
          if( !passesPFJetID(ijet) )               continue;
          if( vjet.pt() < 30. )                    continue;

          vpfjets_p4.push_back( vjet );
        }

        float sumjetpt_jets_p4 = 0.;
        float sumjetpt_jpts_p4 = 0.;
        float sumjetpt_hyp_jets_p4 = 0.;
        float sumjetpt_pfjets_p4 = 0.;
        float meff_jets_p4 = 0.;
        float meff_jpts_p4 = 0.;
        float meff_hyp_jets_p4 = 0.;
        float meff_pfjets_p4 = 0.;
        for(unsigned int ijet = 0; ijet < vjets_p4.size(); ijet++) {
          //USE CORRECTED JET PT
          sumjetpt_jets_p4 += vjets_p4.at(ijet).Pt() * vjets_cor.at(ijet);
          meff_jets_p4 +=     vjets_p4.at(ijet).Pt() * vjets_cor.at(ijet);
        }
        for(unsigned int ijet = 0; ijet < vjpts_p4.size(); ijet++) {
          sumjetpt_jpts_p4 += vjpts_p4.at(ijet).Pt();
          meff_jpts_p4     += vjpts_p4.at(ijet).Pt();
        }
        for(unsigned int ijet = 0; ijet < vhyp_jets_p4.size(); ijet++) {
          sumjetpt_hyp_jets_p4 += vhyp_jets_p4.at(ijet).Pt();
          meff_hyp_jets_p4     += vhyp_jets_p4.at(ijet).Pt();
        }
        for(unsigned int ijet = 0; ijet < vpfjets_p4.size(); ijet++) {
          sumjetpt_pfjets_p4 += vpfjets_p4.at(ijet).Pt();
          meff_pfjets_p4     += vpfjets_p4.at(ijet).Pt();
        }

        float pt_lt  = hyp_lt_p4()[hypIdx].pt();
        float pt_ll  = hyp_ll_p4()[hypIdx].pt();
     
        //genmet stuff
        float genmet    = -9999;
        float gensumet  = -9999;
        float genmetphi = -9999;
     
        if( !isData ){
          genmet    = gen_met();
          gensumet  = gen_sumEt();
          genmetphi = gen_metPhi();
        }
    
        //tcmet stuff
        float tcmet    = -9999;
        float tcsumet  = -9999;
        float tcmetphi = -9999;
        
        if( calculateTCMET ) {
          //on-the-fly tcmet calculation
          metStruct tcmetStruct = correctedTCMET();
          tcmet    = tcmetStruct.met;
          tcsumet  = tcmetStruct.sumet;
          tcmetphi = tcmetStruct.metphi;
        }
        else{
          //take tcmet from event
          tcmet    = evt_tcmet();
          tcsumet  = evt_tcsumet();
          tcmetphi = evt_tcmetPhi();
        }

        //cout<<"caloMET "<<evt_met()<<" tcMET "<<tcmet<<" tcMET (new) "<<tcmetStruct.met<<endl;
        meff_jets_p4      += hyp_ll_p4()[hypIdx].Pt()+hyp_lt_p4()[hypIdx].Pt();
        meff_jets_p4      += tcmet;
        meff_jpts_p4      += hyp_ll_p4()[hypIdx].Pt()+hyp_lt_p4()[hypIdx].Pt();
        meff_jpts_p4      += tcmet;
        meff_hyp_jets_p4  += hyp_ll_p4()[hypIdx].Pt()+hyp_lt_p4()[hypIdx].Pt();
        meff_hyp_jets_p4  += tcmet;
        meff_pfjets_p4    += hyp_ll_p4()[hypIdx].Pt()+hyp_lt_p4()[hypIdx].Pt();
        meff_pfjets_p4    += tcmet;

        // choose which met type to use
        float theMet = -999999.;
        if      (metType == e_tcmet)    theMet = tcmet;
        else if (metType == e_muon)     theMet = evt_metMuonCorr();
        else if (metType == e_muonjes)  theMet = evt_metMuonJESCorr();
        else if (metType == e_pfmet)    theMet = evt_pfmet();
        else {
          std::cout << "UNRECOGNIZED METTYPE!  Learn to program..." << std::endl;
          exit(1);
        }

        // choose which jet type to use
        float theSumJetPt = -999999.;
        int theNJets      = -999999;

        if (jetType == e_JPT) {
          theSumJetPt = sumjetpt_jpts_p4;
          theNJets = vjpts_p4.size();
        } else if (jetType == e_calo) {
          theSumJetPt = sumjetpt_jets_p4;
          theNJets = vjets_p4.size();
        } else if (jetType == e_pfjet) {
          theSumJetPt = sumjetpt_pfjets_p4;
          theNJets = vpfjets_p4.size();
        }
        else {
          std::cout << "UNRECOGNIZED JETTYPE!  Learn to program..." << std::endl;
          exit(1);
        }

     
        //Probably should switch all of this stuff to pfjets!!!!!!!!!!!!!!!!!!!

        // Mt2 with all the jets, no cut on eta
        float thisMt2    = 99999.;
        float thisMet    = tcmet;
        float thisMetPhi = tcmetphi;

        // Using code in the CORE
        float mt2core = MT2(thisMet, thisMetPhi, hyp_ll_p4()[hypIdx], hyp_lt_p4()[hypIdx], 0., false);

        float mt2jcore = -1.;
        if (vjets_noetacut_p4.size() > 1)
          mt2jcore = MT2J(thisMet, thisMetPhi, hyp_ll_p4()[hypIdx], hyp_lt_p4()[hypIdx], vjets_noetacut_p4);

        // Custom Mt2
        for (unsigned int i=0; i<vjets_noetacut_p4.size(); i++) {
          LorentzVector vj1 = vjets_noetacut_p4.at(i);
          for (unsigned int j=i+1; j<vjets_noetacut_p4.size(); j++) {
            LorentzVector vj2 = vjets_noetacut_p4.at(j);

            LorentzVector v1 = vj1 + hyp_lt_p4()[hypIdx];
            LorentzVector v2 = vj2 + hyp_ll_p4()[hypIdx];
            float mt2 = getMt2(v1, v2, thisMet, thisMetPhi);
            if (mt2 < thisMt2) thisMt2 = mt2;

            v1 = vj1 + hyp_ll_p4()[hypIdx];
            v2 = vj2 + hyp_lt_p4()[hypIdx];
            mt2 = getMt2(v1, v2, thisMet, thisMetPhi);
            if (mt2 < thisMt2) thisMt2 = mt2;
          }
        }
        float mt2j = thisMt2;

        m_events.insert(pair<int,int>(evt_event(), 1));

        // The event weight including the kFactor (scaled to 100 pb-1)
        float weight = -1.;
        if(strcmp(prefix,"LMscan") == 0){
          //weight = kFactor * sparm_xsec() * 100. / 10000.; //xsec * lumi (100/pb) / nevents (10000)
          weight = 1;
          cout << "CURRENTLY NOT SET UP TO READ IN SUSY SCAN XSEC, SETTING WEIGHT = 1" << endl;
        }else if( isData ){
          weight = 1;
        }else{
          weight = kFactor * evt_scale1fb() * 0.1;
        }


        // This isn't quite right, and works only if both em and ppmux are in play
        // ibl: is this still applicable? Check! 100302
        if ((! strcmp(prefix, "ppMuX") || ! strcmp(prefix,"EM")) && 
            (hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2)) weight *= 0.5;

        int fakeRateSet = 0;

        if(doFakeApp) {
          double fr_ll = 0;
          // muon truth&ID/iso tags are NOT already in global cuts - need to check here
          //           const double eta = cms2.hyp_ll_p4()[hypIdx].eta();
          //           const double pt = cms2.hyp_ll_p4()[hypIdx].pt();
          switch (abs(cms2.hyp_ll_id()[hypIdx])) {
          case 11:
            if (
 		pass_electronSelection (cms2.hyp_ll_index()[hypIdx], electronSelectionFO_el_ttbarV1_v2) &&
 		//		isFakeableElectron(cms2.hyp_ll_index()[hypIdx],el_v2_cand01) &&
                not isNumeratorElectron_ll &&
                (isNumeratorMuon_lt || isNumeratorElectron_lt)
                ) {
 	      //              fr_ll = elFakeProb(cms2.hyp_ll_index()[hypIdx], el_v2_cand01);
              fr_ll = fr_el->getFR(hyp_ll_p4()[hypIdx].Pt(), hyp_ll_p4()[hypIdx].eta());

              fakeRateSet+=1;
            }
            break;
          case 13:
            if (
 		muonId (cms2.hyp_ll_index()[hypIdx], muonSelectionFO_mu_ttbar) &&
 		//		isFakeableMuon(cms2.hyp_ll_index()[hypIdx], mu_v1) &&
                not isNumeratorMuon_ll  &&
                (isNumeratorMuon_lt || isNumeratorElectron_lt) ) {
 	      //              fr_ll = muFakeProb(cms2.hyp_ll_index()[hypIdx], mu_v1);
              fr_ll = fr_mu->getFR(hyp_ll_p4()[hypIdx].Pt(), hyp_ll_p4()[hypIdx].eta());
              fakeRateSet+=1;
            }
            break;
          default:
            assert(0);
          }

          double fr_lt = 0;
          //           const double eta = cms2.hyp_lt_p4()[hypIdx].eta();
          //           const double pt = cms2.hyp_lt_p4()[hypIdx].pt();
          switch (abs(cms2.hyp_lt_id()[hypIdx])) {
          case 11:
            if (
 		pass_electronSelection (cms2.hyp_lt_index()[hypIdx], electronSelectionFO_el_ttbarV1_v2) &&
 		//isFakeableElectron(cms2.hyp_lt_index()[hypIdx],el_v2_cand01) &&
                not isNumeratorElectron_lt &&
                (isNumeratorMuon_ll || isNumeratorElectron_ll)
                ) {
 	      //              fr_lt = elFakeProb(cms2.hyp_lt_index()[hypIdx], el_v2_cand01);
              fr_ll = fr_el->getFR(hyp_lt_p4()[hypIdx].Pt(), hyp_lt_p4()[hypIdx].eta());
              fakeRateSet+=1;
            }
            break;
          case 13:
            if (
 		muonId (cms2.hyp_lt_index()[hypIdx], muonSelectionFO_mu_ttbar) &&
 		//		isFakeableMuon(cms2.hyp_lt_index()[hypIdx], mu_v1) &&
                not isNumeratorMuon_lt  &&
                (isNumeratorMuon_ll || isNumeratorElectron_ll) ) {
 	      //              fr_lt = muFakeProb(cms2.hyp_lt_index()[hypIdx], mu_v1);
              fr_ll = fr_mu->getFR(hyp_lt_p4()[hypIdx].Pt(), hyp_lt_p4()[hypIdx].eta());

              fakeRateSet+=1;
            }
            break;
          default:
            assert(0);
          }

          weight *= fr_lt / (1 - fr_lt) + fr_ll / (1 - fr_ll);

          // just to make sure: only one of fr_lt and fr_ll should ever be non-zero
          //  printf("fr_lt = %g\tfr_ll = %g\n", fr_lt, fr_ll);
          assert(not ((fr_lt != 0) && (fr_ll != 0)));
        }

        if(fakeRateSet>1)  cout<<"Setting FR multiple times - is this a double fake? Fake Weight was set as often as: "<<fakeRateSet<<endl;

        VofP4 *new_hyp_jets_p4 = &vhyp_jets_p4;
        VofP4 *new_jets_p4 =  &vjets_p4;
        //VofP4 *new_jpts_p4 =  &vjpts_p4;
        int new_hyp_njets =   vhyp_jets_p4.size();
        int new_njets =       vjets_p4.size();
        int new_njpts =       vjpts_p4.size();
        int nHypJetsIdx =     min(new_hyp_njets, 2);
        int nJetsIdx =        min(new_njets, 2);
        //int nJptsIdx =        min(new_njpts, 2);

        //extra variables for baby ntuple
        int pass   = ( theSumJetPt > 200 && theNJets >= 2 && theMet > 50. && id_lt * id_ll < 0 ) ? 1 : 0;
        int passz  = (passZSelection ( hypIdx ) ) ? 1 : 0;
        float etaZ = hyp_p4()[hypIdx].eta();
        float m0   = -9999.;
        float m12  = -9999.;
        
        if(strcmp(prefix,"LMscan") == 0){
          //m0  = sparm_m0();
          //m12 = sparm_m12();
          m0  = -1;
          m12 = -1;
          cout << "CURRENTLY NOT SET UP TO READ IN SUSY PARMS, SETTING M0 = M1/2 = -1" << endl;
        }

        //hyp lepton pt
        float ptll    = hyp_ll_trk_p4()[hypIdx].pt();
        float ptlt    = hyp_lt_trk_p4()[hypIdx].pt();
        float ptl1    = ( ptlt > ptll) ? ptlt : ptll; 
        float ptl2    = ( ptlt > ptll) ? ptll : ptlt; 
        float philep  = ( ptlt > ptll) ? hyp_lt_trk_p4()[hypIdx].phi() : hyp_ll_trk_p4()[hypIdx].phi() ;

        //find leading jet pT
        int   imax  = -9999;
        float ptmax = -9999;
              
        for(unsigned int ijet = 0 ; ijet < vjets_p4.size() ; ++ijet){
          if(vjets_p4.at(ijet).pt() * vjets_cor.at(ijet) > ptmax){
            ptmax = vjets_p4.at(ijet).pt() * vjets_cor.at(ijet);
            imax  = ijet;
          }
        } 

        //find 2nd leading jet pT
        int imax2    = -9999;
        float ptmax2 = -9999;
              
        for(unsigned int ijet = 0 ; ijet < vjets_p4.size() ; ++ijet){
          if(vjets_p4.at(ijet).pt() * vjets_cor.at(ijet) > ptmax2 && ijet!=imax){
            ptmax2 = vjets_p4.at(ijet).pt() * vjets_cor.at(ijet);
            imax2  = ijet;
          }
        } 

        //tranverse mass leading lepton & met
        float deltaPhiLepMet = fabs( philep - tcmetphi );
        if( deltaPhiLepMet > TMath::Pi() ) deltaPhiLepMet = TMath::TwoPi() - deltaPhiLepMet;
        float mt = sqrt( 2 * ( ptl1 * tcmet * (1 - cos(deltaPhiLepMet) ) ) );

        //delta phi btw leptons
        double dphilep = fabs(hyp_lt_p4()[hypIdx].phi() - hyp_ll_p4()[hypIdx].phi());
        if (dphilep > TMath::Pi()) dphilep = TMath::TwoPi() - dphilep;


        //fill tree for baby ntuple 
        if(g_createTree){

          weight_      = weight;                       //event weight
          proc_        = getProcessType(prefix);       //integer specifying sample
          dilmass_     = hyp_p4()[hypIdx].mass();      //dilepton mass
          dilpt_       = hyp_p4()[hypIdx].pt();        //dilepton pT
          dileta_      = hyp_p4()[hypIdx].eta();       //dilepton eta
          dildphi_     = dphilep;                      //dilepton delta phi
          tcmet_       = tcmet;                        //tcmet
          tcsumet_     = tcsumet;                      //tcsumet
          tcmetphi_    = tcmetphi;                     //tcmetphi
          sumjetpt_    = theSumJetPt;                  //scalar sum jet pt
          mt2_         = mt2core;                      //mt2 leptonic
          mt2j_        = mt2j;                         //mt2 with jets
          njets_       = theNJets;                     //njets w pt>30 and |eta|<2.4
          vecjetpt_    = vecjetpt;                     //vector sum jet pt
          pass_        = pass;                         //pass kinematic cuts
          passz_       = passz;                        //pass Z selection
          m0_          = m0;                           //mSUGRA m0
          m12_         = m12;                          //mSUGRA m1/2
          ptl1_        = ptl1;                         //highest pT lepton
          ptl2_        = ptl2;                         //2nd highest pT lepton
          ptj1_        = ptmax;                        //leading jet
          ptj2_        = ptmax2;                       //2nd leading jet
          meff_        = meff_jets_p4;                 //effective mass
          mt_          = mt;                           //transverse mass of leading lepton+met

          leptype_ = -1;
          if (hyp_type()[hypIdx] == 3) leptype_ = 0; // ee
          if (hyp_type()[hypIdx] == 0) leptype_ = 1; // mm
          if (hyp_type()[hypIdx] == 1) leptype_ = 2; // em
          if( hyp_type()[hypIdx] == 2) leptype_ = 2; // em
                
          outTree->Fill();
        }

        //selection (continue statements)-------------------------------
        if(!g_useBitMask){
        
          // all these cuts are ok for the FR application

          if (theSumJetPt < 200.)   continue;
          if (theNJets < 2)         continue;
          if (theMet < 50.)         continue;
      
          if (zveto == e_standard) {
            //veto same-flavor OS dileptons in Z mass window
            if ( ( hyp_type()[hypIdx] == 3 || hyp_type()[hypIdx] == 0 ) 
                 && hyp_p4()[hypIdx].mass() > 76. && hyp_p4()[hypIdx].mass() < 106.) continue;
          }
        
          else if(zveto == e_allzveto){
            //veto dileptons in Z mass window regardless of flavor
            if (hyp_p4()[hypIdx].mass() > 76. && hyp_p4()[hypIdx].mass() < 106.)   continue;
          }
        
          else if(zveto == e_nozveto){
            //no Z veto
          }
        
          else if(zveto == e_selectz){
            //require same-flavor OS dileptons in Z mass window
            if ( hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2 || 
                 hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.) continue;
          }
        
          else{
            cout<<"UNRECOGNIZED ZVETO"<<endl;
            exit(0);
          }
        


          fillHistos(hdilMass,  hyp_p4()[hypIdx].mass()  , weight, myType, nJetsIdx);
          fillHistos(htcmet  ,  tcmet            , weight, myType, nJetsIdx);
          fillHistos(hmetmuon,  evt_metMuonCorr(), weight, myType, nJetsIdx);
          fillHistos(hsumJetPt, sumjetpt_jets_p4 , weight, myType, nJetsIdx);
          fillHistos(hsumJptPt, sumjetpt_jpts_p4 , weight, myType, nJetsIdx);
          hnJet[myType]   ->Fill(new_njets,     weight);
          hnJet[3]        ->Fill(new_njets,     weight);
          hnJpt[myType]   ->Fill(new_njpts,     weight);
          hnJpt[3]        ->Fill(new_njpts,     weight);
          hnHypJet[myType]->Fill(new_hyp_njets, weight);
          hnHypJet[3]     ->Fill(new_hyp_njets, weight);



          //Fill susyscan histos
          if(strcmp(prefix,"LMscan") == 0){
                    
            //cout << " m0 "     << sparm_m0()  << " index " << getIndexFromM0(m0) 
            //     << " m1/2 "   << sparm_m12() << " index " << getIndexFromM12(m12)
            //     << " weight " << weight << endl;
                    
            fillUnderOverFlow( hsusydilPt[ getIndexFromM0(m0) ][ getIndexFromM12(m12)] , hyp_p4()[hypIdx].pt(), weight);
            fillUnderOverFlow( hsusytcmet[ getIndexFromM0(m0) ][ getIndexFromM12(m12)] , tcmet, weight);
            hsusy_met_sumjetpt[ getIndexFromM0(m0) ][ getIndexFromM12(m12)]->Fill(sumjetpt_jets_p4 , tcmet/sqrt(tcsumet), weight);
                    
          }
        }

        //selection (bitmask)-------------------------------------------

        if(g_useBitMask){
          const int ncut = 5;
          bool cut[ncut];
          for(int ic=0;ic<ncut;ic++)cut[ic]=false;
          unsigned cutbit=0;

          cut[0] = (id_lt * id_ll < 0); 
        
          if (zveto == e_standard) {
            //veto same-flavor OS dileptons in Z mass window
            cut[1] = (hyp_type()[hypIdx]==3 || hyp_type()[hypIdx]==0) ? 
              (hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.) : true; 
          }
          else if(zveto == e_allzveto){
            //veto OS dileptons in Z mass window regardless of flavor
            cut[1] = (hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.);
          }
          else if(zveto == e_nozveto){
            //no Z-veto
            cut[1] = true;
          }
          else if(zveto == e_selectz){
            //require same-flavor OS dileptons in Z mass window
            cut[1] = ( ( hyp_type()[hypIdx] == 3 || hyp_type()[hypIdx] == 0 ) && hyp_p4()[hypIdx].mass() > 76. && hyp_p4()[hypIdx].mass() < 106. );
          }

          cut[2] = theMet > 50.;
          cut[3] = theSumJetPt > 200.;
          cut[4] = theNJets    > 1;

          for (int icut=0;icut<ncut;++icut) {
            if (cut[icut]) cutbit+=(1<<icut);
          }

          //Fill N-1 histos
          //NOTE: if(nkcut(cutbit,ncut,X)) means: if event passes all cuts *except* cut X,
          //where X refers to cut[X] = 'some requirement' above

          if(nkcut(cutbit,ncut,1)){//dilepton pt
            fillHistos(hdilMass,  hyp_p4()[hypIdx].mass()  , weight, myType, nJetsIdx);
          }
          if(nkcut(cutbit,ncut,2)){//met
            fillHistos(htcmet  ,  tcmet            , weight, myType, nJetsIdx);
            fillHistos(hmetmuon,  evt_metMuonCorr(), weight, myType, nJetsIdx);

            fillHistos(hetaZ_tcmet,          fabs(etaZ) , tcmet,               weight, myType, nJetsIdx);
            fillHistos(hetaZ_tcmetsqrtsumet, fabs(etaZ) , tcmet/sqrt(tcsumet), weight, myType, nJetsIdx);
            fillHistos(hetaZ_tcmetsumet,     fabs(etaZ) , tcmet/tcsumet,       weight, myType, nJetsIdx);

            float x = tcmet;
            float y = hyp_p4()[hypIdx].pt();

            //full region
            fillHistos(hmt2j_all,      mt2jcore, weight, myType, nJetsIdx);
            fillHistos(hmet_dilpt_all, x, y,     weight, myType, nJetsIdx);

            //signal region
            if( x > 100 && y > 50 && x > y){
              fillHistos(hmt2j_signal,      mt2jcore, weight, myType, nJetsIdx);
              fillHistos(hmet_dilpt_signal, x, y,     weight, myType, nJetsIdx);
            }
            //control region
            if( y > 100 && x > 50 && y > x){
              fillHistos(hmt2j_control,      mt2jcore, weight, myType, nJetsIdx);
              fillHistos(hmet_dilpt_control, x, y,     weight, myType, nJetsIdx);
            }

          }
          if(nkcut(cutbit,ncut,3)){//sumjetpt
            fillHistos(hsumJetPt, sumjetpt_jets_p4 , weight, myType, nJetsIdx);
            fillHistos(hsumJptPt, sumjetpt_jpts_p4 , weight, myType, nJetsIdx);
          }

          if(nkcut(cutbit,ncut,4)){//njets
            hnJet[myType]   ->Fill(new_njets,     weight);
            hnJet[3]        ->Fill(new_njets,     weight);
            hnJpt[myType]   ->Fill(new_njpts,     weight);
            hnJpt[3]        ->Fill(new_njpts,     weight);
            hnHypJet[myType]->Fill(new_hyp_njets, weight);
            hnHypJet[3]     ->Fill(new_hyp_njets, weight);
          }
                  
          if(nkcut(cutbit,ncut,2,3)){
                    
            fillHistos(hsumJetPt_tcmet,          sumjetpt_jets_p4 , tcmet,               weight, myType, nJetsIdx);
            fillHistos(hsumJetPt_tcmetsqrtsumet, sumjetpt_jets_p4 , tcmet/sqrt(tcsumet), weight, myType, nJetsIdx);
            fillHistos(hsumJetPt_tcmetsumet,     sumjetpt_jets_p4 , tcmet/tcsumet,       weight, myType, nJetsIdx);

            fillHistos(hsumJetPt_tcmetsqrtsumet_prof, sumjetpt_jets_p4 , tcmet/sqrt(tcsumet), myType, nJetsIdx);
            fillHistos(htcsumet_tcmet_prof,           tcsumet ,          tcmet,               myType, nJetsIdx);
            fillHistos(hgensumet_genmet_prof,         gensumet   ,       genmet,              myType, nJetsIdx);

          }
                  
          //Fill susyscan histos
          if(strcmp(prefix,"LMscan") == 0){
          
            //cout << " m0 "     << sparm_m0()  << " index " << getIndexFromM0(m0) 
            //     << " m1/2 "   << sparm_m12() << " index " << getIndexFromM12(m12)
            //     << " weight " << weight << endl;
          
            if(nkcut(cutbit,ncut)){
              fillUnderOverFlow( hsusydilPt[ getIndexFromM0(m0) ][ getIndexFromM12(m12)] , hyp_p4()[hypIdx].pt(), weight);
              fillUnderOverFlow( hsusytcmet[ getIndexFromM0(m0) ][ getIndexFromM12(m12)] , tcmet, weight);
            }

            if(nkcut(cutbit,ncut,2,3)){
              hsusy_met_sumjetpt[ getIndexFromM0(m0) ][ getIndexFromM12(m12)]->Fill(sumjetpt_jets_p4 , tcmet/sqrt(tcsumet), weight);
            }
          }
        
          if(!nkcut(cutbit,ncut)) continue; //continue if event doesn't pass ALL cuts
        }
              
              
        //-------------------------------------------------------------
        // Lots of histograms
        //-------------------------------------------------------------
              
        fillHistos(hetaz, fabs(etaZ) , weight, myType, nJetsIdx);
        fillHistos(hmt2jcore, mt2jcore, weight, myType, nJetsIdx);
        fillHistos(hmt2core,  mt2core,  weight, myType, nJetsIdx);
        fillHistos(hmt2j, mt2j, weight, myType, nJetsIdx);
        fillHistos(hmt,   mt,   weight, myType, nJetsIdx);
        //fillHistos(hsumJetPt, sumjetpt_jets_p4, weight, myType, nJetsIdx);
        fillHistos(hDtcmetgenmetVsumJetPt, sumjetpt_jets_p4, tcmet-genmet, weight, myType, nJetsIdx);
        fillHistos(hDmetmuonjesgenmetVsumJetPt, sumjetpt_jets_p4, evt_metMuonJESCorr()-genmet, weight, myType, nJetsIdx);
        fillHistos(hmeffJet, meff_jets_p4, weight, myType, nJetsIdx);
        //fillHistos(hsumJptPt, sumjetpt_jpts_p4, weight, myType, nJetsIdx);
        fillHistos(hmeffJPT, meff_jpts_p4, weight, myType, nJetsIdx);
        fillHistos(hsumHypPt, sumjetpt_hyp_jets_p4, weight, myType, nHypJetsIdx);
        fillHistos(hmeffHyp, meff_hyp_jets_p4, weight, myType, nHypJetsIdx);

        //fillHistos(hetaZ_tcmet, etaZ, tcmet, weight, myType, nJetsIdx);
        //fillHistos(hetaZ_tcmetsqrtsumet, etaZ, tcmet*sqrt(theSumJetPt), weight, myType, nJetsIdx);
        //fillHistos(hetaZ_tcmetsumet, etaZ, tcmet*theSumJetPt, weight, myType, nJetsIdx);
                

        // jet count
        //hnJet[myType]->Fill(new_njets, weight);
        //hnJet[3]->Fill(new_njets, weight);
        //hnJpt[myType]->Fill(new_njpts, weight);
        //hnJpt[3]->Fill(new_njpts, weight);
        //hnHypJet[myType]->Fill(new_hyp_njets, weight);
        //hnHypJet[3]->Fill(new_hyp_njets, weight);

        // lepton Pt
        if (abs(id_lt) == 11) helePt[myType][nJetsIdx]->Fill(pt_lt, weight);
        if (abs(id_ll) == 11) helePt[myType][nJetsIdx]->Fill(pt_ll, weight);
        if (abs(id_lt) == 13) hmuPt[myType][nJetsIdx]->Fill(pt_lt, weight);
        if (abs(id_ll) == 13) hmuPt[myType][nJetsIdx]->Fill(pt_ll, weight);
        hminLepPt[myType][nJetsIdx]->Fill(min(pt_ll, pt_lt), weight);
        hmaxLepPt[myType][nJetsIdx]->Fill(max(pt_ll, pt_lt), weight );

        if (abs(id_lt) == 11) helePt[myType][3]->Fill(pt_lt, weight);
        if (abs(id_ll) == 11) helePt[myType][3]->Fill(pt_ll, weight);
        if (abs(id_lt) == 13) hmuPt[myType][3]->Fill(pt_lt, weight);
        if (abs(id_ll) == 13) hmuPt[myType][3]->Fill(pt_ll, weight);
        hminLepPt[myType][3]->Fill(min(pt_ll, pt_lt), weight);
        hmaxLepPt[myType][3]->Fill(max(pt_ll, pt_lt), weight );

        if (abs(id_lt) == 11) helePt[3][nJetsIdx]->Fill(pt_lt, weight);
        if (abs(id_ll) == 11) helePt[3][nJetsIdx]->Fill(pt_ll, weight);
        if (abs(id_lt) == 13) hmuPt[3][nJetsIdx]->Fill(pt_lt, weight);
        if (abs(id_ll) == 13) hmuPt[3][nJetsIdx]->Fill(pt_ll, weight);
        hminLepPt[3][nJetsIdx]->Fill(min(pt_ll, pt_lt), weight);
        hmaxLepPt[3][nJetsIdx]->Fill(max(pt_ll, pt_lt), weight );

        if (abs(id_lt) == 11) helePt[3][3]->Fill(pt_lt, weight);
        if (abs(id_ll) == 11) helePt[3][3]->Fill(pt_ll, weight);
        if (abs(id_lt) == 13) hmuPt[3][3]->Fill(pt_lt, weight);
        if (abs(id_ll) == 13) hmuPt[3][3]->Fill(pt_ll, weight);
        hminLepPt[3][3]->Fill(min(pt_ll, pt_lt), weight);
        hmaxLepPt[3][3]->Fill(max(pt_ll, pt_lt), weight );

        // lepton Phi
        if (abs(id_lt) == 11) helePhi[myType][nJetsIdx]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 11) helePhi[myType][nJetsIdx]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 13) hmuPhi[myType][nJetsIdx]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 13) hmuPhi[myType][nJetsIdx]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 11) helePhi[myType][3]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 11) helePhi[myType][3]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 13) hmuPhi[myType][3]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 13) hmuPhi[myType][3]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 11) helePhi[3][nJetsIdx]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 11) helePhi[3][nJetsIdx]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 13) hmuPhi[3][nJetsIdx]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 13) hmuPhi[3][nJetsIdx]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 11) helePhi[3][3]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 11) helePhi[3][3]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);
        if (abs(id_lt) == 13) hmuPhi[3][3]->Fill(hyp_lt_p4()[hypIdx].phi(), weight);
        if (abs(id_ll) == 13) hmuPhi[3][3]->Fill(hyp_ll_p4()[hypIdx].phi(), weight);

        // dilepton mass
        //fillHistos(hdilMass, hyp_p4()[hypIdx].mass(), weight, myType, nJetsIdx);

        // delta phi btw leptons
        double dphi = fabs(hyp_lt_p4()[hypIdx].phi() - hyp_ll_p4()[hypIdx].phi());
        if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
        fillHistos(hdphiLep, dphi, weight, myType, nJetsIdx);

        // dphill vs mll, i.e. the 2d correlation between the previous two variables
        fillHistos(hdphillvsmll, hyp_p4()[hypIdx].mass(), dphi, weight, myType, nJetsIdx);

        // lepton Eta
        if (abs(id_lt) == 11) heleEta[myType][nJetsIdx]->Fill(hyp_lt_p4()[hypIdx].eta(), weight);
        if (abs(id_ll) == 11) heleEta[myType][nJetsIdx]->Fill(hyp_ll_p4()[hypIdx].eta(), weight);
        if (abs(id_lt) == 13) hmuEta[myType][nJetsIdx]->Fill(hyp_lt_p4()[hypIdx].eta(), weight);
        if (abs(id_ll) == 13) hmuEta[myType][nJetsIdx]->Fill(hyp_ll_p4()[hypIdx].eta(), weight);
        if (abs(id_lt) == 11) heleEta[myType][3]->Fill(hyp_lt_p4()[hypIdx].eta(), weight);
        if (abs(id_ll) == 11) heleEta[myType][3]->Fill(hyp_ll_p4()[hypIdx].eta(), weight);
        if (abs(id_lt) == 13) hmuEta[myType][3]->Fill(hyp_lt_p4()[hypIdx].eta(), weight);
        if (abs(id_ll) == 13) hmuEta[myType][3]->Fill(hyp_ll_p4()[hypIdx].eta(), weight);
        if (abs(id_lt) == 11) heleEta[3][nJetsIdx]->Fill(hyp_lt_p4()[hypIdx].eta(), weight);
        if (abs(id_ll) == 11) heleEta[3][nJetsIdx]->Fill(hyp_ll_p4()[hypIdx].eta(), weight);
        if (abs(id_lt) == 13) hmuEta[3][nJetsIdx]->Fill(hyp_lt_p4()[hypIdx].eta(), weight);
        if (abs(id_ll) == 13) hmuEta[3][nJetsIdx]->Fill(hyp_ll_p4()[hypIdx].eta(), weight);
        if (abs(id_lt) == 11) heleEta[3][3]->Fill(hyp_lt_p4()[hypIdx].eta(), weight);
        if (abs(id_ll) == 11) heleEta[3][3]->Fill(hyp_ll_p4()[hypIdx].eta(), weight);
        if (abs(id_lt) == 13) hmuEta[3][3]->Fill(hyp_lt_p4()[hypIdx].eta(), weight);
        if (abs(id_ll) == 13) hmuEta[3][3]->Fill(hyp_ll_p4()[hypIdx].eta(), weight);

        // dilepton pt
        fillHistos(hdilPt, hyp_p4()[hypIdx].pt(), weight, myType, nJetsIdx);

        //dilepton pt with z-veto applied
        if(!passZSelection(hypIdx))
          fillHistos(hdilPt_zveto, hyp_p4()[hypIdx].pt(), weight, myType, nJetsIdx);
    
        // smeared dilepton pt
        float sigma = returnSigma(theSumJetPt, metType);
        float bias  = returnBias(theSumJetPt, metType);
        float smear = random3_->Gaus(0, sigma);
        float dilPtSmeared = bias * hyp_p4()[hypIdx].pt() + smear;
        fillHistos(hdilPtSmeared, dilPtSmeared, weight, myType, nJetsIdx);

        // Gen Met and Met Phi
        fillHistos(hgenmet, genmet, weight, myType, nJetsIdx);
        fillHistos(hgenmetPhi, genmetphi , weight, myType, nJetsIdx);

        // Met and Met phi
        //fillHistos(hmetmuon, evt_metMuonCorr(), weight, myType, nJetsIdx);
        fillHistos(hmetmuonPhi, evt_metMuonCorrPhi(), weight, myType, nJetsIdx);

        // pat Met and Met phi
        fillHistos(hmetmuonjes, evt_metMuonJESCorr(), weight, myType, nJetsIdx);
        fillHistos(hmetmuonjesPhi, evt_metMuonJESCorrPhi(), weight, myType, nJetsIdx);

        // tc Met and Met phi
        //fillHistos(htcmet, tcmet, weight, myType, nJetsIdx);
        fillHistos(htcmetPhi, tcmetphi, weight, myType, nJetsIdx);

        // pf Met and Met phi
        fillHistos(hpfmet, evt_pfmet(), weight, myType, nJetsIdx);
        fillHistos(hpfmetPhi, evt_pfmetPhi(), weight, myType, nJetsIdx);

        // Met vs dilepton Pt
        fillHistos(hmetmuonVsDilepPt, evt_metMuonCorr(), hyp_p4()[hypIdx].pt(), weight, myType, nJetsIdx);

        // pat  Met vs dilepton Pt
        fillHistos(hmetmuonjesVsDilepPt, evt_metMuonJESCorr(), hyp_p4()[hypIdx].pt(), weight, myType, nJetsIdx);

        // tc  Met vs dilepton Pt
        fillHistos(htcmetVsDilepPt, tcmet, hyp_p4()[hypIdx].pt(), weight, myType, nJetsIdx);

        // Met over dilepton Pt vs deltaphi btw the two
        double dphi2 = fabs(hyp_p4()[hypIdx].phi() - evt_metMuonCorrPhi() );
        if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
        dphi2 = TMath::Pi() - dphi2;  // changed the definition CC 28 March 08
        fillHistos(hmetmuonOverPtVsDphi, evt_metMuonCorr()/hyp_p4()[hypIdx].pt(), dphi2, weight, myType, nJetsIdx);

        // pat Met over dilepton Pt vs deltaphi btw the two
        dphi2 = fabs(hyp_p4()[hypIdx].phi() - evt_metMuonJESCorrPhi());
        if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
        dphi2 = TMath::Pi() - dphi2;  // changed the definition CC 28 March 08
        fillHistos(hmetmuonjesOverPtVsDphi, evt_metMuonJESCorr()/hyp_p4()[hypIdx].pt(), dphi2, weight, myType, nJetsIdx);

        // tc Met over dilepton Pt vs deltaphi btw the two
        dphi2 = fabs(hyp_p4()[hypIdx].phi() - tcmetphi);
        if (dphi2 > TMath::Pi()) dphi2 = TMath::TwoPi() - dphi2;
        dphi2 = TMath::Pi() - dphi2;  // changed the definition CC 28 March 08
        fillHistos(htcmetOverPtVsDphi, tcmet/hyp_p4()[hypIdx].pt(), dphi2, weight, myType, nJetsIdx);

        // Make a vector of jets sorted by pt and fill jet histograms
        if (new_njets > 0) {
          //vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > my_jets_p4(*new_jets_p4);
          vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > my_jets_p4(*new_jets_p4);
          sort(my_jets_p4.begin(), my_jets_p4.end(), sortByPt);

          fillHistos(hptJet1, my_jets_p4[0].Pt(), weight, myType, nJetsIdx);
          fillHistos(hetaJet1, my_jets_p4[0].Eta(), weight, myType, nJetsIdx);

          if (new_njets > 1) {
            fillHistos(hptJet2, my_jets_p4[1].Pt(), weight, myType, nJetsIdx);
            fillHistos(hetaJet2, my_jets_p4[1].Eta(), weight, myType, nJetsIdx);
          }
          if (new_njets > 2) {
            fillHistos(hptJet3, my_jets_p4[2].Pt(), weight, myType, nJetsIdx);
            fillHistos(hetaJet3, my_jets_p4[2].Eta(), weight, myType, nJetsIdx);
          }
          if (new_njets > 3) {
            fillHistos(hptJet4, my_jets_p4[3].Pt(), weight, myType, nJetsIdx);
            fillHistos(hetaJet4, my_jets_p4[3].Eta(), weight, myType, nJetsIdx);
          }
        }

        // Make a vector of hyp jets sorted by pt and fill jet histograms
        if (new_hyp_njets > 0) {
          //vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > my_hyp_jets_p4(*new_hyp_jets_p4);
          vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > my_hyp_jets_p4(*new_hyp_jets_p4);
          sort(my_hyp_jets_p4.begin(), my_hyp_jets_p4.end(), sortByPt);

          fillHistos(hptHypJet1, my_hyp_jets_p4[0].Pt(), weight, myType, nHypJetsIdx);
          fillHistos(hetaHypJet1, my_hyp_jets_p4[0].Eta(), weight, myType, nHypJetsIdx);

          if (new_hyp_njets > 1) {
            fillHistos(hptHypJet2, my_hyp_jets_p4[1].Pt(), weight, myType, nHypJetsIdx);
            fillHistos(hetaHypJet2, my_hyp_jets_p4[1].Eta(), weight, myType, nHypJetsIdx);
          }
          if (new_hyp_njets > 2) {
            fillHistos(hptHypJet3, my_hyp_jets_p4[2].Pt(), weight, myType, nHypJetsIdx);
            fillHistos(hetaHypJet3, my_hyp_jets_p4[2].Eta(), weight, myType, nHypJetsIdx);
          }
          if (new_hyp_njets > 3) {
            fillHistos(hptHypJet4, my_hyp_jets_p4[3].Pt(), weight, myType, nHypJetsIdx);
            fillHistos(hetaHypJet4, my_hyp_jets_p4[3].Eta(), weight, myType, nHypJetsIdx);
          }
        }

        // Make a vector of jets sorted by pt and fill jet histograms
        if (new_njpts > 0) {
          //vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > my_jpts_p4(*new_jets_p4);
          vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > my_jpts_p4(*new_jets_p4);
          sort(my_jpts_p4.begin(), my_jpts_p4.end(), sortByPt);

          fillHistos(hptJpt1, my_jpts_p4[0].Pt(), weight, myType, nJetsIdx);
          fillHistos(hetaJpt1, my_jpts_p4[0].Eta(), weight, myType, nJetsIdx);

          if (new_njets > 1) {
            fillHistos(hptJpt2, my_jpts_p4[1].Pt(), weight, myType, nJetsIdx);
            fillHistos(hetaJpt2, my_jpts_p4[1].Eta(), weight, myType, nJetsIdx);
          }
          if (new_njets > 2) {
            fillHistos(hptJpt3, my_jpts_p4[2].Pt(), weight, myType, nJetsIdx);
            fillHistos(hetaJpt3, my_jpts_p4[2].Eta(), weight, myType, nJetsIdx);
          }
          if (new_njets > 3) {
            fillHistos(hptJpt4, my_jpts_p4[3].Pt(), weight, myType, nJetsIdx);
            fillHistos(hetaJpt4, my_jpts_p4[3].Eta(), weight, myType, nJetsIdx);
          }
        }
      }
    } // entries
  } // currentFile

  cout << nGoodEl << " Zee events" << endl;
  cout << nGoodMu << " Zmm events" << endl;
  if( nSkip_els_conv_dist > 0 )
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist" << endl;
  
  if(g_createTree) closeTree();
  
  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  
  return 0;
}

//--------------------------------------------------------------------

bool ossusy_looper::passZSelection(int hypIdx)
{
  //require OS leptons, same flavor, in Z mass window
  if(hyp_lt_id()[hypIdx] * hyp_ll_id()[hypIdx] > 0 )                   return false;
  if(hyp_type()[hypIdx]==1 || hyp_type()[hypIdx]==2)                   return false;
  if(hyp_p4()[hypIdx].mass() < 76. || hyp_p4()[hypIdx].mass() > 106.)  return false;

  return true;
}

//--------------------------------------------------------------------

bool ossusy_looper::passTrigger(int dilType)
{
  //  bool hlt_ele15_lw_l1r = passHLTTrigger("HLT_Ele15_SW_L1R");
  bool hlt_ele15_lw_l1r = passHLTTrigger("HLT_Ele15_LW_L1R");
  bool hltMu9           = passHLTTrigger("HLT_Mu9");

  if (dilType == 0 && ! (hltMu9) ) return false;
  if ((dilType == 1 || dilType == 2) && ! (hltMu9 || hlt_ele15_lw_l1r)) return false;
  if (dilType == 3 && ! hlt_ele15_lw_l1r) return false;

  return true;
}

//--------------------------------------------------------------------
 
void ossusy_looper::BookHistos(char *prefix)
{
  // Prefix comes from the sample and it is passed to the scanning function
  // Suffix is "ee" "em" "em" "all" which depends on the final state
  // For example: histogram named tt_hnJet_ee would be the Njet distribution
  // for the ee final state in the ttbar sample.
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  cout << "Begin book histos..." << endl;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  char jetbins[5][7]    = {"0", "1", "2", "3", "#geq 4"};
  char suffixall[4][4]  = {"ee", "mm", "em", "all"};
  char njetCh[4][5]     = {"0j", "1j", "2j", "allj"};
  //char btagpoints[3][7] = {"loose", "medium", "tight"};

  //double binedges1000[11] = {0.,  20.,  40.,  60.,  80., 100., 120., 150.,  200.,  300.,  500.};
  //double binedges1500[11] = {0., 100., 200., 300., 400., 500., 600., 800., 1000., 1200., 1500.};
  double binedges1500[6] = {0., 100., 200., 400., 800., 1500.};
  //double binedges2000[11] = {0., 100., 200., 300., 400., 500., 600., 800., 1000., 1500., 2000.};

  if(strcmp("LMscan",prefix)==0){
    for(int im0 = 0 ; im0 < nm0points ; im0++){
      for(int im12 = 0 ; im12 < nm12points ; im12++){
          
        hsusydilPt[im0][im12] = new TH1F(Form("susy_hdilPt_m0_%i_m12_%i",im0,im12),
                                         Form("susy_hdilPt_m0_%i_m12_%i",im0,im12),60,0.,300.);
        hsusydilPt[im0][im12] -> GetXaxis()->SetTitle("Pt (GeV)");
          
        hsusytcmet[im0][im12] = new TH1F(Form("susy_htcmet_m0_%i_m12_%i",im0,im12),
                                         Form("susy_htcmet_m0_%i_m12_%i",im0,im12),60,0.,300.);
        hsusytcmet[im0][im12] -> GetXaxis()->SetTitle("tcmet (GeV)");

        hsusy_met_sumjetpt[im0][im12] = new TH2F(Form("susy_hmet_sumjetpt_m0_%i_m12_%i",im0,im12),
                                                 Form("susy_hmet_sumjetpt_m0_%i_m12_%i",im0,im12),200,0,2000,500,0,50);

        hsusy_met_sumjetpt[im0][im12] -> GetXaxis()->SetTitle("sumJetPt (GeV)");
        hsusy_met_sumjetpt[im0][im12] -> GetYaxis()->SetTitle("tcmet / #sqrt{tcsumet} (GeV^{1/2})");
      }
    }
  }


  for (int i = 0; i < 4; i++) {
    hnJet[i] = new TH1F(Form("%s_hnJet_%s",prefix,suffixall[i]),Form("%s_nJet_%s",prefix,suffixall[i]),11,-0.5,10.5);	
    hnJet[i]->GetXaxis()->SetTitle("nJets");

    hnJpt[i] = new TH1F(Form("%s_hnJpt_%s",prefix,suffixall[i]),Form("%s_nJpt_%s",prefix,suffixall[i]),11,-0.5,10.5);	
    hnJpt[i]->GetXaxis()->SetTitle("nJpts");

    hnHypJet[i] = new TH1F(Form("%s_hnHypJet_%s",prefix,suffixall[i]),Form("%s_nHypJet_%s",prefix,suffixall[i]),11,-0.5,10.5);	
    hnHypJet[i]->GetXaxis()->SetTitle("nHypJets");

    for(int k = 0; k < 5; k++) {
      hnJet[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
      hnJet[i]->GetXaxis()->SetLabelSize(0.07);

      hnJpt[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
      hnJpt[i]->GetXaxis()->SetLabelSize(0.07);

      hnHypJet[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
      hnHypJet[i]->GetXaxis()->SetLabelSize(0.07);
    }

    for (int j = 0; j < 4; j++) {
      char suffix[7];
      sprintf(suffix, "%s_%s", njetCh[j], suffixall[i]);

      hetaz[i][j]   = new TH1F(Form("%s_hetaz_%s",prefix,suffix),
                               Form("%s_etaz_%s" ,prefix,suffix),1000,0,5);
            
      hmt[i][j]   = new TH1F(Form("%s_hmt_%s",prefix,suffix),
                             Form("%s_mt_%s" ,prefix,suffix),1000,0,5000);
            
      hmt2core[i][j] = new TH1F(Form("%s_hmt2core_%s",prefix,suffix),
                                Form("%s_mt2core_%s" ,prefix,suffix),1000,0,1000);
  
      hmt2jcore[i][j] = new TH1F(Form("%s_hmt2jcore_%s",prefix,suffix),
                                 Form("%s_mt2jcore_%s" ,prefix,suffix),1000,0,1000);
            
      hmt2j[i][j] = new TH1F(Form("%s_hmt2j_%s",prefix,suffix),
                             Form("%s_mt2j_%s" ,prefix,suffix),1000,0,1000);

      hmet_dilpt_all[i][j] = new TH2F(Form("%s_met_dilpt_all_%s",prefix,suffix),
                                      Form("%s_met_dilpt_all_%s" ,prefix,suffix),500,0,500,500,0,500);

      hmet_dilpt_signal[i][j] = new TH2F(Form("%s_met_dilpt_signal_%s",prefix,suffix),
                                         Form("%s_met_dilpt_signal_%s" ,prefix,suffix),500,0,500,500,0,500);

      hmet_dilpt_control[i][j] = new TH2F(Form("%s_met_dilpt_control_%s",prefix,suffix),
                                          Form("%s_met_dilpt_control_%s" ,prefix,suffix),500,0,500,500,0,500);
 
      hmt2j_all[i][j] = new TH1F(Form("%s_hmt2j_all_%s",prefix,suffix),
                                 Form("%s_mt2j_all_%s" ,prefix,suffix),1000,0,1000);

      hmt2j_signal[i][j] = new TH1F(Form("%s_hmt2j_signal_%s",prefix,suffix),
                                    Form("%s_mt2j_signal_%s" ,prefix,suffix),1000,0,1000);

      hmt2j_control[i][j] = new TH1F(Form("%s_hmt2j_control_%s",prefix,suffix),
                                     Form("%s_mt2j_control_%s" ,prefix,suffix),1000,0,1000);
           
      hdilMass_tcmet[i][j] = new TH2F(Form("%s_dilMass_tcmet_%s",prefix,suffix),
                                      Form("%s_dilMass_tcmet_%s",prefix,suffix),500,0,500,500,0,500);
            
      hsumJetPt_tcmet[i][j] = new TH2F(Form("%s_sumJetPt_tcmet_%s",prefix,suffix),
                                       Form("%s_sumJetPt_tcmet_%s",prefix,suffix),200,0,2000,500,0,500);
            
      hsumJetPt_tcmetsqrtsumet[i][j] = new TH2F(Form("%s_sumJetPt_tcmetsqrtsumet_%s",prefix,suffix),
                                                Form("%s_sumJetPt_tcmetsqrtsumet_%s",prefix,suffix),200,0,2000,500,0,50);

      hsumJetPt_tcmetsqrtsumet_prof[i][j] = new TProfile(Form("%s_sumJetPt_tcmetsqrtsumet_prof_%s",prefix,suffix),
                                                         Form("%s_sumJetPt_tcmetsqrtsumet_prof_%s",prefix,suffix),200,0,2000,0,50);
      
      hsumJetPt_tcmetsqrtsumet_prof[i][j]->GetXaxis()->SetTitle("sumJetPt (GeV)");
      hsumJetPt_tcmetsqrtsumet_prof[i][j]->GetYaxis()->SetTitle("tcmet / #sqrt{tcsumet} (GeV^{1/2})");

      htcsumet_tcmet_prof[i][j] = new TProfile(Form("%s_tcsumet_tcmet_prof_%s",prefix,suffix),
                                               Form("%s_tcsumet_tcmet_prof_%s",prefix,suffix),200,0,2000,0,500);
      
      htcsumet_tcmet_prof[i][j]->GetXaxis()->SetTitle("tcsumet (GeV)");
      htcsumet_tcmet_prof[i][j]->GetYaxis()->SetTitle("tcmet (GeV)");
      
      hsumJetPt_tcmetsumet[i][j] = new TH2F(Form("%s_sumJetPt_tcmetsumet_%s",prefix,suffix),
                                            Form("%s_sumJetPt_tcmetsumet_%s",prefix,suffix),200,0,2000,500,0,1);
            
      hetaZ_tcmet[i][j] = new TH2F(Form("%s_etaZ_tcmet_%s",prefix,suffix),
                                   Form("%s_etaZ_tcmet_%s",prefix,suffix),200,0,5,500,0,500);
            
      hetaZ_tcmetsqrtsumet[i][j] = new TH2F(Form("%s_etaZ_tcmetsqrtsumet_%s",prefix,suffix),
                                            Form("%s_etaZ_tcmetsqrtsumet_%s",prefix,suffix),200,0,5,500,0,50);
            
      hetaZ_tcmetsumet[i][j] = new TH2F(Form("%s_etaZ_tcmetsumet_%s",prefix,suffix),
                                        Form("%s_etaZ_tcmetsumet_%s",prefix,suffix),200,0,5,500,0,1);

      hgensumet_genmet_prof[i][j] = new TProfile(Form("%s_gensumet_genmet_prof_%s",prefix,suffix),
                                                 Form("%s_gensumet_genmet_prof_%s",prefix,suffix),200,0,2000,0,500);
            
      //hsumJetPt[i][j] = new TH1F(Form("%s_hsumJetPt_%s",prefix,suffix),Form("%s_sumJetPt_%s",prefix,suffix),5,binedges1500);
      hsumJetPt[i][j] = new TH1F(Form("%s_hsumJetPt_%s",prefix,suffix),Form("%s_sumJetPt_%s",prefix,suffix),200,0,2000);
      //hmeffJet[i][j] = new TH1F(Form("%s_hmeffJet_%s",prefix,suffix),Form("%s_meffJet_%s",prefix,suffix),5,binedges1500);
      hmeffJet[i][j] = new TH1F(Form("%s_hmeffJet_%s",prefix,suffix),Form("%s_meffJet_%s",prefix,suffix),200,0,2000);
      hDtcmetgenmetVsumJetPt[i][j] = new TH2F(Form("%s_hDtcmetgenmetVsumJetPt_%s",prefix,suffix),Form("%s_DtcmetgenmetVsumJetPt_%s",prefix,suffix),5,binedges1500,40,-50.,50.);
      hDmetmuonjesgenmetVsumJetPt[i][j] = new TH2F(Form("%s_hDmetmuonjesgenmetVsumJetPt_%s",prefix,suffix),Form("%s_DmetmuonjesgenmetVsumJetPt_%s",prefix,suffix),5,binedges1500,40,-50.,50.);

      hsumJptPt[i][j] = new TH1F(Form("%s_hsumJptPt_%s",prefix,suffix),Form("%s_sumJptPt_%s",prefix,suffix),5,binedges1500);
      hmeffJPT[i][j] = new TH1F(Form("%s_hmeffJPT_%s",prefix,suffix),Form("%s_meffJPT_%s",prefix,suffix),5,binedges1500);

      hsumHypPt[i][j] = new TH1F(Form("%s_hsumHypPt_%s",prefix,suffix),Form("%s_sumHypPt_%s",prefix,suffix),5,binedges1500);
      hmeffHyp[i][j] = new TH1F(Form("%s_hmeffHyp_%s",prefix,suffix),Form("%s_meffJHyp_%s",prefix,suffix),5,binedges1500);

      //
      // HISTOGRAMS INHERITED FROM SLAVA'S TTDIL LOOPER
      //

      helePt[i][j] = new TH1F(Form("%s_helePt_%s",prefix,suffix),Form("%s_elePt_%s",prefix,suffix),60,0.,300.);
      helePt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hmuPt[i][j]  = new TH1F(Form("%s_hmuPt_%s",prefix,suffix),Form("%s_muPt_%s",prefix,suffix),60,0.,300.);
      hmuPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hminLepPt[i][j]  = new TH1F(Form("%s_hminLepPt_%s",prefix,suffix),Form("%s_minLepPt_%s",prefix,suffix),60,0.,300.);
      hminLepPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hmaxLepPt[i][j]  = new TH1F(Form("%s_hmaxLepPt_%s",prefix,suffix),Form("%s_maxLepPt_%s",prefix,suffix),60,0.,300.);
      hmaxLepPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      helePhi[i][j] = new TH1F(Form("%s_helePhi_%s",prefix,suffix),Form("%s_elePhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      helePhi[i][j]->GetXaxis()->SetTitle("#phi");

      hmuPhi[i][j]  = new TH1F(Form("%s_hmuPhi_%s",prefix,suffix),Form("%s_muPhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      hmuPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hdphiLep[i][j]  = new TH1F(Form("%s_hdphiLep_%s",prefix,suffix),Form("%s_dphiLep_%s",prefix,suffix),50,0.,TMath::Pi());
      hdphiLep[i][j]->GetXaxis()->SetTitle("#delta#phi_{ll}");

      heleEta[i][j] = new TH1F(Form("%s_heleEta_%s",prefix,suffix),Form("%s_eleEta_%s",prefix,suffix),60,-3.,3.);
      heleEta[i][j]->GetXaxis()->SetTitle("#eta");

      hmuEta[i][j]  = new TH1F(Form("%s_hmuEta_%s",prefix,suffix),Form("%s_muEta_%s",prefix,suffix),60,-3.,3.);
      hmuEta[i][j]->GetXaxis()->SetTitle("#eta");

      hdilMass[i][j] = new TH1F(Form("%s_hdilMass_%s",prefix,suffix),Form("%s_dilMass_%s",prefix,suffix),60,0.,300.);
      hdilMass[i][j]->GetXaxis()->SetTitle("Mass_{ll} (GeV)");

      hdilPt[i][j] = new TH1F(Form("%s_hdilPt_%s",prefix,suffix),Form("%s_dilPt_%s",prefix,suffix),60,0.,300.);
      hdilPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      //dilepton pT with z-veto applied
      hdilPt_zveto[i][j] = new TH1F(Form("%s_hdilPt_zveto_%s",prefix,suffix),Form("%s_dilPt_zveto_%s",prefix,suffix),60,0.,300.);
      hdilPt_zveto[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hdilPtSmeared[i][j] = new TH1F(Form("%s_hdilPtSmeared_%s",prefix,suffix),Form("%s_dilPtSmeared_%s",prefix,suffix),60,0.,300.);
      hdilPtSmeared[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      // changed binning from 2 GeV to 10 GeV
      hgenmet[i][j] = new TH1F(Form("%s_hgenmet_%s",prefix,suffix),Form("%s_genmet_%s",prefix,suffix),60,0.,300.);
      hgenmet[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      hgenmetPhi[i][j] = new TH1F(Form("%s_hgenmetPhi_%s",prefix,suffix),Form("%s_genmetPhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      hgenmetPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hpfmet[i][j] = new TH1F(Form("%s_hpfmet_%s",prefix,suffix),Form("%s_pfmet_%s",prefix,suffix),60,0.,300.);
      hpfmet[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      hpfmetPhi[i][j] = new TH1F(Form("%s_hpfmetPhi_%s",prefix,suffix),Form("%s_pfmetPhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      hpfmetPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hmetmuon[i][j] = new TH1F(Form("%s_hmetmuon_%s",prefix,suffix),Form("%s_metmuon_%s",prefix,suffix),60,0.,300.);
      hmetmuon[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      hmetmuonPhi[i][j] = new TH1F(Form("%s_hmetmuonPhi_%s",prefix,suffix),Form("%s_metmuonPhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      hmetmuonPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hmetmuonVsDilepPt[i][j] = new TH2F(Form("%s_hmetmuonVsDilepPt_%s",prefix,suffix),Form("%s_metmuonVsDilepPt_%s",prefix,suffix),60,0.,300.,60,0.,300.);
      hmetmuonVsDilepPt[i][j]->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
      hmetmuonVsDilepPt[i][j]->GetYaxis()->SetTitle("Met (GeV)");

      hmetmuonOverPtVsDphi[i][j] = new TH2F(Form("%s_hmetmuonOverPtVsDphi_%s",prefix,suffix),Form("%s_metmuonOverPtVsDphi_%s",prefix,suffix),30,0.,3.,25,0.,TMath::Pi());
      hmetmuonVsDilepPt[i][j]->GetXaxis()->SetTitle("#Delta#Phi");
      hmetmuonVsDilepPt[i][j]->GetYaxis()->SetTitle("MET/Pt_{ll}");

      // met with muon corr and jes
      hmetmuonjes[i][j] = new TH1F(Form("%s_hmetmuonjes_%s",prefix,suffix),Form("%s_metmuonjes_%s",prefix,suffix),60,0.,300.);
      hmetmuonjes[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      hmetmuonjesPhi[i][j] = new TH1F(Form("%s_hmetmuonjesPhi_%s",prefix,suffix),Form("%s_metmuonjesPhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      hmetmuonjesPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hmetmuonjesVsDilepPt[i][j] = new TH2F(Form("%s_hmetmuonjesVsDilepPt_%s",prefix,suffix),Form("%s_metmuonjesVsDilepPt_%s",prefix,suffix),60,0.,300.,60,0.,300.);
      hmetmuonjesVsDilepPt[i][j]->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
      hmetmuonjesVsDilepPt[i][j]->GetYaxis()->SetTitle("Met (GeV)");

      hmetmuonjesOverPtVsDphi[i][j] = new TH2F(Form("%s_hmetmuonjesOverPtVsDphi_%s",prefix,suffix),Form("%s_metmuonjesOverPtVsDphi_%s",prefix,suffix),30,0.,3.,25,0.,TMath::Pi());
      hmetmuonjesVsDilepPt[i][j]->GetXaxis()->SetTitle("#Delta#Phi");
      hmetmuonjesVsDilepPt[i][j]->GetYaxis()->SetTitle("MET/Pt_{ll}");

      // tc
      htcmet[i][j] = new TH1F(Form("%s_htcmet_%s",prefix,suffix),Form("%s_tcmet_%s",prefix,suffix),60,0.,300.);
      htcmet[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      htcmetPhi[i][j] = new TH1F(Form("%s_htcmetPhi_%s",prefix,suffix),Form("%s_tcmetPhi_%s",prefix,suffix),50,-1*TMath::Pi(),TMath::Pi());
      htcmetPhi[i][j]->GetXaxis()->SetTitle("#phi");

      htcmetVsDilepPt[i][j] = new TH2F(Form("%s_htcmetVsDilepPt_%s",prefix,suffix),Form("%s_tcmetVsDilepPt_%s",prefix,suffix),60,0.,300.,60,0.,300.);
      htcmetVsDilepPt[i][j]->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
      htcmetVsDilepPt[i][j]->GetYaxis()->SetTitle("Met (GeV)");

      htcmetOverPtVsDphi[i][j] = new TH2F(Form("%s_htcmetOverPtVsDphi_%s",prefix,suffix),Form("%s_tcmetOverPtVsDphi_%s",prefix,suffix),30,0.,3.,25,0.,TMath::Pi());
      htcmetOverPtVsDphi[i][j]->GetXaxis()->SetTitle("#Delta#Phi");
      htcmetOverPtVsDphi[i][j]->GetYaxis()->SetTitle("MET/Pt_{ll}");

      hdphillvsmll[i][j] = new TH2F(Form("%s_dphillvsmll_%s",prefix,suffix),Form("%s_dphillvsmll_%s",prefix,suffix),100,10.,210.,50,0.,TMath::Pi());
      hdphillvsmll[i][j]->GetXaxis()->SetTitle("Mass_{ll} (GeV)");
      hdphillvsmll[i][j]->GetYaxis()->SetTitle("#delta#phi_{ll}");

      hptJet1[i][j] = new TH1F(Form("%s_hptJet1_%s",prefix,suffix),Form("%s_ptJet1_%s",prefix,suffix),60,0.,300.);
      hptJet1[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJet2[i][j] = new TH1F(Form("%s_hptJet2_%s",prefix,suffix),Form("%s_ptJet2_%s",prefix,suffix),60,0.,300.);
      hptJet2[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJet3[i][j] = new TH1F(Form("%s_hptJet3_%s",prefix,suffix),Form("%s_ptJet3_%s",prefix,suffix),60,0.,300.);
      hptJet3[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJet4[i][j] = new TH1F(Form("%s_hptJet4_%s",prefix,suffix),Form("%s_ptJet4_%s",prefix,suffix),60,0.,300.);
      hptJet4[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hetaJet1[i][j] = new TH1F(Form("%s_hetaJet1_%s",prefix,suffix),Form("%s_etaJet1_%s",prefix,suffix),50,-4.,4.);
      hetaJet1[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJet2[i][j] = new TH1F(Form("%s_hetaJet2_%s",prefix,suffix),Form("%s_etaJet2_%s",prefix,suffix),50,-4.,4.);
      hetaJet2[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJet3[i][j] = new TH1F(Form("%s_hetaJet3_%s",prefix,suffix),Form("%s_etaJet3_%s",prefix,suffix),50,-4.,4.);
      hetaJet3[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJet4[i][j] = new TH1F(Form("%s_hetaJet4_%s",prefix,suffix),Form("%s_etaJet4_%s",prefix,suffix),50,-4.,4.);
      hetaJet4[i][j]->GetXaxis()->SetTitle("#eta");

      hptJpt1[i][j] = new TH1F(Form("%s_hptJpt1_%s",prefix,suffix),Form("%s_ptJpt1_%s",prefix,suffix),60,0.,300.);
      hptJpt1[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJpt2[i][j] = new TH1F(Form("%s_hptJpt2_%s",prefix,suffix),Form("%s_ptJpt2_%s",prefix,suffix),60,0.,300.);
      hptJpt2[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJpt3[i][j] = new TH1F(Form("%s_hptJpt3_%s",prefix,suffix),Form("%s_ptJpt3_%s",prefix,suffix),60,0.,300.);
      hptJpt3[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJpt4[i][j] = new TH1F(Form("%s_hptJpt4_%s",prefix,suffix),Form("%s_ptJpt4_%s",prefix,suffix),60,0.,300.);
      hptJpt4[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hetaJpt1[i][j] = new TH1F(Form("%s_hetaJpt1_%s",prefix,suffix),Form("%s_etaJpt1_%s",prefix,suffix),50,-4.,4.);
      hetaJpt1[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJpt2[i][j] = new TH1F(Form("%s_hetaJpt2_%s",prefix,suffix),Form("%s_etaJpt2_%s",prefix,suffix),50,-4.,4.);
      hetaJpt2[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJpt3[i][j] = new TH1F(Form("%s_hetaJpt3_%s",prefix,suffix),Form("%s_etaJpt3_%s",prefix,suffix),50,-4.,4.);
      hetaJpt3[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJpt4[i][j] = new TH1F(Form("%s_hetaJpt4_%s",prefix,suffix),Form("%s_etaJpt4_%s",prefix,suffix),50,-4.,4.);
      hetaJpt4[i][j]->GetXaxis()->SetTitle("#eta");

      hptHypJet1[i][j] = new TH1F(Form("%s_hptHypJet1_%s",prefix,suffix),Form("%s_ptHypJet1_%s",prefix,suffix),60,0.,300.);
      hptHypJet1[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptHypJet2[i][j] = new TH1F(Form("%s_hptHypJet2_%s",prefix,suffix),Form("%s_ptHypJet2_%s",prefix,suffix),60,0.,300.);
      hptHypJet2[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptHypJet3[i][j] = new TH1F(Form("%s_hptHypJet3_%s",prefix,suffix),Form("%s_ptHypJet3_%s",prefix,suffix),60,0.,300.);
      hptHypJet3[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptHypJet4[i][j] = new TH1F(Form("%s_hptHypJet4_%s",prefix,suffix),Form("%s_ptHypJet4_%s",prefix,suffix),60,0.,300.);
      hptHypJet4[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hetaHypJet1[i][j] = new TH1F(Form("%s_hetaHypJet1_%s",prefix,suffix),Form("%s_etaHypJet1_%s",prefix,suffix),50,-4.,4.);
      hetaHypJet1[i][j]->GetXaxis()->SetTitle("#eta");

      hetaHypJet2[i][j] = new TH1F(Form("%s_hetaHypJet2_%s",prefix,suffix),Form("%s_etaHypJet2_%s",prefix,suffix),50,-4.,4.);
      hetaHypJet2[i][j]->GetXaxis()->SetTitle("#eta");

      hetaHypJet3[i][j] = new TH1F(Form("%s_hetaHypJet3_%s",prefix,suffix),Form("%s_etaHypJet3_%s",prefix,suffix),50,-4.,4.);
      hetaHypJet3[i][j]->GetXaxis()->SetTitle("#eta");

      hetaHypJet4[i][j] = new TH1F(Form("%s_hetaHypJet4_%s",prefix,suffix),Form("%s_etaHypJet4_%s",prefix,suffix),50,-4.,4.);
      hetaHypJet4[i][j]->GetXaxis()->SetTitle("#eta");

      /*************************/

      hsumJetPt[i][j]->Sumw2(); 
      hmeffJet[i][j]->Sumw2(); 

      hsumJptPt[i][j]->Sumw2(); 
      hmeffJPT[i][j]->Sumw2(); 

      hsumHypPt[i][j]->Sumw2(); 
      hmeffHyp[i][j]->Sumw2(); 

      if (j == 0)  {
        hnJet[i]->Sumw2();
        hnJpt[i]->Sumw2();
        hnHypJet[i]->Sumw2();
      }
      helePt[i][j]->Sumw2();
      hmuPt[i][j]->Sumw2();
      hminLepPt[i][j]->Sumw2();
      hmaxLepPt[i][j]->Sumw2();
      helePhi[i][j]->Sumw2();
      hmuPhi[i][j]->Sumw2();
      hdphiLep[i][j]->Sumw2();
      heleEta[i][j]->Sumw2();
      hmuEta[i][j]->Sumw2();
      hdilMass[i][j]->Sumw2();
      hdilPt[i][j]->Sumw2();
      hdilPt_zveto[i][j]->Sumw2();
      hdilPtSmeared[i][j]->Sumw2();
      hgenmet[i][j]->Sumw2();
      hgenmetPhi[i][j]->Sumw2();
      hmetmuon[i][j]->Sumw2();
      hmetmuonPhi[i][j]->Sumw2();
      hmetmuonjes[i][j]->Sumw2();
      hmetmuonjesPhi[i][j]->Sumw2();
      htcmet[i][j]->Sumw2();
      htcmetPhi[i][j]->Sumw2();
      hpfmet[i][j]->Sumw2();
      hpfmetPhi[i][j]->Sumw2();
      hptJet1[i][j]->Sumw2();
      hptJet2[i][j]->Sumw2();
      hptJet3[i][j]->Sumw2();
      hptJet4[i][j]->Sumw2();
      hetaJet1[i][j]->Sumw2();
      hetaJet2[i][j]->Sumw2();
      hetaJet3[i][j]->Sumw2();
      hetaJet4[i][j]->Sumw2();
      hptJpt1[i][j]->Sumw2();
      hptJpt2[i][j]->Sumw2();
      hptJpt3[i][j]->Sumw2();
      hptJpt4[i][j]->Sumw2();
      hetaJpt1[i][j]->Sumw2();
      hetaJpt2[i][j]->Sumw2();
      hetaJpt3[i][j]->Sumw2();
      hetaJpt4[i][j]->Sumw2();
      hptHypJet1[i][j]->Sumw2();
      hptHypJet2[i][j]->Sumw2();
      hptHypJet3[i][j]->Sumw2();
      hptHypJet4[i][j]->Sumw2();
      hetaHypJet1[i][j]->Sumw2();
      hetaHypJet2[i][j]->Sumw2();
      hetaHypJet3[i][j]->Sumw2();
      hetaHypJet4[i][j]->Sumw2();
    } // njet loop
  } // channel loop

  cout << "End book histos..." << endl;
}// CMS2::BookHistos()

//--------------------------------------------------------------------

void fillUnderOverFlow(TH1F *h1, float value, float weight)
{
  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void fillUnderOverFlow(TH2F *h2, float xvalue, float yvalue, float weight)
{
  float maxx = h2->GetXaxis()->GetXmax();
  float minx = h2->GetXaxis()->GetXmin();
  float maxy = h2->GetYaxis()->GetXmax();
  float miny = h2->GetYaxis()->GetXmin();

  if (xvalue > maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
  if (xvalue < minx) xvalue = h2->GetXaxis()->GetBinCenter(1);
  if (yvalue > maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());
  if (yvalue < miny) yvalue = h2->GetYaxis()->GetBinCenter(1);

  h2->Fill(xvalue, yvalue, weight);
}

//--------------------------------------------------------------------

// void fillUnderOverFlow(TProfile *h2, float xvalue, float yvalue)
// {
//   float maxx = h2->GetXaxis()->GetXmax();
//   float minx = h2->GetXaxis()->GetXmin();
//   float maxy = h2->GetYaxis()->GetXmax();
//   float miny = h2->GetYaxis()->GetXmin();

//   if (xvalue > maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
//   if (xvalue < minx) xvalue = h2->GetXaxis()->GetBinCenter(1);
//   if (yvalue > maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());
//   if (yvalue < miny) yvalue = h2->GetYaxis()->GetBinCenter(1);

//   h2->Fill(xvalue, yvalue);
// }

//--------------------------------------------------------------------

void fillOverFlow(TH1F *h1, float value, float weight)
{
  float max = h1->GetXaxis()->GetXmax();
  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void fillOverFlow(TH2F *h2, float xvalue, float yvalue, float weight)
{
  float maxx = h2->GetXaxis()->GetXmax();
  float maxy = h2->GetYaxis()->GetXmax();

  if (xvalue > maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
  if (yvalue > maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());

  h2->Fill(xvalue, yvalue, weight);
}

//--------------------------------------------------------------------

void fillHistos(TH1F *h1[4][4],float value, float weight, int myType, int nJetsIdx)
{
  fillUnderOverFlow(h1[myType][nJetsIdx], value, weight);      
  fillUnderOverFlow(h1[myType][3],        value, weight);      
  fillUnderOverFlow(h1[3][nJetsIdx],      value, weight);      
  fillUnderOverFlow(h1[3][3],             value, weight);      
}

//--------------------------------------------------------------------

void fillHistos(TH2F *h2[4][4],float xvalue, float yvalue, float weight, int myType, int nJetsIdx)
{
  fillUnderOverFlow(h2[myType][nJetsIdx], xvalue, yvalue, weight);      
  fillUnderOverFlow(h2[myType][3],        xvalue, yvalue, weight);      
  fillUnderOverFlow(h2[3][nJetsIdx],      xvalue, yvalue, weight);      
  fillUnderOverFlow(h2[3][3],             xvalue, yvalue, weight);      
}

//--------------------------------------------------------------------

void fillHistos(TProfile *h2[4][4],float xvalue, float yvalue, int myType, int nJetsIdx)
{
  h2[myType][nJetsIdx] -> Fill(xvalue, yvalue);      
  h2[myType][3]        -> Fill(xvalue, yvalue);      
  h2[3][nJetsIdx]      -> Fill(xvalue, yvalue);      
  h2[3][3]             -> Fill(xvalue, yvalue);      
}

//--------------------------------------------------------------------

float returnSigma (float sumJetPt, ossusy_looper::MetTypeEnum metType)
{
  const float xbins[] = { 100, 200, 400, 800 };
  const float tcmet_sigma[] = { 13.02, 13.44, 14.64, 17.19, 21.95 };
  const float muonjes_sigma[] = { 21.13, 20.3, 21.08, 24.06, 29.08 };

  if (metType == ossusy_looper::e_tcmet) {
    for( int j = 0; j < 4; j++ ) {
      if( sumJetPt < xbins[j] )
        return tcmet_sigma[j];
    }

    return tcmet_sigma[4];
  }
  else if (metType == ossusy_looper::e_muonjes) {
    for( int j = 0; j < 4; j++ ) {
      if( sumJetPt < xbins[j] )
        return muonjes_sigma[j];
    }

    return muonjes_sigma[4];
  }
}

//--------------------------------------------------------------------

float returnBias(float sumJetPt, ossusy_looper::MetTypeEnum metType)
{
  const float xbins[] = { 100, 200, 400, 800 };
  const float tcmet_bias[]  = { 0.87, 0.88, 0.91, 0.95, 0.99 };
  const float muonjes_bias[]  = { 0.94, 1., 1.03, 1.03, 1.07 };

  if (metType == ossusy_looper::e_tcmet) {
    for( int j = 0; j < 4; j++ ) {
      if( sumJetPt < xbins[j] )
        return tcmet_bias[j];
    }

    return tcmet_bias[4];
  }
  else if (metType == ossusy_looper::e_muonjes) {
    for( int j = 0; j < 4; j++ ) {
      if( sumJetPt < xbins[j] )
        return muonjes_bias[j];
    }

    return muonjes_bias[4];
  }
}

//--------------------------------------------------------------------
