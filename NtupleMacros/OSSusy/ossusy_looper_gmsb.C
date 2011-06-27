

//-----------------------------------------------------------------------------------------------------
// Updates since last revision
// 1) Update to Spring10 MC samples
// 2) Update to new implementation of electron ID
//-----------------------------------------------------------------------------------------------------

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

#include "ossusy_looper_gmsb.h"
#include "getMt2.C"
//#include "CORE/CMS2.cc"
#include "CMS2.cc" //include SParm branches
#include "CORE/electronSelections.cc"
#include "CORE/electronSelectionsParameters.cc"
#include "CORE/SimpleFakeRate.cc"
#include "CORE/mcSelections.cc"
#include "CORE/metSelections.cc"
#include "CORE/MT2/MT2.cc"
#include "CORE/muonSelections.cc"
#include "CORE/trackSelections.cc"
#include "CORE/utilities.cc"

using namespace std;
using namespace tas;

//---------------------------------------------------

void fillUnderOverFlow(TH1F *h1, float value, float weight = 1.);

void ossusy_looper_gmsb::makeTree(char *prefix)
{
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  //Super compressed ntuple here
  outFile   = new TFile(Form("ntp/%s_smallTree.root",prefix), "RECREATE");
  outFile->cd();
  outTree = new TTree("t","Tree");

//   //Set branch addresses
//   //variables must be declared in ossusy_looper_gmsb.h
//   outTree->Branch("weight",      &weight_,     "weight/F");
//   outTree->Branch("proc",        &proc_,       "proc/I");
//   outTree->Branch("leptype",     &leptype_,    "leptype/I");
//   outTree->Branch("dilmass",     &dilmass_,    "dilmass/F");
//   outTree->Branch("tcmet",       &tcmet_,      "tcmet/F");
//   outTree->Branch("tcmetphi",    &tcmetphi_,   "tcmetphi/F");
//   outTree->Branch("tcsumet",     &tcsumet_,    "tcsumet/F");
//   outTree->Branch("mt2",         &mt2_,        "mt2/F");  
//   outTree->Branch("mt2j",        &mt2j_,       "mt2j/F");  
//   outTree->Branch("sumjetpt",    &sumjetpt_,   "sumjetpt/F");
//   outTree->Branch("dileta",      &dileta_,     "dileta/F");
//   outTree->Branch("dilpt",       &dilpt_,      "dilpt/F");
//   outTree->Branch("dildphi",     &dildphi_,    "dildphi/F");
//   outTree->Branch("njets",       &njets_,      "njets/I");
//   outTree->Branch("vecjetpt",    &vecjetpt_,   "vecjetpt/F");
//   outTree->Branch("pass",        &pass_,       "pass/I");
//   outTree->Branch("passz",       &passz_,      "passz/I");
//   outTree->Branch("m0",          &m0_,         "m0/F");
//   outTree->Branch("m12",         &m12_,        "m12/F");
//   outTree->Branch("ptl1",        &ptl1_,       "ptl1/F");
//   outTree->Branch("ptl2",        &ptl2_,       "ptl2/F");
//   outTree->Branch("ptj1",        &ptj1_,       "ptj1/F");
//   outTree->Branch("ptj2",        &ptj2_,       "ptj2/F");
//   outTree->Branch("meff",        &meff_,       "meff/F");
//   outTree->Branch("mt",          &mt_,         "mt/F");


}

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

void ossusy_looper_gmsb::closeTree()
{
  outFile->cd();
  outTree->Write();
  outFile->Close();
  delete outFile;
}

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

ossusy_looper_gmsb::ossusy_looper_gmsb()
{
  g_susybaseline = false;
  g_createTree   = false;
  g_useBitMask   = false;
  random3_ = new TRandom3(1);
}

int ossusy_looper_gmsb::ScanChain(TChain* chain, char *prefix, float kFactor, int prescale,
                             JetTypeEnum jetType, MetTypeEnum metType, ZVetoEnum zveto, bool doFakeApp, bool calculateTCMET)
{
  
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  DorkyEventIdentifier dei;

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

  if(g_createTree) makeTree(prefix);

  while((currentFile = (TChainElement*)fileIter.Next())) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);

    unsigned int nEntries = tree->GetEntries();

    for(unsigned int z = 0; z < nEntries; ++z) {
      ++nEventsTotal;
      if (!((nEventsTotal+1)%10000))
        cout << "Processing event " << nEventsTotal+1 << " of " << nEventsChain << " in " << prefix << endl;

      cms2.GetEntry(z);

      // duplicate event?
      bool myduplicate = dei.is_duplicate(DorkyEvent());
      if( myduplicate ) continue;

      bool hlt_ele15_lw_l1r = passHLTTrigger("HLT_Ele15_LW_L1R");
      bool hltMu9           = passHLTTrigger("HLT_Mu9");
      
      if( !hltMu9 && !hlt_ele15_lw_l1r ) continue;

      VofP4 leptons_p4;
      vector<int> leptons_pdgid;
      vector<int> leptons_idx;

      leptons_p4.clear();
      leptons_pdgid.clear();
      leptons_idx.clear();

      int nLep = 0;
      int nEl  = 0;
      int nMu  = 0;

      //add electrons and muons to leptons_p4
      for( unsigned int iel = 0 ; iel < els_p4().size() ; ++iel ){
        
        if( els_p4().at(iel).pt() < 10 )                               continue;
        if( !pass_electronSelection( iel , electronSelection_ttbarV1 ) ) continue;
      
        leptons_p4.push_back( els_p4().at(iel) );
        leptons_pdgid.push_back( els_charge().at(iel) * -11 );
        leptons_idx.push_back( iel );
      
        ++nLep;
        ++nEl;
      }  

      for( unsigned int imu = 0 ; imu < mus_p4().size() ; ++imu ){
        
        if( mus_p4().at(imu).pt() < 10 )      continue;
        if( !muonId( imu , NominalTTbar ) )   continue;
        
        leptons_p4.push_back( mus_p4().at(imu) );
        leptons_pdgid.push_back( mus_charge().at(imu) * -13 );
        leptons_idx.push_back( imu );
  
        ++nLep;
        ++nMu;
      } 

      VofP4 mypfjets_p4;
      mypfjets_p4.clear();

      float sumJetPt = 0;
      int nJets = 0;

      for (unsigned int ijet = 0; ijet < pfjets_p4().size(); ijet++) {

        LorentzVector vjet = pfjets_p4().at(ijet);
        if ( vjet.pt()  < 30   ) continue;
        if ( vjet.eta() > 2.4  ) continue;

        bool skipJet = false;

        for( unsigned int ilep = 0 ; ilep < leptons_p4.size() ; ++ilep){
          LorentzVector vlep  = leptons_p4.at(ilep);
          if (dRbetweenVectors(vjet, vlep) < 0.4) skipJet = true;
        }
        if( skipJet ) continue;
        
        mypfjets_p4.push_back( vjet );
        ++nJets;
        sumJetPt += vjet.pt();
        
      }
      
      float pfmet = evt_pfmet();

      float pt1 = -1.;
      float pt2 = -1.;
      float pt3 = -1.;
      float pt4 = -1.;

      int ilep1 = -1;
      int ilep2 = -1;
      int ilep3 = -1;
      int ilep4 = -1;

      if( nLep > 0 ){
        for( unsigned int ilep = 0 ; ilep < leptons_p4.size() ; ++ilep ){
          float pt = leptons_p4.at(ilep).pt();
          if( pt > pt1 ){
            pt1   = pt;
            ilep1 = ilep;
          }        
        }
      }
    
      if( nLep > 1 ){
        for( unsigned int ilep = 0 ; ilep < leptons_p4.size() ; ++ilep ){
          if( ilep == ilep1 ) continue;
          float pt = leptons_p4.at(ilep).pt();
          if( pt > pt2 ){
            pt2   = pt;
            ilep2 = ilep;
          }        
        }
      }
     
      if( nLep > 2 ){
        for( unsigned int ilep = 0 ; ilep < leptons_p4.size() ; ++ilep ){
          if( ilep == ilep1 || ilep == ilep2 ) continue;
          float pt = leptons_p4.at(ilep).pt();
          if( pt > pt3 ){
            pt3   = pt;
            ilep3 = ilep;
          }        
        }
      }

      if( nLep > 3 ){
        for( unsigned int ilep = 0 ; ilep < leptons_p4.size() ; ++ilep ){
          if( ilep == ilep1 || ilep == ilep2 || ilep == ilep3 ) continue;
          float pt = leptons_p4.at(ilep).pt();
          if( pt > pt4 ){
            pt4   = pt;
            ilep4 = ilep;
          }        
        }
      }
      
      float weight = evt_scale1fb() * 0.1;

      fillUnderOverFlow( hnleptons, nLep ,weight );
      fillUnderOverFlow( hnelectrons, nEl ,weight );
      fillUnderOverFlow( hnmuons, nMu ,weight );
      fillUnderOverFlow( hmet, pfmet ,weight );
      fillUnderOverFlow( hnjets, nJets ,weight );
      fillUnderOverFlow( hsumjetpt, sumJetPt ,weight );

      if( nLep > 0 ) fillUnderOverFlow( hpt1 , pt1 , weight );
      if( nLep > 1 ) fillUnderOverFlow( hpt2 , pt2 , weight );
      if( nLep > 2 ) fillUnderOverFlow( hpt3 , pt3 , weight );
      if( nLep > 3 ) fillUnderOverFlow( hpt4 , pt4 , weight );
      
      for( unsigned int ilep = 0 ; ilep < leptons_p4.size() ; ++ilep ){

        fillUnderOverFlow( hpt , leptons_p4.at(ilep).pt() , weight );
        LorentzVector vi = leptons_p4.at( ilep );

        for( unsigned int jlep = ilep+1 ; jlep < leptons_p4.size() ; ++jlep ){

          LorentzVector vj = leptons_p4.at( jlep );
          float dR = dRbetweenVectors( vi , vj);
          fillUnderOverFlow( hdr_all , dR , weight );
          
          float dilmass2 = ( vi + vj ).M2();
          float dilmass = dilmass2 > 0 ? sqrt(dilmass2) : sqrt(-dilmass2);
          fillUnderOverFlow(hdilmass_all ,  dilmass , weight );
        }
      }

      if ( nLep > 1 ){
        LorentzVector v1 = leptons_p4.at( ilep1 );
        LorentzVector v2 = leptons_p4.at( ilep2 );
        float dR = dRbetweenVectors( v1 , v2);
        fillUnderOverFlow(hdr_ptmax ,  dR , weight );

        float dilmass2 = ( v1 + v2 ).M2();
        float dilmass = dilmass2 > 0 ? sqrt(dilmass2) : sqrt(-dilmass2);
        fillUnderOverFlow(hdilmass_ptmax ,  dilmass , weight );

        hsign->Fill(  leptons_pdgid.at(ilep1) * leptons_pdgid.at(ilep2) > 0 ? 1.5 : 0.5 , weight );
      }
      
      if( sumJetPt < 200. )  continue;
      if( pfmet    < 80.  )  continue;
      if( nLep     < 2    )  continue;
      if( nJets    < 3    )  continue;

      bool samesign = false;
      if( nLep > 2 ) samesign = true;
      else if( nLep > 1 ){
        if( leptons_pdgid.at(ilep1) * leptons_pdgid.at(ilep2) > 0 ) samesign = true;
      }
      if( !samesign ) continue;
      
      fillUnderOverFlow( hnleptons_all_pass , nLep , weight );
      
      if( nLep > 1 ){
        if( leptons_pdgid.at(ilep1) * leptons_pdgid.at(ilep2) > 0 )
          fillUnderOverFlow( hnleptons_ss_pass , nLep , weight );
        
        if( leptons_pdgid.at(ilep1) * leptons_pdgid.at(ilep2) < 0 )
          fillUnderOverFlow( hnleptons_os_pass , nLep , weight );
      }
      

    } // entries
  } // currentFile
      
  if(g_createTree) closeTree();
      
  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    
  return 0;
}

 
void ossusy_looper_gmsb::BookHistos(char *prefix)
{
  // Prefix comes from the sample and it is passed to the scanning function
  // Suffix is "ee" "em" "em" "all" which depends on the final state
  // For example: histogram named tt_hnJet_ee would be the Njet distribution
  // for the ee final state in the ttbar sample.
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  cout << "Begin book histos..." << endl;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  hnleptons_all_pass   = new TH1F(Form("%s_hnleptons_all_pass",prefix), "Num Leptons Passing Cuts + Event Selection",      11,-0.5,10.5);
  hnleptons_ss_pass    = new TH1F(Form("%s_hnleptons_ss_pass",prefix),  "Num SS Leptons Passing Cuts + Event Selection",   11,-0.5,10.5);
  hnleptons_os_pass    = new TH1F(Form("%s_hnleptons_os_pass",prefix),  "Num OS Leptons Passing Cuts + Event Selection",   11,-0.5,10.5);

  hsign           = new TH1F(Form("%s_hsign",prefix),"OS(0) SS (1)",                      2,0,2);
  hnleptons       = new TH1F(Form("%s_hnleptons",prefix),"Num Leptons Passing Cuts",      11,-0.5,10.5);
  hnelectrons     = new TH1F(Form("%s_hnelectrons",prefix),"Num Electrons Passing Cuts",  11,-0.5,10.5);
  hnmuons         = new TH1F(Form("%s_hnmuons",prefix),"Num Muons Passing Cuts",          11,-0.5,10.5);
  hmet            = new TH1F(Form("%s_hmet",prefix),"pfmet",                              100,0,500);
  hnjets          = new TH1F(Form("%s_hnjets",prefix),"nJets",                            11,-0.5,10.5);
  hsumjetpt       = new TH1F(Form("%s_hsumjetpt",prefix),"sumJetPt",                      100,0,1000);
  hpt             = new TH1F(Form("%s_hpt",prefix),"p_{T} all leptons",                   100,0,500);
  hpt1            = new TH1F(Form("%s_hpt1",prefix),"p_{T} of 1^{st} lepton",             100,0,500);
  hpt2            = new TH1F(Form("%s_hpt2",prefix),"p_{T} of 2^{nd} lepton",             100,0,500);
  hpt3            = new TH1F(Form("%s_hpt3",prefix),"p_{T} of 3^{rd} lepton",             100,0,500);
  hpt4            = new TH1F(Form("%s_hpt4",prefix),"p_{T} of 4^{th} lepton",             100,0,500);
  hdr_all         = new TH1F(Form("%s_hdr_all",prefix),"#DeltaR between all selected leptons",               100,0,5);
  hdilmass_all    = new TH1F(Form("%s_hdilmass_all",prefix),"M_{ll} for all selected leptons",               100,0,500);
  hdr_ptmax       = new TH1F(Form("%s_hdr_ptmax",prefix),"#DeltaR between 2 highest p_{T} selected leptons", 100,0,5);
  hdilmass_ptmax  = new TH1F(Form("%s_hdilmass_ptmax",prefix),"M_{ll} for 2 highest p_{T} selected leptons", 100,0,500);

  cout << "End book histos..." << endl;
}// CMS2::BookHistos()


void fillUnderOverFlow(TH1F *h1, float value, float weight)
{
  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();
  
  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}
