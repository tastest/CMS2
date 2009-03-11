#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <set>
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "Math/GenVector/LorentzVector.h"

using namespace std;

#ifndef __CINT__
#include "CMS2_V00_04_00.h"
CMS2 cms2;
#endif

#include "../Tools/selections.C"
#include "../Tools/utilities.C"
#include "warren_functions.C"

//wandrews: This was originally Oli's looper. I've changed it to be a dilepton Generator level looper for use in my susy stuff.

//idea behind buckets: buckets are for final state (e/mu +/-) but also for other particles in event like LSP,nu,b,t,etc (see below)
//we use buckets for all variables: if there is a variable worth looking at for one particle type, it may as well be plotted for all types.

//all other function definitions are in "warren_functions.C"
char* name_from_cuts(char* prefix, char* tag);
const char* print_cuts();

//global vars for CUTS so i don't have to pass them all over the place
bool EtaCut = true;
double EtaCutValue = 2.5;
bool TausAreLeptons = true;
bool PtCut = false;  
double PtCutValue = 50;  //for jets, not necessary for leptons

//other global vars
//LSP's mcid
const Int_t LSP_mcid = 1000022; //chi_10

//define bucket vars--indicies in hist arrays
const unsigned int lepBuckets = 12;//10e,u, 2tau
const unsigned int othBuckets = 7; //tau, lsp, neutrino, zjet, b, t, zlep
const unsigned int allBuckets = lepBuckets + othBuckets + 1;
const unsigned int sst = 10; //same-sign tau
const unsigned int ost = 11; //opposite--
const unsigned int tau = 12;
const unsigned int nu = 13;
const unsigned int LSP = 14;
const unsigned int zjet = 15;
const unsigned int zlep = 16;//ZERO LEPTON ARE NOT IN 'ALL' BUCKET
const unsigned int b = 17;
const unsigned int t = 18;

char *suffix[allBuckets];
//this function's only purpose in life is to declare suffix
void init_suffix() {
  //di-lepton:
  suffix[0]  = "mpmp";
  suffix[1]  = "mpmm";
  suffix[2]  = "mpep";
  suffix[3]  = "mpem";
  suffix[4]  = "mmmm";
  suffix[5]  = "mmep";
  suffix[6]  = "mmem";
  suffix[7] = "epep";
  suffix[8] = "epem";
  suffix[9] = "emem";
  //new:
  suffix[sst] = "sst";
  suffix[ost] = "ost";
  suffix[tau] = "tau";
  suffix[nu] = "nu";
  suffix[LSP] = "LSP";
  suffix[zjet] = "zero-jet";
  suffix[zlep] = "zero-lep"; //do not include zero-lep in 'all' bucket
  suffix[b] = "b";
  suffix[t] = "t"; // top bucket
  suffix[allBuckets - 1] = "all";
}

char* ScanChain( TChain* chain, char * prefix="", char* tag="", int specDY=-1, float kFactor=1.0) {

  char *name = name_from_cuts(prefix, tag);

  // Make sure the specDY flags is kosher
  if (specDY < -1 || specDY > 2) {
    std::cout << "specDY flag is not allowed...quit" << std::endl;
    return "";
  }

  // clear list of duplicates
  already_seen.clear();
  int duplicates_total_n = 0;
  double duplicates_total_weight = 0;
  int i_permille_old = 0;

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;
  unsigned int nCandidatesSelected = 0;
  int filenum = 0;
  float weight=0;

  //particle count variables: taus
  double numtotaltau = 0; //taus which aren't in dilep selection
  double numdileptau = 0; //taus which are in dilep events selected by e,mu
  double numevttau = 0;
  double numevtdileptau = 0;
  //b,t,glu -- vectors for b,t to count number of events with 1,2,3 b,t
  const unsigned int npartcat = 6; //num categories--see count_particles func
  double* numdilepb = new double[npartcat];//was ntotalb-now dilep and tau 
  double* numtaub = new double[npartcat];
  double* numdilept = new double[npartcat];
  double* numtaut = new double[npartcat];
  double* numdilepz = new double[npartcat];
  double* numtauz = new double[npartcat];
  double* numdileph = new double[npartcat];
  double* numtauh = new double[npartcat];
  double numtotalgluon = 0;
  double numevtgluon = 0;    
  //int nbaddecay = 0; //was for gluino->gaugino+glu, which is allowed at 1loop
  int nHardlsp = 0;
  double in_quark_typ[5][2];//up,d-bar,d,u-bar,glu for counting charge
  for(unsigned int i=0; i<npartcat; i++){
	numdilepb[i] = 0;
	numtaub[i] = 0;
	numdilept[i] = 0;
	numtaut[i] = 0;
	numdilepz[i] = 0;
	numtauz[i] = 0;
	numdileph[i] = 0;
	numtauh[i] = 0;
	if( i < 5 ) {
	  for(int j=0;j<2;j++){ in_quark_typ[i][j]=0; } //pp vs mm 
	}
  }
  
  //chgflvs: sssf,ssof,ossf,osof,sst,ost--last two for taus
  unsigned int numchgflv;
  if( TausAreLeptons )
	numchgflv = 6;
  else
	numchgflv = 4;
  //just number of pair produced types, see warren_functions.C
  const unsigned int susy_types = 10;//gg,qg,qq,qX,qs,XX,sX,ss,gX,hh
  double* hard_type = new double[susy_types];
  double* zjet_hard_type = new double[susy_types];
  //this is for charge/flavor of susy pair produced
  //first index is sign/flavor (same/opp), second susy type
  double hard_chgflv[numchgflv][susy_types];
  double zjet_hard_chgflv[numchgflv][susy_types];
  for(unsigned int i=0; i<susy_types;i++) {
	hard_type[i] = 0; zjet_hard_type[i] = 0;
	for(unsigned int j=0;j<numchgflv;j++) {
	  hard_chgflv[j][i] = 0;
	  zjet_hard_chgflv[j][i] = 0;
	}
  }
  init_suffix();

  ////Lepton hists
  TH1F* hPtLepFirst[allBuckets]; //pt max
  TH1F* hPtLepSecond[allBuckets];
  TH1F* hPtLepThird[allBuckets];
  //scalar sum of pt of all leptons in the event
  TH1F* hPtLep[allBuckets];
  //dilepton invariant mass--for now just of two highest pt leptons
  TH1F* hMassLep[allBuckets];
  //num leptons!=tau,!=nu for lepton buckets. For nu bucket, num(nu)
  TH1F* hNLep[allBuckets];

  ////Hadron hists
  //scalar sum of pt of all g,q!=t in the event
  TH1F* hSumEt[allBuckets];
  TH1F* hPtHadFirst[allBuckets]; //pt max of g,q!=t
  TH1F* hPtHadSecond[allBuckets];
  TH1F* hPtHadThird[allBuckets];
  //TH1F* hMassHad[allBuckets]; //mass of leading q, and b's separately
  TH1F* hNjets[allBuckets]; //num g,q!=t
  TH1F* halpha[allBuckets]; //alpha = Etj2/Mjj
  //TH1F* halphat[allBuckets]; //alphat = Etj2/Mtjj
  // leave commented b'c redundant with alpha--very little change noticed, and not obvious that this change is in the 'right' direction.

  //Et, MET hists
  //MET = magnitude Et vector:
  //  LSP bucket is filled with vector sum of Et of LSPs
  TH1F* hEtFirst[allBuckets];
  TH1F* hEtSecond[allBuckets];
  TH1F* hPhi[allBuckets];

  //magnitude of vector sum of Et of LSP+nu
  TH1F* hMET[allBuckets];
  //MEff = mag(met) + SumEt
  TH1F* hMEff[allBuckets];
  //Ht = SumEt + mag(met) + ptLep
  TH1F* hHt[allBuckets];
  
  //chrgflv plots: zjets filled with non-chgflv
  TH1F* hchgflv_PtLepFirst_zjet[numchgflv];
  TH1F* hchgflv_PtLepSecond_zjet[numchgflv];
  TH1F* hchgflv_PtLepThird_zjet[numchgflv];
  TH1F* hchgflv_PtLep_zjet[numchgflv];
  TH1F* hchgflv_MassLep_zjet[numchgflv];
  TH1F* hchgflv_MET_zjet[numchgflv];
  TH1F* hchgflv_SumEt_zjet[numchgflv];
  TH1F* hchgflv_MEff_zjet[numchgflv];
  TH1F* hchgflv_Ht_zjet[numchgflv];
  TH1F* hchgflv_Njets_zjet[numchgflv];
  TH1F* hchgflv_alpha_zjet[numchgflv];

  //b:
  TH1F* hchgflv_PtHadFirst_b[numchgflv];
  TH1F* hchgflv_PtHadSecond_b[numchgflv];
  TH1F* hchgflv_PtHadThird_b[numchgflv];
  TH1F* hchgflv_SumEt_b[numchgflv];
  TH1F* hchgflv_Njets_b[numchgflv];

  for (unsigned int i=0; i < allBuckets; ++i) {
    //Lepton Pt
    const float lowBin = 0.;

    const float highBinF = 400;
    const int nBinsF = 80;

	const float highBinS = 200;
    const int nBinsS = 100;

	hPtLep[i]     = book1DHist(Form("%s_hPtLep_%s",prefix,suffix[i]),
							   Form("%s_hPtLep_%s",prefix,suffix[i]),
							   nBinsF,lowBin,highBinF,"p_{T}lep [GeV]","");
    hPtLepFirst[i]= book1DHist(Form("%s_hPtLepFirst_%s",prefix,suffix[i]),
							   Form("%s_hPtLepFirst_%s",prefix,suffix[i]),
							   nBinsS,lowBin,highBinS,"p_{T} [GeV]","");   
    hPtLepSecond[i]=book1DHist(Form("%s_hPtLepSecond_%s",prefix,suffix[i]),
							   Form("%s_hPtLepSecond_%s",prefix,suffix[i]),
							   nBinsS,lowBin,highBinS,"p_{T} [GeV]","");   
    hPtLepThird[i]= book1DHist(Form("%s_hPtLepThird_%s",prefix,suffix[i]),
							   Form("%s_hPtLepThird_%s",prefix,suffix[i]),
							   nBinsS,lowBin,highBinS,"p_{T} [GeV]","");   
	hMassLep[i]   = book1DHist(Form("%s_hMassLep_%s",prefix,suffix[i]),
							   Form("%s_hMassLep_%s",prefix,suffix[i]),
							   nBinsF,lowBin,highBinF,"Di-lep Mass [GeV]","");
	
	//Hadron Pt
    const float highBinG = 500;
    const int nBinsG = 100;
    const float highBinH = 1000;
    const int nBinsH = 100;
    const float highBinB = 1200;
    const int nBinsB = 120;

	//new notation: sumet is just hadron (jet) pt
	hSumEt[i]      = book1DHist(Form("%s_hSumEt_%s",prefix,suffix[i]),
								Form("%s_hSumEt_%s",prefix,suffix[i]),
								nBinsB,lowBin,highBinB,"SumEt (quark p_{T}) [GeV]","");
    hPtHadFirst[i] = book1DHist(Form("%s_hPtHadFirst_%s",prefix,suffix[i]),
								Form("%s_hPtHadFirst_%s",prefix,suffix[i]),
								nBinsG,lowBin,highBinG,"p_{T} [GeV]","");   
    hPtHadSecond[i]= book1DHist(Form("%s_hPtHadSecond_%s",prefix,suffix[i]),
								Form("%s_hPtHadSecond_%s",prefix,suffix[i]),
								nBinsF,lowBin,highBinF,"p_{T} [GeV]","");   
    hPtHadThird[i] = book1DHist(Form("%s_hPtHadThird_%s",prefix,suffix[i]),
								Form("%s_hPtHadThird_%s",prefix,suffix[i]),
								nBinsS,lowBin,highBinS,"p_{T} [GeV]","");   
	halpha[i]      = book1DHist(Form("%s_hAlpha_%s",prefix,suffix[i]), 
								Form("%s_hAlpha_%s",prefix,suffix[i]),
								20,0,1,"alpha","");
	//halphat[i]     = book1DHist(Form("%s_hAlphat_%s",prefix,suffix[i]), 
	//							Form("%s_hAlphat_%s",prefix,suffix[i]),
	//							20,0,1,"alpha-t","");
	
    // Et, MET 
    const float highBinM = 600;
    const int nBinsM = 120;
    const float highBinT = 2000;
    const int nBinsT = 100;
	
    hEtFirst[i] = book1DHist(Form("%s_hEtFirst_%s",prefix,suffix[i]),
							 Form("%s_hEtFirst_%s",prefix,suffix[i]),
							 nBinsM,lowBin,highBinM,"Et 1st [GeV]","");
    hEtSecond[i]= book1DHist(Form("%s_hEtSecond_%s",prefix,suffix[i]),
							 Form("%s_hEtSecond_%s",prefix,suffix[i]),
							 nBinsG,lowBin,highBinG,"Et 2nd [GeV]","");
    hPhi[i]     = book1DHist(Form("%s_hPhi_%s",prefix,suffix[i]),
							 Form("%s_hPhi_%s",prefix,suffix[i]),
							 10,-3.14,3.14,"Et phi [GeV]","");

	//MET = magnitude of Et(lsp + nu, vector sum)
    hMET[i]     = book1DHist(Form("%s_hMET_%s",prefix,suffix[i]),
							 Form("%s_hMET_%s",prefix,suffix[i]),
							 nBinsH,lowBin,highBinH,"MET [GeV]","");
	//MEff = mag(MET) + SumEt
    hMEff[i]    = book1DHist(Form("%s_hMEff_%s",prefix,suffix[i]),
							 Form("%s_hMEff_%s",prefix,suffix[i]),
							 nBinsT,lowBin,highBinT,"MEff [GeV]","");
	//Ht = mag(MET) + PtLep + SumEt (scalar sum)
    hHt[i]      = book1DHist(Form("%s_hHt_%s",prefix,suffix[i]),
							 Form("%s_hHt_%s",prefix,suffix[i]),
							 nBinsT,lowBin,highBinT,"Ht [GeV]","");

    // nJets, nLep
    const float highBinNjets = 15.;
    const int nBinsNjets = 15;
    const float highBinNlep = 10.;
    const int nBinsNlep = 10;

	hNjets[i] = book1DHist(Form("%s_hNjets_%s",prefix,suffix[i]),
						   Form("%s_hNjets_%s",prefix,suffix[i]),
						   nBinsNjets,lowBin,highBinNjets,"Njets","");
	hNLep[i] = book1DHist(Form("%s_hNLep_%s",prefix,suffix[i]),
						  Form("%s_hNLep_%s",prefix,suffix[i]),
						  nBinsNlep,lowBin,highBinNlep,"Nleptons","");
  }
  //end for i < allBuckets
  
  for(unsigned int i=0;i<numchgflv;i++ ) {
	//note: for zjet, no sumEt plot, no MEff (b'c = Ht), and no Njets
	hchgflv_PtLepFirst_zjet[i] = bookChgflv(hPtLepFirst[zjet], i);
	hchgflv_PtLepSecond_zjet[i] = bookChgflv(hPtLepSecond[zjet], i);
	hchgflv_PtLepThird_zjet[i] = bookChgflv(hPtLepThird[zjet], i);
	hchgflv_PtLep_zjet[i] = bookChgflv(hPtLep[zjet], i);
	hchgflv_MassLep_zjet[i] = bookChgflv(hMassLep[zjet], i);
	hchgflv_MET_zjet[i] = bookChgflv(hMET[zjet], i);
	hchgflv_SumEt_zjet[i] = bookChgflv(hSumEt[zjet], i);
	hchgflv_MEff_zjet[i] = bookChgflv(hMEff[zjet], i);
	hchgflv_Ht_zjet[i] = bookChgflv(hHt[zjet], i);
	hchgflv_Njets_zjet[i] = bookChgflv(hNjets[zjet], i);
	hchgflv_alpha_zjet[i] = bookChgflv(halpha[zjet], i);
	//b:
	hchgflv_PtHadFirst_b[i] = bookChgflv(hPtHadFirst[b], i);
	hchgflv_PtHadSecond_b[i] = bookChgflv(hPtHadSecond[b], i);
	hchgflv_PtHadThird_b[i] = bookChgflv(hPtHadThird[b], i);
	hchgflv_SumEt_b[i] = bookChgflv(hSumEt[b], i);
	hchgflv_Njets_b[i] = bookChgflv(hNjets[b], i);
  }


  //cout << "\nBegin File Loop\n";
  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
	filenum++;

    //cout << "\nBeginning Event Loop\n";
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
	//unsigned int nEvents = 20;
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;

      DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.trks_d0()[0],
								  cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
      if (is_duplicate(id)) {
		duplicates_total_n++;
		duplicates_total_weight += cms2.evt_scale1fb();
		continue;
      }

      // Progress feedback to the user
      //progressBar(i_permille_old, nEventsTotal, nEventsChain);

      // The event weight including the kFactor (scaled to 1 fb-1)
      weight = cms2.evt_scale1fb() * kFactor;
	  //if( cms2.evt_scale1fb() * kFactor != weight )
	  //cout << "\n\nWARNING: Weight changed from original value\n\n";

      // special handling for DY (Check if ntuples run on needs it)
      bool processEvent=true;
      if (specDY == 0) {
		if ( !isDYee() ) processEvent = false;
      } else if (specDY == 1) {
		if ( !isDYmm() ) processEvent = false;
      } else if (specDY == 2) {
		if ( !isDYtt() ) processEvent = false;
      }
      if (!processEvent) { cout<<"\n\nDYERROR\n\n"; continue; }

	  vector<Int_t> list_id ;
	  //for(int i=0;i<100;i++) { list_id[i]=0; } //array implementation

	  double* id_lep = new double[3];
	  double* idx_lep = new double[3]; //for 4-vectors
	  double* id_tau = new double[3];
	  double* idx_tau = new double[3]; //for 4-vectors
	  double* idx_had = new double[3]; //for 4-vectors
	  double* pt_rank_lep = new double[3];
	  double* pt_rank_tau = new double[3];
	  double* pt_rank_had = new double[3]; 
	  double* pt_rank_b = new double[3]; 
      double* metlsp = new double[2]; //guaranteed to have 2 lsp in event
	  double* metnu = new double[3]; //forget about the rest (if > 3)
	  //vector<double> idx_alllep;
	  for(int i=0;i<3;i++) {
		id_lep[i]=0;
		idx_lep[i]=0;
		id_tau[i]=0;
		idx_tau[i]=0;
		idx_had[i]=0;
		pt_rank_lep[i]=0;
		pt_rank_tau[i]=0;
		pt_rank_had[i]=0;
		pt_rank_b[i]=0;
		if(i<2) metlsp[i]=0;
		metnu[i]=0;
	  }

      unsigned int bucket = 0;
      int numlep = 0; //NOT including taus
	  int numnu = 0; int numtau = 0;	  
	  int numhad = 0;//numhad includes b, not t--use for njets
	  int numhadeta = 0; //includes outside eta range--use for zjet
	  int numb = 0; int numbfromt = 0; 
	  int numh = 0; int numteta = 0; int numt = 0;
	  int numz = 0; int numlepfromz = 0; int numtaufromz = 0;
	  int numgluon = 0; int numnufromz = 0;
      double pthad = 0, ptlep = 0, pttau = 0, ptb = 0; 
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vecEt;
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vecEtlsp;
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vecEtnu;
	  vecEt.SetPx(0); vecEt.SetPy(0); vecEt.SetPz(0); vecEt.SetE(0);
	  vecEtlsp.SetPx(0); vecEtlsp.SetPy(0);
	  vecEtlsp.SetPz(0); vecEtlsp.SetE(0);
	  vecEtnu.SetPx(0); vecEtnu.SetPy(0); vecEtnu.SetPz(0); vecEtnu.SetE(0);
	  //lepton 4vectors
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lep1;
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lep2;
	  lep1.SetPx(0); lep1.SetPy(0); lep1.SetPz(0); lep1.SetE(0);
	  lep2.SetPx(0); lep2.SetPy(0); lep2.SetPz(0); lep2.SetE(0);
	  //hadron 4vectors
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > had1;
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > had2;
	  had1.SetPx(0); had1.SetPy(0); had1.SetPz(0); had1.SetE(0);
	  had2.SetPx(0); had2.SetPy(0); had2.SetPz(0); had2.SetE(0);
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > b1;
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > b2;
	  b1.SetPx(0); b1.SetPy(0); b1.SetPz(0); b1.SetE(0);
	  b2.SetPx(0); b2.SetPy(0); b2.SetPz(0); b2.SetE(0);

	  for( int i=0;i<6;i++) {
		//list_id[i]=cms2.genps_id()[i];
		list_id.push_back( cms2.genps_id()[i] );
	  }
	  //loop starts at 6 b'c 0,1 are protons, 2-5 are initial state in proton
      for ( unsigned int par = 6; 
			par < cms2.genps_id().size(); ++par ) {

		double id = cms2.genps_id()[par];
		//cms2.genps_idmother()[par];
		list_id.push_back( (Int_t)id );
		
		//quark!=t, t=6
		if( (abs(id) >= 1 && abs(id) <= 5) || abs(id) == 21  ) { 
		  numhadeta++;//any eta range, any pt
		  if( abs(cms2.genps_id_mother()[par]) == 6 && abs(id) == 5 )
			  numbfromt++;
		  if( EtaCut && abs(cms2.genps_p4()[par].Eta()) > EtaCutValue )
			continue;
		  if( PtCut && abs(cms2.genps_p4()[par].pt()) < PtCutValue )
			continue;

		  //scalar sum of g,q!=t pt
		  pthad += cms2.genps_p4()[par].pt();
		  numhad++;//only eta < EtaCutValue
		  if( abs(id) == 5 ) {
			numb++;
			ptb += cms2.genps_p4()[par].pt();
			pt_rank_b = get_pt_rank(pt_rank_b, cms2.genps_p4()[par].pt() );
			//add idx_b
		  }
		  else if( abs(id) == 21 ) {
			numgluon++;//only gluons found are from h0->gg, h0 from X20->h0X10
		  }
		  //see below (comment above idx_lep = get_id_rank(...) )
		  idx_had = get_id_rank(par, pt_rank_had, idx_had,
								cms2.genps_p4()[par].pt() );
		  pt_rank_had = get_pt_rank(pt_rank_had, cms2.genps_p4()[par].pt() );
		  
		}
		else if( abs(id) == 6 ) { //t quark
		  numteta++;
		  if( EtaCut && abs(cms2.genps_p4()[par].Eta()) > EtaCutValue )
			continue;
		  if( PtCut && abs(cms2.genps_p4()[par].pt()) < PtCutValue )
			continue;
		  numt++;
		}
		else if( abs(id) == 11 || abs(id) == 13 ) { //e,mu
		  if( EtaCut && abs(cms2.genps_p4()[par].Eta()) > EtaCutValue )
			continue;
		
		  ptlep += cms2.genps_p4()[par].pt() ;
		  numlep++;

		  //NOTE: these two lines MUST be before the get_pt_rank call
		  //because pt_rank_lep should be OLD for these two lines
		  id_lep = get_id_rank(id, pt_rank_lep, id_lep,
							   cms2.genps_p4()[par].pt() );
		  //Get index in genps block of 3 highest pt leptons
		  //use id rank, only send it the index in the genps block,
		  // do NOT give it the mcid. 
		  idx_lep = get_id_rank(par, pt_rank_lep, idx_lep,
								cms2.genps_p4()[par].pt() );

		  pt_rank_lep = get_pt_rank(pt_rank_lep, cms2.genps_p4()[par].pt() );
		  //leptau useless b'c need to use e,mu separately
		  //This one always filled b'c it teats tausarelep=true always
		  //pt_rank_leptau = get_pt_rank(pt_rank_leptau, cms2.genps_p4()[par].pt() );
		  //idx_alllep.push_back( par );
		  if( abs(cms2.genps_id_mother()[par]) == 23 )
			numlepfromz++;
		}
		else if( abs(id) == 15 ) {  //tau
		  if( EtaCut && abs(cms2.genps_p4()[par].Eta()) > EtaCutValue )
			continue;

		  pttau += cms2.genps_p4()[par].pt();
		  numtau++;

		  id_tau = get_id_rank(id, pt_rank_tau, id_tau,
							   cms2.genps_p4()[par].pt() );
		  //will use this to merge lep+tau after loop
		  idx_tau = get_id_rank(par, pt_rank_tau, idx_tau,
								cms2.genps_p4()[par].pt() );

		  pt_rank_tau = get_pt_rank(pt_rank_tau, cms2.genps_p4()[par].pt() );
		  //pt_rank_leptau = get_pt_rank(pt_rank_leptau, cms2.genps_p4()[par].pt());
		  //if( TausAreLeptons ) 
		  //idx_alllep.push_back( par );

		  if( abs(cms2.genps_id_mother()[par]) == 23 )
			numtaufromz++;
		}
		else if( abs(id) == 12 || abs(id) == 14 || abs(id) == 16 ) {
		  numnu++;
		  metnu = get_pt_rank(metnu, cms2.genps_p4()[par].Et() );
		  vecEt += cms2.genps_p4()[par];
		  vecEtnu += cms2.genps_p4()[par];
		  if( abs(cms2.genps_id_mother()[par]) == 23 )
			numnufromz++;
		  //hPhi[nu]->Fill(cms2.genps_p4()[par].Phi() , weight);
		}
		else if( abs(id) == 23 ) {
		  numz++;
		}
		else if( abs(id) == 25 ) {
		  numh++;
		}
		else if( abs(id) == LSP_mcid ) {
		  metlsp = get_Et_rank(metlsp, cms2.genps_p4()[par].Et() );
		  vecEt += cms2.genps_p4()[par];
		  vecEtlsp += cms2.genps_p4()[par];
		  hPhi[LSP]->Fill(cms2.genps_p4()[par].Phi() , weight);
		}
		//in case I want to find tree w/ dave's cmssw filter
		//cout << "event: " << cms2.evt_event() << " run: " << cms2.evt_run() << endl;
      }
	  //end genps loop

	  //CUT: must have number leptons >= 2
      if(numlep < 2){
		bool die = false;
		if( numlep == 0 || (numtau+numlep) < 2 ) //lines 2,4
		  die = true;
		else if( !TausAreLeptons )
		  die = true;

		//count taus which could be used if 1 e,mu, 1 tau: category B
		if( numtau + numlep >= 2 && numlep == 1 ) {
		  numtotaltau += numtau*weight;
		  numevttau += weight; //events in category B
		}

		if( die ) {
		  //Fill z-lep bucket.
		  //NOTE: FOR NOW, ZLEP MEANS ZERO OR ONE LEPTON
		  hPtLep[zlep]->Fill(ptlep, weight );
		  hPtLepFirst[zlep]->Fill(pt_rank_lep[0], weight );
		  hPtLepSecond[zlep]->Fill(pt_rank_lep[1], weight );
		  hPtLepThird[zlep]->Fill(pt_rank_lep[2], weight );
		  hSumEt[zlep]->Fill( pthad, weight );
		  hPtHadFirst[zlep]->Fill(pt_rank_had[0], weight );
		  hPtHadSecond[zlep]->Fill(pt_rank_had[1], weight );
		  hPtHadThird[zlep]->Fill(pt_rank_had[2], weight );
		  hNjets[zlep]->Fill(numhad, weight );
		  hMET[zlep]->Fill(vecEt.Et(), weight );
		  hMEff[zlep]->Fill(vecEt.Et() + pthad);
		  hHt[zlep]->Fill(vecEt.Et() + ptlep + pthad, weight);
		  
		  continue;
		}
	  }

	  int chgflv = 999;
      //determine bucket for lepton vars--first 10 buckets if >= 2 e or mu
	  if( !TausAreLeptons || pt_rank_lep[1] != 0 )// pt_rank_tau[0] ) 
		bucket = get_lep_bucket(id_lep[0], id_lep[1]);
	  else {
		//returns either sst or ost
		bucket = get_lep_bucket(id_lep[0], id_tau[0], sst, ost);
		if( bucket == sst ) 
		  chgflv = 4;
		else if( bucket == ost )
		  chgflv = 5;
		else
		  chgflv = 999;
	  }

	  //chgflv is for sign,flavor combo of final state
	  if( bucket <= 9 )
		chgflv = buck_to_chgflv(bucket);
	  else if( chgflv > 5 || chgflv < 0  )
		cout << "\nBAD CHGFLV\n";
	  
	  ++nCandidatesSelected;      	  

	  //count particles of interest--decay products
	  //int nother = numtaufromz + numlepfromz - numz;
	  if( chgflv < 4 ) {
		numdileptau += numtau*weight; //category A
		if( numtau > 0 ) numevtdileptau += weight; //A
		numdilepb = count_particles( numdilepb, numb, 0, npartcat, weight );
		int nhadz =  numz - (numlepfromz + numnufromz)/2;//TAU ARE HADRONIC
		numdilepz = count_particles( numdilepz, numlepfromz/2, nhadz, npartcat, weight);
		numdileph = count_particles( numdileph, numh, 0, npartcat, weight );
		//if(  != numz) //moved into count_particles
		//numdilepz[3] += weight; //3rd indx is other--non-leptonic z

		//NOTE: t->W->tau IS COUNTED AS HADRONIC--see warren_functions.C
		int hd = num_hadronic_top( list_id );
		if( hd >= 0 ) {
		  numdilept = count_particles(numdilept, numt-hd, hd, npartcat, weight);
		}
		else if( hd == -999 ) {
		  cout << "\n\nBad top daughter id\n\n";
		}
	  }
	  else {
		numtaub = count_particles( numtaub, numb, 0, npartcat,  weight );
		int nhadz = numz - (numlepfromz - numnufromz)/2; //tau are leptonic
		numtauz = count_particles( numtauz, numtaufromz/2, nhadz, npartcat, weight);
		numtauh = count_particles( numtauh, numh, 0, npartcat,  weight );
		//if( (numtaufromz + numlepfromz + numnufromz)/2 != numz)
		//numtauz[3] += weight; //3rd indx is other--non-leptonic z

		int hd = num_hadronic_top( list_id );
		if( hd >= 0 ) {
		  numtaut = count_particles( numtaut, numt-hd, hd, npartcat, weight );
		  //if( hd > 0 ) 
		  //numtaut[3] += weight;
		}
		else if( hd == -999 ) {
		  cout << "\n\nBad top daughter id\n\n";
		}
	  }
	  numtotalgluon += numgluon*weight;
	  if( numgluon > 0 ) numevtgluon += weight;

	  //printing for tree test--version 2 only prints on error
	  //test_lists2( filenum, cms2.evt_run(), cms2.evt_event(), list_id);

	  if( cms2.genps_id()[6] == LSP_mcid ) nHardlsp++;
	  if( cms2.genps_id()[7] == LSP_mcid ) nHardlsp++;

	  //susy_pair_type returns index of hard_type array
	  int susy_pair_idx = susy_pair_type(Int_t(cms2.genps_id()[6]),
										 Int_t(cms2.genps_id()[7]) ); 	  
	  if( susy_pair_idx != 999 && chgflv != 999) {
		//weight is now a parameter of print_hard_type(...)
		hard_type[susy_pair_idx] ++;
		hard_chgflv[chgflv][susy_pair_idx] ++;
		if( numhadeta == 0 ) {
		  zjet_hard_type[susy_pair_idx] ++;
		  zjet_hard_chgflv[chgflv][susy_pair_idx] ++;
		}
	  }
	  else
		cout << "\n\nERROR: unknown susy hard scatter type"<<"\n\n";

	  //The index of genps is in idx_lep/had array
	  lep1 = cms2.genps_p4()[(int)idx_lep[0]]; //lep1 is always e|mu
	  if( numlep == 1 && pt_rank_lep[1] < pt_rank_tau[0] && TausAreLeptons)
		lep2 = cms2.genps_p4()[(int)idx_tau[0]];
	  else 
		lep2 = cms2.genps_p4()[(int)idx_lep[1]];
	  had1 = cms2.genps_p4()[(int)idx_had[0]];
	  had2 = cms2.genps_p4()[(int)idx_had[1]];

	  if( numbfromt != numteta )
		cout << "\n\nnumb from t not equal to numt\n\n";

	  //counting initial state quarks (indicies 4,5), only for ++,--:
	  //code for in_quark_typ:first index is:
	  //0:up-type  :2,4
	  //1:dbar-type:-1,-3
	  //2:do-type  :1,3
	  //3:ubar-type:-2,-4
	  //4:gluons   :21
	  //second index is 0 for pp, 1 for mm--same as buck_to_sign
	  unsigned int sign = buck_to_sign(bucket);
	  //if( sign != 999 && sign != 2) { 
	  if( sign < 2 ) { //exclude plus/minus, tau
		for( int i=4;i<6;i++) {
		  if(cms2.genps_id()[i]==2 || cms2.genps_id()[i]==4 )
			in_quark_typ[0][sign] += weight;
		  else if(cms2.genps_id()[i]==-1 || cms2.genps_id()[i]==-3 ||
				  cms2.genps_id()[i]==-5 )
			in_quark_typ[1][sign] += weight;
		  else if(cms2.genps_id()[i]==1 || cms2.genps_id()[i]==3 ||
				  cms2.genps_id()[i]==5 )
			in_quark_typ[2][sign] += weight;
		  else if(cms2.genps_id()[i]==-2 || cms2.genps_id()[i]==-4 )
			in_quark_typ[3][sign] += weight;
		  else if(cms2.genps_id()[i]==21 )
			in_quark_typ[4][sign] += weight;
		  else
			cout << "\nunknown initial state: " << cms2.genps_id()[i] << endl;
		}
	  }
	  else if( sign == 999 ) //tau handeled now
		cout << "\n\nERROR: unknown intitial state sign"<<"\n\n";

      //Lepton hists
	  double sumptlep = 0;
	  if( TausAreLeptons )
		sumptlep = pttau + ptlep;
	  else
		sumptlep = ptlep;

	  //these bucket checks should NO LONGER BE NECESSARY
	  //if( bucket != 999 ) //contents of buckets shouldn't change!
	  hPtLep[bucket]->Fill(sumptlep, weight );
	  if( numtau > 0 )
		hPtLep[tau]->Fill(pttau, weight );
	  if( numhadeta == 0 ) {
		hPtLep[zjet]->Fill(sumptlep, weight );
		hchgflv_PtLep_zjet[chgflv]->Fill(sumptlep, weight);
	  }
	  hPtLep[allBuckets-1]->Fill(sumptlep, weight );

	  //always put lep(e,mu) in Fist hist
	  hPtLepFirst[bucket]->Fill(pt_rank_lep[0], weight );
	  if( numtau > 0 )
		hPtLepFirst[tau]->Fill(pt_rank_tau[0], weight );
	  if( numhadeta == 0 ) {
		hPtLepFirst[zjet]->Fill(pt_rank_lep[0], weight );
		hchgflv_PtLepFirst_zjet[chgflv]->Fill(pt_rank_lep[0], weight );
	  }
	  hPtLepFirst[allBuckets-1]->Fill(pt_rank_lep[0], weight );

	  double* ptsecthr = pt_sec_thr(pt_rank_lep, pt_rank_tau, TausAreLeptons);

	  hPtLepSecond[bucket]->Fill(ptsecthr[0], weight );
	  if( pt_rank_tau[1] > 0 ) //tau bucket
		hPtLepSecond[tau]->Fill(pt_rank_tau[1], weight );
	  if( numhadeta == 0 ) {
		hPtLepSecond[zjet]->Fill(ptsecthr[0], weight );
		hchgflv_PtLepSecond_zjet[chgflv]->Fill(ptsecthr[0], weight );
	  }
	  hPtLepSecond[allBuckets-1]->Fill(ptsecthr[0], weight );

	  //don't fill if no third lepton
	  //if( pt_rank_lep[2] > 0 || (TausAreLeptons && numlep + numtau > 2) ) {
	  if( ptsecthr[1] > 0 ) {
		//if ( bucket != 999 && pt_rank_lep[2] > 0 ) 
		hPtLepThird[bucket]->Fill(ptsecthr[1], weight );
		if( pt_rank_tau[2] > 0 )
		  hPtLepThird[tau]->Fill(pt_rank_tau[2], weight );
		if( numhadeta == 0 ) {
		  hPtLepThird[zjet]->Fill(ptsecthr[1], weight );
		  hchgflv_PtLepThird_zjet[chgflv]->Fill(ptsecthr[1], weight );
		}
		hPtLepThird[allBuckets-1]->Fill(ptsecthr[1], weight );
	  }

	  //Note: check on mass for taus...
	  double dilmass = (lep1+lep2).M();
	  if( bucket != 999 )
		hMassLep[bucket]->Fill(dilmass, weight );
	  if( numhadeta == 0 )
		hMassLep[zjet]->Fill(dilmass, weight );
	  if(  numhadeta == 0 && bucket != 999 )
		hchgflv_MassLep_zjet[chgflv]->Fill(dilmass, weight);
	  hMassLep[allBuckets-1]->Fill(dilmass, weight );

	  if ( bucket != 999 ) 
		hNLep[bucket]->Fill(numlep, weight );
	  hNLep[nu]->Fill(numnu, weight );
	  hNLep[tau]->Fill(numtau, weight );
	  if( numhadeta == 0 )
		hNLep[zjet]->Fill(numlep, weight );
	  if(TausAreLeptons)
		hNLep[allBuckets-1]->Fill(numtau + numlep, weight );
	  else 
		hNLep[allBuckets-1]->Fill(numlep, weight );

      //Hadron hists
	  if( numhad > 0 ) { //don't fill for zero-jet
		//if( bucket != 999 ) 
		hSumEt[bucket]->Fill( pthad, weight );
		if( numb > 0 ) {
		  hSumEt[b]->Fill( ptb, weight );
		  hchgflv_SumEt_b[chgflv]->Fill( ptb, weight );
		}
		hSumEt[allBuckets-1]->Fill( pthad, weight );

		//if( bucket != 999 ) 
		hPtHadFirst[bucket]->Fill(pt_rank_had[0], weight );
		if( pt_rank_b[0] > 0 ) {
		  hPtHadFirst[b]->Fill(pt_rank_b[0], weight);
		  hchgflv_PtHadFirst_b[chgflv]->Fill(pt_rank_b[0], weight);
		}
		hPtHadFirst[allBuckets-1]->Fill(pt_rank_had[0], weight );

		//if( bucket != 999 ) 
		hPtHadSecond[bucket]->Fill(pt_rank_had[1], weight );
		if( pt_rank_b[1] > 0 ) {
		  hPtHadSecond[b]->Fill(pt_rank_b[1], weight);
		  hchgflv_PtHadSecond_b[chgflv]->Fill(pt_rank_b[1], weight);
		}
		hPtHadSecond[allBuckets-1]->Fill(pt_rank_had[1], weight );

		//if( bucket != 999 ) 
		hPtHadThird[bucket]->Fill(pt_rank_had[2], weight );
		if( pt_rank_b[2] > 0 ) {
		  hPtHadThird[b]->Fill(pt_rank_b[2], weight);
		  hchgflv_PtHadThird_b[chgflv]->Fill(pt_rank_b[2], weight);
		}
		hPtHadThird[allBuckets-1]->Fill(pt_rank_had[2], weight );
	  }

	  //if( bucket != 999 ) 
	  hNjets[bucket]->Fill(numhad, weight );
	  hNjets[b]->Fill(numb, weight );
	  hchgflv_Njets_b[chgflv]->Fill(numb, weight );
	  hNjets[t]->Fill(numt, weight );
	  hNjets[allBuckets-1]->Fill(numhad, weight );
	  
	  // alpha === Et2/Mjj, where Et2 is the Et of the jet closest to the MET
	  //alpha variable only defined for njets >= 2,
	  //but i invented a new variable, alphal, which is same but using
	  //the leptons instead of jets--only for zjet for now
	  if( numhad >= 2 ) {
		double alpha = get_alpha(had1, had2, vecEt);
		//double alphat = get_alphat(had1, had2, vecEt);
		if ( bucket != 999 ) {
		  halpha[bucket]->Fill(alpha, weight );
		  //halphat[bucket]->Fill(alphat, weight );
		}
		halpha[allBuckets-1]->Fill(alpha, weight );
		//halphat[allBuckets-1]->Fill(alphat, weight );
	  }
	  else if( numhadeta == 0 ) {
		double alphal = get_alpha(lep1, lep2, vecEt);
		halpha[zjet]->Fill(alphal, weight);
		hchgflv_alpha_zjet[chgflv]->Fill(alphal, weight);
	  }
	  //Next Iteration: Alphab --Alpha For B Jets
	  //in next iteration can calso do mass for bjets if desired

	  //MET hists
      hMET[LSP]->Fill(vecEtlsp.Et() , weight );
      hMET[nu]->Fill(vecEtnu.Et() , weight );
	  if( numhadeta == 0 )
		hMET[zjet]->Fill(vecEt.Et(), weight );
	  if ( bucket != 999 ) 
		hMET[bucket]->Fill(vecEt.Et(), weight );
	  if ( bucket != 999 && numhadeta == 0 )
		hchgflv_MET_zjet[chgflv]->Fill(vecEt.Et(), weight );
	  hMET[allBuckets-1]->Fill(vecEt.Et(), weight );

	  //Et hists
	  hEtFirst[LSP]->Fill(metlsp[0], weight );
	  hEtFirst[nu]->Fill(metnu[0], weight );
	  hEtSecond[LSP]->Fill(metlsp[1], weight );
	  hEtSecond[nu]->Fill(metnu[1], weight );

	  //MEff = mag(MET) + sumEt
	  if( bucket != 999 )
		hMEff[bucket]->Fill(vecEt.Et() + pthad, weight);
	  if( numhadeta == 0 )
		hMEff[zjet]->Fill(vecEt.Et() + pthad, weight);
	  hMEff[allBuckets-1]->Fill(vecEt.Et() + pthad, weight);

	  //Ht = mag(MET) + sumEt + ptLep
	  if( bucket != 999 && TausAreLeptons) 
		hHt[bucket]->Fill(vecEt.Et() + ptlep + pthad + pttau, weight);
	  else if( bucket != 999 )
		hHt[bucket]->Fill(vecEt.Et() + ptlep + pthad, weight);
	  if( numhadeta == 0 ) {
		hHt[zjet]->Fill(vecEt.Et() + ptlep, weight);
		hchgflv_Ht_zjet[chgflv]->Fill(vecEt.Et() + ptlep, weight);
	  }
	  if (TausAreLeptons) 
		hHt[allBuckets-1]->Fill(vecEt.Et() + ptlep + pthad + pttau, weight);
	  else
		hHt[allBuckets-1]->Fill(vecEt.Et() + ptlep + pthad, weight);

	  //cout << "Done fill hists\n";
    }
	//end event loop
  }
  //end file loop

  cout << "\n\n";

  std::cout << std::endl;
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  
  std::cout << "Prefix: " << prefix << " processed: " << nEventsTotal
	    << " Events, found: " << duplicates_total_n 
	    << " Duplicates and selected: " << nCandidatesSelected
	    << " Candidates." << endl << endl;
  //cout << "evt_scale1fb = " << cms2.evt_scale1fb() << endl;
  //if( nbaddecay > 0 ) cout << "\n\nBad decays found : " << nbaddecay*weight << "\n\n";
  
  //function chgflv_plots(...) takes a vector of hists for a single variable and returns a single hist for chargeflav separated
  TH1F* hchgflv_PtLepFirst[numchgflv];
  TH1F* hchgflv_PtLepSecond[numchgflv];
  TH1F* hchgflv_PtLepThird[numchgflv];
  TH1F* hchgflv_PtLep[numchgflv];
  TH1F* hchgflv_MassLep[numchgflv];
  TH1F* hchgflv_MET[numchgflv];
  TH1F* hchgflv_SumEt[numchgflv];
  TH1F* hchgflv_MEff[numchgflv];
  TH1F* hchgflv_Ht[numchgflv];
  TH1F* hchgflv_Njets[numchgflv];
  TH1F* hchgflv_alpha[numchgflv];
  
  //loop over numchgflv=(4||6 depending on tau) chgflv: 
  for(unsigned int i=0;i<numchgflv;i++) {
	hchgflv_PtLepFirst[i] = chgflv_plots(hPtLepFirst, allBuckets, i);
	hchgflv_PtLepSecond[i] = chgflv_plots(hPtLepSecond, allBuckets, i);
	hchgflv_PtLepThird[i] = chgflv_plots(hPtLepThird, allBuckets, i);
	hchgflv_PtLep[i] = chgflv_plots(hPtLep, allBuckets, i);
	hchgflv_MassLep[i] = chgflv_plots(hMassLep, allBuckets, i);
	hchgflv_MET[i] = chgflv_plots(hMET, allBuckets, i);
	hchgflv_SumEt[i] = chgflv_plots(hSumEt, allBuckets, i);
	hchgflv_MEff[i] = chgflv_plots(hMEff, allBuckets, i);
	hchgflv_Ht[i] = chgflv_plots(hHt, allBuckets, i);
	hchgflv_Njets[i] = chgflv_plots(hNjets, allBuckets, i);
	hchgflv_alpha[i] = chgflv_plots(halpha, allBuckets, i);
  }


  //write all hists to root file
  TFile outfile(Form("%s.root",name), "RECREATE");
  cout << "\nWriting hists, name: " << name << "\n\n";

  for(unsigned int j = 0; j < allBuckets; j++) {
    hPtLepFirst[j]->Write();
    hPtLepSecond[j]->Write();
    hPtLepThird[j]->Write();
    hPtLep[j]->Write();
	hMassLep[j]->Write();
	hNLep[j]->Write();

	hSumEt[j]->Write();
    hPtHadFirst[j]->Write();
    hPtHadSecond[j]->Write();
    hPtHadThird[j]->Write();
	hNjets[j]->Write();
	halpha[j]->Write();
	//halphat[j]->Write();

	hMET[j]->Write();
	hEtFirst[j]->Write();
	hEtSecond[j]->Write();

	hMEff[j]->Write();
	hHt[j]->Write();

	if(j < 4) {
	  hchgflv_PtLepFirst[j]->Write();
	  hchgflv_PtLepSecond[j]->Write();
	  hchgflv_PtLepThird[j]->Write();
	  hchgflv_PtLep[j]->Write();
	  hchgflv_MassLep[j]->Write();
	  hchgflv_MET[j]->Write();
	  hchgflv_SumEt[j]->Write();
	  hchgflv_MEff[j]->Write();
	  hchgflv_Ht[j]->Write();
	  hchgflv_Njets[j]->Write();
	  hchgflv_alpha[j]->Write();

	  hchgflv_PtLepFirst_zjet[j]->Write();
	  hchgflv_PtLepSecond_zjet[j]->Write();
	  hchgflv_PtLepThird_zjet[j]->Write();
	  hchgflv_PtLep_zjet[j]->Write();
	  hchgflv_MassLep_zjet[j]->Write();
	  hchgflv_MET_zjet[j]->Write();
	  hchgflv_Ht_zjet[j]->Write();
	  hchgflv_alpha_zjet[j]->Write();

	  hchgflv_PtHadFirst_b[j]->Write();
	  hchgflv_PtHadSecond_b[j]->Write();
	  hchgflv_PtHadThird_b[j]->Write();
	  hchgflv_SumEt_b[j]->Write();
	  hchgflv_Njets_b[j]->Write();
	}
  }

  //if( nHardlsp > 0 ) cout << "Number hard scatter lsp: " << nHardlsp << "\n\n";

  //Write tables
  ofstream file( Form("%s_table",name) );
  file << "---+ " << prefix << "\n\n";
  file << print_cuts();

  file << "\nProduction X-section (fb), All events:\n";
  file << print_hard_type(hard_type, hard_chgflv, susy_types, weight, numchgflv);
  file << "\nProduction X-section (fb), Zero Jet Bin:\n";
  file << print_hard_type(zjet_hard_type, zjet_hard_chgflv, susy_types, weight, numchgflv);

  //put 'all' bucket of variables for which I want 90% efficiency points
  int numvars = 11;
  TH1F allhists[numvars];
  allhists[0] = dilep_sum(hchgflv_PtLepFirst);
  allhists[1] = dilep_sum(hchgflv_PtLepSecond);
  allhists[2] = dilep_sum(hchgflv_PtLepThird);
  allhists[3] = dilep_sum(hchgflv_PtLep);
  allhists[4] = dilep_sum(hchgflv_MassLep);
  allhists[5] = dilep_sum(hchgflv_MET);
  allhists[6] = dilep_sum(hchgflv_SumEt);
  allhists[7] = dilep_sum(hchgflv_MEff);
  allhists[8] = dilep_sum(hchgflv_Ht);
  allhists[9] = dilep_sum(hchgflv_Njets);
  allhists[10]= dilep_sum(hchgflv_alpha);

  TH1F tauhists[numvars];
  if( TausAreLeptons ) {
	tauhists[0] = tau_sum(hchgflv_PtLepFirst);
	tauhists[1] = tau_sum(hchgflv_PtLepSecond);
	tauhists[2] = tau_sum(hchgflv_PtLepThird);
	tauhists[3] = tau_sum(hchgflv_PtLep);
	tauhists[4] = tau_sum(hchgflv_MassLep);
	tauhists[5] = tau_sum(hchgflv_MET);
	tauhists[6] = tau_sum(hchgflv_SumEt);
	tauhists[7] = tau_sum(hchgflv_MEff);
	tauhists[8] = tau_sum(hchgflv_Ht);
	tauhists[9] = tau_sum(hchgflv_Njets);
	tauhists[10]= tau_sum(hchgflv_alpha);
  }
  
  TH1F zjethists[numvars];
  zjethists[0] = dilep_sum(hchgflv_PtLepFirst_zjet);
  zjethists[1] = dilep_sum(hchgflv_PtLepSecond_zjet);
  zjethists[2] = dilep_sum(hchgflv_PtLepThird_zjet);
  zjethists[3] = dilep_sum(hchgflv_PtLep_zjet);
  zjethists[4] = dilep_sum(hchgflv_MassLep_zjet);
  zjethists[5] = dilep_sum(hchgflv_MET_zjet);
  zjethists[6] = dilep_sum(hchgflv_SumEt_zjet);
  zjethists[7] = dilep_sum(hchgflv_MEff_zjet);
  zjethists[8] = dilep_sum(hchgflv_Ht_zjet);
  zjethists[9] = dilep_sum(hchgflv_Njets_zjet);
  zjethists[10]= dilep_sum(hchgflv_alpha_zjet);

  TH1F zjettauhists[numvars];
  if( TausAreLeptons ) {
	zjettauhists[0] = tau_sum(hchgflv_PtLepFirst_zjet);
	zjettauhists[1] = tau_sum(hchgflv_PtLepSecond_zjet);
	zjettauhists[2] = tau_sum(hchgflv_PtLepThird_zjet);
	zjettauhists[3] = tau_sum(hchgflv_PtLep_zjet);
	zjettauhists[4] = tau_sum(hchgflv_MassLep_zjet);
	zjettauhists[5] = tau_sum(hchgflv_MET_zjet);
	zjettauhists[6] = tau_sum(hchgflv_SumEt_zjet);
	zjettauhists[7] = tau_sum(hchgflv_MEff_zjet);
	zjettauhists[8] = tau_sum(hchgflv_Ht_zjet);
	zjettauhists[9] = tau_sum(hchgflv_Njets_zjet);
	zjettauhists[10]= tau_sum(hchgflv_alpha_zjet);
  }

  int numbvars = 5;
  TH1F bhists[numbvars];
  bhists[0] = dilep_sum(hchgflv_PtHadFirst_b);
  bhists[1] = dilep_sum(hchgflv_PtHadSecond_b);
  bhists[2] = dilep_sum(hchgflv_PtHadThird_b);
  bhists[3] = dilep_sum(hchgflv_SumEt_b);
  bhists[4] = dilep_sum(hchgflv_Njets_b);

  TH1F btauhists[numbvars];
  if( TausAreLeptons ) {
	btauhists[0] = tau_sum(hchgflv_PtHadFirst_b);
	btauhists[1] = tau_sum(hchgflv_PtHadSecond_b);
	btauhists[2] = tau_sum(hchgflv_PtHadThird_b);
	btauhists[3] = tau_sum(hchgflv_SumEt_b);
	btauhists[4] = tau_sum(hchgflv_Njets_b);
  }
  
  //print table of 90% Efficiency Points for all and for b's
  file << "\nKinematic plots with 90\% Efficiency Points:";
  if( !TausAreLeptons ) {
	file << print_90eff(allhists, zjethists, numvars );
	file << "\nb-quark hists with 90% Efficiency Points:";
    file << print_90eff(bhists, numbvars );
  }
  else {
	file << print_90eff(allhists, tauhists, zjethists, zjettauhists, numvars );
	file << "\nb-quark hists with 90% Efficiency Points:";
	file << print_90eff(bhists, btauhists, numbvars );	
  }

  //print numbers, num events for certain particles:
  file << print_particle_table(npartcat, numdilepb, numtaub, numdilept, numtaut, numdilepz, numtauz, numdileph, numtauh);

  file << "\nSelected particle info:\n";
  file << "|**|*tau in dilep*|*tau not dilep*|*gluon*|\n";
  file << "|num particles|" << round_char(numdileptau) << "|" << round_char(numtotaltau) << "|" << round_char(numtotalgluon) << "|\n";
  file << "|num events|" << round_char(numevtdileptau) << "|" << round_char(numevttau) << "|" << round_char(numevtgluon) << "|";
  file << "\n\n";

  //print initial state quark charge table
  file << "\nInitial State Quark Charge Table:\n";
  file << "|*lep sign*|*up-type*|*anti-<br>down-type*|*down-type*|*anti-<br>up-type*|*gluon*|*Total*|";
  for(int j=0;j<2;j++){
	if(j==0) file << "\n|plus-plus|";
	else file << "\n|minus-minus|";

	double total = 0;
	for(int i=0;i<5;i++){
	  file << round_char(in_quark_typ[i][j]) << "|";
	  total += in_quark_typ[i][j];
	}
	file << round_char(total) << "|";
  }
  file << "\n\n";

  //ntuple cross-section: not all events make it into ntuples, not all events in ntuples are dilepton.
  file << "\nTotal ntuple cross-section (before dilepton selection): " << nEventsTotal*weight << "\n\n";

  THStack* hs_PtLepFirst[2];
  THStack* hs_PtLepSecond[2];
  THStack* hs_PtLepThird[2];
  THStack* hs_PtLep[2];
  THStack* hs_MassLep[2];
  THStack* hs_MET[2]; 
  THStack* hs_SumEt[2];
  THStack* hs_MEff[2];
  THStack* hs_Ht[2];
  THStack* hs_Njets[2];
  THStack* hs_alpha[2];
  
  THStack* hs_PtLepFirst_zjet[2];
  THStack* hs_PtLepSecond_zjet[2];
  THStack* hs_PtLepThird_zjet[2];
  THStack* hs_PtLep_zjet[2];
  THStack* hs_MassLep_zjet[2];
  THStack* hs_MET_zjet[2];
  THStack* hs_Ht_zjet[2];
  THStack* hs_alpha_zjet[2];
  
  THStack* hs_PtHadFirst_b[2];
  THStack* hs_PtHadSecond_b[2];
  THStack* hs_PtHadThird_b[2];
  THStack* hs_SumEt_b[2];
  THStack* hs_Njets_b[2];

  for(int i=0;i<2;i++) {
	hs_PtLepFirst[i] = make_stack(hchgflv_PtLepFirst, numchgflv,i);
	hs_PtLepSecond[i] = make_stack(hchgflv_PtLepSecond, numchgflv, i);
	hs_PtLepThird[i] = make_stack(hchgflv_PtLepThird, numchgflv, i);
	hs_PtLep[i] = make_stack(hchgflv_PtLep, numchgflv, i);
	hs_MassLep[i] = make_stack(hchgflv_MassLep, numchgflv, i);
	hs_MET[i] = make_stack(hchgflv_MET, numchgflv, i);
	hs_SumEt[i] = make_stack(hchgflv_SumEt, numchgflv, i);
	hs_MEff[i] = make_stack(hchgflv_MEff, numchgflv, i); 
	hs_Ht[i] = make_stack(hchgflv_Ht, numchgflv, i);	  
	hs_Njets[i] = make_stack(hchgflv_Njets, numchgflv, i);
	hs_alpha[i] = make_stack(hchgflv_alpha, numchgflv, i);

	hs_PtLepFirst_zjet[i] = make_stack(hchgflv_PtLepFirst_zjet, numchgflv, i);
	hs_PtLepSecond_zjet[i] = make_stack(hchgflv_PtLepSecond_zjet, numchgflv,i);
	hs_PtLepThird_zjet[i] = make_stack(hchgflv_PtLepThird_zjet, numchgflv, i);
	hs_PtLep_zjet[i] = make_stack(hchgflv_PtLep_zjet, numchgflv, i);
	hs_MassLep_zjet[i] = make_stack(hchgflv_MassLep_zjet, numchgflv, i);
	hs_MET_zjet[i] = make_stack(hchgflv_MET_zjet, numchgflv, i);
	hs_Ht_zjet[i] = make_stack(hchgflv_Ht_zjet, numchgflv, i);	  
	hs_alpha_zjet[i] = make_stack(hchgflv_alpha_zjet, numchgflv, i);

	hs_PtHadFirst_b[i] = make_stack(hchgflv_PtHadFirst_b, numchgflv, i);
	hs_PtHadSecond_b[i] = make_stack(hchgflv_PtHadSecond_b, numchgflv, i);
	hs_PtHadThird_b[i] = make_stack(hchgflv_PtHadThird_b, numchgflv, i);
	hs_SumEt_b[i] = make_stack(hchgflv_SumEt_b, numchgflv, i);
	hs_Njets_b[i] = make_stack(hchgflv_Njets_b, numchgflv, i);
  }
  for(int i=0;i<2;i++) {
	hs_PtLepFirst[i]->Write();
	hs_PtLepSecond[i]->Write();
	hs_PtLepThird[i]->Write();
	hs_PtLep[i]->Write();
	hs_MassLep[i]->Write();
	hs_MET[i]->Write();
	hs_SumEt[i]->Write();
	hs_MEff[i]->Write();
	hs_Ht[i]->Write();
	hs_Njets[i]->Write();
	hs_alpha[i]->Write();

	hs_PtLepFirst_zjet[i]->Write();
	hs_PtLepSecond_zjet[i]->Write();
	hs_PtLepThird_zjet[i]->Write();
	hs_PtLep_zjet[i]->Write();
	hs_MassLep_zjet[i]->Write();
	hs_MET_zjet[i]->Write();
	hs_Ht_zjet[i]->Write();
	hs_alpha_zjet[i]->Write();

	hs_PtHadFirst_b[i]->Write();
	hs_PtHadSecond_b[i]->Write();
	hs_PtHadThird_b[i]->Write();
	hs_SumEt_b[i]->Write();
	hs_Njets_b[i]->Write();
  }
  //i'm ditching the overlays for now b'c tau messed them up
  //if i restore them, will have to restore display in print_90eff (3 versions)
  //save_overlay(hPtLepFirst[allBuckets-1], hPtLepSecond[allBuckets-1], hPtLepThird[allBuckets-1]);
  //save_overlay(hPtLepFirst[zjet], hPtLepSecond[zjet], hPtLepThird[zjet]);
  //save_overlay(hPtHadFirst[b], hPtHadSecond[b], hPtHadThird[b]);

  outfile.Close();
  file.close();

  return name;
}
//end char* ScanChain


const char* print_cuts() {
  ostringstream stream;
  stream.setf(ios::fixed);

  stream << "Cuts: ";
  if( EtaCut ) stream << "Eta < " << round_char(EtaCutValue) << ",  ";
  if( PtCut ) stream << "Pt(jets) > " << PtCutValue << ",  ";
  if( TausAreLeptons ) stream << " including taus \n";
  else stream << " excluding taus \n";

  return stream.str().c_str();
}
//end const char* print_cuts()

char* name_from_cuts(char* prefix, char* tag){ 

  char* name = new char[100];
  strcpy(name, prefix);
  strcat(name, "_");
  strcat(name, tag);
  if( EtaCut ) {
	if( EtaCutValue == 2.5 )
	  strcat(name, "_eta25");
	else if( EtaCutValue == 3.0 )
	  strcat(name, "_eta30");
	else
	  strcat(name, "_eta");
  }
  if( PtCut ) {
	if( PtCutValue == 50 )
	  strcat(name, "_pt50");
	else
	  strcat(name, "_pt");
  }
  if( TausAreLeptons ) strcat(name, "_tau");
  //strcat(name, ".root");//put .filetype in looper

  return name;
}
//end char* name_from_cuts

