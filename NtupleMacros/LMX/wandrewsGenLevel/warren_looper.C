
//now make the source file
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
//#include <cstring>

using namespace std;

#ifndef __CINT__
#include "../../WZ/CMS2_V00_04_00.h"
CMS2 cms2;
#endif

#include "../../CORE/selections.cc"
#include "../../CORE/utilities.cc"
#include "../../Tools/tools.cc"
#include "warren_functions.C"
//#include "warren_hists.C"
//gROOT->ProcessLine(".L warren_hists.C+");

//wandrews: This was originally Oli's looper. I'm butchering (changing) it to be a dilepton Generator level looper for use in my susy stuff

//notes to myself:
//***Make sure twiki is up to date with all detail, definitions, etc

//FKW: Lets stick to dilepton-->toss all events which have < 2 leptons. > 2 leptons is no problem--we'll mainly focus on the greatest two Pt leptons

//idea behind buckets: buckets are for final state (e/mu +/-) but also for other particles in event like LSP,nu,b,t,q
//we should use buckets for all variables: if there is a variable worth looking at for one particle type, it may as well be plotted for all types.
//  g,q != t (had or jet)
//  LSP
//  neutrinos
//  b,t

//function definitions in "warren_functions.C"
char* name_from_cuts(char* prefix, char* tag);
void print_entries_table(TH1F* hist[], int start, int stop);

//global vars for cuts so i don't have to pass them all over the place
bool EtaCut = false;
double EtaCutValue = 2.5;
bool TausAreLeptons = false;

//  bool useMetCut = true;
//  bool usePtCut = false;
//  bool useIdCut = true;

//LSP's mcid. This will depend on which LM I use
//for now, just LM1
const Int_t LSP_mcid = 1000022; //chi_10

//see notes above for bucket explanation
const unsigned int lepBuckets = 10;
const unsigned int othBuckets = 8; //tau, lsp, neutrino, zjet, jet, b, t, zlep
const unsigned int allBuckets = lepBuckets + othBuckets + 1;
const unsigned int tau = 10;
const unsigned int nu = 11;
const unsigned int LSP = 12;
const unsigned int zjet = 13;
const unsigned int zlep = 14;//ZERO LEPTON ARE NOT IN 'ALL' BUCKET
const unsigned int jet = 15;
const unsigned int b = 16;
const unsigned int t = 17;

char *suffix[allBuckets];
//MAKE NEW BUCKET: hard scatter susy--pair produced sparticle
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
  suffix[tau] = "tau";
  suffix[nu] = "nu";
  suffix[LSP] = "LSP";
  suffix[zjet] = "zero-jet";
  suffix[zlep] = "zero-lep"; //do not include zero-lep in 'all' bucket
  suffix[jet] = "jet"; //jet=g,q!=t 
  suffix[b] = "b";
  suffix[t] = "t"; // top bucket
  suffix[allBuckets - 1] = "all";
}


//Start Oli's old code:

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
  float weight=0;

  int bad_id = 0;
  int numidtau = 0; //count taus like counting bad ids, compare to numtau
  int numtotaltau = 0;
  int ngluons = 0;
  int nzerojet = 0;
  int nHardlsp = 0;
  //int dislep_ss = 0; int dislep_os = 0;
  //just number of pair produced types, see warren_functions.C
  const int susy_types = 10;//gg,qg,qq,qX,qs,XX,sX,ss,gX,hh
  double* hard_type = new double[susy_types];
  double* zjet_hard_type = new double[susy_types];
  //this is for charge/flavor of susy pair produced
  //first index is sign/flavor (same/opp), second susy type
  //double** hard_chgflv = new double[4][susy_types];
  double hard_chgflv[4][susy_types];
  double zjet_hard_chgflv[4][susy_types];
  for(int i=0; i<susy_types;i++) {
	hard_type[i] = 0; zjet_hard_type[i] = 0;
	for(int j=0;j<4;j++) {
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
  TH1F* hNjets[allBuckets]; //num g,q!=t

  //Et, MET hists
  //MET = magnitude Et vector:
  //  LSP bucket is filled with vector sum of Et of LSPs
  //No. Changed my mind. 2Et hists:first and second. No scalar sums.
  TH1F* hEtFirst[allBuckets];
  TH1F* hEtSecond[allBuckets];
  TH1F* hPhi[allBuckets];

  //magnitude of vector sum of Et of LSP+nu
  TH1F* hMET[allBuckets];
  //just for checking lsp met:Et minus Pt, just for LSP
  TH1F* hEtmPt[allBuckets];
  //MEff = mag(met) + SumEt
  TH1F* hMEff[allBuckets];
  //Ht = SumEt + mag(met) + ptLep
  TH1F* hHt[allBuckets];
  


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
							   nBinsF,lowBin,highBinF,"p_{T} [GeV]","");   
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

    const float highBinB = 800;
    const int nBinsB = 80;

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
	
    // Et, MET 
    const float highBinM = 600;
    const int nBinsM = 120;
    const float highBinH = 1000;
    const int nBinsH = 100;
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

    hEtmPt[i] = book1DHist(Form("%s_hEtmPt_%s",prefix,suffix[i]),Form("%s_hEtmPt_%s",prefix,suffix[i]),nBinsS,lowBin,highBinS,"Et minus Pt [GeV]","");	
  }
  //end for i < allBuckets = ?

  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
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

	  Int_t* list_id = new Int_t[100];
	  for(int i=0;i<100;i++) { list_id[i]=0; }

	  double* id_lep = new double[3];
	  double* idx_lep = new double[3]; //for 4-vectors
	  double* pt_rank_lep = new double[3];
	  double* pt_rank_leptau = new double[3];
	  double* pt_rank_had = new double[3]; 
	  double* pt_rank_b = new double[3]; 
	  double* pt_rank_tau = new double[3]; 
      double* metlsp = new double[2]; //guaranteed to have 2 lsp in event
	  double* metnu = new double[3]; //forget about the rest (if > 3)
	  for(int i=0;i<3;i++) {
		id_lep[i]=0;
		idx_lep[i]=0;
		pt_rank_lep[i]=0;
		pt_rank_leptau[i]=0;
		pt_rank_had[i]=0;
		pt_rank_b[i]=0;
		pt_rank_tau[i]=0;
		if(i<2) metlsp[i]=0;
		metnu[i]=0;
	  }

      int bucket = 0;
      int numlep = 0; //NOT including taus
	  int numnu = 0; int numtau = 0;	  
	  int numhad = 0;//numhad includes b, not t
	  int numb = 0; int numt = 0;
      double pthad = 0; double ptlep = 0; double pttau = 0;
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vecEt;
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vecEtlsp;
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vecEtnu;
	  vecEt.SetPxPyPzE(0,0,0,0);
	  vecEtlsp.SetPxPyPzE(0,0,0,0);
	  vecEtnu.SetPxPyPzE(0,0,0,0);
	  //lepton 4vectors
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lep1;
	  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lep2;
	  lep1.SetPxPyPzE(0,0,0,0);
	  lep2.SetPxPyPzE(0,0,0,0);

	  for( int i=0;i<6;i++) {
		list_id[i]=cms2.genps_id()[i];
	  }
	  //loop starts at 6 b'c 0,1 are protons, 2-5 are hard scatter in proton
      for ( unsigned int par = 6; 
			par < cms2.genps_id().size(); ++par ) {

		double id = cms2.genps_id()[par];
		//if(event <= 4) { cout << par << ":  "  << id << endl; }
		//if( numhad == 0 )
		list_id[par] = Int_t(id);

		//quark!=t, t=6
		if( (abs(id) >= 1 && abs(id) <= 5) || abs(id) == 21  ) { 
		  if( EtaCut && abs(cms2.genps_p4()[par].Eta()) > EtaCutValue )
			continue;

		  //scalar sum of g,q!=t pt
		  pthad += cms2.genps_p4()[par].pt();
		  numhad++;
		  if( abs(id) == 5) {
			numb++;
			pt_rank_b = get_pt_rank(pt_rank_b, cms2.genps_p4()[par].pt() );
		  }

		  pt_rank_had = get_pt_rank(pt_rank_had, cms2.genps_p4()[par].pt() );
		}
		else if( abs(id) == 6 ) { //t quark
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
		  //check results of idx_lep???
		  if( event < 5 ) { cout << ""; }

		  pt_rank_lep = get_pt_rank(pt_rank_lep, cms2.genps_p4()[par].pt() );
		  //This one always filled b'c it teats tausarelep=true always
		  pt_rank_leptau = get_pt_rank(pt_rank_leptau,
									   cms2.genps_p4()[par].pt() );
		}
		else if( abs(id) == 15 ) {  //tau
		  if( EtaCut && abs(cms2.genps_p4()[par].Eta()) > EtaCutValue )
			continue;

		  pttau += cms2.genps_p4()[par].pt();
		  numtau++;

		  pt_rank_tau = get_pt_rank(pt_rank_tau, cms2.genps_p4()[par].pt() );
		  pt_rank_leptau = get_pt_rank(pt_rank_leptau,
									   cms2.genps_p4()[par].pt());
		}
		else if( abs(id) == 12 || abs(id) == 14 || abs(id) == 16 ) {
		  numnu++;
		  metnu = get_pt_rank(metnu, cms2.genps_p4()[par].Et() );
		  vecEt += cms2.genps_p4()[par];
		  vecEtnu += cms2.genps_p4()[par];
		  //hPhi[nu]->Fill(cms2.genps_p4()[par].Phi() , weight);
		}
		else if( abs(id) == LSP_mcid ) {
		  metlsp = get_Et_rank(metlsp, cms2.genps_p4()[par].Et() );
		  vecEt += cms2.genps_p4()[par];
		  vecEtlsp += cms2.genps_p4()[par];
		  hPhi[LSP]->Fill(cms2.genps_p4()[par].Phi() , weight);
		  hEtmPt[LSP]->Fill(cms2.genps_p4()[par].Et() - cms2.genps_p4()[par].pt() , weight);
		}
		//only gluons found so far are from h0->gg, h0 from X20->h0X10 I thnk.
		else if( abs(id) == 21 ) {
		  //ngluons++;
		  //in case I want to find tree w/ dave's cmssw filter
		  //cout << "event: " << cms2.evt_event() << " run: " << cms2.evt_run() << endl;
		}
      }
	  //end genps loop

	  //CUT: must have number leptons >= 2
      if(numlep < 2){
		bool die = false;
		//cout << "numlep cut\n";
		if( TausAreLeptons && ( numtau + numlep ) < 2 )
		  die = true; 
		else if( !TausAreLeptons )
		  die = true;

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

	  //check results of idx_lep
	  //if( event < 5 ) {
	  //	cout << endl;
	  //	for(int i=0;i<3;i++){
	  //	  cout << pt_rank_lep[i] << "  " << idx_lep[i] << endl;
	  //	}
	  //	for(int i=0; i < 100; i++) {
	  //	  if( list_id[i] != 0) 
	  //		cout << i << " : " << list_id[i] << endl;
	  //	}
	  //}

	  //The index of genps is in idx_lep array
	  lep1 = cms2.genps_p4()[(int)idx_lep[0]];
	  lep2 = cms2.genps_p4()[(int)idx_lep[1]];

	  //cout << "get bucket\n";
      //determine bucket for lepton vars
	  bucket = get_lep_bucket(id_lep[0], id_lep[1]);
	  //************
	  //currently, events with <2 e,mu and with taus will get bucket 999
	  if( bucket == 999 ) {
		//cout << "ERROR: Can't determine lepton bucket. Continue;\n";
		if(event < 5) {
		  //cout << "ids failed: " << id_lep[0] << "," << id_lep[1] << endl;
		}

		//i'm assuming that the zeros are taus--will want to check later
		if( abs(id_lep[0]) == 0 || abs(id_lep[1]) == 0 )
		  numidtau++;
		else
		  bad_id++;
		//continue;
	  }

	  //this part is just for testing the zero jet bin mcid
	  if( numhad == 0 ) {
		nzerojet++;
		//if(nzerojet < 5 ) {
		//  cout << endl;
		//  for(int i=0; i < 100; i++) {
		//	if( list_id[i] != 0) 
		//	  cout << i+6 << " : " << list_id[i] << endl;
		//  }
		//}
	  }//end zero jet bin test
	  //uncomment below to require zero-jet
	  //else
		//continue;

	  ++nCandidatesSelected;      	  
	  
	  numtotaltau += numtau;

	  if( cms2.genps_id()[6] == LSP_mcid ) nHardlsp++;
	  if( cms2.genps_id()[7] == LSP_mcid ) nHardlsp++;

	  //susy_pair_type returns index of hard_type array
	  int susy_pair_idx = susy_pair_type(Int_t(cms2.genps_id()[6]),
								 Int_t(cms2.genps_id()[7]) );
	  //chgflv is for sign,flavor combo of final state
	  int chgflv = buck_to_chgflv(bucket);
 	  
	  if( susy_pair_idx != 999 && chgflv != 999) {
		//weight is now a parameter of print_hard_type(...)
		hard_type[susy_pair_idx] ++;
		hard_chgflv[chgflv][susy_pair_idx] ++;
		if( numhad == 0 ) {
		  zjet_hard_type[susy_pair_idx] ++;
		  zjet_hard_chgflv[chgflv][susy_pair_idx] ++;
		}

	  }
	  else
		cout << "\n\nERROR: unknown susy hard scatter type"<<"\n\n";
	  
	  //cout << "fill hists\n";
      //Lepton hists
	  if ( bucket != 999 ) 
		hPtLep[bucket]->Fill(ptlep, weight );
      hPtLep[tau]->Fill(pttau, weight );
	  if( numhad == 0 )
		hPtLep[zjet]->Fill(ptlep, weight );
	  if( TausAreLeptons )
		hPtLep[allBuckets-1]->Fill(pttau + ptlep, weight );
	  else 
		hPtLep[allBuckets-1]->Fill(ptlep, weight );
      
	  if ( bucket != 999 ) 
		hPtLepFirst[bucket]->Fill(pt_rank_lep[0], weight );
	  hPtLepFirst[tau]->Fill(pt_rank_tau[0], weight );
	  if( numhad == 0 )
		hPtLepFirst[zjet]->Fill(pt_rank_lep[0], weight );
	  if(TausAreLeptons)
		hPtLepFirst[allBuckets-1]->Fill(pt_rank_leptau[0], weight );
	  else 
		hPtLepFirst[allBuckets-1]->Fill(pt_rank_lep[0], weight );

	  if ( bucket != 999 ) 
		hPtLepSecond[bucket]->Fill(pt_rank_lep[1], weight );
	  hPtLepSecond[tau]->Fill(pt_rank_tau[1], weight );
	  if( numhad == 0 )
		hPtLepSecond[zjet]->Fill(pt_rank_lep[1], weight );
	  if(TausAreLeptons)
		hPtLepSecond[allBuckets-1]->Fill(pt_rank_leptau[1], weight );
	  else 
		hPtLepSecond[allBuckets-1]->Fill(pt_rank_lep[1], weight );

	  //don't fill if no third lepton
	  if(pt_rank_lep[2] > 0) { 
		if ( bucket != 999 ) 
		  hPtLepThird[bucket]->Fill(pt_rank_lep[2], weight );
		hPtLepThird[tau]->Fill(pt_rank_tau[2], weight );
		if( numhad == 0 )
		  hPtLepThird[zjet]->Fill(pt_rank_lep[2], weight );
		if(TausAreLeptons)
		  hPtLepThird[allBuckets-1]->Fill(pt_rank_leptau[2], weight );
		else 
		  hPtLepThird[allBuckets-1]->Fill(pt_rank_lep[2], weight );
	  }

	  //Note: NO mass for taus so far...
	  double dilmass = (lep1+lep2).M();
	  if( bucket != 999 )
		hMassLep[bucket]->Fill(dilmass, weight );
	  if( numhad == 0 )
		hMassLep[zjet]->Fill(dilmass, weight );
	  hMassLep[allBuckets-1]->Fill(dilmass, weight );

	  if ( bucket != 999 ) 
		hNLep[bucket]->Fill(numlep, weight );
	  hNLep[nu]->Fill(numnu, weight );
	  hNLep[tau]->Fill(numtau, weight );
	  if( numhad == 0 )
		hNLep[zjet]->Fill(numlep, weight );
	  if(TausAreLeptons)
		hNLep[allBuckets-1]->Fill(numtau + numlep, weight );
	  else 
		hNLep[allBuckets-1]->Fill(numlep, weight );

      //Hadron hists
	  hSumEt[jet]->Fill( pthad, weight ); //all had
	  if ( bucket != 999 ) 
		hSumEt[bucket]->Fill( pthad, weight ); 
	  hSumEt[allBuckets-1]->Fill( pthad, weight );

	  if ( bucket != 999 ) 
		hPtHadFirst[bucket]->Fill(pt_rank_had[0], weight );
	  hPtHadFirst[b]->Fill(pt_rank_b[0], weight);
	  hPtHadFirst[allBuckets-1]->Fill(pt_rank_had[0], weight );

	  if ( bucket != 999 ) 
		hPtHadSecond[bucket]->Fill(pt_rank_had[1], weight );
      hPtHadSecond[allBuckets-1]->Fill(pt_rank_had[1], weight );
	  hPtHadSecond[b]->Fill(pt_rank_b[1], weight);

	  if ( bucket != 999 ) 
		hPtHadThird[bucket]->Fill(pt_rank_had[2], weight );
      hPtHadThird[allBuckets-1]->Fill(pt_rank_had[2], weight );
	  hPtHadThird[b]->Fill(pt_rank_b[2], weight);

	  hNjets[jet]->Fill(numhad, weight );
	  if ( bucket != 999 ) 
		hNjets[bucket]->Fill(numhad, weight );
	  hNjets[b]->Fill(numb, weight );
	  hNjets[t]->Fill(numt, weight );
	  hNjets[allBuckets-1]->Fill(numhad, weight );

	  //MET hists
      hMET[LSP]->Fill(vecEtlsp.Et() , weight );
      hMET[nu]->Fill(vecEtnu.Et() , weight );
	  if( numhad == 0 )
		hMET[zjet]->Fill(vecEt.Et(), weight );
	  if ( bucket != 999 ) 
		hMET[bucket]->Fill(vecEt.Et(), weight );
	  hMET[allBuckets-1]->Fill(vecEt.Et(), weight );

	  //Et hists
	  hEtFirst[LSP]->Fill(metlsp[0], weight );
	  hEtFirst[nu]->Fill(metnu[0], weight );
	  hEtSecond[LSP]->Fill(metlsp[1], weight );
	  hEtSecond[nu]->Fill(metnu[1], weight );

	  //MEff = mag(MET) + sumEt
	  if( bucket != 999 )
		hMEff[bucket]->Fill(vecEt.Et() + pthad);
	  if( numhad == 0 )
		hMEff[zjet]->Fill(vecEt.Et() + pthad);
	  hMEff[allBuckets-1]->Fill(vecEt.Et() + pthad, weight);

	  //Ht = mag(MET) + sumEt + ptLep
	  if ( bucket != 999 && TausAreLeptons) 
		hHt[bucket]->Fill(vecEt.Et() + ptlep + pthad + pttau, weight);
	  else if( bucket != 999 )
		hHt[bucket]->Fill(vecEt.Et() + ptlep + pthad, weight);
	  if( numhad == 0 )
		hHt[zjet]->Fill(vecEt.Et() + ptlep + pthad, weight);
	  if (TausAreLeptons) 
		hHt[allBuckets-1]->Fill(vecEt.Et() + ptlep + pthad + pttau, weight);
	  else
		hHt[allBuckets-1]->Fill(vecEt.Et() + ptlep + pthad, weight);

	  //cout << "Done fill hists\n";
    }//end event loop
	if(bad_id != 0) {
	  cout << "\nBad id: " << bad_id << endl;
	}
	else if( ngluons != 0) {
	  //cout << "\nGluons found. ngluons: " << ngluons << endl;
	}
	
  }
  //end file loop

  cout << "\n\n";
  if(numidtau > 0 ) cout << "\nNum events tau first/second: " << numidtau << endl;
  if(numtotaltau > 0 ) cout << "Total num tau: " << numtotaltau << endl;
  if(nzerojet > 0 ) cout << "Num zero jet: " << nzerojet << endl;

  std::cout << std::endl;
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
  
  std::cout << "Prefix: " << prefix << " processed: " << nEventsTotal
	    << " Events, found: " << duplicates_total_n 
	    << " Duplicates and selected: " << nCandidatesSelected
	    << " Candidates." << endl << endl;
  cout << "evt_scale1fb = " << cms2.evt_scale1fb() << endl;

  //see commented code in 'warren_hists.C' above print_mcid_table

  cout << "\nWriting hists, name: " << name << "\n\n";

  //write all hists to root file
  //TFile outfile("lm1_04_hist_eta25_taus.root","RECREATE");
  //TFile outfile("lm1_05-002_eta25.root","RECREATE");
  if( name == "" ) cout << "ERROR: IMPROPER INIT: name=" << name << " \n\n\n";

  TFile outfile(name, "RECREATE");
  //TFile outfile(name, "RECREATE");
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

	hMET[j]->Write();
	hEtFirst[j]->Write();
	hEtSecond[j]->Write();

	hEtmPt[j]->Write();
	hMEff[j]->Write();
	hHt[j]->Write();

//    hist_id[j]->Write();
//    hist_motherid[j]->Write();
//    hist_genid[j]->Write();
//    hist_genmotherid[j]->Write();
  }

  if( nHardlsp > 0 )
	cout << "Number hard scatter lsp: " << nHardlsp << "\n\n";

  //same for all lepton vars in first 10 buckets b'c of how buckets are filled
  print_entries_table(hPtLep, 0, lepBuckets-1);

  cout << "\nAll events:\n";
  //print_hard_type(hard_type, susy_types);//old version
  print_hard_type(hard_type, hard_chgflv, susy_types, weight);
  cout << "\nZero Jet Bin:\n";
  print_hard_type(zjet_hard_type, zjet_hard_chgflv, susy_types, weight);

  //put 'all' bucket of variables for which I want 90% efficiency points
  int numvars = 9;
  TH1F allhists[numvars];
  allhists[0] = *hPtLepFirst[allBuckets-1];
  allhists[1] = *hPtLepSecond[allBuckets-1];
  allhists[2] = *hPtLepThird[allBuckets-1];
  allhists[3] = *hPtLep[allBuckets-1];
  allhists[4] = *hMET[allBuckets-1];
  allhists[5] = *hSumEt[allBuckets-1];
  allhists[6] = *hMEff[allBuckets-1];	  
  allhists[7] = *hHt[allBuckets-1];
  allhists[8] = *hNjets[allBuckets-1];

  //print table of 90% 
  print_90eff(allhists, numvars );

  outfile.Close();

  return name;
}
//end char* ScanChain


void print_entries_table(TH1F* hist[], int start, int stop) {
//one call prints table for that hist which is passed in.
//start and stop are buckets, and inclusive
//WARNING: the way I added the sign/flavor counting,
//  it includes ALL buckets ALWAYS (even if not displayed)

  int sssf=0,ssof=0,ossf=0,osof=0;
  //use global to enforce hist vector size
  if( stop > int(allBuckets)-1 ) { return; }

  const int numbuckets = stop - start +1;//+1 for inclusive
  double entries[numbuckets]; //display purposes
  double total_entries = 0;

  cout << "Cuts: ";
  if( EtaCut ) cout << "Eta < " << EtaCutValue << "  ";
  if( TausAreLeptons ) cout << " including taus \n";
  else cout << " excluding taus \n";

  cout << "Var: " << hist[0]->GetTitle() << "\n\n";
  cout << "bucket: \t";
  for( int i = start; i <= stop; i++) {
	cout << suffix[i] << "\t";
	entries[i-start] = hist[i]->GetEntries();
	total_entries += hist[i]->GetEntries();
	//entries[i-start] = hist[i]->Integral();
  }
  cout << "\nentries:\t";
  for( int i = 0; i < numbuckets; i++) {
	//NOTE: change this to use function 'buck_to_chgflv'
	cout << entries[i] ;
	if(i==0 || i==4 || i==7 || i==9 )
	  sssf += (int)entries[i];
	else if( i==1 || i==8)
	  ossf += (int)entries[i];
	else if( i==2 || i==6)
	  ssof += (int)entries[i];
	else if( i==3 || i==5)
	  osof += (int)entries[i];
	
	//trying to get the tabs right...
	if(entries[i] < 10) cout << "   \t";
	else if(entries[i] < 100) cout << "  \t";
	else if(entries[i] < 1000) cout << " \t";
	else cout << "\t";
  }

  //cout << "\n\n\tsf\tof" << "\nss:\t" << sssf << "\t" << ssof
  //   << "\nos:\t" << ossf << "\t" << osof;

  //cout << "\n\nTotal entries:  " << total_entries << "   Total events:  "
  //   << total_entries*cms2.evt_scale1fb();
  cout << "\n\n";
}
//end void print_entries_table(...)


char* name_from_cuts(char* prefix, char* tag){ 

  char* name = new char[100];
  //note, below here, paste code back in to revert (minus return)
  //strcpy(name, namebase);
  strcpy(name, prefix);
  strcat(name, "_");
  strcat(name, tag);
  //strcat(name, namebase);
  if( EtaCut ) {
	if( EtaCutValue == 2.5 )
	  strcat(name, "_eta25");
	else if( EtaCutValue == 3.0 )
	  strcat(name, "_eta30");
	else
	  strcat(name, "_eta");
  }
  if( TausAreLeptons ) strcat(name, "_tau");
  strcat(name, ".root");

  return name;
}
//end char* name_from_cuts

