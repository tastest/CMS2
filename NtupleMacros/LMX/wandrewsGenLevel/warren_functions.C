using namespace std;
#include <iostream>
//#include <sstream>
//#include <vector>
//#include <set>
//#include "TChain.h"
//#include "TFile.h"
#include "TH1F.h"
//#include "TH2F.h"


//my functions (old hists functions in warren_hists.C)
//definition is at end of this file
//int get_lep_bucket(double id_1, double id_2);
//double* get_id_rank(int par, double id, double pt_rank_lep[], double id_lep[]);
//double* get_pt_rank(int par, double pt_rank_lep[]);
//double* get_Et_rank(int par, double pt_lep[]);
//void print_hard_type(int* Hard_type, const int susy_types );
//int* susy_pair_type(Int_t id1, Int_t id2);
//bool* susy_flavor(Int_t id);


//actualy, keeping first and second separate doesn't help becuase i need to classify each event, not just each particle
//need function: input is 2 mcids of susy pair production
//return an array whose elements are final state:
//0=2gluino:
//1=gluino-squark:
//8=gluino-gaugino;
//2=2squark  ;
//3=squark-gaugino  ;
//4=squark-slepton  ;
//5=2gaugino  ;
//6=slepton-gaugino  ;
//7=2slepton  ;
//9=higgs;

//this function is slightly redundant with below, but this is for single particle, above is for pair
bool* susy_flavor(Int_t id) {

  const Int_t mill = 1000000;
  enum { squark, slepton, gaugino, higgs };
  bool* id_flav = new bool[higgs+1];
  for(int i=0;i<higgs+1;i++) id_flav[i]=false;

  if( (id > mill && id <= mill+6) ||
	  (id > 2*mill && id <= 2*mill+6) )
	id_flav[squark] = true;
  else if( (id >= mill+11 && id <= mill+16) ||
		   (id >= 2*mill+11 && id <= 2*mill+15) )
	id_flav[slepton] = true;
  else if( id >= mill+22 && id <= mill+37 )
	id_flav[gaugino] = true;
  else if( id >=35 && id <= 37 )
	id_flav[higgs] = true;

  if(!id_flav[squark] && !id_flav[slepton] && !id_flav[gaugino] && !id_flav[higgs] && id != mill+21)
	cout << "susy_flavor->Bad id: " << id << endl;
  return id_flav;
}
//end susy_flvaor

//overload this function: this one returns vector
int* susy_pair_type(Int_t id1, Int_t id2, int idlist[]) {
  id1 = abs(id1);
  id2 = abs(id2);
  if ( id1 < 0 || id2 < 0 )
	cout << "bad id: " << id1 << " " << id2 << endl;
  const int susy_types = 10;
  int* type = new int[susy_types];
  for(int i = 0;i < susy_types;i++) {
	type[i] = idlist[i];
  }
  
  const Int_t gluino = 1000021;
  
  enum { squark, slepton, gaugino, higgs };
  bool* id1_flav;// = new bool[gaugino+1];//not necessary
  bool* id2_flav;// = new bool[gaugino+1];
  
  id1_flav = susy_flavor(id1);
  id2_flav = susy_flavor(id2);

  if(id1==gluino && id2==gluino )
	type[0]++;
  else if( (id1 == gluino && id2_flav[squark]) ||
		   (id2 == gluino && id1_flav[squark]) )
	type[1]++;
  else if( id1_flav[squark] && id2_flav[squark] )
	type[2]++;
  else if( (id1_flav[squark] && id2_flav[gaugino]) ||
		   (id2_flav[squark] && id1_flav[gaugino]))
	type[3]++;
  else if( (id1_flav[squark] && id2_flav[slepton]) ||
		   (id2_flav[squark] && id1_flav[slepton]))
  	type[4]++;
  else if(id1_flav[gaugino] && id2_flav[gaugino] )
	type[5]++;
  else if((id1_flav[slepton] && id2_flav[gaugino]) ||
		  (id2_flav[slepton] && id1_flav[gaugino]))
	type[6]++;
  else if(id1_flav[slepton] && id2_flav[slepton] )
	type[7]++;
  else if( (id1 == gluino && id2_flav[gaugino]) ||
		   (id2 == gluino && id1_flav[gaugino]) )
	type[8]++;
  else if(id1_flav[higgs] && id2_flav[higgs] )
	type[9]++;
  else
	cout << "susy_pair_type not found: " << id1 << " : " << id2 << endl;

  return type;
}
//end susy_pair_type(Int_t, Int_t, int [])

//this version returns the index of the susy particle
int susy_pair_type(Int_t id1, Int_t id2) {

  id1 = abs(id1);
  id2 = abs(id2);
  if ( id1 < 0 || id2 < 0 )
	cout << "bad id: " << id1 << " " << id2 << endl;
  
  const Int_t gluino = 1000021;
  
  enum { squark, slepton, gaugino, higgs };
  bool* id1_flav;// = new bool[gaugino+1];
  bool* id2_flav;// = new bool[gaugino+1];
  
  id1_flav = susy_flavor(id1);
  id2_flav = susy_flavor(id2);

  if(id1==gluino && id2==gluino )
	return 0;
  else if( (id1 == gluino && id2_flav[squark]) ||
		   (id2 == gluino && id1_flav[squark]) )
	return 1;
  else if( id1_flav[squark] && id2_flav[squark] )
	return 2;
  else if( (id1_flav[squark] && id2_flav[gaugino]) ||
		   (id2_flav[squark] && id1_flav[gaugino]))
	return 3;
  else if( (id1_flav[squark] && id2_flav[slepton]) ||
		   (id2_flav[squark] && id1_flav[slepton]))
  	return 4;
  else if(id1_flav[gaugino] && id2_flav[gaugino] )
	return 5;
  else if((id1_flav[slepton] && id2_flav[gaugino]) ||
		  (id2_flav[slepton] && id1_flav[gaugino]))
	return 6;
  else if(id1_flav[slepton] && id2_flav[slepton] )
	return 7;
  else if( (id1 == gluino && id2_flav[gaugino]) ||
		   (id2 == gluino && id1_flav[gaugino]) )
	return 8;
  else if(id1_flav[higgs] && id2_flav[higgs] )
	return 9;
  else
	cout << "susy_pair_type not found: " << id1 << " : " << id2 << endl;

  return 999;
}
//end susy_pair_type(Int_t,Int_t)


void print_90eff(TH1F hist[], int vars) {

  double total;
  int j;
  //double tenpv[vars];//ten percent values

  cout << "\nTen Percent values:\nVariable\t\t\tValue\n";
  for( int i = 0; i< vars; i++ ) {
	j=0;
	total = hist[i].GetBinContent(0);
	while( total < 0.1*(hist[i].Integral()) ) {
	  j++;
	  total += hist[i].GetBinContent(j);
	}
	//tenpv[i] = hist[i].GetBinLowEdge(j+1);
	cout << hist[i].GetTitle() << "\t\t\t" << hist[i].GetBinLowEdge(j+1) << endl;
  }
  //cout << "\nValue:\t\t";
  //for(int k=0;k<vars;k++) cout << tenpv[k] << "\t\t";
  cout << endl;
}
//end void print_90eff(TH1F* hist[])


//this one is for just total column
void print_hard_type(int hard_type[], const int susy_types) {

  char* typ_names[19] = {"2gluino        ",
						 "gluino-squark  ",
						 "2squark        ",
						 "squark-gaugino ",
						 "squark-slepton ",
						 "2gaugino       ",
						 "slepton-gaugino",
						 "2slepton       ",
						 "gluino-gaugino ",
						 "2higgsino      "};
  
  for(int i=0;i<susy_types;i++) {
	std::cout << typ_names[i] << " : " << hard_type[i] << std::endl;
  }

}
//end void print_hard_type(int [], int)

//this one prints totals, plus separate by charge,flavor of final state
void print_hard_type(double hard_type[], double hard_chgflv[][10], const int susy_types, double weight) {

  double totals[5];//totals, for all, and for 4 charg/flav categories
  for(int i=0;i<5;i++) { totals[i] = 0; }

  char* typ_names[19] = {"2-gluino       ",
						 "gluino-squark  ",
						 "2-squark       ",
						 "squark-gaugino ",
						 "squark-slepton ",
						 "2-gaugino      ",
						 "slepton-gaugino",
						 "2-slepton      ",
						 "gluino-gaugino ",
						 "2-higgsino     ",
						 "Totals         "};
  
  cout << "S-Particles\tTotal\t\tsssf\t\tssof\t\tossf\t\tosof\n";
  for(int i=0;i<susy_types;i++) {
	cout << typ_names[i] << "  " << hard_type[i]*weight << "\t\t";
	totals[0] += hard_type[i]*weight;
	//if( hard_type[i] == 0 ) cout << "\t\t";
	for(int j=0;j<4;j++) {
	  cout << hard_chgflv[j][i]*weight << "\t\t";
	  totals[j+1] += hard_chgflv[j][i]*weight; //j+1 bc 0 is total
	  //if( hard_chgflv[i] == 0 ) cout << "\t";
	}
	cout << endl;
  }

  cout << typ_names[susy_types] << "  " ;
  for(int i=0;i<5;i++) { cout << totals[i] << "\t\t" ; }
  cout << "\n\n";
}
//end void print_hard_type(int [], int [][], int)


int buck_to_chgflv(int bucket) {

  if(bucket==0 || bucket==4 || bucket==7 || bucket==9 )
	return 0;
	//sssf += (int)entries[i];
  else if( bucket==1 || bucket==8)
	return 2;
	//ossf += (int)entries[i];
  else if( bucket==2 || bucket==6)
	return 1;
	//ssof += (int)entries[i];
  else if( bucket==3 || bucket==5)
	return 3;
  //osof += (int)entries[i];

  return 999;
}
//end int buck_to_chgflv


double* get_pt_rank(double pt_lep[], double pt) {
  //cout << "pt: "; for(int i=0;i<3;i++) { cout << pt_lep[i] << " "; }
  //cout << endl;
  double* pt_rank_lep = new double[3];
  for(int i=0;i<3;i++) { pt_rank_lep[i] = pt_lep[i]; }
  
  //if ( cms2.genps_p4()[par].pt() > pt_lep[0] ) {
  if ( pt > pt_lep[0] ) {
	pt_rank_lep[2] = pt_lep[1];
	pt_rank_lep[1] = pt_lep[0];
	//pt_rank_lep[0] = cms2.genps_p4()[par].pt();
	pt_rank_lep[0] = pt;
  }
  else if ( pt > pt_lep[1] ) {
	pt_rank_lep[2] = pt_lep[1];
	pt_rank_lep[1] = pt;
  }
  else if ( pt > pt_lep[2] ) {
	pt_rank_lep[2] = pt;
  }
  return pt_rank_lep;
  //return;
}
//end double* get_pt_rank

//this one is for lsp--use et, only two
double* get_Et_rank(double pt_lep[], double Et) {
  double* pt_rank_lep = new double[2];
  for(int i=0;i<2;i++) { pt_rank_lep[i] = pt_lep[i]; }
  
  //if ( cms2.genps_p4()[par].Et() > pt_lep[0] ) {
  if ( Et > pt_lep[0] ) {
	pt_rank_lep[1] = pt_lep[0];
	//pt_rank_lep[0] = cms2.genps_p4()[par].Et();
	pt_rank_lep[0] = Et;
  }
  else if ( Et > pt_lep[1] ) {
	pt_rank_lep[1] = Et;
  }
  return pt_rank_lep;
}


double* get_id_rank(double id_new, double pt_rank_lep[], double id[], double pt) {
  //cout << "id: "; for(int i=0;i<3;i++) { cout << id[i] << " "; }
  //cout << "pt: "; for(int i=0;i<3;i++) { cout << pt_rank_lep[i] << " "; }
//  cout << ". new id: " << id << " new pt: " << cms2.genps_p4()[par].pt() << endl;

  double* id_lep = new double[3];
  for(int i=0;i<3;i++){ id_lep[i] = id[i]; }

  //if ( cms2.genps_p4()[par].pt() > pt_rank_lep[0] ) {
  if ( pt > pt_rank_lep[0] ) {
	id_lep[2] = id[1];
	id_lep[1] = id[0];
	id_lep[0] = id_new;
  }
  else if ( pt > pt_rank_lep[1] ) {
	id_lep[2] = id[1];
	id_lep[1] = id_new;
  }
  else if ( pt > pt_rank_lep[2] ) {
	id_lep[2] = id_new;
  }
  return id_lep;
  //return;
}
//end double* get_id_rank


int get_lep_bucket(double id_1, double id_2) {
  int bucket = 999;
  //these '*= -1' flip sign so that e- = -11, just for my sanity
  id_1 *= -1; id_2 *= -1;
  
  if( id_1 == 13 && id_2 == 13 ) {
    bucket = 0;
  }
  else if( (id_1 == 13 && id_2 == -13)
	   || (id_1 == -13 && id_2 == 13) ) {
    bucket = 1;
  }
  else if( (id_1 == 13 && id_2 == 11)
	   || (id_1 == 11 && id_2 == 13) ) {
    bucket = 2;
  }
  else if( (id_1 == 13 && id_2 == -11)
	   || (id_1 == -11 && id_2 == 13) ) {
    bucket = 3;
  }
  else if( id_1 == -13 && id_2 == -13 ) {
    bucket = 4;
  }
  else if( (id_1 == -13 && id_2 == 11)
	   || (id_1 == 11 && id_2 == -13) ) {
    bucket = 5;
  }
  else if( (id_1 == -13 && id_2 == -11)
	   || (id_1 == -11 && id_2 == -13) ) {
    bucket = 6;
  }
  else if( id_1 == 11 && id_2 == 11 ) {
    bucket = 7;
  }  
  else if( (id_1 == -11 && id_2 == 11)
	   || (id_1 == 11 && id_2 == -11) ) {
    bucket = 8;
  }
  else if( id_1 == -11 && id_2 == -11 ) {
    bucket = 9;
  }
  
  return bucket;
}
//end int getbucket      



//////////////////////////////////////
//this was at the end of the hist declaration in ScanChain

  //need to replace this ridiculous number of bins with 2 hists--one for SM particles, one for SUSY--
  //make sure that susy mcid work here
  //ie, if 1e6 <= mcid < 2e6, plot mcid/1e6
  //  if mcid >= 2e6, plot (on same hist) mcid/2e6 + 1000 (if above line never goes over 1000, or even 100 if that works)
  //make 2 hist vectors for id: one for SM, one for susy

  /*
  //  const Int_t mc_bins = 3e6;
  Int_t gen_bins = 1e7;
  int fakes[allBuckets];

  TH1F* hist_id[allBuckets];
  TH1F* hist_motherid[allBuckets];
  TH1F* hist_genid[allBuckets];
  TH1F* hist_genmotherid[allBuckets];

  //my old hist init
  */

///////////////////////////////
//this used to be right before the file loop in ScanChain
  //CONSTANTS
  //const float zmass = 91.19; //just making sure that I use the same Z mass everywhere!

  // CUTS
  /*const int   goodLeptonsCut        = 3;     // good lepton cut
  const float triggerLeptonMinPtCut = 20.;   // one of the leptons of the trilepton canidate has to have at least this pt
  const float leptonMinPtCut        = 10.;   // all leptons of the trilepton canidate have to have at least this pt
  const float electronMETCut        = 40.;   // cut on MET if lepton not belonging to primary Z is an electron
  const float muonMETCut            = 20.;   // cut on MET if lepton not belonging to primary Z is an muon
  bool        useMETAll             = true; // use metAll corrected for all muons instead of met corrected for muons from cand.
  */
  
