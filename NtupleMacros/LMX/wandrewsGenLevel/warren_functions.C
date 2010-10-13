#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "TH1F.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/VectorUtil.h"
using namespace std;

//code for susy_pair_type function
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

//this function is for single particle, above is for pair
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
  else if( (id >=35 && id <= 37) || id == 25 )//25=SM higgs produced with other
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


char* round_char(double in){
  if( in < 1 ) 
	return Form("%4.2f",in);
  else if ( in < 10 )
	return Form("%3.1f",in);
  else { //in > 1
	if( in - (int)in >= 0.5 )
	  return Form("%i", int(in)+1);
	else
	  return Form("%i", int(in));
  }
}
//end round_char(double)


int get_90_bin(TH1F hist ) {
  int j=0;
  double total = hist.GetBinContent(0);
  //arguments of Integral() includes overflow bin b'c '+1'
  while( total < 0.1*( hist.Integral(0, hist.GetNbinsX()+1) ) ) {
	j++;
	total += hist.GetBinContent(j);
  }
  return j;
}
//end get_90_bin

TString* init_var_desc(const int vars) {
  TString* desc = new TString[vars];
  if( vars != 12 ) { //dummy proof b'c root won't tell me idx out of bounds
	cout << "\nERROR: BAD INPUT TO init_var_desc. see warren_functions.C\n";
	return desc;
  }

  desc[0] = "Pt of leading lepton";
  desc[1] = "Pt of second lepton";
  desc[2] = "Pt of third lepton";
  desc[3] = "Scalar sum of all lepton pts";
  desc[4] = "Mass of two leading leptons";
  desc[5] = "Et of LSP + neutrinos";
  desc[6] = "Scalar sum of quark, gulon pts";
  desc[7] = "!SumEt + MET"; //twiki comment !
  desc[8] = "MEff + lepton pt";
  desc[9] = "Number jets";
  desc[10] = "Number leptons";
  desc[11] = "Etj(closest in phi to MET)/Mjj";

  return desc;
}
//end init_var_desc

TString* init_bvar_desc(const int vars) {
  TString* desc = new TString[vars];
  if( vars != 5 ) { //dummy proof b'c root won't tell me idx out of bounds
	cout << "\nERROR: BAD INPUT TO init_bvar_desc\n";
	return desc;
  }

  desc[0] = "Pt of leading b quark";
  desc[1] = "Pt of leading b quark";
  desc[2] = "Pt of leading b quark";
  desc[3] = "Scalar sum of all b Pt";
  desc[4] = "Number b-jets";

 return desc;
}
//end init_bvar_desc

TString get_path(TString title) {
  TString path = "http://uaf-2.t2.ucsd.edu/~wandrews/SUSY/";
  TString prefix = title;
  //TString col = title; //was: (TString)hist1[0].GetTitle();
  int idx = title.Index("_",0); //find first underscore
  //prefix = col;
  prefix.Resize( idx );
  path += prefix + "/";
  return path;
}
//end get_path(Tstring)

//single column version: using for b-quark table
const char* print_90eff(TH1F hist[], int vars) {

  ostringstream stream;
  stream.setf(ios::fixed);
  int j;
  //const unsigned int suffixlen = 1; //suffix = 'b'
  TString path = get_path((TString)hist[0].GetTitle());
  TString* desc = init_bvar_desc(vars);

  stream << "\n|*Variable*|*Value*|*Variable Definition*|\n";
  for( int i = 0; i< vars; i++ ) {
	j = get_90_bin( hist[i] );
	//total = hist[i].GetBinContent(0);
	//while( total < 0.1*(hist[i].Integral()) ) {
	//  j++;
	//  total += hist[i].GetBinContent(j);
	//}
	TString col= (TString)hist[i].GetTitle();
	//col.resize( strlen(hist[i].GetTitle()) - suffixlen);//names: ..._b_stack

	stream << "|" << col << "|" ;
	stream << "<a href=\"" << path << col << "_stack.png\">";
	stream << round_char(hist[i].GetBinLowEdge(j+1)) << "</a>";
	/*if(i==0) { //getting rid of overlay altogether b'c taus mess it up
	  stream << " <a href=\"" << path << hist[0].GetTitle() << ".png\">";
	  stream << "over</a>";
	  }*/

	stream << "|" << desc[i] << "|\n";

  }
  stream << endl;
  return stream.str().c_str();
}
//end void print_90eff(TH1F* hist[])

//version for two vectors, display in two columns
const char* print_90eff(TH1F hist1[], TH1F hist2[], int vars) {
  ostringstream stream;
  stream.setf(ios::fixed);
  int j1=0, j2=0;
  const unsigned int suffixlen = 4; //suffix = 'all'
  TString* desc;
  if( vars > 10 ) // use vars to distinguish between all vars, b vars
	desc = init_var_desc(vars);
  else if( vars <= 10 )
	desc = init_bvar_desc(vars);
  TString path = get_path((TString)hist1[0].GetTitle());

  //stream << "\nKinematic plots with 90\% Efficiency Points:";
  stream << "\n|*Variable*|* All *|";
  if( vars > 10 )
	stream << "*zero-jet*|"; //again use vars for zjet(all vars), tau(b)
  else
	stream << "*tau*|";
  stream << "*Variable Definition*|\n";
  for( int i = 0; i< vars; i++ ) {
	j1 = get_90_bin( hist1[i] );
	j2 = get_90_bin( hist2[i] );
	TString col = (TString)hist1[i].GetTitle();
	col.Resize( strlen(hist1[i].GetTitle())-suffixlen ); //suffixlen = 4

	stream << "|" << col << "|";
	stream << "<a href=\"" << path << col << "_stack.png\">";
	stream << round_char(hist1[i].GetBinLowEdge(j1+1)) << "</a>"; 
	/*if(i==0) {//getting rid of overlay altogether b'c taus mess it up
	  stream << " <a href=\"" << path << col << "_all.png\">";
	  stream << "over</a>";
	  }*/

	//***second column, now same as first***
	TString nm2 = (TString)hist2[i].GetTitle();
	nm2.Resize( strlen(hist2[i].GetTitle())-4 );
	stream << "|<a href=\"" << path << nm2 << "_stacktau.png\"> ";
	stream << round_char(hist2[i].GetBinLowEdge(j2+1)) << "</a> "; 
	/*if( i==0 ) { 
	  stream << "<a href=\"" << path << hist2[i].GetTitle() <<".png\">";
	  stream << "over</a>";
	  }*/
	stream << "|" << desc[i] << "|\n";

  }//end for vars
  stream << endl;
  return stream.str().c_str();
}
//end void print_90eff(TH1F hist[], TH1F hist2[])

//version for four vectors, display in four columns, separate taus
const char* print_90eff(TH1F hist1[], TH1F hist2[], TH1F hist3[], TH1F hist4[], int vars) {
  //cout << "\nstarting 90 eff\n";
  ostringstream stream;
  stream.setf(ios::fixed);
  int j1=0, j2=0, j3=0, j4=0;
  const unsigned int suffixlen = 4; //suffix = 'all'
  TString* desc = init_var_desc(vars);
  TString path = get_path((TString)hist1[0].GetTitle());
  
  //stream << "\nKinematic plots with 90\% Efficiency Points:";
  stream << "\n|*Variable*|*All dilep*|*All tau*|*zero-jet<br>dilep*|*zero-jet<br>tau*|*Variable Definition*|\n";
  for( int i = 0; i< vars; i++ ) {
	j1 = get_90_bin( hist1[i] );
	j2 = get_90_bin( hist2[i] );
	j3 = get_90_bin( hist3[i] );
	j4 = get_90_bin( hist4[i] );
	TString col = (TString)hist1[i].GetTitle();
	col.Resize( strlen(hist1[i].GetTitle())-suffixlen ); //suffixlen = 4

	//first column: all dilep
	stream << "|" << col;
	stream << "|<a href=\"" << path << col << "_stack.png\">";
	stream << round_char(hist1[i].GetBinLowEdge(j1+1)) << "</a>"; 
	/*if(i==0) { //getting rid of overlay altogether b'c taus mess it up
	  stream << " <a href=\"" << path << col << "_all.png\">";
	  stream << "over</a>";
	  }*/

	//second column: all tau
	stream << "|<a href=\"" << path << col << "_stacktau.png\">";
	stream << round_char(hist2[i].GetBinLowEdge(j2+1)) << "</a>"; 

	//third column, zjet dilep
	TString nm2 = (TString)hist3[i].GetTitle();
	nm2.Resize( strlen(hist3[i].GetTitle())-4 );
	stream << "|<a href=\"" << path << nm2 << "_stack.png\"> ";
	stream << round_char(hist3[i].GetBinLowEdge(j3+1)) << "</a> "; 
	/*if( i==0 ) { //getting rid of overlay altogether b'c taus mess it up
	  stream << "<a href=\"" << path << nm2 <<".png\">";
	  stream << "over</a>";
	  }*/
	
	//forth column, zjet tau
	stream << "|<a href=\"" << path << nm2 << "_stacktau.png\"> ";
	stream << round_char(hist4[i].GetBinLowEdge(j4+1)) << "</a> "; 	

	stream << "|" << desc[i] << "|\n";

  }
  stream << endl;
  return stream.str().c_str();
}
//end void print_90eff(TH1F hist[], TH1F hist2[], TH1F hist3, TH1F hist4)


TH1F* bookChgflv(TH1F* hist, int idx) {
  TString namet = hist->GetTitle();
  
  int titlelen = strlen(hist->GetTitle());
  if( namet(titlelen-3, titlelen) == "all" )
	namet.Resize( titlelen - 3 ); //remove "all"
  else
	namet += "_"; //special for zero-jet
  
  if( idx == 0 )
	namet += "sssf";
  else if( idx == 1 )
	namet += "ssof";
  else if( idx == 2 )
	namet += "ossf";
  else if( idx == 3 )
	namet += "osof";
  else if( idx == 4 )
	namet += "sst_"; //tau
  else if( idx == 5 )
	namet += "ost_";
  
  TH1F* histout = new TH1F(namet, namet, hist->GetNbinsX(), 0, hist->GetXaxis()->GetXmax() );
  histout->Sumw2();

  return histout;
}
//end bookChglfv(TH1F*, int)

//new plan for z-jet chgflv: don't use this function which relies on buckets
//instead, fill the chgflv plots during filling of rest of histos in looper
TH1F* chgflv_plots(TH1F* hvar[], int buckets, int idx) {
  //buckets = allBuckets in warren_looper.C
  //idx is index in loop (from which this fn is called) over chgflv of returned hist

  TH1F* hist = bookChgflv( hvar[buckets-1], idx );
  
  if(idx == 0) { //sssf
	hist->Add( hvar[0], hvar[4] );
	hist->Add( hvar[7] );
	hist->Add( hvar[9] );
  }
  else if( idx == 1 ) { //ssof
	hist->Add( hvar[2], hvar[6] );
  }
  else if( idx == 2 ) { //ossf
	hist->Add( hvar[1], hvar[8] );
  }
  else if( idx == 3 ) { //osof
	hist->Add( hvar[3], hvar[5] );
  }
  else if( idx == 4 ) { //sst
	hist->Add( hvar[10] ); 
  }
  else if( idx == 5 ) { //ost
	hist->Add( hvar[11] );
  }

  return hist;
}
//end chgflv_plots(TH1F[],int,int)

//original version
THStack* make_stack(TH1F* hists[], const unsigned int num) {
  TString name = (TString)hists[0]->GetTitle();
  int start = strlen(hists[0]->GetTitle())-4;

  //cout << "\n suffix: " << suffix << endl;
  name.Resize( start ); //remove "xsxf"
  TH1F* htotal = new TH1F(name+"all_",name+"all_", hists[0]->GetNbinsX(),
						  0, hists[0]->GetXaxis()->GetXmax() );//for stats
  htotal->Sumw2();
  name += "stack";

  THStack *hstack = new THStack( name, name ); //
  TLegend* leg = new TLegend(0.85, 0.65, 0.95, 0.95);
  leg->SetFillColor(0);

  double integrals[num];//these two arrays are for sorting
  bool added[num];
  for( unsigned int i=0; i < num; i++ ) {
	added[i]=false;
	integrals[i] = hists[i]->Integral();	
	hists[i]->SetFillColor(i+2);
	hists[i]->SetLineColor(i+2);
	hists[i]->SetMarkerColor(i+2);
	htotal->Add( hists[i] );
  }
  sort(integrals, integrals + num);

  for( unsigned int i=0; i < num; i++ ){
	for( unsigned int j=0; j < num; j++ ) {
	  if( hists[j]->Integral() == integrals[i] && !added[j] ) {
		added[j] = true;
		hstack->Add(hists[j]);//
		//legend stuff:
		TString iname = (TString)hists[j]->GetTitle();
		TString suffix;
		if( j > 3 )  //for ost_, sst_ -> sst,ost
		  suffix = iname(start, 3 );//second argument is length
		else
		  suffix = iname(start, 4 );
		leg->AddEntry(hists[j],suffix);
	  }
	}
  }
  
  TCanvas c;
  //NOTE: the "hist" option seems to make the optstat box disappear
  //htotal->Draw(); //overwritten 
  hstack->Draw("hist"); //
  //htotal->GetPainter()->PaintStat(1111,0); //nothing at all???
  //htotal->GetPainter()->Draw(); //draws whole hist
  c.SaveAs(name+".png");
  return hstack;
}
//end make_stack(th1f[], int)

//new version: last argument is for dilep vs tau
THStack* make_stack(TH1F* hists[], const unsigned int num, int idx) {
  TString name = (TString)hists[0]->GetTitle();
  int start = strlen(hists[0]->GetTitle())-4;
  name.Resize( start ); //remove "xsxf"
  //TH1F* htotal = new TH1F(name+"all_",name+"all_", hists[0]->GetNbinsX(), 0, hists[0]->GetXaxis()->GetXmax() );//for stats
  //htotal->Sumw2();
  name += "stack";
  if( idx == 1) name += "tau";
  THStack *hstack = new THStack( name, name ); 
  TLegend* leg = new TLegend(0.85, 0.65, 0.95, 0.95);
  leg->SetFillColor(0);

  double integrals[num];//these two arrays are for sorting
  bool added[num];
  for( unsigned int i=0; i < num; i++ ) {
	added[i]=false;
	integrals[i] = hists[i]->Integral();	
	hists[i]->SetFillColor(i+2);
	hists[i]->SetLineColor(i+2);
	hists[i]->SetMarkerColor(i+2);
	//htotal->Add( hists[i] );
  }
  sort(integrals, integrals + num);

  for( unsigned int i=0; i < num; i++ ){
	for( unsigned int j=0; j < num; j++ ) {
	  if( hists[j]->Integral() == integrals[i] && !added[j] ) {
		added[j] = true;
		//legend stuff:
		TString iname = (TString)hists[j]->GetTitle();
		TString suffix;
		if( j > 3 )  //for ost_, sst_ -> sst,ost
		  suffix = iname(start, 3 );//second argument is length
		else
		  suffix = iname(start, 4 );
		if( (idx == 1 && j > 3) || (idx == 0 && j < 4) ) {
		  hstack->Add(hists[j]);
		  leg->AddEntry(hists[j],suffix);
		}
	  }
	}
  }
  
  TCanvas c;
  hstack->Draw("hist"); 
  leg->Draw();
  c.SaveAs(name+".png");
  return hstack;
}
//end make_stack(th1f[], int, int)

//code removed from revised make_stack--reverting
//THStack *hstack[2]; //one dilep, one tau
//  hstack[0] = new THStack( name, name );
//  TString tname = name + "tau";
//  hstack[1] = new THStack( tname, tname );
//hstack[0]->Draw("hist");
//  leg->Draw();
//
//  hstack[1]->Draw("hist");
//  leg->Draw();
//  c.SaveAs(tname+".png");
//		if( j > 3 )
//		  hstack[1]->Add(hists[j]);
//		else
//		  hstack[0]->Add(hists[j]);


//dave's stacking function:
//THStack *getStack(TFile &file, TString nJets, TString hyp_type) {
//  THStack *st_temp = new THStack("st_temp", "");
//  TString histNameSuffix = "_mll_" + nJets + "-N-1_" + hyp_type;
//  std::cout << histNameSuffix.Data() << std::endl;
//  for (Int_t i = 8; i >= 0; --i)
//	{
//	  TString temp = sources[i] + histNameSuffix;
//	  std::cout << "getting " << temp.Data() << std::endl;
//	  TH1F *h1_temp = ((TH1F*)(file.Get(sources[i] + histNameSuffix)->Clone()))->Rebin(5);
//	  h1_temp->Scale(lumiNorm_);
//	  st_temp->Add(h1_temp);
//	  //st_temp->Add(((TH1F*)(file.Get(sources[i] + histNameSuffix)->Clone()))->Rebin(5));
//                         
//	}
//  return st_temp;
//}
//end dave's stacking function


//this fn adds the dilep chgflvs for input into 90_eff--needed to separate dilep from tau
TH1F dilep_sum(TH1F* hists[]) {
  TString name = (TString)hists[0]->GetTitle();
  //int start = strlen(hists[0]->GetTitle())-4;
  name.Resize( strlen(hists[0]->GetTitle())-4 ); //remove "xsxf"
  name += "sum"; //for print_90eff function
  
  TH1F hout(name, name, hists[0]->GetNbinsX(), 0, hists[0]->GetXaxis()->GetXmax() );
  hout.Sumw2();

  hout.Add(hists[0],hists[1]);
  hout.Add(hists[2]);
  hout.Add(hists[3]);//4 chgflv categories
  return hout;
}
//end TH1F dilep_sum(TH1F[])

TH1F tau_sum(TH1F* hists[]) { //repeat above for taus
  TString name = (TString)hists[0]->GetTitle();
  //int start = strlen(hists[0]->GetTitle())-4;
  name.Resize( strlen(hists[0]->GetTitle())-4 ); //remove "xsxf"
  //name += "all"; //for print_90eff function--may cause name duplicate problem
  name += "tsm"; //changed...
  //cout << name << endl;
  
  TH1F hout(name, name, hists[0]->GetNbinsX(), 0, hists[0]->GetXaxis()->GetXmax() );
  hout.Sumw2();

  hout.Add(hists[4],hists[5]);//two tau categories
  return hout;
}
//end TH1F tau_sum(TH1F[])

void save_overlay(TH1F* h1, TH1F* h2, TH1F* h3) {

  h1->SetLineColor( 2 ); //red
  h2->SetLineColor( 3 ); //green
  h3->SetLineColor( 4 ); //blue
  //set rangey of h1 based on max of h2
  if( h2->GetMaximum() > h1->GetMaximum() )
	h1->GetYaxis()->SetRangeUser(0, h2->GetMaximum()*1.1 );

  TCanvas c;
  h1->Draw();
  h2->Draw( "same" );
  h3->Draw( "same" );
  c.SaveAs( (TString)h1->GetTitle() + ".png" );
}
//end save_overlay(TH1F*, TH1F*, TH1F* )


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
const char* print_hard_type(double hard_type[], double hard_chgflv[][10], const int susy_types, double weight, unsigned int numchgflv) {

  ostringstream stream;
  stream.setf(ios::fixed);
  double coltotal[numchgflv+1];//totals, for all, and for charg/flav categories
  double diltotal;//[numchgflv+1];
  double tautotal;//[numchgflv+1];
  //double 
  for(unsigned int i=0;i<numchgflv+1;i++) {
	coltotal[i] = 0; 
	//diltotal[i] = 0; 
	//tautotal[i] = 0;
  }

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
  
  stream << "|*S-Particles*|*Total*|*sssf*|*ssof*|*ossf*|*osof*|";
  if( numchgflv == 6 )
	stream << "*dilep tot*|*sst*|*ost*|*t tot*|";
  stream << "\n";
	  
  //char* out1 = new char[50];
  
  for(int i=0;i<susy_types;i++) {
	if( hard_type[i] == 0 ) continue; //don't print empty rows
	stream << "|" << typ_names[i] << "|*"
		   << round_char(hard_type[i]*weight) << "*|";
	coltotal[0] += hard_type[i]*weight;
	diltotal = 0; tautotal = 0;
	for(unsigned int j=0;j<numchgflv;j++) {
	  stream << round_char(hard_chgflv[j][i]*weight) << "|";
	  coltotal[j+1] += hard_chgflv[j][i]*weight; //j+1 bc 0 is total
	  if( j == 3 && numchgflv > 4 ) { //print total only if taus after
		diltotal += hard_chgflv[j][i]*weight; //j+1 bc 0 is total
		stream << "*" << round_char(diltotal) << "*|";
	  }
	  else if( j < 4 ) //dilep
		diltotal += hard_chgflv[j][i]*weight; //j+1 bc 0 is total
	  else if( j == 5 ) {
		tautotal += hard_chgflv[j][i]*weight; //j+1 bc 0 is total
		stream << "*" << round_char(tautotal) << "*|";
	  }
	  else if( j >= 4 ) //tau
		tautotal += hard_chgflv[j][i]*weight; //j+1 bc 0 is total
	}
	stream << endl;
  }

  diltotal = 0; tautotal = 0;
  stream << "|*" << typ_names[susy_types] << "*|" ;
  for(unsigned int i=0;i<numchgflv+1;i++) {
	if( i == 0 ) stream << "*";
	stream << round_char(coltotal[i]);
	if( i == 0 ) stream << "*";
	stream << "|" ;
	if( i == 4 && numchgflv > 4 ) {
	  diltotal += coltotal[i];
	  stream << "*" << round_char(diltotal) << "*|";
	}
	else if( i < 4 && i != 0 && numchgflv > 4) //dilep
	  diltotal += coltotal[i];
	else if( i == 6 ) {
	  tautotal += coltotal[i];
	  stream << "*" << round_char(tautotal) << "*|";
	}
	else if( i >= 5 ) //tau
	  tautotal += coltotal[i];
  }
  stream << "\n\n";

  return stream.str().c_str();
}
//end void print_hard_type(double [], double [][], int, double)


double get_alpha(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > had1,
			   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > had2,
			   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vecEt)
{
  //built in root function for deltaphi:
  double dif1 = ROOT::Math::VectorUtil::DeltaPhi(vecEt, had1);
  double dif2 = ROOT::Math::VectorUtil::DeltaPhi(vecEt, had2);

  //cout << "Et:" << vecEt.phi() << " h1:" << had1.phi()
  //	   << " h2:" << had2.phi() << endl << " dif1:" << odif1 << "  " << dif1
  //	   << " dif2:" << odif2 << "  " << dif2 << endl;
  if( abs(dif1) < abs(dif2) )  
	return had1.Et()/(had1 + had2).M();
  else
	return had2.Et()/(had1 + had2).M();
}
//end double get_alpha(...)

//alphat is transverse-alpha: instead of mass, use transverse mass
double get_alphat(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > had1,
			   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > had2,
			   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vecEt)
{
  //built in root function for deltaphi:
  double dif1 = ROOT::Math::VectorUtil::DeltaPhi(vecEt, had1);
  double dif2 = ROOT::Math::VectorUtil::DeltaPhi(vecEt, had2);

  //cout << "Et:" << vecEt.phi() << " h1:" << had1.phi()
  //	   << " h2:" << had2.phi() << endl << " dif1:" << odif1 << "  " << dif1
  //	   << " dif2:" << odif2 << "  " << dif2 << endl;
  if( abs(dif1) < abs(dif2) )  
	return had1.Et()/(had1 + had2).Mt();
  else
	return had2.Et()/(had1 + had2).Mt();
}
//end double get_alpha(...)


double* pt_sec_thr(double* pt_rank_lep, double* pt_rank_tau, bool TausAreLeptons){
  //first lepton is pt_rank_lep[0]
  //second is greater of pt_rank_lep[1] and pt_rank_tau[0]
  //third is next greatest

  double* pt = new double[2];

  if( TausAreLeptons ) {
	if( pt_rank_tau[0] > pt_rank_lep[1] ) {
	  pt[0] = pt_rank_tau[0];
	  if( pt_rank_lep[1] > pt_rank_tau[1] ) 
		pt[1] = pt_rank_lep[1];
	  else
		pt[1] = pt_rank_tau[1];
	}
	else {
	  pt[0] = pt_rank_lep[1];
	  if( pt_rank_lep[2] > pt_rank_tau[0] ) 
		pt[1] = pt_rank_lep[2];
	  else
		pt[1] = pt_rank_tau[0];
	}
  }
  else {
	pt[0] = pt_rank_lep[1];
	pt[1] = pt_rank_lep[2];
  }
  
  return pt;
}
//end pt_sec_thr(double*, double*)

//original version: requires only e,mu buckets
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
//end int buck_to_chgflv(int)

int buck_to_chgflv(double idlep, double idtau) {

  if( (idlep < 0 && idtau < 0) || (idlep > 0 && idtau > 0) )
	return 4; //sst
  else 
	return 5; //ost
  
}
//end int buck_to_chgflv(double, double)


int buck_to_sign(int bucket) {
  //pp
  if(bucket==0 || bucket==2 || bucket==7 )
	return 0;
  //mm
  else if( bucket==4 || bucket==6 || bucket==9 )
	return 1;
  //pm (=mp)
  else if( bucket==1 || bucket==3 || bucket==5 || bucket==8 )
	return 2;
  else if( bucket==10 || bucket == 11 )
	return 3;//tau--can't tell pp vs mm, but know ss vs os...

  return 999;
}
//end int buck_to_sign


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

double* count_particles(double count[], int numlep, int numhad,
						int npartcat, double weight){
  // npartcat = 6 currently
  // index: 0, 1, 2  , 3  , 4  , 5
  // row:   1, 2, >=3, lep, had, total(once per event)
  if( numlep <= 0 && numhad <= 0 )
	return count;
  else if( numlep < 0 )
	numlep = 0;

  double* newcount = new double[npartcat];
  for(int i=0;i<npartcat;i++) { newcount[i] = count[i]; }

  int idx = -1;
  if( numlep >= 3 )
	idx = 2; // >=3 row
  else if( numlep > 0 && numlep < 3 )
	idx = numlep-1; //-1 because arrays start at zero

  if( idx >= 0 )
	newcount[idx] += weight; //rows 0-2
  if( numlep > 0 )
	newcount[npartcat-3] += weight; //lep row
  if( numhad > 0 ) 
	newcount[npartcat-2] += weight; //had row
  newcount[npartcat-1] += weight; //total row
  return newcount;
}
//end double* count_particles

const char* print_particle_table( int npartcat, double part1_1[], double part1_2[], double part2_1[], double part2_2[], double part3_1[], double part3_2[],  double part4_1[], double part4_2[]) {
  ostringstream stream;
  stream.setf(ios::fixed);
  
  stream << "\nb, t quark, Z, h event counts:\n";
  stream << "|*Num*|*b in dilep*|*b in tau*|*t in dilep*|*t in tau*|*Z in dilep*|*Z in tau*|*h in dilep*|*h in tau*|\n";
  TString rownms[] = {"|*1*|", "|*2*|", "|*>=3*|",  "|*leptonic*|", "|*hadronic*|", "|*n events*|"};
  for( int i=0; i<npartcat; i++ ) {
	stream << rownms[i];
	stream << round_char(part1_1[i]) << "|" << round_char(part1_2[i]) << "|"
		   << round_char(part2_1[i]) << "|" << round_char(part2_2[i]) << "|"
		   << round_char(part3_1[i]) << "|" << round_char(part3_2[i]) << "|"
		   << round_char(part4_1[i]) << "|" << round_char(part4_2[i]) << "|\n";
  }
  
  stream << "\n";
  return stream.str().c_str();  
}
//end print_particle_table

//Int_t LSP_mcid = 1000022; //mSUGRA
Int_t LSP_mcid = 1000039; //GMSB

bool is_stable( Int_t id ) {
  if( abs(id) >= 1 && abs(id) <= 5 ) //stable quark
	return true;
  else if( abs(id) >= 11 && abs(id) <= 16 ) //any lepton
	return true;
  else if( abs(id) == 22 || abs(id) == 21 || abs(id) == LSP_mcid ) 
	//photon, gluon (both only from h decay), LSP
	return true;
  else
	return false;
}
//end is_stable

bool daughter_charge_mismatch( Int_t mother, Int_t daughter1, Int_t daughter2){
  //this function is NOT sensitive to LSP b'c it is based on charge
  //for now, only chi0s, chi+s, gluino
  mother = abs(mother);
  daughter1 = abs(daughter1);
  daughter2 = abs(daughter2);
  //neutral gauginos: chi(2,3,4)0
  if( mother == 1000023 || mother == 1000025 || mother == 1000035 ) {
	if( daughter1 == 1000022 || daughter1 == 1000023 ) { //neutral
	  //if daughter not Z,h,gamma, then mismatch
	  if( daughter2 != 23 && daughter2 != 25 && daughter2 != 22 ) { 
		return true;
	  }
	}
	else if( daughter1 == 1000024 || daughter1 == 1000037 ) { //charged
	  if( daughter2 != 24 ) {
		return true;
	  }
	}
  }
  //charged gauginos: chi(1,2)+
  else if( mother == 1000024 || mother == 1000037 ) {
	if( daughter1 == 1000022 || daughter1 == 1000023 ) { //neutral
	  if( daughter2 != 24 ) { //daughter must be W or mismatch
		return true;
	  }
	}
  }
  else if( mother == 1000021 ) { //gluino
	//if any gaugino results, must come with 2 others
	if( daughter1 >= 1000022 && daughter1 <= 1000037 ) { 
	  if( daughter2 == 21 ) { 
		return false;
	  }
	  return true;
	}
  }
  
  return false;
}
//end daughter_charge_mismatch

vector<Int_t> make_daughter_list( vector<Int_t> list_id ) {
  //returns vector of indicies of the first daughter of each unstable particle
  vector<Int_t> list_daughter;
  for( unsigned int i=0; i<list_id.size(); i++ ) {
	list_daughter.push_back(0); //if stable, index is zero
  }

  Int_t beg_row = 6; //check for LSP, 3body decay of pair produced sparticle
  Int_t end_row = 7; 
  Int_t rownum = 1;
  Int_t rowunstable = 0;
  Int_t rowextra = 0;

  for( unsigned int i=beg_row; i<list_id.size(); i++ ) {
	if( !is_stable(list_id[i]) ) {
	  list_daughter[i] = end_row + 1 + rowunstable*2 + rowextra;
	  //if daughter+(daugter+1) don't match charge, rowextra++
	  if( daughter_charge_mismatch( list_id[i], list_id[list_daughter[i]], list_id[list_daughter[i]+1]) )
		rowextra++;
	  rowunstable++;
	}
	if( (Int_t)i == end_row ) {
	  rownum++;
	  beg_row = i+1;
	  end_row = end_row + rowunstable*2 + rowextra;
	  rowunstable = 0; 
	  rowextra = 0;
	}
  }

  return list_daughter;
}
//end make_daughter_list

//this one to be called from looper, it calls above
int num_hadronic_top( vector<Int_t> list_id ) {
  vector<Int_t> list_daughter = make_daughter_list( list_id );
  int leptop = 0; int ntop = 0;

  for( unsigned int i=0; i<list_id.size(); i++ ) {
	if( abs(list_id[i])==6 ) {
	  ntop++;
	  int gd = list_daughter[ list_daughter[i] ]; //grand-daughter of top
	  //TAUS ARE NOT COUNTED AS LEPTONS B'C ID <= 14
	  if( abs(list_id[gd]) >= 11 && abs(list_id[gd]) <= 14 ) //isa lepton
		leptop++;
	  else if( list_id[gd] == 0 || abs(list_id[ list_daughter[i] ]) != 24 ) {
		return -999;
	  }
	}
	//this was for finding gluino->gaugino+gluon, which I thought was wrong, but it is ok at 1 loop, but not tree.
	//else if( abs(list_id[i])==1000021 && abs(list_id[ list_daughter[i] ]) >= 1000022 && abs(list_id[ list_daughter[i] ]) <= 1000037 && abs(list_id[ list_daughter[i]+1 ]) == 21 ) {
	//return -998;
	//}
  }
  return ntop-leptop;
}
//end num_hadronic_top

void test_lists( int file, int run, int event, vector<Int_t> list_id ) {
  if( event < 5 ) {
	//if( run == 534 || run == 1056) {
	//if( run >= 1140 ) {
	  cout << "File: " << file << " Run: " << run << " Event: " << event << endl;
	  for( unsigned int i=0; i<list_id.size(); i++ ) {
		cout << i << ":  " << list_id[i] << endl;
	  }
	  //printing for hadronic_top test                                        

	  cout << "hadtop: " << num_hadronic_top( list_id ) << endl;;
	  //make_daughter_list( list_id );
	  cout << endl;
	}
	//}
}
//end test_lists(...)

void test_lists2( int file, int run, int event, vector<Int_t> list_id ) {
  int ht = num_hadronic_top( list_id );
  if( ht == -999 ) {
	cout << "\nFound wrong top daughter--no granddaughter\n\n";
	cout << "File: " << file << " Run: " << run << " Event: " << event << endl;
	vector<Int_t> list_daughter = make_daughter_list( list_id );
	for( unsigned int i=0; i<list_id.size(); i++ ) {
	  cout << i << ":  " << list_id[i] << endl;
	}
	cout << endl;
	for( unsigned int i=0; i<list_daughter.size(); i++ ) {
	  cout << i << ":  " << list_daughter[i] << endl;
	}
	cout << endl;
  }
  //else if( ht == -998 ) { //was error code for bad decay, but there is no problem--see num_hadronic_top
  //cout << "\nFound bad decay: gluino->gaugino+gluon\n";
  //cout << "File: " << file << " Run: " << run << " Event: " << event << endl;
  //}
}
//end test_lists2


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

//original version for 2 e,mu
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
//end int get_lep_bucket(double, double)

//use this version for one e,mu and one tau. Doen't distinguish et from mt.
int get_lep_bucket(double idlep, double idtau, const unsigned int sst,
				   const unsigned int ost) {

  if( (idlep < 0 && idtau < 0) || (idlep > 0 && idtau > 0) )
	return sst; //sst
  else 
	return ost; //ost
}
//end int get_lep_bucket(double, double, int, int)


//this function i no longer use...it's now useless b'c this info--the number of entries in each chgflv bucket--should be on the stack plot
//plus i broke it by commenting out all the global vars in warren_looper.C
void print_entries_table(TH1F* hist[], int start, int stop) {
  //one call prints table for that hist which is passed in.
  //start and stop are buckets, and inclusive
  //WARNING: the way I added the sign/flavor counting,
  //  it includes ALL buckets ALWAYS (even if not displayed)
  
  int sssf=0,ssof=0,ossf=0,osof=0;
  //use global to enforce hist vector size
  //if( stop > int(allBuckets)-1 ) { return; }

  const int numbuckets = stop - start +1;//+1 for inclusive
  double entries[numbuckets]; //display purposes
  double total_entries = 0;

  cout << "Cuts: ";
  //if( EtaCut ) cout << "Eta < " << EtaCutValue << "  ";
  //if( TausAreLeptons ) cout << " including taus \n";
  //else cout << " excluding taus \n";

  cout << "Var: " << hist[0]->GetTitle() << "\n\n";
  cout << "bucket: \t";
  for( int i = start; i <= stop; i++) {
    //cout << suffix[i] << "\t";
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

  //cout << "\n\n\tsf\tof" << "\nss:\t" << sssf << "\t" << ssof ;
  //   << "\nos:\t" << ossf << "\t" << osof;
  //cout << "\n\nTotal entries:  " << total_entries << "   Total events:  "
  //   << total_entries*cms2.evt_scale1fb();
  cout << "\n\n";
}
//end void print_entries_table(...)


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
  
