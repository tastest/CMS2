// ***************************************************************************
// This macro is made to reproduce Sync cut flow at smurf ntuple level 
// One can check following things
//  1) Bits are implemented correctly  
//  2) variables are stored correclty 
// 
// - Cut flow is made by combination of bits(variable : bit) and cuts(variable : addcut)
// - One can add bits and cuts successively to get yields 
// 
// Usage : 
// 		$ root 
// 	 	root[] .L testbit.C++
//  	root[] testbit("qqww.root")
//
// ***************************************************************************

#include "../../../Smurf/Core/SmurfTree.h"
#include "TChain.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1F.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include <iostream> 

using namespace std;

// as a reference
enum Type {
	mm,
	me,
	em,
	ee
};

/*
unsigned int FullSelection	= 	SmurfTree::BaseLine |
								SmurfTree::ChargeMatch |
								SmurfTree::Lep1FullSelection |
								SmurfTree::Lep2FullSelection |
								SmurfTree::FullMET |
								SmurfTree::ZVeto |
								SmurfTree::TopVeto |
								SmurfTree::ExtraLeptonVeto;
*/

void printYield(TChain *ch, unsigned int bit, TString addcut) {

	TH1F *hmm 	= new TH1F("hmm", 	"hmm", 	10, -0.5, 9.5);
	TH1F *hme 	= new TH1F("hme", 	"hme", 	10, -0.5, 9.5);
	TH1F *hem 	= new TH1F("hem", 	"hem", 	10, -0.5, 9.5);
	TH1F *hee 	= new TH1F("hee", 	"hee", 	10, -0.5, 9.5);
	TH1F *hall 	= new TH1F("hall", 	"hall", 10, -0.5, 9.5);

	TString bitmm = Form("((cuts & %u)==%u) && type==0", bit, bit);
	TString bitme = Form("((cuts & %u)==%u) && type==1", bit, bit);
	TString bitem = Form("((cuts & %u)==%u) && type==2", bit, bit);
	TString bitee = Form("((cuts & %u)==%u) && type==3", bit, bit);
	TString bitall = Form("((cuts & %u)==%u)", bit, bit);

	ch->Draw("min(njets,9.4999)>>hmm", 	bitmm+"&&"+addcut, "goff");
	ch->Draw("min(njets,9.4999)>>hme", 	bitme+"&&"+addcut, "goff");
	ch->Draw("min(njets,9.4999)>>hem", 	bitem+"&&"+addcut, "goff");
	ch->Draw("min(njets,9.4999)>>hee", 	bitee+"&&"+addcut, "goff");
	ch->Draw("min(njets,9.4999)>>hall", bitall+"&&"+addcut, "goff");

	cout << hmm->Integral()  << " \t" 
		 << hee->Integral()  << " \t" 
		 << hem->Integral()  << " \t" 
		 << hme->Integral()  << " \t" 
		 << hall->Integral() << endl;

	delete hmm;
	delete hme;
	delete hem;
	delete hee;
	delete hall;
}

void testbit(TString inputsmurf) {

  	TChain *ch = new TChain("tree");
    ch->Add(inputsmurf);

	// header
	cout << "mm \tee \tem \tme \tall" << endl;

	// dummy additional cut
	TString addcut="1";

	//
	unsigned int bit = SmurfTree::BaseLine;
	printYield(ch, 	bit, addcut);

	//
	bit = bit | SmurfTree::ChargeMatch;
	printYield(ch, 	bit, addcut);

	//
	bit = bit | SmurfTree::Lep1FullSelection | SmurfTree::Lep2FullSelection;
	printYield(ch, 	bit, addcut);
	
	//
 	addcut = "met>20";
	printYield(ch, 	bit, addcut);
	
	//
	bit = bit | SmurfTree::ZVeto;
	printYield(ch, 	bit, addcut);
	
	//
 	addcut = addcut+"&&min(pmet,pTrackMet)>20";
	printYield(ch, 	bit, addcut);
	
	// needs some aliases for calculate dPhi(jet1+jet2, dilep)
	// because ROOT does not support (jet1+jet2).pt() in Draw
	ch->SetAlias("jpx", 	"jet1.px()+jet2.px()");
	ch->SetAlias("jpy", 	"jet1.py()+jet2.py()");
	ch->SetAlias("dilpx", 	"dilep.px()");
	ch->SetAlias("dilpy", 	"dilep.py()");
	ch->SetAlias("jp", 		"TMath::Sqrt(jpx*jpx+jpy*jpy)");
	ch->SetAlias("dilp", 	"TMath::Sqrt(dilpx*dilpx+dilpy*dilpy)");
 	addcut = addcut+"&&(type==1 || type==2 || njets<2 || (abs(TMath::ACos((jpx*dilpx+jpy*dilpy)/jp/dilp))<(165.*TMath::Pi()/180.)))";
	printYield(ch, 	bit, addcut);

	//
 	addcut = addcut+"&&nSoftMuons==0";
	printYield(ch, 	bit, addcut);
	
	//
	bit = bit | SmurfTree::ExtraLeptonVeto;
	printYield(ch, 	bit, addcut);

	//
	bit = bit | SmurfTree::TopVeto;
	printYield(ch, 	bit, addcut);

	//
 	addcut = addcut+"&&dilep.pt()>45";
	printYield(ch, 	bit, addcut);

	//
	bit = bit | SmurfTree::FullMET;
	printYield(ch, 	bit, addcut);

	// 
	TString finalcut = addcut;

	// 0 jet bin
	cout << "---------- 0 jet ---------" << endl;
 	addcut = finalcut+"&&njets==0";
	printYield(ch, 	bit, addcut);
 	
	//
	addcut = addcut+"&&max(lep1.pt(),lep2.pt())>30";
	printYield(ch, 	bit, addcut);

	//
	addcut = addcut+"&&min(lep1.pt(),lep2.pt())>25";
	printYield(ch, 	bit, addcut);
	
	// 1 jet bin
	cout << "---------- 1 jet ---------" << endl;
 	addcut = finalcut+"&&njets==1";
	printYield(ch, 	bit, addcut);

	// 2 jet bin
	cout << "---------- 2 jet ---------" << endl;
 	addcut = finalcut+"&&(njets==2 || njets==3)";
	printYield(ch, 	bit, addcut);
 	
	//
	addcut = addcut+"&&abs(jet1.eta())<4.7 && abs(jet2.eta())<4.7";
	printYield(ch, 	bit, addcut);
 	
	//
	addcut = addcut+"&& (njets==2 || jet3.pt()<30 || (((jet1.eta()-jet3.eta())<0 && (jet2.eta()-jet3.eta())<0 ) || ((jet2.eta()-jet3.eta())>0 && (jet1.eta()-jet3.eta())>0)))";
	printYield(ch, 	bit, addcut);

}
