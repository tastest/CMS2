#include <algorithm>
#include <set>
#include "Math/VectorUtil.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TList.h"
#include "TRegexp.h"
#include "TDirectory.h"
#include "tools.h"
#include "DileptonHypType.h"
#include "CORE/CMS2.h"

// use this directory to store our histograms (so they can be written to file in one fell swoop)
TDirectory *histo_directory = new TDirectory("histo_directory", "directory for all our histograms");


unsigned int encodeTriLeptonCand(unsigned int bucket,unsigned int first, unsigned int second, unsigned int third) {
  // encode trilepton candidate according to
  //
  // trilepton candidate is identified by coded unsigned integer: AABBCCDD
  //
  // AA: bucket enumerator
  // BB: index of first lepton from els and mus collection
  // CC: index of second lepton from els and mus collection
  // DD: index of third lepton from els and mus collection
  //
  // flavors are ordered: first mu, then el
  // same flavor leptons in candidate are ordered by increasing index
  //

  return bucket * 1000000 + first * 10000 + second * 100 + third;
}

unsigned int decodeBucket(unsigned int cand) {
  // decode bucket enumerator
  // 
  // trilepton candidate is identified by coded unsigned integer: AABBCCDD
  //
  // AA: bucket enumerator
  // BB: index of first lepton from els and mus collection
  // CC: index of second lepton from els and mus collection
  // DD: index of third lepton from els and mus collection
  //
  // flavors are ordered: first mu, then el
  // same flavor leptons in candidate are ordered by increasing index
  //

  return cand/1000000;
}

unsigned int decodeFirst(unsigned int cand) {
  // decode index of first lepton
  // 
  // trilepton candidate is identified by coded unsigned integer: AABBCCDD
  //
  // AA: bucket enumerator
  // BB: index of first lepton from els and mus collection
  // CC: index of second lepton from els and mus collection
  // DD: index of third lepton from els and mus collection
  //
  // flavors are ordered: first mu, then el
  // same flavor leptons in candidate are ordered by increasing index
  //

  return (cand/10000 - (cand/1000000*100));
}

unsigned int decodeSecond(unsigned int cand) {
  // decode index of first lepton
  // 
  // trilepton candidate is identified by coded unsigned integer: AABBCCDD
  //
  // AA: bucket enumerator
  // BB: index of first lepton from els and mus collection
  // CC: index of second lepton from els and mus collection
  // DD: index of third lepton from els and mus collection
  //
  // flavors are ordered: first mu, then el
  // same flavor leptons in candidate are ordered by increasing index
  //

  return (cand/100 - (cand/10000*100));
}

unsigned int decodeThird(unsigned int cand) {
  // decode index of first lepton
  // 
  // trilepton candidate is identified by coded unsigned integer: AABBCCDD
  //
  // AA: bucket enumerator
  // BB: index of first lepton from els and mus collection
  // CC: index of second lepton from els and mus collection
  // DD: index of third lepton from els and mus collection
  //
  // flavors are ordered: first mu, then el
  // same flavor leptons in candidate are ordered by increasing index
  //

  return (cand - (cand/100*100));
}

TH1F* book1DHist(const char* name, const char* title, unsigned int nbins, float low, float high, const char* xtitle, const char* ytitle, int color) {
  // return histogram instance with called Sumw2
  TH1F *hist = new TH1F(name,title,nbins,low,high);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->Sumw2();
  hist->SetFillColor(color);
  hist->SetLineColor(color);
   
  return hist;   
}

TH1F* book1DVarHist(const char* name, const char* title, unsigned int nbins, float* bins, const char* xtitle, const char* ytitle, int color) {
  // return histogram instance with called Sumw2
  TH1F *hist = new TH1F(name,title,nbins,bins);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->Sumw2();
  hist->SetFillColor(color);
  hist->SetLineColor(color);
   
  return hist;   
}

TH1F* book1DVarHist(const char* name, const char* title, vector<float> &bins, const char* xtitle, const char* ytitle, int color) {
  // return histogram instance with called Sumw2
  const unsigned int nBins = bins.size()-1;
  float binArray[nBins+1];
  for (unsigned int i = 0;
       i < nBins+1;
       ++i) {
    binArray[i] = bins[i];
  }

  TH1F *hist = new TH1F(name,title,nBins,binArray);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->Sumw2();
  hist->SetFillColor(color);
  hist->SetLineColor(color);
   
  return hist;
}

TH2F* book2DHist(const char* name, const char* title, unsigned int nxbins, float xlow, float xhigh, unsigned int nybins, float ylow, float yhigh, const char* xtitle, const char* ytitle, const char* ztitle, int color) {
  // return histogram instance with called Sumw2
  TH2F *hist = new TH2F(name,title,nxbins,xlow,xhigh,nybins,ylow,yhigh);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->SetZTitle(ztitle);
  hist->Sumw2();
  hist->SetFillColor(color);
  hist->SetLineColor(color);
   
  return hist;   
}

TH2F* book2DVarHist(const char* name, const char* title, unsigned int nxbins, float* xbins, unsigned int nybins, float* ybins, const char* xtitle, const char* ytitle, const char* ztitle, int color) {
  // return histogram instance with called Sumw2
  TH2F *hist = new TH2F(name,title,nxbins,xbins,nybins,ybins);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->SetZTitle(ztitle);
  hist->Sumw2();
  hist->SetFillColor(color);
  hist->SetLineColor(color);
   
  return hist;   
}

TH2F* book2DVarHist(const char* name, const char* title, vector<float> &xbins, vector<float> &ybins, const char* xtitle, const char* ytitle, const char* ztitle, int color) {
  // return histogram instance with called Sumw2
  const unsigned int nxBins = xbins.size()-1;
  float xbinArray[nxBins+1];
  for (unsigned int i = 0;
       i < nxBins+1;
       ++i) {
    xbinArray[i] = xbins[i];
  }
  const unsigned int nyBins = ybins.size()-1;
  float ybinArray[nyBins+1];
  for (unsigned int i = 0;
       i < nyBins+1;
       ++i) {
    ybinArray[i] = ybins[i];
  }

  TH2F *hist = new TH2F(name,title,nxBins,xbinArray,nyBins,ybinArray);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->SetZTitle(ztitle);
  hist->Sumw2();
  hist->SetFillColor(color);
  hist->SetLineColor(color);
   
  return hist;
}

TH3F* book3DHist(const char* name, const char* title, 
		 unsigned int nxbins, float xlow, float xhigh, 
		 unsigned int nybins, float ylow, float yhigh, 
		 unsigned int nzbins, float zlow, float zhigh, 
		 const char* xtitle, const char* ytitle, const char* ztitle, int color) {
  // return histogram instance with called Sumw2
  TH3F *hist = new TH3F(name,title,nxbins,xlow,xhigh,nybins,ylow,yhigh,nzbins,zlow,zhigh);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->SetZTitle(ztitle);
  hist->Sumw2();
  hist->SetFillColor(color);
  hist->SetLineColor(color);
   
  return hist;   
}

TH3F* book3DVarHist(const char* name, const char* title, 
		    unsigned int nxbins, float* xbins, 
		    unsigned int nybins, float* ybins, 
		    unsigned int nzbins, float* zbins, 
		    const char* xtitle, const char* ytitle, const char* ztitle, int color) {
  // return histogram instance with called Sumw2
  TH3F *hist = new TH3F(name,title,nxbins,xbins,nybins,ybins,nzbins,zbins);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->SetZTitle(ztitle);
  hist->Sumw2();
  hist->SetFillColor(color);
  hist->SetLineColor(color);
   
  return hist;   
}

TH3F* book3DVarHist(const char* name, const char* title, 
		    vector<float> &xbins, 
		    vector<float> &ybins, 
		    vector<float> &zbins, 
		    const char* xtitle, const char* ytitle, const char* ztitle, int color) {
  // return histogram instance with called Sumw2
  const unsigned int nxBins = xbins.size()-1;
  float xbinArray[nxBins+1];
  for (unsigned int i = 0;
       i < nxBins+1;
       ++i) {
    xbinArray[i] = xbins[i];
  }
  const unsigned int nyBins = ybins.size()-1;
  float ybinArray[nyBins+1];
  for (unsigned int i = 0;
       i < nyBins+1;
       ++i) {
    ybinArray[i] = ybins[i];
  }

  const unsigned int nzBins = zbins.size()-1;
  float zbinArray[nzBins+1];
  for (unsigned int i = 0;
       i < nzBins+1;
       ++i) {
    zbinArray[i] = zbins[i];
  }

  TH3F *hist = new TH3F(name,title,nxBins,xbinArray,nyBins,ybinArray,nzBins,zbinArray);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->SetZTitle(ztitle);
  hist->Sumw2();
  hist->SetFillColor(color);
  hist->SetLineColor(color);
   
  return hist;
}

float mee(int i, int j)
{
  if ( cms2.els_charge()[i] * cms2.els_charge()[j] < 0 ) {
    LorentzVector vec = cms2.els_p4()[i] + cms2.els_p4()[j];
    return vec.mass();
  } else {
    cout << "ERROR: mee tries to calculate Z mass from same charge leptons" << endl;
    return 999999.;
  }
}

float mmm(int i, int j)
{
  if ( cms2.mus_charge()[i] * cms2.mus_charge()[j] < 0 ) {
    LorentzVector vec = cms2.mus_p4()[i] + cms2.mus_p4()[j];
    return vec.mass();
  } else {
    cout << "ERROR: mmm tries to calculate Z mass from same charge leptons" << endl;
    return 999999.;
  }
}

#if 0
bool goodLeptonIsolated(int bucket, int first, int second, int third) {
  // check if all 3 leptons are good isolated leptons
  //

  bool result = false;
  if ( bucket == 16 ||
       bucket == 17 ||
       bucket == 18 ||
       bucket == 19 ) {
    result = goodElectronIsolated(first) && goodElectronIsolated(second) && goodElectronIsolated(third);
  } else if ( bucket == 0 ||
	      bucket == 1 ||
	      bucket == 2 ||
	      bucket == 3 ) {
    result = goodMuonIsolated(first) && goodMuonIsolated(second) && goodMuonIsolated(third);
  } else if ( bucket == 10 ||
	      bucket == 13 ||
	      bucket == 14 ||
	      bucket == 12 ||
	      bucket == 15 ||
	      bucket == 11 ) {
    result = goodMuonIsolated(first) && goodElectronIsolated(second) && goodElectronIsolated(third);
  } else if ( bucket == 4 ||
	      bucket == 5 ||
	      bucket == 7 ||
	      bucket == 6 ||
	      bucket == 9 ||
	      bucket == 8 ) {
    result = goodMuonIsolated(first) && goodMuonIsolated(second) && goodElectronIsolated(third);
  } else {
    cout << "ERROR: goodIsolatedLepton tried to use not existing bucket!" << endl;
  }

  return result;
}
#endif

float ptLowestPtLepton(int bucket, int first, int second, int third) {
  // determine pt of lowest pt lepton
  float array[3] = {0,0,0};

  if ( bucket == 16 ||
       bucket == 17 ||
       bucket == 18 ||
       bucket == 19 ) {
    array[0] = cms2.els_p4()[first].pt();
    array[1] = cms2.els_p4()[second].pt();
    array[2] = cms2.els_p4()[third].pt();
  } else if ( bucket == 0 ||
	      bucket == 1 ||
	      bucket == 2 ||
	      bucket == 3 ) {
    array[0] = cms2.mus_p4()[first].pt();
    array[1] = cms2.mus_p4()[second].pt();
    array[2] = cms2.mus_p4()[third].pt();
  } else if ( bucket == 10 ||
	      bucket == 13 ||
	      bucket == 14 ||
	      bucket == 12 ||
	      bucket == 15 ||
	      bucket == 11 ) {
    array[0] = cms2.mus_p4()[first].pt();
    array[1] = cms2.els_p4()[second].pt();
    array[2] = cms2.els_p4()[third].pt();
  } else if ( bucket == 4 ||
	      bucket == 5 ||
	      bucket == 7 ||
	      bucket == 6 ||
	      bucket == 9 ||
	      bucket == 8 ) {
    array[0] = cms2.mus_p4()[first].pt();
    array[1] = cms2.mus_p4()[second].pt();
    array[2] = cms2.els_p4()[third].pt();
  } else {
    cout << "ERROR: ptLowestPtLepton tried to use not existing bucket!" << endl;
  }

  // sort array
  sort(array,array+3);

  if ( array[0] > 29. ) {
    return 25.;
  }

  return array[0];

}

float* triLeptonPtArray(int bucket, int first, int second, int third) {
  //
  // returns pt of leptons in trilepton cand.
  //
  float *array = new float[3];

  if ( bucket == 16 ||
       bucket == 17 ||
       bucket == 18 ||
       bucket == 19 ) {
    array[0] = cms2.els_p4()[first].pt();
    array[1] = cms2.els_p4()[second].pt();
    array[2] = cms2.els_p4()[third].pt();
  } else if ( bucket == 0 ||
	      bucket == 1 ||
	      bucket == 2 ||
	      bucket == 3 ) {
    array[0] = cms2.mus_p4()[first].pt();
    array[1] = cms2.mus_p4()[second].pt();
    array[2] = cms2.mus_p4()[third].pt();
  } else if ( bucket == 10 ||
	      bucket == 13 ||
	      bucket == 14 ||
	      bucket == 12 ||
	      bucket == 15 ||
	      bucket == 11 ) {
    array[0] = cms2.mus_p4()[first].pt();
    array[1] = cms2.els_p4()[second].pt();
    array[2] = cms2.els_p4()[third].pt();
  } else if ( bucket == 4 ||
	      bucket == 5 ||
	      bucket == 7 ||
	      bucket == 6 ||
	      bucket == 9 ||
	      bucket == 8 ) {
    array[0] = cms2.mus_p4()[first].pt();
    array[1] = cms2.mus_p4()[second].pt();
    array[2] = cms2.els_p4()[third].pt();
  } else {
    cout << "ERROR: triLeptonPtArray tried to use not existing bucket!" << endl;
  }

  return array;
}


bool passTriggerLeptonMinPtCut(int bucket, int first, int second, int third, float triggerLeptonMinPtCut) {
  // return true if at least one lepton in trilepton candidate passes pt > triggerLeptonMinPtCut

  float* array = triLeptonPtArray(bucket,first,second,third);

  // sort array
  sort(array,array+3);

  return array[2] > triggerLeptonMinPtCut;

}
TString printCand(int bucket, int first, int second, int third) {
  // print candidate based on bucket
  //

  TString result = "";

  if ( bucket == 16 ||
       bucket == 17 ||
       bucket == 18 ||
       bucket == 19 ) {

    result.Append("electron 1: " );
    result += cms2.els_p4()[first].pt();
    result.Append(" electron 2: " );
    result += cms2.els_p4()[second].pt();
    result.Append(" electron 3: " );
    result += cms2.els_p4()[third].pt();
  } else if ( bucket == 0 ||
	      bucket == 1 ||
	      bucket == 2 ||
	      bucket == 3 ) {
    result.Append("muon 1: " );
    result += cms2.mus_p4()[first].pt();
    result.Append(" muon 2: " );
    result += cms2.mus_p4()[second].pt();
    result.Append(" muon 3: " );
    result += cms2.mus_p4()[third].pt();
  } else if ( bucket == 10 ||
	      bucket == 13 ||
	      bucket == 14 ||
	      bucket == 12 ||
	      bucket == 15 ||
	      bucket == 11 ) {
    result.Append("muon 1: " );
    result += cms2.mus_p4()[first].pt();
    result.Append(" electron 2: " );
    result += cms2.els_p4()[second].pt();
    result.Append(" electron 3: " );
    result += cms2.els_p4()[third].pt();
  } else if ( bucket == 4 ||
	      bucket == 5 ||
	      bucket == 7 ||
	      bucket == 6 ||
	      bucket == 9 ||
	      bucket == 8 ) {
    result.Append("muon 1: " );
    result += cms2.mus_p4()[first].pt();
    result.Append(" muon 2: " );
    result += cms2.mus_p4()[second].pt();
    result.Append(" electron 3: " );
    result += cms2.els_p4()[third].pt();
  } else {
    cout << "ERROR: goodIsolatedLepton tried to use not existing bucket!" << endl;
  }

  return result;
}
 
void calcPrimZ(int bucket, int first, int second, int third, float zmass, unsigned int* array) {
  // caluclate primary Z candidate for trilepton candidate closest to Z mass
  // determine array holding indiced (first, second and third)
  // first two positions leptons belonging to primary Z candidate
  // third position is the not used lepton
  // 999999 on any of the three positions denotes that no primary Z candidate was found

  float mz1 = 0;
  float  mz2 = 0;
  if ( bucket == 17 && cms2.els_charge()[first] == -1 ){ //found the odd one out
    mz1 = mee(first, second);
    mz2 = mee(first, third);
    array[0] = first;
    array[1] = second;
    array[2] = third;    
    if( TMath::Abs(mz2-zmass) < TMath::Abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
      float tmp = mz1;
      mz1 = mz2;
      mz2 = tmp;
      array[0] = first;
      array[1] = third;
      array[2] = second;    
    }
  } else if ( bucket == 17 && cms2.els_charge()[second] == -1 ){ //found the odd one out
    mz1 = mee(second, first);
    mz2 = mee(second, third);
    array[0] = second;
    array[1] = first;
    array[2] = third;    
    if( TMath::Abs(mz2-zmass) < TMath::Abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
      float tmp = mz1;
      mz1 = mz2;
      mz2 = tmp;
      array[0] = second;
      array[1] = third;
      array[2] = first;    
    }
  } else if ( bucket == 17 && cms2.els_charge()[third] == -1 ){ //found the odd one out
    mz1 = mee(third, first);
    mz2 = mee(third, second);
    array[0] = third;
    array[1] = first;
    array[2] = second;    
    if( TMath::Abs(mz2-zmass) < TMath::Abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
      float tmp = mz1;
      mz1 = mz2;
      mz2 = tmp;
      array[0] = third;
      array[1] = second;
      array[2] = first;    
    }
  } else if ( bucket == 18 && cms2.els_charge()[first] == 1 ){ //found the odd one out
    mz1 = mee(first, second);
    mz2 = mee(first, third);
    array[0] = first;
    array[1] = second;
    array[2] = third;    
    if( TMath::Abs(mz2-zmass) < TMath::Abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
      float tmp = mz1;
      mz1 = mz2;
      mz2 = tmp;
      array[0] = first;
      array[1] = third;
      array[2] = second;    
    }
  } else if ( bucket == 18 && cms2.els_charge()[second] == 1 ){ //found the odd one out
    mz1 = mee(second, first);
    mz2 = mee(second, third);
    array[0] = second;
    array[1] = first;
    array[2] = third;    
    if( TMath::Abs(mz2-zmass) < TMath::Abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
      float tmp = mz1;
      mz1 = mz2;
      mz2 = tmp;
      array[0] = second;
      array[1] = third;
      array[2] = first;    
    }
  } else if ( bucket == 18 && cms2.els_charge()[third] == 1 ){ //found the odd one out
    mz1 = mee(third, first);
    mz2 = mee(third, second);
    array[0] = third;
    array[1] = first;
    array[2] = second;    
    if( TMath::Abs(mz2-zmass) < TMath::Abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
      float tmp = mz1;
      mz1 = mz2;
      mz2 = tmp;
      array[0] = third;
      array[1] = second;
      array[2] = first;    
    }
  } else if ( bucket == 1 && cms2.mus_charge()[first] == -1 ){ //found the odd one out
    mz1 = mmm(first, second);
    mz2 = mmm(first, third);
    array[0] = first;
    array[1] = second;
    array[2] = third;    
    if( TMath::Abs(mz2-zmass) < TMath::Abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
      float tmp = mz1;
      mz1 = mz2;
      mz2 = tmp;
      array[0] = first;
      array[1] = third;
      array[2] = second;    
    }
  } else if ( bucket == 1 && cms2.mus_charge()[second] == -1 ){ //found the odd one out
    mz1 = mmm(second, first);
    mz2 = mmm(second, third);
    array[0] = second;
    array[1] = first;
    array[2] = third;    
    if( TMath::Abs(mz2-zmass) < TMath::Abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
      float tmp = mz1;
      mz1 = mz2;
      mz2 = tmp;
      array[0] = second;
      array[1] = third;
      array[2] = first;    
    }
  } else if ( bucket == 1 && cms2.mus_charge()[third] == -1 ){ //found the odd one out
    mz1 = mmm(third, first);
    mz2 = mmm(third, second);
    array[0] = third;
    array[1] = first;
    array[2] = second;    
    if( TMath::Abs(mz2-zmass) < TMath::Abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
      float tmp = mz1;
      mz1 = mz2;
      mz2 = tmp;
      array[0] = third;
      array[1] = second;
      array[2] = first;    
    }
  } else if ( bucket == 2 && cms2.mus_charge()[first] == 1 ){ //found the odd one out
    mz1 = mmm(first, second);
    mz2 = mmm(first, third);
    array[0] = first;
    array[1] = second;
    array[2] = third;    
    if( TMath::Abs(mz2-zmass) < TMath::Abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
      float tmp = mz1;
      mz1 = mz2;
      mz2 = tmp;
      array[0] = first;
      array[1] = third;
      array[2] = second;    
    }
  } else if ( bucket == 2 && cms2.mus_charge()[second] == 1 ){ //found the odd one out
    mz1 = mmm(second, first);
    mz2 = mmm(second, third);
    array[0] = second;
    array[1] = first;
    array[2] = third;    
    if( TMath::Abs(mz2-zmass) < TMath::Abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
      float tmp = mz1;
      mz1 = mz2;
      mz2 = tmp;
      array[0] = second;
      array[1] = third;
      array[2] = first;    
    }
  } else if ( bucket == 2 && cms2.mus_charge()[third] == 1 ){ //found the odd one out
    mz1 = mmm(third, first);
    mz2 = mmm(third, second);
    array[0] = third;
    array[1] = first;
    array[2] = second;    
    if( TMath::Abs(mz2-zmass) < TMath::Abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
      float tmp = mz1;
      mz1 = mz2;
      mz2 = tmp;
      array[0] = third;
      array[1] = second;
      array[2] = first;    
    }
  } else if ( bucket == 14 || bucket == 11 ) {
    // float mz1 = mee(second, third);
    array[0] = second;
    array[1] = third;
    array[2] = first;    
  } else if ( bucket == 7 || bucket == 6 ) {
    // float mz1 = mmm(first, second);
    array[0] = first;
    array[1] = second;
    array[2] = third;    
  }
}

float calcPrimZMass(int bucket, int first, int second) {
  // calculate mass of primary Z candidate
  float result = 999999.;

  if ( bucket == 17 ||
       bucket == 18 ||
       bucket == 14 ||
       bucket == 11 ) {
    result = mee(first,second);
  } else if (bucket == 1 ||
	     bucket == 2 ||
	     bucket == 7 ||
	     bucket == 6 )  {
    result = mmm(first,second);
  }

  return result;

}

bool passMETCut(int bucket, float metAll, float electronMETCut, float muonMETCut) {
  // cut on MET 
  // use electronMETCut if lepton not belonging to the primary Z is an electron
  // use muonMETCut if lepton not belonging to the primary Z is an muon

  bool result = false;
  
  if ( bucket == 17 ||
       bucket == 18 ||
       bucket == 7 ||
       bucket == 6 ) {
    result = metAll > electronMETCut;
  } else if ( bucket == 14 ||
	      bucket == 11 ||
	      bucket == 1 ||
	      bucket == 2 ) {
    result = metAll > muonMETCut;
  }

  return result;
}

#if 0
float calcSecZMass(int bucket, unsigned int* leptons, float zmass) {
  // combine lepton not belonging to primary Z candidate with all other good leptons without isolation to form Z candidates (OC, Same flavor)
  // return mass closest to zmass
  // leptons is array from calcPrimZ

  float result = 999999.;
  
  if ( bucket == 17 ||
       bucket == 18 ) {
    if ( cms2.evt_nels() > 3 ) {
      for ( unsigned int i = 0; i < cms2.evt_nels(); ++i ) {
	if ( i != leptons[0] && i != leptons[1] && i != leptons[2]) {
	  if ( goodElectronWithoutIsolation(i) ) {
	    if ( cms2.els_charge()[i] * cms2.els_charge()[leptons[2]] < 0 ) {
	      float mass = mee(i,leptons[2]);
	      // store Z mass closest to true Z mass
	      if ( TMath::Abs(zmass-mass) < TMath::Abs(zmass-result) ) {
		result = mass;
	      }
	    }
	  }
	}
      }
    }
  } else if ( bucket == 1 ||
	      bucket == 2 ) {
    if ( cms2.mus_p4().size() > 3 ) {
      for ( unsigned int i = 0; i < cms2.mus_p4().size(); ++i ) {
	if ( i != leptons[0] && i != leptons[1] && i != leptons[2]) {
	  if ( goodMuonWithoutIsolation(i) ) {
	    if ( cms2.mus_charge()[i] * cms2.mus_charge()[leptons[2]] < 0 ) {
	      float mass = mmm(i,leptons[2]);
	      // store Z mass closest to true Z mass
	      if ( TMath::Abs(zmass-mass) < TMath::Abs(zmass-result) ) {
		result = mass;
	      }
	    }
	  }
	}
      }
    }
  } else if ( bucket == 14 ||
	      bucket == 11) {
    if ( cms2.mus_p4().size() > 1 ) {
      for ( unsigned int i = 0; i < cms2.mus_p4().size(); ++i ) {
	if ( i != leptons[2] ) {
	  if ( goodMuonWithoutIsolation(i) ) {
	    if ( cms2.mus_charge()[i] * cms2.mus_charge()[leptons[2]] < 0 ) {
	      float mass = mmm(i,leptons[2]);
	      // store Z mass closest to true Z mass
	      if ( TMath::Abs(zmass-mass) < TMath::Abs(zmass-result) ) {
		result = mass;
	      }
	    }
	  }
	}
      }
    }
  } else if ( bucket == 7 ||
	      bucket == 6 ) {
    if ( cms2.evt_nels() > 1 ) {
      for ( unsigned int i = 0; i < cms2.evt_nels(); ++i ) {
	if ( i != leptons[2] ) {
	  if ( goodElectronWithoutIsolation(i) ) {
	    if ( cms2.els_charge()[i] * cms2.els_charge()[leptons[2]] < 0 ) {
	      float mass = mee(i,leptons[2]);
	      // store Z mass closest to true Z mass
	      if ( TMath::Abs(zmass-mass) < TMath::Abs(zmass-result) ) {
		result = mass;
	      }
	    }
	  }
	}
      }
    }
  }
  
  return result;

}
#endif

bool testJetForElectrons(const LorentzVector& jetP4, const LorentzVector& elP4) {
  
  
  bool matched = false;
  float elphi  = elP4.Phi();
  float jetphi = jetP4.Phi();
   
  float eleta  = elP4.Eta();
  float jeteta = jetP4.Eta();
   
  float dphi = elphi - jetphi;
  float deta = eleta - jeteta;
  if(TMath::Abs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - TMath::Abs(dphi);
   
  double dR = sqrt(dphi*dphi + deta*deta);
  if (dR < 0.4) 
    matched = true;
  
  return !matched;
}

std::vector<LorentzVector> correctJetsForElectrons(int bucket, int first, int second, int third) {

  std::vector<LorentzVector> result;
#if 0
  double hypJetMinEtaCut = -3.0;
  double hypJetMaxEtaCut = 3.0;
  double hypJetMinPtCut  = 15;

  for ( unsigned int jet = 0;
	jet < cms2.jets_p4().size();
	++jet ) {
    if ( bucket == 4 ||
	 bucket == 5 ||
	 bucket == 6 ||
	 bucket == 7 ||
	 bucket == 8 ||
	 bucket == 9 ) {
      if ( !testJetForElectrons(cms2.jets_p4()[jet],cms2.els_p4()[third]) ) continue;
    } else if ( bucket == 10 ||
		bucket == 11 ||
		bucket == 12 ||
		bucket == 13 ||
		bucket == 14 ||
		bucket == 15 ) {
      if ( !testJetForElectrons(cms2.jets_p4()[jet],cms2.els_p4()[second]) ) continue;
      if ( !testJetForElectrons(cms2.jets_p4()[jet],cms2.els_p4()[third]) ) continue;
    } else if ( bucket == 16 ||
		bucket == 17 ||
		bucket == 18 ||
		bucket == 19 ) {
      if ( !testJetForElectrons(cms2.jets_p4()[jet],cms2.els_p4()[first]) ) continue;
      if ( !testJetForElectrons(cms2.jets_p4()[jet],cms2.els_p4()[second]) ) continue;
      if ( !testJetForElectrons(cms2.jets_p4()[jet],cms2.els_p4()[third]) ) continue;
    }

    if ( cms2.jets_p4()[jet].eta() >= hypJetMaxEtaCut ) continue;
    if ( cms2.jets_p4()[jet].eta() <= hypJetMinEtaCut ) continue;
    if ( cms2.jets_p4()[jet].Pt() <= hypJetMinPtCut/cms2.jets_tq_noCorrF()[jet] ) continue;

    result.push_back(cms2.jets_p4()[jet]);

  }
#endif
  return result;

}

void progressBar(int& i_permille_old, int nEventsTotal, int nEventsChain) {
  // 
  // progress bar
  //
  int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
  if (i_permille != i_permille_old) {
    // xterm magic from L. Vacavant and A. Cerri
    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
	   "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
    fflush(stdout);
    i_permille_old = i_permille;
  }
}


void correctMETmuons_crossedE(double& met, double& metPhi, 
			      double muon_pt, double muon_phi,
			      double muon_track_theta, double muon_track_phi,
			      double mu_crossedem_dep, double mu_crossedhad_dep, double mu_crossedho_dep ) {

  // first, account for muon momentum
  double metx =  met*cos(metPhi);
  double mety =  met*sin(metPhi);
  double pt0  =  muon_pt; 
  double phi0 =  muon_phi; 
  metx -= pt0*cos(phi0);
  mety -= pt0*sin(phi0);
  
  
  met = sqrt(metx*metx+mety*mety);
  metPhi = atan2(mety, metx);
   
   double muEx = 0.0;
   double muEy = 0.0;
   
   
   // use muon position at the outer most state of the silicon track if 
   // TrackExtra is available and momentum direction at the origin 
   // otherwise. Both should be fine.
   // NOTICE: MET is built out of towers, which are 5x5 ECAL crystals + one 
   // element of HCAL and HO. Muon energy is reported for individual crossed 
   // elements of all the detectors and 3x3 elements of each time as an 
   // alternative way of energy calculation.
   double theta = muon_track_theta;
   double phi   = muon_track_phi;
	 
   muEx += ( mu_crossedem_dep + mu_crossedhad_dep + mu_crossedho_dep )*sin(theta)*cos( phi );
   muEy += ( mu_crossedem_dep + mu_crossedhad_dep + mu_crossedho_dep )*sin(theta)*sin( phi );
   
   
   metx = met*cos(metPhi) + muEx;
   mety = met*sin(metPhi) + muEy;
   met   = sqrt(metx*metx + mety*mety);
   metPhi = atan2(mety, metx);
}

DorkyEventIdentifier::DorkyEventIdentifier (CMS2 &cms2)
     : run(cms2.evt_run()),
       event(cms2.evt_event()),
       lumi_section(cms2.evt_lumiBlock()),
       trks_d0(cms2.trks_d0().at(0)), // use at() because that way we
				      // get an exception if we are
				      // out of bounds
       trks_pt (cms2.trks_trk_p4().at(0).pt()),
       trks_eta(cms2.trks_trk_p4().at(0).eta()),
       trks_phi(cms2.trks_trk_p4().at(0).phi())
{
     
}

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
     if (run != other.run)
          return run < other.run;
     if (event != other.event)
          return event < other.event;
     if (lumi_section != other.lumi_section)
          return lumi_section < other.lumi_section;
     if (TMath::Abs(trks_d0 - other.trks_d0) > 1e-6 * fabs(trks_d0))
       return trks_d0 < other.trks_d0;
     if (TMath::Abs(trks_pt - other.trks_pt) > 1e-6 * fabs(trks_pt))
       return trks_pt < other.trks_pt;
     if (TMath::Abs(trks_eta - other.trks_eta) > 1e-6 * fabs(trks_eta))
       return trks_eta < other.trks_eta;
     if (TMath::Abs(trks_phi - other.trks_phi) > 1e-6 * fabs(trks_phi))
       return trks_phi < other.trks_phi;
     return false;
}

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
     if (run != other.run)
          return false;
     if (event != other.event)
          return false;
     if (lumi_section != other.lumi_section)
          return false;
     if (TMath::Abs(trks_d0 - other.trks_d0) > 1e-6 * fabs(trks_d0))
          return false;
     if (TMath::Abs(trks_pt - other.trks_pt) > 1e-6 * fabs(trks_pt))
          return false;
     if (TMath::Abs(trks_eta - other.trks_eta) > 1e-6 * fabs(trks_eta))
          return false;
     if (TMath::Abs(trks_phi - other.trks_phi) > 1e-6 * fabs(trks_phi))
          return false;
     return true;
}

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id)
{
     std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
          already_seen.insert(id);
     return !ret.second;
}

Status3Identifier::Status3Identifier (CMS2 &cms2)
     : run(cms2.evt_run()),
       event(cms2.evt_event()),
       lumi_section(cms2.evt_lumiBlock())
{
//      pts.reserve(cms2.genps_p4().size());
     for (unsigned int i = 0; i < std::min(size_t(4), cms2.genps_p4().size()); ++i) {
	  pts.insert(cms2.genps_p4()[i].pt());
	  // only go up to the W for this exercise
	  if (i + 1 < cms2.genps_p4().size() && 
	      abs(cms2.genps_id()[i + 1]) == 24)
	       break;
     }
}

bool Status3Identifier::operator < (const Status3Identifier &other) const
{
     std::multiset<float>::const_iterator ipts = pts.begin();
     std::multiset<float>::const_iterator iother = other.pts.begin();
     for (unsigned int i = 0; i < std::max(pts.size(), other.pts.size()); 
	  ++i, ++ipts, ++iother) {
	  if (i >= pts.size())
	       return true;
	  if (i >= other.pts.size())
	       return false;
// 	  if (fabs(*ipts - *iother) > 1e-6 * fabs(*ipts)) 
// 	       return *ipts < *iother;
	  if (*ipts != *iother)
	       return *ipts < *iother;
     }
     return false;
}

std::set<Status3Identifier> already_seen_stat3;
bool is_duplicate (const Status3Identifier &id)
{
     std::pair<std::set<Status3Identifier>::const_iterator, bool> ret =
          already_seen_stat3.insert(id);
     if (!ret.second) {
	  std::cout << std::setprecision(8);
	  std::cout << "new entry (" << id.run << " " << id.lumi_section << " " << id.event << ") ";
	  for (std::multiset<float>::const_iterator i = id.pts.begin();
	       i != id.pts.end(); ++i) {
	       std::cout << *i << " ";
	  }
	  std::cout << std::endl;
	  std::cout << "old entry (" << ret.first->run << " " << ret.first->lumi_section << " " << ret.first->event << ") ";
	  for (std::multiset<float>::const_iterator i = ret.first->pts.begin();
	       i != ret.first->pts.end(); ++i) {
	       std::cout << *i << " ";
	  }
	  std::cout << std::endl;
	  return true;
     }
     return false;
}

enum DileptonHypType hyp_typeToHypType (int hyp_type)
{
     switch (hyp_type) {
     case 0:
	  return DILEPTON_MUMU;
     case 1: case 2:
	  return DILEPTON_EMU;
     case 3:
	  return DILEPTON_EE;
     default:
	  assert(hyp_type < 4);
     }
     return DILEPTON_ALL;
}

// save histograms to outfile
void saveHist(const char* filename, const char* pat)
{
     TList* list = histo_directory->GetList() ;
     TIterator* iter = list->MakeIterator();
     
     //TRegexp re(pat,kTRUE) ;
     
     TRegexp re(pat, 1) ;
     
     TFile outf(filename,"RECREATE") ;
     while(TObject *obj=iter->Next()) {
	  if (TString(obj->GetName()).Index(re)>=0) {
	       obj->Write() ;
	       cout << "." ;
	       cout.flush() ;
	  }
     }
     cout << endl ;
     outf.Close() ;
     
     delete iter ;
}

double correctd0Phi(int trk_idx) 
{
     return 0;
//   return atan2( -1 * trks_d0corr * sin( trks_trk_p4->phi() ), trks_d0corr * cos( trks_trk_p4->phi() ) );
}


//port of code from CaloTowerDetId.h
int CaloTwr_ieta( int detid ) {
  int zside = (detid & 0x2000) ? 1 : -1 ;
  //warning: multiplication has higher precendnce than bitwise anding
  return zside * ((detid >> 7) & 0x3f);
}

int CaloTwr_iphi( int detid ) {
  return detid & 0x7F;
}


