#include <algorithm>
#include "Math/VectorUtil.h"

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

TH1F* book1DHist(const char* name, const char* title, unsigned int nbins, float low, float high, const char* xtitle, const char* ytitle) {
  // return histogram instance with called Sumw2
  TH1F *hist = new TH1F(name,title,nbins,low,high);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->Sumw2();
   
  return hist;   
}

TH1F* book1DVarHist(const char* name, const char* title, unsigned int nbins, float* bins, const char* xtitle, const char* ytitle) {
  // return histogram instance with called Sumw2
  TH1F *hist = new TH1F(name,title,nbins,bins);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->Sumw2();
   
  return hist;   
}

float mee(int i, int j){
  if ( els_charge[i] * els_charge[j] < 0 ) {
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
      vec = els_p4[i] + els_p4[j];
    return vec.mass();
  } else {
    cout << "ERROR: mee tries to calculate Z mass from same charge leptons" << endl;
    return 999999.;
  }
}

float mmm(int i, int j){
  if ( mus_charge[i] * mus_charge[j] < 0 ) {
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
      vec = mus_p4[i] + mus_p4[j];
    return vec.mass();
  } else {
    cout << "ERROR: mmm tries to calculate Z mass from same charge leptons" << endl;
    return 999999.;
  }
}

bool goodLeptonIsolated(int bucket, int first, int second, int third) {
  // check if all 3 leptons are good isolated leptons
  //

  bool result = false;
  if ( bucket == 0 ||
       bucket == 1 ||
       bucket == 2 ||
       bucket == 3 ) {
    result = goodElectronIsolated(first) && goodElectronIsolated(second) && goodElectronIsolated(third);
  } else if ( bucket == 4 ||
	      bucket == 5 ||
	      bucket == 6 ||
	      bucket == 7 ) {
    result = goodMuonIsolated(first) && goodMuonIsolated(second) && goodMuonIsolated(third);
  } else if ( bucket == 8 ||
	      bucket == 9 ||
	      bucket == 10 ||
	      bucket == 11 ||
	      bucket == 12 ||
	      bucket == 13 ) {
    result = goodMuonIsolated(first) && goodElectronIsolated(second) && goodElectronIsolated(third);
  } else if ( bucket == 14 ||
	      bucket == 15 ||
	      bucket == 16 ||
	      bucket == 17 ||
	      bucket == 18 ||
	      bucket == 19 ) {
    result = goodMuonIsolated(first) && goodMuonIsolated(second) && goodElectronIsolated(third);
  } else {
    cout << "ERROR: goodIsolatedLepton tried to use not existing bucket!" << endl;
  }

  return result;
}


float ptLowestPtLepton(int bucket, int first, int second, int third) {
  // determine pt of lowest pt lepton
  float array[3] = {0,0,0};

  if ( bucket == 0 ||
       bucket == 1 ||
       bucket == 2 ||
       bucket == 3 ) {
    array[0] = els_p4[first].pt();
    array[1] = els_p4[second].pt();
    array[2] = els_p4[third].pt();
  } else if ( bucket == 4 ||
	      bucket == 5 ||
	      bucket == 6 ||
	      bucket == 7 ) {
    array[0] = mus_p4[first].pt();
    array[1] = mus_p4[second].pt();
    array[2] = mus_p4[third].pt();
  } else if ( bucket == 8 ||
	      bucket == 9 ||
	      bucket == 10 ||
	      bucket == 11 ||
	      bucket == 12 ||
	      bucket == 13 ) {
    array[0] = mus_p4[first].pt();
    array[1] = els_p4[second].pt();
    array[2] = els_p4[third].pt();
  } else if ( bucket == 14 ||
	      bucket == 15 ||
	      bucket == 16 ||
	      bucket == 17 ||
	      bucket == 18 ||
	      bucket == 19 ) {
    array[0] = mus_p4[first].pt();
    array[1] = mus_p4[second].pt();
    array[2] = els_p4[third].pt();
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

bool passTriggerLeptonMinPtCut(int bucket, int first, int second, int third, float triggerLeptonMinPtCut) {
  // return true if at least one lepton in trilepton candidate passes pt > triggerLeptonMinPtCut

  float array[3] = {0,0,0};

  if ( bucket == 0 ||
       bucket == 1 ||
       bucket == 2 ||
       bucket == 3 ) {
    array[0] = els_p4[first].pt();
    array[1] = els_p4[second].pt();
    array[2] = els_p4[third].pt();
  } else if ( bucket == 4 ||
	      bucket == 5 ||
	      bucket == 6 ||
	      bucket == 7 ) {
    array[0] = mus_p4[first].pt();
    array[1] = mus_p4[second].pt();
    array[2] = mus_p4[third].pt();
  } else if ( bucket == 8 ||
	      bucket == 9 ||
	      bucket == 10 ||
	      bucket == 11 ||
	      bucket == 12 ||
	      bucket == 13 ) {
    array[0] = mus_p4[first].pt();
    array[1] = els_p4[second].pt();
    array[2] = els_p4[third].pt();
  } else if ( bucket == 14 ||
	      bucket == 15 ||
	      bucket == 16 ||
	      bucket == 17 ||
	      bucket == 18 ||
	      bucket == 19 ) {
    array[0] = mus_p4[first].pt();
    array[1] = mus_p4[second].pt();
    array[2] = els_p4[third].pt();
  } else {
    cout << "ERROR: ptLowestPtLepton tried to use not existing bucket!" << endl;
  }

  // sort array
  sort(array,array+3);

  return array[2] > triggerLeptonMinPtCut;

}
