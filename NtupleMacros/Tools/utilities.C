#include "TH1F.h"
#include <algorithm>
#include "../Tools/selections.h"
#include "Math/VectorUtil.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
#include "TMath.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

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
  if ( cms2.els_charge()[i] * cms2.els_charge()[j] < 0 ) {
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
      vec = cms2.els_p4()[i] + cms2.els_p4()[j];
    return vec.mass();
  } else {
    cout << "ERROR: mee tries to calculate Z mass from same charge leptons" << endl;
    return 999999.;
  }
}

float mmm(int i, int j){
  if ( cms2.mus_charge()[i] * cms2.mus_charge()[j] < 0 ) {
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 
      vec = cms2.mus_p4()[i] + cms2.mus_p4()[j];
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

bool passTriggerLeptonMinPtCut(int bucket, int first, int second, int third, float triggerLeptonMinPtCut) {
  // return true if at least one lepton in trilepton candidate passes pt > triggerLeptonMinPtCut

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
    if( abs(mz2-zmass) < abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
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
    if( abs(mz2-zmass) < abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
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
    if( abs(mz2-zmass) < abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
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
    if( abs(mz2-zmass) < abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
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
    if( abs(mz2-zmass) < abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
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
    if( abs(mz2-zmass) < abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
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
    if( abs(mz2-zmass) < abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
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
    if( abs(mz2-zmass) < abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
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
    if( abs(mz2-zmass) < abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
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
    if( abs(mz2-zmass) < abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
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
    if( abs(mz2-zmass) < abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
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
    if( abs(mz2-zmass) < abs(mz1-zmass) ){ //swap the masses to always make mz1 the one thats closer to Z
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

bool passMETAllCut(int bucket, float metAll, float electronMETAllCut, float muonMETAllCut) {
  // cut on MET all
  // use electronMETAllCut if lepton not belonging to the primary Z is an electron
  // use muonMETAllCut if lepton not belonging to the primary Z is an muon

  bool result = false;
  
  if ( bucket == 17 ||
       bucket == 18 ||
       bucket == 7 ||
       bucket == 6 ) {
    result = metAll > electronMETAllCut;
  } else if ( bucket == 14 ||
	      bucket == 11 ||
	      bucket == 1 ||
	      bucket == 2 ) {
    result = metAll > muonMETAllCut;
  }

  return result;
}

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
	      if ( abs(zmass-mass) < abs(zmass-result) ) {
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
	      if ( abs(zmass-mass) < abs(zmass-result) ) {
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
	      if ( abs(zmass-mass) < abs(zmass-result) ) {
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
	      if ( abs(zmass-mass) < abs(zmass-result) ) {
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

bool testJetForElectrons(const LorentzVector& jetP4, const LorentzVector& elP4) {
  
  
  bool matched = false;
  float elphi  = elP4.Phi();
  float jetphi = jetP4.Phi();
   
  float eleta  = elP4.Eta();
  float jeteta = jetP4.Eta();
   
  float dphi = elphi - jetphi;
  float deta = eleta - jeteta;
  if(fabs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - fabs(dphi);
   
  double dR = sqrt(dphi*dphi + deta*deta);
  if (dR < 0.4) 
    matched = true;
  
  return !matched;
}

std::vector<LorentzVector> correctJetsForElectrons(int bucket, int first, int second, int third) {

  std::vector<LorentzVector> result;

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

  return result;

}
