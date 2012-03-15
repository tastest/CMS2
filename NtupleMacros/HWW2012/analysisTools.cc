
#include "analysisTools.h"

#include <iostream>

//
// Process utilities
//

bool isDYee() {
    if (getDrellYanType() == 0) return true;
    return false;
}
bool isDYmm() {
    if (getDrellYanType() == 1) return true;
    return false;
}
bool isDYtt() {
    if (getDrellYanType() == 2) return true;
    return false;
}

bool isWW() {
    if (getVVType() == 0) return true;
    return false;
}

bool isWZ() {
    if (getVVType() == 1) return true;
    return false;
}

bool isZZ() {
    if (getVVType() == 2) return true;
    return false;
}

//-------------------------------------------------
// Auxiliary function to scan the doc line and 
// identify DY-> ee vs mm vs tt
//-------------------------------------------------
unsigned int getDrellYanType() {
    bool foundEP = false;
    bool foundEM = false;
    bool foundMP = false;
    bool foundMM = false;
    bool foundTP = false;
    bool foundTM = false;
    for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
        if ( cms2.genps_id_mother().at(i) == 23 ){
            switch ( TMath::Abs(cms2.genps_id().at(i)) ){
                case 11:
                    return 0;
                    break;
                case 13:
                    return 1;
                    break;
                case 15:
                    return 2;
                    break;
                default:
                    break;
            }
        }
        switch ( cms2.genps_id().at(i) ){
            case 11:
                foundEM = true;
                break;
            case -11:
                foundEP = true;
                break;
            case 13:
                foundMM = true;
                break;
            case -13:
                foundMP = true;
                break;
            case 15:
                foundTM = true;
                break;
            case -15:
                foundTP = true;
                break;
            default:
                break;
        }
    }

    if ( foundEP && foundEM ) return 0;  //DY->ee
    if ( foundMP && foundMM ) return 1;  //DY->mm
    if ( foundTP && foundTM ) return 2;  //DY->tautau
    std::cout << "Does not look like a DY event" << std::endl;
    return 999;
}

unsigned int getVVType() {
    // types:
    //   0 - WW
    //   1 - WZ
    //   2 - ZZ
    unsigned int nZ(0);
    unsigned int nW(0);
    std::vector<std::vector<int> > leptons;
    std::vector<int> mothers;

    bool verbose = false;

    for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
        int pid = cms2.genps_id().at(i);
        int mid = cms2.genps_id_mother().at(i);
        if ( verbose ) std::cout << "Gen particle id: " << pid << ",\t mother id: " << mid <<std::endl;
        if ( abs(pid)<11 || abs(pid)>16 ) continue;
        if ( mid == 23 ) ++nZ;
        if ( abs(mid) == 24 ) ++nW;
        // now we need to really understand the pattern.
        unsigned int mIndex = 0;
        while ( mIndex < mothers.size() && mid != mothers[mIndex] ) ++mIndex;
        if ( mIndex == mothers.size() ) {
            mothers.push_back(mid);
            leptons.push_back(std::vector<int>());
        }
        leptons[mIndex].push_back(pid);
        if (mothers.size()>3){
            if (verbose) std::cout << "WARNING: failed to identify event (too many mothers)" << std::endl;
            return 999;
        }
    }

    if ( nZ == 4 ) {
        if ( verbose ) std::cout << "Event type ZZ" << std::endl;
        return 2;
    }
    if ( nW == 4 ) {
        if ( verbose ) std::cout << "Event type WW" << std::endl;
        return 0;
    }
    if ( nW == 2 && nZ == 2 ) {
        if ( verbose ) std::cout << "Event type WZ" << std::endl;
        return 1;
    }
    unsigned int nNus(0);
    for ( unsigned int i=0; i<mothers.size(); ++i ){
        nNus += leptons[i].size();
    }
    if ( mothers.size() < 3 && nNus == 4){
        for ( unsigned int i=0; i<mothers.size(); ++i ){
            if ( mothers[i] != 23 && abs(mothers[i]) != 24 ){
                if( leptons[i].size() != 2 && leptons[i].size() != 4){
                    if (verbose) std::cout << "WARNING: failed to identify event (unexpected number of daughters)" << std::endl;
                    if (verbose) std::cout << "\tnumber of daughters for first mother: " <<  leptons[0].size() << std::endl;
                    if (verbose) std::cout << "\tnumber of daughters for second mother: " <<  leptons[1].size() << std::endl;
                    return 999;
                }
                if ( abs(leptons[i][0]) == abs(leptons[i][1]) )
                    nZ += 2;
                else
                    nW += 2;
                if ( leptons[i].size()==4 ){
                    // now it's a wild guess, it's fraction should be small
                    if ( abs(leptons[i][2]) == abs(leptons[i][3]) )
                        nZ += 2;
                    else
                        nW += 2;
                }
            }
        }
    } else {
        // here be dragons

        // if we have 2 leptons and 3 neutrinos and they all of the same
        // generation, we assume it's ZZ (can be WZ also), but if
        // neutrinos are from different generations, than we conclude it's
        // WZ. 

        std::set<int> nus;
        for ( unsigned int i=0; i<mothers.size(); ++i )
            for ( unsigned int j=0; j<leptons[i].size(); ++j )
                if ( abs(leptons[i][j]) == 12 ||
                        abs(leptons[i][j]) == 14 ||
                        abs(leptons[i][j]) == 16 )
                    nus.insert(abs(leptons[i][j]));

        if ( nNus == 5 ){
            if ( nus.size() == 1 ) return 2;
            if ( nus.size() == 2 ) return 1;
        }

        if ( verbose ) std::cout << "WARNING: failed to identify event" << std::endl;
        return 999;
    }

    if ( nZ+nW != 4 ){
        if (verbose) std::cout << "WARNING: failed to identify event (wrong number of bosons)" << std::endl;
        if (verbose) std::cout << "\tfirst mother id: " << mothers[0] << std::endl;
        if (verbose) std::cout << "\tsecond mother id: " << mothers[1] << std::endl;
        if (verbose) std::cout << "\tnumber of daughters for first mother: " << leptons[0].size() << std::endl;
        if (verbose) std::cout << "\tnumber of daughters for second mother: " << leptons[1].size() << std::endl;
        if (verbose) std::cout << "\tnumber of Zs: " << nZ << std::endl;
        if (verbose) std::cout << "\tnumber of Ws: " << nW << std::endl;
        return 999;
    }

    if ( nZ == 4 ) {
        if ( verbose ) std::cout << "Event type ZZ" << std::endl;
        return 2;
    }
    if ( nW == 4 ) {
        if ( verbose ) std::cout << "Event type WW" << std::endl;
        return 0;
    }
    // this covers screws in logic, i.e. most hard to identify events end up being WZ
    if ( verbose ) std::cout << "Event type WZ (can be wrong)" << std::endl;
    return 1;
}

// filter events by process
bool filterByProcess( enum SmurfTree::DataType sample ) {
    switch (sample) {
        case SmurfTree::dyee:
            return isDYee();
        case SmurfTree::dymm:
            return isDYmm();
        case SmurfTree::dytt:
            return isDYtt();
        case SmurfTree::qqww:
            return isWW();
        case SmurfTree::wz:
            return isWZ();
        case SmurfTree::zz:
            return isZZ();
        default:
            return true;
    }
}

bool isIdentified( enum SmurfTree::DataType sample ) {
    switch (sample) {
        case SmurfTree::dyee:
        case SmurfTree::dymm:
        case SmurfTree::dytt:
            return getDrellYanType()!=999;
        case SmurfTree::qqww:
        case SmurfTree::wz:
        case SmurfTree::zz:
            return getVVType()!=999;
        default:
            return true;
    }
}

//
// event utilities
//

EventIdentifier::EventIdentifier(CMS2& cms2, bool isData){
    data = isData;
    run = cms2.evt_run();
    event = cms2.evt_event();
    lumi = cms2.evt_lumiBlock();
    trks_d0 = cms2.trks_d0().empty() ? 0 : cms2.trks_d0().front();
}

bool EventIdentifier::operator < (const EventIdentifier &other) const
{
    if (run != other.run)  return run < other.run;
    if (event != other.event) return event < other.event;
    if (data) return false;
    if (fabs(trks_d0 - other.trks_d0) > 1e-6 * trks_d0) return trks_d0 < other.trks_d0;
    return false;
}

bool EventIdentifier::operator == (const EventIdentifier &other) const
{
    if (run != other.run) return false;
    if (event != other.event) return false;
    if (data) return true;
    if (fabs(trks_d0 - other.trks_d0) > 1e-6 * trks_d0)  return false;
    return true;
}

bool is_duplicate (const EventIdentifier &id)
{
    std::pair<std::set<EventIdentifier>::const_iterator, bool> ret =
        already_seen.insert(id);
    return !ret.second;
}

//
// cut utilities
//

bool CheckCutsNM1(wwcuts_t apply, wwcuts_t remove, wwcuts_t passed)
{
  if ((passed & (apply & (~remove))) == (apply & (~remove))) return true;
  return false;
}

bool CheckCuts(wwcuts_t apply, wwcuts_t passed)
{
  if ((apply & passed) == apply) return true;
  return false;
}

//
// kinematics and weights
//

float getHiggsPt() {
  for (unsigned int i=0; i<cms2.genps_id().size(); ++i) {
    if (cms2.genps_status().at(i) == 3 && cms2.genps_id().at(i) == 25) {
      return cms2.genps_p4().at(i).pt();
    }
  }
  return -1.;
}

float getHiggsPtWeight(float pt, int higgsMass, TH1D *HiggsPtKFactor) {
  if (pt<0) return 1.;
  if (HiggsPtKFactor!=0 && TString(HiggsPtKFactor->GetName()).Contains(Form("%i",higgsMass)))
    return HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(pt));
  TFile *fHiggsPtKFactorFile = TFile::Open("./files/ggHWW_KFactors_PowhegToHQT_WithAdditionalMassPoints.root");
  assert(fHiggsPtKFactorFile);
  char kfactorHistName[100];
  sprintf(kfactorHistName, "KFactor_PowhegToHQT_mH%d", higgsMass);
  HiggsPtKFactor = (TH1D*)(fHiggsPtKFactorFile->Get(kfactorHistName));
  if (HiggsPtKFactor) HiggsPtKFactor->SetDirectory(0);
  assert(HiggsPtKFactor);
  fHiggsPtKFactorFile->Close();
  delete fHiggsPtKFactorFile;
  //cout << "Using new hist for higgs pT k-factors: " << HiggsPtKFactor->GetName() << endl;
  return HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(pt));
}

float getZPt() {
  for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
    if ( cms2.genps_id().at(i) == 23 ){
      return cms2.genps_p4().at(i).pt();
    }
  }
  return -999.;
}

float getZRapidity() {
  for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
    if ( cms2.genps_id().at(i) == 23 ){
      return cms2.genps_p4().at(i).Rapidity();
    }
  }
  return -999.;
}

float getZMass() {
  for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
    if ( cms2.genps_id().at(i) == 23 ){
      return cms2.genps_p4().at(i).mass();
    }
  }
  return -999.;
}

float getDYNNLOWeight(float pt, float rap, float mass, vector<TH2D*> fDYNNLOKFactorHists) {

  float DYNNLOKFactor = 1.;
  //Find Mass Bin
  const UInt_t nMassBins = 41;
  const Double_t massBinLimits[nMassBins+1] = {15,   20,  25,  30,  35,  40,  45,  50,  55,  60,
                           64,   68,  72,  76,  81,  86,  91,  96, 101, 106,
                           110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
                           200, 220, 243, 273, 320, 380, 440, 510, 600, 1000,
                           1500}; //41 bins
  //const UInt_t nMassBins = 13;
  //const Double_t massBinLimits[nMassBins+1] = {0,20,30,40,50,60,76,86,96,106,120,150,200,600}; // 13 bins
  Int_t massBinIndex = -1 ;
  for(UInt_t binIndex=0; binIndex < nMassBins; ++binIndex){
    if( mass >= massBinLimits[binIndex] && mass < massBinLimits[binIndex+1]) {
      massBinIndex = binIndex;
      break;
    }
  }
  //Found the mass bin
  if (massBinIndex >= 0 && massBinIndex < Int_t(nMassBins) ) {
    UInt_t ptBin = fDYNNLOKFactorHists[massBinIndex]->GetXaxis()->FindFixBin(pt);
    UInt_t yBin = fDYNNLOKFactorHists[massBinIndex]->GetYaxis()->FindFixBin(rap);

    if(Int_t(ptBin) == fDYNNLOKFactorHists[massBinIndex]->GetNbinsX() + 1)
      ptBin = fDYNNLOKFactorHists[massBinIndex]->GetNbinsX();
    if(ptBin == 0)
      ptBin = 1;
    if(Int_t(yBin) == fDYNNLOKFactorHists[massBinIndex]->GetNbinsY() + 1)
      yBin = fDYNNLOKFactorHists[massBinIndex]->GetNbinsY();
    if(yBin == 0)
      yBin = 1;
    DYNNLOKFactor = fDYNNLOKFactorHists[massBinIndex]->GetBinContent( ptBin, yBin);
  }
  return DYNNLOKFactor;

}

double mt(double pt1, double pt2, double dphi){
  return 2*sqrt(pt1*pt2)*fabs(sin(dphi/2));
}




