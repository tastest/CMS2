#include "DYMVA.h"

#include "analysisObjects.h"
#include "analysisTools.h"
#include "analysisSelections.h"
#include "../CORE/metSelections.h"

// need to set proper weight files 
#include "files/TMVA_BDTG_0j_mll12.class.C"
#include "files/TMVA_BDTG_1j_mll12.class.C"

#include "../CORE/CMS2.h"

using namespace std;

float recoilvar(float met, float metPhi, LorentzVector* dilep)
{
  float px = met*cos(metPhi) + dilep->px();       
  float py = met*sin(metPhi) + dilep->py();
  return sqrt(px*px+py*py);
}

//###################
//# main function
//###################
float DYMVA(unsigned int ihyp, unsigned int njets, std::vector<JetPair> jets) {

  std::vector<std::string> theInputVars;
  const char* inputVars[] = { "pmet", "pTrackMet","nvtx", "dilpt", "jet1pt", "metSig", "dPhiDiLepJet1", "dPhiDiLepMET", "dPhiJet1MET", "recoil", "mt" };
  for (int i=0;i<11;++i) theInputVars.push_back(inputVars[i]);
  dymva_0j_Zveto::ReadBDTG* rbdtgDy_0j = new dymva_0j_Zveto::ReadBDTG(theInputVars);
  dymva_1j_Zveto::ReadBDTG* rbdtgDy_1j = new dymva_1j_Zveto::ReadBDTG(theInputVars);
 
  //==========================================
  // Loop All Events
  //==========================================  
  //cout << smurfFDir + fileName << " has " << ch->GetEntries() << " entries; \n";

  LorentzVector*  dilep_	= 0;
  float pmet_ 				= 0.;
  float pTrackMet_ 			= 0.;
  float mt_ 				= 0.;
  LorentzVector*  jet1_		= 0;
  float jet1pt_ 			= 0.;
  float dPhiDiLepJet1_ 		= -999.;
  float dPhiDiLepMET_ 		= -999.;
  float met_ 				= 0.;
  unsigned int nvtx_ 		= 0;
  float sumet_ 				= 0.; 
  float metPhi_ 			= 0.;
  unsigned int njets_ 		= 0;
  float dymva_				= -999.;
  float recoil_				= -999.; 
  float dPhiJet1MET_		= -999.;
 
  // input variables
  dilep_ = &cms2.hyp_p4().at(ihyp);
  if(jets.size()>0) {
  	jet1_ = &jets.at(0).first;
	jet1pt_ = jet1_->pt(); 
  } 
  met_ 	= cms2.evt_pfmet();
  metPhi_ = cms2.evt_pfmetPhi();
  
  pmet_ 			= projectedMet(ihyp, met_, metPhi_); 
  metStruct trkMET 	= trackerMET(ihyp,0.1);  
  pTrackMet_ 		= projectedMet(ihyp, trkMET.met, trkMET.metphi);
  nvtx_ 			= nGoodVertex(); 
  sumet_ 			= cms2.evt_pfsumet(); 
  if(jets.size()>0) 
  dPhiDiLepJet1_  	= fabs(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_p4().at(ihyp), *jet1_));
  dPhiDiLepMET_ 	= acos(cos(cms2.hyp_p4().at(ihyp).phi()-metPhi_));
  mt_ 				= mt(dilep_->pt(),met_,dPhiDiLepMET_);  
  njets_ 			= njets;

  LorentzVector metlv( met_*cos(metPhi_), met_*sin(metPhi_), 0, met_ );
  assert((metlv.phi()-metPhi_)<0.001);
  assert((metlv.pt()-met_)<0.001);

  recoil_ = recoilvar(met_, metPhi_, dilep_);
  if(jets.size()>0) 
  dPhiJet1MET_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(*jet1_,metlv));

  //const char* inputVars[] = { "pmet","pTrackMet","nvtx","dilpt","jet1pt","metSig","dPhiDiLepJet1","dPhiDiLepMET","dPhiJet1MET","recoil","mt" };
  std::vector<double> theInputVals;
  const double inputVals[] = { pmet_,pTrackMet_,nvtx_,dilep_->pt(),
				 std::max((float)15.,jet1pt_),
				 met_/sqrt(sumet_),
				 (jet1pt_<15. ? -0.1 : dPhiDiLepJet1_ ),
				 dPhiDiLepMET_,
				 (jet1pt_<15. ? -0.1 : dPhiJet1MET_ ),
				 recoil_,
				 mt_};


  for (int i=0;i<11;++i) theInputVals.push_back(inputVals[i]);
    
  if (njets_==0) dymva_ = rbdtgDy_0j->GetMvaValue(theInputVals);
  else if (njets_==1) dymva_ = rbdtgDy_1j->GetMvaValue(theInputVals);
  else dymva_= -999.;

  // debug
  if(0) 
  {
	  std::cout << cms2.evt_event() << " : " <<  cms2.hyp_type().at(ihyp) << " : ";
	  for(int i=0; i<11; i++) { 
		  std::cout << inputVals[i] << " : ";
	  } 
	  std::cout << njets_ << " :::: " << dymva_ << std::endl;
  }

  delete rbdtgDy_0j;
  delete rbdtgDy_1j;

  return dymva_;

}  


