#include <TFile.h>
#include "TVector3.h"                   
#include <cmath>
#include <vector>
using namespace std;

#include "EGammaMvaEleEstimator.h"

// cms2
#include "../CORE/CMS2.h"
#include "../CORE/electronSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/eventSelections.h"

double electron_d0PV_wwV1_local(unsigned int index) { 
    if ( cms2.vtxs_sumpt().empty() ) return 9999.;
    //double sumPtMax = -1;
    //int iMax = -1;
    //for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
    //    if (!isGoodVertex(i)) continue; 
    //    if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
    //        iMax = i;
    //        sumPtMax = cms2.vtxs_sumpt().at(i);
    //    }
   // }
    //if (iMax<0) return 9999.;
    int iMax = 0; // try the first vertex
	double dxyPV=0;
	if(cms2.els_gsftrkidx()[index]>=0)
	{
    	dxyPV = cms2.els_d0()[index]-
        	cms2.vtxs_position()[iMax].x()*sin(cms2.gsftrks_p4()[cms2.els_gsftrkidx()[index]].phi())+
        	cms2.vtxs_position()[iMax].y()*cos(cms2.gsftrks_p4()[cms2.els_gsftrkidx()[index]].phi());
	}
	else 
	{
    	dxyPV = cms2.els_d0()[index]-
        	cms2.vtxs_position()[iMax].x()*sin(cms2.els_trk_p4()[index].phi())+
        	cms2.vtxs_position()[iMax].y()*cos(cms2.els_trk_p4()[index].phi());
	}

    return dxyPV;
}


//--------------------------------------------------------------------------------------------------
EGammaMvaEleEstimator::EGammaMvaEleEstimator() :
fMethodname("BDTG method"),
fisInitialized(kFALSE),
fPrintMVADebug(kFALSE),
fMVAType(kTrig),
fUseBinnedVersion(kTRUE),
fNMVABins(0)
{
  // Constructor.  
}

//--------------------------------------------------------------------------------------------------
EGammaMvaEleEstimator::~EGammaMvaEleEstimator()
{
  for (unsigned int i=0;i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void EGammaMvaEleEstimator::initialize( std::string methodName,
                                       	std::string weightsfile,
                                       	EGammaMvaEleEstimator::MVAType type)
{
  
  std::vector<std::string> tempWeightFileVector;
  tempWeightFileVector.push_back(weightsfile);
  initialize(methodName,type,kFALSE,tempWeightFileVector);
}


//--------------------------------------------------------------------------------------------------
void EGammaMvaEleEstimator::initialize( std::string methodName,
                                       	EGammaMvaEleEstimator::MVAType type,
                                       	Bool_t useBinnedVersion,
				       					std::vector<std::string> weightsfiles
  ) {

  //clean up first
  for (unsigned int i=0;i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
  fTMVAReader.clear();

  //initialize
  fisInitialized = kTRUE;
  fMVAType = type;
  fMethodname = methodName;
  fUseBinnedVersion = useBinnedVersion;

  //Define expected number of bins
  UInt_t ExpectedNBins = 0;
  if (!fUseBinnedVersion) {
    ExpectedNBins = 1;
  } else if (type == kTrig) {
    ExpectedNBins = 6;
  } else if (type == kNonTrig) {
    ExpectedNBins = 6;
  } else if (type == kIsoRings) {
    ExpectedNBins = 4;
  }
  fNMVABins = ExpectedNBins;
  
  //Check number of weight files given
  if (fNMVABins != weightsfiles.size() ) {
    std::cout << "Error: Expected Number of bins = " << fNMVABins << " does not equal to weightsfiles.size() = " 
              << weightsfiles.size() << std::endl; 
  }

  //Loop over all bins
  for (unsigned int i=0;i<fNMVABins; ++i) {
  
    TMVA::Reader *tmpTMVAReader = new TMVA::Reader( "!Color:!Silent:Error" );  
    tmpTMVAReader->SetVerbose(kTRUE);
  
    if (type == kTrig) {
      // Pure tracking variables
      tmpTMVAReader->AddVariable("fbrem",           &fMVAVar_fbrem);
      tmpTMVAReader->AddVariable("kfchi2",          &fMVAVar_kfchi2);
      tmpTMVAReader->AddVariable("kfhits",          &fMVAVar_kfhits);
      tmpTMVAReader->AddVariable("gsfchi2",         &fMVAVar_gsfchi2);

      // Geometrical matchings
      tmpTMVAReader->AddVariable("deta",            &fMVAVar_deta);
      tmpTMVAReader->AddVariable("dphi",            &fMVAVar_dphi);
      tmpTMVAReader->AddVariable("detacalo",        &fMVAVar_detacalo);
      // tmpTMVAReader->AddVariable("dphicalo",        &fMVAVar_dphicalo);   // Pruned but save in your ntuple. 
    
      // Pure ECAL -> shower shapes
      tmpTMVAReader->AddVariable("see",             &fMVAVar_see);
      tmpTMVAReader->AddVariable("spp",             &fMVAVar_spp);
      tmpTMVAReader->AddVariable("etawidth",        &fMVAVar_etawidth);
      tmpTMVAReader->AddVariable("phiwidth",        &fMVAVar_phiwidth);
      tmpTMVAReader->AddVariable("e1x5e5x5",        &fMVAVar_e1x5e5x5);
      tmpTMVAReader->AddVariable("R9",              &fMVAVar_R9);
      // tmpTMVAReader->AddVariable("nbrems",          &fMVAVar_nbrems); // Pruned but save in your ntuple. 

      // Energy matching
      tmpTMVAReader->AddVariable("HoE",             &fMVAVar_HoE);
      tmpTMVAReader->AddVariable("EoP",             &fMVAVar_EoP); 
      tmpTMVAReader->AddVariable("IoEmIoP",         &fMVAVar_IoEmIoP);
      tmpTMVAReader->AddVariable("eleEoPout",       &fMVAVar_eleEoPout);
      //  tmpTMVAReader->AddVariable("EoPout",          &fMVAVar_EoPout); // Pruned but save in your ntuple.    
      if(i == 2 || i == 5) 
	tmpTMVAReader->AddVariable("PreShowerOverRaw",&fMVAVar_PreShowerOverRaw);
      
      if(!fUseBinnedVersion)
	tmpTMVAReader->AddVariable("PreShowerOverRaw",&fMVAVar_PreShowerOverRaw);

      // IP
      tmpTMVAReader->AddVariable("d0",              &fMVAVar_d0);
      tmpTMVAReader->AddVariable("ip3d",            &fMVAVar_ip3d);
    
      tmpTMVAReader->AddSpectator("eta",            &fMVAVar_eta);
      tmpTMVAReader->AddSpectator("pt",             &fMVAVar_pt);
    }
  
    if (type == kNonTrig) {
      // Pure tracking variables
      tmpTMVAReader->AddVariable("fbrem",           &fMVAVar_fbrem);
      tmpTMVAReader->AddVariable("kfchi2",          &fMVAVar_kfchi2);
      tmpTMVAReader->AddVariable("kfhits",          &fMVAVar_kfhits);
      tmpTMVAReader->AddVariable("gsfchi2",         &fMVAVar_gsfchi2);

      // Geometrical matchings
      tmpTMVAReader->AddVariable("deta",            &fMVAVar_deta);
      tmpTMVAReader->AddVariable("dphi",            &fMVAVar_dphi);
      tmpTMVAReader->AddVariable("detacalo",        &fMVAVar_detacalo);
      // tmpTMVAReader->AddVariable("dphicalo",        &fMVAVar_dphicalo);   // Pruned but save in your ntuple. 
    
      // Pure ECAL -> shower shapes
      tmpTMVAReader->AddVariable("see",             &fMVAVar_see);
      tmpTMVAReader->AddVariable("spp",             &fMVAVar_spp);
      tmpTMVAReader->AddVariable("etawidth",        &fMVAVar_etawidth);
      tmpTMVAReader->AddVariable("phiwidth",        &fMVAVar_phiwidth);
      tmpTMVAReader->AddVariable("e1x5e5x5",        &fMVAVar_e1x5e5x5);
      tmpTMVAReader->AddVariable("R9",              &fMVAVar_R9);
      // tmpTMVAReader->AddVariable("nbrems",          &fMVAVar_nbrems); // Pruned but save in your ntuple. 

      // Energy matching
      tmpTMVAReader->AddVariable("HoE",             &fMVAVar_HoE);
      tmpTMVAReader->AddVariable("EoP",             &fMVAVar_EoP); 
      tmpTMVAReader->AddVariable("IoEmIoP",         &fMVAVar_IoEmIoP);
      tmpTMVAReader->AddVariable("eleEoPout",       &fMVAVar_eleEoPout);
      //  tmpTMVAReader->AddVariable("EoPout",          &fMVAVar_EoPout); // Pruned but save in your ntuple. 
      if(i == 2 || i == 5) 
	tmpTMVAReader->AddVariable("PreShowerOverRaw",&fMVAVar_PreShowerOverRaw);
    
      if(!fUseBinnedVersion)
	tmpTMVAReader->AddVariable("PreShowerOverRaw",&fMVAVar_PreShowerOverRaw);

      tmpTMVAReader->AddSpectator("eta",            &fMVAVar_eta);
      tmpTMVAReader->AddSpectator("pt",             &fMVAVar_pt);
    }

    if (type == kIsoRings) {
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p0To0p1",         &fMVAVar_ChargedIso_DR0p0To0p1        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p1To0p2",         &fMVAVar_ChargedIso_DR0p1To0p2        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p2To0p3",         &fMVAVar_ChargedIso_DR0p2To0p3        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p3To0p4",         &fMVAVar_ChargedIso_DR0p3To0p4        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p4To0p5",         &fMVAVar_ChargedIso_DR0p4To0p5        );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p0To0p1",           &fMVAVar_GammaIso_DR0p0To0p1          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p1To0p2",           &fMVAVar_GammaIso_DR0p1To0p2          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p2To0p3",           &fMVAVar_GammaIso_DR0p2To0p3          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p3To0p4",           &fMVAVar_GammaIso_DR0p3To0p4          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p4To0p5",           &fMVAVar_GammaIso_DR0p4To0p5          );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p0To0p1",   &fMVAVar_NeutralHadronIso_DR0p0To0p1  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p1To0p2",   &fMVAVar_NeutralHadronIso_DR0p1To0p2  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p2To0p3",   &fMVAVar_NeutralHadronIso_DR0p2To0p3  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p3To0p4",   &fMVAVar_NeutralHadronIso_DR0p3To0p4  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p4To0p5",   &fMVAVar_NeutralHadronIso_DR0p4To0p5  );
      tmpTMVAReader->AddSpectator("eta",            &fMVAVar_eta);
      tmpTMVAReader->AddSpectator("pt",             &fMVAVar_pt);
    }
  
    tmpTMVAReader->BookMVA(fMethodname , weightsfiles[i]);
    std::cout << "MVABin " << i << " : MethodName = " << fMethodname 
              << " , type == " << type << " , "
              << "Load weights file : " << weightsfiles[i] 
              << std::endl;
    fTMVAReader.push_back(tmpTMVAReader);
  }
  std::cout << "Electron ID MVA Completed\n";

}


//--------------------------------------------------------------------------------------------------
UInt_t EGammaMvaEleEstimator::GetMVABin( double eta, double pt) const {
  
    //Default is to return the first bin
    unsigned int bin = 0;

    if (fMVAType == EGammaMvaEleEstimator::kIsoRings) {
      if (pt < 10 && fabs(eta) < 1.479) bin = 0;
      if (pt < 10 && fabs(eta) >= 1.479) bin = 1;
      if (pt >= 10 && fabs(eta) < 1.479) bin = 2;
      if (pt >= 10 && fabs(eta) >= 1.479) bin = 3;
    }

    if (fMVAType == EGammaMvaEleEstimator::kNonTrig ) {
      bin = 0;
      if (pt < 10 && fabs(eta) < 0.8) bin = 0;
      if (pt < 10 && fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) bin = 1;
      if (pt < 10 && fabs(eta) >= 1.479) bin = 2;
      if (pt >= 10 && fabs(eta) < 0.8) bin = 3;
      if (pt >= 10 && fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) bin = 4;
      if (pt >= 10 && fabs(eta) >= 1.479) bin = 5;
    }


    if (fMVAType == EGammaMvaEleEstimator::kTrig) {
      bin = 0;
      if (pt < 20 && fabs(eta) < 0.8) bin = 0;
      if (pt < 20 && fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) bin = 1;
      if (pt < 20 && fabs(eta) >= 1.479) bin = 2;
      if (pt >= 20 && fabs(eta) < 0.8) bin = 3;
      if (pt >= 20 && fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) bin = 4;
      if (pt >= 20 && fabs(eta) >= 1.479) bin = 5;
    }

 

    return bin;
}

Double_t EGammaMvaEleEstimator::mvaValue(Int_t ele, Bool_t printDebug) {

	Double_t mvavalue = -999.;

	Double_t fbrem 				=	cms2.els_fbrem()[ele]; 
	Double_t kfchi2				=	cms2.els_trkidx()[ele]>=0 ? cms2.trks_chi2()[cms2.els_trkidx()[ele]]/cms2.trks_ndof()[cms2.els_trkidx()[ele]] : 0.;
	Int_t    kfhits				= 	cms2.els_trkidx()[ele]>=0 ? cms2.trks_nlayers()[cms2.els_trkidx()[ele]] : -1;
	Double_t gsfchi2			= 	cms2.els_chi2()[ele] / cms2.els_ndof()[ele];
	Double_t deta				=	cms2.els_dEtaIn()[ele];
	Double_t dphi				=	cms2.els_dPhiIn()[ele]; 
	Double_t detacalo			= 	cms2.els_dEtaOut()[ele];
	Double_t see				= 	cms2.els_sigmaIEtaIEta()[ele];
	Double_t spp				=	cms2.els_sigmaIPhiIPhi()[ele]; // FIXME : check the case where it's 0 
	Double_t etawidth			=	cms2.els_etaSCwidth()[ele];
	Double_t phiwidth			= 	cms2.els_phiSCwidth()[ele];
	Double_t e1x5e5x5			=	cms2.els_e5x5()[ele]!=0. ? 1. - cms2.els_e1x5()[ele]/cms2.els_e5x5()[ele] : -1; 
	Double_t R9					= 	cms2.els_e3x3()[ele] / cms2.els_eSCRaw()[ele];
	Double_t HoE				=	cms2.els_hOverE()[ele];
	Double_t EoP				=	cms2.els_eOverPIn()[ele];
	//Double_t IoEmIoP			=	1./cms2.els_eSC()[ele] - 1./cms2.els_p4()[ele].P(); 
	Double_t IoEmIoP			=	1./cms2.els_ecalEnergy()[ele] - 1./cms2.els_p4()[ele].P(); // this is consistent with CMSSW 
	Double_t eleEoPout			=	cms2.els_eOverPOut()[ele];
	Double_t PreShowerOverRaw	=	cms2.els_eSCPresh()[ele] / cms2.els_eSCRaw()[ele];
	Double_t d0					=	electron_d0PV_wwV1_local(ele);
	const double gsfsign = ( (gsftrks_d0_pv(cms2.els_gsftrkidx().at(ele),0).first)   >=0 ) ? 1. : -1.;
	Double_t ip3d				=	cms2.els_ip3d()[ele]*gsfsign; 
	Double_t eta				= 	cms2.els_etaSC()[ele];
	Double_t pt					= 	cms2.els_p4()[ele].pt();


	mvavalue =  EGammaMvaEleEstimator::mvaValue(
					 fbrem,
					 kfchi2,
					 kfhits,
					 gsfchi2,
					 deta,
					 dphi,
					 detacalo,
					// dphicalo,
					 see,
					 spp,
					 etawidth,
					 phiwidth,
					 e1x5e5x5,
					 R9,
					//Int_t    nbrems,
					 HoE,
					 EoP,
					 IoEmIoP,
					 eleEoPout,
					 PreShowerOverRaw,
					// EoPout,
					 d0,
					 ip3d,
					 eta,
					 pt,
					 printDebug); 


	return  mvavalue;
}

//--------------------------------------------------------------------------------------------------
Double_t EGammaMvaEleEstimator::mvaValue(Double_t fbrem, 
					Double_t kfchi2,
					Int_t    kfhits,
					Double_t gsfchi2,
					Double_t deta,
					Double_t dphi,
					Double_t detacalo,
					//Double_t dphicalo,
					Double_t see,
					Double_t spp,
					Double_t etawidth,
					Double_t phiwidth,
					Double_t e1x5e5x5,
					Double_t R9,
					//Int_t    nbrems,
					Double_t HoE,
					Double_t EoP,
					Double_t IoEmIoP,
					Double_t eleEoPout,
					Double_t PreShowerOverRaw,
					//Double_t EoPout,
					Double_t d0,
					Double_t ip3d,
					Double_t eta,
					Double_t pt,
					Bool_t printDebug) {
  
  if (!fisInitialized) { 
    std::cout << "Error: EGammaMvaEleEstimator not properly initialized.\n"; 
    return -9999;
  }

  fMVAVar_fbrem           = fbrem; 
  fMVAVar_kfchi2          = kfchi2;
  fMVAVar_kfhits          = float(kfhits);   // BTD does not support int variables
  fMVAVar_gsfchi2         = gsfchi2;

  fMVAVar_deta            = deta;
  fMVAVar_dphi            = dphi;
  fMVAVar_detacalo        = detacalo;
  // fMVAVar_dphicalo        = dphicalo;


  fMVAVar_see             = see;
  fMVAVar_spp             = spp;
  fMVAVar_etawidth        = etawidth;
  fMVAVar_phiwidth        = phiwidth;
  fMVAVar_e1x5e5x5        = e1x5e5x5;
  fMVAVar_R9              = R9;
  //fMVAVar_nbrems          = float(nbrems);   // BTD does not support int variables


  fMVAVar_HoE             = HoE;
  fMVAVar_EoP             = EoP;
  fMVAVar_IoEmIoP         = IoEmIoP;
  fMVAVar_eleEoPout       = eleEoPout;
  fMVAVar_PreShowerOverRaw= PreShowerOverRaw;
  //fMVAVar_EoPout          = EoPout; 

  fMVAVar_d0              = d0;
  fMVAVar_ip3d            = ip3d;
  fMVAVar_eta             = eta;
  fMVAVar_pt              = pt;


  bindVariables();
  Double_t mva = -9999;  
  if (fUseBinnedVersion) {
    mva = fTMVAReader[GetMVABin(fMVAVar_eta,fMVAVar_pt)]->EvaluateMVA(fMethodname);
  } else {
    mva = fTMVAReader[0]->EvaluateMVA(fMethodname);
  }

  if(printDebug) {
    cout << " *** Inside the class fMethodname " << fMethodname << endl;
	cout << " bin " << GetMVABin(fMVAVar_eta,fMVAVar_pt);
	cout << " fbrem " <<  fMVAVar_fbrem  
      	 << " kfchi2 " << fMVAVar_kfchi2  
	 	 << " kfhits " << fMVAVar_kfhits  
	 	 << " gsfchi2 " << fMVAVar_gsfchi2  
	 	 << " deta " <<  fMVAVar_deta  
	 	 << " dphi " << fMVAVar_dphi  
      	 << " detacalo " << fMVAVar_detacalo  
      // << " dphicalo " << fMVAVar_dphicalo  
	 	 << " see " << fMVAVar_see  
	 << " spp " << fMVAVar_spp  
	 << " etawidth " << fMVAVar_etawidth  
	 << " phiwidth " << fMVAVar_phiwidth  
	 << " e1x5e5x5 " << fMVAVar_e1x5e5x5  
	 << " R9 " << fMVAVar_R9  
      // << " mynbrems " << fMVAVar_nbrems  
	 << " HoE " << fMVAVar_HoE  
	 << " EoP " << fMVAVar_EoP  
	 << " IoEmIoP " << fMVAVar_IoEmIoP  
	 << " eleEoPout " << fMVAVar_eleEoPout  
      //<< " EoPout " << fMVAVar_EoPout  
	 << " PreShowerOverRaw " << fMVAVar_PreShowerOverRaw  
	 << " d0 " << fMVAVar_d0  
	 << " ip3d " << fMVAVar_ip3d  
	 << " eta " << fMVAVar_eta  
	 << " pt " << fMVAVar_pt << endl;
    cout << " ### MVA " << mva << endl;
  }


  return mva;
}
//--------------------------------------------------------------------------------------------------
Double_t EGammaMvaEleEstimator::mvaValue(Double_t fbrem, 
					Double_t kfchi2,
					Int_t    kfhits,
					Double_t gsfchi2,
					Double_t deta,
					Double_t dphi,
					Double_t detacalo,
					//Double_t dphicalo,
					Double_t see,
					Double_t spp,
					Double_t etawidth,
					Double_t phiwidth,
					Double_t e1x5e5x5,
					Double_t R9,
					//Int_t    nbrems,
					Double_t HoE,
					Double_t EoP,
					Double_t IoEmIoP,
					Double_t eleEoPout,
					Double_t PreShowerOverRaw,
					//Double_t EoPout,
					Double_t eta,
					Double_t pt,
					Bool_t printDebug) {
  
  if (!fisInitialized) { 
    std::cout << "Error: EGammaMvaEleEstimator not properly initialized.\n"; 
    return -9999;
  }

  fMVAVar_fbrem           = fbrem; 
  fMVAVar_kfchi2          = kfchi2;
  fMVAVar_kfhits          = float(kfhits);   // BTD does not support int variables
  fMVAVar_gsfchi2         = gsfchi2;

  fMVAVar_deta            = deta;
  fMVAVar_dphi            = dphi;
  fMVAVar_detacalo        = detacalo;
  // fMVAVar_dphicalo        = dphicalo;


  fMVAVar_see             = see;
  fMVAVar_spp             = spp;
  fMVAVar_etawidth        = etawidth;
  fMVAVar_phiwidth        = phiwidth;
  fMVAVar_e1x5e5x5        = e1x5e5x5;
  fMVAVar_R9              = R9;
  //fMVAVar_nbrems          = float(nbrems);   // BTD does not support int variables


  fMVAVar_HoE             = HoE;
  fMVAVar_EoP             = EoP;
  fMVAVar_IoEmIoP         = IoEmIoP;
  fMVAVar_eleEoPout       = eleEoPout;
  fMVAVar_PreShowerOverRaw= PreShowerOverRaw;
  //fMVAVar_EoPout          = EoPout; 

  fMVAVar_eta             = eta;
  fMVAVar_pt              = pt;


  bindVariables();
  Double_t mva = -9999;  
  if (fUseBinnedVersion) {
    mva = fTMVAReader[GetMVABin(fMVAVar_eta,fMVAVar_pt)]->EvaluateMVA(fMethodname);
  } else {
    mva = fTMVAReader[0]->EvaluateMVA(fMethodname);
  }



  if(printDebug) {
    cout << " *** Inside the class fMethodname " << fMethodname << endl;
    cout << " bin " <<  GetMVABin(fMVAVar_eta,fMVAVar_pt);  
    cout << " fbrem " <<  fMVAVar_fbrem  
      	 << " kfchi2 " << fMVAVar_kfchi2  
	 << " mykfhits " << fMVAVar_kfhits  
	 << " gsfchi2 " << fMVAVar_gsfchi2  
	 << " deta " <<  fMVAVar_deta  
	 << " dphi " << fMVAVar_dphi  
      	 << " detacalo " << fMVAVar_detacalo  
      // << " dphicalo " << fMVAVar_dphicalo  
	 << " see " << fMVAVar_see  
	 << " spp " << fMVAVar_spp  
	 << " etawidth " << fMVAVar_etawidth  
	 << " phiwidth " << fMVAVar_phiwidth  
	 << " e1x5e5x5 " << fMVAVar_e1x5e5x5  
	 << " R9 " << fMVAVar_R9  
      // << " mynbrems " << fMVAVar_nbrems  
	 << " HoE " << fMVAVar_HoE  
	 << " EoP " << fMVAVar_EoP  
	 << " IoEmIoP " << fMVAVar_IoEmIoP  
	 << " eleEoPout " << fMVAVar_eleEoPout  
      //<< " EoPout " << fMVAVar_EoPout  
	 << " PreShowerOverRaw " << fMVAVar_PreShowerOverRaw  
	 << " eta " << fMVAVar_eta  
	 << " pt " << fMVAVar_pt << endl;
    cout << " ### MVA " << mva << endl;
  }


  return mva;
}

void EGammaMvaEleEstimator::bindVariables() {

  // this binding is needed for variables that sometime diverge. 


  if(fMVAVar_fbrem < -1.)
    fMVAVar_fbrem = -1.;	
  
  fMVAVar_deta = fabs(fMVAVar_deta);
  if(fMVAVar_deta > 0.06)
    fMVAVar_deta = 0.06;
  
  
  fMVAVar_dphi = fabs(fMVAVar_dphi);
  if(fMVAVar_dphi > 0.6)
    fMVAVar_dphi = 0.6;
  
  
//   if(fMVAVar_EoPout > 20.)
//     fMVAVar_EoPout = 20.;
  
  if(fMVAVar_EoP > 20.)
    fMVAVar_EoP = 20.;
  
  if(fMVAVar_eleEoPout > 20.)
    fMVAVar_eleEoPout = 20.;
  
  
  fMVAVar_detacalo = fabs(fMVAVar_detacalo);
  if(fMVAVar_detacalo > 0.2)
    fMVAVar_detacalo = 0.2;
  
  
//   fMVAVar_dphicalo = fabs(fMVAVar_dphicalo);
//   if(fMVAVar_dphicalo > 0.4)
//     fMVAVar_dphicalo = 0.4;
  
  
  if(fMVAVar_e1x5e5x5 < -1.)
    fMVAVar_e1x5e5x5 = -1;
  
  if(fMVAVar_e1x5e5x5 > 2.)
    fMVAVar_e1x5e5x5 = 2.; 
  
  
  
  if(fMVAVar_R9 > 5)
    fMVAVar_R9 = 5;
  
  if(fMVAVar_gsfchi2 > 200.)
    fMVAVar_gsfchi2 = 200;
  
  
  if(fMVAVar_kfchi2 > 10.)
    fMVAVar_kfchi2 = 10.;
  
  
  // Needed for a bug in CMSSW_420, fixed in more recent CMSSW versions
  if(std::isnan(fMVAVar_spp))
    fMVAVar_spp = 0.;	
  
  
  return;
}








