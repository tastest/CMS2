// -*- C++ -*-
//
// Package:    SParmMaker
// Class:      SParmMaker
// 
/**\class SParmMaker SParmMaker.cc CMS2/NtupleMaker/src/SParmMaker.cc

   Description: copy SUSY mSUGRA parameters into the EDM event tree

   Implementation:
   - extract and fill variables
*/
//
// Original Ben Hooberman
// Created:  Wed Mar  24 12:23:38 CDT 2010
// 
//
//


// system include files
#include <memory>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "CMS2/NtupleMaker/interface/SParmMaker.h"

#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"

//
// class declaration
//

//
// constructors and destructor
//

SParmMaker::SParmMaker(const edm::ParameterSet& iConfig) {

    aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    // product of this EDProducer
    produces<TString> (branchprefix+"comment"      ).setBranchAlias(aliasprefix_+"_comment"      );
    produces<float>   (branchprefix+"m0"           ).setBranchAlias(aliasprefix_+"_m0"           );
    produces<float>   (branchprefix+"m12"          ).setBranchAlias(aliasprefix_+"_m12"           );
    produces<float>   (branchprefix+"tanb"         ).setBranchAlias(aliasprefix_+"_tanb"           );
    produces<float>   (branchprefix+"a0"           ).setBranchAlias(aliasprefix_+"_a0"           );
    produces<float>   (branchprefix+"mu"           ).setBranchAlias(aliasprefix_+"_mu"           );
    produces<int>     (branchprefix+"subProcessId" ).setBranchAlias(aliasprefix_+"_subProcessId" );
    produces<float>   (branchprefix+"weight"       ).setBranchAlias(aliasprefix_+"_weight"       );
    produces<float>   (branchprefix+"pdfWeight1"   ).setBranchAlias(aliasprefix_+"_pdfWeight1"   );
    produces<float>   (branchprefix+"pdfWeight2"   ).setBranchAlias(aliasprefix_+"_pdfWeight2"   );
    produces<float>   (branchprefix+"pdfScale"     ).setBranchAlias(aliasprefix_+"_pdfScale"     );
  
  
    // parameters from configuration
    sparm_inputTag = iConfig.getParameter<edm::InputTag>("sparm_inputTag");
}

SParmMaker::~SParmMaker() {}

//
// member functions
//

// ------------ method called to produce the data  ------------
void SParmMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   
    std::auto_ptr<TString> sparm_comment      (new TString);
    std::auto_ptr<int>     sparm_subProcessId (new int);
    std::auto_ptr<float>   sparm_m0           (new float);
    std::auto_ptr<float>   sparm_m12           (new float);
    std::auto_ptr<float>   sparm_tanb           (new float);
    std::auto_ptr<float>   sparm_a0           (new float);
    std::auto_ptr<float>   sparm_mu           (new float);
    std::auto_ptr<float>   sparm_weight       (new float);
    std::auto_ptr<float>   sparm_pdfWeight1   (new float);
    std::auto_ptr<float>   sparm_pdfWeight2   (new float);  
    std::auto_ptr<float>   sparm_pdfScale     (new float);

    edm::Handle<LHEEventProduct> sparm_handle;  

    iEvent.getByLabel(sparm_inputTag, sparm_handle); 

    if( !sparm_handle.isValid() ){
        *sparm_m0		= -9999.;
        *sparm_m12      = -9999.;
        *sparm_tanb     = -9999.;
        *sparm_a0       = -9999.;
        *sparm_mu       = -9999.;
    }

    for (std::vector<std::string>::const_iterator it = sparm_handle->comments_begin(); it != sparm_handle->comments_end(); it++) {      
        TString model_params(*it);
		//std::cout << "params = " << model_params;// << std::endl;
        if (!model_params.Contains("model"))
            continue;
        *sparm_comment = TString(*it);
        TObjArray* tokens = model_params.Tokenize("_");
		const int npar = tokens->GetEntries() - 1; //the crap before the first _ is not a param

		*sparm_m0      = ((TObjString*)tokens->At(1))->GetString().Atof();
		*sparm_m12     = ((TObjString*)tokens->At(2))->GetString().Atof();
		*sparm_tanb    = ((TObjString*)tokens->At(3))->GetString().Atof();
		*sparm_a0      = ((TObjString*)tokens->At(4))->GetString().ReplaceAll("m", "").Atof(); // need to remove "m" if exists
		*sparm_mu      = ((TObjString*)tokens->At(5))->GetString().Atof(); 
		
		delete tokens;
    }

	std::cout << *sparm_m0 << "  " << *sparm_m12 << "  " << *sparm_tanb << "  " << *sparm_a0 << "  " << *sparm_mu << std::endl;

    // now, get info about this event
    const lhef::HEPEUP lhe_info = sparm_handle->hepeup();
    *sparm_subProcessId = lhe_info.IDPRUP;
    *sparm_weight       = lhe_info.XWGTUP;
    *sparm_pdfWeight1   = lhe_info.XPDWUP.first;
    *sparm_pdfWeight2   = lhe_info.XPDWUP.second;

    const gen::PdfInfo* pdf_info = sparm_handle->pdf();
    if (pdf_info != 0) {
        *sparm_pdfScale     = pdf_info->scalePDF;
    }
    else {
        *sparm_pdfScale = lhe_info.SCALUP;
    }

    // put containers into event
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    iEvent.put(sparm_m0           ,branchprefix+"m0"           );
    iEvent.put(sparm_m12          ,branchprefix+"m12"           );
    iEvent.put(sparm_tanb         ,branchprefix+"tanb"           );
    iEvent.put(sparm_a0           ,branchprefix+"a0"           );
    iEvent.put(sparm_mu           ,branchprefix+"mu"           );
    iEvent.put(sparm_weight       ,branchprefix+"weight"       );
    iEvent.put(sparm_pdfWeight1   ,branchprefix+"pdfWeight1"   );
    iEvent.put(sparm_pdfWeight2   ,branchprefix+"pdfWeight2"   );
    iEvent.put(sparm_pdfScale     ,branchprefix+"pdfScale"     );
    iEvent.put(sparm_subProcessId ,branchprefix+"subProcessId" );
    iEvent.put(sparm_comment      ,branchprefix+"comment"      );
}

// ------------ method called once each job just before starting event loop  ------------
void SParmMaker::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void SParmMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SParmMaker);
