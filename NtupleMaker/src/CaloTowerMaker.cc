// -*- C++ -*-
//
// Package:    CaloTowerMaker
// Class:      CaloTowerMaker
// 
/**\class CaloTowerMaker CaloTowerMaker.cc CMS2/NtupleMaker/src/CaloTowerMaker.cc

Description: <produce TaS collection of CaloTowers>

Implementation:
<Currently a bare copy of SCMaker>
 */
//
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "CMS2/NtupleMaker/interface/CaloTowerMaker.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/Ref.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;

//
// class decleration
//

//
// constructors and destructor
//
CaloTowerMaker::CaloTowerMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

	// number of towers in the event
	produces<unsigned int>("evtntwrs").setBranchAlias("evt_ntwrs");

	produces<std::vector<float> >(branchprefix+"eta").setBranchAlias(aliasprefix_+"_eta");
	produces<std::vector<float> >(branchprefix+"phi").setBranchAlias(aliasprefix_+"_phi");
	produces<std::vector<uint32_t> >(branchprefix+"detid").setBranchAlias(aliasprefix_+"_detid");

	// energy contributions from different detectors
	// energy in HO ("outerEnergy")is not included in "hadEnergy"
	//   double emEnergy() const { return emE_ ; }
	produces<std::vector<float> >(branchprefix+"emEnergy").setBranchAlias(aliasprefix_+"_emEnergy");
	//   double hadEnergy() const { return hadE_ ; }
	produces<std::vector<float> >(branchprefix+"hadEnergy").setBranchAlias(aliasprefix_+"_hadEnergy");
	//   double outerEnergy() const { return (id_.ietaAbs()<16)? outerE_ : 0.0; }
	produces<std::vector<float> >(branchprefix+"outerEnergy").setBranchAlias(aliasprefix_+"_outerEnergy");

	// transverse energies wrt to vtx (0,0,0)
	//   double emEt() const { return emE_ * sin( theta() ); }
	produces<std::vector<float> >(branchprefix+"emEt").setBranchAlias(aliasprefix_+"_emEt");
	//   double hadEt() const { return hadE_ * sin( theta() ); }
	produces<std::vector<float> >(branchprefix+"hadEt").setBranchAlias(aliasprefix_+"_hadEt");
	//   double outerEt() const { return (id_.ietaAbs()<16)? outerE_ * sin( theta() ) : 0.0; }
	produces<std::vector<float> >(branchprefix+"outerEt").setBranchAlias(aliasprefix_+"_outerEt");

	// recalculated wrt vertex provided as 3D point
	//   math::PtEtaPhiMLorentzVector p4(Point v) const;
	//   double p (Point v) const { return p4(v).P(); }
	produces<std::vector<float> >(branchprefix+"pcorr").setBranchAlias(aliasprefix_+"_pcorr");
	//   double et(Point v) const { return p4(v).Et(); }
	produces<std::vector<float> >(branchprefix+"etcorr").setBranchAlias(aliasprefix_+"_etcorr");
	//   double emEt(Point v)  const { return  emE_ * sin(p4(v).theta()); }
	produces<std::vector<float> >(branchprefix+"emEtcorr").setBranchAlias(aliasprefix_+"_emEtcorr");
	//   double hadEt(Point v) const { return  hadE_ * sin(p4(v).theta()); }
	produces<std::vector<float> >(branchprefix+"hadEtcorr").setBranchAlias(aliasprefix_+"_hadEtcorr");
	//   double outerEt(Point v) const { return (id_.ietaAbs()<16)? outerE_ * sin(p4(v).theta()) : 0.0; }
	produces<std::vector<float> >(branchprefix+"outerEtcorr").setBranchAlias(aliasprefix_+"_outerEtcorr");
	produces<std::vector<float> >(branchprefix+"etacorr").setBranchAlias(aliasprefix_+"_etacorr");
	produces<std::vector<float> >(branchprefix+"phicorr").setBranchAlias(aliasprefix_+"_phicorr");


	//    // the reference poins in ECAL and HCAL for direction determination
	//    // algorithm and parameters for selecting these points are set in the CaloTowersCreator
	//    const GlobalPoint& emPosition()  const { return emPosition_ ; }
	//    const GlobalPoint& hadPosition() const { return hadPosition_ ; }

	// time (ns) in ECAL/HCAL components of the tower based on weigted sum of the times in the contributing RecHits
	//   float ecalTime() const { return float(ecalTime_) * 0.01; }
	produces<std::vector<float> >(branchprefix+"ecalTime").setBranchAlias(aliasprefix_+"_ecalTime");
	//   float hcalTime() const { return float(hcalTime_) * 0.01; }
	produces<std::vector<float> >(branchprefix+"hcalTime").setBranchAlias(aliasprefix_+"_hcalTime");

	// methods to retrieve status information from the CaloTower:
	// number of bad/recovered/problematic cells in the tower
	// separately for ECAL and HCAL
	//  uint numBadEcalCells() const { return (twrStatusWord_ & 0x1F); }
	produces<std::vector<unsigned int> >(branchprefix+"numBadEcalCells").setBranchAlias(aliasprefix_+"_numBadEcalCells");
	//  uint numRecoveredEcalCells() const { return ((twrStatusWord_ >> 5) & 0x1F); }
	produces<std::vector<unsigned int> >(branchprefix+"numRecoveredEcalCells").setBranchAlias(aliasprefix_+"_numRecoveredEcalCells");
	//  uint numProblematicEcalCells() const { return ((twrStatusWord_ >> 10) & 0x1F); }
	produces<std::vector<unsigned int> >(branchprefix+"numProblematicEcalCells").setBranchAlias(aliasprefix_+"_numProblematicEcalCells");

	//  uint numBadHcalCells() const { return ( (twrStatusWord_ >> 15)& 0x7); }
	produces<std::vector<unsigned int> >(branchprefix+"numBadHcalCells").setBranchAlias(aliasprefix_+"_numBadHcalCells");
	//  uint numRecoveredHcalCells() const { return ((twrStatusWord_ >> 18) & 0x7); }
	produces<std::vector<unsigned int> >(branchprefix+"numRecoveredHcalCells").setBranchAlias(aliasprefix_+"_numRecoveredHcalCells");
	//  uint numProblematicHcalCells() const { return ((twrStatusWord_ >> 21) & 0x7); }
	produces<std::vector<unsigned int> >(branchprefix+"numProblematicHcalCells").setBranchAlias(aliasprefix_+"_numProblematicHcalCells");

	// chi2 prob
        produces<std::vector<float> >(branchprefix+"emMaxChi2Prob").setBranchAlias(aliasprefix_+"_emMaxChi2Prob");
	// vector of 10 samples per max crystal in the calo tower
	produces<std::vector<std::vector<int> > >(branchprefix+"emMaxEcalMGPASampleADC").setBranchAlias(aliasprefix_+"_emMaxEcalMGPASampleADC");
	// time of crystal with highest em energy
	produces<std::vector<float> >(branchprefix+"emMaxTime").setBranchAlias(aliasprefix_+"_emMaxTime");
	// the recoflag of the max energy crystal in the tower
	produces<std::vector<int> >(branchprefix+"emMaxRecoFlag").setBranchAlias(aliasprefix_+"_emMaxRecoFlag");
	// the energy of the max energy crystal in the tower
	produces<std::vector<float> >(branchprefix+"emMax").setBranchAlias(aliasprefix_+"_emMax");
	// the energy in 3x3 crystals centred on the max energy crystal
	produces<std::vector<float> >(branchprefix+"em3x3").setBranchAlias(aliasprefix_+"_em3x3");
	// as above for 5x5 crystals
	produces<std::vector<float> >(branchprefix+"em5x5").setBranchAlias(aliasprefix_+"_em5x5");
	// swiss cross
	produces<std::vector<float> >(branchprefix+"emSwiss").setBranchAlias(aliasprefix_+"_emSwiss");
	//number of crystals 
	produces<std::vector<int> >(branchprefix+"numCrystals").setBranchAlias(aliasprefix_+"_numCrystals");

	// add superclusters to the ntuple if they have ET > scEtMin_
	//   scEtMin_ = iConfig.getParameter<double>("scEtMin");

	// input Tags
	primaryVertexInputTag_ = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
	caloTowersInputTag_ = iConfig.getParameter<edm::InputTag>("caloTowersInputTag");
	ecalRecHitsInputTag_EE_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EE");
	ecalRecHitsInputTag_EB_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EB");
	ecalDigiProducerEE_     = iConfig.getParameter<edm::InputTag>("ecalDigiProducerEE");
	ecalDigiProducerEB_     = iConfig.getParameter<edm::InputTag>("ecalDigiProducerEB");

	// initialise this
	//  cachedCaloGeometryID_ = 0;

}

void CaloTowerMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	// calo topology
	edm::ESHandle<CaloTopology> pTopology;
	iSetup.get<CaloTopologyRecord>().get(pTopology);
	topology_ = pTopology.product();

	//   // get the calo geometry
	//   if (cachedCaloGeometryID_ != iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
	//     cachedCaloGeometryID_ = iSetup.get<CaloGeometryRecord>().cacheIdentifier();
	//     iSetup.get<CaloGeometryRecord>().get(caloGeometry_);
	//   }

	// get the primary vertices
	edm::Handle<reco::VertexCollection> vertexHandle;
	try {
		iEvent.getByLabel(primaryVertexInputTag_, vertexHandle);
	}
	catch ( cms::Exception& ex ) {
		edm::LogError("CaloTowerMakerError") << "Error! can't get the primary vertex";
	}
	const reco::VertexCollection *vertexCollection = vertexHandle.product();
	Point pv(0.0, 0.0, 0.0);
	if (vertexCollection->size() > 0) {
		pv = vertexCollection->at(0).position();
	}

	edm::Handle<CaloTowerCollection> calotower;
	iEvent.getByLabel(caloTowersInputTag_,calotower);

	if(!calotower.isValid()) {
		edm::LogError("CaloTowerMakerError") << "Error! Can't get calotowers!" << std::endl;
		//    return ;
	}

	// get hits
	edm::Handle<EcalRecHitCollection> rhcHandleEE;
	iEvent.getByLabel(ecalRecHitsInputTag_EE_, rhcHandleEE);
	const EcalRecHitCollection *recHitsEE = rhcHandleEE.product();

	edm::Handle<EcalRecHitCollection> rhcHandleEB;
	iEvent.getByLabel(ecalRecHitsInputTag_EB_, rhcHandleEB);
	const EcalRecHitCollection *recHitsEB = rhcHandleEB.product();

	// ECAL Digis
	edm::Handle<EBDigiCollection> ebDigiHandle;
	iEvent.getByLabel(ecalDigiProducerEB_, ebDigiHandle);
	edm::Handle<EEDigiCollection> eeDigiHandle;
	iEvent.getByLabel(ecalDigiProducerEE_, eeDigiHandle);

	const EBDigiCollection *ebDigis = 0;
	const EEDigiCollection *eeDigis = 0;
	digi_ = false;
	if (ebDigiHandle.isValid() && eeDigiHandle.isValid()) {
		digi_ = true;
		ebDigis = ebDigiHandle.product();
		eeDigis = eeDigiHandle.product();
	}

	// ecal cluster shape variables
	// do not use the lazy tools because need to get the hits anyway
	EcalClusterTools clusterTools;

	//   // get hoe variable
	//   HoECalculator hoeCalc(caloGeometry_);

	std::auto_ptr<unsigned int> evt_ntwrs (new unsigned int);

	std::auto_ptr<std::vector<float> > vector_twrs_eta (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_phi (new std::vector<float>);
	std::auto_ptr<std::vector<uint32_t> > vector_twrs_detid (new std::vector<uint32_t>);

	std::auto_ptr<std::vector<float> > vector_twrs_emEnergy (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_hadEnergy (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_outerEnergy (new std::vector<float>);

	std::auto_ptr<std::vector<float> > vector_twrs_emEt (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_hadEt (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_outerEt (new std::vector<float>);

	std::auto_ptr<std::vector<float> > vector_twrs_pcorr (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_etcorr (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_emEtcorr (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_hadEtcorr (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_outerEtcorr (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_etacorr (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_phicorr (new std::vector<float>);

	std::auto_ptr<std::vector<float> > vector_twrs_ecalTime (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_hcalTime (new std::vector<float>);

	std::auto_ptr<std::vector<unsigned int> > vector_twrs_numBadEcalCells (new std::vector<unsigned int>);
	std::auto_ptr<std::vector<unsigned int> > vector_twrs_numRecoveredEcalCells (new std::vector<unsigned int>);
	std::auto_ptr<std::vector<unsigned int> > vector_twrs_numProblematicEcalCells (new std::vector<unsigned int>);

	std::auto_ptr<std::vector<unsigned int> > vector_twrs_numBadHcalCells (new std::vector<unsigned int>);
	std::auto_ptr<std::vector<unsigned int> > vector_twrs_numRecoveredHcalCells (new std::vector<unsigned int>);
	std::auto_ptr<std::vector<unsigned int> > vector_twrs_numProblematicHcalCells (new std::vector<unsigned int>);

	std::auto_ptr<std::vector<float> > vector_twrs_emMaxChi2Prob (new std::vector<float>);
	std::auto_ptr<std::vector<std::vector<int> > > vector_twrs_emMaxEcalMGPASampleADC (new std::vector<std::vector<int> >);
	std::auto_ptr<std::vector<float> > vector_twrs_emMaxTime      (new std::vector<float>);
	std::auto_ptr<std::vector<int> > vector_twrs_emMaxRecoFlag      (new std::vector<int>);
	std::auto_ptr<std::vector<float> > vector_twrs_emMax      (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_em3x3      (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_em5x5      (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_twrs_emSwiss	(new std::vector<float>);	
	std::auto_ptr<std::vector<int>   > vector_twrs_numCrystals(new std::vector<int>);

	*evt_ntwrs = 0;

	for (CaloTowerCollection::const_iterator j = calotower->begin();j !=calotower->end(); j++) {

		*evt_ntwrs +=1;

		//     std::cout << *j << std::endl;
		//     std::cout << "ENERGY HAD " << j->hadEnergy()<< " ENERGY EM " <<j->emEnergy() 
		//          << " ETA " <<j->eta() << " PHI " <<j->phi() << std::endl;

		vector_twrs_eta->push_back(j->eta());
		vector_twrs_phi->push_back(j->phi());
		vector_twrs_detid->push_back(j->id().rawId());

		vector_twrs_emEnergy->push_back(j->emEnergy());
		vector_twrs_hadEnergy->push_back(j->hadEnergy());
		vector_twrs_outerEnergy->push_back(j->outerEnergy());

		vector_twrs_emEt->push_back(j->emEt());
		vector_twrs_hadEt->push_back(j->hadEt());
		vector_twrs_outerEt->push_back(j->outerEt());

		vector_twrs_pcorr->push_back(j->p(pv));
		vector_twrs_etcorr->push_back(j->et(pv));
		vector_twrs_emEtcorr->push_back(j->emEt(pv));
		vector_twrs_hadEtcorr->push_back(j->hadEt(pv));
		vector_twrs_outerEtcorr->push_back(j->outerEt(pv));
		vector_twrs_etacorr->push_back(j->p4(pv).eta());
		vector_twrs_phicorr->push_back(j->p4(pv).phi());

		vector_twrs_ecalTime->push_back(j->ecalTime());
		vector_twrs_hcalTime->push_back(j->hcalTime());

		vector_twrs_numBadEcalCells->push_back(j->numBadEcalCells());
		vector_twrs_numRecoveredEcalCells->push_back(j->numRecoveredEcalCells());
		vector_twrs_numProblematicEcalCells->push_back(j->numProblematicEcalCells());

		vector_twrs_numBadHcalCells->push_back(j->numBadHcalCells());
		vector_twrs_numRecoveredHcalCells->push_back(j->numRecoveredHcalCells());
		vector_twrs_numProblematicHcalCells->push_back(j->numProblematicHcalCells());

		// find the detids of the crystals in this calo tower
		const std::vector<DetId> &towerDetIds = j->constituents(); 

		// get variables for highest em energy crystal in tower
		float emE = 0.0;
		float emMax = 0.0;
		float em3x3 = 0.0;
		float em5x5 = 0.0;
		float emSwiss = 0.0;
		float chi2Prob = -9999.99;
		float emMaxTime = -9999.99;
		DetId emMaxId(0);
		int recoFlag = -1;
		std::vector<int> ecalMGPASampleADC;

		// loop on detids in the tower
		for (size_t i = 0; i < towerDetIds.size(); ++i) {
			// find the energy of this detId if it is in the ecal
			if (towerDetIds[i].subdetId() == EcalEndcap)
				emE = clusterTools.recHitEnergy(towerDetIds[i], recHitsEE);
			if (towerDetIds[i].subdetId() == EcalBarrel)
				emE = clusterTools.recHitEnergy(towerDetIds[i], recHitsEB);
			// compare with previous highest energy
			if (emE > emMax) {
				emMax = emE; 
				emMaxId = towerDetIds[i];
			}
		}

		// find the relevant quantities for the identified crystal
		if (emMaxId != DetId(0)) {
			reco::BasicCluster dummyCluster;
			if (emMaxId.subdetId() == EcalEndcap) {
				chi2Prob = recHitChi2Prob(emMaxId, recHitsEE);
				emMaxTime = recHitTime(emMaxId, recHitsEE);
				recoFlag = recHitFlag(emMaxId, recHitsEE);
				recHitSamples(emMaxId, eeDigis, ecalMGPASampleADC);
				em3x3 = clusterTools.matrixEnergy(dummyCluster, recHitsEE, topology_, emMaxId, -1, 1, -1, 1);
				em5x5 = clusterTools.matrixEnergy(dummyCluster, recHitsEE, topology_, emMaxId, -2, 2, -2, 2);

				// make the swiss cross
				emSwiss = clusterTools.matrixEnergy(dummyCluster, recHitsEE, topology_, emMaxId, 0, 0, -1, 1);
				emSwiss += clusterTools.matrixEnergy(dummyCluster, recHitsEE, topology_, emMaxId, -1, 1, 0, 0);
				emSwiss -= emMax;
			}

			if (emMaxId.subdetId() == EcalBarrel) { 
                                chi2Prob = recHitChi2Prob(emMaxId, recHitsEB);
				emMaxTime = recHitTime(emMaxId, recHitsEB);
				recoFlag = recHitFlag(emMaxId, recHitsEB);
				recHitSamples(emMaxId, ebDigis, ecalMGPASampleADC);
				em3x3 = clusterTools.matrixEnergy(dummyCluster, recHitsEB, topology_, emMaxId, -1, 1, -1, 1);
				em5x5 = clusterTools.matrixEnergy(dummyCluster, recHitsEB, topology_, emMaxId, -2, 2, -2, 2);

                                // make the swiss cross
                                emSwiss = clusterTools.matrixEnergy(dummyCluster, recHitsEB, topology_, emMaxId, 0, 0, -1, 1);
                                emSwiss += clusterTools.matrixEnergy(dummyCluster, recHitsEB, topology_, emMaxId, -1, 1, 0, 0);
                                emSwiss -= emMax;

			}
		}

		vector_twrs_emMaxChi2Prob->push_back(chi2Prob);
		vector_twrs_emMaxEcalMGPASampleADC->push_back(ecalMGPASampleADC);
		vector_twrs_emMaxTime->push_back(emMaxTime);
		vector_twrs_emMaxRecoFlag->push_back(recoFlag);
		vector_twrs_emMax->push_back(emMax);
		vector_twrs_em3x3->push_back(em3x3);
		vector_twrs_em5x5->push_back(em5x5);
		vector_twrs_emSwiss->push_back(emSwiss);
		vector_twrs_numCrystals->push_back(j->numCrystals());

	}

	// put results into the event
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

	iEvent.put(evt_ntwrs, "evtntwrs");
	iEvent.put(vector_twrs_eta, branchprefix+"eta");
	iEvent.put(vector_twrs_phi, branchprefix+"phi");
	iEvent.put(vector_twrs_detid, branchprefix+"detid");
	iEvent.put(vector_twrs_emEnergy, branchprefix+"emEnergy");
	iEvent.put(vector_twrs_hadEnergy, branchprefix+"hadEnergy");
	iEvent.put(vector_twrs_outerEnergy, branchprefix+"outerEnergy");

	iEvent.put(vector_twrs_emEt, branchprefix+"emEt");
	iEvent.put(vector_twrs_hadEt, branchprefix+"hadEt");
	iEvent.put(vector_twrs_outerEt, branchprefix+"outerEt");

	iEvent.put(vector_twrs_pcorr, branchprefix+"pcorr");
	iEvent.put(vector_twrs_etcorr, branchprefix+"etcorr");
	iEvent.put(vector_twrs_emEtcorr, branchprefix+"emEtcorr");
	iEvent.put(vector_twrs_hadEtcorr, branchprefix+"hadEtcorr");
	iEvent.put(vector_twrs_outerEtcorr, branchprefix+"outerEtcorr");
	iEvent.put(vector_twrs_etacorr, branchprefix+"etacorr");
	iEvent.put(vector_twrs_phicorr, branchprefix+"phicorr");

	iEvent.put(vector_twrs_ecalTime, branchprefix+"ecalTime");
	iEvent.put(vector_twrs_hcalTime, branchprefix+"hcalTime");

	iEvent.put(vector_twrs_numBadEcalCells, branchprefix+"numBadEcalCells");
	iEvent.put(vector_twrs_numRecoveredEcalCells, branchprefix+"numRecoveredEcalCells");
	iEvent.put(vector_twrs_numProblematicEcalCells, branchprefix+"numProblematicEcalCells");

	iEvent.put(vector_twrs_numBadHcalCells, branchprefix+"numBadHcalCells");
	iEvent.put(vector_twrs_numRecoveredHcalCells, branchprefix+"numRecoveredHcalCells");
	iEvent.put(vector_twrs_numProblematicHcalCells, branchprefix+"numProblematicHcalCells");

	iEvent.put(vector_twrs_emMaxChi2Prob, branchprefix+"emMaxChi2Prob");
	iEvent.put(vector_twrs_emMaxEcalMGPASampleADC, branchprefix+"emMaxEcalMGPASampleADC");
	iEvent.put(vector_twrs_emMaxTime, branchprefix+"emMaxTime");
	iEvent.put(vector_twrs_emMaxRecoFlag, branchprefix+"emMaxRecoFlag");
	iEvent.put(vector_twrs_emMax, branchprefix+"emMax");
	iEvent.put(vector_twrs_em3x3, branchprefix+"em3x3");
	iEvent.put(vector_twrs_em5x5, branchprefix+"em5x5");
	iEvent.put(vector_twrs_emSwiss, branchprefix+"emSwiss");
	iEvent.put(vector_twrs_numCrystals, branchprefix+"numCrystals");

}

float CaloTowerMaker::recHitChi2Prob(DetId emMaxId, const EcalRecHitCollection *recHits)
{               
        EcalRecHitCollection::const_iterator it = recHits->find(emMaxId);
        if (it != recHits->end()) return it->chi2Prob();
        return -9999.99;
}

float CaloTowerMaker::recHitTime(DetId emMaxId, const EcalRecHitCollection *recHits)
{
	EcalRecHitCollection::const_iterator it = recHits->find(emMaxId);
	if (it != recHits->end()) return it->time();
	return -9999.99;
}

int CaloTowerMaker::recHitFlag(DetId emMaxId, const EcalRecHitCollection *recHits)
{
	EcalRecHitCollection::const_iterator it = recHits->find(emMaxId);
	if (it != recHits->end()) return it->recoFlag();
	return -1;
}

void CaloTowerMaker::recHitSamples(DetId emMaxId, const EcalDigiCollection *digis, std::vector<int> &samples)
{


	samples.clear();
	if (digi_) {
		EcalDigiCollection::const_iterator it = digis->find(emMaxId);
		EcalDataFrame frame;
		if (it != digis->end()) {
			frame = (*it);
			for (size_t i = 0; i < 10; ++i) samples.push_back(frame.sample(i).adc());
		}
	}
	else {
		for (size_t i = 0; i < 10; ++i) samples.push_back(0);
	}

}

// ------------ method called once each job just before starting event loop  ------------
	void 
CaloTowerMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
CaloTowerMaker::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(CaloTowerMaker);

