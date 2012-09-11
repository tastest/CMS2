//now make the source file
#include "LeptonTreeMaker.h"

#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTreeCache.h"
#include "TPRegexp.h"

#include "../HWW2012CORE/analysisObjects.h"
#include "../HWW2012CORE/analysisTools.h"
#include "../HWW2012CORE/analysisSelections.h"

using namespace std;
#include "math.h"

#include "../Tools/MuonIDMVA.h"
#include "../Tools/goodrun.h"
#include "../Tools/EGammaMvaEleEstimator.h"
#include "../Tools/MuonMVAEstimator.h"
#include "../Tools/MuonEffectiveArea.h"

//#ifndef __CINT__
#include "../CORE/CMS2.h"
#include "../CORE/electronSelections.h"
#include "../CORE/MITConversionUtilities.h"
#include "../CORE/muonSelections.h"
//#endif

MuonIDMVA* muonIdMVA_leptree = 0;
EGammaMvaEleEstimator* egammaMvaEleEstimator_leptree = 0;
MuonMVAEstimator* muonMVAEstimator_leptree = 0;

double electron_d0PV_wwV1_local(unsigned int index) { 
    
	if ( cms2.vtxs_sumpt().empty() ) return 9999.;
    
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


//
// constructor and destructor
//

LeptonTreeMaker::LeptonTreeMaker(bool lockToCoreSelectors, bool useLHeleId, bool useMVAeleId, bool doDYNNLOw, unsigned int prescale, bool realData) 
: lockToCoreSelectors_(lockToCoreSelectors), useLHeleId_(useLHeleId), useMVAeleId_(useMVAeleId), 
    doDYNNLOw_(doDYNNLOw), prescale_(prescale), realData_(realData)
{
    std::cout << "LeptonTreeMaker::LeptonTreeMaker" << std::endl;
    std::cout << "\t lockToCoreSelectors: " << lockToCoreSelectors << std::endl;
    std::cout << "\t useLHeleId         : " << useLHeleId << std::endl;
    std::cout << "\t useMVAeleId        : " << useMVAeleId << std::endl;
    std::cout << "\t doDYNNLOw          : " << doDYNNLOw << std::endl;
    std::cout << "\t prescale           : " << prescale << std::endl;
    std::cout << "\t realData           : " << realData << std::endl;


	// --------------- EGamma Id MVA  --------------------------
	vector<std::string> egammaweights;
	egammaweights.push_back("../HWW2012CORE/files/Electrons_BDTG_TrigV0_Cat1.weights.xml"); 
	egammaweights.push_back("../HWW2012CORE/files/Electrons_BDTG_TrigV0_Cat2.weights.xml"); 
	egammaweights.push_back("../HWW2012CORE/files/Electrons_BDTG_TrigV0_Cat3.weights.xml"); 
	egammaweights.push_back("../HWW2012CORE/files/Electrons_BDTG_TrigV0_Cat4.weights.xml"); 
	egammaweights.push_back("../HWW2012CORE/files/Electrons_BDTG_TrigV0_Cat5.weights.xml"); 
	egammaweights.push_back("../HWW2012CORE/files/Electrons_BDTG_TrigV0_Cat6.weights.xml"); 
	egammaMvaEleEstimator_leptree = new EGammaMvaEleEstimator();
	egammaMvaEleEstimator_leptree->initialize("BDT", EGammaMvaEleEstimator::kTrig, true, egammaweights );

	// --------------- Muon RingIso MVA  --------------------------
	vector<std::string> muonisoweights;
	muonisoweights.push_back("../HWW2012CORE/files/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
	muonisoweights.push_back("../HWW2012CORE/files/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
	muonisoweights.push_back("../HWW2012CORE/files/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
	muonisoweights.push_back("../HWW2012CORE/files/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
	muonisoweights.push_back("../HWW2012CORE/files/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml");
	muonisoweights.push_back("../HWW2012CORE/files/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml");
	muonMVAEstimator_leptree =new MuonMVAEstimator();
	muonMVAEstimator_leptree->initialize( "MuonIso_BDTG_IsoRings", MuonMVAEstimator::kIsoRings, true, muonisoweights );

    //
    // set up jet corrections
    //
    std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
    try {
        jetcorr_filenames_pfL1FastJetL2L3.clear();
        if (realData) {
            jetcorr_filenames_pfL1FastJetL2L3.push_back("../HWW2012CORE/files/GR_R_52_V7_L1FastJet_AK5PF.txt"); 
            jetcorr_filenames_pfL1FastJetL2L3.push_back("../HWW2012CORE/files/GR_R_52_V7_L2Relative_AK5PF.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("../HWW2012CORE/files/GR_R_52_V7_L3Absolute_AK5PF.txt"); 
            jetcorr_filenames_pfL1FastJetL2L3.push_back("../HWW2012CORE/files/GR_R_52_V7_L2L3Residual_AK5PF.txt");
        } else {
			jetcorr_filenames_pfL1FastJetL2L3.push_back("../HWW2012CORE/files/START52_V9_L1FastJet_AK5PF.txt");
			jetcorr_filenames_pfL1FastJetL2L3.push_back("../HWW2012CORE/files/START52_V9_L2Relative_AK5PF.txt");
			jetcorr_filenames_pfL1FastJetL2L3.push_back("../HWW2012CORE/files/START52_V9_L3Absolute_AK5PF.txt");
        }
        jet_corrector_pfL1FastJetL2L3_ = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);
    } catch (...){
        cout << "\nFailed to setup correctors needed to get Jet Enetry Scale. Abort\n" << endl;
        assert(0);
    }

    //
    // other parameters
    //

    batchMode_ = false;

}

LeptonTreeMaker::~LeptonTreeMaker() {
    std::cout << "LeptonTreeMaker::~LeptonTreeMaker" << std::endl;
    //if (electronIdMVA_ != 0) delete electronIdMVA_;
    if (jet_corrector_pfL1FastJetL2L3_ != 0) delete jet_corrector_pfL1FastJetL2L3_;
}

//
// process data 
//

void LeptonTreeMaker::ScanChain(TString outfileid,
        TChain* chain, 
		SmurfTree::DataType sample,
        double integratedLumi, // in unit of pb^-1, if negative the weight is 1.
        double xsec,           // in unit of pb, if negative take it from evt_xsec_excl*evt_kfactor
        int nProcessedEvents,  // if negative, take it from evt_nEvts
        bool identifyEvents, 
        bool realData,
        TString cms2_json_file)
{

    std::string prefix = SmurfTree::name(sample);
    std::cout << sample << ", " << prefix << std::endl;

    if ( chain->GetListOfFiles()->GetEntries()==0 ){
        printf("\nERROR: chain is empty for sample: %s\n\n",prefix.c_str());
        assert(0);
    }

    unsigned int nEventsChain = chain->GetEntries();  // number of entries in chain --> number of events from all files
    gErrorIgnoreLevel = -1;
    unsigned int nEventsTotal = 0;

    // make smurf ntuples
    gSystem->MakeDirectory("smurf");
    TFile* fSmurf = TFile::Open(Form("smurf/%s_%s.root",prefix.c_str(), outfileid.Data()),"RECREATE");
    assert(fSmurf);
    LeptonTree leptonTree;
    leptonTree.CreateTree();
    leptonTree.tree_->SetDirectory(fSmurf);

	// triggers whose efficiency to be measured
	// electrons
	HLT_Ele17_Ele8_tag_=0;		HLT_Ele17_Ele8_probe_=0;
	HLT_Ele27_WP80_tag_=0;		HLT_Ele27_WP80_probe_=0;
	HLT_Ele17_Ele8_tag_=0;				HLT_Ele17_Ele8_probe_=0;
	HLT_Ele17_Ele8_LeadingLeg_tag_=0; 	HLT_Ele17_Ele8_LeadingLeg_probe_=0;
	HLT_Ele17_Ele8_TrailingLeg_tag_=0;	HLT_Ele17_Ele8_TrailingLeg_probe_=0;

	// muons
	HLT_IsoMu24_eta2p1_tag_=0;			HLT_IsoMu24_eta2p1_probe_=0;
	HLT_Mu17_Mu8_tag_=0;				HLT_Mu17_Mu8_probe_=0;		
	HLT_Mu17_Mu8_LeadingLeg_tag_=0;		HLT_Mu17_Mu8_LeadingLeg_probe_=0;		
	HLT_Mu17_Mu8_TrailingLeg_tag_=0;	HLT_Mu17_Mu8_TrailingLeg_probe_=0;
	HLT_Mu17_TkMu8_tag_=0;				HLT_Mu17_TkMu8_probe_=0;	
	HLT_Mu17_TkMu8_LeadingLeg_tag_=0;	HLT_Mu17_TkMu8_LeadingLeg_probe_=0;	
	HLT_Mu17_TkMu8_TrailingLeg_tag_=0;	HLT_Mu17_TkMu8_TrailingLeg_probe_=0;

	leptonTree.tree_->Branch("HLT_Ele17_Ele8_tag"					, &HLT_Ele17_Ele8_tag_				,"HLT_Ele17_Ele8_tag/i");
	leptonTree.tree_->Branch("HLT_Ele17_Ele8_probe"					, &HLT_Ele17_Ele8_probe_			,"HLT_Ele17_Ele8_probe/i");
	leptonTree.tree_->Branch("HLT_Ele27_WP80_tag"					, &HLT_Ele27_WP80_tag_				,"HLT_Ele27_WP80_tag/i");
	leptonTree.tree_->Branch("HLT_Ele27_WP80_probe"					, &HLT_Ele27_WP80_probe_			,"HLT_Ele27_WP80_probe/i");
	leptonTree.tree_->Branch("HLT_Ele17_Ele8_tag"					, &HLT_Ele17_Ele8_tag_				,"HLT_Ele17_Ele8_tag/i");
	leptonTree.tree_->Branch("HLT_Ele17_Ele8_probe"					, &HLT_Ele17_Ele8_probe_			,"HLT_Ele17_Ele8_probe/i");
	leptonTree.tree_->Branch("HLT_Ele17_Ele8_LeadingLeg_tag"			, &HLT_Ele17_Ele8_LeadingLeg_tag_		,"HLT_Ele17_Ele8_LeadingLeg_tag/i");
	leptonTree.tree_->Branch("HLT_Ele17_Ele8_LeadingLeg_probe"			, &HLT_Ele17_Ele8_LeadingLeg_probe_		,"HLT_Ele17_Ele8_LeadingLeg_probe/i");
	leptonTree.tree_->Branch("HLT_Ele17_Ele8_TrailingLeg_tag"			, &HLT_Ele17_Ele8_TrailingLeg_tag_		,"HLT_Ele17_Ele8_TrailingLeg_tag/i");
	leptonTree.tree_->Branch("HLT_Ele17_Ele8_TrailingLeg_probe"			, &HLT_Ele17_Ele8_TrailingLeg_probe_	,"HLT_Ele17_Ele8_TrailingLeg_probe/i");

	leptonTree.tree_->Branch("HLT_IsoMu24_eta2p1_tag"				, &HLT_IsoMu24_eta2p1_tag_			,"HLT_IsoMu24_eta2p1_tag/i");
	leptonTree.tree_->Branch("HLT_IsoMu24_eta2p1_probe"				, &HLT_IsoMu24_eta2p1_probe_		,"HLT_IsoMu24_eta2p1_probe/i");
	leptonTree.tree_->Branch("HLT_Mu17_TkMu8_tag"					, &HLT_Mu17_TkMu8_tag_				,"HLT_Mu17_TkMu8_tag/i");
	leptonTree.tree_->Branch("HLT_Mu17_TkMu8_probe"					, &HLT_Mu17_TkMu8_probe_			,"HLT_Mu17_TkMu8_probe/i");
	leptonTree.tree_->Branch("HLT_Mu17_TkMu8_LeadingLeg_tag"		, &HLT_Mu17_TkMu8_LeadingLeg_tag_		,"HLT_Mu17_TkMu8_LeadingLeg_tag/i");
	leptonTree.tree_->Branch("HLT_Mu17_TkMu8_LeadingLeg_probe"		, &HLT_Mu17_TkMu8_LeadingLeg_probe_		,"HLT_Mu17_TkMu8_LeadingLeg_probe/i");
	leptonTree.tree_->Branch("HLT_Mu17_TkMu8_TrailingLeg_tag"		, &HLT_Mu17_TkMu8_TrailingLeg_tag_		,"HLT_Mu17_TkMu8_TrailingLeg_tag/i");
	leptonTree.tree_->Branch("HLT_Mu17_TkMu8_TrailingLeg_probe"		, &HLT_Mu17_TkMu8_TrailingLeg_probe_	,"HLT_Mu17_TkMu8_TrailingLeg_probe/i");
	leptonTree.tree_->Branch("HLT_Mu17_Mu8_tag"						, &HLT_Mu17_Mu8_tag_				,"HLT_Mu17_Mu8_tag/i");
	leptonTree.tree_->Branch("HLT_Mu17_Mu8_probe"					, &HLT_Mu17_Mu8_probe_				,"HLT_Mu17_Mu8_probe/i");
	leptonTree.tree_->Branch("HLT_Mu17_Mu8_LeadingLeg_tag"			, &HLT_Mu17_Mu8_LeadingLeg_tag_		,"HLT_Mu17_Mu8_LeadingLeg_tag/i");
	leptonTree.tree_->Branch("HLT_Mu17_Mu8_LeadingLeg_probe"		, &HLT_Mu17_Mu8_LeadingLeg_probe_	,"HLT_Mu17_Mu8_LeadingLeg_probe/i");
	leptonTree.tree_->Branch("HLT_Mu17_Mu8_TrailingLeg_tag"			, &HLT_Mu17_Mu8_TrailingLeg_tag_	,"HLT_Mu17_Mu8_TrailingLeg_tag/i");
	leptonTree.tree_->Branch("HLT_Mu17_Mu8_TrailingLeg_probe"		, &HLT_Mu17_Mu8_TrailingLeg_probe_	,"HLT_Mu17_Mu8_TrailingLeg_probe/i");

	// FR triggers
	HLT_Ele8_CaloIdL_TrkIdVL_=0;
	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_=0;
	HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_=0;
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_=0;
	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_=0;
	HLT_Mu8_=0;
	HLT_Mu17_=0;

	leptonTree.tree_->Branch(	"HLT_Ele8_CaloIdL_TrkIdVL", 
								&HLT_Ele8_CaloIdL_TrkIdVL_,	
								"HLT_Ele8_CaloIdL_TrkIdVL/i");
	leptonTree.tree_->Branch(	"HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 
								&HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_,
								"HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL/i");
	leptonTree.tree_->Branch(	"HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30",
	 							&HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_,
								"HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30/i");
	leptonTree.tree_->Branch(	"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 
								&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_,
								"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL/i");
	leptonTree.tree_->Branch(	"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30", 
								&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_,
								"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30/i");
	leptonTree.tree_->Branch("HLT_Mu8"					, &HLT_Mu8_				,	"HLT_Mu8/i");
	leptonTree.tree_->Branch("HLT_Mu17"					, &HLT_Mu17_			,	"HLT_Mu17/i");

	// Initialization :  TEST variables
	mva_=-999.;
	iso_=-999.;
	fbrem_=-999.; 
	kfchi2_=-999.;
	kfhits_=0;
	gsfchi2_=-999.;
	deta_=-999.;
	dphi_=-999.;
	detacalo_=-999.;
	see_=-999.;
	spp_=-999.;
	etawidth_=-999.;
	phiwidth_=-999.;
	e1x5e5x5_=-999.;
	R9_=-999.;
	HoE_=-999.;
	EoP_=-999.;
	IoEmIoP_=-999.;
	eleEoPout_=-999.;
	PreShowerOverRaw_=-999.;
	d0_=-999.;
	ip3d_=-999.;
	eta_=-999.;
	pt_=-999;

	leptonTree.tree_->Branch("mva"					, &mva_				,	"mva/F");
	leptonTree.tree_->Branch("iso"					, &iso_				,	"iso/F");
	leptonTree.tree_->Branch("fbrem"				, &fbrem_			,	"fbrem/F");
	leptonTree.tree_->Branch("kfchi2_"				, &kfchi2_			,	"kfchi2/F");
	leptonTree.tree_->Branch("kfhits"				, &kfhits_			,	"kfhits/I");
	leptonTree.tree_->Branch("gsfchi2"				, &gsfchi2_			,	"gsfchi2/F");
	leptonTree.tree_->Branch("deta"					, &deta_			,	"deta/F");
	leptonTree.tree_->Branch("dphi"					, &dphi_			,	"dphi/F");
	leptonTree.tree_->Branch("detacalo"				, &detacalo_		,	"detacalo/F");
	leptonTree.tree_->Branch("see"					, &see_				,	"see/F");
	leptonTree.tree_->Branch("spp"					, &spp_				,	"spp/F");
	leptonTree.tree_->Branch("etawidth"				, &etawidth_		,	"etawidth/F");
	leptonTree.tree_->Branch("phiwidth"				, &phiwidth_		,	"phiwidth/F");
	leptonTree.tree_->Branch("e1x5e5x5"				, &e1x5e5x5_		,	"e1x5e5x5/F");
	leptonTree.tree_->Branch("R9"					, &R9_				,	"R9/F");
	leptonTree.tree_->Branch("HoE"					, &HoE_				,	"HoE/F");
	leptonTree.tree_->Branch("EoP"					, &EoP_				,	"EoP/F");
	leptonTree.tree_->Branch("IoEmIoP"				, &IoEmIoP_			,	"IoEmIoP/F");
	leptonTree.tree_->Branch("eleEoPout"			, &eleEoPout_		,	"eleEoPout/F");
	leptonTree.tree_->Branch("PreShowerOverRaw"		, &PreShowerOverRaw_,	"PreShowerOverRaw/F");
	leptonTree.tree_->Branch("d0"					, &d0_				,	"d0/F");
	leptonTree.tree_->Branch("ip3d"					, &ip3d_			,	"ip3d/F");
	leptonTree.tree_->Branch("eta"					, &eta_				,	"eta/F");
	leptonTree.tree_->Branch("pt"					, &pt_				,	"pt/F");

    std::map<unsigned int, std::set<unsigned int> > runList;

    // clear list of duplicates
    already_seen.clear();
    int duplicates_total_n = 0;
    double duplicates_total_weight = 0;
  	int nFailedIdentification = 0;
    int nFilteredOut = 0;

    int i_permille_old = 0;

    if(realData) {
        if(cms2_json_file=="") {  
            std::cout<<"\n WARNING: Running on the real Data, but not JSON file!"<<std::endl;
            // assert(0);
        } else {
            //set_goodrun_file_json(cms2_json_file);
            set_goodrun_file(cms2_json_file);

        }
    }

    vector<TH2D*>     fDYNNLOKFactorHists;           //vector of hist for Drell-Yan NNLO Kfactor
    if (doDYNNLOw_ && (sample==SmurfTree::dyee || sample==SmurfTree::dymm || sample==SmurfTree::dytt) ) {
        TFile *tmpFile = new TFile( "../HWW2012CORE/files/fewz_powheg_weights_stepwise_2011_fine7.root", "READ");
        const int nMassBins = 41;
        for(int i=0; i<nMassBins; i++){         
            TString hname = TString::Format("weight_%02d",i+1);
            TH2D *tmpHist = (TH2D*)tmpFile->Get(hname);
            tmpHist->SetDirectory(0);
            fDYNNLOKFactorHists.push_back(tmpHist);
        }
        tmpFile->Close();
        delete tmpFile;
    }

    // file loop
    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);

    while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
        //        printf("current file: %s (%s), %s\n", currentFile->GetName(), 
        // 	      currentFile->GetTitle(), currentFile->IsA()->GetName());

        TFile *f = TFile::Open(currentFile->GetTitle()); 
        assert(f);
        TTree *tree = (TTree*)f->Get("Events");
        assert(tree);
        TTreeCache::SetLearnEntries(10);
        tree->SetCacheSize(128*1024*1024);
        cms2.Init(tree);  // set branch addresses for TTree tree
        TStopwatch t;

        //Event Loop
        unsigned int nEvents = tree->GetEntries();

        for( unsigned int event = 0; event < nEvents; ++event) {
            tree->LoadTree(event);
            cms2.GetEntry(event);  // get entries for Event number event from branches of TTree tree
            if (cms2.evt_event() % prescale_ !=0) continue;

            // Select the good runs from the json file
            if(realData && cms2_json_file!="") {
                if( !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()) ) continue;
            }
            runList[cms2.evt_run()].insert(cms2.evt_lumiBlock());

            double weight = 1.0;
            if ( !realData && integratedLumi>0 ){
                double mcweight = cms2.genps_weight() > 0.0 ? 1.0 : -1.0;
                weight = integratedLumi * mcweight * (xsec>0?xsec:cms2.evt_xsec_excl()*cms2.evt_kfactor()*cms2.evt_filt_eff()) /
                    (nProcessedEvents>0?nProcessedEvents:cms2.evt_nEvts());
            }       
            ++nEventsTotal;

            //
            // basic event preselection
            //

            if (cms2.trks_d0().size() == 0) continue;  // needed to get rid of back Monte Carlo events in CMSSW_2_X analysis
            if (cms2.hyp_type().size() == 0) continue; // skip events without hypothesis

            //
            // duplicate removal
            //

            EventIdentifier id(cms2,realData);
            if (is_duplicate(id)) {
                duplicates_total_n++;
                if(!realData) duplicates_total_weight += cms2.evt_scale1fb();
                //cout << "Duplicate event found. Run: " << cms2.evt_run() << ", Event:" << cms2.evt_event() << ", Lumi: " << cms2.evt_lumiBlock() << endl;
                continue;
            }

            //
            // counter
            //

//            if (!batchMode_) {
                int i_permille = (int)floor(10 * nEventsTotal / float(nEventsChain));
                if (i_permille != i_permille_old) {
                    // xterm magic from L. Vacavant and A. Cerri
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                            "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
                    fflush(stdout);
                    i_permille_old = i_permille;
                }
//            }

            //
            // filter out invalid contributions
            // to the sample being analyzed
            //

            if ( identifyEvents ){
                // check if we know what we are looking at
                if ( ! isIdentified(sample) ) nFailedIdentification++;

                // filter by process
                if ( ! filterByProcess(sample) ) {
                    nFilteredOut++;
                    continue;
                }
            }

            //
            // fill lepton tree
            //

			unsigned int eventSelection = 0;

			//
			// electron FR
			//
			TPMERegexp regexp_ele8					(	"HLT_Ele8_CaloIdL_TrkIdVL_v(\\d+)", "o"								);
			TPMERegexp regexp_ele8trkidiso			(	"HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v(\\d+)", "o"			);
			TPMERegexp regexp_ele8trkidiso_jet30	(	"HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v(\\d+)", "o"	);
			TPMERegexp regexp_ele17trkidiso			(	"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v(\\d+)", "o"			);
			TPMERegexp regexp_ele17trkidiso_jet30	(	"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v(\\d+)", "o"	);

			unsigned int eleFRMask = LeptonTree::QCDFakeEle;
			if (passedTriggerRegExp(regexp_ele8) 				||            		
					passedTriggerRegExp(regexp_ele8trkidiso) 		||            
					passedTriggerRegExp(regexp_ele8trkidiso_jet30) 	||      
					passedTriggerRegExp(regexp_ele17trkidiso) 		||           
					passedTriggerRegExp(regexp_ele17trkidiso_jet30) )     	eventSelection |= LeptonTree::QCDFakeEle;
			if (eventSelection & eleFRMask)                 			MakeElectronFakeRateTree(leptonTree, weight, sample, eventSelection);

			//
			// muon FR
            //
            TPMERegexp regexp_mu8	("HLT_Mu8_v(\\d+)", "o"	);
            TPMERegexp regexp_mu17	("HLT_Mu17_v(\\d+)", "o");

            unsigned int muonFRMask = LeptonTree::QCDFakeMu;
            if (	passedTriggerRegExp(regexp_mu8) || 
					passedTriggerRegExp(regexp_mu17) ) 		eventSelection |= LeptonTree::QCDFakeMu;
            if (eventSelection & muonFRMask)               	MakeMuonFakeRateTree(leptonTree, weight, sample, eventSelection);

            //
            // electron tag and probe
            //

            if (cms2.els_p4().size() > 1)    MakeElectronTagAndProbeTree(leptonTree, weight, sample);

            //
            // muon tag and probe tree
            //

            if (cms2.mus_p4().size() > 1)    MakeMuonTagAndProbeTree(leptonTree, weight, sample);


        }

        //
        // timing information
        //

        t.Stop();
        printf("Finished processing file: %s\n",currentFile->GetTitle());
        printf("Real time: %u events / %f s = %e event/s\n", nEvents, 
                t.RealTime(), nEvents / t.RealTime());
        printf("CPU time: %u events / %f s = %e event/s\n", nEvents, 
                t.CpuTime(), nEvents / t.CpuTime());
        printf("Total duplicate count: %d.  Total weight %f\n",   
                duplicates_total_n, duplicates_total_weight);

        //
        // tidy up
        //

        delete tree;
        delete f;
    }

    //
    // print monitoring information
    //

    if ( nEventsChain != nEventsTotal ) {
        printf("ERROR: number of events from files (%d) is not equal to total number"
                " of events (%d)\n", nEventsChain, nEventsTotal);
    }
    printf("Total number of skipped events due to bad identification: %d (%0.0f %%)\n",   
            nFailedIdentification, nFailedIdentification*100.0/(nEventsChain+1e-5));
    printf("Total number of filtered out events: %d (%0.0f %%)\n",   
            nFilteredOut, nFilteredOut*100.0/(nEventsChain+1e-5));

    //
    // write smurf tree
    // and tidy up
    //

    fSmurf->cd(); 
    leptonTree.tree_->Write();
    fSmurf->Close();

}


void LeptonTreeMaker::MakeElectronRecoTagAndProbeTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample)
{
/*
    // loop on tags
    for (int el = 0; el < cms2.els_p4().size(); ++el) {

        // loop on probes
        for (int sc = 0; sc < cms2.scs_energy().size(); ++sc) {

            // check probe and tag do not overlap
            // this is based on comparing the edm::Ref 
            // in the supercluster and electron supercluster
            if (cms2.scs_elsidx()[sc] == el) continue;

            // basic kinematics
            if (fabs(cms2.scs_eta()[sc]) > 2.5) continue;
            if (cms2.scs_energy()[sc]/cosh(cms2.scs_eta()[sc]) < 10.0) continue;

            // fill the leptonTree
            SetCommonTreeVariables(leptonTree, weight, sample);
            leptonTree.pt_          = cms2.scs_energy()[sc]/cosh(cms2.scs_eta()[sc]);
            leptonTree.eta_         = cms2.scs_eta()[sc];
            leptonTree.phi_         = cms2.scs_phi()[sc];
            leptonTree.q_           = 0;

            // away jet
            //const std::vector<JetPair> &jets = getJets(jetType(), cms2.els_p4()[el], cms2.scs_p4()[sc], 0, 5.0, true, false);
            //leptonTree.leadingAwayJetPt_ =  -999.;
            //if (jets.size() > 0) leptonTree.leadingAwayJetPt_ = jets[0].first.Pt();

            // tag and probe mass
            LorentzVector tag = cms2.els_p4()[el];
            LorentzVector probe = cms2.scs_p4()[sc];
            LorentzVector tp = tag + probe;
            leptonTree.tagAndProbeMass_ = tp.M();

            // lepton selection
            leptonTree.leptonSelection_                                       |= LeptonTree::PassEleSC;
            if (cms2.scs_elsidx()[sc] != -1) leptonTree.leptonSelection_      |= LeptonTree::PassEleReco;

            // event selection
            leptonTree.eventSelection_                                        |= LeptonTree::ZeeTagAndProbe;

            // fill it
            leptonTree.tree_->Fill();
        }

    }
*/
}

void LeptonTreeMaker::MakeElectronTagAndProbeTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample)
{

    // loop on tags
    for (unsigned int tag = 0; tag < cms2.els_p4().size(); ++tag) {
            
		if (cms2.els_p4()[tag].Pt() < 20.)      continue;
        if (fabs(cms2.els_etaSC()[tag]) > 2.5)  continue;

        // check ID criteria on tag
        //if (!goodElectronIsolated(tag, useLHeleId_, useMVAeleId_, electronIdMVA_, lockToCoreSelectors_)) continue;
        if (! goodElectronIsolated(tag, useLHeleId_, useMVAeleId_, egammaMvaEleEstimator_leptree, lockToCoreSelectors_)) continue;

        // loop on probes
        for (unsigned int probe = 0; probe < cms2.els_p4().size(); ++probe) {

            // check probe and tag do not overlap
            if (tag == probe) continue;

            // basic probe denominator
            if (cms2.els_p4()[probe].Pt() < 10.)      continue;
            if (fabs(cms2.els_etaSC()[probe]) > 2.5)  continue;

            // fill the tree - event general variables
            SetCommonTreeVariables(leptonTree, weight, sample);
            leptonTree.eventSelection_ = LeptonTree::ZeeTagAndProbe;

            // fill the tree - probe specific variables
            leptonTree.probe_       = cms2.els_p4()[probe];
            leptonTree.qProbe_      = cms2.els_charge()[probe];
            leptonTree.tag_         = cms2.els_p4()[tag];
            leptonTree.qTag_        = cms2.els_charge()[tag];
            leptonTree.met_         = cms2.evt_pfmet();
            leptonTree.metPhi_      = cms2.evt_pfmetPhi();

            leptonTree.tagAndProbeMass_ = (cms2.els_p4()[probe] + cms2.els_p4()[tag]).M();
           
            HLT_Ele17_Ele8_tag_        = cms2.els_HLT_Ele17_Ele8()[tag];
            HLT_Ele17_Ele8_probe_      = cms2.els_HLT_Ele17_Ele8()[probe];
			HLT_Ele27_WP80_tag_		   = cms2.els_HLT_Ele27_WP80()[tag];
			HLT_Ele27_WP80_probe_	   = cms2.els_HLT_Ele27_WP80()[probe];
			HLT_Ele17_Ele8_tag_					= cms2.els_HLT_Ele17_Ele8()[tag];				
		 	HLT_Ele17_Ele8_probe_				= cms2.els_HLT_Ele17_Ele8()[probe];
		 	HLT_Ele17_Ele8_LeadingLeg_tag_		= cms2.els_HLT_Ele17_Ele8_LeadingLeg()[tag]; 	
		 	HLT_Ele17_Ele8_LeadingLeg_probe_	= cms2.els_HLT_Ele17_Ele8_LeadingLeg()[probe];
		 	HLT_Ele17_Ele8_TrailingLeg_tag_		= cms2.els_HLT_Ele17_Ele8_TrailingLeg()[tag];	
		 	HLT_Ele17_Ele8_TrailingLeg_probe_	= cms2.els_HLT_Ele17_Ele8_TrailingLeg()[probe];

			// fill the tree - criteria the probe passed 
			const std::vector<JetPair> &jets = getJets(jetType(), cms2.els_p4()[tag], cms2.els_p4()[probe], 0, 4.7, true, jet_corrector_pfL1FastJetL2L3_);
			if (jets.size()>0)	leptonTree.jet1_ = jets.at(0).first;
			if (jets.size()>1)	leptonTree.jet2_ = jets.at(1).first;
			if (jets.size()>2)	leptonTree.jet3_ = jets.at(2).first;
			
			// FO
			if (ElectronFOV4(probe))
                leptonTree.leptonSelection_     |= LeptonTree::PassEleFO;
			
            // ID
			if (goodElectronWithoutIsolation(probe, useLHeleId_, useMVAeleId_, egammaMvaEleEstimator_leptree))
                leptonTree.leptonSelection_     |= LeptonTree::PassEleID;

            // ISO
			if (ww_elIso(probe))
                leptonTree.leptonSelection_     |= LeptonTree::PassEleIso;

			// TEST 
			mva_	=	egammaMvaEleEstimator_leptree->mvaValue(probe,false);	
			iso_	=	ww_elIsoVal(probe);	

			fbrem_ 		=	cms2.els_fbrem()[probe]; 
			kfchi2_		=	cms2.els_trkidx()[probe]>=0 ? cms2.trks_chi2()[cms2.els_trkidx()[probe]]/cms2.trks_ndof()[cms2.els_trkidx()[probe]] : 0.;
			kfhits_		= 	cms2.els_trkidx()[probe]>=0 ? cms2.trks_nlayers()[cms2.els_trkidx()[probe]] : -1;
			gsfchi2_	= 	cms2.els_chi2()[probe] / cms2.els_ndof()[probe];
			deta_		=	cms2.els_dEtaIn()[probe];
			dphi_		=	cms2.els_dPhiIn()[probe]; 
			detacalo_	= 	cms2.els_dEtaOut()[probe];
			see_		= 	cms2.els_sigmaIEtaIEta()[probe];
			spp_		=	cms2.els_sigmaIPhiIPhi()[probe]; // FIXME : check the case where it's 0 
			etawidth_	=	cms2.els_etaSCwidth()[probe];
			phiwidth_	= 	cms2.els_phiSCwidth()[probe];
			e1x5e5x5_	=	cms2.els_e5x5()[probe]!=0. ? 1. - cms2.els_e1x5()[probe]/cms2.els_e5x5()[probe] : -1; 
			R9_			= 	cms2.els_e3x3()[probe] / cms2.els_eSCRaw()[probe];
			HoE_	 			=	cms2.els_hOverE()[probe];
			EoP_				=	cms2.els_eOverPIn()[probe];
			IoEmIoP_			=	1./cms2.els_ecalEnergy()[probe] - 1./cms2.els_p4()[probe].P(); // this is consistent with CMSSW 
			eleEoPout_			=	cms2.els_eOverPOut()[probe];
			PreShowerOverRaw_	=	cms2.els_eSCPresh()[probe] / cms2.els_eSCRaw()[probe];
			d0_					=	electron_d0PV_wwV1_local(probe);
			const double gsfsign = ( (gsftrks_d0_pv(cms2.els_gsftrkidx().at(probe),0).first)   >=0 ) ? 1. : -1.;
			ip3d_				=	cms2.els_ip3d()[probe]*gsfsign; 
			eta_				= 	cms2.els_etaSC()[probe];
			pt_					= 	cms2.els_p4()[probe].pt();
          	 
			// fill it
            leptonTree.tree_->Fill();

        }

    }
}

void LeptonTreeMaker::MakeMuonTagAndProbeTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample) {
	
	std::vector<Int_t> nullMu; // null identified muons // FIXME
	std::vector<Int_t> nullEle; // null identified electrons  // FIXME

    // loop on tags
    for (unsigned int tag = 0; tag < cms2.mus_p4().size(); ++tag) {
        
		if (cms2.mus_p4()[tag].Pt() < 20.0)                 continue; // pt cut
        if (fabs(cms2.mus_p4()[tag].Eta()) > 2.4)           continue; // eta cut

        // check ID criteria on tag
        if (!goodMuonIsolated(tag, lockToCoreSelectors_, false, muonIdMVA_leptree, muonMVAEstimator_leptree,  nullMu, nullEle) ) continue;

        // loop on probes
        for (unsigned int probe = 0; probe < cms2.mus_p4().size(); ++probe) {

            // check probe and tag do not overlap
            if (tag == probe) continue;

            // basic probe denominator
            if (cms2.mus_p4()[probe].Pt() < 10.0)                 continue; // pt cut
            if (fabs(cms2.mus_p4()[probe].Eta()) > 2.4)          continue; // eta cut
            if (((cms2.mus_type()[probe]) & (1<<2)) == 0)        continue; // tracker muon

            // fill the tree - event general variables
            SetCommonTreeVariables(leptonTree, weight, sample);
            leptonTree.eventSelection_ = LeptonTree::ZmmTagAndProbe;

            // fill the tree - probe specific variables
            leptonTree.probe_       = cms2.mus_p4()[probe];
            leptonTree.qProbe_      = cms2.mus_charge()[probe];
            leptonTree.tag_         = cms2.mus_p4()[tag];
            leptonTree.qTag_        = cms2.mus_charge()[tag];
            leptonTree.met_         = cms2.evt_pfmet();
            leptonTree.metPhi_      = cms2.evt_pfmetPhi();

            leptonTree.tagAndProbeMass_ = (cms2.mus_p4()[probe] + cms2.mus_p4()[tag]).M();

			HLT_IsoMu24_eta2p1_tag_ 	=	cms2.mus_HLT_IsoMu24_eta2p1()[tag];
			HLT_IsoMu24_eta2p1_probe_ 	=	cms2.mus_HLT_IsoMu24_eta2p1()[probe];
			HLT_Mu17_Mu8_tag_			= 	cms2.mus_HLT_Mu17_Mu8()[tag];
			HLT_Mu17_Mu8_probe_			= 	cms2.mus_HLT_Mu17_Mu8()[probe];
			HLT_Mu17_Mu8_LeadingLeg_tag_		=   cms2.mus_HLT_Mu17_Mu8_LeadingLeg()[tag];		
		 	HLT_Mu17_Mu8_LeadingLeg_probe_		=   cms2.mus_HLT_Mu17_Mu8_LeadingLeg()[probe];		
		 	HLT_Mu17_Mu8_TrailingLeg_tag_		=   cms2.mus_HLT_Mu17_Mu8_TrailingLeg()[tag];	
		 	HLT_Mu17_Mu8_TrailingLeg_probe_		=   cms2.mus_HLT_Mu17_Mu8_TrailingLeg()[probe];
			HLT_Mu17_TkMu8_tag_					=	cms2.mus_HLT_Mu17_TkMu8()[tag];
			HLT_Mu17_TkMu8_probe_				=	cms2.mus_HLT_Mu17_TkMu8()[probe];
		 	HLT_Mu17_TkMu8_LeadingLeg_tag_		=   cms2.mus_HLT_Mu17_TkMu8_LeadingLeg()[tag];	
		 	HLT_Mu17_TkMu8_LeadingLeg_probe_	=   cms2.mus_HLT_Mu17_TkMu8_LeadingLeg()[probe];	
		 	HLT_Mu17_TkMu8_TrailingLeg_tag_		=   cms2.mus_HLT_Mu17_TkMu8_TrailingLeg()[tag];	
		 	HLT_Mu17_TkMu8_TrailingLeg_probe_	=   cms2.mus_HLT_Mu17_TkMu8_TrailingLeg()[probe];

            // fill the tree - criteria the probe passed 
			const std::vector<JetPair> &jets = getJets(jetType(), cms2.mus_p4()[tag], cms2.mus_p4()[probe], 0, 4.7, true, jet_corrector_pfL1FastJetL2L3_);
			if (jets.size()>0)	leptonTree.jet1_ = jets.at(0).first;
			if (jets.size()>1)	leptonTree.jet2_ = jets.at(1).first;
			if (jets.size()>2)	leptonTree.jet3_ = jets.at(2).first;

            // FO
        	if ( fakableMuon(probe, MuFOV2, muonMVAEstimator_leptree,  nullMu, nullEle) ) 
                	leptonTree.leptonSelection_     |= LeptonTree::PassMuFO;

            // ID
        	if ( goodMuonWithoutIsolation(probe,false, muonIdMVA_leptree, muonMVAEstimator_leptree, nullMu, nullEle) )
                	leptonTree.leptonSelection_     |= LeptonTree::PassMuID;

            // ISO
        	if ( ww_muIso(probe, muonMVAEstimator_leptree, nullMu, nullEle) )
                	leptonTree.leptonSelection_     |= LeptonTree::PassMuIso;

			mva_ 	=	muonMVAEstimator_leptree->mvaValueIso( probe, cms2.evt_ww_rho(), MuonEffectiveArea::kMuEAFall11MC,
			                                                     nullMu, nullEle, false );
            // fill it
            leptonTree.tree_->Fill();

        }

    }
}

//
// fake rates
//

void LeptonTreeMaker::MakeMuonFakeRateTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample, const unsigned int eventSelection)
{
	std::vector<Int_t> nullMu; // null identified muons // FIXME
	std::vector<Int_t> nullEle; // null identified electrons  // FIXME

    //
    // count FO
    //

    unsigned int fo = 0;
    unsigned int nfo = 0;
    for (unsigned int i = 0; i < cms2.mus_p4().size(); ++i) {
        if ( fakableMuon(i, MuFOV2, muonMVAEstimator_leptree,  nullMu, nullEle) ) { ++nfo; fo = i; }
    }

    //
    // require exactly one FO
    //

    if (nfo == 1) {

        //
        // check event properties
        //

        if (metValue() > 20.0) return;

        // 
        // fill the tree
        //

        SetCommonTreeVariables(leptonTree, weight, sample);

        // fill the tree - probe specific variables
        leptonTree.probe_       = cms2.mus_p4()[fo];
        leptonTree.qProbe_      = cms2.mus_charge()[fo];
        leptonTree.met_         = cms2.evt_pfmet();
        leptonTree.metPhi_      = cms2.evt_pfmetPhi();

        LorentzVector nullvector(0.0, 0.0, 0.0, 0.0);

        leptonTree.eventSelection_ = eventSelection;

		// FR triggers
        TPMERegexp regexp_mu8	("HLT_Mu8_v(\\d+)", "o"	);
        TPMERegexp regexp_mu17	("HLT_Mu17_v(\\d+)", "o");

        if(	passedTriggerRegExp(regexp_mu8)	)  	HLT_Mu8_ = true;
		if(	passedTriggerRegExp(regexp_mu17)) 	HLT_Mu17_ = true;

       	if ( fakableMuon(fo, MuFOV2, muonMVAEstimator_leptree,  nullMu, nullEle) ) 
          	leptonTree.leptonSelection_     |= LeptonTree::PassMuFO;
        if (goodMuonWithoutIsolation(fo,false, muonIdMVA_leptree, muonMVAEstimator_leptree, nullMu, nullEle))
            leptonTree.leptonSelection_     |= LeptonTree::PassMuID;
        if (ww_muIso(fo, muonMVAEstimator_leptree, nullMu,  nullEle))
            leptonTree.leptonSelection_     |= LeptonTree::PassMuIso;

    	const std::vector<JetPair> &jets = getJets(jetType(), cms2.mus_p4()[fo], nullvector, 0, 4.7, true, jet_corrector_pfL1FastJetL2L3_);
		if (jets.size()>0)	leptonTree.jet1_ = jets.at(0).first;
		if (jets.size()>1)	leptonTree.jet2_ = jets.at(1).first;
		if (jets.size()>2)	leptonTree.jet3_ = jets.at(2).first;

        leptonTree.tree_->Fill(); 

    }

}

void LeptonTreeMaker::MakeElectronFakeRateTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample, const unsigned int eventSelection)
{

	//
	// count FO
	//

	// 1 and only 1 very loose
	unsigned int fo = 0;
	unsigned int nfo = 0;
	for (unsigned int i = 0; i < cms2.els_p4().size(); ++i) {
		if (ElectronFOV4(i)) { 	++nfo; fo = i; }
	}

	//
	// require exactly one FO
	//

	if (nfo == 1) {

			//
			// check event properties
			//

			if (metValue() > 20.0) return;

			// 
			// fill the tree
			//

			SetCommonTreeVariables(leptonTree, weight, sample);

			// fill the tree - probe specific variables
			leptonTree.probe_       = cms2.els_p4()[fo];
			leptonTree.qProbe_      = cms2.els_charge()[fo];
        	leptonTree.met_         = cms2.evt_pfmet();
        	leptonTree.metPhi_      = cms2.evt_pfmetPhi();
			LorentzVector nullvector(0.0, 0.0, 0.0, 0.0);
			leptonTree.eventSelection_ = eventSelection;

			// FR triggers
            TPMERegexp regexp_ele8					(	"HLT_Ele8_CaloIdL_TrkIdVL_v(\\d+)", "o"								);
            TPMERegexp regexp_ele8trkidiso			(	"HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v(\\d+)", "o"			);
            TPMERegexp regexp_ele8trkidiso_jet30	(	"HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v(\\d+)", "o"	);
            TPMERegexp regexp_ele17trkidiso			(	"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v(\\d+)", "o"			);
            TPMERegexp regexp_ele17trkidiso_jet30	(	"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v(\\d+)", "o"	);

            if( passedTriggerRegExp(regexp_ele8) )			           	HLT_Ele8_CaloIdL_TrkIdVL_ = true;
            if(	passedTriggerRegExp(regexp_ele8trkidiso) )		        HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_ = true;
            if(	passedTriggerRegExp(regexp_ele8trkidiso_jet30)) 		HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_ = true;
           	if(	passedTriggerRegExp(regexp_ele17trkidiso) )				HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_ = true;
            if(	passedTriggerRegExp(regexp_ele17trkidiso_jet30) )     	HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_ = true;

			const std::vector<JetPair> &jets = getJets(jetType(), cms2.els_p4()[fo], nullvector, 0, 4.7, true, jet_corrector_pfL1FastJetL2L3_);
			if (jets.size()>0)	leptonTree.jet1_ = jets.at(0).first;
			if (jets.size()>1)  leptonTree.jet2_ = jets.at(1).first;
			if (jets.size()>2)	leptonTree.jet3_ = jets.at(2).first;
			
			if (ElectronFOV4(fo))
			    leptonTree.leptonSelection_     |= LeptonTree::PassEleFO;
			if (goodElectronWithoutIsolation(fo, useLHeleId_, useMVAeleId_, egammaMvaEleEstimator_leptree))
				leptonTree.leptonSelection_     |= LeptonTree::PassEleID;
			if (ww_elIso(fo))
				leptonTree.leptonSelection_     |= LeptonTree::PassEleIso;
			leptonTree.tree_->Fill();

		}

}


void LeptonTreeMaker::SetCommonTreeVariables(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample)
{

    leptonTree.InitVariables();
    leptonTree.run_         = cms2.evt_run();
    leptonTree.event_       = cms2.evt_event();
    leptonTree.lumi_        = cms2.evt_lumiBlock();
    leptonTree.nvtx_        = nGoodVertex();
    leptonTree.scale1fb_    = weight;
    //leptonTree.rho_         = cms2.evt_ww_rho_vor();
  	if (sample!=SmurfTree::data && sample!=SmurfTree::dyttDataDriven){
    	for (unsigned int nbc=0;nbc<cms2.puInfo_nPUvertices().size();++nbc) {
		      if (cms2.puInfo_bunchCrossing().at(nbc)==0) leptonTree.npu_ = cms2.puInfo_nPUvertices().at(nbc);
    	} 
	}
}

float LeptonTreeMaker::GetAwayJetPt(LorentzVector lep1, LorentzVector lep2)
{
    const std::vector<JetPair> &jets = getJets(jetType(), lep1, lep2, 0, 4.7, true, jet_corrector_pfL1FastJetL2L3_);
    if (jets.size() > 0 && TMath::Abs(ROOT::Math::VectorUtil::DeltaR(lep1,jets[0].first))>1.0) return jets[0].first.Pt();
    else return -999.9;
}

