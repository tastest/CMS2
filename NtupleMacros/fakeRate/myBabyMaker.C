#include <iostream>
#include <set>

// ROOT includes
#include "TSystem.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TChainElement.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Math/VectorUtil.h"

// TAS includes
#include "./CMS2.cc"
#include "../CORE/trackSelections.cc"
#include "../CORE/eventSelections.cc"

#include "../CORE/muonSelections.cc"
#include "../CORE/electronSelections.cc"
#include "../CORE/MITConversionUtilities.cc"
#include "../CORE/electronSelectionsParameters.cc"
#include "../CORE/metSelections.cc"

#include "../CORE/triggerUtils.cc"
#include "../Tools/goodrun.cc"
#include "./myBabyMaker.h"
#include "../CORE/jetSelections.cc"

using namespace std;
using namespace tas;

// function for dR matching offline letpon to trigger object 
pair<int, float> TriggerMatch( LorentzVector lepton_p4, const char* trigString, double dR_cut = 0.4 ){
  float dR_min = numeric_limits<float>::max();
  dR_min = 99.0;
  int nTrig = nHLTObjects( trigString );
  if (nTrig > 0) {
    bool match = false;
    for (int itrg=0; itrg<nTrig; itrg++) {
      LorentzVector p4tr = p4HLTObject( trigString, itrg );
      double dr = ROOT::Math::VectorUtil::DeltaR( lepton_p4, p4tr);
      if (dr < dR_cut) match = true;
      if (dr < dR_min) dR_min = dr;
    }
    if (match) {
      nTrig = 2;
    } 
    else {
      nTrig = 1;
    }
  }
  pair<int, float> answer;
  answer.first  = nTrig;
  answer.second = dR_min;
  return answer;
}



/* THIS NEEDS TO BE IN CORE */

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

// transverse mass
float Mt( LorentzVector p4, float met, float met_phi ){
  return sqrt( 2*met*( p4.pt() - ( p4.Px()*cos(met_phi) + p4.Py()*sin(met_phi) ) ) );
}

//-----------------------------------
// Looper code starts here
// eormu=-1 do both e and mu
//      =11 do electrons
//      =13 do muons
//-----------------------------------
void myBabyMaker::ScanChain( TChain* chain, const char *babyFilename, bool isData, int eormu) {

  already_seen.clear();

  // Make a baby ntuple
  MakeBabyNtuple(babyFilename);

  // Set the JSON file
  if(isData){
    set_goodrun_file_json("json/Cert_TopNov5_Merged_135821-149442_allPVT.txt");
    //set_goodrun_file_json("json/json_135821_148058_15.21pb.txt");
    //set_goodrun_file_json("json/Cert_TopAug13_Merged_135059-142664.txt");
  }

  std::vector<std::string> jetcorr_filenames;
  jetcorr_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L2Relative_AK5PF.txt");
  jetcorr_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L3Absolute_AK5PF.txt");
  FactorizedJetCorrector *jet_corrector = makeJetCorrector(jetcorr_filenames);

  // The deltaR requirement between objects and jets to remove the jet trigger dependence
  float deltaRCut   = 1.0;
  float deltaPhiCut = 2.5;

  //--------------------------
  // File and Event Loop
  //---------------------------
  int i_permilleOld = 0;
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = 0;
  int nEvents = -1;
  if (nEvents==-1){
    nEventsChain = chain->GetEntries();
  } else {
    nEventsChain = nEvents;
  }
  nEventsChain = chain->GetEntries();
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  map<int,int> m_events;
  while(TChainElement *currentFile = (TChainElement*)fileIter.Next() ) {
    TString filename = currentFile->GetTitle();
    
    TFile f(filename.Data());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    unsigned int nEntries = tree->GetEntries();
    unsigned int nLoop = nEntries;
    unsigned int z;
    for( z = 0; z < nLoop; z++) { // Event Loop
      cms2.GetEntry(z);

      if(isData){
        // Good  Runs
        //if(!goodrun( evt_run(), evt_lumiBlock() )) continue;
        if(!goodrun_json( evt_run(), evt_lumiBlock() )) continue;

        // check for duplicated
        DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
        if (is_duplicate(id) ){ 
          cout << "\t! ERROR: found duplicate." << endl;
          continue;
        }
      }

      // looper progress
      ++nEventsTotal;
      int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
      if (i_permille != i_permilleOld) {
        printf("  \015\033[32m ---> \033[1m\033[31m%4.1f%%" "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
        fflush(stdout);
        i_permilleOld = i_permille;
      }
      
      // Event cleaning (careful, it requires technical bits)
      //if (!cleaning_BPTX(isData))   continue;
      //if (!cleaning_beamHalo())   continue;
      if (!cleaning_goodVertexAugust2010()) continue;
      if (!cleaning_goodTracks()) continue;

      // Loop over jets and see what is btagged
      // Medium operating point from https://twiki.cern.ch/twiki/bin/view/CMS/BTagPerformanceOP
      int this_nbjet = 0;
      vector<unsigned int> bindex;
      for (unsigned int iJet = 0; iJet < jets_p4().size(); iJet++) {
        if (jets_p4().at(iJet).pt() < 15.) continue;
        if (jets_simpleSecondaryVertexHighEffBJetTag().at(iJet) < 1.74) continue;
        this_nbjet++;
        bindex.push_back(iJet);
      }

// PF Jets

      int this_nbpfjet = 0;
      vector<unsigned int> bpfindex;
      for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
        if ( !passesPFJetID(iJet)) continue;
        LorentzVector jp4 = pfjets_p4()[iJet];
        float jet_cor = jetCorrection(jp4, jet_corrector);
        LorentzVector jp4cor = jp4 * jet_cor;
        if (jp4cor.pt() < 15) continue;
        if (pfjets_simpleSecondaryVertexHighEffBJetTag().at(iJet) < 1.74) continue;
        this_nbpfjet++;
        bpfindex.push_back(iJet);
      }

      
/* Electrons */
      
      if (eormu == -1 || eormu==11) {
      for (unsigned int iEl = 0 ; iEl < els_p4().size(); iEl++) {

        // Apply a pt cut (Changed it from 5 GeV to 10 GeV...Claudio 10 July 2010)
        if ( els_p4().at(iEl).pt() < 10.) continue;

        // Initialize baby ntuple
        InitBabyNtuple();

        // Add spike veto
        num_ = pass_electronSelection( iEl, electronSelection_ttbarV1 ) && (!isSpikeElectron(iEl));
        numv1_ = pass_electronSelection( iEl, electronSelection_ttbarV1 );
        v1_  = pass_electronSelection( iEl, electronSelectionFO_el_ttbarV1_v1 );
        v2_  = pass_electronSelection( iEl, electronSelectionFO_el_ttbarV1_v2 );
        v3_  = pass_electronSelection( iEl, electronSelectionFO_el_ttbarV1_v3 );

        numSS_ = pass_electronSelection(iEl, electronSelection_ss);
        v1SS_  = pass_electronSelection(iEl, electronSelectionFO_ssVBTF80_v1);
        v2SS_  = pass_electronSelection(iEl, electronSelectionFO_ssVBTF80_v2);
        v3SS_  = pass_electronSelection(iEl, electronSelectionFO_ssVBTF80_v3);

        numAug9_ = pass_electronSelection( iEl, electronSelection_ttbarV1, isData, true ) && (!isSpikeElectron(iEl));
        v1Aug9_  = v1_;
        v2Aug9_  = v2_;
        v3Aug9_  = pass_electronSelection( iEl, electronSelectionFO_el_ttbarV1_v3, isData, true );

        numSSAug9_ = pass_electronSelection(iEl, electronSelection_ss, isData, true);
        v1SSAug9_  = pass_electronSelection(iEl, electronSelectionFO_ssVBTF80_v1, isData, true);
        v2SSAug9_  = pass_electronSelection(iEl, electronSelectionFO_ssVBTF80_v2, isData, true);
        v3SSAug9_  = pass_electronSelection(iEl, electronSelectionFO_ssVBTF80_v3, isData, true);

	numOct6_ = pass_electronSelection( iEl, electronSelection_ttbarV1_pass5);
	v1Oct6_  = pass_electronSelection( iEl, electronSelectionFO_el_ttbarV1_v1_pass5);
	v2Oct6_  = pass_electronSelection( iEl, electronSelectionFO_el_ttbarV1_v2_pass5);
	v3Oct6_  = pass_electronSelection( iEl, electronSelectionFO_el_ttbarV1_v3_pass5);

	numOSOct18_ = pass_electronSelection( iEl, electronSelection_el_OSV1);
	v1OSOct18_  = pass_electronSelection( iEl, electronSelectionFO_el_OSV1_v1);
	v2OSOct18_  = pass_electronSelection( iEl, electronSelectionFO_el_OSV1_v2);
	v3OSOct18_  = pass_electronSelection( iEl, electronSelectionFO_el_OSV1_v3);

	numSSOct18_ = pass_electronSelection( iEl, electronSelection_ss, false, false) && (!isSpikeElectron(iEl));
	v1SSOct18_  = pass_electronSelection( iEl, electronSelectionFO_ssVBTF80_v1, false, false);
	v2SSOct18_  = pass_electronSelection( iEl, electronSelectionFO_ssVBTF80_v2, false, false);
	v3SSOct18_  = pass_electronSelection( iEl, electronSelectionFO_ssVBTF80_v3, false, false);

        v1_wwV1_  = pass_electronSelection( iEl, electronSelectionFO_el_wwV1_v1);
        v2_wwV1_  = pass_electronSelection( iEl, electronSelectionFO_el_wwV1_v2);
        v3_wwV1_  = pass_electronSelection( iEl, electronSelectionFO_el_wwV1_v3);
        v4_wwV1_  = pass_electronSelection( iEl, electronSelectionFO_el_wwV1_v4);
        num_wwV1_ = pass_electronSelection( iEl, electronSelection_wwV1);

        numSSV2_ = pass_electronSelection( iEl, electronSelection_ssV2, false, false) && (!isSpikeElectron(iEl));
        v1SSV2_  = pass_electronSelection( iEl, electronSelectionFOV2_ssVBTF80_v1, false, false);
        v2SSV2_  = pass_electronSelection( iEl, electronSelectionFOV2_ssVBTF80_v2, false, false);
        v3SSV2_  = pass_electronSelection( iEl, electronSelectionFOV2_ssVBTF80_v3, false, false);


        // Sanity
        if (numOSOct18_ && (!v1OSOct18_)) cout << "bad v1OSOct18_" << endl;
        if (numOSOct18_ && (!v2OSOct18_)) cout << "bad v2OSOct18_" << endl;
        if (numOSOct18_ && (!v3OSOct18_)) cout << "bad v3OSOct18_" << endl;

        if (numSSOct18_ && (!v1SSOct18_)) cout << "bad v1SSOct18_" << endl;
        if (numSSOct18_ && (!v2SSOct18_)) cout << "bad v2SSOct18_" << endl;
        if (numSSOct18_ && (!v3SSOct18_)) cout << "bad v3SSOct18_" << endl;

        if (num_ && (!v1_)) cout << "bad v1_" << endl;
        if (num_ && (!v2_)) cout << "bad v2_" << endl;
        if (num_ && (!v3_)) cout << "bad v3_" << endl;

        if (numSS_ && (!v1SS_)) cout << "bad v1SS_" << endl;
        if (numSS_ && (!v2SS_)) cout << "bad v2SS_" << endl;
        if (numSS_ && (!v3SS_)) cout << "bad v3SS_" << endl;

        if (numSSAug9_ && (!v1SSAug9_)) cout << "bad v1SSAug9_" << endl;
        if (numSSAug9_ && (!v2SSAug9_)) cout << "bad v2SSAug9_" << endl;
        if (numSSAug9_ && (!v3SSAug9_)) cout << "bad v3SSAug9_" << endl;

        if (numAug9_ && (!v1Aug9_)) cout << "bad v1Aug9_" << endl;
        if (numAug9_ && (!v2Aug9_)) cout << "bad v2Aug9_" << endl;
        if (numAug9_ && (!v3Aug9_)) cout << "bad v3Aug9_" << endl;

        if (numOct6_ && (!v1Oct6_)) cout << "bad v1Oct6_" << endl;
        if (numOct6_ && (!v2Oct6_)) cout << "bad v2Oct6_" << endl;
        if (numOct6_ && (!v3Oct6_)) cout << "bad v3Oct6_" << endl;

        if (num_wwV1_ && (!v1_wwV1_)) cout << "bad v1_wwV1_" << endl;
        if (num_wwV1_ && (!v2_wwV1_)) cout << "bad v2_wwV1_" << endl;
        if (num_wwV1_ && (!v3_wwV1_)) cout << "bad v3_wwV1_" << endl;
        if (num_wwV1_ && (!v4_wwV1_)) cout << "bad v4_wwV1_" << endl;

        if (numSSV2_ && (!v1SSV2_)) cout << "bad v1SSV2_" << endl;
        if (numSSV2_ && (!v2SSV2_)) cout << "bad v2SSV2_" << endl;
        if (numSSV2_ && (!v3SSV2_)) cout << "bad v3SSV2_" << endl;



        // If there is no v1/v2/v3 lepton quit
        if (  (!v1_) && (!v2_) && (!v3_) && 
              (!v1SS_) && (!v2SS_) && (!v3SS_) && 
              (!v1Aug9_) && (!v2Aug9_) && (!v3Aug9_) &&
              (!v1SSAug9_) && (!v2SSAug9_) && (!v3SSAug9_) &&
              (!v1OSOct18_) && (!v2OSOct18_) && (!v3OSOct18_) &&
              (!v1SSOct18_) && (!v2SSOct18_) && (!v3SSOct18_) &&
              (!v1SSV2_) && (!v2SSV2_) && (!v3SSV2_) &&
              (!v1Oct6_) && (!v2Oct6_) && (!v3Oct6_) &&
              (!v1_wwV1_) && (!v2_wwV1_) && (!v3_wwV1_) && (!v4_wwV1_)
        ) continue;
        
        // If it is above 20 GeV see if we can make a 
        // Z with another pt>20 FO.  Will use the v1 FO since 
        // these are the loosest
        bool isaZ = false;
        if (els_p4().at(iEl).pt() > 20.) {
          for (unsigned int jEl = 0 ; jEl < els_p4().size(); jEl++) {
            if (iEl == jEl)                             continue;
            if (els_p4().at(jEl).pt() < 20.)            continue;
            if ( ! pass_electronSelection( jEl, electronSelectionFO_el_ttbarV1_v1 ) ) continue;
            if ( ! v1_ ) continue;
            LorentzVector w = els_p4().at(iEl) + els_p4().at(jEl);
            if (abs(w.mass()-91.) > 20.) continue;
            isaZ = true;
          }
        }
        if (isaZ) continue;
    
        // Load the electron and event quantities
        run_   = evt_run();
        ls_    = evt_lumiBlock();
        evt_   = evt_event();
		weight_ = isData ? 1. : evt_scale1fb();
        pt_    = els_p4().at(iEl).pt();
        eta_   = els_p4().at(iEl).eta();
        phi_   = els_p4().at(iEl).phi();
        scet_  = els_eSC()[iEl] / cosh( els_etaSC()[iEl] );
        id_    = 11*els_charge().at(iEl);
        tcmet_ = evt_tcmet();
        tcmetphi_ = evt_tcmetPhi();
        pfmet_ = evt_pfmet();
        pfmetphi_ = evt_pfmetPhi();
        iso_ = electronIsolation_rel(iEl, true);

        // Pileup - PUSummaryInfoMaker
        pu_nPUvertices_ = puInfo_nPUvertices();
        pu_zpositions_  = puInfo_zpositions();
        pu_sumptlowpt_  = puInfo_sumpt_lowpt();
        pu_sumpthighpt_ = puInfo_sump_highpt();
        pu_instLumi_    = puInfo_instLumi();
        pu_ntrkslowpt_  = puInfo_ntrks_lowpt();
        pu_ntrkshighpt_ = puInfo_ntrks_highpt();

        // Pileup - VertexMaker
        evt_nvtxs_       = evt_nvtxs();
        vtxs_xError_     = vtxs_xError();
        vtxs_yError_     = vtxs_yError();
        vtxs_zError_     = vtxs_zError();
        vtxs_chi2_       = vtxs_chi2();
        vtxs_ndof_       = vtxs_ndof();
        vtxs_sumpt_      = vtxs_sumpt();
        vtxs_isFake_     = vtxs_isFake();
        vtxs_isValid_    = vtxs_isValid();
        vtxs_tracksSize_ = vtxs_tracksSize();
        vtxs_covMatrix_  = vtxs_covMatrix();
        vtxs_position_   = vtxs_position();

        // Pileup - VertexMaker
        evt_ndavtxs_       = evt_ndavtxs();
        davtxs_xError_     = davtxs_xError();
        davtxs_yError_     = davtxs_yError();
        davtxs_zError_     = davtxs_zError();
        davtxs_chi2_       = davtxs_chi2();
        davtxs_ndof_       = davtxs_ndof();
        davtxs_sumpt_      = davtxs_sumpt();
        davtxs_isFake_     = davtxs_isFake();
        davtxs_isValid_    = davtxs_isValid();
        davtxs_tracksSize_ = davtxs_tracksSize();
        davtxs_covMatrix_  = davtxs_covMatrix();
        davtxs_position_   = davtxs_position();



        if (! isData) {
            mcid_       = els_mc_id().at(iEl);
            mcmotherid_ = els_mc_motherid().at(iEl);
        }
    
        // do the 3 electron charges agree?
        int iCTF = els_trkidx().at(iEl);
        if( iCTF >= 0 ){
          int qCTF = trks_charge().at( iCTF );
          int qGSF = els_trk_charge().at(iEl);
          int qPIX = els_sccharge().at(iEl);
          if( qCTF == qGSF && qCTF == qPIX && qGSF == qPIX ) q3_ = true;
        }

       // Missing hits info
       // Warning els_exp_innerlayers39X_ is set to 999 if this branch doesn't exits

        els_exp_innerlayers_ = els_exp_innerlayers().at(iEl);
        //els_exp_innerlayers39X_ = els_exp_innerlayers39X().at(iEl);

        // W transverse mass
        mt_ = Mt( els_p4().at(iEl), tcmet_, tcmetphi_ );
        pfmt_ = Mt( els_p4().at(iEl), pfmet_, pfmetphi_ );
        
        // The btag information
        nbjet_ = this_nbjet;
        dRbNear_ = 99.;
        dRbFar_  = -99.;
        for (int ii=0; ii<nbjet_; ii++) {
          unsigned int iJet = bindex[ii];
          float dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), jets_p4().at(iJet));
          if (dr < dRbNear_) dRbNear_ = dr;
          if (dr > dRbFar_)   dRbFar_  = dr;
        }

        // btag info for corrected pfjet

        nbpfcjet_ = this_nbpfjet;
        dRbpfcNear_ = 99.;
        dRbpfcFar_  = -99.;
        for (int ii=0; ii<nbpfcjet_; ii++) {
          unsigned int iJet = bpfindex[ii];
          LorentzVector jp4 = pfjets_p4()[iJet];
          float jet_cor = jetCorrection(jp4, jet_corrector);
          LorentzVector jp4cor = jp4 * jet_cor;
          float dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), jp4cor);
          if (dr < dRbpfcNear_) dRbpfcNear_ = dr;
          if (dr > dRbpfcFar_)   dRbpfcFar_  = dr; 
        }
          
 
        ///////////////////////  
        // 2011 Triggers     //
        ///////////////////////

          // Electrons
          pair<int, float> pair_ele8_v2                                           = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_v2"                                         );
          pair<int, float> pair_ele8_CaloIdL_TrkIdVL_v2                           = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_CaloIdL_TrkIdVL_v2"                         );
          pair<int, float> pair_ele8_CaloIdL_CaloIsoVL_Jet40_v2                   = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2"                 );
          pair<int, float> pair_ele8_CaloIdL_CaloIsoVL_v2                         = TriggerMatch( els_p4().at(iEl), "HLT_Ele8_CaloIdL_CaloIsoVL_v2"                       );
          pair<int, float> pair_ele17_CaloIdL_CaloIsoVL_v2                        = TriggerMatch( els_p4().at(iEl), "HLT_Ele17_CaloIdL_CaloIsoVL_v2"                      );  
          //pair<int, float> pair_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2  = TriggerMatch( els_p4().at(iEl), "HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2");
  
          ele8_v2_                                              = pair_ele8_v2.first;
          ele8_CaloIdL_TrkIdVL_v2_                              = pair_ele8_CaloIdL_TrkIdVL_v2.first; 
          ele8_CaloIdL_CaloIsoVL_Jet40_v2_                      = pair_ele8_CaloIdL_CaloIsoVL_Jet40_v2.first;
          ele8_CaloIdL_CaloIsoVL_v2_                            = pair_ele8_CaloIdL_CaloIsoVL_v2.first;
          ele17_CaloIdL_CaloIsoVL_v2_                           = pair_ele17_CaloIdL_CaloIsoVL_v2.first; 
          //photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_     = pair_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2.first;
          
          dr_ele8_v2_                                           = pair_ele8_v2.second;
          dr_ele8_CaloIdL_TrkIdVL_v2_                           = pair_ele8_CaloIdL_TrkIdVL_v2.second;
          dr_ele8_CaloIdL_CaloIsoVL_Jet40_v2_                   = pair_ele8_CaloIdL_CaloIsoVL_Jet40_v2.second;
          dr_ele8_CaloIdL_CaloIsoVL_v2_                         = pair_ele8_CaloIdL_CaloIsoVL_v2.second;
          dr_ele17_CaloIdL_CaloIsoVL_v2_                        = pair_ele17_CaloIdL_CaloIsoVL_v2.second;
          //dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2_  = pair_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2.second; 

          // Muons
          pair<int, float> pair_mu3_v3       = TriggerMatch( els_p4().at(iEl), "HLT_Mu3_v3"       );
          pair<int, float> pair_mu5_v3       = TriggerMatch( els_p4().at(iEl), "HLT_Mu5_v3"       );
          pair<int, float> pair_mu8_v1       = TriggerMatch( els_p4().at(iEl), "HLT_Mu8_v1"       );
          pair<int, float> pair_mu12_v1      = TriggerMatch( els_p4().at(iEl), "HLT_Mu12_v1"      );
          pair<int, float> pair_mu15_v2      = TriggerMatch( els_p4().at(iEl), "HLT_Mu15_v2"      );
          pair<int, float> pair_mu20_v1      = TriggerMatch( els_p4().at(iEl), "HLT_Mu20_v1"      );
          pair<int, float> pair_mu24_v1      = TriggerMatch( els_p4().at(iEl), "HLT_Mu24_v1"      );
          pair<int, float> pair_mu30_v1      = TriggerMatch( els_p4().at(iEl), "HLT_Mu30_v1"      );
          pair<int, float> pair_mu8_Jet40_v3 = TriggerMatch( els_p4().at(iEl), "HLT_Mu8_Jet40_v3" );

          mu3_v3_         = pair_mu3_v3.first;
          mu5_v3_         = pair_mu5_v3.first;
          mu8_v1_         = pair_mu8_v1.first;
          mu12_v1_        = pair_mu12_v1.first;
          mu15_v2_        = pair_mu15_v2.first;
          mu20_v1_        = pair_mu20_v1.first;
          mu24_v1_        = pair_mu24_v1.first;
          mu30_v1_        = pair_mu30_v1.first;
          mu8_Jet40_v3_   = pair_mu8_Jet40_v3.first;

          dr_mu3_v3_       = pair_mu3_v3.second;                                                
          dr_mu5_v3_       = pair_mu5_v3.second;                                                
          dr_mu8_v1_       = pair_mu8_v1.second;                                                
          dr_mu12_v1_      = pair_mu12_v1.second;                                               
          dr_mu15_v2_      = pair_mu15_v2.second;                                               
          dr_mu20_v1_      = pair_mu20_v1.second;  
          dr_mu24_v1_      = pair_mu24_v1.second;
          dr_mu30_v1_      = pair_mu30_v1.second;                                               
          dr_mu8_Jet40_v3_ = pair_mu8_Jet40_v3.second;

        ///////////////////////  
        // 2011 Triggers     //
        ///////////////////////



        ///////////////////////  
        // 2010 Triggers     //
        ///////////////////////

        // Our jet trigger flags
        hlt15u_ = min(2,nHLTObjects("HLT_Jet15U")); 
        hlt30u_ = min(2,nHLTObjects("HLT_Jet30U")); 
        hlt50u_ = min(2,nHLTObjects("HLT_Jet50U")); 
        l16u_   = min(2,nHLTObjects("HLT_L1Jet6U"));
        l110u_  = min(2,nHLTObjects("HLT_L1Jet10U"));

        // If only one jet triggered, see if it is far enough away 
        if (hlt15u_ == 1) {
          LorentzVector p4j = p4HLTObject("HLT_Jet15U",0);
          double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4j);
          if (dr > deltaRCut) hlt15u_ = 2;
        }
        if (hlt30u_ == 1) {
          LorentzVector p4j = p4HLTObject("HLT_Jet30U",0);
          double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4j);
          if (dr > deltaRCut) hlt30u_ = 2;
        }
        if (hlt50u_ == 1) {
          LorentzVector p4j = p4HLTObject("HLT_Jet50U",0);
          double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4j);
          if (dr > deltaRCut) hlt50u_ = 2;
        }
        if (l16u_ == 1) {
          LorentzVector p4j = p4HLTObject("HLT_L1Jet6U",0);
          double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4j);
          if (dr > deltaRCut) l16u_ = 2;
        }
        if (l110u_ == 1) {
          LorentzVector p4j = p4HLTObject("HLT_L1Jet10U",0);
          double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), p4j);
          if (dr > deltaRCut) l110u_ = 2;
        }
        
        // Fill electron & photon triggers
        pair<int, float> pair_el10_lw     = TriggerMatch( els_p4().at(iEl), "HLT_Ele10_LW_L1R");
        pair<int, float> pair_el10_sw     = TriggerMatch( els_p4().at(iEl), "HLT_Ele10_SW_L1R");
        pair<int, float> pair_el10_sw_v2  = TriggerMatch( els_p4().at(iEl), "HLT_Ele10_SW_L1R_v2");

        pair<int, float> pair_el10_lw_id  = TriggerMatch( els_p4().at(iEl), "HLT_Ele10_LW_EleId_L1R");
        pair<int, float> pair_el10_sw_id  = TriggerMatch( els_p4().at(iEl), "HLT_Ele10_SW_EleId_L1R");

        pair<int, float> pair_el15_lw     = TriggerMatch( els_p4().at(iEl), "HLT_Ele15_LW_L1R");
        pair<int, float> pair_el15_sw     = TriggerMatch( els_p4().at(iEl), "HLT_Ele15_SW_L1R");

        pair<int, float> pair_el15_lw_id  = TriggerMatch( els_p4().at(iEl), "HLT_Ele15_LW_EleId_L1R");
        pair<int, float> pair_el15_sw_id  = TriggerMatch( els_p4().at(iEl), "HLT_Ele15_SW_EleId_L1R");

        pair<int, float> pair_el15_sw_cid = TriggerMatch( els_p4().at(iEl), "HLT_Ele15_SW_CaloEleId_L1R");

        pair<int, float> pair_el20_sw     = TriggerMatch( els_p4().at(iEl), "HLT_Ele20_SW_L1R");
        pair<int, float> pair_el25_sw     = TriggerMatch( els_p4().at(iEl), "HLT_Ele25_SW_L1R");

        pair<int, float> pair_Del10_sw    = TriggerMatch( els_p4().at(iEl), "HLT_DoubleEle10_SW_L1R");

        pair<int, float> pair_ph10    = TriggerMatch( els_p4().at(iEl), "HLT_Photon10_L1R");
        pair<int, float> pair_ph10C   = TriggerMatch( els_p4().at(iEl), "HLT_Photon10_Cleaned_L1R");
        pair<int, float> pair_ph15    = TriggerMatch( els_p4().at(iEl), "HLT_Photon15_L1R");
        pair<int, float> pair_ph15C   = TriggerMatch( els_p4().at(iEl), "HLT_Photon15_Cleaned_L1R");
        pair<int, float> pair_ph20C   = TriggerMatch( els_p4().at(iEl), "HLT_Photon20_Cleaned_L1R");

        pair<int, float> pair_el17_sw     = TriggerMatch( els_p4().at(iEl), "HLT_Ele17_SW_L1R");
        pair<int, float> pair_el17_sw_v2  = TriggerMatch( els_p4().at(iEl), "HLT_Ele17_SW_L1R_v2");
        pair<int, float> pair_el17_iso    = TriggerMatch( els_p4().at(iEl), "HLT_Ele17_Isol_L1R");
        pair<int, float> pair_el17_loose  = TriggerMatch( els_p4().at(iEl), "HLT_Ele17_SW_LooseEleId_L1R");
        pair<int, float> pair_el17_sw_cid = TriggerMatch( els_p4().at(iEl), "HLT_Ele17_SW_CaloEleId_L1R");
        pair<int, float> pair_el17_sw_id  = TriggerMatch( els_p4().at(iEl), "HLT_Ele17_SW_EleId_L1R");
        pair<int, float> pair_el17_tiso   = TriggerMatch( els_p4().at(iEl), "HLT_Ele17_SW_TightEleIdIsol_L1R");
        pair<int, float> pair_el17_tiso_v1  = TriggerMatch( els_p4().at(iEl), "HLT_Ele17_SW_TightEleIdIsol_L1R_v1");


        int   ph10    = max( pair_ph10.first, pair_ph10C.first );
        int   ph15    = max( pair_ph15.first, pair_ph15C.first );
        float drph10  = ( pair_ph10.first > pair_ph10C.first ? pair_ph10.second : pair_ph10C.second );
        float drph15  = ( pair_ph15.first > pair_ph15C.first ? pair_ph15.second : pair_ph15C.second );

        // trigger matching
        el17_sw_       = pair_el17_sw.first;
        el17_sw_v2_    = pair_el17_sw_v2.first;
        el17_iso_      = pair_el17_iso.first;
        el17_loose_    = pair_el17_loose.first;
        el17_sw_cid_   = pair_el17_sw_cid.first;
        el17_sw_id_    = pair_el17_sw_id.first;
        el17_tiso_     = pair_el17_tiso.first;
        el17_tiso_v1_  = pair_el17_tiso_v1.first;

        el10_lw_      = pair_el10_lw.first;
        el10_sw_      = pair_el10_sw.first;
        el10_sw_v2_   = pair_el10_sw_v2.first;

        el10_lw_id_   = pair_el10_lw_id.first;
        el10_sw_id_   = pair_el10_sw_id.first;

        el15_lw_      = pair_el15_lw.first;
        el15_sw_      = pair_el15_sw.first;

        el15_lw_id_   = pair_el15_lw_id.first;
        el15_sw_id_   = pair_el15_sw_id.first;

        el15_sw_cid_  = pair_el15_sw_cid.first;

        el20_sw_      = pair_el20_sw.first;
        el25_sw_      = pair_el25_sw.first;
        Del10_sw_     = pair_Del10_sw.first;


        ph10_     = ph10;
        ph15_     = ph15;
        ph20_     = pair_ph20C.first;

        // dr between lepton and closest jet
        drel17_sw_       = pair_el17_sw.second;
        drel17_sw_v2_    = pair_el17_sw_v2.second;
        drel17_iso_      = pair_el17_iso.second;
        drel17_loose_    = pair_el17_loose.second;
        drel17_sw_cid_   = pair_el17_sw_cid.second;
        drel17_sw_id_    = pair_el17_sw_id.second;
        drel17_tiso_     = pair_el17_tiso.second;
        drel17_tiso_v1_  = pair_el17_tiso_v1.second;

        drel10_lw_      = pair_el10_lw.second;
        drel10_sw_      = pair_el10_sw.second;
        drel10_sw_v2_   = pair_el10_sw_v2.second;

        drel10_lw_id_   = pair_el10_lw_id.second;
        drel10_sw_id_   = pair_el10_sw_id.second;

        drel15_lw_      = pair_el15_lw.second;
        drel15_sw_      = pair_el15_sw.second;

        drel15_lw_id_   = pair_el15_lw_id.second;
        drel15_sw_id_   = pair_el15_sw_id.second;

        drel15_sw_cid_  = pair_el15_sw_cid.second;

        drel20_sw_      = pair_el20_sw.second;
        drel25_sw_      = pair_el25_sw.second;
        drDel10_sw_     = pair_Del10_sw.second;

        drph10_     = drph10;
        drph15_     = drph15;
        drph20_     = pair_ph20C.second;

        ///////////////////////  
        // End 2010 Triggers //
        ///////////////////////





        // Find the highest Pt jet separated by at least dRcut from this lepton and fill the jet Pt
        ptj1_       = -999.0;
        ptj1_b2b_   = -999.0;
        dphij1_b2b_ = -999.0;
        nj1_        = 0;
        for (unsigned int iJet = 0; iJet < jets_p4().size(); iJet++) {
          double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), jets_p4().at(iJet) );
          if( dr > deltaRCut && jets_p4().at(iJet).pt() > 10 ) nj1_++;
          if ( dr > deltaRCut && jets_p4().at(iJet).pt() > ptj1_ ){
            ptj1_ = jets_p4().at(iJet).pt();
      
            // back to back in phi
            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iEl), jets_p4().at(iJet) ) );
            if( dphi > deltaPhiCut && jets_p4().at(iJet).pt() > ptj1_b2b_ ){ 
              ptj1_b2b_   = jets_p4().at(iJet).pt();
              dphij1_b2b_ = dphi;
            }
          }
        }

        // Find the highest Pt pfjet separated by at least dRcut from this lepton and fill the pfjet Pt
        ptpfj1_       = -999.0;
        ptpfj1_b2b_   = -999.0;
        dphipfj1_b2b_ = -999.0;
        npfj1_        = 0;
        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
          double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), pfjets_p4().at(iJet) );
          if( dr > deltaRCut && pfjets_p4().at(iJet).pt() > 10 ) npfj1_++;
          if ( dr > deltaRCut && pfjets_p4().at(iJet).pt() > ptpfj1_ ){
            ptpfj1_ = pfjets_p4().at(iJet).pt();
      
            // back to back in phi
            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iEl), pfjets_p4().at(iJet) ) );
            if( dphi > deltaPhiCut && pfjets_p4().at(iJet).pt() > ptpfj1_b2b_ ){ 
              ptpfj1_b2b_   = pfjets_p4().at(iJet).pt();
              dphij1_b2b_ = dphi;
            }
          }
        }

      // Find the highest Pt PF corrected jet separated by at least dRcut from this lepton and fill the jet Pt

        ptpfcj1_       = -999.0; 
        ptpfcj1_b2b_   = -999.0;
        dphipfcj1_b2b_ = -999.0;
        npfcj1_        = 0;
        btagpfc_       = false;

        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
          if ( !passesPFJetID(iJet)) continue;
          LorentzVector jp4 = pfjets_p4()[iJet];
          float jet_cor = jetCorrection(jp4, jet_corrector);
          LorentzVector jp4cor = jp4 * jet_cor;
          if (jp4cor.pt() > 15 && pfjets_simpleSecondaryVertexHighEffBJetTag().at(iJet) > 1.74 ) btagpfc_ = true;
          double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iEl), jp4cor );
          if( dr > deltaRCut && jp4cor.pt() > 10 ) npfcj1_++;
          if ( dr > deltaRCut && jp4cor.pt() > ptpfcj1_ ){
            ptpfcj1_ = jp4cor.pt();

            // back to back in phi
            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iEl), jp4cor ) );
            if( dphi > deltaPhiCut && jp4cor.pt() > ptpfcj1_b2b_ ){
              ptpfcj1_b2b_   = jp4cor.pt();
              dphipfcj1_b2b_ = dphi;
            } 
          }
        }


	// Time to fill the baby for the electrons
        FillBabyNtuple();

      } // closes loop over electrons
      } // closes if statements about whether we want to fill electrons

/* Muons */

      if (eormu == -1 || eormu==13) {
      for ( unsigned int iMu = 0; iMu < mus_p4().size(); iMu++) {
        
        // Apply a pt cut --- moved the cut to 10 GeV (Claudio, 10 Jul 2010)
        if ( mus_p4().at(iMu).pt() < 10.) continue;
        
        // If it is above 20 GeV see if we can make a 
        // Z with another pt>20 FO.  
        bool isaZ = false;
        if (mus_p4().at(iMu).pt() > 20.) {
          for (unsigned int jMu = 0 ; jMu < mus_p4().size(); jMu++) {
            if (iMu == jMu)                             continue;
            if (mus_p4().at(jMu).pt() < 20.)            continue;
            if ( ! muonId( jMu, muonSelectionFO_mu_ttbar) ) continue;
            if ( ! muonId( iMu, muonSelectionFO_mu_ttbar) ) continue;
            LorentzVector w = mus_p4().at(iMu) + mus_p4().at(jMu);
            if (abs(w.mass()-91.) > 20.) continue;
            isaZ = true;
          }
        }
        if (isaZ) continue;
        
        // Initialize baby ntuple
        InitBabyNtuple();
        
        // Load the muon and event quantities
        run_  = evt_run();
        ls_   = evt_lumiBlock();
        evt_  = evt_event();
		weight_ = isData ? 1. : evt_scale1fb();
        pt_   = mus_p4().at(iMu).pt();
        eta_  = mus_p4().at(iMu).eta();
        phi_  = mus_p4().at(iMu).phi();
        id_   = 13*mus_charge().at(iMu);
        tcmet_ = evt_tcmet();
        tcmetphi_ = evt_tcmetPhi();
        pfmet_ = evt_pfmet();
        pfmetphi_ = evt_pfmetPhi();
        iso_ = muonIsoValue(iMu);

        // Pileup - PUSummaryInfoMaker
        pu_nPUvertices_ = puInfo_nPUvertices();
        pu_zpositions_  = puInfo_zpositions();
        pu_sumptlowpt_  = puInfo_sumpt_lowpt();
        pu_sumpthighpt_ = puInfo_sump_highpt();
        pu_instLumi_    = puInfo_instLumi();
        pu_ntrkslowpt_  = puInfo_ntrks_lowpt();
        pu_ntrkshighpt_ = puInfo_ntrks_highpt();

        // Pileup - VertexMaker
        evt_nvtxs_       = evt_nvtxs();
        vtxs_xError_     = vtxs_xError();
        vtxs_yError_     = vtxs_yError();
        vtxs_zError_     = vtxs_zError();
        vtxs_chi2_       = vtxs_chi2();
        vtxs_ndof_       = vtxs_ndof();
        vtxs_sumpt_      = vtxs_sumpt();
        vtxs_isFake_     = vtxs_isFake();
        vtxs_isValid_    = vtxs_isValid();
        vtxs_tracksSize_ = vtxs_tracksSize();
        vtxs_covMatrix_  = vtxs_covMatrix();
        vtxs_position_   = vtxs_position();

        // Pileup - VertexMaker
        evt_ndavtxs_       = evt_ndavtxs();
        davtxs_xError_     = davtxs_xError();
        davtxs_yError_     = davtxs_yError();
        davtxs_zError_     = davtxs_zError();
        davtxs_chi2_       = davtxs_chi2();
        davtxs_ndof_       = davtxs_ndof();
        davtxs_sumpt_      = davtxs_sumpt();
        davtxs_isFake_     = davtxs_isFake();
        davtxs_isValid_    = davtxs_isValid();
        davtxs_tracksSize_ = davtxs_tracksSize();
        davtxs_covMatrix_  = davtxs_covMatrix();
        davtxs_position_   = davtxs_position();



        if (! isData) {
            mcid_       = mus_mc_id().at(iMu);
            mcmotherid_ = mus_mc_motherid().at(iMu);
        }

        //
        num_    = muonId(iMu, NominalTTbarV2);
        numv1_  = muonId(iMu, NominalTTbar);
        numSS_  = muonId(iMu, Nominal);
        numNomSS_  = muonId(iMu, NominalSS);
		numNomSSv2_ = muonId(iMu, NominalSSv2);
        num_wwV1_  = muonId(iMu, NominalWWV1);
        num_OSGv1_ = muonId(iMu, OSGeneric_v1);
        num_OSZv1_ = muonId(iMu, OSZ_v1);
    
        fo_04_  = muonId(iMu, muonSelectionFO_mu_ttbar);
        fo_10_  = muonId(iMu, muonSelectionFO_mu_ttbar_iso10);

        fo_muss04_ = muonId(iMu, muonSelectionFO_mu_ss);
        fo_muss10_ = muonId(iMu, muonSelectionFO_mu_ss_iso10);
        fo_mussV2_04_ = muonId(iMu, muonSelectionFO_mu_ssV2);
        fo_mussV2_10_ = muonId(iMu, muonSelectionFO_mu_ssV2_iso10);

        fo_wwV1_04_  = muonId(iMu, muonSelectionFO_mu_wwV1);
        fo_wwV1_10_  = muonId(iMu, muonSelectionFO_mu_wwV1_iso10);
        fo_wwV1_10_d0_  = muonId(iMu, muonSelectionFO_mu_wwV1_iso10_d0);

        numAug9_ = num_;

        if( !fo_04_ && !fo_10_ &&
            !fo_wwV1_04_ && !fo_wwV1_10_ && !fo_wwV1_10_d0_ &&
            !fo_muss04_ && !fo_muss10_ &&
            !fo_mussV2_04_ && !fo_mussV2_10_
           ) continue;

        // Now REALLY fix it (July 14, 2010)
        if (pt_ > 10.) {
          if (!wasMetCorrectedForThisMuon(iMu, usingTcMet)) {
            float metX = tcmet_ * cos(evt_tcmetPhi());
            float metY = tcmet_ * sin(evt_tcmetPhi());
            fixMetForThisMuon(iMu, metX, metY, usingTcMet);
            tcmet_ = sqrt(metX*metX + metY*metY);
          }
        }
        
        // W transverse mass
        mt_ = Mt( mus_p4().at(iMu), tcmet_, tcmetphi_ );
        pfmt_ = Mt( mus_p4().at(iMu), pfmet_, pfmetphi_ );
        
        // The btag information
        nbjet_ = this_nbjet;
        dRbNear_ =  99.;
        dRbFar_  = -99.;
        for (int ii=0; ii<nbjet_; ii++) {
          unsigned int iJet = bindex[ii];
          float dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), jets_p4().at(iJet));
          if (dr < dRbNear_) dRbNear_ = dr;
          if (dr > dRbFar_)  dRbFar_  = dr;
        }
         // The btag information for pfjets

        nbpfcjet_ = this_nbpfjet;
        dRbpfcNear_ = 99.;
        dRbpfcFar_  = -99.;
        for (int ii=0; ii<nbpfcjet_; ii++) {
          unsigned int iJet = bpfindex[ii];
          LorentzVector jp4 = pfjets_p4()[iJet];
          float jet_cor = jetCorrection(jp4, jet_corrector);
          LorentzVector jp4cor = jp4 * jet_cor;
          float dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), jp4cor);
          if (dr < dRbpfcNear_) dRbpfcNear_ = dr;
          if (dr > dRbpfcFar_)   dRbpfcFar_  = dr;
        }

        
        // Our jet trigger flags
        hlt15u_ = min(2,nHLTObjects("HLT_Jet15U")); 
        hlt30u_ = min(2,nHLTObjects("HLT_Jet30U")); 
        hlt50u_ = min(2,nHLTObjects("HLT_Jet50U")); 
        l16u_   = min(2,nHLTObjects("HLT_L1Jet6U"));
        l110u_  = min(2,nHLTObjects("HLT_L1Jet10U"));
        
        // If only one jet triggered, see if it is far enough away 
        if (hlt15u_ == 1) {
          LorentzVector p4j = p4HLTObject("HLT_Jet15U",0);
          double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4j);
          if (dr > deltaRCut) hlt15u_ = 2;
        }
        if (hlt30u_ == 1) {
          LorentzVector p4j = p4HLTObject("HLT_Jet30U",0);
          double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4j);
          if (dr > deltaRCut) hlt30u_ = 2;
        }
        if (hlt50u_ == 1) {
          LorentzVector p4j = p4HLTObject("HLT_Jet50U",0);
          double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4j);
          if (dr > deltaRCut) hlt50u_ = 2;
        }
        if (l16u_ == 1) {
          LorentzVector p4j = p4HLTObject("HLT_L1Jet6U",0);
          double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4j);
          if (dr > deltaRCut) l16u_ = 2;
        }
        if (l110u_ == 1) {
          LorentzVector p4j = p4HLTObject("HLT_L1Jet10U",0);
          double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), p4j);
          if (dr > deltaRCut) l110u_ = 2;
        }
        
        // Now fill the muon trigger flags
        pair<int, float> pair_mu5 = TriggerMatch( mus_p4().at(iMu), "HLT_Mu5");
        pair<int, float> pair_mu7 = TriggerMatch( mus_p4().at(iMu), "HLT_Mu7");
        pair<int, float> pair_mu9 = TriggerMatch( mus_p4().at(iMu), "HLT_Mu9");
        pair<int, float> pair_mu11 = TriggerMatch( mus_p4().at(iMu), "HLT_Mu11");
        pair<int, float> pair_mu13 = TriggerMatch( mus_p4().at(iMu), "HLT_Mu13_v1");
        pair<int, float> pair_mu15 = TriggerMatch( mus_p4().at(iMu), "HLT_Mu15_v1");
        pair<int, float> pair_mu17 = TriggerMatch( mus_p4().at(iMu), "HLT_Mu17_v1");
        mu5_     = pair_mu5.first;
        drmu5_   = pair_mu5.second;
        mu7_     = pair_mu7.first;
        drmu7_   = pair_mu7.second;
        mu9_     = pair_mu9.first;
        drmu9_   = pair_mu9.second;
        mu11_    = pair_mu11.first;
        drmu11_  = pair_mu11.second;
        mu13_    = pair_mu13.first;
        drmu13_  = pair_mu13.second;
        mu15_    = pair_mu15.first;
        drmu15_  = pair_mu15.second;
        mu17_    = pair_mu17.first;
        drmu17_  = pair_mu17.second;

        // Find the highest Pt jet separated by at least dRcut from this lepton and fill the jet Pt
        ptj1_       = -999.0;
        ptj1_b2b_   = -999.0;
        dphij1_b2b_ = -999.0;
        nj1_        = 0;
        for (unsigned int iJet = 0; iJet < jets_p4().size(); iJet++) {
          double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), jets_p4().at(iJet) );
          if( dr > deltaRCut && jets_p4().at(iJet).pt() > 10 ) nj1_++;
          if ( dr > deltaRCut && jets_p4().at(iJet).pt() > ptj1_ ){
            ptj1_ = jets_p4().at(iJet).pt();
    
            // back to back in phi
            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iMu), jets_p4().at(iJet) ) );
            if( dphi > deltaPhiCut && jets_p4().at(iJet).pt() > ptj1_b2b_ ){        
              ptj1_b2b_   = jets_p4().at(iJet).pt();
              dphij1_b2b_ = dphi;
            }
          }
        }

        // Find the highest Pt pfjet separated by at least dRcut from this lepton and fill the pfjet Pt
	ptpfj1_       = -999.0;
        ptpfj1_b2b_   = -999.0;
        dphipfj1_b2b_ = -999.0;
        npfj1_        = 0;
        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
          double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), pfjets_p4().at(iJet) );
          if( dr > deltaRCut && pfjets_p4().at(iJet).pt() > 10 ) npfj1_++;
          if ( dr > deltaRCut && pfjets_p4().at(iJet).pt() > ptpfj1_ ){
            ptpfj1_ = pfjets_p4().at(iJet).pt();
    
            // back to back in phi
            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iMu), pfjets_p4().at(iJet) ) );
            if( dphi > deltaPhiCut && pfjets_p4().at(iJet).pt() > ptpfj1_b2b_ ){        
              ptpfj1_b2b_   = pfjets_p4().at(iJet).pt();
              dphipfj1_b2b_ = dphi;
            }
          }
        }

       // Find the highest Pt PF corrected jet separated by at least dRcut from this lepton and fill the jet Pt

        ptpfcj1_       = -999.0;
        ptpfcj1_b2b_   = -999.0;
        dphipfcj1_b2b_ = -999.0;
        npfcj1_        = 0;
        btagpfc_       = false;

        for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
       // JetID
          if ( !passesPFJetID(iJet)) continue;
          LorentzVector jp4 = pfjets_p4()[iJet];
          float jet_cor = jetCorrection(jp4, jet_corrector);
          LorentzVector jp4cor = jp4 * jet_cor;
          if (jp4cor.pt() > 15 && pfjets_simpleSecondaryVertexHighEffBJetTag().at(iJet) > 1.74 ) btagpfc_ = true;
          double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iMu), jp4cor );
          if( dr > deltaRCut && jp4cor.pt() > 10 ) npfcj1_++;
          if ( dr > deltaRCut && jp4cor.pt() > ptpfcj1_ ){
            ptpfcj1_ = jp4cor.pt();

            // back to back in phi
            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iMu), jp4cor ) );
            if( dphi > deltaPhiCut && jp4cor.pt() > ptpfcj1_b2b_ ){
              ptpfcj1_b2b_   = jp4cor.pt();
              dphipfcj1_b2b_ = dphi;
            }
          }
        }


        // Time to fill the baby for the muons
        FillBabyNtuple();
      
      }// closes loop over muons
      } // closes if statements about whether we want to fill muons

    }// closes loop over events
  }  // closes loop over files

  cout << "   " <<endl;
  CloseBabyNtuple();
  return;

} // closes myLooper function  

//------------------------------------------
// Initialize baby ntuple variables
//------------------------------------------
void myBabyMaker::InitBabyNtuple () {
  run_ = -1;
  ls_  = -1;
  evt_ = -1;
  weight_ = 1.;

  // Pileup - PUSummaryInfoMaker
  pu_nPUvertices_ = -1;
  pu_sumptlowpt_.clear();
  pu_sumpthighpt_.clear();
  pu_instLumi_.clear();
  pu_ntrkslowpt_.clear();
  pu_ntrkshighpt_.clear();

  // Pileup - VertexMaker
  evt_nvtxs_ = -1;
  vtxs_xError_.clear();
  vtxs_yError_.clear();
  vtxs_zError_.clear();
  vtxs_chi2_.clear();
  vtxs_ndof_.clear();
  vtxs_sumpt_.clear();
  vtxs_isFake_.clear();
  vtxs_isValid_.clear();
  vtxs_tracksSize_.clear();
  vtxs_covMatrix_.clear();
  vtxs_position_.clear();

  // Pileup - VertexMaker
  evt_ndavtxs_ = -1;
  davtxs_xError_.clear();
  davtxs_yError_.clear();
  davtxs_zError_.clear();
  davtxs_chi2_.clear();
  davtxs_ndof_.clear();
  davtxs_sumpt_.clear();
  davtxs_isFake_.clear();
  davtxs_isValid_.clear();
  davtxs_tracksSize_.clear();
  davtxs_covMatrix_.clear();
  davtxs_position_.clear();



  id_  = -1;
  pt_  = -999.;
  eta_ = -999.;
  phi_ = -999.;
  scet_ = -999.;
  tcmet_ = -999.;
  tcmetphi_ = -999.;
  pfmet_ = -999.;
  pfmetphi_ = -999.;
  iso_ = -999.;
  hlt15u_ = 0;
  hlt30u_ = 0;
  hlt50u_ = 0;
  l16u_   = 0;
  l110u_  = 0;
  fo_04_ = false;
  fo_10_ = false;

  fo_muss04_ = false;
  fo_muss10_ = false;
  fo_mussV2_04_ = false;
  fo_mussV2_10_ = false;

  fo_wwV1_04_ = false;
  fo_wwV1_10_ = false;
  fo_wwV1_10_d0_ = false;

  v1_  = false;
  v2_  = false;
  v3_  = false;
  num_ = false;

  v1SS_  = false;
  v2SS_  = false;
  v3SS_  = false;
  numSS_ = false;
  numNomSS_ = false;
  numNomSSv2_ = false;

  v1SSAug9_  = false;
  v2SSAug9_  = false;
  v3SSAug9_  = false;
  numSSAug9_ = false;

  numv1_ = false;
  numAug9_ = false;
  v1Aug9_  = false;
  v2Aug9_  = false;
  v3Aug9_  = false;

  numOct6_ = false;
  v1Oct6_  = false;
  v2Oct6_  = false;
  v3Oct6_  = false;

  numOSOct18_ = false;
  v1OSOct18_  = false;
  v2OSOct18_  = false;
  v3OSOct18_  = false;

  numSSOct18_ = false;
  v1SSOct18_  = false;
  v2SSOct18_  = false;
  v3SSOct18_  = false;

  numSSV2_ = false;
  v1SSV2_  = false;
  v2SSV2_  = false;
  v3SSV2_  = false;

  v1_wwV1_  = false;
  v2_wwV1_  = false;
  v3_wwV1_  = false;
  v4_wwV1_  = false;
  num_wwV1_ = false;

  num_OSGv1_ = false;
  num_OSZv1_ = false;

  //
  ph10_ = 0;
  ph15_ = 0;
  ph20_ = 0;
  el10_lw_ = 0;
  el10_sw_ = 0;
  el10_sw_v2_ = 0;
  el10_lw_id_ = 0;
  el10_sw_id_ = 0;
  el15_lw_ = 0;
  el15_sw_ = 0;
  el15_lw_id_ = 0;
  el15_sw_id_ = 0;
  el15_sw_cid_ = 0;
  el20_sw_ = 0;
  el25_sw_ = 0;
  Del10_sw_ = 0;

  el17_sw_ = 0;
  el17_sw_v2_ = 0;
  el17_iso_ =0;
  el17_loose_ =0;
  el17_sw_cid_ =0;
  el17_sw_id_ =0;
  el17_tiso_ =0;
  el17_tiso_v1_ =0;

  //
  drph10_ = 99.;
  drph15_ = 99.;
  drph20_ = 99.;
  drel10_lw_ = 99.;
  drel10_sw_ = 99.;
  drel10_sw_v2_ = 99.;
  drel10_lw_id_ = 99.;
  drel10_sw_id_ = 99.;
  drel15_lw_ = 99.;
  drel15_sw_ = 99.;
  drel15_lw_id_ = 99.;
  drel15_sw_id_ = 99.;
  drel15_sw_cid_ = 99.;
  drel20_sw_ = 99.;
  drel25_sw_ = 99.;
  drDel10_sw_ = 99.;
  drel17_sw_ = 99;
  drel17_sw_v2_ = 99;
  drel17_iso_ =99;
  drel17_loose_ =99;
  drel17_sw_cid_ =99;
  drel17_sw_id_ =99;
  drel17_tiso_ =99;
  drel17_tiso_v1_ =99;


  //
  mu5_  = 0;
  mu7_  = 0;
  mu9_  = 0;
  mu11_  = 0;
  mu13_  = 0;
  mu15_  = 0;
  mu17_  = 0;

  //
  drmu5_  = 99.;
  drmu7_  = 99.;
  drmu9_  = 99.;
  drmu11_  = 99.;
  drmu13_  = 99.;
  drmu15_  = 99.;
  drmu17_  = 99.;

  //
  nbjet_  = 0;
  dRbNear_ = 99.;
  dRbFar_ = -99.;
  nbpfcjet_  = 0;
  dRbpfcNear_ = 99.;
  dRbpfcFar_ = -99.;


  ptj1_   = 0.;
  nj1_    = 0;
  ptj1_b2b_ = -999.;
  dphij1_b2b_ = -999.;
  ptpfj1_   = 0.;
  npfj1_    = 0;
  ptpfj1_b2b_ = -999.;
  dphipfj1_b2b_ = -999.;

  ptpfcj1_   = 0.;
  npfcj1_    = 0;
  ptpfcj1_b2b_ = -999.;
  dphipfcj1_b2b_ = -999.;
  btagpfc_ = false;

  mt_ = -999;
  pfmt_ = -999;

  //
  q3_ = false;
  els_exp_innerlayers_ = 999;
  els_exp_innerlayers39X_ = 999;


  mcid_ = 0;
  mcmotherid_ = 0;

}
//-------------------------------------
// Book the baby ntuple
//-------------------------------------
void myBabyMaker::MakeBabyNtuple(const char *babyFilename)
{
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();
    babyFile_ = new TFile(Form("%s", babyFilename), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");

    babyTree_->Branch("run",          &run_,         "run/I"         );
    babyTree_->Branch("ls",           &ls_,          "ls/I"          );
    babyTree_->Branch("evt",          &evt_,         "evt/I"         );
	  babyTree_->Branch("weight",       &weight_,      "weight/F"      );

    // Pileup - PUSummaryInfoMaker
    babyTree_->Branch("pu_nPUvertices", &pu_nPUvertices_ );
    babyTree_->Branch("pu_zpositions" , &pu_zpositions_  );
    babyTree_->Branch("pu_sumptlowpt" , &pu_sumptlowpt_  );
    babyTree_->Branch("pu_sumpthighpt", &pu_sumpthighpt_ );
    babyTree_->Branch("pu_instLumi"   , &pu_instLumi_    );
    babyTree_->Branch("pu_ntrkslowpt" , &pu_ntrkslowpt_  );
    babyTree_->Branch("pu_ntrkshighpt", &pu_ntrkshighpt_ );

    // Pileup - VertexMaker
    babyTree_->Branch("evt_nvtxs"      , &evt_nvtxs_      );
    babyTree_->Branch("vtxs_xError"    , &vtxs_xError_    );
    babyTree_->Branch("vtxs_yError"    , &vtxs_yError_    );
    babyTree_->Branch("vtxs_zError"    , &vtxs_zError_    );
    babyTree_->Branch("vtxs_chi2"      , &vtxs_chi2_      );
    babyTree_->Branch("vtxs_ndof"      , &vtxs_ndof_      );
    babyTree_->Branch("vtxs_sumpt"     , &vtxs_sumpt_     );
    babyTree_->Branch("vtxs_isFake"    , &vtxs_isFake_    );
    babyTree_->Branch("vtxs_isValid"   , &vtxs_isValid_   );
    babyTree_->Branch("vtxs_tracksSize", &vtxs_tracksSize_);
    babyTree_->Branch("vtxs_covMatrix" , &vtxs_covMatrix_ );
    babyTree_->Branch("vtxs_position"  , &vtxs_position_  );

    // Pileup - VertexMaker
    babyTree_->Branch("evt_ndavtxs"      , &evt_ndavtxs_      );
    babyTree_->Branch("davtxs_xError"    , &davtxs_xError_    );
    babyTree_->Branch("davtxs_yError"    , &davtxs_yError_    );
    babyTree_->Branch("davtxs_zError"    , &davtxs_zError_    );
    babyTree_->Branch("davtxs_chi2"      , &davtxs_chi2_      );
    babyTree_->Branch("davtxs_ndof"      , &davtxs_ndof_      );
    babyTree_->Branch("davtxs_sumpt"     , &davtxs_sumpt_     );
    babyTree_->Branch("davtxs_isFake"    , &davtxs_isFake_    );
    babyTree_->Branch("davtxs_isValid"   , &davtxs_isValid_   );
    babyTree_->Branch("davtxs_tracksSize", &davtxs_tracksSize_);
    babyTree_->Branch("davtxs_covMatrix" , &davtxs_covMatrix_ );
    babyTree_->Branch("davtxs_position"  , &davtxs_position_  );



    babyTree_->Branch("pt",           &pt_,          "pt/F"         );
    babyTree_->Branch("eta",          &eta_,         "eta/F"         );
    babyTree_->Branch("phi",          &phi_,         "phi/F"         );
    babyTree_->Branch("scet",          &scet_,         "scet/F"         );
    babyTree_->Branch("tcmet",          &tcmet_,         "tcmet/F"         );
    babyTree_->Branch("tcmetphi",          &tcmetphi_,         "tcmetphi/F"         );
    babyTree_->Branch("pfmet",          &pfmet_,         "pfmet/F"         );
    babyTree_->Branch("pfmetphi",          &pfmetphi_,         "pfmetphi/F"         );
    babyTree_->Branch("iso",          &iso_,         "iso/F"         );
    babyTree_->Branch("id",          &id_,         "id/I"         );

    babyTree_->Branch("hlt15u",       &hlt15u_,       "hlt15u/I"      );
    babyTree_->Branch("hlt30u",       &hlt30u_,       "hlt30u/I"      );
    babyTree_->Branch("hlt50u",       &hlt50u_,       "hlt50u/I"      );
    babyTree_->Branch("l16u",         &l16u_,         "l16uu/I"      );
    babyTree_->Branch("l110",         &l110u_,        "l110u/I"      );

    babyTree_->Branch("fo_04",         &fo_04_,        "fo_04/O"      );
    babyTree_->Branch("fo_10",         &fo_10_,        "fo_10/O"      );
    babyTree_->Branch("fo_muss04",         &fo_muss04_,        "fo_muss04/O"      );
    babyTree_->Branch("fo_muss10",         &fo_muss10_,        "fo_muss10/O"      );
    babyTree_->Branch("fo_mussV2_04",         &fo_mussV2_04_,        "fo_mussV2_04/O"      );
    babyTree_->Branch("fo_mussV2_10",         &fo_mussV2_10_,        "fo_mussV2_10/O"      );

    babyTree_->Branch("fo_wwV1_04",         &fo_wwV1_04_,        "fo_wwV1_04/O"      );
    babyTree_->Branch("fo_wwV1_10",         &fo_wwV1_10_,        "fo_wwV1_10/O"      );
    babyTree_->Branch("fo_wwV1_10_d0",         &fo_wwV1_10_d0_,        "fo_wwV1_10_d0/O"      );
    babyTree_->Branch("v1",         &v1_,        "v1/O"      );
    babyTree_->Branch("v2",         &v2_,        "v2/O"      );
    babyTree_->Branch("v3",         &v3_,        "v3/O"      );
    babyTree_->Branch("num",         &num_,        "num/O"      );
    babyTree_->Branch("numv1",         &numv1_,        "numv1/O"      );

    babyTree_->Branch("v1SS",         &v1SS_,        "v1SS/O"      );
    babyTree_->Branch("v2SS",         &v2SS_,        "v2SS/O"      );
    babyTree_->Branch("v3SS",         &v3SS_,        "v3SS/O"      );
    babyTree_->Branch("numSS",         &numSS_,        "numSS/O"      );
    babyTree_->Branch("numNomSS",         &numNomSS_,        "numNomSS/O"      );
    babyTree_->Branch("numNomSSv2",         &numNomSSv2_,        "numNomSSv2/O"      );

    babyTree_->Branch("v1SSAug9",         &v1SSAug9_,        "v1SSAug9/O"      );
    babyTree_->Branch("v2SSAug9",         &v2SSAug9_,        "v2SSAug9/O"      );
    babyTree_->Branch("v3SSAug9",         &v3SSAug9_,        "v3SSAug9/O"      );
    babyTree_->Branch("numSSAug9",         &numSSAug9_,        "numSSAug9/O"      );

    babyTree_->Branch("v1Aug9",         &v1Aug9_,        "v1Aug9/O"      );
    babyTree_->Branch("v2Aug9",         &v2Aug9_,        "v2Aug9/O"      );
    babyTree_->Branch("v3Aug9",         &v3Aug9_,        "v3Aug9/O"      );
    babyTree_->Branch("numAug9",         &numAug9_,        "numAug9/O"      );

    babyTree_->Branch("v1Oct6",         &v1Oct6_,        "v1Oct6/O"      );
    babyTree_->Branch("v2Oct6",         &v2Oct6_,        "v2Oct6/O"      );
    babyTree_->Branch("v3Oct6",         &v3Oct6_,        "v3Oct6/O"      );
    babyTree_->Branch("numOct6",         &numOct6_,        "numOct6/O"      );

    babyTree_->Branch("v1SSOct18",         &v1SSOct18_,        "v1SSOct18/O"      );
    babyTree_->Branch("v2SSOct18",         &v2SSOct18_,        "v2SSOct18/O"      );
    babyTree_->Branch("v3SSOct18",         &v3SSOct18_,        "v3SSOct18/O"      );
    babyTree_->Branch("numSSOct18",         &numSSOct18_,        "numSSOct18/O"      );

    babyTree_->Branch("v1SSV2",         &v1SSV2_,        "v1SSV2/O"      );
    babyTree_->Branch("v2SSV2",         &v2SSV2_,        "v2SSV2/O"      );
    babyTree_->Branch("v3SSV2",         &v3SSV2_,        "v3SSV2/O"      );
    babyTree_->Branch("numSSV2",         &numSSV2_,        "numSSV2/O"      );


    babyTree_->Branch("v1OSOct18",         &v1OSOct18_,        "v1OSOct18/O"      );
    babyTree_->Branch("v2OSOct18",         &v2OSOct18_,        "v2OSOct18/O"      );
    babyTree_->Branch("v3OSOct18",         &v3OSOct18_,        "v3OSOct18/O"      );
    babyTree_->Branch("numOSOct18",         &numOSOct18_,        "numOSOct18/O"      );

    babyTree_->Branch("v1_wwV1",         &v1_wwV1_,        "v1_wwV1/O"      );
    babyTree_->Branch("v2_wwV1",         &v2_wwV1_,        "v2_wwV1/O"      );
    babyTree_->Branch("v3_wwV1",         &v3_wwV1_,        "v3_wwV1/O"      );
    babyTree_->Branch("v4_wwV1",         &v4_wwV1_,        "v4_wwV1/O"      );
    babyTree_->Branch("num_wwV1",         &num_wwV1_,        "num_wwV1/O"      );

    babyTree_->Branch("num_OSGv1",         &num_OSGv1_,        "num_OSGv1/O"      );
    babyTree_->Branch("num_OSZv1",         &num_OSZv1_,        "num_OSZv1/O"      );



    //
    babyTree_->Branch("ph10",       &ph10_,       "ph10/I"      );
    babyTree_->Branch("ph15",       &ph15_,       "ph15/I"      );
    babyTree_->Branch("ph20",       &ph20_,       "ph20/I"      );
    babyTree_->Branch("el10_lw",         &el10_lw_,         "el10_lw/I"      );
    babyTree_->Branch("el10_sw",         &el10_sw_,         "el10_sw/I"      );
    babyTree_->Branch("el10_sw_v2",         &el10_sw_v2_,         "el10_sw_v2/I"      );
    babyTree_->Branch("el10_lw_id",         &el10_lw_id_,         "el10_lw_id/I"      );
    babyTree_->Branch("el10_sw_id",         &el10_sw_id_,         "el10_sw_id/I"      );
    babyTree_->Branch("el15_lw",         &el15_lw_,         "el15_lw/I"      );
    babyTree_->Branch("el15_sw",         &el15_sw_,         "el15_sw/I"      );
    babyTree_->Branch("el15_lw_id",         &el15_lw_id_,         "el15_lw_id/I"      );
    babyTree_->Branch("el15_sw_id",         &el15_sw_id_,         "el15_sw_id/I"      );
    babyTree_->Branch("el15_sw_cid",         &el15_sw_cid_,         "el15_sw_cid/I"      );
    babyTree_->Branch("el20_sw",         &el20_sw_,         "el20_sw/I"      );
    babyTree_->Branch("el25_sw",         &el25_sw_,         "el25_sw/I"      );
    babyTree_->Branch("Del10_sw",         &Del10_sw_,         "Del10_sw/I"      );

    babyTree_->Branch("el17_sw",         &el17_sw_,         "el17_sw/I"      );
    babyTree_->Branch("el17_sw_v2",         &el17_sw_v2_,         "el17_sw_v2/I"      );
    babyTree_->Branch("el17_iso",         &el17_iso_,         "el17_iso/I"      );
    babyTree_->Branch("el17_loose",         &el17_loose_,         "el17_loose/I"      );
    babyTree_->Branch("el17_sw_cid",         &el17_sw_cid_,         "el17_sw_cid/I"      );
    babyTree_->Branch("el17_sw_id",         &el17_sw_id_,         "el17_sw_id/I"      );
    babyTree_->Branch("el17_tiso",         &el17_tiso_,         "el17_tiso/I"      );
    babyTree_->Branch("el17_tiso_v1",         &el17_tiso_v1_,         "el17_tiso_v1/I"      );


    //
    babyTree_->Branch("drph10",       &drph10_,       "drph10/F"      );
    babyTree_->Branch("drph15",       &drph15_,       "drph15/F"      );
    babyTree_->Branch("drph20",       &drph20_,       "drph20/F"      );
    babyTree_->Branch("drel10_lw",         &drel10_lw_,         "drel10_lw/F"      );
    babyTree_->Branch("drel10_sw",         &drel10_sw_,         "drel10_sw/F"      );
    babyTree_->Branch("drel10_sw_v2",         &drel10_sw_v2_,         "drel10_sw_v2/F"      );
    babyTree_->Branch("drel10_lw_id",         &drel10_lw_id_,         "drel10_lw_id/F"      );
    babyTree_->Branch("drel10_sw_id",         &drel10_sw_id_,         "drel10_sw_id/F"      );
    babyTree_->Branch("drel15_lw",         &drel15_lw_,         "drel15_lw/F"      );
    babyTree_->Branch("drel15_sw",         &drel15_sw_,         "drel15_sw/F"      );
    babyTree_->Branch("drel15_lw_id",         &drel15_lw_id_,         "drel15_lw_id/F"      );
    babyTree_->Branch("drel15_sw_id",         &drel15_sw_id_,         "drel15_sw_id/F"      );
    babyTree_->Branch("drel15_sw_cid",         &drel15_sw_cid_,         "drel15_sw_cid/F"      );
    babyTree_->Branch("drel20_sw",         &drel20_sw_,         "drel20_sw/F"      );
    babyTree_->Branch("drel25_sw",         &drel25_sw_,         "drel25_sw/F"      );
    babyTree_->Branch("drDel10_sw",         &drDel10_sw_,         "drDel10_sw/F"      );

    babyTree_->Branch("drel17_sw",         &drel17_sw_,         "drel17_sw/F"      );
    babyTree_->Branch("drel17_sw_v2",         &drel17_sw_v2_,         "drel17_sw_v2/F"      );
    babyTree_->Branch("drel17_iso",         &drel17_iso_,         "drel17_iso/F"      );
    babyTree_->Branch("drel17_loose",         &drel17_loose_,         "drel17_loose/F"      );
    babyTree_->Branch("drel17_sw_cid",         &drel17_sw_cid_,         "drel17_sw_cid/F"      );
    babyTree_->Branch("drel17_sw_id",         &drel17_sw_id_,         "drel17_sw_id/F"      );
    babyTree_->Branch("drel17_tiso",         &drel17_tiso_,         "drel17_tiso/F"      );
    babyTree_->Branch("drel17_tiso_v1",         &drel17_tiso_v1_,         "drel17_tiso_v1/F"      );

    //
    babyTree_->Branch("mu17",       &mu17_,       "mu17/I"      );
    babyTree_->Branch("mu15",       &mu15_,       "mu15/I"      );
    babyTree_->Branch("mu13",       &mu13_,       "mu13/I"      );
    babyTree_->Branch("mu11",       &mu11_,       "mu11/I"      );
    babyTree_->Branch("mu9",       &mu9_,       "mu9/I"      );
    babyTree_->Branch("mu7",       &mu7_,       "mu7/I"      );
    babyTree_->Branch("mu5",       &mu5_,       "mu5/I"      );

    //
    babyTree_->Branch("drmu17",       &drmu17_,       "drmu17/F"      );
    babyTree_->Branch("drmu15",       &drmu15_,       "drmu15/F"      );
    babyTree_->Branch("drmu13",       &drmu13_,       "drmu13/F"      );
    babyTree_->Branch("drmu11",       &drmu11_,       "drmu11/F"      );
    babyTree_->Branch("drmu9",       &drmu9_,       "drmu9/F"      );
    babyTree_->Branch("drmu7",       &drmu7_,       "drmu7/F"      );
    babyTree_->Branch("drmu5",       &drmu5_,       "drmu5/F"      );

    babyTree_->Branch("nbjet",       &nbjet_,       "nbjet/I"      );
    babyTree_->Branch("dRNear",       &dRbNear_,       "dRbNear/F"      );
    babyTree_->Branch("dRFar",       &dRbFar_,       "dRbFar/F"      );

    babyTree_->Branch("nbpfcjet",       &nbpfcjet_,       "nbpfcjet/I"      );
    babyTree_->Branch("dRpfcNear",       &dRbpfcNear_,       "dRbpfcNear/F"      );   
    babyTree_->Branch("dRpfcFar",       &dRbpfcFar_,       "dRbpfcFar/F"      );


    babyTree_->Branch("ptj1",       &ptj1_,       "ptj1/F"      );
    babyTree_->Branch("nj1",       &nj1_,       "nj1/I"      );
    babyTree_->Branch("ptj1_b2b",       &ptj1_b2b_,       "ptj1_b2b/F"      );
    babyTree_->Branch("dphij1_b2b",       &dphij1_b2b_,       "dphij1_b2b/F"      );
    babyTree_->Branch("ptpfj1",       &ptpfj1_,       "ptpfj1/F"      );
    babyTree_->Branch("npfj1",       &npfj1_,       "npfj1/I"      );
    babyTree_->Branch("ptpfj1_b2b",       &ptpfj1_b2b_,       "ptpfj1_b2b/F"      );
    babyTree_->Branch("dphipfj1_b2b",       &dphipfj1_b2b_,       "dphipfj1_b2b/F"      );

    babyTree_->Branch("ptpfcj1",       &ptpfcj1_,       "ptpfcj1/F"      );
    babyTree_->Branch("npfcj1",       &npfcj1_,       "npfcj1/I"      );
    babyTree_->Branch("ptpfcj1_b2b",       &ptpfcj1_b2b_,       "ptpfcj1_b2b/F"      );
    babyTree_->Branch("dphipfcj1_b2b",       &dphipfcj1_b2b_,       "dphipfcj1_b2b/F"      );
    babyTree_->Branch("btagpfc",       &btagpfc_,       "btagpfc/O"      );

    babyTree_->Branch("mt",          &mt_,         "mt/F"         );
    babyTree_->Branch("pfmt",          &pfmt_,         "pfmt/F"         );

    babyTree_->Branch("q3",          &q3_,         "q3/O"         );

    babyTree_->Branch("els_exp_innerlayers", &els_exp_innerlayers_, "els_exp_innerlayers/I" );
    babyTree_->Branch("els_exp_innerlayers39X", &els_exp_innerlayers39X_, "els_exp_innerlayers39X/I" );
    babyTree_->Branch("mcid",       &mcid_,       "mcid/I"      );
    babyTree_->Branch("mcmotherid", &mcmotherid_, "mcmotherid/I"      );


}
//----------------------------------
// Fill the baby
//----------------------------------
void myBabyMaker::FillBabyNtuple()
{
    babyTree_->Fill();
}
//--------------------------------
// Close the baby
//--------------------------------
void myBabyMaker::CloseBabyNtuple()
{
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();
}




