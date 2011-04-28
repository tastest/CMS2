// C++ includes
#include <iostream>
#include <set>
#include <string>

// ROOT includes
#include "TSystem.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TChainElement.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Math/VectorUtil.h"
#include "TObjArray.h"

#if 1
#define private public
#include "TString.h"
#endif

#include "TRegexp.h"

// TAS includes
#include "myBabyMaker.h"
#include "../CORE/electronSelections.cc"
#include "../CORE/electronSelectionsParameters.cc"
#include "../CORE/eventSelections.cc"
#include "../CORE/jetSelections.cc"
#include "../CORE/metSelections.cc"
#include "../CORE/MITConversionUtilities.cc"
#include "../CORE/muonSelections.cc"
#include "../CORE/trackSelections.cc"
#include "../CORE/triggerUtils.cc"
#include "../Tools/goodrun.cc"

// namespaces
using namespace std;
using namespace tas;
bool header1 = false;
bool header2 = false;

void PrintTriggerDebugHeader(string outfileName){

    int width  = 7;
    ofstream outfile( Form("triggerStudy/%s", outfileName.c_str() ), ios::app );
    outfile  
        <<  setw(width) << "itrg"
        <<  setw(width) << "id"
        <<  setw(width) << "match"
        <<  setw(width) << "matchId"
        <<  setw(width) << "dr"
        <<  setw(width) << "lep pt"
        <<  setw(width) << "trg pt"
        <<  setw(width) << "lep eta"
        <<  setw(width) << "trg eta"
        <<  setw(width) << "lep phi"
        <<  setw(width) << "trg phi" 
        << "\t"         << "trigString" << endl << endl;
    outfile.close();
}
void PrintTriggerDebugLine(int itrg, int id, bool match, bool matchId, double dr, LorentzVector lepton_p4, LorentzVector p4tr, string trigString, int nTrig, string outfileName ){

    ofstream outfile( Form("triggerStudy/%s", outfileName.c_str() ), ios::app );

    int precis = 2;
    int width  = 7;
    outfile.setf( ios::fixed, ios::floatfield );
    outfile << setprecision(precis) << setw(width) << setfill(' ') << itrg
            << setprecision(precis) << setw(width) << setfill(' ') << id
            << setprecision(precis) << setw(width) << setfill(' ') << match
            << setprecision(precis) << setw(width) << setfill(' ') << matchId
            << setprecision(precis) << setw(width) << setfill(' ') << dr
            << setprecision(precis) << setw(width) << setfill(' ') << lepton_p4.pt()
            << setprecision(precis) << setw(width) << setfill(' ') << p4tr.pt()
            << setprecision(precis) << setw(width) << setfill(' ') << lepton_p4.eta()
            << setprecision(precis) << setw(width) << setfill(' ') << p4tr.eta()
            << setprecision(precis) << setw(width) << setfill(' ') << lepton_p4.phi()
            << setprecision(precis) << setw(width) << setfill(' ') << p4tr.phi()
            << "\t" << trigString << endl;
    if( itrg == nTrig-1 ) outfile << endl;
    outfile.close();
}

bool found_ele8 = false;
bool found_ele8_CaloIdL_TrkIdVL = false;
bool found_ele8_CaloIdL_CaloIsoVL = false;
bool found_ele17_CaloIdL_CaloIsoVL = false;
bool found_ele8_CaloIdL_CaloIsoVL_Jet40 = false;
// function for dR matching offline letpon to trigger object 
pair<int, float> TriggerMatch( LorentzVector lepton_p4, const char* trigString, double dR_cut = 0.4, int pid = 11 ){
    float dR_min = numeric_limits<float>::max();
    dR_min = 99.0;
    int nTrig = nHLTObjects( trigString );
    if (nTrig > 0) {
        bool match   = false;
        bool matchId = false;


        for (int itrg=0; itrg<nTrig; itrg++) {
            LorentzVector p4tr = p4HLTObject( trigString, itrg );
            int id             = idHLTObject( trigString, itrg );
            double dr = ROOT::Math::VectorUtil::DeltaR( lepton_p4, p4tr);
            if ( dr < dR_cut ){
                match = true;
                if( abs(id) == abs(pid) ) matchId = true;
            }
            if (dr < dR_min) dR_min = dr;


            //////////////////////////
            // Debug Mixed Triggers //
            //////////////////////////
            vector<string> triggers;
            triggers.push_back("HLT_Ele8");
            triggers.push_back("HLT_Ele8_CaloIdL_TrkIdVL");
            triggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL");
            triggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL");
            triggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40");
            for (unsigned int i=0; i < triggers.size(); i++) {
                if (
                    ( strcmp( trigString, Form( "%s_%s", triggers.at(i).c_str(), "v1" ) ) == 0 ) || 
                    ( strcmp( trigString, Form( "%s_%s", triggers.at(i).c_str(), "v2" ) ) == 0 )
                    ) {

                    //
                    if(
                        ( strcmp( trigString, "HLT_Ele8_v1") == 0 ) || 
                        ( strcmp( trigString, "HLT_Ele8_v2") == 0 )
                        ){
                        if(found_ele8 == false){
                            found_ele8 = true;
                            PrintTriggerDebugHeader( Form("%s.txt",triggers.at(i).c_str()) );
                        }
                    }

                    //
                    if(
                        ( strcmp( trigString, "HLT_Ele8_CaloIdL_TrkIdVL_v1") == 0 ) || 
                        ( strcmp( trigString, "HLT_Ele8_CaloIdL_TrkIdVL_v2") == 0 )
                        ){
                        if(found_ele8_CaloIdL_TrkIdVL == false){
                            found_ele8_CaloIdL_TrkIdVL = true;
                            PrintTriggerDebugHeader( Form("%s.txt",triggers.at(i).c_str()) );
                        }
                    }

                    //
                    if(
                        ( strcmp( trigString, "HLT_Ele8_CaloIdL_CaloIsoVL_v1") == 0 ) || 
                        ( strcmp( trigString, "HLT_Ele8_CaloIdL_CaloIsoVL_v2") == 0 )
                        ){
                        if(found_ele8_CaloIdL_CaloIsoVL == false){
                            found_ele8_CaloIdL_CaloIsoVL = true;
                            PrintTriggerDebugHeader( Form("%s.txt",triggers.at(i).c_str()) );
                        }
                    }

                    //
                    if(
                        ( strcmp( trigString, "HLT_Ele17_CaloIdL_CaloIsoVL_v1") == 0 ) ||
                        ( strcmp( trigString, "HLT_Ele17_CaloIdL_CaloIsoVL_v2") == 0 )
                        ){
                        if(found_ele17_CaloIdL_CaloIsoVL == false){
                            found_ele17_CaloIdL_CaloIsoVL = true;
                            PrintTriggerDebugHeader( Form("%s.txt",triggers.at(i).c_str()) );
                        }
                    }


                    //
                    if(
                        ( strcmp( trigString, "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1") == 0 ) || 
                        ( strcmp( trigString, "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2") == 0 )
                        ){
                        if(found_ele8_CaloIdL_CaloIsoVL_Jet40 == false){
                            found_ele8_CaloIdL_CaloIsoVL_Jet40 = true;
                            PrintTriggerDebugHeader( Form("%s.txt",triggers.at(i).c_str()) );
                        }
                    }





                    PrintTriggerDebugLine( itrg, id, match, matchId, dr, lepton_p4, p4tr, trigString, nTrig, Form("%s.txt", triggers.at(i).c_str() ) );
                }
            }

        } // end loop on triggers
        if(matchId){
            nTrig = 3;
        }
        else {
            if (match) {
                nTrig = 2;
            } 
            else {
                nTrig = 1;
            }
        }
    }
    pair<int, float> answer;
    answer.first  = nTrig;
    answer.second = dR_min;
    return answer;
}

// struct for trigger matching
struct triggerMatchStruct {
    triggerMatchStruct (int NumHLTObjs_, float delta_R_, int vers_) {
        nHLTObjects_ = NumHLTObjs_;
        dR_ = delta_R_;
        version_ = vers_;
    }

    int nHLTObjects_;
    float dR_;
    int version_;
};

// wrapper around TriggerMatch that takes a TRegExp for matching a class of triggers
triggerMatchStruct MatchTriggerClass(LorentzVector lepton_p4, const TRegexp trigClass, double dR_cut = 0.4, int pid = 11)
{
    std::pair<int, float> triggerMatchValues = make_pair (0, 99.);
    triggerMatchStruct triggerMatchInfo = triggerMatchStruct(triggerMatchValues.first, triggerMatchValues.second, -1);
    
    unsigned int loopCounts = 0;
    for (unsigned int tidx = 0; tidx < cms2.hlt_trigNames().size(); tidx++) {
        int length = -1;
        if (trigClass.Index(cms2.hlt_trigNames().at(tidx), &length) != 0)
            continue;

        ++loopCounts;

        // get lepton-trigger matching information
        triggerMatchValues = TriggerMatch(lepton_p4, cms2.hlt_trigNames().at(tidx).Data(), dR_cut, pid);        

        // now, figure out which version of the trigger we're dealing with
        const TString delimeter = "_";
        TObjArray* tokens = cms2.hlt_trigNames().at(tidx).Tokenize(delimeter);

        TObjString* last_token = (TObjString*)tokens->Last();
        TString last_stoken(last_token->GetString());
        TSubString last_subtoken = last_stoken.Strip(TString::kLeading, 'v');
        TString subtoken = TString(last_subtoken);

        int version = -1;
        if (subtoken.IsDigit())
            version = subtoken.Atoi();

        triggerMatchInfo = triggerMatchStruct(triggerMatchValues.first, triggerMatchValues.second, version);
    }

    assert (loopCounts < 2);

    return triggerMatchInfo;
}

//////////////////////////////
// THIS NEEDS TO BE IN CORE //
//////////////////////////////

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
        //set_goodrun_file_json("json/Cert_TopNov5_Merged_135821-149442_allPVT.txt");
    }

    // Jet Corrections
    std::vector<std::string> jetcorr_pf_L2L3_filenames;
    jetcorr_pf_L2L3_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L2Relative_AK5PF.txt");
    jetcorr_pf_L2L3_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L3Absolute_AK5PF.txt");
    FactorizedJetCorrector *jet_pf_L2L3corrector = makeJetCorrector(jetcorr_pf_L2L3_filenames);

    std::vector<std::string> jetcorr_pf_L1FastL2L3_filenames;
    //jetcorr_pf_L1FastL2L3_filenames.push_back("../CondFormats/JetMETObjects/data/Jec10V1_L1FastJet_AK5PF.txt");
    jetcorr_pf_L1FastL2L3_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L2Relative_AK5PF.txt");
    jetcorr_pf_L1FastL2L3_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L3Absolute_AK5PF.txt");
    FactorizedJetCorrector *jet_pf_L1FastL2L3corrector = makeJetCorrector(jetcorr_pf_L1FastL2L3_filenames);

    std::vector<std::string> jetcorr_jpt_L2L3_filenames;
    jetcorr_jpt_L2L3_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L2Relative_AK5JPT.txt");
    jetcorr_jpt_L2L3_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L3Absolute_AK5JPT.txt");
    FactorizedJetCorrector *jet_jpt_L2L3corrector = makeJetCorrector(jetcorr_jpt_L2L3_filenames);



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
                //if(!goodrun_json( evt_run(), evt_lumiBlock() )) continue;


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
                float jet_cor = jetCorrection(jp4, jet_pf_L2L3corrector);
                LorentzVector jp4cor = jp4 * jet_cor;
                if (jp4cor.pt() < 15) continue;
                if (pfjets_simpleSecondaryVertexHighEffBJetTag().at(iJet) < 1.74) continue;
                this_nbpfjet++;
                bpfindex.push_back(iJet);
            }
      
// Electrons
      
            if (eormu == -1 || eormu==11) {
                for (unsigned int iLep = 0 ; iLep < els_p4().size(); iLep++) {

                    // Apply a pt cut (Changed it from 5 GeV to 10 GeV...Claudio 10 July 2010)
                    if ( els_p4().at(iLep).pt() < 10.) continue;

                    // Initialize baby ntuple
                    InitBabyNtuple();

                    //////////////////////////////////////////////////////
                    // Fake Rate Numerator & Denominator Selections     //
                    //////////////////////////////////////////////////////

                    //////////
                    // 2011 //
                    //////////

                    // variable naming convention:  [num or vX]_[el or mu]_[analysis][version]

                    // SS
                    num_el_ssV3_    = pass_electronSelection( iLep, electronSelection_ssV3             );
                    v1_el_ssV3_     = pass_electronSelection( iLep, electronSelectionFOV3_ssVBTF80_v1  );
                    v2_el_ssV3_     = pass_electronSelection( iLep, electronSelectionFOV3_ssVBTF80_v2  );
                    v3_el_ssV3_     = pass_electronSelection( iLep, electronSelectionFOV3_ssVBTF80_v3  );

                    // WW
                    num_el_smurfV3_ = pass_electronSelection( iLep, electronSelection_smurfV3          );
                    v1_el_smurfV1_  = pass_electronSelection( iLep, electronSelectionFO_el_smurf_v1    );
                    v2_el_smurfV1_  = pass_electronSelection( iLep, electronSelectionFO_el_smurf_v2    );
                    v3_el_smurfV1_  = pass_electronSelection( iLep, electronSelectionFO_el_smurf_v3    );
                    v4_el_smurfV1_  = pass_electronSelection( iLep, electronSelectionFO_el_smurf_v4    );

                    //////////
                    // 2010 //
                    //////////

                    // ttbar
                    numOct6_ = pass_electronSelection( iLep, electronSelection_ttbarV1_pass5         );
                    v1Oct6_  = pass_electronSelection( iLep, electronSelectionFO_el_ttbarV1_v1_pass5 );
                    v2Oct6_  = pass_electronSelection( iLep, electronSelectionFO_el_ttbarV1_v2_pass5 );
                    v3Oct6_  = pass_electronSelection( iLep, electronSelectionFO_el_ttbarV1_v3_pass5 );

                    // Same Sign Susy
                    numSSV2_ = pass_electronSelection( iLep, electronSelection_ssV2           , false, false ) && (!isSpikeElectron(iLep));
                    v1SSV2_  = pass_electronSelection( iLep, electronSelectionFOV2_ssVBTF80_v1, false, false );
                    v2SSV2_  = pass_electronSelection( iLep, electronSelectionFOV2_ssVBTF80_v2, false, false );
                    v3SSV2_  = pass_electronSelection( iLep, electronSelectionFOV2_ssVBTF80_v3, false, false );

                    // Opposite Sign Susy
                    numOSOct18_ = pass_electronSelection( iLep, electronSelection_el_OSV1);
                    v1OSOct18_  = pass_electronSelection( iLep, electronSelectionFO_el_OSV1_v1);
                    v2OSOct18_  = pass_electronSelection( iLep, electronSelectionFO_el_OSV1_v2);
                    v3OSOct18_  = pass_electronSelection( iLep, electronSelectionFO_el_OSV1_v3);

                    // WW
                    num_wwV1_ = pass_electronSelection( iLep, electronSelection_wwV1);
                    v1_wwV1_  = pass_electronSelection( iLep, electronSelectionFO_el_wwV1_v1);
                    v2_wwV1_  = pass_electronSelection( iLep, electronSelectionFO_el_wwV1_v2);
                    v3_wwV1_  = pass_electronSelection( iLep, electronSelectionFO_el_wwV1_v3);
                    v4_wwV1_  = pass_electronSelection( iLep, electronSelectionFO_el_wwV1_v4);

                    // loosest denominator - used for Z veto
                    // bool loosest = 
                    //   v1Oct6_        | v2Oct6_     | v3Oct6_        |                   // ttbar
                    //   v1SSV2_        | v2SSV2_     | v3SSV2_        |                   // SS 2010
                    //   v1_el_ssV3_    | v2_el_ssV3_ | v3_el_ssV3_    |                   // SS 2011
                    //   v1OSOct18_     | v2OSOct18_  | v3OSOct18_     |                   // OS
                    //   v1_wwV1_       | v2_wwV1_    | v3_wwV1_       | v4_wwV1_        | // WW 2010
                    //   v1_el_smurfV1_ |               v3_el_smurfV1_ | v4_el_smurfV1_  ; // WW 2011


                    ////////////////////////////////////////////////////////////
                    // Skip this electron if it fails the loosest denominator //
                    ////////////////////////////////////////////////////////////

                    // ttbar
                    if (numOct6_ && (!v1Oct6_)) cout << "bad v1Oct6_" << endl;
                    if (numOct6_ && (!v2Oct6_)) cout << "bad v2Oct6_" << endl;
                    if (numOct6_ && (!v3Oct6_)) cout << "bad v3Oct6_" << endl;
  
                    // SS
                    if (numSSV2_ && (!v1SSV2_)) cout << "bad v1SSV2_" << endl;
                    if (numSSV2_ && (!v2SSV2_)) cout << "bad v2SSV2_" << endl;
                    if (numSSV2_ && (!v3SSV2_)) cout << "bad v3SSV2_" << endl;
  
                    // OS
                    if (numOSOct18_ && (!v1OSOct18_)) cout << "bad v1OSOct18_" << endl;
                    if (numOSOct18_ && (!v2OSOct18_)) cout << "bad v2OSOct18_" << endl;
                    if (numOSOct18_ && (!v3OSOct18_)) cout << "bad v3OSOct18_" << endl;
    
                    // WW
                    if (num_wwV1_ && (!v1_wwV1_)) cout << "bad v1_wwV1_" << endl;
                    if (num_wwV1_ && (!v2_wwV1_)) cout << "bad v2_wwV1_" << endl;
                    if (num_wwV1_ && (!v3_wwV1_)) cout << "bad v3_wwV1_" << endl;
                    if (num_wwV1_ && (!v4_wwV1_)) cout << "bad v4_wwV1_" << endl;

                    // 
                    if (  
                        (!v1Oct6_)        && (!v2Oct6_)        && (!v3Oct6_)        &&                      // ttbar
                        (!v1SSV2_)        && (!v2SSV2_)        && (!v3SSV2_)        &&                      // SS 2010
                        (!v1_el_ssV3_)    && (!v2_el_ssV3_)    && (!v3_el_ssV3_)    &&                      // SS 2011
                        (!v1OSOct18_)     && (!v2OSOct18_)     && (!v3OSOct18_)     &&                      // OS
                        (!v1_wwV1_)       && (!v2_wwV1_)       && (!v3_wwV1_)       && (!v4_wwV1_)       && // WW 2010
                        (!v1_el_smurfV1_) && (!v1_el_smurfV1_) && (!v3_el_smurfV1_) && (!v4_el_smurfV1_)    // WW 2011
                        ){ 
                        //continue;
                    }
 
                    //////////////////////////////////////////////////////
                    // End Fake Rate Numerator & Denominator Selections //
                    //////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
// NEED TO THINK ABOUT THIS... Z'S ARE VETOED BASED ON TOP SELECTIONS //
////////////////////////////////////////////////////////////////////////
       
                    // If it is above 20 GeV see if we can make a 
                    // Z with another pt>20 FO.  Will use the v1 FO since 
                    // these are the loosest
                    bool isaZ = false;
                    if (els_p4().at(iLep).pt() > 20.) {
                        for (unsigned int jEl = 0 ; jEl < els_p4().size(); jEl++) {
                            if (iLep == jEl)                             continue;
                            if (els_p4().at(jEl).pt() < 20.)            continue;
                            if ( ! pass_electronSelection( jEl, electronSelectionFO_el_ttbarV1_v1 ) ) continue;
                            if ( ! v1Oct6_ ) continue;
                            LorentzVector w = els_p4().at(iLep) + els_p4().at(jEl);
                            if (abs(w.mass()-91.) > 20.) continue;
                            isaZ = true;
                        }
                    }
                    if (isaZ) continue;
    

                    /////////////////////////// 
                    // Event Information     //
                    ///////////////////////////

                    // Load the electron and event quantities
                    run_   = evt_run();
                    ls_    = evt_lumiBlock();
                    evt_   = evt_event();
                    weight_ = isData ? 1. : evt_scale1fb();
  
                    if(!isData){
                        // Pileup - PUSummaryInfoMaker
                        pu_nPUvertices_ = puInfo_nPUvertices();
                    } else {
  
                        // Pileup - VertexMaker
                        evt_nvtxs_       = evt_nvtxs();
    
                        // Pileup - VertexMaker
                        evt_ndavtxs_       = evt_ndavtxs();

                    }

                    /////////////////////////// 
                    // End Event Information //
                    ///////////////////////////



                    //////////////////////////// 
                    // Lepton Information     //
                    ////////////////////////////

                    // Basic Quantities
                    pt_             = els_p4().at(iLep).pt();
                    eta_            = els_p4().at(iLep).eta();
                    phi_            = els_p4().at(iLep).phi();
                    scet_           = els_eSC()[iLep] / cosh( els_etaSC()[iLep] );
                    id_             = 11*els_charge().at(iLep);
                    tcmet_          = evt_tcmet();
                    tcmetphi_       = evt_tcmetPhi();
                    pfmet_          = evt_pfmet();
                    pfmetphi_       = evt_pfmetPhi();

                    // Isolation
                    //iso_              = electronIsolation_rel         ( iLep, true               );
                    //iso_nps_          = electronIsolation_rel         ( iLep, true, 0.0          );
                    //nt_iso_           = electronIsolation_rel_v1      ( iLep, true, false        );
                    //nt_iso_nps_       = electronIsolation_rel_v1      ( iLep, true,  false, 0.0  );
                    //trck_iso_         = electronIsolation_rel         ( iLep, false              );
                    //trck_nt_iso_      = electronIsolation_rel_v1      ( iLep, false, false       );
                    //ecal_iso_         = electronIsolation_ECAL_rel    ( iLep                     );
                    //ecal_iso_nps_     = electronIsolation_ECAL_rel    ( iLep, 0.0                );
                    //ecal_nt_iso_      = electronIsolation_ECAL_rel_v1 ( iLep, false              );
                    //ecal_nt_iso_nps_  = electronIsolation_ECAL_rel_v1 ( iLep, false, 0.0         );
                    //hcal_iso_         = electronIsolation_HCAL_rel    ( iLep                     );
                    //hcal_nt_iso_      = electronIsolation_HCAL_rel_v1 ( iLep, false              );

                    // Isolation
                    iso_              = electronIsolation_rel         ( iLep, true  );
                    iso_nps_          = electronIsolation_rel         ( iLep, true  );  // wrong
                    nt_iso_           = electronIsolation_rel_v1      ( iLep, true  );
                    nt_iso_nps_       = electronIsolation_rel_v1      ( iLep, true  );

                    trck_iso_         = electronIsolation_rel         ( iLep, false );
                    trck_nt_iso_      = electronIsolation_rel_v1      ( iLep, false );

                    ecal_iso_         = electronIsolation_ECAL_rel    ( iLep        );
                    ecal_iso_nps_     = electronIsolation_ECAL_rel    ( iLep        );  // wrong
                    ecal_nt_iso_      = electronIsolation_ECAL_rel_v1 ( iLep, true  );
                    ecal_nt_iso_nps_  = electronIsolation_ECAL_rel_v1 ( iLep, false );
                    hcal_iso_         = electronIsolation_HCAL_rel    ( iLep );
                    hcal_nt_iso_      = electronIsolation_HCAL_rel_v1 ( iLep );
          

                    // safety
                    /*
                      Float_t x1 = electronIsolation_relOriginal(iLep, true);    // Original isolation, truncated
                      Float_t x2 = electronIsolation_rel_v1(iLep, true, true);   // new function, truncated
                      if( fabs(x1 - x2) > .000001 ){
                      //cout << "T Iso: " << x1 << " != " << x2 << endl;
                      }

                      Float_t y1 = electronIsolation_rel_v1Original(iLep, true);  // Original isolation, not truncated
                      Float_t y2 = electronIsolation_rel_v1(iLep, true, false);   // new function, non truncated
                      if( fabs(y1 - y2) > .000001 ){
                      //cout << "NT Iso: " << y1 << " != " << y2 << endl;
                      }
                    */

                    el_id_smurfV3_  = pass_electronSelection( iLep, electronSelection_smurfV3_id );
                    el_id_vbtf80_   = electronId_VBTF(iLep, VBTF_35X_80, false, false);
                    el_id_vbtf90_   = electronId_VBTF(iLep, VBTF_35X_90, false, false);
                    if( els_closestMuon().at(iLep) == -1 ) closestMuon_ = true;
                    if (! isData) {
                        mcid_       = els_mc_id().at(iLep);
                        mcmotherid_ = els_mc_motherid().at(iLep);
                    }
      
                    // PV
                    d0PV_wwV1_ = electron_d0PV_wwV1(iLep);
                    dzPV_wwV1_ = electron_dzPV_wwV1(iLep);

                    // W transverse mass
                    mt_   = Mt( els_p4().at(iLep), tcmet_, tcmetphi_ );
                    pfmt_ = Mt( els_p4().at(iLep), pfmet_, pfmetphi_ );

                    // Do the 3 electron charges agree?
                    int iCTF = els_trkidx().at(iLep);
                    if( iCTF >= 0 ){
                        int qCTF = trks_charge().at( iCTF );
                        int qGSF = els_trk_charge().at(iLep);
                        int qPIX = els_sccharge().at(iLep);
                        if( qCTF == qGSF && qCTF == qPIX && qGSF == qPIX ) q3_ = true;
                    }
  
                    // Missing hits info
                    els_exp_innerlayers_ = els_exp_innerlayers().at(iLep);

                    // Conversion Rejection  
                    convHitPattern_   = isFromConversionHitPattern(iLep);
                    convPartnerTrack_ = isFromConversionPartnerTrack(iLep);
                    convMIT_          = isFromConversionMIT(iLep);
                    if( els_exp_innerlayers().at(iLep) == 0 ) conv0MissHits_ = true;

                    // HT
                    ht_calo_           = (float) sumPt (iLep, JETS_TYPE_CALO_UNCORR  , JETS_CLEAN_SINGLE_E );
                    ht_calo_L2L3_      = (float) sumPt (iLep, JETS_TYPE_CALO_CORR    , JETS_CLEAN_SINGLE_E );
                    ht_jpt_L2L3_       = (float) sumPt (iLep, JETS_TYPE_JPT          , JETS_CLEAN_SINGLE_E );
                    ht_pf_             = (float) sumPt (iLep, JETS_TYPE_PF_UNCORR    , JETS_CLEAN_SINGLE_E );
                    ht_pf_L2L3_        = (float) sumPt (iLep, JETS_TYPE_PF_CORR      , JETS_CLEAN_SINGLE_E );
                    //ht_pf_L1FastL2L3_  = (float) sumPt (iLep, JETS_TYPE_PF_FAST_CORR , JETS_CLEAN_SINGLE_E );

                    //////////////////////////// 
                    // End Lepton Information //
                    ////////////////////////////



                    ///////////////////////  
                    // 2011 Triggers     //
                    ///////////////////////

                    // Electrons
                    triggerMatchStruct struct_ele8_vstar                                           = MatchTriggerClass( els_p4().at(iLep), TRegexp("^HLT_Ele8_v[0-9]+$"                                         ));
                    triggerMatchStruct struct_ele8_CaloIdL_TrkIdVL_vstar                           = MatchTriggerClass( els_p4().at(iLep), TRegexp("^HLT_Ele8_CaloIdL_TrkIdVL_v[0-9]+$"                         ));
                    triggerMatchStruct struct_ele8_CaloIdL_CaloIsoVL_Jet40_vstar                   = MatchTriggerClass( els_p4().at(iLep), TRegexp("^HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v[0-9]+$"                 ));
                    triggerMatchStruct struct_ele8_CaloIdL_CaloIsoVL_vstar                         = MatchTriggerClass( els_p4().at(iLep), TRegexp("^HLT_Ele8_CaloIdL_CaloIsoVL_v[0-9]+$"                       ));
                    triggerMatchStruct struct_ele17_CaloIdL_CaloIsoVL_vstar                        = MatchTriggerClass( els_p4().at(iLep), TRegexp("^HLT_Ele17_CaloIdL_CaloIsoVL_v[0-9]+$"                      ));  
                    //struct<int, float> struct_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar  = MatchTriggerClass( els_p4().at(iLep), TRegexp("^HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v[0-9]+$"));
  
                    ele8_vstar_                                              = struct_ele8_vstar.nHLTObjects_;
                    ele8_CaloIdL_TrkIdVL_vstar_                              = struct_ele8_CaloIdL_TrkIdVL_vstar.nHLTObjects_; 
                    ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                      = struct_ele8_CaloIdL_CaloIsoVL_Jet40_vstar.nHLTObjects_;
                    ele8_CaloIdL_CaloIsoVL_vstar_                            = struct_ele8_CaloIdL_CaloIsoVL_vstar.nHLTObjects_;
                    ele17_CaloIdL_CaloIsoVL_vstar_                           = struct_ele17_CaloIdL_CaloIsoVL_vstar.nHLTObjects_; 
                    //photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_     = struct_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar.nHLTObjects_;

                    ele8_version_                                              = struct_ele8_vstar.nHLTObjects_;
                    ele8_CaloIdL_TrkIdVL_version_                              = struct_ele8_CaloIdL_TrkIdVL_vstar.nHLTObjects_; 
                    ele8_CaloIdL_CaloIsoVL_Jet40_version_                      = struct_ele8_CaloIdL_CaloIsoVL_Jet40_vstar.nHLTObjects_;
                    ele8_CaloIdL_CaloIsoVL_version_                            = struct_ele8_CaloIdL_CaloIsoVL_vstar.nHLTObjects_;
                    ele17_CaloIdL_CaloIsoVL_version_                           = struct_ele17_CaloIdL_CaloIsoVL_vstar.nHLTObjects_; 
                    //photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_version_     = struct_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar.nHLTObjects_;
          
                    dr_ele8_vstar_                                           = struct_ele8_vstar.dR_;
                    dr_ele8_CaloIdL_TrkIdVL_vstar_                           = struct_ele8_CaloIdL_TrkIdVL_vstar.dR_;
                    dr_ele8_CaloIdL_CaloIsoVL_Jet40_vstar_                   = struct_ele8_CaloIdL_CaloIsoVL_Jet40_vstar.dR_;
                    dr_ele8_CaloIdL_CaloIsoVL_vstar_                         = struct_ele8_CaloIdL_CaloIsoVL_vstar.dR_;
                    dr_ele17_CaloIdL_CaloIsoVL_vstar_                        = struct_ele17_CaloIdL_CaloIsoVL_vstar.dR_;
                    //dr_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar_  = struct_photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_vstar.dR_; 


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
                        double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), p4j);
                        if (dr > deltaRCut) hlt15u_ = 2;
                    }
                    if (hlt30u_ == 1) {
                        LorentzVector p4j = p4HLTObject("HLT_Jet30U",0);
                        double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), p4j);
                        if (dr > deltaRCut) hlt30u_ = 2;
                    }
                    if (hlt50u_ == 1) {
                        LorentzVector p4j = p4HLTObject("HLT_Jet50U",0);
                        double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), p4j);
                        if (dr > deltaRCut) hlt50u_ = 2;
                    }
                    if (l16u_ == 1) {
                        LorentzVector p4j = p4HLTObject("HLT_L1Jet6U",0);
                        double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), p4j);
                        if (dr > deltaRCut) l16u_ = 2;
                    }
                    if (l110u_ == 1) {
                        LorentzVector p4j = p4HLTObject("HLT_L1Jet10U",0);
                        double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), p4j);
                        if (dr > deltaRCut) l110u_ = 2;
                    }
        
                    // Fill electron & photon triggers
                    pair<int, float> pair_el10_lw       = TriggerMatch( els_p4().at(iLep), "HLT_Ele10_LW_L1R"                  );
                    pair<int, float> pair_el10_sw       = TriggerMatch( els_p4().at(iLep), "HLT_Ele10_SW_L1R"                  );
                    pair<int, float> pair_el10_sw_v2    = TriggerMatch( els_p4().at(iLep), "HLT_Ele10_SW_L1R_v2"               );
                    pair<int, float> pair_el10_lw_id    = TriggerMatch( els_p4().at(iLep), "HLT_Ele10_LW_EleId_L1R"            );
                    pair<int, float> pair_el10_sw_id    = TriggerMatch( els_p4().at(iLep), "HLT_Ele10_SW_EleId_L1R"            );
                    pair<int, float> pair_el15_lw       = TriggerMatch( els_p4().at(iLep), "HLT_Ele15_LW_L1R"                  );
                    pair<int, float> pair_el15_sw       = TriggerMatch( els_p4().at(iLep), "HLT_Ele15_SW_L1R"                  );
                    pair<int, float> pair_el15_lw_id    = TriggerMatch( els_p4().at(iLep), "HLT_Ele15_LW_EleId_L1R"            );
                    pair<int, float> pair_el15_sw_id    = TriggerMatch( els_p4().at(iLep), "HLT_Ele15_SW_EleId_L1R"            );
                    pair<int, float> pair_el15_sw_cid   = TriggerMatch( els_p4().at(iLep), "HLT_Ele15_SW_CaloEleId_L1R"        );
                    pair<int, float> pair_el20_sw       = TriggerMatch( els_p4().at(iLep), "HLT_Ele20_SW_L1R"                  );
                    pair<int, float> pair_el25_sw       = TriggerMatch( els_p4().at(iLep), "HLT_Ele25_SW_L1R"                  );
                    pair<int, float> pair_Del10_sw      = TriggerMatch( els_p4().at(iLep), "HLT_DoubleEle10_SW_L1R"            );
                    pair<int, float> pair_ph10          = TriggerMatch( els_p4().at(iLep), "HLT_Photon10_L1R"                  );
                    pair<int, float> pair_ph10C         = TriggerMatch( els_p4().at(iLep), "HLT_Photon10_Cleaned_L1R"          );
                    pair<int, float> pair_ph15          = TriggerMatch( els_p4().at(iLep), "HLT_Photon15_L1R"                  );
                    pair<int, float> pair_ph15C         = TriggerMatch( els_p4().at(iLep), "HLT_Photon15_Cleaned_L1R"          );
                    pair<int, float> pair_ph20C         = TriggerMatch( els_p4().at(iLep), "HLT_Photon20_Cleaned_L1R"          );
                    pair<int, float> pair_el17_sw       = TriggerMatch( els_p4().at(iLep), "HLT_Ele17_SW_L1R"                  );
                    pair<int, float> pair_el17_sw_v2    = TriggerMatch( els_p4().at(iLep), "HLT_Ele17_SW_L1R_v2"               );
                    pair<int, float> pair_el17_iso      = TriggerMatch( els_p4().at(iLep), "HLT_Ele17_Isol_L1R"                );
                    pair<int, float> pair_el17_loose    = TriggerMatch( els_p4().at(iLep), "HLT_Ele17_SW_LooseEleId_L1R"       );
                    pair<int, float> pair_el17_sw_cid   = TriggerMatch( els_p4().at(iLep), "HLT_Ele17_SW_CaloEleId_L1R"        );
                    pair<int, float> pair_el17_sw_id    = TriggerMatch( els_p4().at(iLep), "HLT_Ele17_SW_EleId_L1R"            );
                    pair<int, float> pair_el17_tiso     = TriggerMatch( els_p4().at(iLep), "HLT_Ele17_SW_TightEleIdIsol_L1R"   );
                    pair<int, float> pair_el17_tiso_v1  = TriggerMatch( els_p4().at(iLep), "HLT_Ele17_SW_TightEleIdIsol_L1R_v1");

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
                    ph10_         = ph10;
                    ph15_         = ph15;
                    ph20_         = pair_ph20C.first;

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
                    drph10_         = drph10;
                    drph15_         = drph15;
                    drph20_         = pair_ph20C.second;

                    ///////////////////////  
                    // End 2010 Triggers //
                    ///////////////////////



                    //////////////
                    // Jets     //
                    //////////////

                    // Calo Jets
                    // Find the highest Pt jet separated by at least dRcut from this lepton and fill the jet Pt
                    ptj1_       = -999.0;
                    ptj1_b2b_   = -999.0;
                    dphij1_b2b_ = -999.0;
                    nj1_        = 0;
                    for (unsigned int iJet = 0; iJet < jets_p4().size(); iJet++) {
                        double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jets_p4().at(iJet) );
                        if( dr > deltaRCut && jets_p4().at(iJet).pt() > 10 ) nj1_++;
                        if ( dr > deltaRCut && jets_p4().at(iJet).pt() > ptj1_ ){
                            ptj1_ = jets_p4().at(iJet).pt();
          
                            // back to back in phi
                            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iLep), jets_p4().at(iJet) ) );
                            if( dphi > deltaPhiCut && jets_p4().at(iJet).pt() > ptj1_b2b_ ){ 
                                ptj1_b2b_   = jets_p4().at(iJet).pt();
                                dphij1_b2b_ = dphi;
                            }
                        }
                    }
    
                    // PF Jets
                    // Find the highest Pt pfjet separated by at least dRcut from this lepton and fill the pfjet Pt
                    ptpfj1_       = -999.0;
                    ptpfj1_b2b_   = -999.0;
                    dphipfj1_b2b_ = -999.0;
                    npfj1_        = 0;
                    for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                        double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), pfjets_p4().at(iJet) );
                        if( dr > deltaRCut && pfjets_p4().at(iJet).pt() > 10 ) npfj1_++;
                        if ( dr > deltaRCut && pfjets_p4().at(iJet).pt() > ptpfj1_ ){
                            ptpfj1_ = pfjets_p4().at(iJet).pt();
          
                            // back to back in phi
                            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iLep), pfjets_p4().at(iJet) ) );
                            if( dphi > deltaPhiCut && pfjets_p4().at(iJet).pt() > ptpfj1_b2b_ ){ 
                                ptpfj1_b2b_   = pfjets_p4().at(iJet).pt();
                                dphij1_b2b_ = dphi;
                            }
                        }
                    }
    
                    // L2L3 PF Jets
                    // Find the highest Pt PF L2L3 corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                    ptpfcj1_       = -999.0; 
                    ptpfcj1_b2b_   = -999.0;
                    dphipfcj1_b2b_ = -999.0;
                    npfcj1_        = 0;
                    btagpfc_       = false;
                    for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                        if ( !passesPFJetID(iJet)) continue;
                        LorentzVector jp4 = pfjets_p4()[iJet];
                        float jet_cor = jetCorrection(jp4, jet_pf_L2L3corrector);
                        LorentzVector jp4cor = jp4 * jet_cor;
                        if (jp4cor.pt() > 15 && pfjets_simpleSecondaryVertexHighEffBJetTag().at(iJet) > 1.74 ) btagpfc_ = true;
                        double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jp4cor );
                        if( dr > deltaRCut && jp4cor.pt() > 10 ) npfcj1_++;
                        if ( dr > deltaRCut && jp4cor.pt() > ptpfcj1_ ){
                            ptpfcj1_ = jp4cor.pt();
    
                            // back to back in phi
                            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iLep), jp4cor ) );
                            if( dphi > deltaPhiCut && jp4cor.pt() > ptpfcj1_b2b_ ){
                                ptpfcj1_b2b_   = jp4cor.pt();
                                dphipfcj1_b2b_ = dphi;
                            } 
                        }
                    }

                    // L1FastL2L3 PF Jets
                    // Find the highest Pt PF L1FastL2L3 corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                    ptpfcL1Fj1_       = -999.0;
                    ptpfcL1Fj1_b2b_   = -999.0;
                    dphipfcL1Fj1_b2b_ = -999.0;
                    npfcL1Fj1_        = 0;
                    btagpfcL1F_       = false;
                    for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                        if ( !passesPFJetID(iJet)) continue;
                        LorentzVector jp4 = pfjets_p4()[iJet];
                        float jet_cor = jetCorrection(jp4, jet_pf_L1FastL2L3corrector);
                        LorentzVector jp4cor = jp4 * jet_cor;
                        if (jp4cor.pt() > 15 && pfjets_simpleSecondaryVertexHighEffBJetTag().at(iJet) > 1.74 ) btagpfcL1F_ = true;
                        double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jp4cor );
                        if( dr > deltaRCut && jp4cor.pt() > 10 ) npfcL1Fj1_++;
                        if ( dr > deltaRCut && jp4cor.pt() > ptpfcL1Fj1_ ){
                            ptpfcL1Fj1_ = jp4cor.pt();

                            // back to back in phi
                            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iLep), jp4cor ) );
                            if( dphi > deltaPhiCut && jp4cor.pt() > ptpfcL1Fj1_b2b_ ){
                                ptpfcL1Fj1_b2b_   = jp4cor.pt();
                                dphipfcL1Fj1_b2b_ = dphi;
                            }
                        }
                    }

                    // L2L3 JPT Jets
                    // Find the highest Pt JPT L2L3 corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                    ptjptcj1_       = -999.0;
                    ptjptcj1_b2b_   = -999.0;
                    dphijptcj1_b2b_ = -999.0;
                    njptcj1_        = 0;
                    btagjptc_       = false;
                    for (unsigned int iJet = 0; iJet < jpts_p4().size(); iJet++) {
                        LorentzVector jp4 = jpts_p4()[iJet];
                        if ( !passesCaloJetID(jp4)) continue;
                        float jet_cor = jetCorrection(jp4, jet_jpt_L2L3corrector);
                        LorentzVector jp4cor = jp4 * jet_cor;
                        if (jp4cor.pt() > 15 && jpts_simpleSecondaryVertexHighEffBJetTag().at(iJet) > 1.74 ) btagjptc_ = true;
                        double dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jp4cor );
                        if( dr > deltaRCut && jp4cor.pt() > 10 ) njptcj1_++;
                        if ( dr > deltaRCut && jp4cor.pt() > ptjptcj1_ ){
                            ptjptcj1_ = jp4cor.pt();
   
                            // back to back in phi
                            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( els_p4().at(iLep), jp4cor ) );
                            if( dphi > deltaPhiCut && jp4cor.pt() > ptjptcj1_b2b_ ){
                                ptjptcj1_b2b_   = jp4cor.pt();
                                dphijptcj1_b2b_ = dphi;
                            }
                        }
                    }

                    //////////////
                    // End Jets //
                    //////////////



                    ///////////////////
                    // B Tagging     //
                    ///////////////////

                    // The btag information
                    nbjet_ = this_nbjet;
                    dRbNear_ = 99.;
                    dRbFar_  = -99.;
                    for (int ii=0; ii<nbjet_; ii++) {
                        unsigned int iJet = bindex[ii];
                        float dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jets_p4().at(iJet));
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
                        float jet_cor = jetCorrection(jp4, jet_pf_L2L3corrector);
                        LorentzVector jp4cor = jp4 * jet_cor;
                        float dr = ROOT::Math::VectorUtil::DeltaR( els_p4().at(iLep), jp4cor);
                        if (dr < dRbpfcNear_) dRbpfcNear_ = dr;
                        if (dr > dRbpfcFar_)   dRbpfcFar_  = dr;
                    }

                    ///////////////////
                    // End B Tagging //
                    ///////////////////



                    // Time to fill the baby for the electrons
                    FillBabyNtuple();

                } // closes loop over electrons
            } // closes if statements about whether we want to fill electrons

// Muons

            if (eormu == -1 || eormu==13) {
                for ( unsigned int iLep = 0; iLep < mus_p4().size(); iLep++) {
        
                    // Apply a pt cut --- moved the cut to 10 GeV (Claudio, 10 Jul 2010)
                    if ( mus_p4().at(iLep).pt() < 5.0) continue;
        
////////////////////////////////////////////////////////////////////////
// NEED TO THINK ABOUT THIS... Z'S ARE VETOED BASED ON TOP SELECTIONS //
////////////////////////////////////////////////////////////////////////

                    // If it is above 20 GeV see if we can make a 
                    // Z with another pt>20 FO.  
                    bool isaZ = false;
                    if (mus_p4().at(iLep).pt() > 20.) {
                        for (unsigned int jMu = 0 ; jMu < mus_p4().size(); jMu++) {
                            if (iLep == jMu)                             continue;
                            if (mus_p4().at(jMu).pt() < 20.)            continue;
                            if ( ! muonId( jMu, muonSelectionFO_mu_ttbar) ) continue;
                            if ( ! muonId( iLep, muonSelectionFO_mu_ttbar) ) continue;
                            LorentzVector w = mus_p4().at(iLep) + mus_p4().at(jMu);
                            if (abs(w.mass()-91.) > 20.) continue;
                            isaZ = true;
                        }
                    }
                    if (isaZ) continue;
        
                    // Initialize baby ntuple
                    InitBabyNtuple();
 
       

                    /////////////////////////// 
                    // Event Information     //
                    ///////////////////////////

                    // Load the electron and event quantities
                    run_   = evt_run();
                    ls_    = evt_lumiBlock();
                    evt_   = evt_event();
                    weight_ = isData ? 1. : evt_scale1fb();

                    if(!isData){
                        // Pileup - PUSummaryInfoMaker
                        pu_nPUvertices_ = puInfo_nPUvertices();
                    } 
                    else {
                        // Pileup - VertexMaker
                        evt_nvtxs_       = evt_nvtxs();
  
                        // Pileup - VertexMaker
                        evt_ndavtxs_       = evt_ndavtxs();
                    }

                    /////////////////////////// 
                    // End Event Information //
                    ///////////////////////////



                    //////////////////////////// 
                    // Lepton Information     //
                    ////////////////////////////

                    // Basic Quantities
                    pt_       = mus_p4().at(iLep).pt();
                    eta_      = mus_p4().at(iLep).eta();
                    phi_      = mus_p4().at(iLep).phi();
                    id_       = 13*mus_charge().at(iLep);
                    tcmet_    = evt_tcmet();
                    tcmetphi_ = evt_tcmetPhi();
                    pfmet_    = evt_pfmet();
                    pfmetphi_ = evt_pfmetPhi();
          
                    // Isolation
                    iso_          = muonIsoValue(iLep);
                    nt_iso_       = muonIsoValue(iLep, false);
                    trck_iso_     = muonIsoValue_TRK(iLep);
                    trck_nt_iso_  = muonIsoValue_TRK(iLep, false);
                    ecal_iso_     = muonIsoValue_ECAL(iLep);
                    ecal_nt_iso_  = muonIsoValue_ECAL(iLep, false);
                    hcal_iso_     = muonIsoValue_HCAL(iLep);
                    hcal_nt_iso_  = muonIsoValue_HCAL(iLep, false);

                    // Safety
                    Float_t x = muonIsoValueOriginal(iLep);
                    if( x != iso_ ){
                        //cout << "Muons " << x << " != " << iso_ << endl;
                    }

                    if (! isData) {
                        mcid_       = mus_mc_id().at(iLep);
                        mcmotherid_ = mus_mc_motherid().at(iLep);
                    }

                    // Correct tcmet (July 14, 2010)
                    if (pt_ > 10.) {
                        if (!wasMetCorrectedForThisMuon(iLep, usingTcMet)) {
                            float metX = tcmet_ * cos(evt_tcmetPhi());
                            float metY = tcmet_ * sin(evt_tcmetPhi());
                            fixMetForThisMuon(iLep, metX, metY, usingTcMet);
                            tcmet_ = sqrt(metX*metX + metY*metY);
                        }
                    }

                    // W transverse mass
                    mt_   = Mt( mus_p4().at(iLep), tcmet_, tcmetphi_ );
                    pfmt_ = Mt( mus_p4().at(iLep), pfmet_, pfmetphi_ );

                    // HT
                    ht_calo_           = (float) sumPt (iLep, JETS_TYPE_CALO_UNCORR  , JETS_CLEAN_SINGLE_MU );
                    ht_calo_L2L3_      = (float) sumPt (iLep, JETS_TYPE_CALO_CORR    , JETS_CLEAN_SINGLE_MU );
                    ht_jpt_L2L3_       = (float) sumPt (iLep, JETS_TYPE_JPT          , JETS_CLEAN_SINGLE_MU );
                    ht_pf_             = (float) sumPt (iLep, JETS_TYPE_PF_UNCORR    , JETS_CLEAN_SINGLE_MU );
                    ht_pf_L2L3_        = (float) sumPt (iLep, JETS_TYPE_PF_CORR      , JETS_CLEAN_SINGLE_MU );
                    //ht_pf_L1FastL2L3_  = (float) sumPt (iLep, JETS_TYPE_PF_FAST_CORR , JETS_CLEAN_SINGLE_MU );

                    //////////////////////////// 
                    // End Lepton Information //
                    ////////////////////////////



                    //////////////////////////////////////////////////////
                    // Fake Rate Numerator & Denominator Selections     //
                    //////////////////////////////////////////////////////


                    //////////
                    // 2011 //
                    //////////
                    numNomSSv3_   = muonId(iLep, NominalSSv3          );
                    fo_mussV3_04_ = muonId(iLep, muonSelectionFO_ssV3 );

                    //////////
                    // 2010 //
                    //////////

                    // ttbar
                    num_            = muonId(iLep, NominalTTbarV2                    );
                    fo_04_          = muonId(iLep, muonSelectionFO_mu_ttbar          );
                    fo_10_          = muonId(iLep, muonSelectionFO_mu_ttbar_iso10    );
  
                    // Same Sign Susy
                    numNomSSv2_     = muonId(iLep, NominalSSv2                       );
                    fo_mussV2_04_   = muonId(iLep, muonSelectionFO_mu_ssV2           );
                    fo_mussV2_10_   = muonId(iLep, muonSelectionFO_mu_ssV2_iso10     );
  
                    // Opposite Sign Susy
                    num_OSGv1_      = muonId(iLep, OSGeneric_v1                      );
                    num_OSZv1_      = muonId(iLep, OSZ_v1                            );
  
                    // WW
                    num_wwV1_       = muonId(iLep, NominalWWV1                       );
                    fo_wwV1_04_     = muonId(iLep, muonSelectionFO_mu_wwV1           );
                    fo_wwV1_10_     = muonId(iLep, muonSelectionFO_mu_wwV1_iso10     );
                    fo_wwV1_10_d0_  = muonId(iLep, muonSelectionFO_mu_wwV1_iso10_d0  );
 
                    // WW
                    num_mu_smurfV3_    = muonId(iLep, NominalSmurfV3                    );
                    fo_mu_smurf_04_    = muonId(iLep, muonSelectionFO_mu_smurf_04       );
                    fo_mu_smurf_10_    = muonId(iLep, muonSelectionFO_mu_smurf_10       );

 
                    // 
                    if( 
                        !fo_04_ && !fo_10_ &&                                      // ttbar
                        !fo_mussV2_04_   && !fo_mussV2_10_   &&                    // SS
                        !fo_wwV1_04_     && !fo_wwV1_10_     && !fo_wwV1_10_d0_ && // WW
                        !fo_mu_smurf_04_ && !fo_mu_smurf_10_                       // WW
                        ){ 
                        //continue;
                    }
  
                    //////////////////////////////////////////////////////
                    // End Fake Rate Numerator & Denominator Selections //
                    //////////////////////////////////////////////////////



                    ///////////////////////  
                    // 2011 Triggers     //
                    ///////////////////////
  
                    // Muons
                    triggerMatchStruct struct_mu3_vstar       = MatchTriggerClass( mus_p4().at(iLep), TRegexp("^HLT_Mu3_v[0-9]+$"       ));
                    triggerMatchStruct struct_mu5_vstar       = MatchTriggerClass( mus_p4().at(iLep), TRegexp("^HLT_Mu5_v[0-9]+$"       ));
                    triggerMatchStruct struct_mu8_vstar       = MatchTriggerClass( mus_p4().at(iLep), TRegexp("^HLT_Mu8_v[0-9]+$"       ));
                    triggerMatchStruct struct_mu12_vstar      = MatchTriggerClass( mus_p4().at(iLep), TRegexp("^HLT_Mu12_v[0-9]+$"      ));
                    triggerMatchStruct struct_mu15_vstar      = MatchTriggerClass( mus_p4().at(iLep), TRegexp("^HLT_Mu15_v[0-9]+$"      ));
                    triggerMatchStruct struct_mu20_vstar      = MatchTriggerClass( mus_p4().at(iLep), TRegexp("^HLT_Mu20_v[0-9]+$"      ));
                    triggerMatchStruct struct_mu24_vstar      = MatchTriggerClass( mus_p4().at(iLep), TRegexp("^HLT_Mu24_v[0-9]+$"      ));
                    triggerMatchStruct struct_mu30_vstar      = MatchTriggerClass( mus_p4().at(iLep), TRegexp("^HLT_Mu30_v[0-9]+$"      ));
                    triggerMatchStruct struct_mu8_Jet40_vstar = MatchTriggerClass( mus_p4().at(iLep), TRegexp("^HLT_Mu8_Jet40_v[0-9]+$" ));

                    mu3_vstar_         = struct_mu3_vstar.nHLTObjects_;
                    mu5_vstar_         = struct_mu5_vstar.nHLTObjects_;
                    mu8_vstar_         = struct_mu8_vstar.nHLTObjects_;
                    mu12_vstar_        = struct_mu12_vstar.nHLTObjects_;
                    mu15_vstar_        = struct_mu15_vstar.nHLTObjects_;
                    mu20_vstar_        = struct_mu20_vstar.nHLTObjects_;
                    mu24_vstar_        = struct_mu24_vstar.nHLTObjects_;
                    mu30_vstar_        = struct_mu30_vstar.nHLTObjects_;
                    mu8_Jet40_vstar_   = struct_mu8_Jet40_vstar.nHLTObjects_;

                    mu3_version_         = struct_mu3_vstar.version_;
                    mu5_version_         = struct_mu5_vstar.version_;
                    mu8_version_         = struct_mu8_vstar.version_;
                    mu12_version_        = struct_mu12_vstar.version_;
                    mu15_version_        = struct_mu15_vstar.version_;
                    mu20_version_        = struct_mu20_vstar.version_;
                    mu24_version_        = struct_mu24_vstar.version_;
                    mu30_version_        = struct_mu30_vstar.version_;
                    mu8_Jet40_version_   = struct_mu8_Jet40_vstar.version_;

                    dr_mu3_vstar_       = struct_mu3_vstar.dR_;
                    dr_mu5_vstar_       = struct_mu5_vstar.dR_;
                    dr_mu8_vstar_       = struct_mu8_vstar.dR_;
                    dr_mu12_vstar_      = struct_mu12_vstar.dR_;
                    dr_mu15_vstar_      = struct_mu15_vstar.dR_;
                    dr_mu20_vstar_      = struct_mu20_vstar.dR_;
                    dr_mu24_vstar_      = struct_mu24_vstar.dR_;
                    dr_mu30_vstar_      = struct_mu30_vstar.dR_;
                    dr_mu8_Jet40_vstar_ = struct_mu8_Jet40_vstar.dR_;

                    ///////////////////////  
                    // End 2011 Triggers //
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
                        double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), p4j);
                        if (dr > deltaRCut) hlt15u_ = 2;
                    }
                    if (hlt30u_ == 1) {
                        LorentzVector p4j = p4HLTObject("HLT_Jet30U",0);
                        double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), p4j);
                        if (dr > deltaRCut) hlt30u_ = 2;
                    }
                    if (hlt50u_ == 1) {
                        LorentzVector p4j = p4HLTObject("HLT_Jet50U",0);
                        double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), p4j);
                        if (dr > deltaRCut) hlt50u_ = 2;
                    }
                    if (l16u_ == 1) {
                        LorentzVector p4j = p4HLTObject("HLT_L1Jet6U",0);
                        double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), p4j);
                        if (dr > deltaRCut) l16u_ = 2;
                    }
                    if (l110u_ == 1) {
                        LorentzVector p4j = p4HLTObject("HLT_L1Jet10U",0);
                        double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), p4j);
                        if (dr > deltaRCut) l110u_ = 2;
                    }
        
                    // Now fill the muon trigger flags
                    pair<int, float> pair_mu5  = TriggerMatch( mus_p4().at(iLep), "HLT_Mu5"    );
                    pair<int, float> pair_mu7  = TriggerMatch( mus_p4().at(iLep), "HLT_Mu7"    );
                    pair<int, float> pair_mu9  = TriggerMatch( mus_p4().at(iLep), "HLT_Mu9"    );
                    pair<int, float> pair_mu11 = TriggerMatch( mus_p4().at(iLep), "HLT_Mu11"   );
                    pair<int, float> pair_mu13 = TriggerMatch( mus_p4().at(iLep), "HLT_Mu13_v1");
                    pair<int, float> pair_mu15 = TriggerMatch( mus_p4().at(iLep), "HLT_Mu15_v1");
                    pair<int, float> pair_mu17 = TriggerMatch( mus_p4().at(iLep), "HLT_Mu17_v1");
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

                    ///////////////////////  
                    // End 2010 Triggers //
                    ///////////////////////



                    //////////////
                    // Jets     //
                    //////////////

                    // Calo Jets
                    // Find the highest Pt jet separated by at least dRcut from this lepton and fill the jet Pt
                    ptj1_       = -999.0;
                    ptj1_b2b_   = -999.0;
                    dphij1_b2b_ = -999.0;
                    nj1_        = 0;
                    for (unsigned int iJet = 0; iJet < jets_p4().size(); iJet++) {
                        double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jets_p4().at(iJet) );
                        if( dr > deltaRCut && jets_p4().at(iJet).pt() > 10 ) nj1_++;
                        if ( dr > deltaRCut && jets_p4().at(iJet).pt() > ptj1_ ){
                            ptj1_ = jets_p4().at(iJet).pt();
        
                            // back to back in phi
                            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iLep), jets_p4().at(iJet) ) );
                            if( dphi > deltaPhiCut && jets_p4().at(iJet).pt() > ptj1_b2b_ ){        
                                ptj1_b2b_   = jets_p4().at(iJet).pt();
                                dphij1_b2b_ = dphi;
                            }
                        }
                    }
  
                    // PF Jets
                    // Find the highest Pt pfjet separated by at least dRcut from this lepton and fill the pfjet Pt
                    ptpfj1_       = -999.0;
                    ptpfj1_b2b_   = -999.0;
                    dphipfj1_b2b_ = -999.0;
                    npfj1_        = 0;
                    for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                        double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), pfjets_p4().at(iJet) );
                        if( dr > deltaRCut && pfjets_p4().at(iJet).pt() > 10 ) npfj1_++;
                        if ( dr > deltaRCut && pfjets_p4().at(iJet).pt() > ptpfj1_ ){
                            ptpfj1_ = pfjets_p4().at(iJet).pt();
        
                            // back to back in phi
                            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iLep), pfjets_p4().at(iJet) ) );
                            if( dphi > deltaPhiCut && pfjets_p4().at(iJet).pt() > ptpfj1_b2b_ ){        
                                ptpfj1_b2b_   = pfjets_p4().at(iJet).pt();
                                dphipfj1_b2b_ = dphi;
                            }
                        }
                    }
  
                    // L2L3 PF Jets
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
                        float jet_cor = jetCorrection(jp4, jet_pf_L2L3corrector);
                        LorentzVector jp4cor = jp4 * jet_cor;
                        if (jp4cor.pt() > 15 && pfjets_simpleSecondaryVertexHighEffBJetTag().at(iJet) > 1.74 ) btagpfc_ = true;
                        double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jp4cor );
                        if( dr > deltaRCut && jp4cor.pt() > 10 ) npfcj1_++;
                        if ( dr > deltaRCut && jp4cor.pt() > ptpfcj1_ ){
                            ptpfcj1_ = jp4cor.pt();
    
                            // back to back in phi
                            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iLep), jp4cor ) );
                            if( dphi > deltaPhiCut && jp4cor.pt() > ptpfcj1_b2b_ ){
                                ptpfcj1_b2b_   = jp4cor.pt();
                                dphipfcj1_b2b_ = dphi;
                            }
                        }
                    }

                    // L1FastL2L3 PF Jets
                    // Find the highest Pt PF corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                    ptpfcL1Fj1_       = -999.0;
                    ptpfcL1Fj1_b2b_   = -999.0;
                    dphipfcL1Fj1_b2b_ = -999.0;
                    npfcL1Fj1_        = 0;
                    btagpfcL1F_       = false;
                    for (unsigned int iJet = 0; iJet < pfjets_p4().size(); iJet++) {
                        // JetID
                        if ( !passesPFJetID(iJet)) continue;
                        LorentzVector jp4 = pfjets_p4()[iJet];
                        float jet_cor = jetCorrection(jp4, jet_pf_L1FastL2L3corrector);
                        LorentzVector jp4cor = jp4 * jet_cor;
                        if (jp4cor.pt() > 15 && pfjets_simpleSecondaryVertexHighEffBJetTag().at(iJet) > 1.74 ) btagpfcL1F_ = true;
                        double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jp4cor );
                        if( dr > deltaRCut && jp4cor.pt() > 10 ) npfcL1Fj1_++;
                        if ( dr > deltaRCut && jp4cor.pt() > ptpfcL1Fj1_ ){
                            ptpfcL1Fj1_ = jp4cor.pt();

                            // back to back in phi
                            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iLep), jp4cor ) );
                            if( dphi > deltaPhiCut && jp4cor.pt() > ptpfcL1Fj1_b2b_ ){
                                ptpfcL1Fj1_b2b_   = jp4cor.pt();
                                dphipfcL1Fj1_b2b_ = dphi;
                            }
                        }
                    }


                    // L2L3 JPT Jets
                    // Find the highest Pt JPT corrected jet separated by at least dRcut from this lepton and fill the jet Pt
                    ptjptcj1_       = -999.0;
                    ptjptcj1_b2b_   = -999.0;
                    dphijptcj1_b2b_ = -999.0;
                    njptcj1_        = 0;
                    btagjptc_       = false;
                    for (unsigned int iJet = 0; iJet < jpts_p4().size(); iJet++) {
                        // JetID
                        LorentzVector jp4 = jpts_p4()[iJet];
                        if ( !passesCaloJetID(jp4)) continue;
                        float jet_cor = jetCorrection(jp4, jet_jpt_L2L3corrector);
                        LorentzVector jp4cor = jp4 * jet_cor;
                        if (jp4cor.pt() > 15 && jpts_simpleSecondaryVertexHighEffBJetTag().at(iJet) > 1.74 ) btagjptc_ = true;
                        double dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jp4cor );
                        if( dr > deltaRCut && jp4cor.pt() > 10 ) njptcj1_++;
                        if ( dr > deltaRCut && jp4cor.pt() > ptjptcj1_ ){
                            ptjptcj1_ = jp4cor.pt();

                            // back to back in phi
                            float dphi = fabs( ROOT::Math::VectorUtil::DeltaPhi( mus_p4().at(iLep), jp4cor ) );
                            if( dphi > deltaPhiCut && jp4cor.pt() > ptjptcj1_b2b_ ){
                                ptjptcj1_b2b_   = jp4cor.pt();
                                dphijptcj1_b2b_ = dphi;
                            }
                        }
                    }

                    //////////////
                    // End Jets //
                    //////////////



                    ///////////////////
                    // B Tagging     //
                    ///////////////////

                    // The btag information
                    nbjet_ = this_nbjet;
                    dRbNear_ =  99.;
                    dRbFar_  = -99.;
                    for (int ii=0; ii<nbjet_; ii++) {
                        unsigned int iJet = bindex[ii];
                        float dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jets_p4().at(iJet));
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
                        float jet_cor = jetCorrection(jp4, jet_pf_L2L3corrector);
                        LorentzVector jp4cor = jp4 * jet_cor;
                        float dr = ROOT::Math::VectorUtil::DeltaR( mus_p4().at(iLep), jp4cor);
                        if (dr < dRbpfcNear_) dRbpfcNear_ = dr;
                        if (dr > dRbpfcFar_)   dRbpfcFar_  = dr;
                    }

                    ///////////////////
                    // End B Tagging //
                    ///////////////////


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

