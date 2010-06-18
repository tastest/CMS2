/* Usage:
   root[0] .L ttDilRefS10_looper.C++
   root [1] TFile *_file0 = TFile::Open("ntuple_file.root")
   root [2] TChain *chain = new TChain("Events")
   root [3] chain->Add("ntuple_file.root")
   root [4] ttDilRefS10_looper a 
   root [5] a.ScanChain(chain) // will give the same results
*/
#include <iostream>
#include <vector>
#include <map>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "TChain.h"
#include "TFile.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TROOT.h"

//from slava
#include "TH1F.h"
#include "TH2F.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include <algorithm>
#include <fstream>

#include "CORE/CMS2.h"
#include "ttDilRefS10_looper.h"
#include "CORE/selections.cc"
#include "CORE/muonSelections.cc"
#include "CORE/electronSelections.cc"
#include "CORE/metSelections.cc"
#include "CORE/utilities.cc"
#include "CORE/MT2.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef std::vector<LorentzVector >  VofLorentzVector;



//------------------------------------------------------------


void ttDilRefS10_looper::fill1D(TH1F* h, double v, double w){
  unsigned int nB = h->GetNbinsX();
  double hMin = h->GetXaxis()->GetBinLowEdge(1);
  double hMax = h->GetXaxis()->GetBinUpEdge(nB);
  double bminw =  h->GetXaxis()->GetBinWidth(1);
  double bmaxw =  h->GetXaxis()->GetBinWidth(nB);

  h->Fill(min(max(v,hMin+bminw*0.01),hMax-bmaxw*0.01),w);
}

int ttDilRefS10_looper::ScanChain (std::string fName, std::string prefix, float kFactor, int prescale, 
				   unsigned long long int cutsMask){
  TChain* chain = new TChain("Events");
  chain->Add(fName.c_str()); //chain->GetEntries("mus_p4.pt()");
  return ScanChain(chain, prefix, kFactor, prescale, cutsMask);
}
int ttDilRefS10_looper::ScanChain (TChain* chain, std::string prefix, float kFactor, int prescale, 
				   unsigned long long int cutsMask){
  ProcDSS pds(prefix,kFactor,prescale); 
  ProcDSChain pdsTmp(chain, pds.name);
  pds.add(pdsTmp);
  return ScanChain(pds, cutsMask);
}
int ttDilRefS10_looper::ScanChain (ProcDSS& pds, unsigned long long int cutsMask){
  TDirectory *rootdir = gROOT->GetDirectory("root:");
  if (rootdir)  rootdir->cd();
  else {
    std::cout<<"Cant find root: . Current dir is "<<gDirectory->GetName()<<std::endl;
    rootdir = gROOT->GetDirectory("Rint:");
    if (rootdir){
      std::cout<<"OK, got Rint: "<<std::endl;
      rootdir->cd();
    } else {
      std::cout<<"Cant find Rint: either . Current dir is "<<gDirectory->GetName()<<std::endl;
    }
  }
  
  using namespace std;
  
  std::string prefix(pds.name);
  int prescale(pds.prescale);
  float kFactor(pds.kFactor);
  
  ofstream outf;
  outf.open(Form("eventsCMS2_%d.txt",cutsMask));
  outf << "Starting " << prefix << " bitmask: " << cutsMask << endl;
  std::cout << "Starting " << prefix << " bitmask: " << cutsMask << std::endl;
  
  
  //book Histograms
  bookHistos(prefix);

  compactConfig = "";

  bool genDilCount =  cutsMask & 1;
  bool genDilCountHyp = genDilCount; //MC truth defines the hyps
  if( genDilCount) {
    cout << "GenP dileptons count" << endl;
    compactConfig = compactConfig+ "_genpDil";
  }

  bool genDilTrigBit = (cutsMask>>1) & 1;
  if( genDilTrigBit ) {
    cout << "genDil type-driven trigger bit selection " << endl;
    compactConfig = compactConfig + "_genDilTrigBit";
  }
  
  bool lepCount = ((cutsMask >> 2)&1);
  if (lepCount){
    std::cout<< "The weight for each nJet is nLeptons passing selection"<<std::endl;
    compactConfig = compactConfig + "_lepCnt";
  }

  //the following cuts need a hyp
  bool needAHyp = (cutsMask>>3) != 0;
  bool refS10IDIso =  (cutsMask>>3) & 1;
  if( refS10IDIso ) {
    cout << "TaS lepton ID and isolation (make sure to turn off other ID and iso)" << endl;
    compactConfig = compactConfig + "_refS10IDIso";
  } 
  
  bool osSelection = ((cutsMask >> 4 ) & 1);
  if (osSelection){
    std::cout<<"Require OS"<<std::endl;
    compactConfig = compactConfig + "_OS";
  }
  
  bool zVetoCut = ((cutsMask>>5)&1);
  if (zVetoCut){
    std::cout<<"Apply Z mass veto on same flavor dils, use TTDil08 selections"<<std::endl;
    compactConfig = compactConfig + "_outZ";
  }
  
  bool ge2Jets = ((cutsMask>>6)&1);
  if (ge2Jets){
    std::cout<<"Apply >=2 jets "<<std::endl;
    compactConfig = compactConfig + "_ge2j";
  }

  bool usePF = (cutsMask>>7) & 1; 
  if( usePF ) {
    cout << "Use PF jets (pt>30, eta<2.4) and PFMet" << endl;
    compactConfig = compactConfig + "_pf";
  } else {
    cout << "Use Calo jets (pt>30, eta<2.4) and CaloMet" << endl;
    compactConfig = compactConfig + "_cal";
  }
  // Now we get counted jets
  JetCollectionType jetCType = CaloJetCorr_jct;
  if (usePF) jetCType = PF_jct;

  bool muJetClean = true;

  bool passMET2030 = ((cutsMask>>8)&1);
  if (passMET2030){
    std::cout<<"Apply MET emu met >20, mm,ee met>30"<<std::endl;
    compactConfig = compactConfig + "_met2030";
  }
   
  METCollectionType metCType = PatMETCor_mct;
  if (usePF) metCType = PFMET_mct;
  
  bool ge2BtagsHiEff1p7 =  ((cutsMask>>9)&1);
  if (ge2BtagsHiEff1p7) {
    std::cout<< "Apply two btags with track-counting high-eff disc>1.7"<<std::endl;
    compactConfig = compactConfig + "_2bHE";
  }
  bool ge1BtagsHiEff1p7 =  ((cutsMask>>10)&1);
  if (ge1BtagsHiEff1p7) {
    std::cout<< "Apply two btags with track-counting high-eff disc>1.7"<<std::endl;
    compactConfig = compactConfig + "_1bHE";
  }

  std::cout<<"Compact config string is "<<compactConfig.c_str()<<std::endl;
  
  //Event Loop
  int nAfterPrescale = 0;
  unsigned int nEventsTotal = 0;
  int nAllEvents = 0;
  map<int,int> m_events;


  for(unsigned int iPDS = 0; iPDS<pds.size(); ++iPDS){
    TChain* chain = pds[iPDS].events;
    bool useWeigtFromBranch = pds[iPDS].useWeigtFromBranch;
    float scale1fb = pds[iPDS].scale1fb;

    TObjArray *listOfFiles = chain->GetListOfFiles();
    //    unsigned int nEventsChain=chain->GetEntries();
    
    std::map<EIDiif,bool> evId;
    bool doCheckDuplicateEvents = pds[iPDS].checkDuplicates;
    
    bool printEvents = false;
    
    // file loop
    TIter fileIter(listOfFiles);
    while(TChainElement *currentFile = (TChainElement*)fileIter.Next() ) {
      TFile f(currentFile->GetTitle());
      TTree *tree = (TTree*)f.Get("Events");
      cms2.Init(tree);
      unsigned int nEntries = tree->GetEntries();
      unsigned int nLoop = nEntries;
      
      unsigned int z;
      
      for( z = 0; z < nLoop; z++) {
	nAllEvents++;
	// Progress feedback to the user
	int iz = nAllEvents/10000;
	if (nAllEvents-10000*iz == 0) cout << "Processing event " << nAllEvents+1 
					   << " of sample " << prefix << endl;
	
	nAfterPrescale++;
	cms2.GetEntry(z);
	++nEventsTotal;


      //check if it's a correct genp-event
      std::string prefixStr(prefix);
      if (prefixStr == "ttdil" && getNumberStatus3Leptons() != 2) continue;
      if (prefixStr == "ttotr" && getNumberStatus3Leptons() ==2) continue;

      if (prefixStr == "DYeemm" && genpCountPDGId(11) != 2 && genpCountPDGId(13) != 2) continue;
      if (prefixStr == "DYee" && genpCountPDGId(11) != 2) continue;
      if (prefixStr == "DYmm" && genpCountPDGId(13) != 2) continue;
      if (prefixStr == "DYtautau" && genpCountPDGId(15) != 2) continue;

      int genDilType = getGenDilType();

      //decide weather or not the event passed
      bool eventPassed = false;
      

      std::vector<unsigned int> goodHyps(0);

      if (doCheckDuplicateEvents){
	EIDiif eid(cms2.evt_run(),cms2.evt_event(),cms2.evt_met());
	if (evId[eid]){
	  std::cout<<prefixStr<<" Duplicate evt: "
		   <<" ::: "<<eid.i0<<" "<<eid.i1<<" "<<eid.f0
		   <<std::endl;
	  continue;
	} else {
	  evId[eid] = true;
	}
	if (evId.size()%100000 == 0)std::cout<<evId.size()<<" unique events read"<<std::endl;
      }

      if (genDilCountHyp && genDilType == -1 ) continue;
      //loop over hyps and select good ones
      unsigned int nHyps = cms2.hyp_p4().size();
      if (genDilCountHyp && !needAHyp && nHyps == 0 ) eventPassed = true;
      for(unsigned int hypIdx = 0; hypIdx < nHyps; hypIdx++) {
	
	unsigned int i_lt = cms2.hyp_lt_index()[hypIdx];
	unsigned int i_ll = cms2.hyp_ll_index()[hypIdx];
	
	int id_lt = cms2.hyp_lt_id()[hypIdx];
	int id_ll = cms2.hyp_ll_id()[hypIdx];
	
	float pt_lt = cms2.hyp_lt_p4()[hypIdx].pt();
	float pt_ll = cms2.hyp_ll_p4()[hypIdx].pt();
	
	int hyp_type = cms2.hyp_type()[hypIdx];
	
	if(zVetoCut) {
	  // Z mass veto using hyp_leptons for ee and mumu final states
	  if (((!genDilCountHyp) && (hyp_type == 0 || hyp_type == 3))
	      || (genDilCountHyp && (genDilType ==0 || genDilType == 1))) {
	    if (inZmassWindow(cms2.hyp_p4()[hypIdx].mass())) continue;
	  }    
	}
	
	if (refS10IDIso){
	  if (abs(id_lt)==11 && ! electronSelectionTTbar_refS10(i_lt) ) continue;
	  if (abs(id_ll)==11 && ! electronSelectionTTbar_refS10(i_ll) ) continue;
	  
	  if (abs(id_lt)==13 && ! muonId(i_lt,TTbarRefS10)) continue;
	  if (abs(id_ll)==13 && ! muonId(i_ll,TTbarRefS10)) continue;
	}
	
	if (osSelection){
	  if ( id_lt * id_ll > 0 ) continue;
	}
	
	VofLorentzVector new_hyp_jets_p4_val(jetsForCounting(hypIdx, 1., muJetClean, jetCType, 30, 2.4, 0.4));
	unsigned int nJetsCounted = new_hyp_jets_p4_val.size();
	using ROOT::Math::VectorUtil::DeltaR;
	/*
	for (unsigned int iJ=0; iJ < nJetsCounted; ++iJ){
	  for (unsigned int iMu = 0; iMu < cms2.mus_type().size(); ++iMu){
	    if (!muonId(iMu, TTbarRefS10)) continue;
	    if (DeltaR(new_hyp_jets_p4_val[iJ], cms2.mus_p4()[iMu]) < 0.4){
	      std::cout<<"Found unclean jet wrt muon, jet pt: "<<new_hyp_jets_p4_val[iJ].pt()<<std::endl;
	    }
	  }
	  for (unsigned int iEl = 0; iEl < cms2.els_type().size(); ++iEl){
	    if (! electronSelectionTTbar_refS10(iEl)) continue;
	    if (DeltaR(new_hyp_jets_p4_val[iJ], cms2.els_p4()[iEl])< 0.4){
	      std::cout<<"Found unclean jet wrt muon, jet pt: "<<new_hyp_jets_p4_val[iJ].pt()<<std::endl;
	    }
	  }
	}
	*/

	// ! for TTDil analysis this should be made for the event-qualifying hyp only
	if (passMET2030){
	  if (!genDilCountHyp){
	    // do this "correctly" here?
	    if (passMET2030 && !passMet_OF_SF(hypIdx,metCType,20,30)) continue;
	  } else {
	    if (genDilType == 0 || genDilType == 1){
	      if (!usePF && cms2.met_pat_metCor() < 30 ) continue;
	      if (usePF && cms2.evt_pfmet() < 30 ) continue;
	    } else {
	      if (!usePF && cms2.met_pat_metCor() < 20 ) continue;
	      if (usePF && cms2.evt_pfmet() < 20 ) continue;
	    }
	  }
	}

	if (ge2Jets && nJetsCounted < 2 ) continue;
	
	if (ge2BtagsHiEff1p7 || ge1BtagsHiEff1p7){
	  int nBtags = 0;
	  if (usePF){
	    //we don't have it properly :(
	    //so, I'm matching to calojets or trk jets
	    for (unsigned int iJ = 0; iJ < nJetsCounted; ++iJ){
	      //lets find the jet index
	      int iJet = -1;
	      // calo
	      for (unsigned int jJ = 0; jJ <cms2.jets_p4().size(); ++jJ ){
		if (DeltaR(new_hyp_jets_p4_val[iJ], cms2.jets_p4()[jJ]) < 0.1){
		  iJet = jJ; break;
		}
	      }
	      if (iJet != -1 && cms2.jets_trackCountingHighEffBJetTag()[iJet] > 1.7) nBtags++;
	      if (iJet == -1){//try trk
		for (unsigned int jJ = 0; jJ <cms2.jets_p4().size(); ++jJ ){
		  if (DeltaR(new_hyp_jets_p4_val[iJ], cms2.trkjets_p4()[jJ]) < 0.1){
		    iJet = jJ; break;
		  }
		}
		if (iJet != -1 && cms2.trkjets_trackCountingHighEffBJetTag()[iJet] > 1.7) nBtags++;
	      }//  if (iJet == -1)
	    }
	  }  else {
	    //calojets
	    for (unsigned int iJ = 0; iJ < nJetsCounted; ++iJ){
	      //lets find the jet index
	      int iJet = -1;
	      for (unsigned int jJ = 0; jJ <cms2.jets_p4().size(); ++jJ ){
		if (DeltaR(new_hyp_jets_p4_val[iJ], cms2.jets_p4()[jJ]) < 0.1){
		  iJet = jJ; break;
		}
	      }
	      if (iJet == -1) continue;
	      if (cms2.jets_trackCountingHighEffBJetTag()[iJet] > 1.7) nBtags++;
	    }
	  }

	  //check now
	  if (ge2BtagsHiEff1p7 && nBtags < 2 ) continue;
	  if (ge1BtagsHiEff1p7 && nBtags == 0 ) continue;
	}

	goodHyps.push_back(hypIdx);
	//done with cuts on hyps		   
	
	eventPassed = true;
	
      }
      if (genDilTrigBit){
	bool haveBit = false;
	if (genDilType == 0 || genDilType == 2){
	  if (cms2.passHLT8E29Trigger("HLT_Ele15_LW_L1R")) haveBit = true;
	}
	if (genDilType == 1 || genDilType == 2){
	  if (cms2.passHLT8E29Trigger("HLT_Mu9")) haveBit = true;
	}
	if (! haveBit) eventPassed = false;
	//	  else std::cout<<"Have bit for "<<genDilType<<std::endl;
      }
      
      if (! eventPassed) continue;
      
      //=============================================================================================

      //now fill the histograms

      //=============================================================================================
      if (genDilCountHyp && genDilType != -1){
	// correspond to 0, 1, ge.2, 2-3, ge.3, ge.4
        bool fillJetSel[6] = {0, 0, 1, 1, 0, 0};
	float weightGL = useWeigtFromBranch? kFactor*cms2.evt_scale1fb()*0.01 : scale1fb*0.01; // /100; //10pb^-1
	if (lepCount ){
	  int nLeps = 0;
	  
	  int nMus = cms2.mus_type().size();
	  for (int iMu=0; iMu<nMus; ++iMu){
	    if (muonId(iMu, TTbarRefS10)){
	      nLeps++;
	    }
	  }
	  int nEls = cms2.els_type().size();
	  for (int iEl=0; iEl<nEls; ++iEl){
	    if (electronSelectionTTbar_refS10(iEl)){
	      nLeps++;
	    }
	  }
	  float nLepsF = nLeps;
	  weightGL *= nLepsF; // step 2 only (this is a terrible way to do this)
	}
	
	fill1D(hnJet[(genDilType&3)], 2, weightGL);
	fill1D(hnJet[3], 2, weightGL);
      }
      
      }//event loop
    }//file loop
  }//pds loop
  //  if ( nEventsChain != nEventsTotal ) {
  //    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  //  }

  std::cout<<"Done with "<<prefix<<std::endl;
  outf.close();
  std::cout<<"Closed outf "<<std::endl;
  rootdir = gROOT->GetDirectory("root:"); 
  if (rootdir) rootdir->cd(); 
  else{
    std::cout<<"Cant find root: . Current dir is "<<gDirectory->GetName()<<std::endl;
    rootdir = gROOT->GetDirectory("Rint:");
    if (rootdir){
      std::cout<<"OK, got Rint: "<<std::endl;
      rootdir->cd();
    } else {
      std::cout<<"Cant find Rint: either . Current dir is "<<gDirectory->GetName()<<std::endl;
    }
  }

  return 0;
}




void ttDilRefS10_looper::bookHistos(std::string& prefixS) {

  //  Book histograms...
  //  Naming Convention:
  //  Prefix comes from the sample and it is passed to the scanning function
  //  Suffix is "ee" "em" "em" "all" which depends on the final state
  //  For example: histogram named tt_hnJet_ee would be the Njet distribution
  //  for the ee final state in the ttbar sample.
  
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  //  const char* prefix = prefixS.c_str();
  std::string jetbins[5] = {"0", "1", "2", "3", "#geq 4"};
  TDirectory *rootdir = gROOT->GetDirectory("root:");
  std::cout<<"Current dir is "<<gDirectory->GetName()<<std::endl;
  if (rootdir == 0){
    std::cout<<"Head directory root: not found. Try Rint: ..."<<std::endl;
    rootdir = gROOT->GetDirectory("Rint:");
    if (rootdir){
      std::cout<<"OK: Got Rint:"<<std::endl;
    } else {
      std::cout<<"ERROR: no root: or Rint: found. Histograms will likely be lost"<<std::endl;
    }
  }

  for (int i=0; i<4; i++) {
    for (int j=0; j<6; j++) { //6 jet regions
      std::string suffixall[4] = { "ee",  "mm", "em",  "all"};
      //      char *suffix[4];
  
      if (j == 0){
	hnJet[i] = new TH1F(Form("%s_hnJet_%s",prefixS.c_str(),suffixall[i].c_str()),
			    Form("%s_nJet_%s",prefixS.c_str(),suffixall[i].c_str()),
			    5,0.,5.);	
	hnJet[i]->SetDirectory(rootdir);
	hnJet[i]->GetXaxis()->SetTitle("nJets");

	for(int k = 0; k<5; k++) {
	  hnJet[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k].c_str());
	  hnJet[i]->GetXaxis()->SetLabelSize(0.07);
	}
      }
    
      if (j==0){
	hnJet[i]->Sumw2();
      }
    }
  }//channel loop
  //  gDirectory->ls();
  //  tmpFile_->ls();
  //  tmpFile_->Print();
}//CMS2::bookHistos()
