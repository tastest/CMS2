/*
  Common EventList creation tools.
  Usage:
     TFile *_file0 = TFile::Open("/home/users/dmytro/ntuples/cms2-V00-04-00/merge_WW.root")
     .L EventListCreator.C+
     EventListCreator creator
     ww=(TTree*)_file0->Get("Events")
     list = creator.getAllUniqueEvents(ww)
     ww->SetEventList(list)
     Events->Draw("hyp_type");

*/
#include <iostream>
#include <vector>
#include <set>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TEventList.h"
#include "EventListCreator.h"
#include "TDatabasePDG.h"
#include "TProfile.h"

int EventListCreator::dumpDocLines ( TTree* tree, int first, int last) {

  Init(tree);


	// figure out how many events to dump
	int nEntries = tree->GetEntries();
        int howMany = nEntries;
	if (last > 0) howMany = min(nEntries, last+1);
	

	// Initialize particle database
	TDatabasePDG *pdg = new TDatabasePDG(); 	

	 //Event Loop
	for( int i = first; i < howMany; i++) {
	  GetEntry(i);

	  cout << "   " << endl;
	  cout << "---------" << " Entry = " << i << " ------------------------------" << endl;
          cout << "ID             Pt      Phi         Eta        Mass         MotherID" << endl;
	  for (unsigned int j=0; j<genps_id().size(); j++) {
	    cout << setw(9) << left << pdg->GetParticle(genps_id().at(j))->GetName() << " "
		 << setw(7) << right << setprecision(4) << genps_p4().at(j).pt() << "  "
                 << setw(7) << right << setprecision(4) << genps_p4().at(j).phi() << "  "
                 << setw(10) << right << setprecision(4) << genps_p4().at(j).eta() << "  "
                 << setw(10) << right << setprecision(4) << genps_p4().at(j).mass() << "  " 
		 << pdg->GetParticle(genps_id_mother().at(j))->GetName() << endl;
	  }
	  /*
	  // get decay tree by plain combinatorics starting from the bottom.
	  std::vector<std::set<int > > children(genps_id().size());
	  std::fill( children.begin(), children.end(), std::set<int>() );
	  std::vector<int> parents(genps_id().size());
	  std::fill( parents.begin(), parents.end(), -1 );
	  // 2Body
	  int iChild1 = genps_id().size()-1;
	  while (iChild1 > 1){
	    int iChild2 = iChild1-1;
	    while (iChild2 > 0){
	      LorentzVector p4 = genps_p4()[iChild1]+genps_p4()[iChild2];
	      //cout << iChild1 << "\t " << iChild2 << "\t " << p4.mass() << "\t " << 
	      // "\t " << p4.P() << "\t " << p4 << endl;
	      int iParent = iChild2-1;
	      while ( iParent >= 0 ){
		LorentzVector pp4 = p4 - genps_p4()[iParent];
		if ( pp4.P() < 1e-3 ) {
		  //matched
		  parents[iChild1] = iParent;
		  parents[iChild2] = iParent;
		  children[iParent].insert( iChild1 );
		  children[iParent].insert( iChild2 );
		  cout << pdg->GetParticle(genps_id().at(iParent))->GetName() << " -> " << 
		    pdg->GetParticle(genps_id().at(iChild1))->GetName() << " + " <<
		    pdg->GetParticle(genps_id().at(iChild2))->GetName() << endl;
		  break;
		}
		--iParent;
	      }
	      if ( parents[iChild2] > -1 ) break;
	      --iChild2;
	    }
	    --iChild1;
	  }
	  for (unsigned int i=0; i<genps_id().size(); ++i) {
	    if ( children[i]->empty() ){
	      cout << pdg->GetParticle(genps_id().at(i))->GetName() << " -> " << endl;
	      for (unsigned int j=0; j<genps_id().size(); j++) {
		if ( parent[j] == -1 && genps_id_mother()[j] == genps_id()[i] )
		
	  */
 	}

	// Clean up ofter ourselves
	delete pdg;
	return 0;
}


TEventList* EventListCreator::getAllUniqueEvents( TTree* tree, 
				      const char* name /*= "AllUniqueEvents"*/, 
				      const char* title /*= "All events without duplicates"*/ )
{
  Init(tree);
  TEventList* list = new TEventList(name,title);
  std::set<ULong64_t> known_events;
  unsigned int nEvents = tree->GetEntries();
  for( unsigned int event = 0; event < nEvents; ++event) {
    if ( tree->GetEventList() && !tree->GetEventList()->Contains(event) ) continue;
    GetEntry(event);
    ULong64_t key = int((evt_met()-int(evt_met()))*10000) + 
      ULong64_t(evt_run()%100000)*10000 + 
      ULong64_t(evt_event()%100000)*1000000000;
    if ( known_events.insert(key).second ) list->Enter( event );
  }
  return list;
}


void EventListCreator::test_wjets_bkg( TTree* tree )
{
  TH1F* h1 = new TH1F("h1","Parton PDG id",30,0,30);
  TH1F* h2 = new TH1F("h2","Parton energy fraction",20,0.,1.);
  TH1F* h3 = new TH1F("h3","Parton Pt",50,0,100);
  TH1F* h4 = new TH1F("h4","Parton Eta",20,0,5.);
  Init(tree);
  TEventList* list = getAllUniqueEvents(tree);
  unsigned int nEvents = tree->GetEntries();
  int index = 0;
  double nIsoSidebandFOIdBarrel = 0;
  double nIsoSidebandFOIdEndcap = 0;
  double nIsoSignalFOIdBarrel = 0;
  double nIsoSignalFOIdEndcap = 0;
  double nIsoSidebandFOD0Barrel = 0;
  double nIsoSidebandFOD0Endcap = 0;
  double nIsoSignalFOD0Barrel = 0;
  double nIsoSignalFOD0Endcap = 0;
  double nIsoSidebandBarrel = 0;
  double nIsoSidebandEndcap = 0;
  double nIsoSignalBarrel = 0;
  double nIsoSignalEndcap = 0;
  for( unsigned int event = 0; event < nEvents; ++event) {
    if ( !list->Contains(event) ) continue;
    GetEntry(event);
    for ( unsigned int hypi = 0; hypi < hyp_p4().size(); ++hypi ){
      bool nominalSelection 
	= (((ww_ltgoodmuiso()[hypi]==1 && ww_llgoodel()[hypi]==1 && hyp_type()[hypi]==1)||
	    (ww_ltgoodel()[hypi]==1 && ww_llgoodmuiso()[hypi]==1 && hyp_type()[hypi]==2)) &&
	   ww_pass4met()[hypi]==1 && ww_passzveto()[hypi]==1 && ww_passaddzveto()[hypi]==1 && 
	   hyp_ll_p4()[hypi].Pt()>20 && ww_oppsign()[hypi]==1 && hyp_njets()[hypi]==0);
      bool relaxedIdSelection
	= (((ww_ltgoodmuiso()[hypi]==1 && (abs(els_d0()[hyp_ll_index()[hypi]]) <= 0.025) && hyp_type()[hypi]==1)||
	    ((abs(els_d0()[hyp_lt_index()[hypi]]) <= 0.025) && ww_llgoodmuiso()[hypi]==1 && hyp_type()[hypi]==2)) &&
	   ww_pass4met()[hypi]==1 && ww_passzveto()[hypi]==1 && ww_passaddzveto()[hypi]==1 && 
	   hyp_ll_p4()[hypi].Pt()>20 && ww_oppsign()[hypi]==1 && hyp_njets()[hypi]==0);

      bool relaxedD0Selection 
	= (((ww_ltgoodmuiso()[hypi]==1 && els_tightId()[hyp_lt_index()[hypi]]==1 && hyp_type()[hypi]==1)||
	    (els_tightId()[hyp_ll_index()[hypi]]==1 && ww_llgoodmuiso()[hypi]==1 && hyp_type()[hypi]==2)) &&
	   ww_pass4met()[hypi]==1 && ww_passzveto()[hypi]==1 && ww_passaddzveto()[hypi]==1 && 
	   hyp_ll_p4()[hypi].Pt()>20 && ww_oppsign()[hypi]==1 && hyp_njets()[hypi]==0);

      // isolation
      double elIso = 0;
      if (  hyp_type()[hypi]==1 )
	elIso = hyp_ll_p4()[hypi].pt()/(hyp_ll_p4()[hypi].pt()+
					els_tkIso()[hyp_ll_index()[hypi]]+
					els_ecalJuraIso()[hyp_ll_index()[hypi]]+
					els_hcalConeIso()[hyp_ll_index()[hypi]]);
      else if ( hyp_type()[hypi]==2 )
	elIso = hyp_lt_p4()[hypi].pt()/(hyp_lt_p4()[hypi].pt()+
					els_tkIso()[hyp_lt_index()[hypi]]+
					els_ecalJuraIso()[hyp_lt_index()[hypi]]+
					els_hcalConeIso()[hyp_lt_index()[hypi]]);
      if ( ! nominalSelection && ! relaxedIdSelection && ! relaxedD0Selection ) continue;
      int elIndex = hyp_type()[hypi]==1 ? hyp_ll_index()[hypi] : hyp_lt_index()[hypi];

      if ( relaxedIdSelection ){
	if ( elIso > 0.7 && elIso < 0.9 ){
	  if ( abs(els_p4()[elIndex].eta()) < 1.5 )
	    nIsoSidebandFOIdBarrel += evt_scale1fb();
	  else
	    nIsoSidebandFOIdEndcap += evt_scale1fb();
	} else {
	  if ( elIso > 0.9 ){
	    if ( abs(els_p4()[elIndex].eta()) < 1.5 )
	    nIsoSignalFOIdBarrel += evt_scale1fb();
	  else
	    nIsoSignalFOIdEndcap += evt_scale1fb();
	  }
	}
      }

      if ( relaxedD0Selection ){
	if ( elIso > 0.7 && elIso < 0.9 ){
	  if ( abs(els_p4()[elIndex].eta()) < 1.5 )
	    nIsoSidebandFOD0Barrel += evt_scale1fb();
	  else
	    nIsoSidebandFOD0Endcap += evt_scale1fb();
	} else {
	  if ( elIso > 0.9 ){
	    if ( abs(els_p4()[elIndex].eta()) < 1.5 )
	    nIsoSignalFOD0Barrel += evt_scale1fb();
	  else
	    nIsoSignalFOD0Endcap += evt_scale1fb();
	  }
	}
      }
      
      if ( nominalSelection ){
	if ( elIso > 0.7 && elIso < 0.9 ){
	  if ( abs(els_p4()[elIndex].eta()) < 1.5 )
	    nIsoSidebandBarrel += evt_scale1fb();
	  else
	    nIsoSidebandEndcap += evt_scale1fb();
	} else {
	  if ( elIso > 0.9 ){
	    if ( abs(els_p4()[elIndex].eta()) < 1.5 )
	    nIsoSignalBarrel += evt_scale1fb();
	  else
	    nIsoSignalEndcap += evt_scale1fb();
	  }
	}
      }
      if ( elIso < 0.9 ) continue;

      printf("%d \t %0.2f\t % d\n",index, elIso, els_mc_id()[elIndex]);
      ++index;
      // loop over gen particles
      double dR = 9999;
      int iMatch = -1;
      for ( unsigned int igen = 0; igen < genps_p4().size(); ++igen ){
	double d = deltaR( els_p4()[elIndex].eta(), els_p4()[elIndex].phi(),
			   genps_p4()[igen].eta(), genps_p4()[igen].phi() );
	if ( d < dR ){
	  dR = d;
	  iMatch = igen;
	}
      }
      if ( iMatch != -1 ){
	printf("\tParton id: % d\t dR: %0.2f\t Energy fraction: %0.2f\n", 
	       genps_id()[iMatch], 
	       dR,
	       els_p4()[elIndex].E()/genps_p4()[iMatch].E() );
	if (dR<0.5){
	  h1->Fill(abs(genps_id()[iMatch])-0.001);
	  h2->Fill(els_p4()[elIndex].E()/genps_p4()[iMatch].E());
	  h3->Fill(genps_p4()[iMatch].pt());
	  h4->Fill(abs(genps_p4()[iMatch].eta()));
	}
      }
    }
  }
  
  printf("\nnIsoSidebandFOIdBarrel: \t%0.2f\n", nIsoSidebandFOIdBarrel);
  printf("nIsoSidebandFOIdEndcap: \t%0.2f\n", nIsoSidebandFOIdEndcap);
  printf("nIsoSignalFOIdBarrel: \t%0.2f\n", nIsoSignalFOIdBarrel);
  printf("nIsoSignalFOIdEndcap: \t%0.2f\n", nIsoSignalFOIdEndcap);
  printf("nIsoSidebandBarrel: \t%0.2f\n", nIsoSidebandBarrel);
  printf("nIsoSidebandEndcap: \t%0.2f\n", nIsoSidebandEndcap);
  printf("nIsoSignalBarrel: \t%0.2f\n", nIsoSignalBarrel);
  printf("nIsoSignalEndcap: \t%0.2f\n", nIsoSignalEndcap);
  double fakeRateIdBarrel = nIsoSidebandBarrel/nIsoSidebandFOIdBarrel;
  printf("fake rate (ID) barrel(nIsoSidebandBarrel/nIsoSidebandFOIdBarrel): \t%0.2f\n", fakeRateIdBarrel);
  double fakeRateIdEndcap = nIsoSidebandEndcap/nIsoSidebandFOIdEndcap;
  printf("fake rate (ID) endcap(nIsoSidebandEndcap/nIsoSidebandFOIdEndcap): \t%0.2f\n", fakeRateIdEndcap);
  printf("background estimate: \t%0.2f\n", fakeRateIdBarrel*nIsoSignalFOIdBarrel+fakeRateIdEndcap*nIsoSignalFOIdEndcap);
  printf("true value: \t%0.2f\n", nIsoSignalBarrel+nIsoSignalEndcap);

  printf("\nnIsoSidebandFOD0Barrel: \t%0.2f\n", nIsoSidebandFOD0Barrel);
  printf("nIsoSidebandFOD0Endcap: \t%0.2f\n", nIsoSidebandFOD0Endcap);
  printf("nIsoSignalFOD0Barrel: \t%0.2f\n", nIsoSignalFOD0Barrel);
  printf("nIsoSignalFOD0Endcap: \t%0.2f\n", nIsoSignalFOD0Endcap);
  printf("nIsoSidebandBarrel: \t%0.2f\n", nIsoSidebandBarrel);
  printf("nIsoSidebandEndcap: \t%0.2f\n", nIsoSidebandEndcap);
  printf("nIsoSignalBarrel: \t%0.2f\n", nIsoSignalBarrel);
  printf("nIsoSignalEndcap: \t%0.2f\n", nIsoSignalEndcap);
  double fakeRateD0Barrel = nIsoSidebandBarrel/nIsoSidebandFOD0Barrel;
  printf("fake rate (d0) barrel(nIsoSidebandBarrel/nIsoSidebandFOD0Barrel): \t%0.2f\n", fakeRateD0Barrel);
  double fakeRateD0Endcap = nIsoSidebandEndcap/nIsoSidebandFOD0Endcap;
  printf("fake rate (d0) endcap(nIsoSidebandEndcap/nIsoSidebandFOD0Endcap): \t%0.2f\n", fakeRateD0Endcap);
  printf("background estimate: \t%0.2f\n", fakeRateD0Barrel*nIsoSignalFOD0Barrel+fakeRateD0Endcap*nIsoSignalFOD0Endcap);



  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);
  c1->cd(1);
  h1->Draw();
  c1->cd(2);
  h2->Draw();
  c1->cd(3);
  h3->Draw();
  c1->cd(4);
  h4->Draw();
}

void EventListCreator::test_sep08_cuts( TTree* tree, double kFactor )
{
  Init(tree);
  TEventList* list = getAllUniqueEvents(tree);
  unsigned int nEvents = tree->GetEntries();
  double nFebSelection(0);
  double nNewNotTagged(0);
  // double nTrkNew(0);
  double nFinalSoftTag(0);
  double nFinal(0);
  double nFinalNoTrkJetsTag(0);
  double nFinalNoTrkJets(0);
  double nBeforeTagSmart(0);
  TH1F* hSoftMuon = new TH1F("hSoftMuon","Soft muon pt",20,0,20);
  TH1F* hTrkJetTrackPt = new TH1F("hTrkJetTrackPt","TrkJet constituent track Pt",20,0,20);
  TH1F* hTrkJetTrackPtMax = new TH1F("hTrkJetTrackPtMax","Maximum of TrkJet constituent track Pt",20,0,20);
  
  // Let's look at trkJetEt/genJetEt as a function of trkJetEt for a few cases:
  // 1) soft
  // 2) medium
  // 3) hard
  TProfile* hTrkJetResponseSoft = new TProfile("hTrkJetResponseSoft","TrkJetEt/GenJetEt for jets made of soft tracks",10,0,100);
  TProfile* hTrkJetResponseMediumSoft = new TProfile("hTrkJetResponseMediumSoft","TrkJetEt/GenJetEt for jets made of medium soft tracks",10,0,100);
  TProfile* hTrkJetResponseMediumHard = new TProfile("hTrkJetResponseMediumHard","TrkJetEt/GenJetEt for jets made of medium hard tracks",10,0,100);
  TProfile* hTrkJetResponseHard = new TProfile("hTrkJetResponseHard","TrkJetEt/GenJetEt for jets made of hard tracks",10,0,100);
  TProfile* hTrkJetResponseNotTagged = new TProfile("hTrkJetResponseNotTagged","TrkJetEt/GenJetEt for jets without muons",10,0,100);
  TProfile* hTrkJetResponseTagged = new TProfile("hTrkJetResponseTagged","TrkJetEt/GenJetEt for jets with muons",10,0,100);

  for( unsigned int event = 0; event < nEvents; ++event) {
    if ( !list->Contains(event) ) continue;
    GetEntry(event);
    for ( unsigned int hypi = 0; hypi < hyp_p4().size(); ++hypi ){
      bool nominalSelection 
	= (((ww_ltgoodmuiso()[hypi]==1 && ww_llgoodel()[hypi]==1 && hyp_type()[hypi]==1)||
	    (ww_ltgoodel()[hypi]==1 && ww_llgoodmuiso()[hypi]==1 && hyp_type()[hypi]==2)) &&
	   ww_pass4met()[hypi]==1 && ww_passzveto()[hypi]==1 && ww_passaddzveto()[hypi]==1 && 
	   hyp_ll_p4()[hypi].Pt()>20 && ww_oppsign()[hypi]==1 );
      if ( ! nominalSelection ) continue;
      // isolation
      double elTrkIso = 0;
      if (  hyp_type()[hypi]==1 )
	elTrkIso = hyp_ll_p4()[hypi].pt()/(hyp_ll_p4()[hypi].pt()+els_tkIso()[hyp_ll_index()[hypi]]);
      else if ( hyp_type()[hypi]==2 )
	elTrkIso = hyp_lt_p4()[hypi].pt()/(hyp_lt_p4()[hypi].pt()+els_tkIso()[hyp_lt_index()[hypi]]);
      
      if ( elTrkIso > 0.92 && hyp_njets()[hypi]==0) nFebSelection += evt_scale1fb();
      double elIso = 0;
      if (  hyp_type()[hypi]==1 )
	elIso = hyp_ll_p4()[hypi].pt()/(hyp_ll_p4()[hypi].pt()+
					els_tkIso()[hyp_ll_index()[hypi]]+
					els_ecalJuraIso()[hyp_ll_index()[hypi]]+
					els_hcalConeIso()[hyp_ll_index()[hypi]]);
      else if ( hyp_type()[hypi]==2 )
	elIso = hyp_lt_p4()[hypi].pt()/(hyp_lt_p4()[hypi].pt()+
					els_tkIso()[hyp_lt_index()[hypi]]+
					els_ecalJuraIso()[hyp_lt_index()[hypi]]+
					els_hcalConeIso()[hyp_lt_index()[hypi]]);
      if ( elIso < 0.9 ) continue;
      if ( hyp_njets()[hypi]>0 ) continue;

      // int elIndex = hyp_type()[hypi]==1 ? hyp_ll_index()[hypi] : hyp_lt_index()[hypi];
      /*
      std::vector<std::vector<std::string> > map;
      for ( unsigned int i = 0; i < 10; ++i )
	map.push_back( std::vector<std::string>(10) );
      printf("---------------------------------------------------------------------\n");
      char buf[1024];
      sprintf(buf, "lt \tpt: %0.2f, eta: %0.2f, phi: %0.2f, id:%d\n", hyp_lt_p4()[hypi].pt(),
	      hyp_lt_p4()[hypi].eta(), hyp_lt_p4()[hypi].phi(), hyp_lt_id()[hypi] );
      addEntryToMap(map, hyp_lt_p4()[hypi].eta(), hyp_lt_p4()[hypi].phi(), buf);
      sprintf(buf, "ll \tpt: %0.2f, eta: %0.2f, phi: %0.2f, id:%d\n", hyp_ll_p4()[hypi].pt(),
	      hyp_ll_p4()[hypi].eta(), hyp_ll_p4()[hypi].phi(), hyp_ll_id()[hypi]);
      addEntryToMap(map, hyp_ll_p4()[hypi].eta(), hyp_ll_p4()[hypi].phi(), buf);
      // print Calo Jets
      for ( unsigned int ijet = 0; ijet < jets_p4().size(); ++ijet ){
	sprintf(buf, "calo jet \tpt: %0.2f, eta: %0.2f, phi: %0.2f\n", jets_p4()[ijet].pt(),
		jets_p4()[ijet].eta(), jets_p4()[ijet].phi() );
	addEntryToMap(map, jets_p4()[ijet].eta(), jets_p4()[ijet].phi(), buf);
      }
      for ( unsigned int ijet = 0; ijet < hyp_jets_p4()[hypi].size(); ++ijet ){
	sprintf(buf, "hyp jet \tpt: %0.2f, eta: %0.2f, phi: %0.2f\n", hyp_jets_p4()[hypi][ijet].pt(),
		hyp_jets_p4()[hypi][ijet].eta(), hyp_jets_p4()[hypi][ijet].phi() );
	addEntryToMap(map, hyp_jets_p4()[hypi][ijet].eta(), hyp_jets_p4()[hypi][ijet].phi(), buf);
      }
      */
      // loop over trkJets
      unsigned int nTrkJets(0);
      unsigned int nNotTaggedTrkJets(0);
      for ( unsigned int ijet = 0; ijet < trkjets_p4().size(); ++ijet ){
	// sprintf(buf, "trk jet \tpt: %0.2f, eta: %0.2f, phi: %0.2f\n", trkjets_p4()[ijet].pt(),
	//       trkjets_p4()[ijet].eta(), trkjets_p4()[ijet].phi() );
	// addEntryToMap( map, trkjets_p4()[ijet].eta(), trkjets_p4()[ijet].phi(), buf );
	if ( trkjets_p4()[ijet].pt() < 15 ) continue;
	double dlt = deltaR( hyp_lt_p4()[hypi].eta(), hyp_lt_p4()[hypi].phi(),
			    trkjets_p4()[ijet].eta(), trkjets_p4()[ijet].phi() );
	if ( dlt < 0.4 ) continue;
	double dll = deltaR( hyp_ll_p4()[hypi].eta(), hyp_ll_p4()[hypi].phi(),
			     trkjets_p4()[ijet].eta(), trkjets_p4()[ijet].phi() );
	if ( dll < 0.4 ) continue;
	// addEntryToMap( map, trkjets_p4()[ijet].eta(), trkjets_p4()[ijet].phi(), "\tpassed\n" );
	++nTrkJets;
	bool muonTag = false;
	bool smartMuonTag = false;
	for ( unsigned int imuon = 0; imuon < mus_p4().size(); ++imuon ){
	  double d = deltaR( mus_p4()[imuon].eta(), mus_p4()[imuon].phi(),
			     trkjets_p4()[ijet].eta(), trkjets_p4()[ijet].phi() );
	  if ( d > 0.4 ) continue;
	  muonTag = true;
	  if ( trkjets_p4()[ijet].pt() > 15 + mus_p4()[imuon].pt() ) continue;
	  smartMuonTag = true;
	}
	if ( !smartMuonTag ) ++nNotTaggedTrkJets;

	// get track pt
	double maxPt = 0;
	for ( unsigned int itrk = 0; itrk < trks_trk_p4().size(); ++itrk ){
	  if ( trks_validHits()[itrk] < 7 ) continue;
	  if ( fabs(trks_d0()[itrk]) > 0.25 ) continue;
	  double d = deltaR( trks_trk_p4()[itrk].eta(), trks_trk_p4()[itrk].phi(),
			     trkjets_p4()[ijet].eta(), trkjets_p4()[ijet].phi() );
	  if ( d > 0.5 ) continue;
	  if ( trks_trk_p4()[itrk].pt() > maxPt ) maxPt = trks_trk_p4()[itrk].pt();
	  if (!muonTag) hTrkJetTrackPt->Fill( trks_trk_p4()[itrk].pt() );
	}
	if (!muonTag) hTrkJetTrackPtMax->Fill(maxPt);

	// find gen jets
	int iMatchedGenJet = -1;
	float matchedGenJetPt = 0;
	for ( unsigned int igenjet = 0; igenjet < genjets_p4().size(); ++igenjet ){
	  double d = deltaR( genjets_p4()[igenjet].eta(), genjets_p4()[igenjet].phi(),
			     trkjets_p4()[ijet].eta(), trkjets_p4()[ijet].phi() );
	  if ( d > 0.5 ) continue;
	  if ( genjets_p4()[igenjet].pt() > matchedGenJetPt ) {
	    matchedGenJetPt = genjets_p4()[igenjet].pt();
	    iMatchedGenJet = igenjet;
	  }
	}
	if ( iMatchedGenJet == -1 ) continue;
	if (maxPt/genjets_p4()[iMatchedGenJet].pt()<0.1) 
	  hTrkJetResponseSoft->Fill(trkjets_p4()[ijet].pt(),trkjets_p4()[ijet].pt()/genjets_p4()[iMatchedGenJet].pt());
	if (maxPt/genjets_p4()[iMatchedGenJet].pt()<0.2 && maxPt/genjets_p4()[iMatchedGenJet].pt()>0.1) 
	  hTrkJetResponseMediumSoft->Fill(trkjets_p4()[ijet].pt(),trkjets_p4()[ijet].pt()/genjets_p4()[iMatchedGenJet].pt());
	if (maxPt/genjets_p4()[iMatchedGenJet].pt()<0.3 && maxPt/genjets_p4()[iMatchedGenJet].pt()>0.2) 
	  hTrkJetResponseMediumHard->Fill(trkjets_p4()[ijet].pt(),trkjets_p4()[ijet].pt()/genjets_p4()[iMatchedGenJet].pt());
	if ( maxPt/genjets_p4()[iMatchedGenJet].pt()>0.3 )
	  hTrkJetResponseHard->Fill(trkjets_p4()[ijet].pt(),trkjets_p4()[ijet].pt()/genjets_p4()[iMatchedGenJet].pt());
	if ( muonTag )
	  hTrkJetResponseTagged->Fill(trkjets_p4()[ijet].pt(),trkjets_p4()[ijet].pt()/genjets_p4()[iMatchedGenJet].pt());
	else
	  hTrkJetResponseNotTagged->Fill(trkjets_p4()[ijet].pt(),trkjets_p4()[ijet].pt()/genjets_p4()[iMatchedGenJet].pt());
      }
    
    if ( hyp_njets()[hypi]>0 ) continue;

      /*
      unsigned int nAllTrkJets(0);
      for ( unsigned int ijet = 0; ijet < alltrkjets_p4().size(); ++ijet ){
	sprintf(buf, "all trk jet \tpt: %0.2f, eta: %0.2f, phi: %0.2f\n", alltrkjets_p4()[ijet].pt(),
	       alltrkjets_p4()[ijet].eta(), alltrkjets_p4()[ijet].phi() );
	addEntryToMap( map, alltrkjets_p4()[ijet].eta(), alltrkjets_p4()[ijet].phi(), buf );
	if ( alltrkjets_p4()[ijet].pt() < 15 ) continue;
	double dlt = deltaR( hyp_lt_p4()[hypi].eta(), hyp_lt_p4()[hypi].phi(),
			    alltrkjets_p4()[ijet].eta(), alltrkjets_p4()[ijet].phi() );
	if ( dlt < 0.5 ) continue;
	double dll = deltaR( hyp_ll_p4()[hypi].eta(), hyp_ll_p4()[hypi].phi(),
			     alltrkjets_p4()[ijet].eta(), alltrkjets_p4()[ijet].phi() );
	if ( dll < 0.5 ) continue;
	addEntryToMap( map, alltrkjets_p4()[ijet].eta(), alltrkjets_p4()[ijet].phi(), "\tpassed\n" );
	++nAllTrkJets;
      }
      */
      // if (
      // if ( nAllTrkJets == 0) nTrkNew += evt_scale1fb();
      if ( nTrkJets == 0 ) nNewNotTagged += evt_scale1fb();
      if ( nNotTaggedTrkJets == 0 ) nBeforeTagSmart += evt_scale1fb();
      // soft muon veto
      bool softMuonFound = false;
      bool softMuonFound2 = false;
      for ( unsigned int imuon = 0; imuon < mus_p4().size(); ++imuon ){
	// sprintf(buf, "muon \tpt: %0.2f, eta: %0.2f, phi: %0.2f\n", mus_p4()[imuon].pt(),
	//       mus_p4()[imuon].eta(), mus_p4()[imuon].phi() );
	// addEntryToMap( map, mus_p4()[imuon].eta(), mus_p4()[imuon].phi(), buf);
	if ( hyp_type()[hypi]==1 && int(imuon) == hyp_lt_index()[hypi] ) continue;
	if ( hyp_type()[hypi]==2 && int(imuon) == hyp_ll_index()[hypi] ) continue;
	if ( mus_p4()[imuon].pt() < 20 ) softMuonFound = true;
	softMuonFound2 = true;
	hSoftMuon->Fill( mus_p4()[imuon].pt() );
      }
      if ( !softMuonFound && nTrkJets == 0 ) nFinalSoftTag += evt_scale1fb();
      if ( !softMuonFound2 && nTrkJets == 0) nFinal += evt_scale1fb();
      if ( !softMuonFound2 ) nFinalNoTrkJetsTag += evt_scale1fb();
      nFinalNoTrkJets += evt_scale1fb();

      // if ( elIso >0.92 && !softMuonFound ) nNew4 += evt_scale1fb();
      /*
      for ( unsigned int i = 0; i < 10; ++i )
	for ( unsigned int j = 0; j < 10; ++j )
	  cout << map[i][j];
      std::cout << std::endl;
      */
    }
  }
  
  printf("February selection: \t%0.2f\n", nFebSelection*kFactor);
  printf("+ calo iso: \t%0.2f\n", nFinalNoTrkJets*kFactor);
  printf("+ calo iso + muon tag: \t%0.2f\n", nFinalNoTrkJetsTag*kFactor);
  printf("+ calo iso + track jet veto - tagged jets: \t%0.2f\n", nBeforeTagSmart*kFactor);
  printf("+ calo iso + track jet veto: \t%0.2f\n", nNewNotTagged*kFactor);
  printf("+ calo iso + track jet veto + soft muon: \t%0.2f\n", nFinalSoftTag*kFactor);
  printf("+ calo iso + track jet veto + extra muon: \t%0.2f\n", nFinal*kFactor);
  printf("muon tagging efficiency no trk jet veto: %0.2f%%\n", (1-nFinalNoTrkJetsTag/nFinalNoTrkJets)*100);
  printf("advanced muon tagging efficiency: %0.2f%%\n", (1-nFinal/nBeforeTagSmart)*100);
  

  TCanvas* c2 = new TCanvas("c2","",1200,400);
  c2->Divide(3,1);
  c2->cd(1);
  hTrkJetTrackPt->Draw();
  c2->cd(2);
  hTrkJetTrackPtMax->Draw();
  c2->cd(3);
  hSoftMuon->Draw();

  TF1* f1 = new TF1("f1","15/x",0,100);
  TCanvas* c1 = new TCanvas("c1","",800,800);
  c1->Divide(2,2);
  c1->cd(1);
  hTrkJetResponseSoft->Draw();
  f1->Draw("same");
  c1->cd(2);
  hTrkJetResponseMediumSoft->Draw();
  f1->Draw("same");
  c1->cd(3);
  hTrkJetResponseMediumHard->Draw();
  f1->Draw("same");
  c1->cd(4);
  hTrkJetResponseHard->Draw();
  f1->Draw("same");

  TCanvas* c3 = new TCanvas("c3","",800,400);
  c3->Divide(2,1);
  c3->cd(1);
  hTrkJetResponseTagged->Draw();
  c3->cd(2);
  hTrkJetResponseNotTagged->Draw();
}

void EventListCreator::addEntryToMap( std::vector<std::vector<std::string> >& map, 
				      const double eta, const double phi,
				      const char* text ){
  unsigned int ieta(0);
  if ( eta > 3.0 ) ieta = map.size()-1;
  if ( eta > -3.0 && eta < 3.0 )
    ieta = int((eta+3.0)/6.0001*map.size());

  unsigned int iphi(0);
  iphi = int((phi+3.1416)/3.1417/2*map[0].size());
  map[ieta][iphi] += text;
}

void EventListCreator::test_mt( TTree* tree, std::string name)
{
  Init(tree);
  TEventList* list = getAllUniqueEvents(tree);
  unsigned int nEvents = tree->GetEntries();
  unsigned int nSelected(0);
  unsigned int nAllCuts(0);
  TH1F* hMt = new TH1F((name+"hMt").c_str(),"Min Mt",50,0,100);
  hMt->SetDirectory(0);
  std::cout << hMt << std::endl;
  TH1F* hPMet = new TH1F((name+"hPMet").c_str(),"Projected MET",50,0,100);
  hPMet->SetDirectory(0);
  std::cout << hPMet << std::endl;
  
  for( unsigned int event = 0; event < nEvents; ++event) {
    if ( !list->Contains(event) ) continue;
    GetEntry(event);
    for ( unsigned int hypi = 0; hypi < hyp_p4().size(); ++hypi ){
      bool nominalSelection 
	= (((ww_ltgoodmuiso()[hypi]==1 && ww_llgoodel()[hypi]==1 && hyp_type()[hypi]==1)||
	    (ww_ltgoodel()[hypi]==1 && ww_llgoodmuiso()[hypi]==1 && hyp_type()[hypi]==2)) &&
	   ww_pass2met()[hypi]==1 &&
	   // ww_pass4met()[hypi]==1 &&
	   ww_passzveto()[hypi]==1 && ww_passaddzveto()[hypi]==1 && 
	   hyp_ll_p4()[hypi].Pt()>20 && ww_oppsign()[hypi]==1 && hyp_njets()[hypi]==0 );
      if ( ! nominalSelection ) continue;
      // isolation
      double elTrkIso = 0;
      if (  hyp_type()[hypi]==1 )
	elTrkIso = hyp_ll_p4()[hypi].pt()/(hyp_ll_p4()[hypi].pt()+els_tkIso()[hyp_ll_index()[hypi]]);
      else if ( hyp_type()[hypi]==2 )
	elTrkIso = hyp_lt_p4()[hypi].pt()/(hyp_lt_p4()[hypi].pt()+els_tkIso()[hyp_lt_index()[hypi]]);
      
      double elIso = 0;
      if (  hyp_type()[hypi]==1 )
	elIso = hyp_ll_p4()[hypi].pt()/(hyp_ll_p4()[hypi].pt()+
					els_tkIso()[hyp_ll_index()[hypi]]+
					els_ecalJuraIso()[hyp_ll_index()[hypi]]+
					els_hcalConeIso()[hyp_ll_index()[hypi]]);
      else if ( hyp_type()[hypi]==2 )
	elIso = hyp_lt_p4()[hypi].pt()/(hyp_lt_p4()[hypi].pt()+
					els_tkIso()[hyp_lt_index()[hypi]]+
					els_ecalJuraIso()[hyp_lt_index()[hypi]]+
					els_hcalConeIso()[hyp_lt_index()[hypi]]);
      if ( elIso < 0.9 ) continue;

      // loop over trkJets
      unsigned int nTrkJets(0);
      for ( unsigned int ijet = 0; ijet < trkjets_p4().size(); ++ijet ){
	if ( trkjets_p4()[ijet].pt() < 15 ) continue;
	double dlt = deltaR( hyp_lt_p4()[hypi].eta(), hyp_lt_p4()[hypi].phi(),
			    trkjets_p4()[ijet].eta(), trkjets_p4()[ijet].phi() );
	if ( dlt < 0.4 ) continue;
	double dll = deltaR( hyp_ll_p4()[hypi].eta(), hyp_ll_p4()[hypi].phi(),
			     trkjets_p4()[ijet].eta(), trkjets_p4()[ijet].phi() );
	if ( dll < 0.4 ) continue;
	++nTrkJets;
      }
      if ( nTrkJets>0 ) continue;
      ++nSelected;
      double metX = hyp_met()[hypi]*cos(hyp_metPhi()[hypi]);
      double metY = hyp_met()[hypi]*sin(hyp_metPhi()[hypi]);
      double mt1 = getMt( metX, metY, hyp_met()[hypi],
			  hyp_lt_p4()[hypi].px(), hyp_lt_p4()[hypi].py(), hyp_lt_p4()[hypi].Et() );
      double mt2 = getMt( metX, metY, hyp_met()[hypi],
			  hyp_ll_p4()[hypi].px(), hyp_ll_p4()[hypi].py(), hyp_ll_p4()[hypi].Et() );
      if ( mt2 > mt1 )
	hMt->Fill( mt1 );
      else
	hMt->Fill( mt2 );
      
      // compute pMET
      double tightDPhi = TMath::Min(TMath::Abs(hyp_lt_p4()[hypi].Phi() - hyp_metPhi()[hypi]), 
				    2*TMath::Pi() - TMath::Abs(hyp_lt_p4()[hypi].Phi() - hyp_metPhi()[hypi]));
      double looseDPhi = TMath::Min(TMath::Abs(hyp_ll_p4()[hypi].Phi() - hyp_metPhi()[hypi]), 
				    2*TMath::Pi() - TMath::Abs(hyp_ll_p4()[hypi].Phi() - hyp_metPhi()[hypi]));
      double minDPhi = TMath::Min(tightDPhi, looseDPhi);
      double pMet = hyp_met()[hypi];
      if (minDPhi < TMath::Pi()/2) pMet = pMet*TMath::Sin(minDPhi);
      hPMet->Fill( pMet );
      // soft muon veto
      bool softMuonFound = false;
      bool softMuonFound2 = false;
      for ( unsigned int imuon = 0; imuon < mus_p4().size(); ++imuon ){
	if ( hyp_type()[hypi]==1 && int(imuon) == hyp_lt_index()[hypi] ) continue;
	if ( hyp_type()[hypi]==2 && int(imuon) == hyp_ll_index()[hypi] ) continue;
	if ( mus_p4()[imuon].pt() < 20 ) softMuonFound = true;
	softMuonFound2 = true;
      }
      if ( ww_pass4met()[hypi]==1 ) ++nAllCuts;
    }
  }
  
  printf("Unique events: \t%d\n",list->GetN());
  printf("Passed seleciton: \t%d\n", nSelected);
  printf("Passed all seleciton cuts include pMET: \t%d\n", nAllCuts);
  TCanvas* c1 = new TCanvas((name+"c").c_str(),"c1",1000,500);
  c1->Divide(2,1);
  c1->cd(1);
  hMt->Draw();
  c1->cd(2);
  hPMet->Draw();
}

double EventListCreator::getMt( double px1, double py1, double et1,
				double px2, double py2, double et2 )
{
  return sqrt( (et1+et2)*(et1+et2) - (px1+px2)*(px1+px2) - (py1+py2)*(py1+py2) );
}
