// C++ includes
#include <iostream>
#include <set>
#include <fstream>

// ROOT includes
#include "TSystem.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TChainElement.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Math/VectorUtil.h"

// TAS includes
#include "../CORE/CMS2.cc"
//#include "./CMS2.cc"
#include "./textFileMaker.h"
using namespace std;
using namespace tas;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

void textFileMaker::ScanChain( TChain* chain) {


  // Open the text file for output
  ofstream myfile ("decay.txt");

  // Sanity Histogram
  TH2F* hm  = new TH2F("xvscsthetam","xvscstheta muon",100,-1.,1.,100,0.,1.2);
  TH2F* he  = new TH2F("xvscsthetae","xvscstheta ele ",100,-1.,1.,100,0.,1.2);
  TH1F* lpt = new TH1F("lpt","lepton pt",100,0.,100.);

  // x is defined as E(lepton)/E(max)
  // E(max) is (m(tau)^2+(m(lep)^2)/(2 m(tau)
  float mtau = 1.778;
  float mmu  = 0.1066;
  float me   = 0.000511;
  float xmumax = (mtau*mtau+mmu*mmu)/(2*mtau);
  float xemax  = (mtau*mtau+me*me)/(2*mtau);

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
    for( z = 0; z < nLoop; z++) {	// Event Loop
      cms2.GetEntry(z);


      // looper progress
      ++nEventsTotal;
      int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
      if (i_permille != i_permilleOld) {
        printf("  \015\033[32m ---> \033[1m\033[31m%4.1f%%" "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
        fflush(stdout);
        i_permilleOld = i_permille;
      }
      
      // loop over status==3 and look for taus and its daughters

      for(unsigned int i = 0; i < genps_id().size(); i++) {//status 3 loop
      	if(abs(genps_id()[i]) != 15 ) continue; 
	//LorentzVector vtau = genps_p4()[i];

	// Because the status=3 and the status=1 particles are in slightly different
	// frames, build the tau 4-vector from all of its daughters
	LorentzVector vtau;
	for(unsigned int j = 0; j < genps_lepdaughter_id()[i].size(); j++) { //loop over the tau's status1 daughters
	  vtau = vtau + genps_lepdaughter_p4()[i][j];
	}
	// Sanity check.  Check the tau mass..  It worked, but not always.  Not sure why.  Ask that it works!!  
	//	cout << "tau mass  " << vtau.mass() << endl;
	if (abs(vtau.mass()-mtau) > 0.1)  continue;


	// We are looking for tau->lep nu_{lep} nu_{tau}.  So ask for three daughters.  
	// This is not perfect because of cases with radiated photons, but let's ignore that for now.
	if (genps_lepdaughter_id()[i].size() != 3) continue;

	for(unsigned int j = 0; j < genps_lepdaughter_id()[i].size(); j++) { //loop over the tau's status1 daughters
	  if (abs(genps_lepdaughter_id()[i][j]) == 15) cout << "We find a tau in the tau daughter list" << endl;
	  if (abs(genps_lepdaughter_id()[i][j]) == 11 || abs(genps_lepdaughter_id()[i][j]) == 13) { // found ele or mu

	    LorentzVector vlep = genps_lepdaughter_p4()[i][j];
	    lpt->Fill(vlep.pt());
   
	    // This is the boost back to the tau rest frame
	    ROOT::Math::Boost boost(vtau.BoostToCM().x(), vtau.BoostToCM().y(), vtau.BoostToCM().z());
 
	    // Now boost the lepton into the tau rest frame
	    LorentzVector vlepCM = boost*vlep;

	    // Sanity check: boost the tau its rest frame...check that momentum=0
	    // CHECKED: THIS WORKS, get typically less than 1 MeV
	    //LorentzVector vtauCM = boost*vtau;
	    //cout << vtauCM.P() << endl;
	  
	    // calculate the cos of the angle between the original direction of 
	    // the tau and the lepton direction in the tau rest frame
	    float costheta = (vlepCM.px()*vtau.px()+vlepCM.py()*vtau.py()+vlepCM.pz()*vtau.Pz())/(vlepCM.P()*vtau.P());
	      
	    // Fill sanity histograms
	    if (abs(genps_lepdaughter_id()[i][j]) == 11) he->Fill(costheta, vlepCM.E()/xemax);
	    if (abs(genps_lepdaughter_id()[i][j]) == 13) hm->Fill(costheta, vlepCM.E()/xmumax);


	    //Now we write out a text file with leptonType : Momentum in Tau Frame : costheta
	    myfile << abs(genps_lepdaughter_id()[i][j]) << " " << vlepCM.P() << " " << costheta << endl; 

	  } //  closes if-block over found ele or mu
	}   //  closes loop of tau daughtters
      }     //  clses loop over status==3
    }       // closes loop over events
  }         // closes loop over files

  return;
}    // closes myLooper function  



