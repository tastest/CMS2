#include "ScanChainUtilities.h"
//don't use this if it's in ScanChain.C, but need it if this file is on its own
//#include "CMS2.cc"

using namespace tas;


bool isGoodTrk(unsigned int idx ){
  //will fill this later
  if(trks_qualityMask().at(idx) & 4){
    return true;
  } else {
    return false;
  }

}

bool passesTrigger(bool isGEN) {
  
  //if(passL1Trigger("L1_SingleHfBitCountsRing1_1"))
  //return true;
  
  //if(passL1Trigger("L1_SingleHfBitCountsRing2_1"))
  //return true;
  
  //Beam Halo triggers
  if(l1_techbits2() & (1<<4) || l1_techbits2() & (1<<5) || l1_techbits2() & (1<<6) || l1_techbits2() & (1<<7) )
    return false;

  //BPTX triggers
  if(!(l1_techbits1() & (1<<0)) && !isGEN)
    return false;

  //BSC triggers
  if(l1_techbits2() & (1<<8) || l1_techbits2() & (1 << 9))
    return true;
  
  return false;
}


bool passesTrackCuts() {
    
  int nGoodVtxs = 0;
  for(unsigned int i = 0; i < vtxs_isFake().size(); i++) {
    if(vtxs_isFake().at(i))
      continue;
    //if(vtxs_tracksSize().at(i) < 4 )
	//  continue;
	if(vtxs_ndof().at(i) < 5 ) //the calomet guys use this
	  continue;
    if( fabs( vtxs_position().at(i).z() ) > 15 )
      continue;
    //if( vtxs_position().at(i).pt() > 2 )
	//continue;
    nGoodVtxs++;
  }
  
  //require that there be at least one good vertex
  if(nGoodVtxs==0)
    return false;
  
  if(trks_trk_p4().size() <= 10){
    return true;
  } else {

    //require that the fraction of highPurity tracks be > 50%
    int nGoodTrks = 0;
    for(unsigned int i = 0; i < trks_trk_p4().size(); i++) {
      if(trks_qualityMask().at(i) & 4)
		nGoodTrks++;
    }
    //if((float)nGoodTrks/trks_trk_p4().size() < 0.2) //20%
	if((float)nGoodTrks/trks_trk_p4().size() < 0.25) //25%
      return false;
  }

  return true;
}

//utility function to get the dist and delta cot theta
std::pair<float, float> getConversionInfo(int idx1, int idx2, float bField){
  
  int trk1_q = trks_charge().at(idx1);
  int trk2_q = trks_charge().at(idx2);
  double trk1_d0 =  trks_d0().at(idx1);
  double trk2_d0 =  trks_d0().at(idx2);
  double trk1_pt = trks_trk_p4().at(idx1).pt();
  double trk2_pt = trks_trk_p4().at(idx2).pt();
  double trk1_phi = trks_trk_p4().at(idx1).phi();
  double trk2_phi = trks_trk_p4().at(idx2).phi();
  double trk1_theta = trks_trk_p4().at(idx1).theta();
  double trk2_theta = trks_trk_p4().at(idx2).theta();

  if( trk1_pt == 0 || trk2_pt == 0 || bField == 0 || trk1_q == 0 || trk2_q == 0 || trk1_theta == 0 || trk2_theta == 0) 
    cout << "about to barf because of division by zero" << endl;
 
  double tk1Curvature = -0.3*bField*(trk1_q/trk1_pt)/100.;
  double rTk1 = fabs(1./tk1Curvature);
  double xTk1 = (1./tk1Curvature - trk1_d0)*cos(trk1_phi);
  double yTk1 = (1./tk1Curvature - trk1_d0)*sin(trk1_phi);
    
  double tk2Curvature = -0.3*bField*(trk2_q/trk2_pt)/100.;
  double rTk2 = fabs(1./tk2Curvature);
  double xTk2 = (1./tk2Curvature - trk2_d0)*cos(trk2_phi);
  double yTk2 = (1./tk2Curvature - trk2_d0)*sin(trk2_phi);
         
  double dist = sqrt(pow(xTk1-xTk2, 2) + pow(yTk1-yTk2 , 2));
  dist = dist - (rTk1 + rTk2);

  double dcot = 1/tan(trk1_theta) - 1/tan(trk2_theta);

  return make_pair(dist, dcot);
  
}


//fn to reproduce what PJ calls 'S9': the (em or had) of a tower plus the em+had of the four immediate neighbors
//input: index in twrs block of seed, bool for em or had

// range of ieta is -41 <= ieta <= 41, EXCLUDING zero
// range of iphi is   1 <= iphi <= 72 in barrel (-20 <= ieta <= 20),
//                    1 <= iphi <= 71, ODD NUMBERS ONLY in endcap      : 21 <= abs(ieta) <= 39
//                    3 <= iphi <= 71, EVERY FOURTH ONLY in far forward: 40 <= abs(ieta) <= 41
//                    ie, 3, 7, 11, 15

float getTwrHFSwiss( int seedidx, bool em ) {

  if( twrs_eta()[seedidx] < 2.964 || twrs_eta()[seedidx] > 4.889 ) //see above
	return -1;

  float result = ( em ? twrs_emEnergy()[seedidx] : twrs_hadEnergy()[seedidx] );
  unsigned int added = 0;
  const int seedieta = CaloTwr_ieta( twrs_detid()[seedidx] );
  const int seediphi = CaloTwr_iphi( twrs_detid()[seedidx] );
  const int maxiphibarrel = 72;
  const int maxiphiendcap = 71;

  //cout << evt_run() << "  " << evt_lumiBlock() << "  " << evt_event() << endl;

  for( unsigned int i=0; i<twrs_eta().size(); i++ ) {
	if( int(i) == seedidx )
	  continue;

	bool phineighbors = false;
	bool etaneighbors = false;
	const int thisieta = CaloTwr_ieta( twrs_detid()[i] );
	const int thisiphi = CaloTwr_iphi( twrs_detid()[i] );
	const int maxiphi = (abs(thisieta) <= 20 ? maxiphibarrel : maxiphiendcap );
	int diffiphi = 1; //barrel
	if( abs(thisieta) > 20 && abs(thisieta) <= 39 )
	  diffiphi = 2;
	else if( abs(thisieta) > 39 )
	  diffiphi = 4;

	if( (seediphi == 1       && thisiphi == maxiphi) || //take account of periodicity of phi
		(seediphi == maxiphi && thisiphi == 1      ) ||
		abs(thisiphi-seediphi) == diffiphi ) { //abs bc doesn't matter if + or - diffiphi
	  phineighbors = true;
	}

	if( abs(thisieta-seedieta) == 1 ||
		(thisieta == 1  && seedieta == -1) ||
		(thisieta == -1 && seedieta ==  1) )
	  etaneighbors = true;
		
	if( ( seediphi == thisiphi && etaneighbors ) || //same phi, neighbors in eta
		( phineighbors && seedieta == thisieta ) ) { //same eta, neighbors in phi
	  result += twrs_emEnergy()[i] + twrs_hadEnergy()[i];
	  //if( added == 0 )
		//cout << "seed " << seedieta << "  " << seediphi << endl;
	  added++;
	  //cout << "added " << thisieta << "  " << thisiphi << endl;
	}

	/* ///////////////////////////////////////////
	//DON"T USE THIS---IT"S BROKEN
	//float dphi = ROOT::Math::VectorUtil::DeltaPhi( LorentzVector(, , , ),  );
	float dphi = deltaPhi( twrs_phi()[seedidx], twrs_phi()[i] );
	float deta = twrs_eta()[seedidx] - twrs_eta()[i];

	if( fabs(dphi) * 180./pi < 10. && fabs(deta) < 0.35 ) { //same phi, neighbors in eta
	  result += twrs_emEnergy()[i] + twrs_hadEnergy()[i]; 
	  added++;
	}
	else if( fabs(deta) < 0.175 && fabs(dphi) * 180./pi < 20. ) { //same eta, neighbors in phi
	  result += twrs_emEnergy()[i] + twrs_hadEnergy()[i];
	  added++;
	}
	*/
  }

  //cout << endl;
  if( added > 4 )
	cout << "added too many towers--fix" << endl;

  return result;

}


float deltaPhi(float phi1,float phi2){

  float deltaphi=phi1-phi2;
  if(deltaphi> acos(-1.)) deltaphi=2*acos(-1.)-deltaphi;
  if(deltaphi<-acos(-1.)) deltaphi=2*acos(-1.)+deltaphi;
  return deltaphi;

}

//port of code from CaloTowerDetId.h
int CaloTwr_ieta( int detid ) {
  int zside = (detid & 0x2000) ? 1 : -1 ;
  //warning: multiplication has higher precendnce than bitwise anding
  return zside * ((detid >> 7) & 0x3f);
}

int CaloTwr_iphi( int detid ) {
  return detid & 0x7F;
}


void progressbar( int nEventsTotal, int nEventsChain ) {
  // Progress feedback to the user
  int i_permille = (int)floor(1000 * float(nEventsTotal) / float(nEventsChain));
  if (i_permille != i_permille_old) {
	// xterm magic from L. Vacavant and A. Cerri
	if (isatty(1)) {
	  printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
			 "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	  fflush(stdout);
	}
	i_permille_old = i_permille;
  }
}



void NewHist(TH1F*& h, const char* name, const char* title, int bins, double min, double max, bool isGEN) {
  h = new TH1F(name, title, bins, min, max);
  h->Sumw2();
  if( isGEN ) {
	int gencolor = 2; //red
	h->SetFillColor(gencolor);
	h->SetLineColor(gencolor);
  }
}

void NewHist(TH2F*& h, const char* name, const char* title, int xbins, double xmin, double xmax, int ybins, double ymin, double ymax, bool isGEN) {
  h = new TH2F(name, title, xbins, xmin, xmax, ybins, ymin, ymax);
  //h->Sumw2();
  if( isGEN ) {
	int gencolor = 2; //red
	h->SetFillColor(gencolor);
	h->SetLineColor(gencolor);
  }
}

TH1F* MakeHist(const char* name, const char* title, int bins, double min, double max, bool isGEN) {
  TH1F* h = new TH1F(name, title, bins, min, max);
  h->Sumw2();
  if( isGEN ) {
	int gencolor = 2; //red
	h->SetFillColor(gencolor);
	h->SetLineColor(gencolor);
  }
  return h;
}

TH2F* MakeHist(const char* name, const char* title, int xbins, double xmin, double xmax, int ybins, double ymin, double ymax, bool isGEN) {
  TH2F* h = new TH2F(name, title, xbins, xmin, xmax, ybins, ymin, ymax);
  h->Sumw2();
  if( isGEN ) {
	int gencolor = 2; //red
	h->SetFillColor(gencolor);
	h->SetLineColor(gencolor);
  }
  return h;
}

