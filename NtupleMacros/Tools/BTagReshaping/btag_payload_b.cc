#include <vector>
// #include <iostream>
#include "btag_payload_b.h"

namespace beff
{

  void setpt()
  {
	//20-30 is also covered for mistag
	// float ptmin[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
	// float ptmax[] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
	// size_t bins=16;
	ptmin.clear();
	ptmin.push_back(20.0);
	ptmin.push_back(30.0);
	ptmin.push_back(40.0);
	ptmin.push_back(50.0);
	ptmin.push_back(60.0);
	ptmin.push_back(70.0);
	ptmin.push_back(80.0);
	ptmin.push_back(100.0);
	ptmin.push_back(120.0);
	ptmin.push_back(160.0);
	ptmin.push_back(210.0);
	ptmin.push_back(260.0);
	ptmin.push_back(320.0);
	ptmin.push_back(400.0);
	ptmin.push_back(500.0);
	ptmin.push_back(600.0);

	ptmax.clear();
	ptmax.push_back(30.0);
	ptmax.push_back(40.0);
	ptmax.push_back(50.0);
	ptmax.push_back(60.0);
	ptmax.push_back(70.0);
	ptmax.push_back(80.0);
	ptmax.push_back(100.0);
	ptmax.push_back(120.0);
	ptmax.push_back(160.0);
	ptmax.push_back(210.0);
	ptmax.push_back(260.0);
	ptmax.push_back(320.0);
	ptmax.push_back(400.0);
	ptmax.push_back(500.0);
	ptmax.push_back(600.0);
	ptmax.push_back(800.0);
	return;
  }

  void setbins()
  {
	bins = 16;
  }

  void seteta()
  {
	etamin.clear();
	etamin.push_back(0.0);
	etamin.push_back(0.5);
	etamin.push_back(0.8);
	etamin.push_back(1.0);
	etamin.push_back(1.5);

	etamax.clear();
	etamax.push_back(0.5);
	etamax.push_back(0.8);
	etamax.push_back(1.0);
	etamax.push_back(1.5);
	etamax.push_back(2.5);
  }

  void setCSVL_SFb_error()
  {
	CSVL_SFb_error.clear();
	CSVL_SFb_error.push_back(0.0484285);
	CSVL_SFb_error.push_back(0.0126178);
	CSVL_SFb_error.push_back(0.0120027);
	CSVL_SFb_error.push_back(0.0141137);
	CSVL_SFb_error.push_back(0.0145441);
	CSVL_SFb_error.push_back(0.0131145);
	CSVL_SFb_error.push_back(0.0168479);
	CSVL_SFb_error.push_back(0.0160836);
	CSVL_SFb_error.push_back(0.0126209);
	CSVL_SFb_error.push_back(0.0136017);
	CSVL_SFb_error.push_back(0.019182);
	CSVL_SFb_error.push_back(0.0198805);
	CSVL_SFb_error.push_back(0.0386531);
	CSVL_SFb_error.push_back(0.0392831);
	CSVL_SFb_error.push_back(0.0481008);
	CSVL_SFb_error.push_back(0.0474291); 
  }

  void setCSVM_SFb_error()
  {
	CSVM_SFb_error.clear();
	CSVM_SFb_error.push_back(0.0554504);
	CSVM_SFb_error.push_back(0.0209663);
	CSVM_SFb_error.push_back(0.0207019);
	CSVM_SFb_error.push_back(0.0230073);
	CSVM_SFb_error.push_back(0.0208719);
	CSVM_SFb_error.push_back(0.0200453);
	CSVM_SFb_error.push_back(0.0264232);
	CSVM_SFb_error.push_back(0.0240102);
	CSVM_SFb_error.push_back(0.0229375);
	CSVM_SFb_error.push_back(0.0184615);
	CSVM_SFb_error.push_back(0.0216242);
	CSVM_SFb_error.push_back(0.0248119);
	CSVM_SFb_error.push_back(0.0465748);
	CSVM_SFb_error.push_back(0.0474666);
	CSVM_SFb_error.push_back(0.0718173);
	CSVM_SFb_error.push_back(0.0717567);
  }

  void setCSVT_SFb_error()
  {
	CSVT_SFb_error.clear();
	CSVT_SFb_error.push_back(0.0567059);
	CSVT_SFb_error.push_back(0.0266907);
	CSVT_SFb_error.push_back(0.0263491);
	CSVT_SFb_error.push_back(0.0342831);
	CSVT_SFb_error.push_back(0.0303327);
	CSVT_SFb_error.push_back(0.024608);
	CSVT_SFb_error.push_back(0.0333786);
	CSVT_SFb_error.push_back(0.0317642);
	CSVT_SFb_error.push_back(0.031102);
	CSVT_SFb_error.push_back(0.0295603);
	CSVT_SFb_error.push_back(0.0474663);
	CSVT_SFb_error.push_back(0.0503182);
	CSVT_SFb_error.push_back(0.0580424);
	CSVT_SFb_error.push_back(0.0575776);
	CSVT_SFb_error.push_back(0.0769779);
	CSVT_SFb_error.push_back(0.0898199);
  }

  // Tagger: CSVL within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
  float CSVL_SFb(float x) 
  { 
  	return  0.981149*((1.+(-0.000713295*x))/(1.+(-0.000703264*x)));
  }

  // Tagger: CSVM within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
  float CSVM_SFb(float x) 
  { 
  	return  0.726981*((1.+(0.253238*x))/(1.+(0.188389*x)));
  }

  // Tagger: CSVT within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
  float CSVT_SFb(float x) 
  { 
  	return 0.869965*((1.+(0.0335062*x))/(1.+(0.0304598*x))); 
  }

}
