#include "TFile.h"
#include "TSystem.h"
#include "TH2.h"
#include "CORE/CMS2.h"
#include "CORE/selections.h"
#include "Tools/chargeflip.h"

static TFile *el_chargeflipfile_v1 = 0;
static TH2F  *el_chargeflip_v1 = 0;

static TFile *el_chargeflipfile_v2 = 0;
static TH2F  *el_chargeflip_v2 = 0;

double elchargeflipprob (int i_el, int add_error_times)
{
     float prob = 0.0;
     float prob_error = 0.0;
     TH2F *theFakeRate = &flipRate();
     // cut definition
     float pt = cms2.els_p4()[i_el].Pt();
     float upperEdge = theFakeRate->GetYaxis()->GetBinLowEdge(theFakeRate->GetYaxis()->GetNbins()) + theFakeRate->GetYaxis()->GetBinWidth(theFakeRate->GetYaxis()->GetNbins()) - 0.001;
     if ( pt > upperEdge )
       pt = upperEdge;
     prob = theFakeRate->GetBinContent(theFakeRate->FindBin(TMath::Abs(cms2.els_p4()[i_el].Eta()),pt));
     prob_error =
	  theFakeRate->GetBinError(theFakeRate->FindBin(TMath::Abs(cms2.els_p4()[i_el].Eta()),pt));
     
     if (prob>1.0 || prob<0.0) {
 	  std::cout<<"ERROR FROM FAKE RATE!!! prob = " << prob << std::endl;
     }
     if (prob==0.0){
 	  std::cout<<"ERROR FROM FAKE RATE!!! prob = " << prob
 		   <<" for Et = " <<cms2.els_p4()[i_el].Pt()
 		   <<" and Eta = " <<cms2.els_p4()[i_el].Eta()
 		   << std::endl;
     }
     return prob+add_error_times*prob_error;
}


#define USE_V1

TH2F &flipRate ()
{
#ifdef USE_V1
     if ( el_chargeflip_v1 == 0 ) {
	  el_chargeflipfile_v1 = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/el_chargeflip_v1.root", "read"); 
	  if ( el_chargeflipfile_v1 == 0 ) {
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/el_chargeflip-v1.root could not be found!!" << std::endl;
	       std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/el_chargeflip-v1.root exists!" << std::endl;
	       gSystem->Exit(1);
	  }
	  el_chargeflip_v1 = dynamic_cast<TH2F *>(el_chargeflipfile_v1->Get("zee_mcFLRPtvsEta"));
	  //el_flipRate_err_v5 = dynamic_cast<TH2F *>(el_flipRateFile_v5->Get("flipRateTemplateError_wo_leading_elt_flipRatesFull"));
     }
     return *el_chargeflip_v1;
#endif

#ifdef USE_V2
     if ( el_chargeflip_v2 == 0 ) {
	  el_chargeflipfile_v2 = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/el_chargeflip_v2.root", "read"); 
	  if ( el_chargeflipfile_v2 == 0 ) {
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/el_chargeflip-v2.root could not be found!!" << std::endl;
	       std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
	       std::cout << "$CMS2_LOCATION/NtupleMacros/data/el_chargeflip-v2.root exists!" << std::endl;
	       gSystem->Exit(1);
	  }
	  el_chargeflip_v2 = dynamic_cast<TH2F *>(el_chargeflipfile_v2->Get("zee_mcFLRPtvsEta"));
	  //el_flipRate_err_v5 = dynamic_cast<TH2F *>(el_flipRateFile_v5->Get("flipRateTemplateError_wo_leading_elt_flipRatesFull"));
     }
     return *el_chargeflip_v2;
#endif

}
//  LocalWords:  flipRateErrorMuon
