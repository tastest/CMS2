#include <assert.h>
#include <algorithm>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "jetSelections.h"

// CMS2 includes                                                                                                                                    
#include "CMS2.h"
#include "utilities.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

bool comparePt (const LorentzVector &lv1,
                const LorentzVector &lv2)
{
  return lv1.pt() > lv2.pt();
}

//   jets_p4   
std::vector<LorentzVector> getCaloJets(int i_hyp) {
  std::vector<LorentzVector> calo_jets;
  calo_jets.clear();
  
  for (unsigned int jj=0; jj < cms2.jets_p4().size(); ++jj) {
    if ((dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.jets_p4()[jj]) < 0.4)||
	(dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.jets_p4()[jj]) < 0.4)
	) continue;
    if (cms2.jets_cor()[jj]*cms2.jets_p4()[jj].pt() < 30) continue;
    if (fabs(cms2.jets_p4()[jj].Eta()) > 2.4) continue;
    //fkw July21 2009 if (cms2.jets_emFrac()[jj] < 0.1) continue;
    calo_jets.push_back(cms2.jets_cor()[jj]*cms2.jets_p4()[jj]);
  }
  
  if (calo_jets.size() > 1) {
       sort(calo_jets.begin(), calo_jets.end(),  comparePt);
  }
  return calo_jets;
}
