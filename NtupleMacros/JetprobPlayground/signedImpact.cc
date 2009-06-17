#include "Math/LorentzVector.h"
#include "Tools/signedImpact.h"
#include "CORE/CMS2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

std::pair<double, double> signedImpact (const LorentzVector &jet, unsigned int trkidx) 
{
     TVector3 pjet(jet.px(), jet.py(), jet.pz());
     TVector3 d0;
     d0.SetPtEtaPhi(fabs(cms2.trks_d0()[trkidx]), 0, cms2.trks_vertexphi()[trkidx]);
     double signedImpact = d0.Dot(pjet) / pjet.Mag();
     double signedImpactErr = fabs(signedImpact / d0.Mag()) * cms2.trks_d0Err()[trkidx];
     return std::pair<double, double>(signedImpact, signedImpactErr);
}

void getSignedImpact(std::pair<LorentzVector, std::vector<unsigned int> > &trkjet,
		     std::vector<double>&          signedImpactVec,
                     std::vector<double>&          signedImpactErrVec,
		     std::vector<unsigned int>&    trackIndexVec)
{
  TLorentzVector trkJet(trkjet.first.Px(), trkjet.first.Py(),
                        trkjet.first.Pz(), trkjet.first.E());

  // Loop over all tracks and match to the track jet
  for (unsigned int trkIter = 0; trkIter < trkjet.second.size(); ++trkIter)
  {
    TLorentzVector trk(cms2.trks_trk_p4()[trkIter].Px(), cms2.trks_trk_p4()[trkIter].Py(),
                       cms2.trks_trk_p4()[trkIter].Pz(), cms2.trks_trk_p4()[trkIter].E() );

    if (trkJet.DeltaR(trk) < 0.5)
    {
      // This track matches the track jet, make some quality cuts
      if (fabs( cms2.trks_d0()[trkIter])   < 0.1 &&
	  cms2.trks_chi2()[trkIter]/cms2.trks_ndof()[trkIter] < 4 &&
          cms2.trks_validHits()[trkIter]   > 6   &&
          cms2.trks_trk_p4()[trkIter].Pt() > 2.0 &&
          cms2.trks_trk_p4()[trkIter].Pt() < 50.0 )
      {
        TVector3 d0Vector(fabs(cms2.trks_d0()[trkIter]) * cos(cms2.trks_vertexphi()[trkIter]),
                          fabs(cms2.trks_d0()[trkIter]) * sin(cms2.trks_vertexphi()[trkIter]),
                          0);

        TVector3 trkJet3Vec(trkJet.Px(), trkJet.Py(), 0);

        // Calculate signed impact parameter between track and track jet
        double signedImpact = d0Vector.Dot(trkJet3Vec) / trkJet3Vec.Mag();
        double signedImpactErr = fabs((signedImpact / d0Vector.Mag())) * cms2.trks_d0Err()[trkIter];

        // Add the results to be returned to the caller
        signedImpactVec.push_back(signedImpact);
        signedImpactErrVec.push_back(signedImpactErr);
        trackIndexVec.push_back(trkIter);
      }
    }
  }
}

double getJetProb(std::vector<double>& signedImpactVec,
                  std::vector<double>& signedImpactErrVec,
                  double               maxSignedImpact,
                  double               minSDSignificance)
{
  double jetProbRet = 0.0;

  for (unsigned int i = 0; i < signedImpactVec.size(); ++i)
  {
    double jetProb = signedImpactVec[i] / signedImpactErrVec[i];

    if (signedImpactVec[i] < maxSignedImpact &&
        jetProb > minSDSignificance) 
    {
      jetProbRet += jetProb;
    } 
  }

  return jetProbRet;
}
