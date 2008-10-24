// -*- C++ -*-
#ifndef WWDRIVER_H
#define WWDRIVER_H

int WW_Results ();		// plots, tables etc (old baseline + track jets + b tag + calo iso)
int WW_Results_softmu ();	// plots, tables etc (old baseline + track jets + b tag (< 20 GeV) + calo iso)
int WW_NoLLPtCut ();		// plots, tables etc (old baseline + track jets + b tag + calo iso), no ll pt cut
int WW_NoTrackJets ();		// plots, tables etc (baseline - track jets)
int WW_NoMuTag ();		// plots, tables etc (baseline - b tag)
int WW_NoCaloIso ();		// plots, tables etc (baseline - calo iso)
int WW_OldResults ();		// plots, tables etc (old baseline)
int WW_OldResults_trkjets ();		// plots, tables etc (old baseline + trkjets)
int WW_OldResults_muveto ();		// plots, tables etc (old baseline + muveto)
int WW_OldResults_caloiso ();		// plots, tables etc (old baseline + calo iso)
int WW_OldResults_caloiso_trkjets ();		// plots, tables etc (old baseline + calo iso + trkjets)
int WW_OldResults_caloiso_muveto ();		// plots, tables etc (old baseline + calo iso + muveto)
int WW_TopEstimate ();
int WW_TopEstimate_noptcut ();
int WW_TopEstimate_oldcuts ();
int WW_SS_Results (); 		// same-sign control sample 
int WW_SS_NoCaloIso (); 		// same-sign control sample 
int WW_Wjets_Dumbo ();		// get W+jets estimate from the Dumbo method
int WW_Wjets_SS_Dumbo ();	// get W+jets estimate in same-sign control sample from the Dumbo method
int WW_Wjets_Fakerates ();	// get W+jets estimate from el fake rate
int WW_Wjets_Fakerates_barrel ();	// get W+jets estimate from el fake rate (barrel only)
int WW_Wjets_Numerator ();	// count W+jets electron numerator objects
int WW_Wjets_Numerator_barrel ();	// count W+jets electron numerator objects, barrel only
int WW_Wjets_Numerator_barrel_supertight ();	// count W+jets electron numerator objects, barrel only + supertight el cuts
int WW_Wjets_Numerator_barrel_supertight_caloiso ();	// count W+jets electron numerator objects, barrel only + supertight el cuts, cutting on calo iso
int WW_Wjets_FOs ();
int WW_Wjets_FOs_supertight ();
int WW_Wjets_FOs_not_numerator ();
int WW_Wjets_SS_Fakerates ();	// get W+jets (same-sign) estimate from el fake rate
int WW_Wjets_SS_Fakerates_barrel ();	// get W+jets (same-sign) estimate from el fake rate (barrel only)
int WW_Wjets_SS_Numerator ();	// count W+jets (same-sign) electron numerator objects
int WW_Wjets_SS_Numerator_barrel ();	// count W+jets electron numerator objects, barrel only
int WW_Wjets_SS_Numerator_barrel_supertight ();	// count W+jets electron numerator objects, barrel only + supertight el cuts
int WW_Wjets_SS_Numerator_barrel_supertight_caloiso ();	// count W+jets electron numerator objects, barrel only + supertight el cuts, cutting on calo iso
int WW_Wjets_SS_FOs ();
int WW_Wjets_SS_FOs_not_numerator ();
int WW_Wjets_SS_FOs_supertight ();
int WW_DY_in_emu ();
int WW_TrkCorrMET_NoMetNoZVeto ();
int WW_TrkCorrMET_NoMet ();

#endif
