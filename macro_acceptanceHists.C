void macro_acceptanceHists(const char* rootfile){
gROOT->ProcessLine(".L acceptanceplots.C+");
acceptanceplots("lepAzimAsym2",false,Form("results/%s",rootfile));
acceptanceplots("lepAzimAsym",false,Form("results/%s",rootfile));
acceptanceplots("lepChargeAsym",false,Form("results/%s",rootfile));
acceptanceplots("topCosTheta",false,Form("results/%s",rootfile));
acceptanceplots("lepPlusCosTheta",false,Form("results/%s",rootfile));
acceptanceplots("lepMinusCosTheta",false,Form("results/%s",rootfile));
acceptanceplots("topSpinCorr",false,Form("results/%s",rootfile));
acceptanceplots("rapiditydiff",false,Form("results/%s",rootfile));
acceptanceplots("pseudorapiditydiff",false,Form("results/%s",rootfile));
acceptanceplots("rapiditydiffMarco",false,Form("results/%s",rootfile));
acceptanceplots("lepCosTheta",false,Form("results/%s",rootfile));
}
