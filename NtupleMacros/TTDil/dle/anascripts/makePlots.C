
{

	gROOT->ProcessLine(".L ../../../Tools/Utilities.cc+");
	gROOT->ProcessLine(".L ../../../Tools/DataSource.cc+");
	gROOT->ProcessLine(".L ../../../Tools/AllDataSources.cc+");
	gROOT->ProcessLine(".L ../../../Tools/HistogramUtilities.cc+");
	gROOT->ProcessLine(".L plotResults.cc+");
	// lumi is set to 1fb in the looper)
	// Scale it to the luminosity group measurement
	// https://twiki.cern.ch/twiki/bin/view/CMS/LumiWiki2010Data

	float luminorm = 1.0/1e+06; // 1.0 nb
	plotResults("all", "../histos_data", "../histos_mc", "1 nb^{-1}", luminorm);

}

