
{

	gROOT->ProcessLine(".L ../../Tools/Utilities.cc+");
	gROOT->ProcessLine(".L ../../Tools/DataSource.cc+");
	gROOT->ProcessLine(".L ../../Tools/AllDataSources.cc+");
	gROOT->ProcessLine(".L ../../Tools/HistogramUtilities.cc+");
	gROOT->ProcessLine(".L plotResults.cc+");
	// lumi is set to 1fb in the looper)
	// Scale it to the luminosity group measurement
	// https://twiki.cern.ch/twiki/bin/view/CMS/LumiWiki2010Data
        float lumi = 0.201;// in unit of nano barn
	float luminorm = 10*lumi/1e+06; // reference being scaled up 10 times
	plotResultsW("all", "histos_data", "histos_reference", luminorm);



}

