
{

	gROOT->ProcessLine(".L ../../../Tools/Utilities.cc+");
	gROOT->ProcessLine(".L ../../../Tools/DataSource.cc+");
	gROOT->ProcessLine(".L ../../../Tools/AllDataSources.cc+");
	gROOT->ProcessLine(".L ../../../Tools/HistogramUtilities.cc+");

	gROOT->ProcessLine(".L ResultsHistograms.cc+");
	gROOT->ProcessLine(".L processDYEstResults.cc+");

	processDYEstResults("../histos_mc_3x.root");

}

