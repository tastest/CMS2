
{

	gROOT->ProcessLine(".L ../../Tools/Utilities.cc+");
	gROOT->ProcessLine(".L ../../Tools/DataSource.cc+");
	gROOT->ProcessLine(".L ../../Tools/AllDataSources.cc+");
	gROOT->ProcessLine(".L ../../Tools/HistogramUtilities.cc+");
	gROOT->ProcessLine(".L plotResults.cc+");

    plotResultsW("eb", "histos_mc");
    plotResultsW("ee", "histos_mc");



}

