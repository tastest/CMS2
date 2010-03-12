
{

	gROOT->ProcessLine(".L ../../Tools/Utilities.cc+");
	gROOT->ProcessLine(".L ../../Tools/DataSource.cc+");
	gROOT->ProcessLine(".L ../../Tools/AllDataSources.cc+");
	gROOT->ProcessLine(".L ../../Tools/HistogramUtilities.cc+");
	gROOT->ProcessLine(".L plotResults.cc+");

/*
    plotResultsW("EB", "histos_eleid_pt20up", "pt20up_v0_");
    plotResultsW("EE", "histos_eleid_pt20up", "pt20up_v0_");

    plotResultsW("EB", "histos_eleid_pt10to20", "pt10to20_v0_");
    plotResultsW("EE", "histos_eleid_pt10to20", "pt10to20_v0_");
*/

    plotResultsW("EB", "histos_eleid_pt10up", "pt10up_v0_");
    plotResultsW("EE", "histos_eleid_pt10up", "pt10up_v0_");



}

