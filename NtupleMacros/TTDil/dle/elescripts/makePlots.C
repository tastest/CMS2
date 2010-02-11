
{

	gROOT->ProcessLine(".L ../../../Tools/Utilities.cc+");
	gROOT->ProcessLine(".L ../../../Tools/DataSource.cc+");
	gROOT->ProcessLine(".L ../../../Tools/AllDataSources.cc+");
	gROOT->ProcessLine(".L ../../../Tools/HistogramUtilities.cc+");
	gROOT->ProcessLine(".L plotResults.cc+");

// first argument, "det"
// second argument, "hyp"
//	plotResultsID("ee", "ee", "../histos_mc_3x", "dyee");
 //   plotResultsID("eb", "ee", "../histos_mc_3x", "dyee");

    plotResultsID("ee", "ee", "../histos_eleid_hypbased", "eleid_v1");
    plotResultsID("eb", "ee", "../histos_eleid_hypbased", "eleid_v1");



}

