
{

	gROOT->ProcessLine(".L ../../../Tools/Utilities.cc+");
	gROOT->ProcessLine(".L ../../../Tools/DataSource.cc+");
	gROOT->ProcessLine(".L ../../../Tools/AllDataSources.cc+");
	gROOT->ProcessLine(".L ../../../Tools/HistogramUtilities.cc+");
	gROOT->ProcessLine(".L plotResults.cc+");

	plotResults("ee", "histos_mc_3x");
	plotResults("em", "histos_mc_3x");
	plotResults("mm", "histos_mc_3x");
    plotResults("all", "histos_mc_3x");

}

