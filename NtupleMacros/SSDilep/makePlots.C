
{

	gROOT->ProcessLine(".L ../Tools/Utilities.cc+");
	gROOT->ProcessLine(".L ../Tools/DataSource.cc+");
	gROOT->ProcessLine(".L ../Tools/AllDataSources.cc+");
	gROOT->ProcessLine(".L ../Tools/HistogramUtilities.cc+");
	gROOT->ProcessLine(".L plotResults.cc+");

    // in looper, norm is to 10pb
    float L = 10.0; // 10 pb
    std::string L_str = "10.0 pb^{-1}";
    float luminorm = 1.0; //L/(1e+03);

    // make the plots
    plotResults("ee", "v00", "histos_data", "histos_mc", L_str, luminorm);
    plotResults("mm", "v00", "histos_data", "histos_mc", L_str, luminorm);
    plotResults("em", "v00", "histos_data", "histos_mc", L_str, luminorm);
    plotResults("all", "v00", "histos_data", "histos_mc", L_str, luminorm);


    // make the tables
    makeTables("histos_data", "histos_mc", L_str, luminorm);


}

