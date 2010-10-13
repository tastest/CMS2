
{

    gROOT->ProcessLine(".L ../../../Tools/Utilities.cc+");
    gROOT->ProcessLine(".L ../../../Tools/DataSource.cc+");
    gROOT->ProcessLine(".L ../../../Tools/AllDataSources.cc+");
    gROOT->ProcessLine(".L ../../../Tools/HistogramUtilities.cc+");

    gROOT->ProcessLine(".L ResultsHistograms.cc+");
    gROOT->ProcessLine(".L processDYEstResults.cc+");

    // lumi is set to 1fb in the looper)
    // Scale it to the luminosity group measurement
    // https://twiki.cern.ch/twiki/bin/view/CMS/LumiWiki2010Data

    float L = 11.52; // 1 nb
    std::string L_str = "11.52 nb^{-1}";
    float luminorm = L/(1e+06*0.01);
    plotResults("all", "histos_data", "histos_mc", L_str, luminorm);


}

