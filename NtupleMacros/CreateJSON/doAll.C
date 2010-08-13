{

	gROOT->ProcessLine(".L ScanChain.C+");

	const string prefix = (getenv("CMS2_NTUPLE_LOCATION") != 0) ? string(getenv("CMS2_NTUPLE_LOCATION")) + "/" : "/data/tmp/";

	TChain *ch = new TChain("Events"); 

	string dataset = prefix + "MinimumBias_Commissioning10-SD_EG-Jun14thSkim_v1_RECO/V03-04-26-02/merge*.root";
	ch->Add(dataset.c_str());
	// dataset = prefix + "MinimumBias_Commissioning10-SD_Mu-Jun14thSkim_v1_RECO/V03-04-26-02/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "EG_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "Mu_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "EG_Run2010A-Jul16thReReco-v2_RECO/V03-04-26-07/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "Mu_Run2010A-Jul16thReReco-v1_RECO/V03-04-26-07/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "EG_Run2010A-PromptReco-v4_RECO/V03-04-25/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "EG_Run2010A-PromptReco-v4_RECO/V03-04-26-01/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "EG_Run2010A-PromptReco-v4_RECO/V03-04-26-07/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "EG_Run2010A-PromptReco-v4_RECO/V03-04-26-12/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "Mu_Run2010A-PromptReco-v4_RECO/V03-04-25/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-01/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-02/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-07/merge*.root";
	// ch->Add(dataset.c_str());
	// dataset = prefix + "Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-12/merge*.root";
	// ch->Add(dataset.c_str());

	ScanChain(ch); 
}