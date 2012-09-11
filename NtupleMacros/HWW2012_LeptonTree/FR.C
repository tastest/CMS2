#include "../Tools/goodrun.cc"
#include "../Tools/vtxreweight.cc"

void printline(TH2F* h2)
{
	for (unsigned int x = 1; x <= h2->GetXaxis()->GetNbins(); ++x) {

		Float_t min = h2->GetXaxis()->GetBinLowEdge(x);
		Float_t max = min + h2->GetXaxis()->GetBinWidth(x);
		printf("  %4.1f - %4.1f  & ", min, max);

		for (unsigned int y = 1; y <= h2->GetYaxis()->GetNbins(); ++y)
		{
			Float_t eff = h2->GetBinContent(x, y);
			Float_t err = h2->GetBinError(x, y);
			if (y == h2->GetYaxis()->GetNbins())
				printf("\t%4.4f $\\pm$ %4.4f \\\\", eff, err);
			else
				printf("\t%4.4f $\\pm$ %4.4f & ", eff, err);
		}
		printf("\n");
	}
}

void FR(bool doele) {

	// 
	// Files
	// 
	TChain *chdata = new TChain("leptons");
	if(doele)  chdata->Add("../BatchSubmitCMS2/merged/merged_LeptonTree_DoubleElectron_Run2012AB-PromptReco-v1_AOD_3p6ifb.root");
	if(!doele) chdata->Add("../BatchSubmitCMS2/merged/merged_LeptonTree_DoubleMu_Run2012AB-PromptReco-v1_AOD_3p6ifb.root");
	//
	// bins 
	//
	float ptbin[] 	= {10., 15., 20., 25., 30., 35.}; 	int nptbin=5;
	float etabin[] 	= {0, 1.0, 1.479, 2.0, 2.5}; 	int netabin=4; 

	//
	// histogram
	//
	//deno
	TH2F *hdata_deno 	= new TH2F("hdata_deno", "hdata_deno", nptbin, ptbin, netabin, etabin);
	hdata_deno->Sumw2();
	//num
	TH2F *hdata_num 	= new TH2F("hdata_num", "hdata_num", nptbin, ptbin, netabin, etabin);
	hdata_num->Sumw2();

	// eff
	TH2F *hdata 	= new TH2F("hdata", "hdata", nptbin, ptbin, netabin, etabin);
	hdata->Sumw2();

	
	TH1F *hdr 	= new TH1F("hdr", "hdr", 100, 0, 5);
	TH1F *hdphi = new TH1F("hdphi", "hdphi", 100, 0, 6.5);

	//
	// cuts
	//
	if(doele) {
		TCut frcut 	= "met<20 && (eventSelection&16)==16 && (probe->eta())<2.5 && jet1->pt()>35";  // electron
		TCut dr		= "sqrt((jet1->eta()-probe->eta())*(jet1->eta()-probe->eta())+(jet1->phi()-probe->phi())*(jet1->phi()-probe->phi()))>1.0";  
		TCut fo 	= "((leptonSelection&4)==4)";  		// ele fo 
		TCut id 	= "((leptonSelection&8)==8)";  		// ele id 
		TCut iso 	= "((leptonSelection&16)==16)";  	// ele iso
		TCut ele8		= "HLT_Ele8_CaloIdL_TrkIdVL>0";
		TCut ele8T		= "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL>0";
		TCut ele8T30	= "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30>0";
		TCut ele17T		= "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL>0";
		TCut ele17T30	= "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30>0";
		TCut trg		= ele8 || ele8T || ele8T30;
	}
	if(!doele) {
		TCut frcut 	= "met<20 && (eventSelection&32)==32 && (probe->eta())<2.4 && jet1->pt()>15";
		TCut dr		= "sqrt((jet1->eta()-probe->eta())*(jet1->eta()-probe->eta())+(jet1->phi()-probe->phi())*(jet1->phi()-probe->phi()))>1.0";  
		TCut fo 	= "((leptonSelection&32768)==32768)";  		// mu fo 
		TCut id 	= "((leptonSelection&65536)==65536)";  		// mu id 
		TCut iso 	= "((leptonSelection&131072)==131072)";  	// mu iso 
		TCut mu8	= "HLT_Mu8>0";
		TCut mu17	= "HLT_Mu17>0";
		//TCut trg	= mu8 || mu17;
		TCut trg	= mu8;
	}
	TCut deno = frcut+dr+trg+fo;
	TCut num = deno+iso+id;

	if(doele) 	cout << "..... doing Electron ..... " << endl;
	if(!doele) 	cout << "..... doing Muon ..... " << endl;

	cout << "Total DATA yields 	: " << chdata->GetEntries(frcut+trg+fo) << endl;
	//cout << "Basic selection : " << (frcut+trg+fo).Data() << endl;
	//
	// Fill histograms
	//
	chdata->Draw("abs(probe->eta()):probe->pt()>>hdata_deno", 	deno,			"goff");
	chdata->Draw("abs(probe->eta()):probe->pt()>>hdata_num", 	num,			"goff");
	chdata->Draw("sqrt((jet1->eta()-probe->eta())*(jet1->eta()-probe->eta())+(jet1->phi()-probe->phi())*(jet1->phi()-probe->phi()))>>hdr", 	deno,			"goff");
	chdata->Draw("sqrt((jet1->phi()-probe->phi())*(jet1->phi()-probe->phi()))>>hdphi", 	deno,			"goff");


	// get efficiencies 
	hdata->Divide(hdata_num,hdata_deno,1,1,"B");
	
	// print table
	cout << " ------ FR ----- " << endl;
	printline(hdata);

}
