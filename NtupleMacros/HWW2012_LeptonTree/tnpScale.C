#include "../Tools/goodrun.cc"
#include "../Tools/vtxreweight.cc"

// PU reweighting
TFile* file = TFile::Open("/home/users/jaehyeok/scratch/puWeights_Summer12_5000ipb.root");
TH1D *hpu = (TH1D*) file->Get("puWeights");  

//float puweight(int nvtx, TH1D* hpu) { 
float puweight(int nvtx) { 
	
	if( nvtx > hpu->GetNbinsX() ) nvtx = hpu->GetNbinsX();

	float weight = 0;
  	weight = hpu->GetBinContent( hpu->FindBin(nvtx) );

	return weight;
}

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

void tnpScale(bool doele) {

	// 
	// Files
	// 
	TChain *chmc = new TChain("leptons");
	chmc->Add("../BatchSubmitCMS2/merged/merged_LeptonTree_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2_AOD.root");
	TChain *chdata = new TChain("leptons");
	if(doele)  chdata->Add("../BatchSubmitCMS2/merged/merged_LeptonTree_SingleElectron_Run2012AB-PromptReco-v1_AOD_3p6ifb.root");
	if(!doele) chdata->Add("../BatchSubmitCMS2/merged/merged_LeptonTree_SingleMu_Run2012AB-PromptReco-v1_AOD_3p6ifb.root");

	//
	// bins 
	//
	float ptbin[] = {10., 15., 20., 30., 40., 50., 7000.}; 	int nptbin=6;
	if(doele) 	{ float etabin[] = {0, 0.8, 1.479, 2.0, 2.5}; 	int netabin=4; }
	if(!doele) 	{ float etabin[] = {0, 0.8, 1.2, 2.4};			int netabin=3; }

	//
	// histogram
	//
	//deno
	TH2F *hmcid_deno 	= new TH2F("hmcid_deno", "hmcid_deno", nptbin, ptbin, netabin, etabin);
	TH2F *hmciso_deno 	= new TH2F("hmciso_deno", "hmciso_deno", nptbin, ptbin, netabin, etabin);
	TH2F *hdataid_deno 	= new TH2F("hdataid_deno", "hdataid_deno", nptbin, ptbin, netabin, etabin);
	TH2F *hdataiso_deno	= new TH2F("hdataiso_deno", "hdataiso_deno", nptbin, ptbin, netabin, etabin);
	hmcid_deno->Sumw2();
	hmciso_deno->Sumw2();
	hdataid_deno->Sumw2();
	hdataiso_deno->Sumw2();
	//num
	TH2F *hmcid_num 	= new TH2F("hmcid_num", "hmcid_num", nptbin, ptbin, netabin, etabin);
	TH2F *hmciso_num 	= new TH2F("hmciso_num", "hmciso_num", nptbin, ptbin, netabin, etabin);
	TH2F *hdataid_num 	= new TH2F("hdataid_num", "hdataid_num", nptbin, ptbin, netabin, etabin);
	TH2F *hdataiso_num 	= new TH2F("hdataiso_num", "hdataiso_num", nptbin, ptbin, netabin, etabin);
	hmcid_num->Sumw2();
	hmciso_num->Sumw2();
	hdataid_num->Sumw2();
	hdataiso_num->Sumw2();
	// eff
	TH2F *hmcid 	= new TH2F("hmcid", "hmcid", nptbin, ptbin, netabin, etabin);
	TH2F *hmciso 	= new TH2F("hmciso", "hmciso", nptbin, ptbin, netabin, etabin);
	TH2F *hdataid 	= new TH2F("hdataid", "hdataid", nptbin, ptbin, netabin, etabin);
	TH2F *hdataiso 	= new TH2F("hdataiso", "hdataiso", nptbin, ptbin, netabin, etabin);
	hmcid->Sumw2();
	hmciso->Sumw2();
	hdataid->Sumw2();
	hdataiso->Sumw2();
	// SF
	TH2F *hsfid 	= new TH2F("hsfid", "hsfid", nptbin, ptbin, netabin, etabin);
	TH2F *hsfiso 	= new TH2F("hsfiso", "hsfiso", nptbin, ptbin, netabin, etabin);
	hsfid->Sumw2();
	hsfiso->Sumw2();

	// nvtx check
	TH1F *hnvtxdata	=	new TH1F("hnvtxdata", 	"hnvtxdata", 100, -0.5, 99.5);
	TH1F *hnvtxmcrew	=	new TH1F("hnvtxmcrew", 	"hnvtxmcrew", 100, -0.5, 99.5);

	//
	// cuts
	//
	if(doele) {
		TString tnpcut 	= "abs(tagAndProbeMass-91)<15&&(eventSelection&1)==1&&qProbe*qTag<0&&tag->pt()>32&&probe->pt()>10&&abs(tag->eta())<2.5&&abs(probe->eta())<2.5";  // electron
		TString elfo 	= "((leptonSelection&4)==4)";  		// ele fo 
		TString elid 	= "((leptonSelection&8)==8)";  		// ele id 
		TString eliso 	= "((leptonSelection&16)==16)";  	// ele iso
		TString foriso_deno 	= tnpcut+"&&"+elid;
		TString forid_deno 		= tnpcut+"&&"+eliso;
		TString num 			= tnpcut+"&&"+elid+"&&"+eliso;
		TString trg	=	"HLT_Ele27_WP80_tag>0";
	}
	if(!doele) {
		TString tnpcut 	= "abs(tagAndProbeMass-91)<30&&(eventSelection&2)==2&&qProbe*qTag<0&&tag->pt()>30&&probe->pt()>10&&abs(tag->eta())<2.4&&abs(probe->eta())<2.4";  // muon
		TString mufo 	= "((leptonSelection&32768)==32768)";  		// mu fo 
		TString muid 	= "((leptonSelection&65536)==65536)";  		// mu id 
		TString muiso 	= "((leptonSelection&131072)==131072)";  	// mu iso 
		TString foriso_deno 	= tnpcut+"&&"+muid;
		TString forid_deno 		= tnpcut+"&&"+muiso;
		TString num 			= tnpcut+"&&"+muid+"&&"+muiso;
		TString trg	=	"HLT_IsoMu24_eta2p1_tag>0";
	}
	TString puWeight	=	"puweight(nvtx)";

	if(doele) 	cout << "..... doing Electron ..... " << endl;
	if(!doele) 	cout << "..... doing Muon ..... " << endl;

	cout << "Total MC yields 	: " << chmc->GetEntries(tnpcut) << endl;
	cout << "Total DATA yields 	: " << chdata->GetEntries(tnpcut) << endl;
	cout << "Basic selection : " << tnpcut << endl;
	//
	// Fill histograms
	//
	chmc->Draw("abs(probe->eta()):probe->pt()>>hmcid_deno", 		"("+forid_deno+")*"+puWeight,	"goff");
	chmc->Draw("abs(probe->eta()):probe->pt()>>hmcid_num", 			"("+num+")*"+puWeight,			"goff");
	chmc->Draw("abs(probe->eta()):probe->pt()>>hmciso_deno", 		"("+foriso_deno+")*"+puWeight,	"goff");
	chmc->Draw("abs(probe->eta()):probe->pt()>>hmciso_num", 		"("+num+")*"+puWeight,			"goff");
	chdata->Draw("abs(probe->eta()):probe->pt()>>hdataid_deno", 	forid_deno+"&&"+trg,			"goff");
	chdata->Draw("abs(probe->eta()):probe->pt()>>hdataid_num", 		num+"&&"+trg,					"goff");
	chdata->Draw("abs(probe->eta()):probe->pt()>>hdataiso_deno", 	foriso_deno+"&&"+trg,			"goff");
	chdata->Draw("abs(probe->eta()):probe->pt()>>hdataiso_num", 	num+"&&"+trg,					"goff");
	// vtx reweighting check	
	chdata->Draw("nvtx>>hnvtxdata",									tnpcut,			"goff");
	//chmc->Draw("nvtx>>hnvtxmcrew",									"("+tnpcut+")*"+puWeight,	"goff");
	chmc->Draw("nvtx>>hnvtxmcrew",									tnpcut,	"goff");


	// get efficiencies 
	hmcid->Divide(hmcid_num,hmcid_deno,1,1,"B");
	hmciso->Divide(hmciso_num,hmciso_deno,1,1,"B");
	hdataid->Divide(hdataid_num,hdataid_deno,1,1,"B");
	hdataiso->Divide(hdataiso_num,hdataiso_deno,1,1,"B");
	
	// get scale factors
	hsfid->Divide(hdataid, hmcid, 1, 1);
	hsfiso->Divide(hdataiso, hmciso, 1, 1);

	// Draw histograms	
	//hmcid->Draw("text");
	
	// print table
	cout << " ------ MC ID ----- " << endl;
	printline(hmcid);
	cout << " ------ MC ISO ----- " << endl;
	printline(hmciso);
	cout << " ------ DATA ID ----- " << endl;
	printline(hdataid);
	cout << " ------ DATA ISO ----- " << endl;
	printline(hdataiso);
	cout << " ------ Scale Factor ID ----- " << endl;
	printline(hsfid);
	cout << " ------ Scale Factor ISO ----- " << endl;
	printline(hsfiso);
	
//	hnvtxdata->DrawNormalized();
//	hnvtxmcrew->DrawNormalized("P");
}
