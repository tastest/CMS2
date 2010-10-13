{
//-----------------------------------
// This is an example macro that can
// be used to compare the yield of leptons
// with the prediction.
// For example, imagine that you made the FR 
// with a pt>15 GeV jet and you wanted to 
// see whether you can predict the number
// of leptons in a sample of jet 
// with pt>30 GeV.   This macro should make it easy.
// Just modify away.
//
// It needs
// (1) a standard 2D FR histogram
// (2) a standard baby ntuple for FR work
//  From (2) you can get the number of numerator
//  leptons and the number of FO taht fail the
//  numerator requirement.  From (1) you get 
//  the FR.  You put the two together and  
//  you have 
//------------------------------------

//---------------------------
// Load some tools
//---------------------------
gROOT->LoadMacro("eff2.C");
gROOT->LoadMacro("histio.cc");
gStyle->SetOptStat(0);

//--------------------------------
// Here we load the baby ntuple
// Change the name as needed
//--------------------------------
TChain *ch2 = new TChain("tree");
ch2->Add("Aug16th/Mu_840nb.root");

//-------------------------------------
// Here we load the FR histogram
// You need to know the name of the 
// FR file and the name of the 
// FR histogram.....
// For this example:
//  * FakeRate30August.root
//  * muFR15u
// You should edit as needed
//------------------------------------
 loadHist("FakeRates30August.root",0,"muFR15u");
 TH2F* thisFR = muFR15u;

//----------------------------------------------------
// These are the cuts that define the new selection
// In this example, I put the jet cut at 30 GeV.
// I also concentrate on the pt(muon)>20 GeV region
//----------------------------------------------------
// A cut against Ws
TCut notWCut = "tcmet<20 && mt<25"; 

// A pt cut...
TCut ptCut   = "pt>20";

// Only consider events with >=1 jet of uncorrected
// calojet pt above a threshold separated by at least
// dR...here the threshold is set to 30 GeV
TCut jetCut  = "ptj1>30";

// The trigger selection
TCut trgCut = "mu9>1";

// The numerator selection
TCut isNum = "num&&abs(id)==13"+trgCut+jetCut+ptCut+notWCut;

//The denominator but not denominator selection
TCut isDenomNotNum = "fo_04&&(!num)&&abs(id)==13"+trgCut+jetCut+ptCut+notWCut;

//-------------------------------------------
// Take the FR histogram and make three empty
// histograms with exactly the same binning
//-------------------------------------------
TH2F* myNum = thisFR->Clone();
myNum->Reset();
myNum->SetTitle("myNum");
myNum->SetName("myNum");
TH2F* myDenNotNum = thisFR->Clone();
myDenNotNum->Reset();
myDenNotNum->SetTitle("myDenNotNum");
myDenNotNum->SetName("myDenNotNum");
TH2F* prediction = thisFR->Clone();
prediction->Reset();
prediction->SetTitle("prediction");
prediction->SetName("prediction");


//----------------------------------------------------------
// Fill two of the three histograms.
// The 1st one is filled with the numerator
// The 2nd one is filled with the denominator but not numerator
// The 3rd one will be filled with the prediction
//----------------------------------------------------------
float ptmax = myNum->GetYaxis()->GetXmax();
cout << "Any lepton with pt > " << ptmax   << " will be put in the last bin" << endl; 
 //ch2->Draw("pt:abs(eta)>>myNum",isNum);
 // ch2->Draw("pt:abs(eta)>>myDenNotNum", isDenomNotNum);
ch2->Draw(Form("min(pt,%f-0.1):abs(eta)>>myNum",ptmax),isNum);
ch2->Draw(Form("min(pt,%f-0.1):abs(eta)>>myDenNotNum",ptmax),isDenomNotNum);

//----------------------------------------------
// Now we need to get a histogram full of ones
//-----------------------------------------------
TH2F* ones = thisFR->Clone();
ones->Reset();
ones->SetTitle("ones");
ones->SetName("ones");
 for (int ix=1; ix<ones->GetNbinsX()+1; ix++) {
   for (int iy=1; iy<ones->GetNbinsY()+1; iy++) {
     ones->SetBinContent(ix,iy,1.);
   }
 }

//----------------------------------
// Now do the actual calculation
//----------------------------------
TH2F* junk1 = thisFR->Clone();
TH2F* junk2 = thisFR->Clone();
junk1->Reset();
junk2->Reset();
junk1->SetTitle("junk1");
junk2->SetTitle("junk2");
junk1->SetName("junk1");
junk2->SetName("junk2");

junk1->Add(ones,thisFR,1,-1);  // = 1-FR
junk2->Divide(thisFR,junk1);   // = FR/(1-FR)
prediction->Multiply(myDenNotNum,junk2);

//----------------------------
// and now output the results
//-----------------------------
 cout << "Observed  = " << myNum->Integral() << endl;
 cout << "Predicted = " << prediction->Integral() << endl;

}
