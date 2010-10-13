//----------------------------------------------------
// This macro is a template that you can use to make
// your own 2D fake rate.
// It all starts from a baby ntuple.  The baby ntuple
// has one entry per lepton and is made using the 
// code in myBabyMaker.C.
//
// The content of the ntuple is documented in the 
// header file.
//
// This macro reads a baby ntuple and makes
// a muon fake rate with some sensible cuts.
// It should be easy to follow modify to do 
// what ever else you want to do
//--------------------------------------------------

{
//---------------------------
// Load some useful tools
//----------------------------
gROOT->LoadMacro("eff2.C");
gStyle->SetOptStat(0);

//-------------------------------------------
// Here you load the lepton data that
// you want to use to make your fake rate.
// It should be a baby ntuple.  
// Make sure you have one available
//-------------------------------------------
TChain *ch2 = new TChain("tree");
ch2->Add("Mu_840nb.root");

//----------------------------------------------------
// Here we define cuts for numerator and denominator
//-----------------------------------------------------
// A cut against Ws
TCut notWCut = "tcmet<20 && mt<25"; 

// A pt cut...
// Remember, we use 35 (at least for muons) to
// minimize the impact of Ws
TCut ptCut   = "pt>10 && pt<35";

// Only consider events with >=1 jet of uncorrected
// calojet pt above a threshold separated by at least
// dR...here the threshold is set to 15 GeV
TCut jetCut  = "ptj1>15";

// The trigger selection
TCut trgCut = "mu9>1";

// The numerator selection
TCut isNum = "num&&abs(id)==13"+trgCut+jetCut+ptCut+notWCut;

//The denominator selection
TCut isDenom = "fo_04&&abs(id)==13"+trgCut+jetCut+ptCut+notWCut;

//-------------------------------------------------------
// Now you define the pt and eta bins for your fake rate
//-------------------------------------------------------
double ybin[6]={10.,15.,20.,25.,30., 35.};
int nbinsy = 5;
double xbin[5]={0.0, 1.0, 1.479, 2.0, 2.5};
int nbinsx = 4;

//--------------------------------------------------------
// Book your numerator and denominator histograms
//--------------------------------------------------------
TH2F* num = new TH2F("num","num",  nbinsx, xbin, nbinsy, ybin);
TH2F* fo  = new TH2F("fo", "fo",   nbinsx, xbin, nbinsy, ybin);

//------------------------------------------
// Fill the Histograms
//-------------------------------------------
 ch2->Draw("pt:abs(eta)>>num",isNum);
 ch2->Draw("pt:abs(eta)>>fo", isDenom);

//------------------------------------------
// Get the fake rate
// The output histogram name is "fr"
//------------------------------------------
 TH2F* fr = eff2(fo,num,"fr");

}
