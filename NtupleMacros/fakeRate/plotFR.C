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

void plotFR( TChain* ch2, TCut numCut, TCut denCut, char* label){
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
//TChain *ch2 = new TChain("tree");
//ch2->Add("Mu.root");

//----------------------------------------------------
// Here we define cuts for numerator and denominator
//-----------------------------------------------------
// A cut against Ws
TCut notWCut = "tcmet<20 && mt<25"; 

// A pt cut...
// Remember, we use 35 (at least for muons) to
// minimize the impact of Ws
TCut ptCut  = "pt>10 && pt<35";

// The numerator selection
TCut isNum  = numCut + ptCut + notWCut;

//The denominator selection
TCut isDen  = denCut + ptCut + notWCut;

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
TH2F* num = new TH2F( Form("%s_num", label), Form("%s_num", label), nbinsx, xbin, nbinsy, ybin);
TH2F* den = new TH2F( Form("%s_den", label), Form("%s_den", label), nbinsx, xbin, nbinsy, ybin);

//------------------------------------------
// Fill the Histograms
//-------------------------------------------
 TCanvas *can = new TCanvas();
 ch2->Draw( Form("pt:abs(eta)>>%s_num", label), isNum);
 ch2->Draw( Form("pt:abs(eta)>>%s_den", label), isDen);
 delete can;

//------------------------------------------
// Get the fake rate
// The output histogram name is "fr"
//------------------------------------------
 TH2F* fr = eff2(den, num, Form("%s_fr", label) );
 return;
}
