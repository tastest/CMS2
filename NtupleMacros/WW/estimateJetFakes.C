{
bool runWjetBackground1 = false;
bool runWjetBackground2 = true;

// Load various tools
gROOT->ProcessLine(".x init.C");

TFile* f1 = TFile::Open("processed_data.root") ;
if ( !f1 ) {
  std::cout << "Failed to open processed_data.root" << std::endl;
  exit(1);
}

RooDataSet* dataset = (RooDataSet*)f1->Get("fulldataset");
if ( !dataset ) {
  std::cout << "Failed to get dataset" << std::endl;
  exit(1);
}

TFile* f2 = TFile::Open("fakeIsoControlSamples.root") ;
if ( !f2 ) {
  std::cout << "Failed to open fakeIsoControlSamples.root" << std::endl;
  exit(1);
}

TH1F* h_electron_qcd30  = (TH1F*)f2->Get("h_electron_qcd30");
assert(h_electron_qcd30);
TH1F* h_muon_qcd30      = (TH1F*)f2->Get("h_muon_qcd30");
assert(h_muon_qcd30);
TH1F* h_electron_qcd80  = (TH1F*)f2->Get("h_electron_qcd80");
assert(h_electron_qcd80);
TH1F* h_muon_qcd80      = (TH1F*)f2->Get("h_muon_qcd80");
assert(h_muon_qcd80);
TH1F* h_electron_qcd170 = (TH1F*)f2->Get("h_electron_qcd170");
assert(h_electron_qcd170);
TH1F* h_muon_qcd170     = (TH1F*)f2->Get("h_muon_qcd170");
assert(h_muon_qcd170);

if ( h_electron_qcd30 && h_muon_qcd30 &&
     h_electron_qcd80 && h_electron_qcd80 &&
     h_electron_qcd170 && h_muon_qcd170 ){
  if ( runWjetBackground1 ){ 
    TCanvas* c1 = new TCanvas("wjetsBackgroundEstimates_sidebandfit","",800,800);
    c1->Divide(2,2);
    c1->cd(1);
    fit_isolation(dataset,0,2,"Wjets e-fake background (pdf2)");
    c1->cd(2);
    fit_isolation(dataset,0,1,"Wjets e-fake background (pdf1)");
    c1->cd(3);
    fit_isolation(dataset,1,2,"Wjets mu-fake background (pdf2)");
    c1->cd(4);
    fit_isolation(dataset,1,1,"Wjets mu-fake background (pdf1)");
  }
  if ( runWjetBackground2 ){ 
    TFile* fcs = TFile::Open("fakeIsoControlSamples.root");
    if ( fcs ){
      TCanvas* c2 = new TCanvas("wjetsBackgroundEstimates_qcd_sideband","",600,900);
      c2->Divide(2,3);
      c2->cd(1);
      fit_isolation(dataset,0,3,"Wjets e-fake background (QCD30)",h_electron_qcd30);
      c2->cd(2);
      fit_isolation(dataset,1,3,"Wjets mu-fake background (QCD30)",h_muon_qcd30);
      c2->cd(3);
      fit_isolation(dataset,0,3,"Wjets e-fake background (QCD80)",h_electron_qcd80);
      c2->cd(4);
      fit_isolation(dataset,1,3,"Wjets mu-fake background (QCD80)",h_muon_qcd80);
      c2->cd(5);
      fit_isolation(dataset,0,3,"Wjets e-fake background (QCD170)",h_electron_qcd170);
      c2->cd(6);
      fit_isolation(dataset,1,3,"Wjets mu-fake background (QCD170)",h_muon_qcd170);
    }
  }
}
}
