TCut tightid   = "els_tightId22XMinMatteo ==1";
TCut pt20      = "els_p4.pt() > =20";
TCut forward   = "abs(els_p4.eta()) > 1.47";
TCut central   = "abs(els_p4.eta()) < 1.47";
TCut d0corr250 = "abs(els_d0corr) < 0.025";
TCut d0corr200 = "abs(els_d0corr) < 0.020";
TCut conv_dcot = "abs(els_conv_dcot)<0.02"; 
TCut conv_dist = "abs(els_conv_dist)<0.02";
TCut conversion_1 =  conv_dcot && conv_dist;
TCut notmuon   = "els_closestMuon == -1";
TCut susy_iso  = "(els_pat_trackIso+max(0., (els_pat_ecalIso-2.)) + els_pat_hcalIso)/max(els_p4.pt(),20.) < 0.1";
TCut ww_iso    = "els_p4.pt()/(els_p4.pt()+els_pat_trackIso+els_pat_ecalIso+els_pat_hcalIso) >0.92";

TCut susy_eta  = "TMath::Abs(els_p4.eta()) <= 2.4";
TCut ww_baseline_0 = pt20 && tightid && notmuon && ww_iso;
TCut ww_baseline_1 = ww_baseline_0 && d0corr250;
TCut ww_baseline_2 = ww_baseline_1 && !conversion_1;

TCut susy_baseline_0 = susy_eta&& pt20 && tightid && notmuon && susy_iso;
TCut susy_baseline_1 = susy_baseline_0 && d0corr200;
TCut susy_baseline_2 = susy_baseline_1 && !conversion_1;


// total electrons = 1.42007e+06
// total electrons in the forward region = 495648
// electrons removed by cutting on the number of pixel hits  = 0
// electrons removed by cutting on the layer  = 6321
// electrons removed by cutting on the charge  = 1842
// the inefficeincy of prompt electrons in forward region = 0.0164694
// the inefficeincy of prompt electrons = 0.0057483

// total conversions = 25410
// total conversions in the forward region = 15717
// conversions removed by cutting on the number of pixel hits  = 738
// conversions removed by cutting on the layer  = 2777
// conversions removed by cutting on the charge  = 4740
// the efficeincy of conversions in forward region = 0.478272
// the efficeincy of conversions = 0.324872



runLooper_electron(){
  
  
  //  gStyle->SetTitleX(0.1f);
  //   gStyle->SetTitleW(0.5f);
  
  TChain *Chain = new TChain("Events");
 
 
  //Chain->Add("/store/disk01/yanjuntu/2_2_10/singleElectron_ntuple_v3.root");
 
  Chain->Add("/store/disk02/gutsche/cms2-V1-02-06/SingleElectronFlatPt5To100/merged_ntuple*.root");
 
  TCanvas* c1 = new TCanvas("c1","c1");
  TH1F hist_singleElectron("singleElectron", "singleElectron; Eta", 100,-3,3);
  Chain->Draw("els_p4.eta() >> singleElectron", susy_baseline_2);
  c1->SaveAs("~/sample_plots/conversions/2_2_x/singleElectron_ele_eta_susy.eps");
  float ntotal = hist_singleElectron.Integral();
 
  
  TH1F hist_singleElectron("singleElectron", "singleElectron; Number of valid pixel hits", 16,0,8);
  Chain->Draw("els_valid_pixelhits >> singleElectron",susy_baseline_2 && forward);
  c1->SaveAs("~/sample_plots/conversions/2_2_x/singleElectron_ele_valid_pixelhits_susy.eps");
  float ntotal_forward = hist_singleElectron.Integral();
  float n_removed_pixelhit = hist_singleElectron.GetBinContent(0);
  
  TH1F hist_singleElectron("singleElectron", "singleElectron; DetID", 10,0,5);
  Chain->Draw("els_layer1_det  >> singleElectron", susy_baseline_2 && forward &&"(els_valid_pixelhits==2||els_valid_pixelhits==1)");
  c1->SaveAs("~/sample_plots/conversions/2_2_x/singleElectron_ele_layer1_detid_susy.eps");

  TH1F hist_singleElectron("singleElectron", "singleElectron; Layer of the first valid pixel hit (Barrel)", 10,0,5);
  Chain->Draw("els_layer1_layer  >> singleElectron", susy_baseline_2 && forward &&"(els_valid_pixelhits==2||els_valid_pixelhits==1) && els_layer1_det==1");
  c1->SaveAs("~/sample_plots/conversions/2_2_x/singleElectron_ele_layer1_layer_barrel_susy.eps");
  float n_removed_layer = hist_singleElectron.Integral(hist_singleElectron.FindBin(2),hist_singleElectron.GetNbinsX()+1);
 

  TH1F hist_singleElectron("singleElectron", "singleElectron; Charge of the first valid pixel hit (Disk)", 100,10000,200000);
  Chain->Draw("els_layer1_charge  >> singleElectron", susy_baseline_2 && forward && "(els_valid_pixelhits==2||els_valid_pixelhits==1) &&els_layer1_det==2");
  c1->SaveAs("~/sample_plots/conversions/2_2_x/singleElectron_ele_layer1_charge_disk_susy.eps");
  float n_removed_charge = hist_singleElectron.Integral(hist_singleElectron.FindBin(40000),hist_singleElectron.GetNbinsX()+1);
 
  
  TH2F hist_Electron("Electron", "Electron; Charge of the first valid pixel hit (Disk); Eta", 100, 10000,100000, 300,1.5,2.8);
  Chain->Draw("abs(els_p4.eta()):els_layer1_charge  >> Electron",susy_baseline_2 && forward && "(els_valid_pixelhits==2||els_valid_pixelhits==1) &&els_layer1_det==2", "box");
  c1->SaveAs("~/sample_plots/conversions/2_2_x/singleElectron_ele_layer1_chargeVsEta_disk_susy.eps");
  // TH1F hist_singleElectron("singleElectron", "singleElectron; Number of expected inner layers; Events", 16,0,8);
//   Chain->Draw("els_ninner >> singleElectron", susy_baseline_2 );
//   c1->SaveAs("~/sample_plots/conversions/singleElectron_ele_ninner_susy.eps");
  

  float ineff_forward = (n_removed_layer+n_removed_charge)/ntotal_forward;
  float ineff = (n_removed_pixelhit+n_removed_layer+n_removed_charge)/ntotal;
  std::cout<< "total electrons = "<<ntotal <<std::endl;
  std::cout<< "total electrons in the forward region = "<<ntotal_forward <<std::endl;
  std::cout<< "electrons removed by cutting on the number of pixel hits  = "<< n_removed_pixelhit<<std::endl;
  std::cout<< "electrons removed by cutting on the layer  = "<< n_removed_layer<<std::endl;
  std::cout<< "electrons removed by cutting on the charge  = "<< n_removed_charge<<std::endl;
  std::cout<< "the inefficeincy of prompt electrons in forward region = "<<ineff_forward <<std::endl;
  std::cout<< "the inefficeincy of prompt electrons = "<<ineff <<std::endl;
}



runLooper_gamma(){
  TChain *Chain = new TChain("Events");
 
  
  //Chain->Add("/store/disk01/yanjuntu/2_2_10/singleGamma_ntuple_v3.root");
 
  Chain->Add("/store/disk02/yanjuntu/cms2-V01-02-06/RelValSingleGammaFlatPt20To100/preprocessing/ntuple*.root");
 
  TCanvas* c1 = new TCanvas("c1","c1");
  TH1F hist_singleGamma("singleGamma", "singleGamma; Eta", 100,-3,3);
  Chain->Draw("els_p4.eta() >> singleGamma", susy_baseline_2);    
  c1->SaveAs("~/sample_plots/conversions/2_2_x/singleGamma_ele_eta_susy.eps");
  float ntotal = hist_singleGamma.Integral();
   
  TH1F hist_singleGamma("singleGamma", "singleGamma; Number of valid pixel hits", 16,0,8);
  Chain->Draw("els_valid_pixelhits >> singleGamma",susy_baseline_2 && forward);
  c1->SaveAs("~/sample_plots/conversions/2_2_x/singleGamma_ele_valid_pixelhits_susy.eps");
  float ntotal_forward = hist_singleGamma.Integral();
  float n_removed_pixelhit = hist_singleGamma.GetBinContent(hist_singleGamma.FindBin(0));

  TH1F hist_singleGamma("singleGamma", "singleGamma; DetID", 10,0,5);
  Chain->Draw("els_layer1_det  >> singleGamma", susy_baseline_2 && forward &&"(els_valid_pixelhits==2||els_valid_pixelhits==1)");
  c1->SaveAs("~/sample_plots/conversions/2_2_x/singleGamma_ele_layer1_detid_susy.eps");

  TH1F hist_singleGamma("singleGamma", "singleGamma; Layer of the first valid pixel hit (Barrel)", 10,0,5);
  Chain->Draw("els_layer1_layer  >> singleGamma", susy_baseline_2 && forward &&"(els_valid_pixelhits==2||els_valid_pixelhits==1) && els_layer1_det==1");
  c1->SaveAs("~/sample_plots/conversions/2_2_x/singleGamma_ele_layer1_layer_barrel_susy.eps");
  float n_removed_layer = hist_singleGamma.Integral(hist_singleGamma.FindBin(2),hist_singleGamma.GetNbinsX()+1);

  TH1F hist_singleGamma("singleGamma", "singleGamma; Charge of the first valid pixel hit (Disk)", 100,10000,200000);
  Chain->Draw("els_layer1_charge  >> singleGamma", susy_baseline_2 && forward && "(els_valid_pixelhits==2||els_valid_pixelhits==1) &&els_layer1_det==2");
  c1->SaveAs("~/sample_plots/conversions/2_2_x/singleGamma_ele_layer1_charge_disk_susy.eps");
  float n_removed_charge = hist_singleGamma.Integral(hist_singleGamma.FindBin(40000),hist_singleGamma.GetNbinsX()+1);

  TH2F hist_Gamma("Gamma", "Gamma; Charge of the first valid pixel hit (Disk); Eta", 100, 10000,100000, 300,1.5,2.8);
  Chain->Draw("abs(els_p4.eta()):els_layer1_charge  >> Gamma",susy_baseline_2 && forward && "(els_valid_pixelhits==2||els_valid_pixelhits==1) &&els_layer1_det==2", "box");
  c1->SaveAs("~/sample_plots/conversions/2_2_x/singleGamma_ele_layer1_chargeVsEta_disk_susy.eps");
  // TH1F hist_singleGamma("singleGamma", "singleGamma; Number of expected inner layers; Events", 16,0,8);
  //   Chain->Draw("els_ninner >> singleGamma", susy_baseline_2 );
  //   c1->SaveAs("~/sample_plots/conversions/singleGamma_ele_ninner_susy.eps");

  float eff_forward = (n_removed_pixelhit+n_removed_layer+n_removed_charge)/ntotal_forward;
  float eff = (n_removed_pixelhit+n_removed_layer+n_removed_charge)/ntotal;
  std::cout<< "total conversions = "<<ntotal <<std::endl;
  std::cout<< "total conversions in the forward region = "<<ntotal_forward <<std::endl;
  std::cout<< "conversions removed by cutting on the number of pixel hits  = "<< n_removed_pixelhit<<std::endl;
  std::cout<< "conversions removed by cutting on the layer  = "<< n_removed_layer<<std::endl;
  std::cout<< "conversions removed by cutting on the charge  = "<< n_removed_charge<<std::endl;
  std::cout<< "the efficeincy of conversions in forward region = "<<eff_forward <<std::endl;
  std::cout<< "the efficeincy of conversions = "<<eff <<std::endl;

}

