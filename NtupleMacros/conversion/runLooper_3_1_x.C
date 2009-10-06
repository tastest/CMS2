//TCut tightid   = "els_tightId22XMinMatteo ==1"; //used in 2_x 
TCut tightid   = "els_egamma_tightId ==1";  //move to 3_x
TCut pt20      = "els_p4.pt() > =20";
TCut forward   = "abs(els_p4.eta()) > 1.47";
TCut central   = "abs(els_p4.eta()) < 1.47";
TCut d0corr250 = "abs(els_d0corr) < 0.025";
TCut d0corr200 = "abs(els_d0corr) < 0.020";
TCut conv_dcot = "abs(els_conv_dcot)<0.02"; 
TCut conv_dist = "abs(els_conv_dist)<0.02";
TCut conversion_1 =  conv_dcot && conv_dist;
TCut notmuon   = "els_closestMuon == -1";
//TCut susy_iso  = "(els_tkIso03+max(0., (els_ecalIso03-2.)) + els_hcalIso03)/max(els_p4.pt(),20.) < 0.1";  //ecal is buggy in 3_1_0
TCut susy_iso  = "(els_tkIso03 + els_hcalIso03)/max(els_p4.pt(),20.) < 0.1";
//TCut ww_iso    = "els_p4.pt()/(els_p4.pt()+els_tkIso03+els_ecalIso+els_hcalIso03) >0.92";
TCut ww_iso    = "els_p4.pt()/(els_p4.pt()+els_tkIso03+els_hcalIso03) >0.92";

TCut susy_eta  = "TMath::Abs(els_p4.eta()) <= 2.4";
TCut ww_baseline_0 = pt20 && tightid && notmuon && ww_iso;
TCut ww_baseline_1 = ww_baseline_0 && d0corr250;
//TCut ww_baseline_2 = ww_baseline_1 && !conversion_1;
TCut ww_baseline_2 = ww_baseline_1 ;

TCut susy_baseline_0 = susy_eta&& pt20 && tightid && notmuon && susy_iso;
TCut susy_baseline_1 = susy_baseline_0 && d0corr200;
//TCut susy_baseline_2 = susy_baseline_1 && !conversion_1;
TCut susy_baseline_2 = susy_baseline_1 ;  

// total electrons = 1.33213e+06
// total electrons in the forward region = 437081
// electrons removed by cutting on the number of pixel hits  = 0
// electrons removed by cutting on the layer  = 3383
// electrons removed by cutting on the charge  = 1315
// the inefficeincy of prompt electrons in forward region = 0.0107486
// the inefficeincy of prompt electrons = 0.00352668
// electrons removed by cutting ninner layers (loose)  = 40541
// electrons removed by cutting ninner layers (tight)  = 2486
// method 2 loose cut: the inefficeincy of prompt electrons = 0.0304332
// method 2 tight cut: the inefficeincy of prompt electrons = 0.00186618


// total conversions = 28470
// total conversions in the forward region = 18243
// conversions removed by cutting on the number of pixel hits  = 4568
// conversions removed by cutting on the layer  = 1882
// conversions removed by cutting on the charge  = 5421
// the efficeincy of conversions in forward region = 0.650715
// the efficeincy of conversions = 0.416965
// conversions removed by cutting ninner layers (loose)  = 22272
// conversions removed by cutting ninner layers (tight)  = 12017
// method 2 loose cut: the inefficeincy of conversions = 0.782297
// method 2 tight cut: the inefficeincy of conversions = 0.422093

//turn off dcot and dist cut

// total electrons = 1.33721e+06
// total electrons in the forward region = 463080
// electrons removed by cutting on the number of pixel hits  = 0
// electrons removed by cutting on the layer  = 3438
// electrons removed by cutting on the charge  = 1746
// the inefficeincy of prompt electrons in forward region = 0.0111946
// the inefficeincy of prompt electrons = 0.00387673
// electrons removed by cutting ninner layers (loose)  = 42212
// electrons removed by cutting ninner layers (tight)  = 2810
// method 2 loose cut: the inefficeincy of prompt electrons = 0.0315673
// method 2 tight cut: the inefficeincy of prompt electrons = 0.00210139


// total conversions = 48004
// total conversions in the forward region = 29798
// conversions removed by cutting on the number of pixel hits  = 6480
// conversions removed by cutting on the layer  = 2409
// conversions removed by cutting on the charge  = 9177
// the efficeincy of conversions in forward region = 0.606282
// the efficeincy of conversions = 0.376344
// conversions removed by cutting ninner layers (loose)  = 34979
// conversions removed by cutting ninner layers (tight)  = 17672
// method 2 loose cut: the inefficeincy of conversions = 0.728668
// method 2 tight cut: the inefficeincy of conversions = 0.368136


runLooper_electron(){
  
  
  //  gStyle->SetTitleX(0.1f);
  //   gStyle->SetTitleW(0.5f);
  
  TChain *Chain = new TChain("Events");
 
  // Chain->Add("/store/disk01/yanjuntu/SingleElectron83055667ceda31da68d5aa532f72919d/preprocessing/ntuple*.root");  //tas03
  Chain->Add("/hadoop/cms/store/user/yanjuntu/conversion_sample_cms2-V01-03-01/SingleElectron83055667ceda31da68d5aa532f72919d/preprocessing/ntuple*.root");   //uaf6
  //Chain->Add("/store/disk01/yanjuntu/2_2_10/singleElectron_ntuple_v3.root");
 
  
 
  TCanvas* c1 = new TCanvas("c1","c1");
  TH1F hist_singleElectron("singleElectron", "; Eta", 100,-3,3);
  Chain->Draw("els_p4.eta() >> singleElectron", susy_baseline_2);
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleElectron_ele_eta_susy.eps");
  float ntotal = hist_singleElectron.Integral();
 
  
  TH1F hist_singleElectron("singleElectron", "; Number of valid pixel hits", 16,0,8);
  Chain->Draw("els_valid_pixelhits >> singleElectron",susy_baseline_2 && forward);
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleElectron_ele_valid_pixelhits_susy.eps");
  float ntotal_forward = hist_singleElectron.Integral();
  float n_removed_pixelhit = hist_singleElectron.GetBinContent(0);
  
  TH1F hist_singleElectron("singleElectron", "; DetID", 10,0,5);
  Chain->Draw("els_layer1_det  >> singleElectron", susy_baseline_2 && forward &&"(els_valid_pixelhits==2||els_valid_pixelhits==1)");
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleElectron_ele_layer1_detid_susy.eps");

  TH1F hist_singleElectron("singleElectron", "; Layer of the first valid pixel hit (Barrel)", 10,0,5);
  Chain->Draw("els_layer1_layer  >> singleElectron", susy_baseline_2 && forward &&"(els_valid_pixelhits==2||els_valid_pixelhits==1) && els_layer1_det==1");
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleElectron_ele_layer1_layer_barrel_susy.eps");
  float n_removed_layer = hist_singleElectron.Integral(hist_singleElectron.FindBin(2),hist_singleElectron.GetNbinsX()+1);
 

  TH1F hist_singleElectron("singleElectron", "; Charge of the first valid pixel hit (Disk)", 100,10000,200000);
  Chain->Draw("els_layer1_charge  >> singleElectron", susy_baseline_2 && forward && "(els_valid_pixelhits==2||els_valid_pixelhits==1) &&els_layer1_det==2");
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleElectron_ele_layer1_charge_disk_susy.eps");
  float n_removed_charge = hist_singleElectron.Integral(hist_singleElectron.FindBin(40000),hist_singleElectron.GetNbinsX()+1);
 
  
  TH2F hist_Electron("Electron", "; Charge of the first valid pixel hit (Disk); Eta", 100, 10000,100000, 300,1.5,2.8);
  Chain->Draw("abs(els_p4.eta()):els_layer1_charge  >> Electron",susy_baseline_2 && forward && "(els_valid_pixelhits==2||els_valid_pixelhits==1) &&els_layer1_det==2", "box");
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleElectron_ele_layer1_chargeVsEta_disk_susy.eps");
  
  TH1F hist_singleElectron("singleElectron", "; Number of expected inner layers", 16,0,8);
  Chain->Draw("els_n_inner_layers >> singleElectron", susy_baseline_2 );
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleElectron_ele_ninner_susy.eps");
  float n_removed_m2_0 = hist_singleElectron.Integral(hist_singleElectron.FindBin(1),hist_singleElectron.GetNbinsX()+1);
  float n_removed_m2_1 = hist_singleElectron.Integral(hist_singleElectron.FindBin(2),hist_singleElectron.GetNbinsX()+1);
 
  float ineff_forward = (n_removed_pixelhit+n_removed_layer+n_removed_charge)/ntotal_forward;
  float ineff = (n_removed_pixelhit+n_removed_layer+n_removed_charge)/ntotal;
  std::cout<< "total electrons = "<<ntotal <<std::endl;
  std::cout<< "total electrons in the forward region = "<<ntotal_forward <<std::endl;
  std::cout<< "electrons removed by cutting on the number of pixel hits  = "<< n_removed_pixelhit<<std::endl;
  std::cout<< "electrons removed by cutting on the layer  = "<< n_removed_layer<<std::endl;
  std::cout<< "electrons removed by cutting on the charge  = "<< n_removed_charge<<std::endl;
  std::cout<< "the inefficeincy of prompt electrons in forward region = "<<ineff_forward <<std::endl;
  std::cout<< "the inefficeincy of prompt electrons = "<<ineff <<std::endl;
  
  float ineff_m2_0 = n_removed_m2_0/ntotal;
  float ineff_m2_1 = n_removed_m2_1/ntotal;
  std::cout<< "electrons removed by cutting ninner layers (loose)  = "<<n_removed_m2_0 <<std::endl;
  std::cout<< "electrons removed by cutting ninner layers (tight)  = "<<n_removed_m2_1 <<std::endl;
  std::cout<< "method 2 loose cut: the inefficeincy of prompt electrons = "<<ineff_m2_0 <<std::endl;
  std::cout<< "method 2 tight cut: the inefficeincy of prompt electrons = "<<ineff_m2_1 <<std::endl;
}



runLooper_gamma(){
  TChain *Chain = new TChain("Events");
 
  //Chain->Add("/store/disk01/yanjuntu/SingleGamma4016bf049eded18e8e9dac4c09773936/preprocessing/ntuple*.root"); //tas03
  Chain->Add("/hadoop/cms/store/user/yanjuntu/conversion_sample_cms2-V01-03-01/SingleGamma4016bf049eded18e8e9dac4c09773936/preprocessing/ntuple*.root"); //uaf6
  //Chain->Add("/store/disk01/yanjuntu/2_2_10/singleGamma_ntuple_v3.root");
 
  
 
  TCanvas* c1 = new TCanvas("c1","c1");
  TH1F hist_singleGamma("singleGamma", "; Eta", 100,-3,3);
  Chain->Draw("els_p4.eta() >> singleGamma", susy_baseline_2);    
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleGamma_ele_eta_susy.eps");
  float ntotal = hist_singleGamma.Integral();
   
  TH1F hist_singleGamma("singleGamma", "; Number of valid pixel hits", 16,0,8);
  Chain->Draw("els_valid_pixelhits >> singleGamma",susy_baseline_2 && forward);
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleGamma_ele_valid_pixelhits_susy.eps");
  float ntotal_forward = hist_singleGamma.Integral();
  float n_removed_pixelhit = hist_singleGamma.GetBinContent(hist_singleGamma.FindBin(0));

  TH1F hist_singleGamma("singleGamma", "; DetID", 10,0,5);
  Chain->Draw("els_layer1_det  >> singleGamma", susy_baseline_2 && forward &&"(els_valid_pixelhits==2||els_valid_pixelhits==1)");
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleGamma_ele_layer1_detid_susy.eps");

  TH1F hist_singleGamma("singleGamma", "; Layer of the first valid pixel hit (Barrel)", 10,0,5);
  Chain->Draw("els_layer1_layer  >> singleGamma", susy_baseline_2 && forward &&"(els_valid_pixelhits==2||els_valid_pixelhits==1) && els_layer1_det==1");
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleGamma_ele_layer1_layer_barrel_susy.eps");
  float n_removed_layer = hist_singleGamma.Integral(hist_singleGamma.FindBin(2),hist_singleGamma.GetNbinsX()+1);

  TH1F hist_singleGamma("singleGamma", "; Charge of the first valid pixel hit (Disk)", 100,10000,200000);
  Chain->Draw("els_layer1_charge  >> singleGamma", susy_baseline_2 && forward && "(els_valid_pixelhits==2||els_valid_pixelhits==1) &&els_layer1_det==2");
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleGamma_ele_layer1_charge_disk_susy.eps");
  float n_removed_charge = hist_singleGamma.Integral(hist_singleGamma.FindBin(40000),hist_singleGamma.GetNbinsX()+1);

  TH2F hist_Gamma("Gamma", "; Charge of the first valid pixel hit (Disk); Eta", 100, 10000,100000, 300,1.5,2.8);
  Chain->Draw("abs(els_p4.eta()):els_layer1_charge  >> Gamma",susy_baseline_2 && forward && "(els_valid_pixelhits==2||els_valid_pixelhits==1) &&els_layer1_det==2", "box");
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleGamma_ele_layer1_chargeVsEta_disk_susy.eps");
  
  TH1F hist_singleGamma("singleGamma", "; Number of expected inner layers", 16,0,8);
  Chain->Draw("els_n_inner_layers >> singleGamma", susy_baseline_2 );
  c1->SaveAs("~/sample_plots/conversions/3_1_x/singleGamma_ele_ninner_susy.eps");
  float n_removed_m2_0 = hist_singleGamma.Integral(hist_singleGamma.FindBin(1),hist_singleGamma.GetNbinsX()+1);
  float n_removed_m2_1 = hist_singleGamma.Integral(hist_singleGamma.FindBin(2),hist_singleGamma.GetNbinsX()+1);
  
  float eff_forward = (n_removed_pixelhit+n_removed_layer+n_removed_charge)/ntotal_forward;
  float eff = (n_removed_pixelhit+n_removed_layer+n_removed_charge)/ntotal;
  std::cout<< "total conversions = "<<ntotal <<std::endl;
  std::cout<< "total conversions in the forward region = "<<ntotal_forward <<std::endl;
  std::cout<< "conversions removed by cutting on the number of pixel hits  = "<< n_removed_pixelhit<<std::endl;
  std::cout<< "conversions removed by cutting on the layer  = "<< n_removed_layer<<std::endl;
  std::cout<< "conversions removed by cutting on the charge  = "<< n_removed_charge<<std::endl;
  std::cout<< "the efficeincy of conversions in forward region = "<<eff_forward <<std::endl;
  std::cout<< "the efficeincy of conversions = "<<eff <<std::endl;

  
  float ineff_m2_0 = n_removed_m2_0/ntotal;
  float ineff_m2_1 = n_removed_m2_1/ntotal;
  std::cout<< "conversions removed by cutting ninner layers (loose)  = "<<n_removed_m2_0 <<std::endl;
  std::cout<< "conversions removed by cutting ninner layers (tight)  = "<<n_removed_m2_1 <<std::endl;
  std::cout<< "method 2 loose cut: the inefficeincy of conversions = "<<ineff_m2_0 <<std::endl;
  std::cout<< "method 2 tight cut: the inefficeincy of conversions = "<<ineff_m2_1 <<std::endl;

}

