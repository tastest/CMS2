{
  TCut goodIsolatedElectron("els_tightId22XMinMatteo==1&&els_closestMuon==-1&&TMath::Abs(els_d0corr)<=0.025&&els_p4.pt()/(els_p4.pt()+els_pat_trackIso+1e-5)>0.92");
  TCut noConvElectron("TMath::Abs(els_conv_dist) >= 0.02 || TMath::Abs(els_conv_dcot) >= 0.02");
  TCut convElectronYanjun("els_trk_p4.eta()>1.5 && (els_valid_pixelhits==0 || ((els_valid_pixelhits==1||els_valid_pixelhits==2)&&((els_layer1_det==1&&els_layer1_layer>1)||(els_layer1_det==2&&els_layer1_charge>40000))))");
  TCut olisGoodElectron = goodIsolatedElectron && noConvElectron && !convElectronYanjun;

  TCut eleTrkInConsistent("els_trkidx >= 0 && els_charge != trks_charge[els_trkidx]");

  TCut eleTrueMatch("els_mc_id!=-999");

  TCut corCharge("(els_charge == -1 && els_mc_id == 11) || (els_charge == 1 && els_mc_id == -11)");

  TCut barrelElectron("TMath::Abs(els_mc_p4.eta()) <= 1.479");

  TCut forwardCor = goodIsolatedElectron && noConvElectron && !convElectronYanjun && !eleTrkInConsistent && eleTrueMatch && corCharge && !barrelElectron;
  TCut forwardIncor = goodIsolatedElectron && noConvElectron && !convElectronYanjun && !eleTrkInConsistent && eleTrueMatch && !corCharge && !barrelElectron;



}
