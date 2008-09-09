{
   // final numbers with projected met and flags
   TCut final_emu_new2  = "((ww_ltgoodmuiso&&ww_llgoodeliso&&hyp_type==1)||(ww_ltgoodeliso&&ww_llgoodmuiso&&hyp_type==2))&&ww_pass4met&&ww_passzveto&&ww_passaddzveto&&hyp_ll_p4.Pt()>20&&ww_oppsign&&hyp_njets==0";
   TCut final_ee_new2  = "ww_ltgoodeliso&&ww_llgoodeliso&&hyp_type==3&&ww_pass4met&&ww_passzveto&&ww_passaddzveto&&hyp_ll_p4.Pt()>20&&ww_oppsign&&hyp_njets==0";
   TCut final_mm_new2  = "ww_ltgoodmuiso&&ww_llgoodmuiso&&hyp_type==0&&ww_pass4met&&ww_passzveto&&ww_passaddzveto&&hyp_ll_p4.Pt()>20&&ww_oppsign&&hyp_njets==0";

   // divide hypotheses
   TCut ee (" abs(hyp_lt_id)==11 && abs(hyp_ll_id)==11");
   TCut em (" abs(hyp_lt_id)==11 && abs(hyp_ll_id)==13");
   TCut me (" abs(hyp_lt_id)==13 && abs(hyp_ll_id)==11");
   TCut mm (" abs(hyp_lt_id)==13 && abs(hyp_ll_id)==13");
   TCut dyel = "genps_id[8] == 11 || genps_id[8] == -11";
   TCut dymu = "genps_id[8] == 13 || genps_id[8] == -13";
   TCut dyta = "genps_id[8] == 15 || genps_id[8] == -15";
   
   //Muon cuts
   TCut muo_lt_id (" mus_validHits[hyp_lt_index] > 6 && mus_gfit_chi2[hyp_lt_index]/mus_gfit_ndof[hyp_lt_index] <= 5. ");
   TCut muo_lt_d0 (" abs(mus_d0[hyp_lt_index]) <= 0.25 "); 
   TCut muo_lt_iso (" mus_p4[hyp_lt_index].pt()/(mus_p4[hyp_lt_index].pt() + (mus_iso03_sumPt[hyp_lt_index]+mus_iso03_emEt[hyp_lt_index]+mus_iso03_hadEt[hyp_lt_index]) ) >= 0.92 ") ;
   TCut gd_lt_muo = muo_lt_id + muo_lt_d0 + muo_lt_iso;
   
   TCut muo_ll_id (" mus_validHits[hyp_ll_index] > 6 && mus_gfit_chi2[hyp_ll_index]/mus_gfit_ndof[hyp_ll_index] <= 5. ");
   TCut muo_ll_d0 (" abs(mus_d0[hyp_ll_index]) <= 0.25 "); 
   TCut muo_ll_iso (" mus_p4[hyp_ll_index].pt()/(mus_p4[hyp_ll_index].pt() + (mus_iso03_sumPt[hyp_ll_index]+mus_iso03_emEt[hyp_ll_index]+mus_iso03_hadEt[hyp_ll_index]) ) >= 0.92 ") ;
   TCut gd_ll_muo = muo_ll_id + muo_ll_d0 + muo_ll_iso;

   //Electron cuts
   TCut ele_lt_id ( "els_tightId[hyp_lt_index]     ==  1" ) ;
   TCut ele_lt_d0 (" abs(els_d0[hyp_lt_index]) <= 0.025 "); 
   TCut ele_lt_iso (" hyp_lt_p4.pt()/(hyp_lt_p4.pt() + els_tkIso[hyp_lt_index]) > 0.92 && els_closestMuon[hyp_lt_index]==-1 ");
   TCut gd_lt_ele = ele_lt_id + ele_lt_d0 + ele_lt_iso;

   TCut ele_ll_id ( "els_tightId[hyp_ll_index]     ==  1" ) ;
   TCut ele_ll_d0 (" abs(els_d0[hyp_ll_index]) <= 0.025 "); 
   TCut ele_ll_iso (" hyp_ll_p4.pt()/(hyp_ll_p4.pt() + els_tkIso[hyp_ll_index]) > 0.92 && els_closestMuon[hyp_ll_index]==-1");
   TCut gd_ll_ele = ele_ll_id + ele_ll_d0 + ele_ll_iso;

   //MET cut
   TCut metcut_emu = (em||me) && "hyp_met>20";
   TCut metcut_ll = (ee||mm) && "hyp_met>30 && abs(hyp_p4.mass()-91)>15 && (hyp_met/hyp_p4.pt()>0.6 || acos(cos(hyp_metPhi-hyp_p4.Phi()-3.1416))>0.25)";
   TCut metcut = metcut_emu || metcut_ll;
   
   // first set of cuts
   TCut final_emu  = ((gd_lt_muo+gd_ll_ele+"hyp_type==1")||(gd_lt_ele+gd_ll_muo+"hyp_type==2"))+metcut+"hyp_ll_p4.Pt()>20&&hyp_lt_id*hyp_ll_id==-143&&hyp_njets==0";
   TCut final_emu_new  = "((ww_ltgoodmuiso&&ww_llgoodeliso&&hyp_type==1)||(ww_ltgoodeliso&&ww_llgoodmuiso&&hyp_type==2))&&ww_pass2met&&ww_passzveto&&hyp_ll_p4.Pt()>20&&ww_oppsign&&hyp_njets==0";
   

}
