#include <iostream>

void xsecCalc(double obs, double expS, double expSErr, double bg, double bgErr){
  double obsS = obs - bg;
  double obsSErr = sqrt( obs + bgErr*bgErr);

  double sigma = 157.5 * (obs - bg)/ expS;

  double sigmaErrStatFrac = sqrt(obs)/obsS;
  double sigmaErrSystFrac = sqrt(bgErr*bgErr/obsS/obsS + expSErr*expSErr/expS/expS);
  double sigmaErrFrac = sqrt(sigmaErrStatFrac*sigmaErrStatFrac + sigmaErrSystFrac*sigmaErrSystFrac);
  
  double sigmaErrFracLum = 0.11;

  double sigmaErrStat = sigma* sigmaErrStatFrac;
  double sigmaErrSyst = sigma* sigmaErrSystFrac;
  double sigmaErr = sigma*sigmaErrFrac;
  double sigmaErrLum = sigma*sigmaErrFracLum;

  double sigmaErrWLum = sigma * sqrt( sigmaErrFrac*sigmaErrFrac + sigmaErrFracLum*sigmaErrFracLum);

  std::cout<<"sigma = "
	   << sigma << " \\pm "<< sigmaErrStat << " (stat) \\pm " << sigmaErrSyst 
	   << " (syst) \\pm "<<sigmaErrLum << " (lum)"<<std::endl
	   <<"\t = "<<sigma << " \\pm " << sigmaErr << " (syst+stat) \\pm "<< sigmaErrLum << " (lum)"
	   <<"\n\t = "<<sigma<<" ( 1 \\pm "<<sigmaErrStat/sigma<<" (stat) \\pm "<<sigmaErrSyst/sigma 
	   <<" (syst) \\pm "<<sigmaErrFracLum<<" (lum) )"
	   <<"\n\t = "<<sigma<<" ( 1 \\pm "
	   <<sigmaErr/sigma << " (syst+stat) \\pm "
	   <<sigmaErrFracLum <<" )" <<std::endl;

}


struct TTxsecStruct {
  double tt_exp; double tt_exp_me;

  double dytt_exp; double dytt_exp_me;
  double vv_exp; double vv_exp_me;
  double tw_exp; double tw_exp_me;

  double data;
  double sr; double sr_stat; double sr_syst;

  double qcd_exp; double qcd_stat; double qcd_syst;
  double wjraw_exp; double wjraw_exp_stat;
  double wjf_systFrac;

  double dy_exp; double dy_stat; double dy_syst; double dy_e;

  double SF_e;//
  double tt_AE_eRel;

};

void xsecCalc_inStruct(TTxsecStruct& tt){

  double mcbg_exp = tt.dytt_exp + tt.vv_exp + tt.tw_exp;
  double mcbg_e   = sqrt(0.25*(tt.dytt_exp*tt.dytt_exp 
				     + tt.vv_exp*tt.vv_exp
				     + tt.tw_exp*tt.tw_exp)
			       + tt.dytt_exp_me*tt.dytt_exp_me
			       + tt.vv_exp_me*tt.vv_exp_me
			       + tt.tw_exp_me*tt.tw_exp_me); 

  double spill_exp = tt.sr * tt.data;
  double spill_stat = sqrt( tt.data* tt.sr*tt.sr );
  double spill_syst = sqrt(  (tt.data * tt.sr_stat)* (tt.data * tt.sr_stat)
				   + (tt.sr_syst * tt.data)*(tt.sr_syst * tt.data)
				   );
  double spill_e = sqrt(spill_stat*spill_stat  +  spill_syst*spill_syst );

  double wjf_exp = tt.wjraw_exp - 2.* tt.qcd_exp - spill_exp;
  if (wjf_exp< 0) std::cout<<"Warning: negative wj estimate: "<<wjf_exp<<std::endl;
  double wjf_syst = tt.wjf_systFrac*wjf_exp;

  double fake_exp  = wjf_exp + tt.qcd_exp;
  double fake_stat = sqrt( tt.wjraw_exp_stat*tt.wjraw_exp_stat 
				 + tt.qcd_stat*tt.qcd_stat
				 + spill_stat*spill_stat);
  double fake_syst = sqrt(wjf_syst*wjf_syst
				+ tt.qcd_syst*tt.qcd_syst
				+ spill_syst*spill_syst );

  double fake_e = sqrt(fake_stat * fake_stat 
			     + fake_syst*fake_syst );
//   std::cout<<" wjf_exp "<<wjf_exp<<" wjfSF "<<tt.wjf_systFrac<<" wjfsyst "<<wjf_syst<<" qcdsyst "<< tt.qcd_syst<<" ssyst "<<spill_syst
// 	   <<" fstat "<<fake_stat<<" fsyst "<<fake_syst<<std::endl;

  double dy_e = sqrt(tt.dy_stat*tt.dy_stat + tt.dy_syst*tt.dy_syst);

  double bg_exp = mcbg_exp + fake_exp + tt.dy_exp;
  double bg_e = sqrt(mcbg_e*mcbg_e + fake_e*fake_e + dy_e*dy_e );

  //this is it:
  std::cout<<"Data: "<<tt.data
	   <<"\n\t   MCBG: "<<mcbg_exp <<" +/- " << mcbg_e
	   <<"\n\t   Fake: "<<fake_exp <<" +/- " << fake_stat <<" +/- "<< fake_syst
	   <<"\n\t     DY: "<<tt.dy_exp <<" +/- " << tt.dy_stat <<" +/- "<< tt.dy_syst
	   <<"\n\t All BG: "<<bg_exp <<" +/- " << bg_e
	   <<"\n\t tt exp: "<<tt.tt_exp*tt.SF_e<<" +/- "<<tt.tt_exp*tt.SF_e*tt.tt_AE_eRel
	   <<"\tS/B: "<<tt.tt_exp*tt.SF_e/bg_exp<<std::endl;
  xsecCalc(tt.data, tt.tt_exp*tt.SF_e, tt.tt_exp*tt.SF_e*tt.tt_AE_eRel, bg_exp, bg_e );
  

}

void xsecCalc_35pb_pass5(){
// Sample  & ee & $\mu\mu$ & e$\mu$ & all \\ \hline
//     ttdil  & 26.552 $\pm$  0.314 & 29.600 $\pm$  0.332   & 55.057 $\pm$  0.452    & 111.208 $\pm$  0.643 \\  \hline
//     ttotr  &  0.818 $\pm$  0.055 &  0.164 $\pm$  0.025   &  1.167 $\pm$  0.066    &  2.148 $\pm$  0.089 \\
//     wjets  &  0.763 $\pm$  0.289 &  0.000 $\pm$  0.000   &  0.763 $\pm$  0.289    &  1.527 $\pm$  0.408 \\
//  DYtautau  &  0.983 $\pm$  0.311 &  1.207 $\pm$  0.342   &  0.787 $\pm$  0.278    &  2.976 $\pm$  0.539 \\
//        VV  &  0.934 $\pm$  0.039 &  0.918 $\pm$  0.039   &  1.050 $\pm$  0.041    &  2.903 $\pm$  0.069 \\
//        tw  &  0.873 $\pm$  0.026 &  0.940 $\pm$  0.027   &  1.804 $\pm$  0.038    &  3.617 $\pm$  0.054 \\
//       QCD  &  0.000 $\pm$  0.000 &  0.000 $\pm$  0.000   &  0.000 $\pm$  0.000    &  0.000 $\pm$  0.000 \\
//    DYeemm  &288.997 $\pm$  5.307 & 331.322 $\pm$ 5.687   &  0.295 $\pm$  0.170    & 620.614 $\pm$  7.780 \\  \hline
//  Total MC  &319.920 $\pm$  5.333 & 364.149 $\pm$ 5.707   & 60.923 $\pm$  0.634    & 744.993 $\pm$  7.837 \\  \hline
//      Data  &383                  & 532                   & 57                     & 972 \\  \hline
// \textrm{$>=$ 2 Jet Bin:} \\
//               Sample & ee & $\mu\mu$ & e$\mu$ & all \\ \hline
//     WJets Prediction  & 27.0031 $\pm$ 2.7892 & 13.6705 $\pm$ 2.1510  & 5.2129 $\pm$ 1.1978  & 45.8865 $\pm$ 3.7203 \\
//       QCD Prediction  & 2.4330  $\pm$ 0.4241 & 0.2810  $\pm$ 0.1750  & 0.7029 $\pm$ 0.2475  & 3.4169 $\pm$ 0.5213 \\ \hline
//        FR Prediction  & 24.5701 $\pm$ 2.8212 & 13.3895 $\pm$ 2.1581  & 4.5100 $\pm$ 1.2231  & 42.4696 $\pm$ 3.7567 \\
// Data Yield in Z mass  & 320                  & 437                   & 57                   & 814 \\
//    Spillage Fraction  & 0.0768 $\pm$ 0.0098 & 0.0306 $\pm$ 0.0052 & 0.0537 $\pm$ 0.0052 &  -  &  \\ \hline

  std::cout<<"Doing JPT/TC"<<std::endl;
  std::cout<<"EMU (no MET)"<<std::endl;
  TTxsecStruct tem;
  tem.tt_exp   = 55.057; tem.tt_exp_me    =    0.452;
  tem.dytt_exp = 0.787 ;  tem.dytt_exp_me =    0.278;
  tem.vv_exp   = 1.050 ;  tem.vv_exp_me   =    0.041;
  tem.tw_exp   = 1.804 ;  tem.tw_exp_me   =    0.038;

  tem.data = 57;

  tem.sr  = 0.0537; tem.sr_stat = 0.0052; tem.sr_syst = (0.0537 - 0.0307)*0.5;


  tem.qcd_exp   = 0.7029; tem.qcd_stat       = 0.2475; tem.qcd_syst = tem.qcd_exp;
  tem.wjraw_exp = 5.2129; tem.wjraw_exp_stat = 1.1978;
  tem.wjf_systFrac = 0.75;

  tem.dy_exp = 0; tem.dy_stat = 0; tem.dy_syst = 0;

  tem.SF_e = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  tem.tt_AE_eRel = 0.062;

  xsecCalc_inStruct(tem);

  std::cout<<"Doing JPT/TC"<<std::endl;
  std::cout<<"EMU (w MET)"<<std::endl;
  TTxsecStruct teme;
  teme.tt_exp   = 50.542; teme.tt_exp_me    =    0.433;
  teme.dytt_exp = 0.787 ;  teme.dytt_exp_me =    0.278;
  teme.vv_exp   = 0.890 ;  teme.vv_exp_me   =    0.038;
  teme.tw_exp   = 1.645 ;  teme.tw_exp_me   =    0.036;

  teme.data = 49;

  teme.sr  = 0.0537; teme.sr_stat = 0.0052; teme.sr_syst = (0.0537 - 0.0307)*0.5;


  teme.qcd_exp   = 0.4480; teme.qcd_stat       = 0.2094; teme.qcd_syst = teme.qcd_exp;
  teme.wjraw_exp = 4.1669 ; teme.wjraw_exp_stat = 1.0667;
  teme.wjf_systFrac = 0.75;

  teme.dy_exp = 0; teme.dy_stat = 0; teme.dy_syst = 0;

  teme.SF_e = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  teme.tt_AE_eRel = 0.062;
 
  xsecCalc_inStruct(teme);
 
  std::cout<<"Doing JPT/TC"<<std::endl;
  std::cout<<"EE (w MET)"<<std::endl;
  TTxsecStruct tee;
  tee.tt_exp   = 16.849; tee.tt_exp_me    =   0.250;
  tee.dytt_exp = 0.492;  tee.dytt_exp_me =   0.220;
  tee.vv_exp   = 0.320;  tee.vv_exp_me   =   0.023;
  tee.tw_exp   = 0.567;  tee.tw_exp_me   =   0.021;

  tee.data = 21;

  tee.sr  = 0.0768; tee.sr_stat = 0.0098; tee.sr_syst = (0.0768 - 0.0425)*0.5;


  tee.qcd_exp   = 0.3948; tee.qcd_stat       = 0.2055; tee.qcd_syst = tee.qcd_exp;
  tee.wjraw_exp =  3.2631; tee.wjraw_exp_stat = 0.9226;
  tee.wjf_systFrac = 0.5;

  tee.dy_exp = 4.033; tee.dy_stat = 1.194; tee.dy_syst = 0.5*tee.dy_exp;

  tee.SF_e = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  tee.tt_AE_eRel = 0.086;
 
  xsecCalc_inStruct(tee);


  std::cout<<"Doing JPT/TC"<<std::endl;
  std::cout<<"MM (w MET)"<<std::endl;
  TTxsecStruct tmm;
  tmm.tt_exp   = 18.886; tmm.tt_exp_me    =   0.265;
  tmm.dytt_exp = 0.787;  tmm.dytt_exp_me =   0.278;
  tmm.vv_exp   = 0.289;  tmm.vv_exp_me   =   0.022;
  tmm.tw_exp   = 0.614;  tmm.tw_exp_me   =   0.022;

  tmm.data = 29;

  tmm.sr  = 0.0306; tmm.sr_stat = 0.0052; tmm.sr_syst = fabs(0.0306 - 0.0188)*0.5;


  tmm.qcd_exp   = 0.0766; tmm.qcd_stat       = 0.0766; tmm.qcd_syst = tmm.qcd_exp;
  tmm.wjraw_exp =  2.4158; tmm.wjraw_exp_stat = 0.9520;
  tmm.wjf_systFrac = 0.75;

  tmm.dy_exp = 7.723; tmm.dy_stat = 1.794; tmm.dy_syst = 0.5*tmm.dy_exp;

  tmm.SF_e = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  tmm.tt_AE_eRel = 0.071;
 
  xsecCalc_inStruct(tmm);

  //=====================================================================================
  // NOW PF
  //=====================================================================================
  
  std::cout<<"Doing PF"<<std::endl;
  std::cout<<"EM (no MET)"<<std::endl;
  TTxsecStruct pem;
  pem.tt_exp   = 54.310; pem.tt_exp_me    =   0.449;
  pem.dytt_exp = 0.885;  pem.dytt_exp_me =   0.295;
  pem.vv_exp   = 0.952;  pem.vv_exp_me   =   0.039;
  pem.tw_exp   = 1.767;  pem.tw_exp_me   =   0.037;

  pem.data = 58;

  pem.sr  = 0.0544; pem.sr_stat = 0.0056; pem.sr_syst = fabs(0.0544 - 0.0307)*0.5;


  pem.qcd_exp   = 0.5322; pem.qcd_stat       = 0.2114; pem.qcd_syst = pem.qcd_exp;
  pem.wjraw_exp =  4.8042; pem.wjraw_exp_stat = 1.1378;
  pem.wjf_systFrac = 0.75;

  pem.dy_exp = 0; pem.dy_stat = 0; pem.dy_syst = 0;

  pem.SF_e = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  pem.tt_AE_eRel = 0.062;
 
  xsecCalc_inStruct(pem);


  std::cout<<"Doing PF"<<std::endl;
  std::cout<<"EM (w MET)"<<std::endl;
  TTxsecStruct peme;
  peme.tt_exp   = 50.181; peme.tt_exp_me    =   0.432;
  peme.dytt_exp = 0.787;  peme.dytt_exp_me =   0.278;
  peme.vv_exp   = 0.828;  peme.vv_exp_me   =   0.037;
  peme.tw_exp   = 1.627;  peme.tw_exp_me   =   0.036;

  peme.data = 50;

  peme.sr  = 0.0544; peme.sr_stat = 0.0056; peme.sr_syst = fabs(0.0544 - 0.0307)*0.5;


  peme.qcd_exp   = 0.2901; peme.qcd_stat       = 0.1687; peme.qcd_syst = peme.qcd_exp;
  peme.wjraw_exp =  3.3574; peme.wjraw_exp_stat = 0.9567;
  peme.wjf_systFrac = 0.75;

  peme.dy_exp = 0; peme.dy_stat = 0; peme.dy_syst = 0;

  peme.SF_e = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  peme.tt_AE_eRel = 0.062;
 
  xsecCalc_inStruct(peme);


  std::cout<<"Doing PF"<<std::endl;
  std::cout<<"EE (w MET)"<<std::endl;
  TTxsecStruct pee;
  pee.tt_exp   = 16.849; pee.tt_exp_me    =   0.250;
  pee.dytt_exp = 0.590;  pee.dytt_exp_me =   0.241;
  pee.vv_exp   = 0.281;  pee.vv_exp_me   =   0.021;
  pee.tw_exp   = 0.553;  pee.tw_exp_me   =   0.021;

  pee.data = 22;

  pee.sr  = 0.0795; pee.sr_stat = 0.0106; pee.sr_syst = fabs(0.0795 - 0.0425)*0.5;


  pee.qcd_exp   = 0.3007; pee.qcd_stat       = 0.1827; pee.qcd_syst = pee.qcd_exp;
  pee.wjraw_exp =  3.2624; pee.wjraw_exp_stat = 0.9341;
  pee.wjf_systFrac = 0.5;

  pee.dy_exp = 2.495; pee.dy_stat = 0.987; pee.dy_syst = 0.5* pee.dy_exp;

  pee.SF_e = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  pee.tt_AE_eRel = 0.086;
 
  xsecCalc_inStruct(pee);

  std::cout<<"Doing PF"<<std::endl;
  std::cout<<"MM (w MET)"<<std::endl;
  TTxsecStruct pmm;
  pmm.tt_exp   = 18.622; pmm.tt_exp_me    =   0.263;
  pmm.dytt_exp = 0.813;  pmm.dytt_exp_me =   0.279;
  pmm.vv_exp   = 0.286;  pmm.vv_exp_me   =   0.022;
  pmm.tw_exp   = 0.587;  pmm.tw_exp_me   =   0.022;

  pmm.data = 22;

  pmm.sr  = 0.0293; pmm.sr_stat = 0.0056; pmm.sr_syst = fabs(0.0293 - 0.0188)*0.5;


  pmm.qcd_exp   = 0.0766; pmm.qcd_stat       = 0.0766; pmm.qcd_syst = pmm.qcd_exp;
  pmm.wjraw_exp =  1.9240; pmm.wjraw_exp_stat = 0.8862;
  pmm.wjf_systFrac = 0.75;

  pmm.dy_exp = 4.636; pmm.dy_stat = 1.451; pmm.dy_syst = 0.5*pmm.dy_exp;

  pmm.SF_e = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  pmm.tt_AE_eRel = 0.071;
 
  xsecCalc_inStruct(pmm);


//   std::cout<<"Doing JPT/TC"<<std::endl;
//   std::cout<<"EE (w MET)"<<std::endl;
//   TTxsecStruct tee;
//   tee.tt_exp   = ; tee.tt_exp_me    =   ;
//   tee.dytt_exp = ;  tee.dytt_exp_me =   ;
//   tee.vv_exp   = ;  tee.vv_exp_me   =   ;
//   tee.tw_exp   = ;  tee.tw_exp_me   =   ;

//   tee.data = ;

//   tee.sr  = ; tee.sr_stat = ; tee.sr_syst = fabs( - )*;


//   tee.qcd_exp   = ; tee.qcd_stat       = ; tee.qcd_syst = tee.qcd_exp;
//   tee.wjraw_exp =  ; tee.wjraw_exp_stat = ;
//   tee.wjf_systFrac = ;

//   tee.dy_exp = ; tee.dy_stat = ; tee.dy_syst = ;

//   tee.SF_e = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
//   tee.tt_AE_eRel = ;
 
//   xsecCalc_inStruct(tee);
 
  
  
}
