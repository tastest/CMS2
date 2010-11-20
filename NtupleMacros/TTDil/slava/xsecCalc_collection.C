#include <iostream>

double oplus(double a, double b){
  return sqrt(a*a + b*b);
}
double oplus(double a, double b, double c){
  return sqrt(a*a + b*b + c*c);
}


void xsecCalc(double obs, double expS, double expSErr, double bg, double bgErr){
  double obsS = obs - bg;
  //notused  double obsSErr = sqrt( obs + bgErr*bgErr);

  double sigma = 157.5 * (obs - bg)/ expS;

  double sigmaErrStatFrac = sqrt(obs)/obsS;
  double sigmaErrSystFrac = oplus(bgErr/obsS, expSErr/expS);
  double sigmaErrFrac =  oplus(sigmaErrStatFrac, sigmaErrSystFrac);
  
  double sigmaErrFracLum = 0.11;

  double sigmaErrStat = sigma* sigmaErrStatFrac;
  double sigmaErrSyst = sigma* sigmaErrSystFrac;
  double sigmaErr = sigma*sigmaErrFrac;
  double sigmaErrLum = sigma*sigmaErrFracLum;

  double sigmaErrWLum = sigma * oplus(sigmaErrFrac, sigmaErrFracLum);

  std::cout<<"sigma = "
	   << sigma << " \\pm "<< sigmaErrStat << " (stat) \\pm " << sigmaErrSyst 
	   << " (syst) \\pm "<<sigmaErrLum << " (lum)"<<std::endl
	   <<"\t = "<<sigma << " \\pm " << sigmaErr << " (syst+stat) \\pm "<< sigmaErrLum << " (lum)"
	   <<"\t = "<<sigma << " \\pm " << sigmaErrWLum
	   <<"\n\t = "<<sigma<<" ( 1 \\pm "<<sigmaErrStat/sigma<<" (stat) \\pm "<<sigmaErrSyst/sigma 
	   <<" (syst) \\pm "<<sigmaErrFracLum<<" (lum) )"
	   <<"\n\t = "<<sigma<<" ( 1 \\pm "
	   <<sigmaErr/sigma << " (syst+stat) \\pm "
	   <<sigmaErrFracLum <<" )" 
	   <<"\t = "<< sigma<<" ( 1 \\pm "<<sigmaErrWLum/sigma << " )"
	   <<std::endl;

}


class TTxsecStruct {
public:
  double data;

  double tt_exp;       double tt_stat;

  double dytt_exp;     double dytt_stat;    double dytt_syst; 
  double vv_exp;       double vv_stat;      double vv_syst;
  double tw_exp;       double tw_stat;      double tw_syst;
  double mcbg_exp;     double mcbg_stat;    double mcbg_syst;  double mcbg_e;

  double sr_exp;       double sr_stat;      double sr_syst;
  double spill_exp;    double spill_stat;   double spill_syst; double spill_e;
  double qcd_exp;      double qcd_stat;     double qcd_syst;
  double wjraw_exp;    double wjraw_stat;
  double wjf_exp;      double wjf_stat;     double wjf_syst;
  double wjf_systFrac;
  double fake_exp;     double fake_stat;    double fake_syst;  double fake_e;

  double dy_exp;       double dy_stat;      double dy_syst;    double dy_e;

  double bg_exp;                                               double bg_e;

  double sf_exp;//
  double tt_AE_eRel;

  void setDependentPars_mcbg_v0(){
    dytt_syst = 0.5*dytt_exp;
    vv_syst   = 0.5*vv_exp;
    tw_syst   = 0.5*tw_exp;

    mcbg_exp = dytt_exp + vv_exp + tw_exp;
    mcbg_stat= oplus(dytt_stat, vv_stat, tw_stat);
    mcbg_syst= oplus(dytt_syst, vv_syst, tw_syst);
    mcbg_e   = oplus(mcbg_stat, mcbg_syst);
  }

  void setDependentPars_fake_v0(){
    spill_exp  = sr_exp * data;
    spill_stat = sqrt( data )*sr_exp;
    spill_syst = oplus(sr_stat, sr_syst)*data;
    spill_e    = oplus(spill_stat, spill_syst);

    wjf_exp = wjraw_exp - 2.* qcd_exp - spill_exp;
    if (wjf_exp< 0) std::cout<<"Warning: negative wj estimate: "<<wjf_exp<<std::endl;
    wjf_syst = wjf_systFrac*wjf_exp;

    fake_exp  = wjf_exp + qcd_exp;
    fake_stat = oplus(wjraw_stat, qcd_stat, spill_stat);
    fake_syst = oplus(wjf_syst, qcd_syst, spill_syst);

    fake_e    = oplus(fake_stat, fake_syst);
  }

  void setDependentPars_v0(){
    setDependentPars_mcbg_v0();
    setDependentPars_fake_v0();

    dy_e = oplus(dy_stat, dy_syst);

    bg_exp = mcbg_exp + fake_exp + dy_exp;
    bg_e = oplus(mcbg_e, fake_e, dy_e);

    //this is it
  }


//   void simpleAdd(TTxsecStruct& t){
//     data;
    
//     tt_exp;       tt_stat;
    
//     dytt_exp;     dytt_stat;    dytt_syst; 
//     vv_exp;       vv_stat;      vv_syst;
//     tw_exp;       tw_stat;      tw_syst;
//     mcbg_exp;     mcbg_stat;    mcbg_syst;  mcbg_e;
    
//     sr_exp;       sr_stat;      sr_syst;
//     spill_exp;    spill_stat;   spill_syst; spill_e;
//     qcd_exp;      qcd_stat;     qcd_syst;
//     wjraw_exp;    wjraw_stat;
//     wjf_exp;      wjf_stat;     wjf_syst;
//     wjf_systFrac;
//     fake_exp;     fake_stat;    fake_syst;  fake_e;
    
//     dy_exp;       dy_stat;      dy_syst;    dy_e;
    
//     bg_exp;                                               bg_e;
    
//     sf_exp;//
//     tt_AE_eRel;


//   }

  void printSummary(){
    std::cout<<"Data: "<<data
	     <<"\n\t   MCBG: "<<mcbg_exp <<" +/- " << mcbg_e
	     <<"\n\t   Fake: "<<fake_exp <<" +/- " << fake_stat <<" +/- "<< fake_syst
	     <<"\n\t     DY: "<<dy_exp <<" +/- " << dy_stat <<" +/- "<< dy_syst
	     <<"\n\t All BG: "<<bg_exp <<" +/- " << bg_e
	     <<"\n\t tt exp: "<<tt_exp*sf_exp<<" +/- "<<tt_exp*sf_exp*tt_AE_eRel
	     <<"\tS/B: "<<tt_exp*sf_exp/bg_exp<<std::endl;    
  }
};

void xsecCalc_inStruct(TTxsecStruct& tt){
  tt.setDependentPars_v0();
  
  tt.printSummary();

  xsecCalc(tt.data, tt.tt_exp*tt.sf_exp, tt.tt_exp*tt.sf_exp*tt.tt_AE_eRel, tt.bg_exp, tt.bg_e );
  
}

void xsecCalc_35pb_pass5(){
  std::cout<<"Doing JPT/TC"<<std::endl;
  std::cout<<"EMU (no MET)"<<std::endl;
  TTxsecStruct tem;
  tem.tt_exp   = 55.057; tem.tt_stat    =    0.452;
  tem.dytt_exp = 0.787 ;  tem.dytt_stat =    0.278;
  tem.vv_exp   = 1.050 ;  tem.vv_stat   =    0.041;
  tem.tw_exp   = 1.804 ;  tem.tw_stat   =    0.038;

  tem.data = 57;

  tem.sr_exp  = 0.0537; tem.sr_stat = 0.0052; tem.sr_syst = (0.0537 - 0.0307)*0.5;


  tem.qcd_exp   = 0.7029; tem.qcd_stat       = 0.2475; tem.qcd_syst = tem.qcd_exp;
  tem.wjraw_exp = 5.2129; tem.wjraw_stat = 1.1978;
  tem.wjf_systFrac = 0.75;

  tem.dy_exp = 0; tem.dy_stat = 0; tem.dy_syst = 0;

  tem.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  tem.tt_AE_eRel = 0.062;

  xsecCalc_inStruct(tem);

  std::cout<<"Doing JPT/TC"<<std::endl;
  std::cout<<"EMU (w MET)"<<std::endl;
  TTxsecStruct teme;
  teme.tt_exp   = 50.542; teme.tt_stat    =    0.433;
  teme.dytt_exp = 0.787 ;  teme.dytt_stat =    0.278;
  teme.vv_exp   = 0.890 ;  teme.vv_stat   =    0.038;
  teme.tw_exp   = 1.645 ;  teme.tw_stat   =    0.036;

  teme.data = 49;

  teme.sr_exp  = 0.0537; teme.sr_stat = 0.0052; teme.sr_syst = (0.0537 - 0.0307)*0.5;


  teme.qcd_exp   = 0.4480; teme.qcd_stat       = 0.2094; teme.qcd_syst = teme.qcd_exp;
  teme.wjraw_exp = 4.1669 ; teme.wjraw_stat = 1.0667;
  teme.wjf_systFrac = 0.75;

  teme.dy_exp = 0; teme.dy_stat = 0; teme.dy_syst = 0;

  teme.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  teme.tt_AE_eRel = 0.062;
 
  xsecCalc_inStruct(teme);
 
  std::cout<<"Doing JPT/TC"<<std::endl;
  std::cout<<"EE (w MET)"<<std::endl;
  TTxsecStruct tee;
  tee.tt_exp   = 16.849; tee.tt_stat    =   0.250;
  tee.dytt_exp = 0.492;  tee.dytt_stat =   0.220;
  tee.vv_exp   = 0.320;  tee.vv_stat   =   0.023;
  tee.tw_exp   = 0.567;  tee.tw_stat   =   0.021;

  tee.data = 21;

  tee.sr_exp  = 0.0768; tee.sr_stat = 0.0098; tee.sr_syst = (0.0768 - 0.0425)*0.5;


  tee.qcd_exp   = 0.3948; tee.qcd_stat       = 0.2055; tee.qcd_syst = tee.qcd_exp;
  tee.wjraw_exp =  3.2631; tee.wjraw_stat = 0.9226;
  tee.wjf_systFrac = 0.5;

  tee.dy_exp = 4.033; tee.dy_stat = 1.194; tee.dy_syst = 0.5*tee.dy_exp;

  tee.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  tee.tt_AE_eRel = 0.086;
 
  xsecCalc_inStruct(tee);


  std::cout<<"Doing JPT/TC"<<std::endl;
  std::cout<<"MM (w MET)"<<std::endl;
  TTxsecStruct tmm;
  tmm.tt_exp   = 18.886; tmm.tt_stat    =   0.265;
  tmm.dytt_exp = 0.787;  tmm.dytt_stat =   0.278;
  tmm.vv_exp   = 0.289;  tmm.vv_stat   =   0.022;
  tmm.tw_exp   = 0.614;  tmm.tw_stat   =   0.022;

  tmm.data = 29;

  tmm.sr_exp  = 0.0306; tmm.sr_stat = 0.0052; tmm.sr_syst = fabs(0.0306 - 0.0188)*0.5;


  tmm.qcd_exp   = 0.0766; tmm.qcd_stat       = 0.0766; tmm.qcd_syst = tmm.qcd_exp;
  tmm.wjraw_exp =  2.4158; tmm.wjraw_stat = 0.9520;
  tmm.wjf_systFrac = 0.75;

  tmm.dy_exp = 7.723; tmm.dy_stat = 1.794; tmm.dy_syst = 0.5*tmm.dy_exp;

  tmm.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  tmm.tt_AE_eRel = 0.071;
 
  xsecCalc_inStruct(tmm);

  //=====================================================================================
  // NOW PF
  //=====================================================================================
  
  std::cout<<"Doing PF"<<std::endl;
  std::cout<<"EM (no MET)"<<std::endl;
  TTxsecStruct pem;
  pem.tt_exp   = 54.310; pem.tt_stat    =   0.449;
  pem.dytt_exp = 0.885;  pem.dytt_stat =   0.295;
  pem.vv_exp   = 0.952;  pem.vv_stat   =   0.039;
  pem.tw_exp   = 1.767;  pem.tw_stat   =   0.037;

  pem.data = 58;

  pem.sr_exp  = 0.0544; pem.sr_stat = 0.0056; pem.sr_syst = fabs(0.0544 - 0.0307)*0.5;


  pem.qcd_exp   = 0.5322; pem.qcd_stat       = 0.2114; pem.qcd_syst = pem.qcd_exp;
  pem.wjraw_exp =  4.8042; pem.wjraw_stat = 1.1378;
  pem.wjf_systFrac = 0.75;

  pem.dy_exp = 0; pem.dy_stat = 0; pem.dy_syst = 0;

  pem.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  pem.tt_AE_eRel = 0.062;
 
  xsecCalc_inStruct(pem);


  std::cout<<"Doing PF"<<std::endl;
  std::cout<<"EM (w MET)"<<std::endl;
  TTxsecStruct peme;
  peme.tt_exp   = 50.181; peme.tt_stat    =   0.432;
  peme.dytt_exp = 0.787;  peme.dytt_stat =   0.278;
  peme.vv_exp   = 0.828;  peme.vv_stat   =   0.037;
  peme.tw_exp   = 1.627;  peme.tw_stat   =   0.036;

  peme.data = 50;

  peme.sr_exp  = 0.0544; peme.sr_stat = 0.0056; peme.sr_syst = fabs(0.0544 - 0.0307)*0.5;


  peme.qcd_exp   = 0.2901; peme.qcd_stat       = 0.1687; peme.qcd_syst = peme.qcd_exp;
  peme.wjraw_exp =  3.3574; peme.wjraw_stat = 0.9567;
  peme.wjf_systFrac = 0.75;

  peme.dy_exp = 0; peme.dy_stat = 0; peme.dy_syst = 0;

  peme.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  peme.tt_AE_eRel = 0.062;
 
  xsecCalc_inStruct(peme);


  std::cout<<"Doing PF"<<std::endl;
  std::cout<<"EE (w MET)"<<std::endl;
  TTxsecStruct pee;
  pee.tt_exp   = 16.849; pee.tt_stat    =   0.250;
  pee.dytt_exp = 0.590;  pee.dytt_stat =   0.241;
  pee.vv_exp   = 0.281;  pee.vv_stat   =   0.021;
  pee.tw_exp   = 0.553;  pee.tw_stat   =   0.021;

  pee.data = 22;

  pee.sr_exp  = 0.0795; pee.sr_stat = 0.0106; pee.sr_syst = fabs(0.0795 - 0.0425)*0.5;


  pee.qcd_exp   = 0.3007; pee.qcd_stat       = 0.1827; pee.qcd_syst = pee.qcd_exp;
  pee.wjraw_exp =  3.2624; pee.wjraw_stat = 0.9341;
  pee.wjf_systFrac = 0.5;

  pee.dy_exp = 2.495; pee.dy_stat = 0.987; pee.dy_syst = 0.5* pee.dy_exp;

  pee.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  pee.tt_AE_eRel = 0.086;
 
  xsecCalc_inStruct(pee);

  std::cout<<"Doing PF"<<std::endl;
  std::cout<<"MM (w MET)"<<std::endl;
  TTxsecStruct pmm;
  pmm.tt_exp   = 18.622; pmm.tt_stat    =   0.263;
  pmm.dytt_exp = 0.813;  pmm.dytt_stat =   0.279;
  pmm.vv_exp   = 0.286;  pmm.vv_stat   =   0.022;
  pmm.tw_exp   = 0.587;  pmm.tw_stat   =   0.022;

  pmm.data = 22;

  pmm.sr_exp  = 0.0293; pmm.sr_stat = 0.0056; pmm.sr_syst = fabs(0.0293 - 0.0188)*0.5;


  pmm.qcd_exp   = 0.0766; pmm.qcd_stat       = 0.0766; pmm.qcd_syst = pmm.qcd_exp;
  pmm.wjraw_exp =  1.9240; pmm.wjraw_stat = 0.8862;
  pmm.wjf_systFrac = 0.75;

  pmm.dy_exp = 4.636; pmm.dy_stat = 1.451; pmm.dy_syst = 0.5*pmm.dy_exp;

  pmm.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  pmm.tt_AE_eRel = 0.071;
 
  xsecCalc_inStruct(pmm);


//   std::cout<<"Doing JPT/TC"<<std::endl;
//   std::cout<<"EE (w MET)"<<std::endl;
//   TTxsecStruct tee;
//   tee.tt_exp   = ; tee.tt_stat    =   ;
//   tee.dytt_exp = ;  tee.dytt_stat =   ;
//   tee.vv_exp   = ;  tee.vv_stat   =   ;
//   tee.tw_exp   = ;  tee.tw_stat   =   ;

//   tee.data = ;

//   tee.sr_exp  = ; tee.sr_stat = ; tee.sr_syst = fabs( - )*; //syst is 0.5 of |sr_exp-sr_inclusive|


//   tee.qcd_exp   = ; tee.qcd_stat       = ; tee.qcd_syst = tee.qcd_exp;
//   tee.wjraw_exp =  ; tee.wjraw_stat = ;
//   tee.wjf_systFrac = ;

//   tee.dy_exp = ; tee.dy_stat = ; tee.dy_syst = ;

//   tee.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
//   tee.tt_AE_eRel = ;
 
//   xsecCalc_inStruct(tee);
 
  
  
}
