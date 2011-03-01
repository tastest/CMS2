#include <iostream>
#include <math.h>
#include <limits>

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

enum TTChannel {
  ee_ch = 1,
  mm_ch = 2,
  em_ch = 4,
  eemm_ch = 3,
  eeem_ch = 5,
  mmem_ch = 6,
  eemmem_ch = 7,
  other_ch = 8
};

class TTxsecStruct {
public:
  TTChannel channel;
  double data;

  double tt_exp;       double tt_stat;

  double dytt_exp;     double dytt_stat;    double dytt_syst_self;  double dytt_syst_corr;  double dytt_syst;  
  double vv_exp;       double vv_stat;      double vv_syst_self;    double vv_syst_corr;    double vv_syst;    
  double tw_exp;       double tw_stat;      double tw_syst_self;    double tw_syst_corr;    double tw_syst;    
  double mcbg_exp;     double mcbg_stat;    double mcbg_syst;  double mcbg_e;

  double ttotr_mc;     double ttotr_mc_stat;
  double wj_mc;        double wj_mc_stat;

  double sr_exp;       double sr_stat;      double sr_syst;
  double spill_exp;    double spill_stat;   double spill_syst; double spill_e;
  double qcd_exp;      double qcd_stat;     double qcd_syst;
  double wjraw_exp;    double wjraw_stat;
  double wjf_exp;      double wjf_stat;     double wjf_syst;
  double wjf_systFrac;
  double fake_exp;     double fake_stat;    double fake_syst;  double fake_e;

  double dy_mc;        double dy_mc_stat;
  double dy_exp;       double dy_stat;      double dy_syst;    double dy_e;
  double dy_roi;       double dy_roi_stat;

  double bg_exp;                                               double bg_e;

  double all_mc;       double all_mc_stat;

  double sf_exp;//
  double tt_AE_eRel;

  double tt_xsecTh_eRel;
  double lum_eRel;

  double tt_sigma;
  double lum_total;

  void setDependentPars_mcbg_v0(){
    dytt_syst_self = 0.3*dytt_exp;
    vv_syst_self   = 0.3*vv_exp;
    tw_syst_self   = 0.3*tw_exp;

//     dytt_syst_corr = oplus(0.11,0.1)*dytt_exp;
//     vv_syst_corr   = oplus(0.11,0.1)*vv_exp;
//     tw_syst_corr   = oplus(0.11,0.1)*tw_exp;

    //this is only interesting by itself, don't use in combinations with the rest
    dytt_syst      = oplus(dytt_syst_self, dytt_syst_corr);
    vv_syst      = oplus(vv_syst_self, vv_syst_corr);
    tw_syst      = oplus(tw_syst_self, tw_syst_corr);

    mcbg_exp = dytt_exp + vv_exp + tw_exp;
    mcbg_stat= oplus(dytt_stat, vv_stat, tw_stat);
    mcbg_syst= oplus(dytt_syst_self, vv_syst_self, tw_syst_self);
    mcbg_syst= oplus(mcbg_syst, dytt_syst_corr+vv_syst_corr+tw_syst_corr);
    mcbg_e   = oplus(mcbg_stat, mcbg_syst);

    all_mc = mcbg_exp + ttotr_mc + wj_mc + dy_mc + tt_exp;
    all_mc_stat = oplus(mcbg_stat, ttotr_mc_stat, wj_mc_stat);
    all_mc_stat = oplus(all_mc_stat, dy_mc_stat, tt_stat);
  }

  void setDependentPars_fake_v0(){
    double real_dils = data;
    // correct for fakes 
//     real_dils = real_dils - ( wjraw_exp - qcd_exp );
//     if (real_dils < 0 ){
//       std::cout<<"Warning: negative real_dils estimate: "<<real_dils<<" truncate to 0"<<std::endl;
//       real_dils = 0;      
//     }
    spill_exp  = sr_exp * real_dils;
    spill_stat = sqrt( real_dils )*sr_exp;
    spill_syst = oplus(sr_stat, sr_syst)*real_dils;
    spill_e    = oplus(spill_stat, spill_syst);

    wjf_exp = wjraw_exp - 2.* qcd_exp - spill_exp;
    if (wjf_exp< 0){ 
      std::cout<<"Warning: negative wj estimate: "<<wjf_exp<<" truncate to 0"<<std::endl;
      wjf_exp = 0;
    }
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


  void simpleAdd(TTxsecStruct& t, bool corrSyst=true){
    TTChannel oldChannel = channel;
    int aChI = (int)channel;
    int bChI = (int)t.channel;
    int chI = aChI | bChI;
    TTChannel newChannel = (TTChannel)chI;
    std::cout<<"Combining "<<aChI<<" and "<<bChI<<" into "<<chI<<std::endl;
    if (((aChI&1) == (bChI&1) && (aChI&1)==1)
	|| ((aChI&2) == (bChI&2) && (aChI&2)==2)
	|| ((aChI&4) == (bChI&4) && (aChI&4)==4)
	){
      std::cout<<"Can't simple-combine the same channels"<<std::endl;
      return;
    }
    channel = newChannel;

//     double dataOld = data;
    double tt_exp_old = tt_exp;
    data      += t.data;
    
    tt_exp    += t.tt_exp;
    tt_stat = oplus(tt_stat, t.tt_stat);
    
    dytt_exp  += t.dytt_exp;
    dytt_stat = oplus(dytt_stat, t.dytt_stat);
    dytt_syst_self += t.dytt_syst_self;
    dytt_syst_corr += t.dytt_syst_corr;
    dytt_syst = oplus(dytt_syst_self,  dytt_syst_corr);

    vv_exp    += t.vv_exp;
    vv_stat   = oplus(vv_stat, t.vv_stat);
    vv_syst_self += t.vv_syst_self;
    vv_syst_corr += t.vv_syst_corr;
    vv_syst = oplus(vv_syst_self,  vv_syst_corr);

    tw_exp    += t.tw_exp;
    tw_stat   = oplus(tw_stat, t.tw_stat);
    tw_syst_self += t.tw_syst_self;
    tw_syst_corr += t.tw_syst_corr;
    tw_syst = oplus(tw_syst_self,  tw_syst_corr);

    mcbg_exp  = dytt_exp + vv_exp + tw_exp;
    mcbg_stat = oplus(dytt_stat, vv_stat, tw_stat);
    mcbg_syst= oplus(dytt_syst_self, vv_syst_self, tw_syst_self);
    mcbg_syst= oplus(mcbg_syst, dytt_syst_corr+vv_syst_corr+tw_syst_corr);
    mcbg_e    = oplus(mcbg_stat, mcbg_syst);
    
    ttotr_mc       += t.ttotr_mc;
    ttotr_mc_stat   = oplus(ttotr_mc_stat, t.ttotr_mc_stat);
    wj_mc          += wj_mc;
    wj_mc_stat      = oplus(wj_mc_stat, t.wj_mc_stat);
    dy_mc          += t.dy_mc;
    dy_mc_stat      = oplus(dy_mc_stat, t.dy_mc_stat);

    //just recalculate from the above
    all_mc = mcbg_exp + ttotr_mc + wj_mc + dy_mc + tt_exp;
    all_mc_stat = oplus(mcbg_stat, ttotr_mc_stat, wj_mc_stat);
    all_mc_stat = oplus(all_mc_stat, dy_mc_stat, tt_stat);

    sr_exp=std::numeric_limits<double>::quiet_NaN();    
    sr_stat=std::numeric_limits<double>::quiet_NaN();      
    sr_syst=std::numeric_limits<double>::quiet_NaN();

    spill_exp += t.spill_exp;
    spill_stat= oplus(spill_stat, t.spill_stat);
    if ((oldChannel == ee_ch && t.channel == mm_ch) 
	|| (oldChannel == ee_ch && t.channel == mm_ch) ){
      spill_syst =  oplus(spill_syst, t.spill_syst);
    } 
    else if (corrSyst) spill_syst += t.spill_syst; 
    else spill_syst = oplus(spill_syst, t.spill_syst);
    spill_e   = oplus(spill_stat, spill_syst);

    qcd_exp   += t.qcd_exp;
    qcd_stat  = oplus(qcd_stat, t.qcd_stat);
    if ((oldChannel == ee_ch && t.channel == mm_ch) 
	|| (oldChannel == ee_ch && t.channel == mm_ch) ){
      qcd_syst = oplus(qcd_syst, t.qcd_syst);
    }    
    else if (corrSyst) qcd_syst += t.qcd_syst; 
    else qcd_syst = oplus(qcd_syst, t.qcd_syst);

    wjraw_exp += t.wjraw_exp;
    wjraw_stat= oplus(wjraw_stat, t.wjraw_stat);

    wjf_exp   = wjraw_exp - 2.* qcd_exp - spill_exp;
    if (wjf_exp< 0) std::cout<<"Warning: negative wj estimate: "<<wjf_exp<<std::endl;
    wjf_stat  = oplus(wjf_stat, t.wjf_stat);
    if ((oldChannel == ee_ch && t.channel == mm_ch) 
	|| (oldChannel == ee_ch && t.channel == mm_ch) ){
      wjf_syst = oplus(wjf_syst, t.wjf_syst);
    }
    else if (corrSyst) wjf_syst += t.wjf_syst; 
    else wjf_syst = oplus(wjf_syst, t.wjf_syst);
    wjf_systFrac = std::numeric_limits<double>::quiet_NaN();

    fake_exp  = wjf_exp + qcd_exp;
    fake_stat = oplus(wjraw_stat, qcd_stat, spill_stat);
    fake_syst = oplus(wjf_syst, qcd_syst, spill_syst);  
    fake_e    = oplus(fake_stat, fake_syst);
    
    dy_exp    += t.dy_exp;
    dy_stat   = oplus(dy_stat, t.dy_stat);
    dy_syst   = oplus(dy_syst, t.dy_syst);    
    dy_e      = oplus(dy_stat, dy_syst);
    
    bg_exp    += t.bg_exp;
    bg_e      = oplus(mcbg_e, fake_e, dy_e);
    
    //use simple average here
    sf_exp    = (sf_exp*tt_exp_old + t.sf_exp*t.tt_exp)/(tt_exp_old+t.tt_exp);
    tt_AE_eRel= (tt_AE_eRel*tt_exp_old + t.tt_AE_eRel*t.tt_exp)/(tt_exp_old+t.tt_exp);


  }

  void printSummary(){
    double l_sig_exp = tt_exp*sf_exp;
    double l_dytt_err = oplus(dytt_stat,dytt_syst);
    double l_vv_err = oplus(vv_stat,vv_syst);
    double l_tw_err = oplus(tw_stat,tw_syst);
    double l_fake_err = oplus(fake_stat, fake_syst);
    double l_dy_err = oplus(dy_stat, dy_syst);
    double l_sig_obs = data - bg_exp;
    std::cout<<"Data: "<<data<<" in "<<channel
	     <<"\n\t   MCBG: "<<mcbg_exp <<" +/- " << mcbg_e
	     <<"\n\t\t  = ( DYtt "<<dytt_exp<<" +/- "<<l_dytt_err << " ( "<<l_dytt_err/l_sig_obs<<" )"
	     <<"\n\t\t      VV   "<<vv_exp<<" +/- "<<l_vv_err << " ( "<<l_vv_err/l_sig_obs<<" )"
	     <<"\n\t\t      tW   "<<tw_exp<<" +/- "<<l_tw_err  << " ( "<<l_tw_err/l_sig_obs<<" )" <<"   )"
	     <<"\n\t   Fake: "<<fake_exp <<" +/- " << fake_stat <<" +/- "<< fake_syst  << " ( "<<l_fake_err/l_sig_obs<<" )"
	     <<"\n\t\t  ( QCDraw " << qcd_exp <<" +/- "<<qcd_stat<<" +/- "<<qcd_syst 
	     <<"\n\t\t\t WJraw "<<wjraw_exp<<" +/- "<<wjraw_stat
	     <<"\n\t\t\t Spill "<<spill_exp<<" +/- "<<spill_stat<<" +/- "<<spill_syst<< " )"
	     <<"\n\tFake MC: "<<ttotr_mc+wj_mc <<" +/- " << oplus(ttotr_mc_stat, wj_mc_stat)
	     <<"\n\t     DY: "<<dy_exp <<" +/- " << dy_stat <<" +/- "<< dy_syst  << " ( "<<l_dy_err/l_sig_obs<<" )"
	     <<"\n\t  DY MC: "<<dy_mc <<" +/- " << dy_mc_stat<<"\t R_{out/in}: "<<dy_roi<<" +/- "<< dy_roi_stat
	     <<"\n\t All BG: "<<bg_exp <<" +/- " << bg_e
	     <<"\n\t All MC: "<<all_mc <<" +/- "<<all_mc_stat

	     <<"\n\t tt exp: "<<l_sig_exp
	     <<" +/- "<<l_sig_exp*tt_AE_eRel<<"(syst) +/- "<<l_sig_exp*tt_xsecTh_eRel
	     <<"(th) +/- "<<l_sig_exp*lum_eRel<<"(lum)  "
	     <<"\n\t\t= "<<l_sig_exp
	     <<" +/- " <<l_sig_exp*tt_AE_eRel 
	     <<"(syst) +/- "
	     <<oplus(l_sig_exp*tt_xsecTh_eRel, l_sig_exp*lum_eRel) << "(th+lum)"
	     <<"\n\tA: "<<l_sig_exp/tt_sigma/lum_total*100.
	     <<" +/- "<<l_sig_exp/tt_sigma/lum_total*100.*tt_AE_eRel
	     <<" \tS/B: "<<l_sig_exp/bg_exp<<std::endl;    
  }
  
  TTxsecStruct(){
    tt_xsecTh_eRel = 0.15;
    lum_eRel       = 0.11;
    tt_sigma       = 157.5;
  }
};

void xsecCalc_inStruct(TTxsecStruct& tt, bool setDepPars = true){
  if (setDepPars) tt.setDependentPars_v0();
  
  tt.printSummary();

  xsecCalc(tt.data, tt.tt_exp*tt.sf_exp, tt.tt_exp*tt.sf_exp*tt.tt_AE_eRel, tt.bg_exp, tt.bg_e );
  std::cout<<"\nExpected performance: "<<std::endl;
  xsecCalc(tt.bg_exp+tt.tt_exp*tt.sf_exp, tt.tt_exp*tt.sf_exp, tt.tt_exp*tt.sf_exp*tt.tt_AE_eRel, tt.bg_exp, tt.bg_e );
  
}

void xsecCalc_35pb_pass5(){
  std::cout<<"Doing JPT/TC"<<std::endl;
  std::cout<<"EMU (no MET)"<<std::endl;
  TTxsecStruct tem;
  tem.channel = em_ch;
  tem.tt_exp   = 55.057;  tem.tt_stat   =    0.452;
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
  teme.channel = em_ch;
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
  tee.channel = ee_ch;
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
  tmm.channel = mm_ch;
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
  pem.channel = em_ch;
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
  peme.channel = em_ch;
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
  pee.channel = ee_ch;
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
  pmm.channel = mm_ch;
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



  bool bgSystIsCorr = false;
  //combine
  std::cout<<"======================================="<<std::endl;
  std::cout<<"======       COMBINATIONS      ========"<<std::endl; 
  std::cout<<"======================================="<<std::endl;
  std::cout<<"EEMM (wMET) PF"<<std::endl;
  TTxsecStruct peemm(pee);
  peemm.simpleAdd(pmm);
  xsecCalc_inStruct(peemm, false);

  std::cout<<"EEMMEM (wMET) PF"<<std::endl;
  TTxsecStruct peemmem(peemm);
  peemmem.simpleAdd(pem);
  xsecCalc_inStruct(peemmem, false);


  std::cout<<"======================================="<<std::endl;
  std::cout<<"======       BTAGGING          ========"<<std::endl; 
  std::cout<<"======================================="<<std::endl;
  /* Current strategy (how soon this comment goes out of sync with what's below?)
   * Use the same SR as in pre-tagged sample
   * simply add 5% to the acceptance systematics, based on +/-10% tag, +/-25% mistag values
   ( note that these will be averaged out in combination, which is not quite correct in 
   combination of tagged and untagged result)
  */

  std::cout<<"Doing PF w BTAG "<<std::endl;
  std::cout<<"EM (no MET)"<<std::endl;
  TTxsecStruct pbem;
  pbem.channel = em_ch;
  pbem.tt_exp   = 50.349;  pbem.tt_stat   =   0.433;
  pbem.dytt_exp =  0.197;  pbem.dytt_stat =   0.139;
  pbem.vv_exp   =  0.222;  pbem.vv_stat   =   0.019;
  pbem.tw_exp   =  1.526;  pbem.tw_stat   =   0.035;

  pbem.data = 52;

  pbem.sr_exp  = 0.0544; pbem.sr_stat = 0.0056; pbem.sr_syst = fabs(0.0544 - 0.0307)*0.5;


  pbem.qcd_exp   = 0.1792; pbem.qcd_stat  =  0.1273; pbem.qcd_syst = pbem.qcd_exp;
  pbem.wjraw_exp = 2.4110; pbem.wjraw_stat= 0.8263;
  pbem.wjf_systFrac = 0.75;

  pbem.dy_exp = 0; pbem.dy_stat = 0; pbem.dy_syst = 0;

  pbem.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  pbem.tt_AE_eRel = oplus(0.062, 0.05);
 
  xsecCalc_inStruct(pbem);


  std::cout<<"Doing PF w BTAG"<<std::endl;
  std::cout<<"EM (w MET)"<<std::endl;
  TTxsecStruct pbeme;
  pbeme.channel = em_ch;
  pbeme.tt_exp   = 46.499;  pbeme.tt_stat   =   0.416;
  pbeme.dytt_exp =  0.197;  pbeme.dytt_stat =   0.139;
  pbeme.vv_exp   =  0.211;  pbeme.vv_stat   =   0.019;
  pbeme.tw_exp   =  1.408;  pbeme.tw_stat   =   0.033;

  pbeme.data = 46;

  pbeme.sr_exp  = 0.0544; pbeme.sr_stat = 0.0056; pbeme.sr_syst = fabs(0.0544 - 0.0307)*0.5;


  pbeme.qcd_exp   = 0.0985; pbeme.qcd_stat  = 0.0985; pbeme.qcd_syst = pbeme.qcd_exp;
  pbeme.wjraw_exp = 2.4110; pbeme.wjraw_stat= 0.8263;
  pbeme.wjf_systFrac = 0.75;

  pbeme.dy_exp = 0; pbeme.dy_stat = 0; pbeme.dy_syst = 0;

  pbeme.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  pbeme.tt_AE_eRel = oplus(0.062, 0.05);
 
  xsecCalc_inStruct(pbeme);


  std::cout<<"Doing PF w BTAG"<<std::endl;
  std::cout<<"EE (w MET)"<<std::endl;
  TTxsecStruct pbee;
  pbee.channel = ee_ch;
  pbee.tt_exp   = 15.623;  pbee.tt_stat   =   0.241;
  pbee.dytt_exp =  0.295;  pbee.dytt_stat =   0.197;
  pbee.vv_exp   =  0.065;  pbee.vv_stat   =   0.082;
  pbee.tw_exp   =  0.487;  pbee.tw_stat   =   0.506;

  pbee.data = 15;

  pbee.sr_exp  = 0.0795; pbee.sr_stat = 0.0106; pbee.sr_syst = fabs(0.0795 - 0.0425)*0.5;


  pbee.qcd_exp   = 0.0953; pbee.qcd_stat  = 0.0953; pbee.qcd_syst = pbee.qcd_exp;
  pbee.wjraw_exp = 1.7113; pbee.wjraw_stat= 0.6696;
  pbee.wjf_systFrac = 0.5;

  pbee.dy_exp = 0.677; pbee.dy_stat = 0.586; pbee.dy_syst = 0.5* pbee.dy_exp;

  pbee.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  pbee.tt_AE_eRel = oplus(0.086, 0.05);
 
  xsecCalc_inStruct(pbee);

  std::cout<<"Doing PF w BTAG"<<std::endl;
  std::cout<<"MM (w MET)"<<std::endl;
  TTxsecStruct pbmm;
  pbmm.channel = mm_ch;
  pbmm.tt_exp   = 17.306;  pbmm.tt_stat   =   0.254;
  pbmm.dytt_exp =  0.197;  pbmm.dytt_stat =   0.197;
  pbmm.vv_exp   =  0.082;  pbmm.vv_stat   =   0.211;
  pbmm.tw_exp   =  0.506;  pbmm.tw_stat   =   1.408;

  pbmm.data = 20;

  pbmm.sr_exp  = 0.0293; pbmm.sr_stat = 0.0056; pbmm.sr_syst = fabs(0.0293 - 0.0188)*0.5;

  //!!!!! asymmetric error would be better
  pbmm.qcd_exp   = 0.0000; pbmm.qcd_stat  = 0.09; pbmm.qcd_syst = pbmm.qcd_exp;
  pbmm.wjraw_exp = 1.0140; pbmm.wjraw_stat = 0.6082;
  pbmm.wjf_systFrac = 0.75;

  pbmm.dy_exp = 1.365; pbmm.dy_stat = 0.820; pbmm.dy_syst = 0.5*pbmm.dy_exp;

  pbmm.sf_exp = (0.108*9.)*(0.108*9.);//multiply MC by this. This is for Madgraph !
  pbmm.tt_AE_eRel = oplus(0.071, 0.05);
 
  xsecCalc_inStruct(pbmm);
  
  //combine
  std::cout<<"======================================="<<std::endl;
  std::cout<<"======       COMBINATIONS      ========"<<std::endl; 
  std::cout<<"======================================="<<std::endl;
  std::cout<<"EEMM (wMET) PF w BTAG"<<std::endl;
  TTxsecStruct pbeemm(pbee);
  pbeemm.simpleAdd(pbmm);
  xsecCalc_inStruct(pbeemm, false);

  std::cout<<"EEMMEM (wMET) PF w BTAG"<<std::endl;
  TTxsecStruct pbeemmem(pbeemm);
  pbeemmem.simpleAdd(pbem);
  xsecCalc_inStruct(pbeemmem, false);

  std::cout<<"EEMMEM PF: eemm wMET,BTAG em noMET"<<std::endl;
  TTxsecStruct pbeemmNoBem(pbeemm);
  pbeemmNoBem.simpleAdd(pem);
  xsecCalc_inStruct(pbeemmNoBem, false);

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


void xsecCalc_36pb_pass6(){
  //=====================================================================================
  // NOW PF
  //=====================================================================================

  double rescaleMuSF = (0.992/0.9889);
  
  std::cout<<"#=============================================================#"<<std::endl;
  std::cout<<"# It's all PF                                                 #"<<std::endl;
  std::cout<<"#-------------------------------------------------------------#"<<std::endl;

  std::cout<<"\n=======\n\t EE 0 jets\n========="<<std::endl;
  TTxsecStruct pee0j;
  pee0j.lum_total = 36.1;
  double lum_eRel = pee0j.lum_eRel;
  pee0j.channel = ee_ch;
  pee0j.tt_exp   = 0.673;  pee0j.tt_stat    =   0.052;
  pee0j.dytt_exp = 0.618;  pee0j.dytt_stat =   0.165; pee0j.dytt_syst_corr = oplus(0.1,lum_eRel)*pee0j.dytt_exp;
  pee0j.vv_exp   = 3.122;  pee0j.vv_stat   =   0.543; pee0j.vv_syst_corr = oplus(0.1,lum_eRel)*pee0j.vv_exp;
  pee0j.tw_exp   = 0.199;  pee0j.tw_stat   =   0.012; pee0j.tw_syst_corr = oplus(0.1,lum_eRel)*pee0j.tw_exp;    

  pee0j.data = 12;

  pee0j.sr_exp  = 0.0331; pee0j.sr_stat = 0.0011; pee0j.sr_syst = 0.0004; 

  pee0j.ttotr_mc = 0.012;  pee0j.ttotr_mc_stat  =   0.007;
  pee0j.wj_mc    = 3.346;  pee0j.wj_mc_stat     =   0.543;

  pee0j.qcd_exp   = 0.3730; pee0j.qcd_stat  =  0.1558; pee0j.qcd_syst = pee0j.qcd_exp;
  pee0j.wjraw_exp = 4.8377; pee0j.wjraw_stat=  1.1990;
  pee0j.wjf_systFrac = 0.5;

  pee0j.dy_mc  =    1.79;  pee0j.dy_mc_stat = 0.36;
  pee0j.dy_exp = 3.58; pee0j.dy_stat = 1.27; pee0j.dy_syst = pee0j.dy_exp*0.5;
  pee0j.dy_roi = (0.2807+0.1630)*0.5; pee0j.dy_roi_stat = (0.0513+0.0358)*0.5;

  pee0j.sf_exp = (0.108*9.)*(0.108*9.) * 0.9231* 1.0126;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pee0j.tt_AE_eRel = 0.113; //use the same as 1 jet
 
  xsecCalc_inStruct(pee0j);

  std::cout<<"\n=======\n\t EE 1 jets\n========="<<std::endl;
  TTxsecStruct pee1j;
  pee1j.lum_total = 36.1;
  pee1j.channel = ee_ch;
  pee1j.tt_exp   = 5.459; pee1j.tt_stat    =   0.149;
  pee1j.dytt_exp = 1.500;  pee1j.dytt_stat =   0.257; pee1j.dytt_syst_corr = oplus(0.1,lum_eRel)*pee1j.dytt_exp;
  pee1j.vv_exp   = 0.975;  pee1j.vv_stat   =   0.035; pee1j.vv_syst_corr = oplus(0.1,lum_eRel)*pee1j.vv_exp;	
  pee1j.tw_exp   = 0.913;  pee1j.tw_stat   =   0.027; pee1j.tw_syst_corr = oplus(0.1,lum_eRel)*pee1j.tw_exp;    

  pee1j.data = 19;

  pee1j.sr_exp   = 0.0340; pee1j.sr_stat = 0.0030; pee1j.sr_syst = 0.0001; 

  pee1j.ttotr_mc = 0.098;  pee1j.ttotr_mc_stat  =   0.020;
  pee1j.wj_mc    = 1.639;  pee1j.wj_mc_stat     =   0.377;

  pee1j.qcd_exp   = 0.5508; pee1j.qcd_stat  =  0.2193; pee1j.qcd_syst = pee1j.qcd_exp;
  pee1j.wjraw_exp = 3.4900; pee1j.wjraw_stat=  0.9813;
  pee1j.wjf_systFrac = 0.5;

  pee1j.dy_mc  =    1.31;  pee1j.dy_mc_stat = 0.285;
  pee1j.dy_exp =    3.23; pee1j.dy_stat =    0.98; pee1j.dy_syst =   pee1j.dy_exp*0.5;
  pee1j.dy_roi = (0.1729+0.0870)*0.5; pee1j.dy_roi_stat = (0.0306+0.0239)*0.5;

  pee1j.sf_exp = (0.108*9.)*(0.108*9.) * 0.9231* 1.0126;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pee1j.tt_AE_eRel = 0.113; //using FOM-based systematics (too lazy to redo it here for these )
 
  xsecCalc_inStruct(pee1j);


  std::cout<<"\n=======\n\t EE 2 jets\n========="<<std::endl;
  TTxsecStruct pee2j;
  pee2j.lum_total = 36.1;
  pee2j.channel = ee_ch;
  pee2j.tt_exp   =16.659; pee2j.tt_stat    =   0.261;
  pee2j.dytt_exp = 0.574;  pee2j.dytt_stat =   0.159; pee2j.dytt_syst_corr = oplus(0.1,lum_eRel)*pee2j.dytt_exp;
  pee2j.vv_exp   = 0.282;  pee2j.vv_stat   =   0.018; pee2j.vv_syst_corr = oplus(0.1,lum_eRel)*pee2j.vv_exp;	
  pee2j.tw_exp   = 0.562;  pee2j.tw_stat   =   0.021; pee2j.tw_syst_corr = oplus(0.1,lum_eRel)*pee2j.tw_exp;    

  pee2j.data = 23;

  pee2j.sr_exp  = 0.063; pee2j.sr_stat = 0.009; pee2j.sr_syst = 0.01455; 

  pee2j.ttotr_mc = 0.379;  pee2j.ttotr_mc_stat  =   0.039;
  pee2j.wj_mc    = 0.178;  pee2j.wj_mc_stat     =   0.126;

  pee2j.qcd_exp   = 0.2587; pee2j.qcd_stat  =  0.1665; pee2j.qcd_syst = pee2j.qcd_exp;
  pee2j.wjraw_exp = 3.6512; pee2j.wjraw_stat=  1.0043;
  pee2j.wjf_systFrac = 0.5;

  pee2j.dy_mc  =    1.72;  pee2j.dy_mc_stat = 0.34;
  pee2j.dy_exp =    3.00; pee2j.dy_stat =  0.97; pee2j.dy_syst =   pee2j.dy_exp*0.5;
  pee2j.dy_roi = (0.1174+0.1610)*0.5; pee2j.dy_roi_stat = (0.0221+0.0345)*0.5;

  pee2j.sf_exp = (0.108*9.)*(0.108*9.) * 0.9231* 1.0126;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pee2j.tt_AE_eRel = 0.082; //mind the correlations
 
  xsecCalc_inStruct(pee2j);


  std::cout<<"\n=======\n\t EE 2 jets combined fakes\n========="<<std::endl;
  TTxsecStruct pee2jFComb;
  pee2jFComb.lum_total = 36.1;
  pee2jFComb.channel = ee_ch;
  pee2jFComb.tt_exp   =16.659; pee2jFComb.tt_stat    =   0.261;
  pee2jFComb.dytt_exp = 0.574;  pee2jFComb.dytt_stat =   0.159; pee2jFComb.dytt_syst_corr = oplus(0.1,lum_eRel)*pee2jFComb.dytt_exp;
  pee2jFComb.vv_exp   = 0.282;  pee2jFComb.vv_stat   =   0.018; pee2jFComb.vv_syst_corr = oplus(0.1,lum_eRel)*pee2jFComb.vv_exp;	
  pee2jFComb.tw_exp   = 0.562;  pee2jFComb.tw_stat   =   0.021; pee2jFComb.tw_syst_corr = oplus(0.1,lum_eRel)*pee2jFComb.tw_exp;    

  pee2jFComb.data = 23;

  pee2jFComb.sr_exp  = 0.0; pee2jFComb.sr_stat = 0.0; pee2jFComb.sr_syst = 0.0; 

  pee2jFComb.ttotr_mc = 0.379;  pee2jFComb.ttotr_mc_stat  =   0.039;
  pee2jFComb.wj_mc    = 0.178;  pee2jFComb.wj_mc_stat     =   0.126;

  pee2jFComb.qcd_exp   = 0.0; pee2jFComb.qcd_stat  =  0.0; pee2jFComb.qcd_syst = pee2jFComb.qcd_exp;
  pee2jFComb.wjraw_exp = 1.06; pee2jFComb.wjraw_stat=  1.44;
  pee2jFComb.wjf_systFrac = 0.0;

  pee2jFComb.dy_mc  =    1.72;  pee2jFComb.dy_mc_stat = 0.34;
  pee2jFComb.dy_exp =    3.00; pee2jFComb.dy_stat =  0.97; pee2jFComb.dy_syst =   pee2jFComb.dy_exp*0.5;
  pee2jFComb.dy_roi = (0.1174+0.1610)*0.5; pee2jFComb.dy_roi_stat = (0.0221+0.0345)*0.5;

  pee2jFComb.sf_exp = (0.108*9.)*(0.108*9.) * 0.9231* 1.0126;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pee2jFComb.tt_AE_eRel = 0.082; //mind the correlations
 
  xsecCalc_inStruct(pee2jFComb);




  std::cout<<"\n=======\n\t MM 0 jets\n========="<<std::endl;
  TTxsecStruct pmm0j;
  pmm0j.lum_total = 36.1;
  pmm0j.channel = mm_ch;
  pmm0j.tt_exp   = 0.754; pmm0j.tt_stat    =   0.055;
  pmm0j.dytt_exp = 0.574;  pmm0j.dytt_stat =   0.159; pmm0j.dytt_syst_corr = oplus(0.1,lum_eRel)*pmm0j.dytt_exp;
  pmm0j.vv_exp   = 3.931;  pmm0j.vv_stat   =   0.073; pmm0j.vv_syst_corr = oplus(0.1,lum_eRel)*pmm0j.vv_exp;	
  pmm0j.tw_exp   = 0.223;  pmm0j.tw_stat   =   0.013; pmm0j.tw_syst_corr = oplus(0.1,lum_eRel)*pmm0j.tw_exp;    

  pmm0j.data = 11;

  pmm0j.sr_exp  = 0.0165; pmm0j.sr_stat = 0.0007; pmm0j.sr_syst = 0.0004; 

  pmm0j.ttotr_mc = 0.000;  pmm0j.ttotr_mc_stat  =   0.004;
  pmm0j.wj_mc    = 0.000;  pmm0j.wj_mc_stat     =   0.111;

  pmm0j.qcd_exp   = 0.0000; pmm0j.qcd_stat  =  0.09; pmm0j.qcd_syst = pmm0j.qcd_exp;
  pmm0j.wjraw_exp = 0.4874; pmm0j.wjraw_stat=  0.4874;
  pmm0j.wjf_systFrac = 0.75;

  pmm0j.dy_mc  =  4.59;  pmm0j.dy_mc_stat = 0.575;
  pmm0j.dy_exp =  14.505; pmm0j.dy_stat = 3.225; pmm0j.dy_syst = pmm0j.dy_exp*0.5;
  pmm0j.dy_roi = (0.6044+0.2058)*0.5; pmm0j.dy_roi_stat = (0.0660+0.0293)*0.5;

  pmm0j.sf_exp = (0.108*9.)*(0.108*9.) * 0.9613* 1.0126* rescaleMuSF*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pmm0j.tt_AE_eRel = 0.112; //use the same as 1 jet
 
  xsecCalc_inStruct(pmm0j);

  std::cout<<"\n=======\n\t MM 1 jets\n========="<<std::endl;
  TTxsecStruct pmm1j;
  pmm1j.lum_total = 36.1;
  pmm1j.channel = mm_ch;
  pmm1j.tt_exp   = 6.458; pmm1j.tt_stat    =   0.162;
  pmm1j.dytt_exp = 1.897;  pmm1j.dytt_stat =   0.289; pmm1j.dytt_syst_corr = oplus(0.1,lum_eRel)*pmm1j.dytt_exp;
  pmm1j.vv_exp   = 1.163;  pmm1j.vv_stat   =   0.038; pmm1j.vv_syst_corr = oplus(0.1,lum_eRel)*pmm1j.vv_exp;	
  pmm1j.tw_exp   = 1.066;  pmm1j.tw_stat   =   0.029; pmm1j.tw_syst_corr = oplus(0.1,lum_eRel)*pmm1j.tw_exp;    

  pmm1j.data = 19;

  pmm1j.sr_exp   = 0.0220; pmm1j.sr_stat = 0.0023; pmm1j.sr_syst = 0.0023; 

  pmm1j.ttotr_mc = 0.012;  pmm1j.ttotr_mc_stat  =   0.007;
  pmm1j.wj_mc    = 0.111;  pmm1j.wj_mc_stat     =   0.111;

  pmm1j.qcd_exp   = 0.2225; pmm1j.qcd_stat  =  0.1573; pmm1j.qcd_syst = pmm1j.qcd_exp;
  pmm1j.wjraw_exp = 0.2892; pmm1j.wjraw_stat=  0.2892;
  pmm1j.wjf_systFrac = 0.75;

  pmm1j.dy_mc  =    6.185;  pmm1j.dy_mc_stat = 0.665;
  pmm1j.dy_exp =    13.18; pmm1j.dy_stat =    2.775; pmm1j.dy_syst =   pmm1j.dy_exp*0.5;
  pmm1j.dy_roi = (0.4493+0.3013)*0.5; pmm1j.dy_roi_stat = (0.0391+0.0385)*0.5;

  pmm1j.sf_exp = (0.108*9.)*(0.108*9.) * 0.9613* 1.0126*rescaleMuSF*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pmm1j.tt_AE_eRel = 0.112; //use FOM-based values, too lazy to recompute for these specific cuts
 
  xsecCalc_inStruct(pmm1j);


  std::cout<<"\n=======\n\t MM 2 jets\n========="<<std::endl;
  TTxsecStruct pmm2j;
  pmm2j.lum_total = 36.1;
  pmm2j.channel = mm_ch;
  pmm2j.tt_exp   =19.925; pmm2j.tt_stat    =   0.285;
  pmm2j.dytt_exp = 0.515;  pmm2j.dytt_stat =   0.149; pmm2j.dytt_syst_corr = oplus(0.1,lum_eRel)*pmm2j.dytt_exp;
  pmm2j.vv_exp   = 0.313;  pmm2j.vv_stat   =   0.018; pmm2j.vv_syst_corr = oplus(0.1,lum_eRel)*pmm2j.vv_exp;	
  pmm2j.tw_exp   = 0.658;  pmm2j.tw_stat   =   0.023; pmm2j.tw_syst_corr = oplus(0.1,lum_eRel)*pmm2j.tw_exp;    

  pmm2j.data = 28;

  pmm2j.sr_exp  = 0.033; pmm2j.sr_stat = 0.005; pmm2j.sr_syst = 0.0060; 

  pmm2j.ttotr_mc = 0.077;  pmm2j.ttotr_mc_stat  =   0.018;
  pmm2j.wj_mc    = 0.000;  pmm2j.wj_mc_stat     =   0.111;

  pmm2j.qcd_exp   = 0.0841; pmm2j.qcd_stat  =  0.0841; pmm2j.qcd_syst = pmm2j.qcd_exp;
  pmm2j.wjraw_exp = 2.0411; pmm2j.wjraw_stat=  0.8397;
  pmm2j.wjf_systFrac = 0.75;

  pmm2j.dy_mc  =    3.345;  pmm2j.dy_mc_stat = 0.49;
  pmm2j.dy_exp =     7.44;  pmm2j.dy_stat =    1.825; pmm2j.dy_syst =   pmm2j.dy_exp*0.5;
  pmm2j.dy_roi = (0.2431+0.1909)*0.5; pmm2j.dy_roi_stat = (0.0285+0.0324)*0.5;

  pmm2j.sf_exp = (0.108*9.)*(0.108*9.) * 0.9613* 1.0126*rescaleMuSF*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pmm2j.tt_AE_eRel = 0.081; //mind the correlations
 
  xsecCalc_inStruct(pmm2j);


  std::cout<<"\n=======\n\t MM 2 jets Fakes combined w AN406\n========="<<std::endl;
  TTxsecStruct pmm2jFComb;
  pmm2jFComb.lum_total = 36.1;
  pmm2jFComb.channel = mm_ch;
  pmm2jFComb.tt_exp   =19.925; pmm2jFComb.tt_stat    =   0.285;
  pmm2jFComb.dytt_exp = 0.515;  pmm2jFComb.dytt_stat =   0.149; pmm2jFComb.dytt_syst_corr = oplus(0.1,lum_eRel)*pmm2jFComb.dytt_exp;
  pmm2jFComb.vv_exp   = 0.313;  pmm2jFComb.vv_stat   =   0.018; pmm2jFComb.vv_syst_corr = oplus(0.1,lum_eRel)*pmm2jFComb.vv_exp;	
  pmm2jFComb.tw_exp   = 0.658;  pmm2jFComb.tw_stat   =   0.023; pmm2jFComb.tw_syst_corr = oplus(0.1,lum_eRel)*pmm2jFComb.tw_exp;    

  pmm2jFComb.data = 28;

  pmm2jFComb.sr_exp  = 0.0; pmm2jFComb.sr_stat = 0.0; pmm2jFComb.sr_syst = 0.0; 

  pmm2jFComb.ttotr_mc = 0.077;  pmm2jFComb.ttotr_mc_stat  =   0.018;
  pmm2jFComb.wj_mc    = 0.000;  pmm2jFComb.wj_mc_stat     =   0.111;

  pmm2jFComb.qcd_exp   = 0.0; pmm2jFComb.qcd_stat  =  0.0; pmm2jFComb.qcd_syst = pmm2jFComb.qcd_exp;
  pmm2jFComb.wjraw_exp = 0.60; pmm2jFComb.wjraw_stat=  1.14;
  pmm2jFComb.wjf_systFrac = 0.0;

  pmm2jFComb.dy_mc  =    3.345;  pmm2jFComb.dy_mc_stat = 0.49;
  pmm2jFComb.dy_exp =     7.44;  pmm2jFComb.dy_stat =    1.825; pmm2jFComb.dy_syst =   pmm2jFComb.dy_exp*0.5;
  pmm2jFComb.dy_roi = (0.2431+0.1909)*0.5; pmm2jFComb.dy_roi_stat = (0.0285+0.0324)*0.5;

  pmm2jFComb.sf_exp = (0.108*9.)*(0.108*9.) * 0.9613* 1.0126*rescaleMuSF*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pmm2jFComb.tt_AE_eRel = 0.081; //mind the correlations
 
  xsecCalc_inStruct(pmm2jFComb);



  std::cout<<"\n=======\n\t EM 0 jets\n========="<<std::endl;
  TTxsecStruct pem0j;
  pem0j.lum_total = 36.1;
  pem0j.channel = em_ch;
  pem0j.tt_exp   = 2.100; pem0j.tt_stat    =   0.093;
  pem0j.dytt_exp = 65.567;  pem0j.dytt_stat =   1.701; pem0j.dytt_syst_corr = oplus(0.1,lum_eRel)*pem0j.dytt_exp;
  pem0j.vv_exp   = 13.943;  pem0j.vv_stat   =   0.140; pem0j.vv_syst_corr = oplus(0.1,lum_eRel)*pem0j.vv_exp;	
  pem0j.tw_exp   = 0.688;  pem0j.tw_stat   =   0.023;  pem0j.tw_syst_corr = oplus(0.1,lum_eRel)*pem0j.tw_exp;    

  pem0j.data = 79;

  pem0j.sr_exp  = 0.0248; pem0j.sr_stat = 0.0007; pem0j.sr_syst = 0.0004; 

  pem0j.ttotr_mc = 0.029;  pem0j.ttotr_mc_stat  =   0.011;
  pem0j.wj_mc    = 13.747;  pem0j.wj_mc_stat    =   1.221;

  pem0j.qcd_exp   = 5.3854; pem0j.qcd_stat  =  0.7161; pem0j.qcd_syst = pem0j.qcd_exp;
  pem0j.wjraw_exp = 23.9090; pem0j.wjraw_stat=  2.7220;
  pem0j.wjf_systFrac = 0.5;

  pem0j.dy_mc  =3.331;  pem0j.dy_mc_stat = 0.375;
  pem0j.dy_exp =  0.0;  pem0j.dy_stat =    0.0; pem0j.dy_syst = pem0j.dy_exp*0.5;
  pem0j.dy_roi = 0;     pem0j.dy_roi_stat = 0;

  pem0j.sf_exp = (0.108*9.)*(0.108*9.) * 0.9444* 1.0126*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pem0j.tt_AE_eRel = 0.138; //use the same as 1jet in FOM (this isn't quite right though)
 
  xsecCalc_inStruct(pem0j);

  std::cout<<"\n=======\n\t EM 1 jets\n========="<<std::endl;
  TTxsecStruct pem1j;
  pem1j.lum_total = 36.1;
  pem1j.channel = em_ch;
  pem1j.tt_exp   = 17.936; pem1j.tt_stat    =   0.270;
  pem1j.dytt_exp = 10.487;  pem1j.dytt_stat =   0.680; pem1j.dytt_syst_corr = oplus(0.1,lum_eRel)*pem1j.dytt_exp;
  pem1j.vv_exp   = 3.852;  pem1j.vv_stat   =   0.069;  pem1j.vv_syst_corr = oplus(0.1,lum_eRel)*pem1j.vv_exp;	
  pem1j.tw_exp   = 3.130;  pem1j.tw_stat   =   0.049;  pem1j.tw_syst_corr = oplus(0.1,lum_eRel)*pem1j.tw_exp;    

  pem1j.data = 37;

  pem1j.sr_exp   = 0.028; pem1j.sr_stat = 0.002; pem1j.sr_syst = 0.0012; 

  pem1j.ttotr_mc = 0.220;  pem1j.ttotr_mc_stat  =   0.030;
  pem1j.wj_mc    = 2.199;  pem1j.wj_mc_stat     =   0.484;

  pem1j.qcd_exp   = 1.6999; pem1j.qcd_stat  =  0.3630; pem1j.qcd_syst = pem1j.qcd_exp;
  pem1j.wjraw_exp = 12.4671; pem1j.wjraw_stat=  1.9140;
  pem1j.wjf_systFrac = 0.5;

  pem1j.dy_mc  =    0.804;  pem1j.dy_mc_stat = 0.185;
  pem1j.dy_exp =    0.0; pem1j.dy_stat =     0.0; pem1j.dy_syst =   pem1j.dy_exp*0.5;
  pem1j.dy_roi = 0;     pem1j.dy_roi_stat = 0;

  pem1j.sf_exp = (0.108*9.)*(0.108*9.) * 0.9444* 1.0126*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pem1j.tt_AE_eRel = 0.138; //use the same as for the FOM Njet=1 cuts
 
  xsecCalc_inStruct(pem1j);


  std::cout<<"\n=======\n\t EM 2 jets\n========="<<std::endl;
  TTxsecStruct pem2j;
  pem2j.lum_total = 36.1;
  pem2j.channel = em_ch;
  pem2j.tt_exp   = 58.203; pem2j.tt_stat    =   0.487;
  pem2j.dytt_exp = 2.500;  pem2j.dytt_stat =   0.331; pem2j.dytt_syst_corr = oplus(0.1,lum_eRel)*pem2j.dytt_exp;
  pem2j.vv_exp   = 0.896;  pem2j.vv_stat   =   0.032; pem2j.vv_syst_corr = oplus(0.1,lum_eRel)*pem2j.vv_exp;	
  pem2j.tw_exp   = 1.862;  pem2j.tw_stat   =   0.038; pem2j.tw_syst_corr = oplus(0.1,lum_eRel)*pem2j.tw_exp;    

  pem2j.data = 60;

  pem2j.sr_exp  = 0.048; pem2j.sr_stat = 0.0052; pem2j.sr_syst = 0.0079; 

  pem2j.ttotr_mc = 0.934;  pem2j.ttotr_mc_stat  =   0.062;
  pem2j.wj_mc    = 0.646;  pem2j.wj_mc_stat     =   0.264;

  pem2j.qcd_exp   = 0.5717; pem2j.qcd_stat  =  0.2541; pem2j.qcd_syst = pem2j.qcd_exp;
  pem2j.wjraw_exp = 4.9342; pem2j.wjraw_stat=  1.2022;
  pem2j.wjf_systFrac = 0.5;

  pem2j.dy_mc  =    0.423;  pem2j.dy_mc_stat = 0.134;
  pem2j.dy_exp =    0.0;  pem2j.dy_stat =    0.0; pem2j.dy_syst =   pem2j.dy_exp*0.5;
  pem2j.dy_roi = 0;     pem2j.dy_roi_stat = 0;

  pem2j.sf_exp = (0.108*9.)*(0.108*9.) * 0.9444* 1.0126*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pem2j.tt_AE_eRel = 0.063; //mind the correlations
 
  xsecCalc_inStruct(pem2j);


  std::cout<<"\n=======\n\t EM 2 jets fakes combined with AN406\n========="<<std::endl;
  TTxsecStruct pem2jFComb;
  pem2jFComb.lum_total = 36.1;
  pem2jFComb.channel = em_ch;
  pem2jFComb.tt_exp   = 58.203; pem2jFComb.tt_stat    =   0.487;
  pem2jFComb.dytt_exp = 2.500;  pem2jFComb.dytt_stat =   0.331; pem2jFComb.dytt_syst_corr = oplus(0.1,lum_eRel)*pem2jFComb.dytt_exp;
  pem2jFComb.vv_exp   = 0.896;  pem2jFComb.vv_stat   =   0.032; pem2jFComb.vv_syst_corr = oplus(0.1,lum_eRel)*pem2jFComb.vv_exp;	
  pem2jFComb.tw_exp   = 1.862;  pem2jFComb.tw_stat   =   0.038; pem2jFComb.tw_syst_corr = oplus(0.1,lum_eRel)*pem2jFComb.tw_exp;    

  pem2jFComb.data = 60;

  pem2jFComb.sr_exp  = 0.0; pem2jFComb.sr_stat = 0.00; pem2jFComb.sr_syst = 0.00; 

  pem2jFComb.ttotr_mc = 0.934;  pem2jFComb.ttotr_mc_stat  =   0.062;
  pem2jFComb.wj_mc    = 0.646;  pem2jFComb.wj_mc_stat     =   0.264;

  pem2jFComb.qcd_exp   = 0.0; pem2jFComb.qcd_stat  =  0.0; pem2jFComb.qcd_syst = pem2jFComb.qcd_exp;
  pem2jFComb.wjraw_exp = 1.44; pem2jFComb.wjraw_stat=  1.58;
  pem2jFComb.wjf_systFrac = 0.0;

  pem2jFComb.dy_mc  =    0.423;  pem2jFComb.dy_mc_stat = 0.134;
  pem2jFComb.dy_exp =    0.0;  pem2jFComb.dy_stat =    0.0; pem2jFComb.dy_syst =   pem2jFComb.dy_exp*0.5;
  pem2jFComb.dy_roi = 0;     pem2jFComb.dy_roi_stat = 0;

  pem2jFComb.sf_exp = (0.108*9.)*(0.108*9.) * 0.9444* 1.0126*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pem2jFComb.tt_AE_eRel = 0.063; //mind the correlations
 
  xsecCalc_inStruct(pem2jFComb);


  std::cout<<"\n\n\n\n"<<std::endl;
  std::cout<<"#=============================================================#"<<std::endl;
  std::cout<<"# It's all PF with b-tags                                     #"<<std::endl;
  std::cout<<"#-------------------------------------------------------------#"<<std::endl;

  std::cout<<"\n=======\n\t EE 1 jets with b-tags\n========="<<std::endl;
  TTxsecStruct pee1j1bj;
  pee1j1bj.lum_total = 36.1;
  pee1j1bj.channel = ee_ch;
  pee1j1bj.tt_exp   = 3.918; pee1j1bj.tt_stat    =   0.126;
  pee1j1bj.dytt_exp = 0.088;  pee1j1bj.dytt_stat =   0.062; pee1j1bj.dytt_syst_corr = oplus(0.1,lum_eRel,0.25)*pee1j1bj.dytt_exp;
  pee1j1bj.vv_exp   = 0.104;  pee1j1bj.vv_stat   =   0.012; pee1j1bj.vv_syst_corr = oplus(0.1,lum_eRel,0.25)*pee1j1bj.vv_exp;
  pee1j1bj.tw_exp   = 0.637;  pee1j1bj.tw_stat   =   0.022; pee1j1bj.tw_syst_corr = oplus(0.1,lum_eRel,0.10)*pee1j1bj.tw_exp;    

  pee1j1bj.data = 8;

  pee1j1bj.sr_exp   = 0.0340; pee1j1bj.sr_stat = 0.0030; pee1j1bj.sr_syst = 0.0001; 

  pee1j1bj.ttotr_mc = 0.061;  pee1j1bj.ttotr_mc_stat  =   0.016;
  pee1j1bj.wj_mc    = 0.161;  pee1j1bj.wj_mc_stat     =   0.115;

  pee1j1bj.qcd_exp   = 0.0990; pee1j1bj.qcd_stat  =  0.0990; pee1j1bj.qcd_syst = pee1j1bj.qcd_exp;
  pee1j1bj.wjraw_exp = 0.7193; pee1j1bj.wjraw_stat=  0.4505;
  pee1j1bj.wjf_systFrac = 0.5;

  pee1j1bj.dy_mc  =    (0.254+0.059)*0.5;  pee1j1bj.dy_mc_stat = (0.104+0.042)*0.5;
  pee1j1bj.dy_exp =    (0.975+0.214)*0.5;  pee1j1bj.dy_stat  =   (0.565+0.176)*0.5; 
  pee1j1bj.dy_syst =   pee1j1bj.dy_exp*0.5;
  pee1j1bj.dy_roi = (0.1364+0.0300)*0.5; pee1j1bj.dy_roi_stat = (0.0566+0.0215)*0.5;

  pee1j1bj.sf_exp = (0.108*9.)*(0.108*9.) * 0.9231* 1.0126;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pee1j1bj.tt_AE_eRel = oplus(0.050,0.113); //use the same as FOM Njet=1
 
  xsecCalc_inStruct(pee1j1bj);


  std::cout<<"\n=======\n\t EE 2 jets with b-tags\n========="<<std::endl;
  TTxsecStruct pee2j1bj;
  pee2j1bj.lum_total = 36.1;
  pee2j1bj.channel = ee_ch;
  pee2j1bj.tt_exp   =15.285; pee2j1bj.tt_stat    =   0.250;
  pee2j1bj.dytt_exp = 0.176;  pee2j1bj.dytt_stat =   0.088; pee2j1bj.dytt_syst_corr = oplus(0.1,lum_eRel,0.25)*pee2j1bj.dytt_exp; 
  pee2j1bj.vv_exp   = 0.083;  pee2j1bj.vv_stat   =   0.010; pee2j1bj.vv_syst_corr = oplus(0.1,lum_eRel,0.25)*pee2j1bj.vv_exp;	 
  pee2j1bj.tw_exp   = 0.481;  pee2j1bj.tw_stat   =   0.019; pee2j1bj.tw_syst_corr = oplus(0.1,lum_eRel,0.10)*pee2j1bj.tw_exp;    

  pee2j1bj.data = 15;

  pee2j1bj.sr_exp  = 0.063; pee2j1bj.sr_stat = 0.009; pee2j1bj.sr_syst = 0.01455; 

  pee2j1bj.ttotr_mc = 0.342;  pee2j1bj.ttotr_mc_stat  =   0.037;
  pee2j1bj.wj_mc    = 0.000;  pee2j1bj.wj_mc_stat     =   0.111;

  pee2j1bj.qcd_exp   = 0.0000; pee2j1bj.qcd_stat  =  0.0900; pee2j1bj.qcd_syst = pee2j1bj.qcd_exp;
  pee2j1bj.wjraw_exp = 2.5821; pee2j1bj.wjraw_stat=  0.8366;
  pee2j1bj.wjf_systFrac = 0.5;

  pee2j1bj.dy_mc  =    (0.417+0.720)*0.5;  pee2j1bj.dy_mc_stat = (0.135+0.273)*0.5;
  pee2j1bj.dy_exp =    (0.589+0.929)*0.5;  pee2j1bj.dy_stat   =  (0.438+0.717)*0.5; 
  pee2j1bj.dy_syst =   pee2j1bj.dy_exp*0.5;
  pee2j1bj.dy_roi = (0.1120+0.1767)*0.5; pee2j1bj.dy_roi_stat = (0.0368+0.0686)*0.5;

  pee2j1bj.sf_exp = (0.108*9.)*(0.108*9.) * 0.9231* 1.0126;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pee2j1bj.tt_AE_eRel = 0.099; //mind the correlations
 
  xsecCalc_inStruct(pee2j1bj);


  std::cout<<"\n=======\n\t MM 1 jets with b-tag\n========="<<std::endl;
  TTxsecStruct pmm1j1bj;
  pmm1j1bj.lum_total = 36.1;
  pmm1j1bj.channel = mm_ch;
  pmm1j1bj.tt_exp   = 4.721; pmm1j1bj.tt_stat    =   0.139;
  pmm1j1bj.dytt_exp = 0.309;  pmm1j1bj.dytt_stat =   0.117; pmm1j1bj.dytt_syst_corr = oplus(0.1,lum_eRel,0.25)*pmm1j1bj.dytt_exp;
  pmm1j1bj.vv_exp   = 0.114;  pmm1j1bj.vv_stat   =   0.012; pmm1j1bj.vv_syst_corr = oplus(0.1,lum_eRel,0.25)*pmm1j1bj.vv_exp;	 
  pmm1j1bj.tw_exp   = 0.795;  pmm1j1bj.tw_stat   =   0.025; pmm1j1bj.tw_syst_corr = oplus(0.1,lum_eRel,0.10)*pmm1j1bj.tw_exp;    

  pmm1j1bj.data = 7;

  pmm1j1bj.sr_exp   = 0.0220; pmm1j1bj.sr_stat = 0.0023; pmm1j1bj.sr_syst = 0.0023; 

  pmm1j1bj.ttotr_mc = 0.004;  pmm1j1bj.ttotr_mc_stat  =   0.004;
  pmm1j1bj.wj_mc    = 0.000;  pmm1j1bj.wj_mc_stat     =   0.111;

  pmm1j1bj.qcd_exp   = 0.0000; pmm1j1bj.qcd_stat  =  0.0900; pmm1j1bj.qcd_syst = pmm1j1bj.qcd_exp;
  pmm1j1bj.wjraw_exp = 0.2892; pmm1j1bj.wjraw_stat=  0.2892;
  pmm1j1bj.wjf_systFrac = 0.75;

  pmm1j1bj.dy_mc  =    (0.683+0.953)*0.5;  pmm1j1bj.dy_mc_stat = (0.173+0.307)*0.5;
  pmm1j1bj.dy_exp =    (2.909+2.842)*0.5;  pmm1j1bj.dy_stat =    (1.212+1.314)*0.5; 
  pmm1j1bj.dy_syst =   pmm1j1bj.dy_exp*0.5;
  pmm1j1bj.dy_roi = (0.2688+0.2626)*0.5; pmm1j1bj.dy_roi_stat = (0.0686+0.0852)*0.5;

  pmm1j1bj.sf_exp = (0.108*9.)*(0.108*9.) * 0.9613* 1.0126*rescaleMuSF*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pmm1j1bj.tt_AE_eRel = oplus(0.050,0.112); //use the same as FOM Njet=1
 
  xsecCalc_inStruct(pmm1j1bj);


  std::cout<<"\n=======\n\t MM 2 jets with b-tags\n========="<<std::endl;
  TTxsecStruct pmm2j1bj;
  pmm2j1bj.lum_total = 36.1;
  pmm2j1bj.channel = mm_ch;
  pmm2j1bj.tt_exp   =18.527; pmm2j1bj.tt_stat    =   0.275;
  pmm2j1bj.dytt_exp = 0.250;  pmm2j1bj.dytt_stat =   0.103; pmm2j1bj.dytt_syst_corr = oplus(0.1,lum_eRel,0.25)*pmm2j1bj.dytt_exp;
  pmm2j1bj.vv_exp   = 0.082;  pmm2j1bj.vv_stat   =   0.009; pmm2j1bj.vv_syst_corr = oplus(0.1,lum_eRel,0.25)*pmm2j1bj.vv_exp;	 
  pmm2j1bj.tw_exp   = 0.562;  pmm2j1bj.tw_stat   =   0.021; pmm2j1bj.tw_syst_corr = oplus(0.1,lum_eRel,0.10)*pmm2j1bj.tw_exp;    

  pmm2j1bj.data = 24;

  pmm2j1bj.sr_exp  = 0.033; pmm2j1bj.sr_stat = 0.005; pmm2j1bj.sr_syst = 0.0060; 

  pmm2j1bj.ttotr_mc = 0.069;  pmm2j1bj.ttotr_mc_stat  =   0.017;
  pmm2j1bj.wj_mc    = 0.000;  pmm2j1bj.wj_mc_stat     =   0.111;

  pmm2j1bj.qcd_exp   = 0.0000; pmm2j1bj.qcd_stat  =  0.0900; pmm2j1bj.qcd_syst = pmm2j1bj.qcd_exp;
  pmm2j1bj.wjraw_exp = 1.2856; pmm2j1bj.wjraw_stat=  0.6478;
  pmm2j1bj.wjf_systFrac = 0.75;

  pmm2j1bj.dy_mc  =    (1.081+1.432)*0.5;  pmm2j1bj.dy_mc_stat = (0.210+0.387)*0.5;
  pmm2j1bj.dy_exp =    (3.251+3.004)*0.5;  pmm2j1bj.dy_stat =    (1.339+1.361)*0.5; 
  pmm2j1bj.dy_syst =   pmm2j1bj.dy_exp*0.5;
  pmm2j1bj.dy_roi = (0.2387+0.2206)*0.5; pmm2j1bj.dy_roi_stat = (0.0467+0.0600)*0.5;

  pmm2j1bj.sf_exp = (0.108*9.)*(0.108*9.) * 0.9613* 1.0126*rescaleMuSF*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pmm2j1bj.tt_AE_eRel = 0.098; //mind the correlations
 
  xsecCalc_inStruct(pmm2j1bj);


  std::cout<<"\n=======\n\t EM 1 jets with b-tag\n========="<<std::endl;
  TTxsecStruct pem1j1bj;
  pem1j1bj.lum_total = 36.1;
  pem1j1bj.channel = em_ch;
  pem1j1bj.tt_exp   =13.496; pem1j1bj.tt_stat    =   0.235;
  pem1j1bj.dytt_exp = 1.235;  pem1j1bj.dytt_stat =   0.233; pem1j1bj.dytt_syst_corr = oplus(0.1,lum_eRel,0.25)*pem1j1bj.dytt_exp;
  pem1j1bj.vv_exp   = 0.371;  pem1j1bj.vv_stat   =   0.022; pem1j1bj.vv_syst_corr = oplus(0.1,lum_eRel,0.25)*pem1j1bj.vv_exp;	 
  pem1j1bj.tw_exp   = 2.291;  pem1j1bj.tw_stat   =   0.042; pem1j1bj.tw_syst_corr = oplus(0.1,lum_eRel,0.10)*pem1j1bj.tw_exp;    

  pem1j1bj.data = 14;

  pem1j1bj.sr_exp   = 0.028; pem1j1bj.sr_stat = 0.002; pem1j1bj.sr_syst = 0.0012; 

  pem1j1bj.ttotr_mc = 0.159;  pem1j1bj.ttotr_mc_stat  =   0.025;
  pem1j1bj.wj_mc    = 0.111;  pem1j1bj.wj_mc_stat     =   0.111;

  pem1j1bj.qcd_exp   = 0.5974; pem1j1bj.qcd_stat  =  0.2278; pem1j1bj.qcd_syst = pem1j1bj.qcd_exp;
  pem1j1bj.wjraw_exp = 2.2868; pem1j1bj.wjraw_stat=  0.8257;
  pem1j1bj.wjf_systFrac = 0.5;

  pem1j1bj.dy_mc  =    0.085;  pem1j1bj.dy_mc_stat = 0.060;
  pem1j1bj.dy_exp =    0.0; pem1j1bj.dy_stat =     0.0; pem1j1bj.dy_syst =   pem1j1bj.dy_exp*0.5;
  pem1j1bj.dy_roi = 0;      pem1j1bj.dy_roi_stat = 0;

  pem1j1bj.sf_exp = (0.108*9.)*(0.108*9.) * 0.9444* 1.0126*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pem1j1bj.tt_AE_eRel = oplus(0.050, 0.138); //use the same as with cuts in FOM=1jet +btagging
 
  xsecCalc_inStruct(pem1j1bj);


  std::cout<<"\n=======\n\t EM 2 jets with b-tags\n========="<<std::endl;
  TTxsecStruct pem2j1bj;
  pem2j1bj.lum_total = 36.1;
  pem2j1bj.channel = em_ch;
  pem2j1bj.tt_exp   =53.825; pem2j1bj.tt_stat    =   0.468;
  pem2j1bj.dytt_exp = 0.647;  pem2j1bj.dytt_stat =   0.168; pem2j1bj.dytt_syst_corr = oplus(0.1,lum_eRel,0.25)*pem2j1bj.dytt_exp;
  pem2j1bj.vv_exp   = 0.226;  pem2j1bj.vv_stat   =   0.017; pem2j1bj.vv_syst_corr = oplus(0.1,lum_eRel,0.25)*pem2j1bj.vv_exp;	 
  pem2j1bj.tw_exp   = 1.563;  pem2j1bj.tw_stat   =   0.035; pem2j1bj.tw_syst_corr = oplus(0.1,lum_eRel,0.10)*pem2j1bj.tw_exp;    

  pem2j1bj.data = 51;

  pem2j1bj.sr_exp  = 0.048; pem2j1bj.sr_stat = 0.0052; pem2j1bj.sr_syst = 0.0079; 

  pem2j1bj.ttotr_mc = 0.856;  pem2j1bj.ttotr_mc_stat  =   0.059;
  pem2j1bj.wj_mc    = 0.111;  pem2j1bj.wj_mc_stat     =   0.111;

  pem2j1bj.qcd_exp   = 0.2191; pem2j1bj.qcd_stat  =  0.1611; pem2j1bj.qcd_syst = pem2j1bj.qcd_exp;
  pem2j1bj.wjraw_exp = 2.9343; pem2j1bj.wjraw_stat=  0.9498;
  pem2j1bj.wjf_systFrac = 0.5;

  pem2j1bj.dy_mc  =    0.085;  pem2j1bj.dy_mc_stat = 0.060;
  pem2j1bj.dy_exp =    0.0;  pem2j1bj.dy_stat =    0.0; pem2j1bj.dy_syst =   pem2j1bj.dy_exp*0.5;
  pem2j1bj.dy_roi = 0;       pem2j1bj.dy_roi_stat = 0;

  pem2j1bj.sf_exp = (0.108*9.)*(0.108*9.) * 0.9444* 1.0126*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pem2j1bj.tt_AE_eRel = 0.083; //mind the correlations
 
  xsecCalc_inStruct(pem2j1bj);



  std::cout<<"\n\n\n\n"<<std::endl;
  std::cout<<"#=============================================================#"<<std::endl;
  std::cout<<"# It's all PF in 1 jet FOM cuts                               #"<<std::endl;
  std::cout<<"#-------------------------------------------------------------#"<<std::endl;

  std::cout<<"\n=======\n\t EE 1 jet with FOM cuts\n========="<<std::endl;
  TTxsecStruct pee1jFOM;
  pee1jFOM.lum_total = 36.1;
  pee1jFOM.channel = ee_ch;
  pee1jFOM.tt_exp   = 3.930;  pee1jFOM.tt_stat   =   0.127;
  pee1jFOM.dytt_exp = 0.662;  pee1jFOM.dytt_stat =   0.171; pee1jFOM.dytt_syst_corr = oplus(0.1,lum_eRel)*pee1jFOM.dytt_exp;
  pee1jFOM.vv_exp   = 0.549;  pee1jFOM.vv_stat   =   0.026; pee1jFOM.vv_syst_corr   = oplus(0.1,lum_eRel)*pee1jFOM.vv_exp;
  pee1jFOM.tw_exp   = 0.648;  pee1jFOM.tw_stat   =   0.022; pee1jFOM.tw_syst_corr   = oplus(0.1,lum_eRel)*pee1jFOM.tw_exp;    

  pee1jFOM.data = 8;

  pee1jFOM.sr_exp   = 0.0340; pee1jFOM.sr_stat = 0.0030; pee1jFOM.sr_syst = 0.0001; 

  pee1jFOM.ttotr_mc = 0.077;  pee1jFOM.ttotr_mc_stat  =   0.018;
  pee1jFOM.wj_mc    = 0.945;  pee1jFOM.wj_mc_stat     =   0.286;

  pee1jFOM.qcd_exp   = 0.2862; pee1jFOM.qcd_stat  =  0.1654; pee1jFOM.qcd_syst = pee1jFOM.qcd_exp;
  pee1jFOM.wjraw_exp = 0.9267; pee1jFOM.wjraw_stat=  0.4827;
  pee1jFOM.wjf_systFrac = 0.5;

  pee1jFOM.dy_mc  =    (0.085+0.220)*0.5;  pee1jFOM.dy_mc_stat = (0.060+0.156)*0.5;
  pee1jFOM.dy_exp =    (0.141+0.299)*0.5;  pee1jFOM.dy_stat  =   (0.251+0.535)*0.5; 
  pee1jFOM.dy_syst =   pee1jFOM.dy_exp*0.5;
  pee1jFOM.dy_roi =    (0.118+0.250)*0.5; pee1jFOM.dy_roi_stat = (0.084+0.182)*0.5;

  pee1jFOM.sf_exp = (0.108*9.)*(0.108*9.) * 0.9231* 1.0126;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pee1jFOM.tt_AE_eRel = 0.113; //
 
  xsecCalc_inStruct(pee1jFOM);


  std::cout<<"\n=======\n\t MM 1 jet with FOM cuts\n========="<<std::endl;
  TTxsecStruct pmm1jFOM;
  pmm1jFOM.lum_total = 36.1;
  pmm1jFOM.channel = mm_ch;
  pmm1jFOM.tt_exp   = 4.778; pmm1jFOM.tt_stat    =   0.140;
  pmm1jFOM.dytt_exp = 0.574;  pmm1jFOM.dytt_stat =   0.159; pmm1jFOM.dytt_syst_corr = oplus(0.1,lum_eRel)*pmm1jFOM.dytt_exp;
  pmm1jFOM.vv_exp   = 0.712;  pmm1jFOM.vv_stat   =   0.030; pmm1jFOM.vv_syst_corr   = oplus(0.1,lum_eRel)*pmm1jFOM.vv_exp;	 
  pmm1jFOM.tw_exp   = 0.743;  pmm1jFOM.tw_stat   =   0.024; pmm1jFOM.tw_syst_corr   = oplus(0.1,lum_eRel)*pmm1jFOM.tw_exp;    

  pmm1jFOM.data = 10;

  pmm1jFOM.sr_exp   = 0.0220; pmm1jFOM.sr_stat = 0.0023; pmm1jFOM.sr_syst = 0.0023; 

  pmm1jFOM.ttotr_mc = 0.008;  pmm1jFOM.ttotr_mc_stat  =   0.006;
  pmm1jFOM.wj_mc    = 0.111;  pmm1jFOM.wj_mc_stat     =   0.111;

  pmm1jFOM.qcd_exp   = 0.0000; pmm1jFOM.qcd_stat  =  0.0900; pmm1jFOM.qcd_syst = pmm1jFOM.qcd_exp;
  pmm1jFOM.wjraw_exp = 0.2892; pmm1jFOM.wjraw_stat=  0.2892;
  pmm1jFOM.wjf_systFrac = 0.75;

  pmm1jFOM.dy_mc  =    (1.451+1.32)*0.5;  pmm1jFOM.dy_mc_stat = (0.246+0.38)*0.5;
  pmm1jFOM.dy_exp =    (6.068+3.825)*0.5;  pmm1jFOM.dy_stat =    (3.841+2.59)*0.5; 
  pmm1jFOM.dy_syst =   pmm1jFOM.dy_exp*0.5;
  pmm1jFOM.dy_roi =   (1.2692+0.800)*0.5; pmm1jFOM.dy_roi_stat = (0.255039+0.249)*0.5;

  pmm1jFOM.sf_exp = (0.108*9.)*(0.108*9.) * 0.9613* 1.0126*rescaleMuSF*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pmm1jFOM.tt_AE_eRel = 0.112; //JES=3.6 from PK (should I be using Pedram's value ~8.7? since it's larger)
 
  xsecCalc_inStruct(pmm1jFOM);




  std::cout<<"\n=======\n\t EM 1 jet with FOM cuts\n========="<<std::endl;
  TTxsecStruct pem1jFOM;
  pem1jFOM.lum_total = 36.1;
  pem1jFOM.channel = em_ch;
  pem1jFOM.tt_exp   =11.999; pem1jFOM.tt_stat    =   0.221;
  pem1jFOM.dytt_exp = 0.000;  pem1jFOM.dytt_stat =   0.044; pem1jFOM.dytt_syst_corr = oplus(0.1,lum_eRel)*pem1jFOM.dytt_exp;
  pem1jFOM.vv_exp   = 1.878;  pem1jFOM.vv_stat   =   0.049; pem1jFOM.vv_syst_corr   = oplus(0.1,lum_eRel)*pem1jFOM.vv_exp;	 
  pem1jFOM.tw_exp   = 2.036;  pem1jFOM.tw_stat   =   0.040; pem1jFOM.tw_syst_corr   = oplus(0.1,lum_eRel)*pem1jFOM.tw_exp;    

  pem1jFOM.data = 18;

  pem1jFOM.sr_exp   = 0.028; pem1jFOM.sr_stat = 0.002; pem1jFOM.sr_syst = 0.0012; 

  pem1jFOM.ttotr_mc = 0.106;  pem1jFOM.ttotr_mc_stat  =   0.021;
  pem1jFOM.wj_mc    = 0.869;  pem1jFOM.wj_mc_stat     =   0.308;

  pem1jFOM.qcd_exp   = 0.3223; pem1jFOM.qcd_stat  =  0.1613; pem1jFOM.qcd_syst = pem1jFOM.qcd_exp;
  pem1jFOM.wjraw_exp = 1.6467; pem1jFOM.wjraw_stat=  0.6898;
  pem1jFOM.wjf_systFrac = 0.5;

  pem1jFOM.dy_mc  =    0.127;  pem1jFOM.dy_mc_stat = 0.073;
  pem1jFOM.dy_exp =    0.0; pem1jFOM.dy_stat =     0.0; pem1jFOM.dy_syst =   pem1jFOM.dy_exp*0.5;
  pem1jFOM.dy_roi = 0;      pem1jFOM.dy_roi_stat = 0;

  pem1jFOM.sf_exp = (0.108*9.)*(0.108*9.) * 0.9444* 1.0126*rescaleMuSF;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pem1jFOM.tt_AE_eRel = 0.138; //
 
  xsecCalc_inStruct(pem1jFOM);





  std::cout<<"\n\n\n\n"<<std::endl;
  std::cout<<"#=============================================================#"<<std::endl;
  std::cout<<"# It's all PF in 1 jet FOM with numbers from LIP              #"<<std::endl;
  std::cout<<"#-------------------------------------------------------------#"<<std::endl;

  std::cout<<"\n=======\n\t EE 1 jet with FOMLIP cuts\n========="<<std::endl;
  TTxsecStruct pee1jFOMLIP;
  pee1jFOMLIP.lum_total = 36.1;
  pee1jFOMLIP.channel = ee_ch;
  pee1jFOMLIP.tt_exp   = 3.4794;  pee1jFOMLIP.tt_stat   =   0.127;
  pee1jFOMLIP.dytt_exp = 0.5000;  pee1jFOMLIP.dytt_stat =   0.1200; pee1jFOMLIP.dytt_syst_corr =  oplus(0.1,lum_eRel)*pee1jFOMLIP.dytt_exp;
  pee1jFOMLIP.vv_exp   = 0.5400;  pee1jFOMLIP.vv_stat   =   0.0277; pee1jFOMLIP.vv_syst_corr   = oplus(0.1,lum_eRel)*pee1jFOMLIP.vv_exp;
  pee1jFOMLIP.tw_exp   = 0.5700;  pee1jFOMLIP.tw_stat   =   0.0177; pee1jFOMLIP.tw_syst_corr   = oplus(0.1,lum_eRel)*pee1jFOMLIP.tw_exp;

  pee1jFOMLIP.data = 8;

  pee1jFOMLIP.sr_exp   = 0.0; pee1jFOMLIP.sr_stat = 0.0; pee1jFOMLIP.sr_syst = 0.0; 

  pee1jFOMLIP.ttotr_mc = 0.077;  pee1jFOMLIP.ttotr_mc_stat  =   0.018;
  pee1jFOMLIP.wj_mc    = 0.945;  pee1jFOMLIP.wj_mc_stat     =   0.286;

  pee1jFOMLIP.qcd_exp   = 0.000; pee1jFOMLIP.qcd_stat  =  0.000; pee1jFOMLIP.qcd_syst = pee1jFOMLIP.qcd_exp;
  pee1jFOMLIP.wjraw_exp = 0.29; pee1jFOMLIP.wjraw_stat=  0.49;
  pee1jFOMLIP.wjf_systFrac = 0.0;

  pee1jFOMLIP.dy_mc  =    (0.085+0.220)*0.5;  pee1jFOMLIP.dy_mc_stat = (0.060+0.156)*0.5;
  pee1jFOMLIP.dy_exp =    0.16;  pee1jFOMLIP.dy_stat  =   0.29; 
  pee1jFOMLIP.dy_syst =   0;
  pee1jFOMLIP.dy_roi =    (0.118+0.250)*0.5; pee1jFOMLIP.dy_roi_stat = (0.084+0.182)*0.5;

  pee1jFOMLIP.sf_exp = 1.0;//(0.108*9.)*(0.108*9.) * 0.9231* 1.0126;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pee1jFOMLIP.tt_AE_eRel = 0.113; //
 
  xsecCalc_inStruct(pee1jFOMLIP);


  std::cout<<"\n=======\n\t MM 1 jet with FOMLIP cuts\n========="<<std::endl;
  TTxsecStruct pmm1jFOMLIP;
  pmm1jFOMLIP.lum_total = 36.1;
  pmm1jFOMLIP.channel = mm_ch;
  pmm1jFOMLIP.tt_exp   = 4.4147; pmm1jFOMLIP.tt_stat    =   0.140;
  pmm1jFOMLIP.dytt_exp = 0.5300; pmm1jFOMLIP.dytt_stat =   0.1190; pmm1jFOMLIP.dytt_syst_corr = oplus(0.1,lum_eRel)*pmm1jFOMLIP.dytt_exp;
  pmm1jFOMLIP.vv_exp   = 0.6400; pmm1jFOMLIP.vv_stat   =   0.0298; pmm1jFOMLIP.vv_syst_corr   = oplus(0.1,lum_eRel)*pmm1jFOMLIP.vv_exp;	 
  pmm1jFOMLIP.tw_exp   = 0.6800; pmm1jFOMLIP.tw_stat   =   0.0184; pmm1jFOMLIP.tw_syst_corr   = oplus(0.1,lum_eRel)*pmm1jFOMLIP.tw_exp;    

  pmm1jFOMLIP.data = 10;

  pmm1jFOMLIP.sr_exp   = 0.0; pmm1jFOMLIP.sr_stat = 0.0; pmm1jFOMLIP.sr_syst = 0.0; 

  pmm1jFOMLIP.ttotr_mc = 0.008;  pmm1jFOMLIP.ttotr_mc_stat  =   0.006;
  pmm1jFOMLIP.wj_mc    = 0.111;  pmm1jFOMLIP.wj_mc_stat     =   0.111;

  pmm1jFOMLIP.qcd_exp   = 0.0000; pmm1jFOMLIP.qcd_stat  =  0.00; pmm1jFOMLIP.qcd_syst = pmm1jFOMLIP.qcd_exp;
  pmm1jFOMLIP.wjraw_exp = 0.07; pmm1jFOMLIP.wjraw_stat=  0.38;
  pmm1jFOMLIP.wjf_systFrac = 0.0;

  pmm1jFOMLIP.dy_mc  =    (1.451+1.32)*0.5;  pmm1jFOMLIP.dy_mc_stat = (0.246+0.38)*0.5;
  pmm1jFOMLIP.dy_exp =    5.19;  pmm1jFOMLIP.dy_stat =    4.28; 
  pmm1jFOMLIP.dy_syst =   0;
  pmm1jFOMLIP.dy_roi =   (1.2692+0.800)*0.5; pmm1jFOMLIP.dy_roi_stat = (0.255039+0.249)*0.5;

  pmm1jFOMLIP.sf_exp = 1.0;//(0.108*9.)*(0.108*9.) * 0.9613* 1.0126;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pmm1jFOMLIP.tt_AE_eRel = 0.112; //JES=3.6 from PK (should I be using Pedram's value ~8.7? since it's larger)
 
  xsecCalc_inStruct(pmm1jFOMLIP);




  std::cout<<"\n=======\n\t EM 1 jet with FOMLIP cuts\n========="<<std::endl;
  TTxsecStruct pem1jFOMLIP;
  pem1jFOMLIP.lum_total = 36.1;
  pem1jFOMLIP.channel = em_ch;
  pem1jFOMLIP.tt_exp   =10.9028; pem1jFOMLIP.tt_stat    =   0.221;
  pem1jFOMLIP.dytt_exp = 0.0500; pem1jFOMLIP.dytt_stat =   0.0500; pem1jFOMLIP.dytt_syst_corr = oplus(0.1,lum_eRel)*pem1jFOMLIP.dytt_exp;
  pem1jFOMLIP.vv_exp   = 1.7700; pem1jFOMLIP.vv_stat   =   0.0472; pem1jFOMLIP.vv_syst_corr   = oplus(0.1,lum_eRel)*pem1jFOMLIP.vv_exp;	 
  pem1jFOMLIP.tw_exp   = 1.8600; pem1jFOMLIP.tw_stat   =   0.0362; pem1jFOMLIP.tw_syst_corr   = oplus(0.1,lum_eRel)*pem1jFOMLIP.tw_exp;    

  pem1jFOMLIP.data = 18;

  pem1jFOMLIP.sr_exp   = 0.0; pem1jFOMLIP.sr_stat = 0.0; pem1jFOMLIP.sr_syst = 0.0; 

  pem1jFOMLIP.ttotr_mc = 0.106;  pem1jFOMLIP.ttotr_mc_stat  =   0.021;
  pem1jFOMLIP.wj_mc    = 0.869;  pem1jFOMLIP.wj_mc_stat     =   0.308;

  pem1jFOMLIP.qcd_exp   = 0.0; pem1jFOMLIP.qcd_stat  =  0.0; pem1jFOMLIP.qcd_syst = pem1jFOMLIP.qcd_exp;
  pem1jFOMLIP.wjraw_exp = 1.25; pem1jFOMLIP.wjraw_stat=  1.28;
  pem1jFOMLIP.wjf_systFrac = 0.0;

  pem1jFOMLIP.dy_mc  =    0.127;  pem1jFOMLIP.dy_mc_stat = 0.073;
  pem1jFOMLIP.dy_exp =    0.0; pem1jFOMLIP.dy_stat =     0.0; pem1jFOMLIP.dy_syst =   pem1jFOMLIP.dy_exp*0.5;
  pem1jFOMLIP.dy_roi = 0;      pem1jFOMLIP.dy_roi_stat = 0;

  pem1jFOMLIP.sf_exp = 1.;//(0.108*9.)*(0.108*9.) * 0.9444* 1.0126;//(0.108*9.) is for Madgraph ! 1.0126 is PU
  pem1jFOMLIP.tt_AE_eRel = 0.138; //
 
  xsecCalc_inStruct(pem1jFOMLIP);


  TTxsecStruct pee;
  TTxsecStruct pmm;
  TTxsecStruct pem;




  //combine
  std::cout<<"====================================================="<<std::endl;
  std::cout<<"======       COMBINATIONS without b-tags     ========"<<std::endl; 
  std::cout<<"====================================================="<<std::endl;
  std::cout<<"EEMM 2jets no tags"<<std::endl;
  TTxsecStruct peemm2j(pee2j);
  peemm2j.lum_total = 36.1;
  peemm2j.simpleAdd(pmm2j);
  xsecCalc_inStruct(peemm2j, false);

  std::cout<<"EE-MM-EM 2 jets no tags"<<std::endl;
  TTxsecStruct peemmem2j(peemm2j);
  peemmem2j.lum_total = 36.1;
  peemmem2j.simpleAdd(pem2j);
  xsecCalc_inStruct(peemmem2j, false);

  std::cout<<"EE-MM-EM 2 jets no tags, uncorr backgrounds"<<std::endl;
  TTxsecStruct peemmem2jUC(peemm2j);
  peemmem2jUC.lum_total = 36.1;
  peemmem2jUC.simpleAdd(pem2j, false);
  xsecCalc_inStruct(peemmem2jUC, false);

  std::cout<<"===================================================="<<std::endl;
  std::cout<<"======       COMBINATIONS  with b-tags      ========"<<std::endl; 
  std::cout<<"===================================================="<<std::endl;
  std::cout<<"EE-MM 2jets with tags"<<std::endl;
  TTxsecStruct peemm2j1bj(pee2j1bj);
  peemm2j1bj.lum_total = 36.1;
  peemm2j1bj.simpleAdd(pmm2j1bj);
  xsecCalc_inStruct(peemm2j1bj, false);

  std::cout<<"EE-MM-EM 2 jets with tags"<<std::endl;
  TTxsecStruct peemmem2j1bj(peemm2j1bj);
  peemmem2j1bj.lum_total = 36.1;
  peemmem2j1bj.simpleAdd(pem2j1bj);
  xsecCalc_inStruct(peemmem2j1bj, false);

  std::cout<<"EE-MM 1jet with a tag"<<std::endl;
  TTxsecStruct peemm1j1bj(pee1j1bj);
  peemm1j1bj.lum_total = 36.1;
  peemm1j1bj.simpleAdd(pmm1j1bj);
  xsecCalc_inStruct(peemm1j1bj, false);

  std::cout<<"EE-MM-EM 1 jet with a tag"<<std::endl;
  TTxsecStruct peemmem1j1bj(peemm1j1bj);
  peemmem1j1bj.lum_total = 36.1;
  peemmem1j1bj.simpleAdd(pem1j1bj);
  xsecCalc_inStruct(peemmem1j1bj, false);



  std::cout<<"\n\n\nFINAL\n(EE-MM with tag)-(EM no tag) 2 jets"<<std::endl;
  TTxsecStruct peemmem2jFinal(peemm2j1bj);
  peemmem2jFinal.lum_total = 36.1;
  peemmem2jFinal.simpleAdd(pem2j); 
  // now a hack: start from the untagged uncertainty
  peemmem2jFinal.tt_AE_eRel = (peemm2j.tt_AE_eRel*peemm2j1bj.tt_exp 
			       + pem2j.tt_AE_eRel*pem2j.tt_exp)/(peemm2j1bj.tt_exp+pem2j.tt_exp);
  peemmem2jFinal.tt_AE_eRel = oplus(peemmem2jFinal.tt_AE_eRel, 0.050*peemm2j1bj.tt_exp/(peemm2j1bj.tt_exp+pem2j.tt_exp));
  xsecCalc_inStruct(peemmem2jFinal, false);


  std::cout<<"\n\n\nFINAL\n(EE-MM with tag)-(EM no tag) 2 jets and EE-MM-EM tag 1 jet"<<std::endl;
  TTxsecStruct peemmem1n2jFinal(peemmem2jFinal);
  peemmem1n2jFinal.lum_total = 36.1;
  peemmem1n2jFinal.channel = other_ch;
  peemmem1n2jFinal.simpleAdd(peemmem1j1bj); 
  xsecCalc_inStruct(peemmem1n2jFinal, false);

  return;


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

void xsecCalc_comb_pass6_normLum(bool forAN410 = true){
  double mean_a = 165.926;//160.314;
  double k_b    = 1./0.9848;
  if (forAN410){
    mean_a = 159.791;
    k_b    = 1./0.9838;
  }
  double mean_b = mean_a*k_b;

  double relStat = 17.5814/mean_a;
  double relSyst = 13.7953/mean_a;
  if (forAN410){
    relStat = 0.117599;
    relSyst = 0.08432;
  }
  std::cout<<"Val b: "<<mean_b<<" \t\\pm "<<relStat*mean_b<<" \t\\pm "<<relSyst*mean_b;
  if (forAN410) std::cout <<" \t\\pm "<<0.052*mean_b<< "(= \t\\pm "<<oplus(relStat*mean_b,relSyst*mean_b);
  else std::cout <<" \t\\pm "<<0.05*mean_b <<"(= \t\\pm "<<oplus(relStat*mean_b,relSyst*mean_b);
  std::cout<<std::endl;
  
  double relComb = oplus(relStat, relSyst);
  double mean_[2] = {mean_a, mean_b};

  double sig_a = oplus(relComb* mean_a, 0.11*mean_a);
  double sig_aa = sig_a*sig_a;
  double sig_b = oplus(relComb* mean_a*k_b, 0.05*mean_a*k_b); 
  if (forAN410){
    sig_b = oplus(relComb* mean_a*k_b, 0.052*mean_a*k_b); 
  }
  double sig_bb = sig_b*sig_b;
  double sig_ab = k_b*relComb* mean_a*relComb* mean_a;
  double det = (sig_aa*sig_bb - sig_ab*sig_ab);
  double h_ab[2][2] = {
    {sig_bb/det,    -sig_ab/det},
    {-sig_ab/det,   sig_aa/det}
  };

  double mean = 0;
  double sumH = 0;
  for (int i =0;i<2; ++i){
    for(int j = 0; j<2; ++j){
      sumH += h_ab[i][j];
      mean += h_ab[i][j]*mean_[j];
    }
  }
  mean /= sumH;

  double sigma = sqrt(1./sumH);
  double sigma_normC = sqrt(sigma*sigma -  mean*relComb*mean*relComb);

  std::cout<<mean<<" \\pm "<<sigma
	   <<"\t = "<<mean<<" \\pm "<<mean*relStat<<"(stat) \\pm "<<mean*relSyst<<"(syst) \\pm "<< sigma_normC<<"(norm)" 
	   <<"\t = "<<mean<<" \\pm "<<mean*relComb<<"(stat+syst) \\pm "<<sigma_normC<<std::endl;
  
}
