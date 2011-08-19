void dimuonTrigEff(TH1F* hSF = 0){
  double frac_b   = 0.15;   double fracE_b  =  0.05;
  double eff_a_mc = 0.949; double effE_a_mc = 0.0003;
  double eff_b_mc = 0.419; double effE_b_mc = 0.003;

  double eff_a_9 = 0.880; double effE_a_9 = 0.008;
  double eff_b_9 = 0.37;  double effE_b_9 = 0.06;

  double eff_a_11 = 0.912; double effE_a_11 = 0.005;
  double eff_b_11 = 0.39;  double effE_b_11 = 0.04;

  double eff_a_15 = 0.922; double effE_a_15 = 0.002;
  double eff_b_15 = 0.43;  double effE_b_15 = 0.02;

  double eff_aa_mc = 1. - (1.-eff_a_mc)*(1.-eff_a_mc);// 2*eff_a_mc - eff_a_mc*eff_a_mc;
  double eff_ab_mc = 1. - (1.-eff_a_mc)*(1.-eff_b_mc);// eff_a_mc + eff_b_mc - eff_a_mc*eff_b_mc;
  double effE_aa_mc = 2.*(1. - eff_a_mc)*effE_a_mc;
  double effE_ab_mc = sqrt((1.-eff_a_mc)*(1.-eff_a_mc)*effE_a_mc*effE_a_mc
			   +(1.-eff_b_mc)*(1.-eff_b_mc)*effE_b_mc*effE_b_mc
			   );

  double eff_aa_9 = 1. - (1.-eff_a_9)*(1.-eff_a_9);
  double eff_ab_9 = 1. - (1.-eff_a_9)*(1.-eff_b_9);
  double effE_aa_9 = 2.*(1. - eff_a_9)*effE_a_9;
  double effE_ab_9 = sqrt((1.-eff_a_9)*(1.-eff_a_9)*effE_a_9*effE_a_9
			   +(1.-eff_b_9)*(1.-eff_b_9)*effE_b_9*effE_b_9
			   );

  double eff_aa_11 = 1. - (1.-eff_a_11)*(1.-eff_a_11);
  double eff_ab_11 = 1. - (1.-eff_a_11)*(1.-eff_b_11);
  double effE_aa_11 = 2.*(1. - eff_a_11)*effE_a_11;
  double effE_ab_11 = sqrt((1.-eff_a_11)*(1.-eff_a_11)*effE_a_11*effE_a_11
			   +(1.-eff_b_11)*(1.-eff_b_11)*effE_b_11*effE_b_11
			   );

  double eff_aa_15 = 1. - (1.-eff_a_15)*(1.-eff_a_15);
  double eff_ab_15 = 1. - (1.-eff_a_15)*(1.-eff_b_15);
  double effE_aa_15 = 2.*(1. - eff_a_15)*effE_a_15;
  double effE_ab_15 = sqrt((1.-eff_a_15)*(1.-eff_a_15)*effE_a_15*effE_a_15
			   +(1.-eff_b_15)*(1.-eff_b_15)*effE_b_15*effE_b_15
			   );

  double sf_aa_9 = eff_aa_9/eff_aa_mc;
  double sf_ab_9 = eff_ab_9/eff_ab_mc;

  double sf_aa_11 = eff_aa_11/eff_aa_mc;
  double sf_ab_11 = eff_ab_11/eff_ab_mc;

  double sf_aa_15 = eff_aa_15/eff_aa_mc;
  double sf_ab_15 = eff_ab_15/eff_ab_mc;

  std::cout<<sf_aa_9<<"\t"<<sf_ab_9<<std::endl;
  std::cout<<sf_aa_11<<"\t"<<sf_ab_11<<std::endl;
  std::cout<<sf_aa_15<<"\t"<<sf_ab_15<<std::endl;

  //this is fixed
  double eff_mc = eff_aa_mc*(1.-frac_b) + eff_ab_mc*frac_b;

  //this is the central value
  double eff_9 = eff_aa_9 + (eff_ab_9 - eff_aa_9)*frac_b;//eff_aa_9*(1.-frac_b) + eff_ab_9*frac_b;
  //the uncertainty. The frac is varied as well. The meaning of it is this fraction may differ
  // in data compared to simulation for a given particular signal
  double effEstat_9 = sqrt((1.-frac_b)*(1.-frac_b)* effE_aa_9*effE_aa_9
			   + frac_b   *   frac_b  * effE_ab_9*effE_ab_9);
  double effEsyst_9 = (eff_ab_9 - eff_aa_9)*fracE_b;
  double effE_9 = sqrt(effEstat_9*effEstat_9 + effEsyst_9*effEsyst_9);

  //this is the central value
  double eff_11 = eff_aa_11 + (eff_ab_11 - eff_aa_11)*frac_b;//eff_aa_11*(1.-frac_b) + eff_ab_11*frac_b;
  //the uncertainty. The frac is varied as well. The meaning of it is this fraction may differ
  // in data compared to simulation for a given particular signal
  double effEstat_11 = sqrt((1.-frac_b)*(1.-frac_b)* effE_aa_11*effE_aa_11
			   + frac_b   *   frac_b  * effE_ab_11*effE_ab_11);
  double effEsyst_11 = (eff_ab_11 - eff_aa_11)*fracE_b;
  double effE_11 = sqrt(effEstat_11*effEstat_11 + effEsyst_11*effEsyst_11);

  //this is the central value
  double eff_15 = eff_aa_15 + (eff_ab_15 - eff_aa_15)*frac_b;//eff_aa_15*(1.-frac_b) + eff_ab_15*frac_b;
  //the uncertainty. The frac is varied as well. The meaning of it is this fraction may differ
  // in data compared to simulation for a given particular signal
  double effEstat_15 = sqrt((1.-frac_b)*(1.-frac_b)* effE_aa_15*effE_aa_15
			   + frac_b   *   frac_b  * effE_ab_15*effE_ab_15);
  double effEsyst_15 = (eff_ab_15 - eff_aa_15)*fracE_b;
  double effE_15 = sqrt(effEstat_15*effEstat_15 + effEsyst_15*effEsyst_15);

  double lumF_9 = 3.1/36.1;
  double lumF_11 = 5/36.1;
  double lumF_15 = 28/36.1;
    
  double eff_av = lumF_9*eff_9 + lumF_11*eff_11 + lumF_15*eff_15;
  double effEstat_av = sqrt(lumF_9*lumF_9 * effEstat_9*effEstat_9
			    + lumF_11*lumF_11 * effEstat_11*effEstat_11
			    + lumF_15*lumF_15 * effEstat_15*effEstat_15);
  double effEsyst_av = fabs(lumF_9*effEsyst_9 + lumF_11*effEsyst_11 + lumF_15*effEsyst_15);
  double effE_av = sqrt(effEstat_av*effEstat_av + effEsyst_av*effEsyst_av);

  std::cout<<"MC\t"<<eff_mc<<"\t  9: "<<eff_9<<"\t \\pm "<<effE_9<<std::endl;
  std::cout<<"MC\t"<<eff_mc<<"\t 11: "<<eff_11<<"\t \\pm "<<effE_11<<std::endl;
  std::cout<<"MC\t"<<eff_mc<<"\t 15: "<<eff_15<<"\t \\pm "<<effE_15<<std::endl;
  std::cout<<"MC\t"<<eff_mc<<"\tAve: "<<eff_av
	   <<"\t \\pm "<<effEstat_av<<"\t \\pm "<<effEsyst_av
	   <<"\t = "<<eff_av<<"\t \\pm "<<effE_av<<std::endl;
  std::cout<<"SF(data/MC)\t  9: "<<eff_9/eff_mc<<"\t \\pm "<<effE_9/eff_mc<<std::endl;
  std::cout<<"SF(data/MC)\t 11: "<<eff_11/eff_mc<<"\t \\pm "<<effE_11/eff_mc<<std::endl;
  std::cout<<"SF(data/MC)\t 15: "<<eff_15/eff_mc<<"\t \\pm "<<effE_15/eff_mc<<std::endl;
  std::cout<<"SF(data/MC)\tAve: "<<eff_av/eff_mc
	   <<"\t \\pm "<<effEstat_av/eff_mc<<"\t \\pm "<<effEsyst_av/eff_mc
	   <<"\t = "<<eff_av/eff_mc<<"\t \\pm "<<effE_av/eff_mc<<std::endl;
			     


}


void dieleTriggerEff(){
  double eff[5] =  {0.983, 0.991, 0.989, 0.963, 0.961 }; 
  double effE[5] = {0.003, 0.005, 0.002, 0.003, 0.002 }; 
  double lum[5]  = {  2.4,   0.7,   5.0,   9.5, 18.5  };

  double effLoose[5] = {0.983, 0.983, 0.992, 0.992, 0.992};

  double lumSum = 0;
  double eff_av = 0;
  double effE_av = 0;
  double effDD = 0;
  for (int i = 0; i< 5; ++i){
    lumSum+= lum[i];
    eff_av+= (1.-(1.-eff[i])*(1.-eff[i]))*lum[i];
    effE_av+= 4.*(1.-eff[i])*(1.-eff[i])*lum[i]*lum[i]*effE[i]*effE[i];
    effDD += (effLoose[i]-eff[i])*(effLoose[i]-eff[i])*lum[i];
  }

  eff_av = eff_av/lumSum;
  effE_av = sqrt(effE_av)/lumSum;
  effDD = effDD/lumSum;
  std::cout<<"Ave: "<<eff_av<<"\t (+ "<< effDD<<" = "<< eff_av+effDD <<" )\t \\pm "<<effE_av<<std::endl;
}

void emuTriggerEff(){
  double eff_ele_mc = 0.970;
  double eff_ele[5] =  {0.983, 0.991, 0.989, 0.963, 0.961 }; 
  double effE_ele[5] = {0.003, 0.005, 0.002, 0.003, 0.002 }; 
  double lum[5]  = {  2.4,   0.7,   5.0,   9.5, 18.5  };
  double lumSum = 0; 
  double lumFrac[5];
  for (int i =0; i< 5; ++i) lumSum+=lum[i];
  for (int i =0; i< 5; ++i) lumFrac[i] = lum[i]/lumSum;

  double frac_b_mu   = 0.10;   double fracE_b_mu  =  0.05;

  double eff_a_mu_mc = 0.949; double effE_a_mu_mc = 0.0003;
  double eff_b_mu_mc = 0.419; double effE_b_mu_mc = 0.003;

  double eff_a_mu[5] =  {0.880, 0.880, 0.912, 0.922, 0.922 }; 
  double effE_a_mu[5] = {0.008, 0.008, 0.005, 0.002, 0.002 };

  double eff_b_mu[5]  = {0.37, 0.37, 0.39, 0.43, 0.43 };
  double effE_b_mu[5] = {0.06, 0.06, 0.04, 0.02, 0.02 };

  double corr_mu[5][5] = {{ 1, 1, 0, 0, 0},
			  { 1, 1, 0, 0, 0},
			  { 0, 0, 1, 0, 0},
			  { 0, 0, 0, 1, 1},
			  { 0, 0, 0, 1, 1}};

  double eff_emu[5];
  double eff_Demu[5];
  double eff_eDamu[5];
  double eff_eDbmu[5];
  double eff_emuDfrac[5];
  double effCov_emu[5][5];
  
  for (int i = 0; i < 5; ++i){
    eff_emu[i] = 1. - (1.-eff_ele[i])*((1.-frac_b_mu)*(1.-eff_a_mu[i])
				       + frac_b_mu*(1.-eff_b_mu[i]));
    eff_Demu[i] = ((1.-frac_b_mu)*(1.-eff_a_mu[i])
		   + frac_b_mu*(1.-eff_b_mu[i]));
    eff_eDamu[i] = (1.-eff_ele[i])*(1.-frac_b_mu);
    eff_eDbmu[i] = (1.-eff_ele[i])*frac_b_mu;
    eff_emuDfrac[i] = -(1.-eff_ele[i])*(eff_a_mu[i] -eff_b_mu[i]);
  }

  
  for (int i =0; i<5; ++i){
    for (int j = 0; j< 5; ++j){
      effCov_emu[i][j] = 0;
      if (i==j) effCov_emu[i][j] += eff_Demu[i]*eff_Demu[i]*effE_ele[i]*effE_ele[i];
      effCov_emu[i][j] += eff_eDamu[i]*eff_eDamu[j]*corr_mu[i][j]*effE_a_mu[i]*effE_a_mu[j];
      effCov_emu[i][j] += eff_eDbmu[i]*eff_eDbmu[j]*corr_mu[i][j]*effE_b_mu[i]*effE_b_mu[j];
      effCov_emu[i][j] += eff_emuDfrac[i]* eff_emuDfrac[j]*fracE_b_mu*fracE_b_mu;	
    }
  }

  double eff_emu_av = 0;
  double effE_emu_av = 0;
  for (int i =0; i<5; ++i){
    std::cout<<i<<" "<<eff_emu[i]<<"\t \\pm "<<sqrt(effCov_emu[i][i])<<" \t in "<<lumFrac[i]<<std::endl;
    eff_emu_av += eff_emu[i]*lumFrac[i];
    for (int j = 0; j< 5; ++j){
      effE_emu_av += effCov_emu[i][j]*lumFrac[i]*lumFrac[j];
    }
  }
  if (effE_emu_av<0)std::cout<<"Weighted err is negative, ooh"<<std::endl;
  effE_emu_av = sqrt(effE_emu_av);
  double eff_emu_mc = 1. - (1.-eff_ele_mc)*((1.-frac_b_mu)*(1.-eff_a_mu_mc)
					    + frac_b_mu*(1.-eff_b_mu_mc));
  std::cout<<"MC "<<eff_emu_mc
	   <<"\t Ave: "<<eff_emu_av<<"\t \\pm "<<effE_emu_av
	   <<"\t SFav: "<< eff_emu_av/eff_emu_mc<<"\t \\pm "<<effE_emu_av/eff_emu_mc
	   <<std::endl;
}


/*

e_ab = e_a*e_b;
ee_ab = 1.- (1.-e_ab)*(1.-e_ab);

AB,AB
ab,ab
  


*/
