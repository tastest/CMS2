void xsecCalc(float obs, float expS, float expSErr, float bg, float bgErr){
  float obsS = obs - bg;
  float obsSErr = sqrt( obs + bgErr*bgErr);

  float sigma = 157.5 * (obs - bg)/ expS;

  float sigmaErrStatFrac = sqrt(obs)/obsS;
  float sigmaErrSystFrac = sqrt(bgErr*bgErr/obsS/obsS + expSErr*expSErr/expS/expS);
  float sigmaErrFrac = sqrt(sigmaErrStatFrac*sigmaErrStatFrac + sigmaErrSystFrac*sigmaErrSystFrac);
  
  float sigmaErrFracLum = 0.11;

  float sigmaErrStat = sigma* sigmaErrStatFrac;
  float sigmaErrSyst = sigma* sigmaErrSystFrac;
  float sigmaErr = sigma*sigmaErrFrac;
  float sigmaErrLum = sigma*sigmaErrFracLum;

  float sigmaErrWLum = sigma * sqrt( sigmaErrFrac*sigmaErrFrac + sigmaErrFracLum*sigmaErrFracLum);

  std::cout<<"sigma = "
	   << sigma << " \\pm "<< sigmaErrStat << " (stat) \\pm " << sigmaErrSyst 
	   << " (syst) \\pm "<<sigmaErrLum << " (lum)"<<std::endl
	   <<"\t = "<<sigma << " \\pm " << sigmaErr << " (syst+stat) \\pm "<< sigmaErrLum << " (lum)"<<std::endl;

}
