double zbi(double n_on, double mu_b_hat, double sigma_b){
 
  //double n_on     = 140.;                         // total events in signal region (S+B)
  //double mu_b_hat = 83.33;                        // mean num of BG events expected in sig. region
  //double sigma_b  = 8.333;                        // uncertainty of mu_b_hat
  
  double tau      = mu_b_hat / (sigma_b*sigma_b); // scale factor to corresp. Noff/Non              
  double n_off    = tau*mu_b_hat;
  double P_Bi     = TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1);
  double Z_Bi     = sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi);           

  cout  <<"  total events in signal region (S+B)               - n_on     " <<n_on      <<endl
        <<"  mean num of BG events expected in sig. region     - mu_b_hat " <<mu_b_hat  <<endl
        <<"  uncertainty of mu_b_hat                           - sigma_b  " <<sigma_b   <<endl
        <<"  scale factor to corresp. Noff/Non                 - tau      " <<tau       <<endl
        <<"  tau*mu_b_hat                                      - n_off    " <<n_off     <<endl
        <<"  TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1) - P_Bi     " <<P_Bi      <<endl
        <<"  sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi)           - Z_Bi     " <<Z_Bi      <<endl;

  return Z_Bi;
}
