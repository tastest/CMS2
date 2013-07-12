#include <fstream>
#include "TH1D.h"
#include "TMatrixD.h"
#include "TMath.h"

double xsection = 154.0*0.06451;

int lineWidth=3;

TString observablename;
TString acceptanceName;
TString xaxislabel;

float xmin=-1.0;
float xmax= 1.0;

const int nbins1D=6;
const int nbins2D=6;

double xbins1D[nbins1D+1];
double xbins2D[nbins2D+1];


Double_t stat_corr  [nbins1D]; //errors include syst error in the unfolding
Double_t stat_uncorr[nbins1D]; //errors do not include syst error in the unfolding
Double_t syst_corr  [nbins1D];

Float_t sign(Float_t t) 
{
  if( t >= 0.0 )
    return 1.0;
  else
    return -1.0;
}


void GetAfb(TH1D* h, Float_t &afb, Float_t  &afberr){

  Int_t nbins = h->GetNbinsX();
  Float_t event_minus;
  Float_t event_plus;
  Float_t event_total;
  Double_t event_plus_err;
  Double_t event_minus_err;

  //event_minus  = h-> IntegralAndError(0, nbins/2, event_plus_err,"");
  event_minus  = h-> IntegralAndError(0, nbins/2, event_minus_err,"");
  //event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_minus_err,"");
  event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_plus_err,"");
  event_total = event_plus + event_minus;

  //cout<<event_minus<<" "<<event_minus_err<<" "<<event_plus<<" "<<event_plus_err<<" "<<event_total<<endl;

  afb = (event_plus-event_minus)/(event_plus+event_minus);
  afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
    + event_minus*event_minus*event_plus_err*event_plus_err)/
    (event_total*event_total*event_total*event_total));

}


//void GetAfbBinByBin(TH1D* h, Float_t &afbbin, Float_t  &afberrbin){
void GetAfbBinByBin(TH1D* h){

  Int_t nbins = h->GetNbinsX();
  const int nbins2 = nbins/2 +1;
  Float_t event_minus[nbins2];
  Float_t event_plus[nbins2];
  Float_t event_total[nbins2];
  Double_t event_plus_err[nbins2];
  Double_t event_minus_err[nbins2];
  Double_t afbbin[nbins2];
  Double_t afberrbin[nbins2];

  Double_t event_plus_total = 0.;
  Double_t event_minus_total = 0.;
  Double_t event_total_total = 0.;

  for(int i=0;i<nbins2;i++){
    //event_minus[i]  = h-> IntegralAndError(i, i, event_plus_err[i],"");
    event_minus[i]  = h-> IntegralAndError(i, i, event_minus_err[i],"");
    event_minus_total += event_minus[i];
    //event_plus[i]   = h-> IntegralAndError(nbins+1-i, nbins+1-i, event_minus_err[i],"");
    event_plus[i]   = h-> IntegralAndError(nbins+1-i, nbins+1-i, event_plus_err[i],"");
    event_plus_total += event_plus[i];
    event_total[i] = event_plus[i] + event_minus[i];
    event_total_total += event_total[i];

    //cout<<event_minus[i]<<" "<<event_minus_err[i]<<" "<<event_plus[i]<<" "<<event_plus_err[i]<<" "<<event_total[i]<<endl;

    afbbin[i] = (event_plus[i]-event_minus[i])/(event_plus[i]+event_minus[i]);
    afberrbin[i]   = sqrt(4*(event_plus[i]*event_plus[i]*event_minus_err[i]*event_minus_err[i] 
      + event_minus[i]*event_minus[i]*event_plus_err[i]*event_plus_err[i])/
      (event_total[i]*event_total[i]*event_total[i]*event_total[i]));
    cout<<i<<" AFB = "<<afbbin[i]<<" +/- "<<afberrbin[i]<<endl;
  }
}


void GetAvsY(TH1* histogram, TMatrixD &covarianceM, vector<double> &myafb, vector<double> &myerr, ofstream& second_output_file){

  myafb.clear();
  myerr.clear();

    //Get info from histogram
  int nbins = histogram->GetNbinsX();
  double n[16];
  for(int i=0;i<nbins;i++){
    n[i] = histogram->GetBinContent(i+1);
  }

    //Output
  double afb[8], err[8];

    //Setup Some Needed Vectors
  double alpha[16], beta[16], dfdn[16];

    //Get Asymmetry for each Y bin
  for(int j=0;j<nbins/2;j++){

    int forBin = nbins/2 + j;
    int bacBin = nbins/2 - j - 1;

    for(int i=0;i<nbins;i++){
      if( i == forBin ){ 
        alpha[i] = 1; 
        beta[i] = 1;
      }else if( i == bacBin ){ 
        alpha[i] = -1; 
        beta[i] = 1;
      }else{
        alpha[i] = 0;
        beta[i] = 0;
      }
    }

    double sum = 0. , diff = 0.;
    for(int i=0;i<nbins;i++){
      sum += beta[i] * n[i];
      diff += alpha[i] * n[i];
    }

      //Calculate Everything
    if(sum > 0){ 

  //Error Calculation
      for(int i=0;i<nbins;i++){
        dfdn[i] = ( alpha[i] * sum - beta[i] * diff ) / pow(sum,2);
      }

      double afberr = 0.;
      for(int i=0;i<nbins;i++){
        for(int k=0;k<nbins;k++){
          afberr += covarianceM(i,k) * dfdn[i] * dfdn[k];
      //if(i==k) cout<<"DAH: "<<n[i]<<" "<<k<<" "<<covarianceM(i,k)<<endl;
        }
      }
      afberr = sqrt(afberr);

      err[j] = afberr;
      afb[j] = diff / sum; 

    }else{ 

      afb[j] = 0.; 
      err[j] = 0.;

    }
    myafb.push_back(afb[j]);
    myerr.push_back(err[j]);

    cout<<j<<" AFB = "<<afb[j]<<" +/- "<<err[j]<<endl;
    second_output_file << acceptanceName << " " << observablename << " AFB" << j << ": " << afb[j] << " +/- " << err[j] << endl;    	
  }
}

void GetCorrectedAfb(TH1D* histogram, TMatrixD &covarianceM, Float_t &afb, Float_t  &afberr){

    //Need to calculate AFB and Error for the fully corrected distribution, m_correctE(j,i)

    //Get histogram info
  int nbins = histogram->GetNbinsX();
  double n[16];
  for(int i=0;i<nbins;i++){
    n[i] = histogram->GetBinContent(i+1);
  }

    //Setup Alpha Vector
  double alpha[16], beta[16];
  for(int i=0;i<nbins;i++) if(i < nbins/2 ){ alpha[i] = -1;}else{ alpha[i] = 1;}

    //Components of the error calculation
  double sum_n = 0.;
  double sum_alpha_n = 0.;
  for(int i=0;i<nbins;i++){
    sum_n += n[i];
    sum_alpha_n += alpha[i] * n[i];
  }

  double dfdn[16];
  for(int i=0;i<nbins;i++){
    dfdn[i] = ( alpha[i] * sum_n - sum_alpha_n ) / pow(sum_n,2);
  }

    //Error Calculation
  afberr = 0.;
  for(int i=0;i<nbins;i++){
    for(int j=0;j<nbins;j++){
      afberr += covarianceM(i,j) * dfdn[i] * dfdn[j];
    }
  }
  afberr = sqrt(afberr);

    //Calculate Afb
  afb = sum_alpha_n / sum_n;

    //    cout<<"AFB = "<<afb<<" "<<afberr<<endl;
}




void GetCorrectedAfbBinByBin(TH1D* histogram, TMatrixD &covarianceM, vector<double> &myafb, vector<double> &myerr, ofstream& second_output_file){

    //Need to calculate AFB and Error for the fully corrected distribution, m_correctE(j,i)

  myafb.clear();
  myerr.clear();

    //Get histogram info
  int nbins = histogram->GetNbinsX();
  const int nbins2 = nbins/2;

  Double_t afbbin[nbins2];
  Double_t afberrbin[nbins2];

  double n[16];
  for(int i=0;i<nbins;i++){
    n[i] = histogram->GetBinContent(i+1);
  }

    //Setup Alpha Vector
  double alpha[16], beta[16];
  for(int i=0;i<nbins;i++) if(i < nbins/2 ){ alpha[i] = -1;}else{ alpha[i] = 1;}

    //Components of the error calculation
  double sum_n[nbins2];
  double sum_alpha_n[nbins2];
  double sum_n_total = 0.;
  double sum_alpha_n_total = 0.;


  for(int i=0;i<nbins2;i++){
    sum_n[i] = n[i] + n[nbins-1-i];
    sum_alpha_n[i] = alpha[i] * n[i] + alpha[nbins-1-i] * n[nbins-1-i];
    sum_n_total += sum_n[i];
    sum_alpha_n_total += sum_alpha_n[i];
  }

  double dfdn[16];
  for(int i=0;i<nbins;i++){
    int k = -999;
    if (i < nbins2) k = i;
    else k = nbins-1-i;
    dfdn[i] = ( alpha[i] * sum_n[k] - sum_alpha_n[k] ) / pow(sum_n[k],2);
  }

    //Error Calculation

  for(int k=0;k<nbins2;k++){
    afberrbin[k] = 0.;
    for(int i=0;i<nbins;i++){
      for(int j=0;j<nbins;j++){
        if( (i==k || i==nbins-1-k ) && (j==k || j==nbins-1-k ) ) {
          afberrbin[k] += covarianceM(i,j) * dfdn[i] * dfdn[j];
              //cout<<covarianceM(i,j)<<" "<<dfdn[i]<<" "<<dfdn[j]<<" "<<endl;
        }
      }
    }
    afberrbin[k] = sqrt(afberrbin[k]);
    afbbin[k] = sum_alpha_n[k] / sum_n[k];
    cout<<k<<" AFB = "<<afbbin[k]<<" +/- "<<afberrbin[k]<<endl;
    second_output_file << acceptanceName << " " << observablename << " AFB" << i << ": " << afbbin[k] << " +/- " << afberrbin[k] << endl;

    myafb.push_back(afbbin[k]);
    myerr.push_back(afberrbin[k]);
  }
  double afb = sum_alpha_n_total / sum_n_total;
    //cout<<"AFB = "<<afb<<endl;
}







void Initialize1DBinning(int iVar){


  switch (iVar)
  {
    //   Lepton Charge Asymmetry
    case 0:
    {
      observablename="lep_charge_asymmetry";
      xaxislabel="|#eta_{l+}|-|#eta_{l-}|";
      acceptanceName="lepChargeAsym";
      xbins1D[0]=-2.0; xbins1D[1]=-0.8; xbins1D[2]=-0.4; xbins1D[3]=0.0; xbins1D[4]=0.4; xbins1D[5]=0.8; xbins1D[6]=2.0;
syst_corr[0] =  0.001794  ; stat_corr[0] =  0.005346  ; stat_uncorr[0] =  0.003457  ;
syst_corr[1] =  0.001924  ; stat_corr[1] =  0.006949  ; stat_uncorr[1] =  0.004681  ;
syst_corr[2] =  0.003580  ; stat_corr[2] =  0.007129  ; stat_uncorr[2] =  0.005290  ;
syst_corr[3] =  0.003847  ; stat_corr[3] =  0.006887  ; stat_uncorr[3] =  0.005316  ;
syst_corr[4] =  0.001992  ; stat_corr[4] =  0.006162  ; stat_uncorr[4] =  0.004696  ;
syst_corr[5] =  0.002427  ; stat_corr[5] =  0.004645  ; stat_uncorr[5] =  0.003491  ;
      xmin=-2.0;
      xmax= 2.0;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="#Delta#phi_{l+l-}";
      acceptanceName="lepAzimAsym2";
      Double_t pi = 3.141592653589793;
      xbins1D[0]=0.0; xbins1D[1]=4.*pi/20.; xbins1D[2]=7.*pi/20.; xbins1D[3]=10.*pi/20.; xbins1D[4]=13.*pi/20.; xbins1D[5]=16.*pi/20.; xbins1D[6]=pi;
syst_corr[0] =  0.005598  ; stat_corr[0] =  0.007941  ; stat_uncorr[0] =  0.005696  ;
syst_corr[1] =  0.004784  ; stat_corr[1] =  0.006734  ; stat_uncorr[1] =  0.004737  ;
syst_corr[2] =  0.002737  ; stat_corr[2] =  0.006163  ; stat_uncorr[2] =  0.004353  ;
syst_corr[3] =  0.002514  ; stat_corr[3] =  0.006491  ; stat_uncorr[3] =  0.004404  ;
syst_corr[4] =  0.005171  ; stat_corr[4] =  0.006872  ; stat_uncorr[4] =  0.004837  ;
syst_corr[5] =  0.006179  ; stat_corr[5] =  0.008617  ; stat_uncorr[5] =  0.006307  ;
      xmin=0.0;
      xmax=pi;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lep_costheta_cms";
      xaxislabel="cos(#theta^{+}_{l})";
      acceptanceName="lepPlusCosTheta";
      xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
syst_corr[0] =  0.026433  ; stat_corr[0] =  0.022191  ; stat_uncorr[0] =  0.017472  ;
syst_corr[1] =  0.021013  ; stat_corr[1] =  0.016553  ; stat_uncorr[1] =  0.011947  ;
syst_corr[2] =  0.013153  ; stat_corr[2] =  0.017841  ; stat_uncorr[2] =  0.012534  ;
syst_corr[3] =  0.010084  ; stat_corr[3] =  0.015056  ; stat_uncorr[3] =  0.011457  ;
syst_corr[4] =  0.019675  ; stat_corr[4] =  0.015144  ; stat_uncorr[4] =  0.010476  ;
syst_corr[5] =  0.028864  ; stat_corr[5] =  0.025112  ; stat_uncorr[5] =  0.016376  ;
      xmin=-1.0;
      xmax= 1.0;
      break;
    }
    //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos(#theta^{-}_{l})";
      acceptanceName="lepMinusCosTheta";
      xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
syst_corr[0] =  0.023301  ; stat_corr[0] =  0.024151  ; stat_uncorr[0] =  0.017495  ;
syst_corr[1] =  0.009467  ; stat_corr[1] =  0.015927  ; stat_uncorr[1] =  0.011962  ;
syst_corr[2] =  0.006665  ; stat_corr[2] =  0.017457  ; stat_uncorr[2] =  0.012514  ;
syst_corr[3] =  0.009335  ; stat_corr[3] =  0.015462  ; stat_uncorr[3] =  0.011421  ;
syst_corr[4] =  0.013100  ; stat_corr[4] =  0.014852  ; stat_uncorr[4] =  0.010626  ;
syst_corr[5] =  0.016321  ; stat_corr[5] =  0.023841  ; stat_uncorr[5] =  0.016631  ;
      xmin=-1.0;
      xmax= 1.0;
      break;
    }
    //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lep_costheta_cms";
      xaxislabel="cos(#theta^{*}_{l})";
      acceptanceName="lepCosTheta";
      xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
syst_corr[0] =  0.024422  ; stat_corr[0] =  0.019040  ; stat_uncorr[0] =  0.012361  ;
syst_corr[1] =  0.014372  ; stat_corr[1] =  0.012204  ; stat_uncorr[1] =  0.008452  ;
syst_corr[2] =  0.008457  ; stat_corr[2] =  0.011424  ; stat_uncorr[2] =  0.008855  ;
syst_corr[3] =  0.008922  ; stat_corr[3] =  0.011026  ; stat_uncorr[3] =  0.008088  ;
syst_corr[4] =  0.016064  ; stat_corr[4] =  0.010472  ; stat_uncorr[4] =  0.007461  ;
syst_corr[5] =  0.021688  ; stat_corr[5] =  0.015978  ; stat_uncorr[5] =  0.011670  ; 
      xmin=-1.0;
      xmax= 1.0;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="cos(#theta_{l+,n})cos(#theta_{l-,n})";
      acceptanceName="topSpinCorr";
      xbins1D[0]=-1.0; xbins1D[1]=-0.5; xbins1D[2]=-0.2; xbins1D[3]=0.0; xbins1D[4]=0.2; xbins1D[5]=0.5; xbins1D[6]=1.0;
syst_corr[0] =  0.009005  ; stat_corr[0] =  0.015109  ; stat_uncorr[0] =  0.009931  ;
syst_corr[1] =  0.016934  ; stat_corr[1] =  0.030011  ; stat_uncorr[1] =  0.019719  ;
syst_corr[2] =  0.021668  ; stat_corr[2] =  0.040640  ; stat_uncorr[2] =  0.030935  ;
syst_corr[3] =  0.023122  ; stat_corr[3] =  0.040636  ; stat_uncorr[3] =  0.027907  ;
syst_corr[4] =  0.017626  ; stat_corr[4] =  0.026800  ; stat_uncorr[4] =  0.019072  ;
syst_corr[5] =  0.005839  ; stat_corr[5] =  0.010713  ; stat_uncorr[5] =  0.007566  ;
      xmin=-1.0;
      xmax= 1.0;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapidtiydiff_Marco";
      xaxislabel="|y_{top}|-|y_{tbar}|";
      acceptanceName="rapiditydiffMarco";
      xbins1D[0]=-2.0; xbins1D[1]=-0.7; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.7; xbins1D[6]=2.0;
syst_corr[0] =  0.003890  ; stat_corr[0] =  0.006153  ; stat_uncorr[0] =  0.004584  ;
syst_corr[1] =  0.003297  ; stat_corr[1] =  0.011308  ; stat_uncorr[1] =  0.008458  ;
syst_corr[2] =  0.008537  ; stat_corr[2] =  0.014555  ; stat_uncorr[2] =  0.010529  ;
syst_corr[3] =  0.012801  ; stat_corr[3] =  0.014194  ; stat_uncorr[3] =  0.010549  ;
syst_corr[4] =  0.003135  ; stat_corr[4] =  0.010855  ; stat_uncorr[4] =  0.008355  ;
syst_corr[5] =  0.002232  ; stat_corr[5] =  0.006124  ; stat_uncorr[5] =  0.004587  ;
      xmin=-2.0;
      xmax= 2.0;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="cos(#theta_{top})";
      acceptanceName="topCosTheta";
      xbins1D[0]=-1.0; xbins1D[1]=-0.7; xbins1D[2]=-0.4; xbins1D[3]=0.0; xbins1D[4]=0.4; xbins1D[5]=0.7; xbins1D[6]=1.0;
syst_corr[0] =  0.016200  ; stat_corr[0] =  0.031149  ; stat_uncorr[0] =  0.023930  ;
syst_corr[1] =  0.009050  ; stat_corr[1] =  0.014825  ; stat_uncorr[1] =  0.010633  ;
syst_corr[2] =  0.006642  ; stat_corr[2] =  0.011610  ; stat_uncorr[2] =  0.008159  ;
syst_corr[3] =  0.005216  ; stat_corr[3] =  0.010612  ; stat_uncorr[3] =  0.008202  ;
syst_corr[4] =  0.007476  ; stat_corr[4] =  0.013658  ; stat_uncorr[4] =  0.010493  ;
syst_corr[5] =  0.014691  ; stat_corr[5] =  0.032764  ; stat_uncorr[5] =  0.023655  ;
      xmin=-1.0;
      xmax= 1.0;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 8:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="cos(#Delta#phi_{l+l-})";
      acceptanceName="lepAzimAsym";
      xbins1D[0]=-1.0; xbins1D[1]=-0.8; xbins1D[2]=-0.4; xbins1D[3]=0.0; xbins1D[4]=0.4; xbins1D[5]=0.8; xbins1D[6]=1.0;
      stat_corr[0] = 0.01; stat_corr [1] = 0.01;  stat_corr [2] = 0.01;  stat_corr [3] = 0.01; stat_corr [4] = 0.01;  stat_corr [5] = 0.01;  stat_corr [6] = 0.01; 
      stat_uncorr[0] = 0.00; stat_uncorr[1] = 0.00; stat_uncorr[2] = 0.00; stat_uncorr[3] = 0.00; stat_uncorr[4] = 0.00; stat_uncorr[5] = 0.00; stat_uncorr[6] = 0.00;
      syst_corr[0] = 0.02; syst_corr [1] = 0.02;  syst_corr [2] = 0.02;  syst_corr [3] = 0.02; syst_corr [4] = 0.02;  syst_corr [5] = 0.02;  syst_corr [6] = 0.02; 
      xmin=-1.0;
      xmax= 1.0;
      break;
    }
    //   Top Asy I
    case 9:
    {
      observablename="top_pseudorapidtiydiff_cms";
      xaxislabel="|#eta_{top}|-|#eta_{tbar}|";
      acceptanceName="pseudorapiditydiff";
      xbins1D[0]=-4.0; xbins1D[1]=-1.0; xbins1D[2]=-0.5; xbins1D[3]=0.0; xbins1D[4]=0.5; xbins1D[5]=1.0; xbins1D[6]=4.0;
      stat_corr[0] = 0.01; stat_corr [1] = 0.01;  stat_corr [2] = 0.01;  stat_corr [3] = 0.01; stat_corr [4] = 0.01;  stat_corr [5] = 0.01;  stat_corr [6] = 0.01; 
      stat_uncorr[0] = 0.00; stat_uncorr[1] = 0.00; stat_uncorr[2] = 0.00; stat_uncorr[3] = 0.00; stat_uncorr[4] = 0.00; stat_uncorr[5] = 0.00; stat_uncorr[6] = 0.00;
      syst_corr[0] = 0.02; syst_corr [1] = 0.02;  syst_corr [2] = 0.02;  syst_corr [3] = 0.02; syst_corr [4] = 0.02;  syst_corr [5] = 0.02;  syst_corr [6] = 0.02; 
      xmin=-4.0;
      xmax= 4.0;
      break;
    }
    //   Top Asy II
    case 10:
    {
      observablename="top_rapidtiydiff_cms";
      xaxislabel="(y_{top}-y_{tbar})(y_{top}+y_{tbar})";
      acceptanceName="rapiditydiff";
      xbins1D[0]=-4.0; xbins1D[1]=-0.8; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.8; xbins1D[6]=4.0;
      stat_corr[0] = 0.01; stat_corr [1] = 0.01;  stat_corr [2] = 0.01;  stat_corr [3] = 0.01; stat_corr [4] = 0.01;  stat_corr [5] = 0.01;  stat_corr [6] = 0.01; 
      stat_uncorr[0] = 0.00; stat_uncorr[1] = 0.00; stat_uncorr[2] = 0.00; stat_uncorr[3] = 0.00; stat_uncorr[4] = 0.00; stat_uncorr[5] = 0.00; stat_uncorr[6] = 0.00;
      syst_corr[0] = 0.02; syst_corr [1] = 0.02;  syst_corr [2] = 0.02;  syst_corr [3] = 0.02; syst_corr [4] = 0.02;  syst_corr [5] = 0.02;  syst_corr [6] = 0.02; 
      xmin=-4.0;
      xmax= 4.0;
      break;
    }
    default:
    {
      cout<<"Set the variable switch";
    }
  }
}

void Initialize2DBinning(int iVar){


  switch (iVar)
  {
    //   Lepton Charge Asymmetry
    case 0:
    {
      observablename="lep_charge_asymmetry";
      xaxislabel="M_{t#bar t}";
      acceptanceName="lepChargeAsym";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.006227  ; stat_corr[0] =  0.007378  ; stat_uncorr[0] =  0.005607  ;
syst_corr[1] =  0.006933  ; stat_corr[1] =  0.016666  ; stat_uncorr[1] =  0.011877  ;
syst_corr[2] =  0.006046  ; stat_corr[2] =  0.023699  ; stat_uncorr[2] =  0.016971  ;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="M_{t#bar t}";
      acceptanceName="lepAzimAsym2";
      Double_t pi = 3.141592653589793;
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.023448  ; stat_corr[0] =  0.008685  ; stat_uncorr[0] =  0.006075  ;
syst_corr[1] =  0.024505  ; stat_corr[1] =  0.018431  ; stat_uncorr[1] =  0.012799  ;
syst_corr[2] =  0.023928  ; stat_corr[2] =  0.023935  ; stat_uncorr[2] =  0.016685  ;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lep_costheta_cms";
      xaxislabel="M_{t#bar t}";
      acceptanceName="lepPlusCosTheta";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.014551  ; stat_corr[0] =  0.009797  ; stat_uncorr[0] =  0.007038  ;
syst_corr[1] =  0.048631  ; stat_corr[1] =  0.026873  ; stat_uncorr[1] =  0.018776  ;
syst_corr[2] =  0.071372  ; stat_corr[2] =  0.038785  ; stat_uncorr[2] =  0.026632  ;
      break;
    }
      //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="M_{t#bar t}";
      acceptanceName="lepMinusCosTheta";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.011335  ; stat_corr[0] =  0.009746  ; stat_uncorr[0] =  0.007026  ;
syst_corr[1] =  0.033818  ; stat_corr[1] =  0.026464  ; stat_uncorr[1] =  0.018820  ;
syst_corr[2] =  0.055506  ; stat_corr[2] =  0.038864  ; stat_uncorr[2] =  0.026844  ;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lep_costheta_cms";
      xaxislabel="M_{t#bar t}";
      acceptanceName="lepCosTheta";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.012416  ; stat_corr[0] =  0.006649  ; stat_uncorr[0] =  0.004972  ;
syst_corr[1] =  0.040286  ; stat_corr[1] =  0.018013  ; stat_uncorr[1] =  0.013292  ;
syst_corr[2] =  0.062468  ; stat_corr[2] =  0.026273  ; stat_uncorr[2] =  0.018905  ;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="M_{t#bar t}";
      acceptanceName="topSpinCorr";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.015438  ; stat_corr[0] =  0.011742  ; stat_uncorr[0] =  0.008875  ;
syst_corr[1] =  0.031088  ; stat_corr[1] =  0.033467  ; stat_uncorr[1] =  0.025140  ;
syst_corr[2] =  0.036884  ; stat_corr[2] =  0.048251  ; stat_uncorr[2] =  0.036018  ;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapidtiydiff_Marco";
      xaxislabel="M_{t#bar t}";
      acceptanceName="rapiditydiffMarco";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.002885  ; stat_corr[0] =  0.008910  ; stat_uncorr[0] =  0.006635  ;
syst_corr[1] =  0.007974  ; stat_corr[1] =  0.024875  ; stat_uncorr[1] =  0.018155  ;
syst_corr[2] =  0.009536  ; stat_corr[2] =  0.036063  ; stat_uncorr[2] =  0.025779  ;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="M_{t#bar t}";
      acceptanceName="topCosTheta";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.003006  ; stat_corr[0] =  0.009327  ; stat_uncorr[0] =  0.006702  ;
syst_corr[1] =  0.006987  ; stat_corr[1] =  0.025379  ; stat_uncorr[1] =  0.018333  ;
syst_corr[2] =  0.010131  ; stat_corr[2] =  0.035790  ; stat_uncorr[2] =  0.026031  ;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 8:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="M_{t#bar t}";
      acceptanceName="lepAzimAsym";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    //   Top Asy I
    case 9:
    {
      observablename="top_pseudorapidtiydiff_cms";
      xaxislabel="M_{t#bar t}";
      acceptanceName="pseudorapiditydiff";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    //   Top Asy II
    case 10:
    {
      observablename="top_rapidtiydiff_cms";
      xaxislabel="M_{t#bar t}";
      acceptanceName="rapiditydiff";
      xbins2D[0]=-800.0; xbins2D[1]=-510.0; xbins2D[2]=-410.0; xbins2D[3]=0.0; xbins2D[4]=410; xbins2D[5]=510.0; xbins2D[6]=800.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    default:
    {
      cout<<"Set the variable switch";
    }
  }
}

void Initialize2DBinningttpt(int iVar){


  switch (iVar)
  {
    //   Lepton Charge Asymmetry
    case 0:
    {
      observablename="lep_charge_asymmetry";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepChargeAsym";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.004836  ; stat_corr[0] =  0.008115  ; stat_uncorr[0] =  0.005153  ;
syst_corr[1] =  0.005043  ; stat_corr[1] =  0.017765  ; stat_uncorr[1] =  0.011602  ;
syst_corr[2] =  0.011328  ; stat_corr[2] =  0.024255  ; stat_uncorr[2] =  0.016729  ;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepAzimAsym2";
      Double_t pi = 3.141592653589793;
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.008098  ; stat_corr[0] =  0.008007  ; stat_uncorr[0] =  0.005079  ;
syst_corr[1] =  0.014371  ; stat_corr[1] =  0.018177  ; stat_uncorr[1] =  0.011628  ;
syst_corr[2] =  0.029878  ; stat_corr[2] =  0.026496  ; stat_uncorr[2] =  0.016883  ;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lep_costheta_cms";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepPlusCosTheta";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.016303  ; stat_corr[0] =  0.009889  ; stat_uncorr[0] =  0.007379  ;
syst_corr[1] =  0.037422  ; stat_corr[1] =  0.028339  ; stat_uncorr[1] =  0.020271  ;
syst_corr[2] =  0.050842  ; stat_corr[2] =  0.042338  ; stat_uncorr[2] =  0.029472  ;
      break;
    }
      //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepMinusCosTheta";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.013284  ; stat_corr[0] =  0.009941  ; stat_uncorr[0] =  0.007402  ;
syst_corr[1] =  0.024501  ; stat_corr[1] =  0.027525  ; stat_uncorr[1] =  0.020346  ;
syst_corr[2] =  0.033907  ; stat_corr[2] =  0.040658  ; stat_uncorr[2] =  0.029570  ;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lep_costheta_cms";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepCosTheta";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.013950  ; stat_corr[0] =  0.007652  ; stat_uncorr[0] =  0.005226  ;
syst_corr[1] =  0.027920  ; stat_corr[1] =  0.020792  ; stat_uncorr[1] =  0.014360  ;
syst_corr[2] =  0.037950  ; stat_corr[2] =  0.029614  ; stat_uncorr[2] =  0.020875  ;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="topSpinCorr";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.007068  ; stat_corr[0] =  0.014116  ; stat_uncorr[0] =  0.009345  ;
syst_corr[1] =  0.019308  ; stat_corr[1] =  0.039766  ; stat_uncorr[1] =  0.026297  ;
syst_corr[2] =  0.023303  ; stat_corr[2] =  0.058066  ; stat_uncorr[2] =  0.038355  ;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapidtiydiff_Marco";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="rapiditydiffMarco";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.002966  ; stat_corr[0] =  0.010176  ; stat_uncorr[0] =  0.007568  ;
syst_corr[1] =  0.010711  ; stat_corr[1] =  0.028293  ; stat_uncorr[1] =  0.020901  ;
syst_corr[2] =  0.010683  ; stat_corr[2] =  0.041044  ; stat_uncorr[2] =  0.030389  ;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="topCosTheta";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.004295  ; stat_corr[0] =  0.010517  ; stat_uncorr[0] =  0.007593  ;
syst_corr[1] =  0.012710  ; stat_corr[1] =  0.029721  ; stat_uncorr[1] =  0.020969  ;
syst_corr[2] =  0.015388  ; stat_corr[2] =  0.043829  ; stat_uncorr[2] =  0.030487  ;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 8:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepAzimAsym";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    //   Top Asy I
    case 9:
    {
      observablename="top_pseudorapidtiydiff_cms";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="pseudorapiditydiff";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    //   Top Asy II
    case 10:
    {
      observablename="top_rapidtiydiff_cms";
      xaxislabel="p_{T,t#bar{t}}";
      acceptanceName="rapiditydiff";
      xbins2D[0]=-100.0; xbins2D[1]=-52.0; xbins2D[2]=-24.0; xbins2D[3]=0.0; xbins2D[4]=24; xbins2D[5]=52.0; xbins2D[6]=100.0;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    default:
    {
      cout<<"Set the variable switch";
    }
  }
}

void Initialize2DBinningttrapidity2(int iVar){


  switch (iVar)
  {
    //   Lepton Charge Asymmetry
    case 0:
    {
      observablename="lep_charge_asymmetry";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="lepChargeAsym";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.003099  ; stat_corr[0] =  0.007396  ; stat_uncorr[0] =  0.005278  ;
syst_corr[1] =  0.005609  ; stat_corr[1] =  0.016024  ; stat_uncorr[1] =  0.011594  ;
syst_corr[2] =  0.012690  ; stat_corr[2] =  0.022065  ; stat_uncorr[2] =  0.016291  ;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="lepAzimAsym2";
      Double_t pi = 3.141592653589793;
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.009396  ; stat_corr[0] =  0.006898  ; stat_uncorr[0] =  0.005202  ;
syst_corr[1] =  0.013361  ; stat_corr[1] =  0.015107  ; stat_uncorr[1] =  0.011635  ;
syst_corr[2] =  0.025226  ; stat_corr[2] =  0.021287  ; stat_uncorr[2] =  0.016371  ;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lep_costheta_cms";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="lepPlusCosTheta";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.015559  ; stat_corr[0] =  0.011571  ; stat_uncorr[0] =  0.007450  ;
syst_corr[1] =  0.035438  ; stat_corr[1] =  0.030569  ; stat_uncorr[1] =  0.019981  ;
syst_corr[2] =  0.049866  ; stat_corr[2] =  0.042447  ; stat_uncorr[2] =  0.028264  ;
      break;
    }
      //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="lepMinusCosTheta";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.012326  ; stat_corr[0] =  0.011079  ; stat_uncorr[0] =  0.007475  ;
syst_corr[1] =  0.026414  ; stat_corr[1] =  0.029408  ; stat_uncorr[1] =  0.020071  ;
syst_corr[2] =  0.040896  ; stat_corr[2] =  0.041478  ; stat_uncorr[2] =  0.028425  ;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lep_costheta_cms";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="lepCosTheta";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.012946  ; stat_corr[0] =  0.007235  ; stat_uncorr[0] =  0.005277  ;
syst_corr[1] =  0.030079  ; stat_corr[1] =  0.019050  ; stat_uncorr[1] =  0.014160  ;
syst_corr[2] =  0.044795  ; stat_corr[2] =  0.026322  ; stat_uncorr[2] =  0.020041  ;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="topSpinCorr";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.011796  ; stat_corr[0] =  0.012633  ; stat_uncorr[0] =  0.009257  ;
syst_corr[1] =  0.029986  ; stat_corr[1] =  0.035685  ; stat_uncorr[1] =  0.025730  ;
syst_corr[2] =  0.036571  ; stat_corr[2] =  0.051251  ; stat_uncorr[2] =  0.036619  ;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapidtiydiff_Marco";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="rapiditydiffMarco";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.007830  ; stat_corr[0] =  0.008845  ; stat_uncorr[0] =  0.006797  ;
syst_corr[1] =  0.009733  ; stat_corr[1] =  0.024428  ; stat_uncorr[1] =  0.018642  ;
syst_corr[2] =  0.012993  ; stat_corr[2] =  0.035239  ; stat_uncorr[2] =  0.026382  ;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="topCosTheta";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] =  0.007489  ; stat_corr[0] =  0.009511  ; stat_uncorr[0] =  0.006897  ;
syst_corr[1] =  0.010966  ; stat_corr[1] =  0.026434  ; stat_uncorr[1] =  0.018885  ;
syst_corr[2] =  0.013571  ; stat_corr[2] =  0.038226  ; stat_uncorr[2] =  0.026731  ;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 8:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="lepAzimAsym";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    //   Top Asy I
    case 9:
    {
      observablename="top_pseudorapidtiydiff_cms";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="pseudorapiditydiff";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    //   Top Asy II
    case 10:
    {
      observablename="top_rapidtiydiff_cms";
      xaxislabel="y_{t#bar{t}}";
      acceptanceName="rapiditydiff";
      xbins2D[0]=-1.5; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=1.5;
      xmin=xbins2D[0];
      xmax=xbins2D[6];
      break;
    }
    default:
    {
      cout<<"Set the variable switch";
    }
  }
}

void fillUnderOverFlow(TH1D *h1, float value, double weight, int Nsolns)
{
  double min = h1->GetXaxis()->GetXmin();
  double max = h1->GetXaxis()->GetXmax();

  if (value >= max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value <= min) value = h1->GetBinCenter(1);

  int bin_number = h1->FindBin(value);
  double orig_content = h1->GetBinContent(bin_number);
  double orig_error = h1->GetBinError(bin_number);

  //h1->Fill(value, weight);
  h1->SetBinContent( bin_number, orig_content+weight );
  h1->SetBinError( bin_number, sqrt( orig_error*orig_error + weight*weight*double(Nsolns) ) );
}

//--------------------------------------------------------------------

void fillUnderOverFlow(TH2D *h2, float xvalue, float yvalue, double weight, int Nsolns)
{
  double maxx = h2->GetXaxis()->GetXmax();
  double minx = h2->GetXaxis()->GetXmin();
  double maxy = h2->GetYaxis()->GetXmax();
  double miny = h2->GetYaxis()->GetXmin();

  if (xvalue >= maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
  if (xvalue <= minx) xvalue = h2->GetXaxis()->GetBinCenter(1);
  if (yvalue >= maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());
  if (yvalue <= miny) yvalue = h2->GetYaxis()->GetBinCenter(1);

  int bin_number = h2->FindBin(xvalue,yvalue);
  double orig_content = h2->GetBinContent(bin_number);
  double orig_error = h2->GetBinError(bin_number);

  //h2->Fill(xvalue, yvalue, weight);
  h2->SetBinContent( bin_number, orig_content+weight );
  h2->SetBinError( bin_number, sqrt( orig_error*orig_error + weight*weight*double(Nsolns) ) );
}
