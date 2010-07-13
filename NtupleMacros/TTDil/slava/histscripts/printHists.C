#include "TMath.h"
#include <iostream>
#include <vector>
#include "TROOT.h"
#include "TH1F.h"
#include "Math/Math.h"
#include "Math/SpecFunc.h"
#include "TFile.h"
#include <cmath>
#include "histtools.h"

#include "specFun.C"


//template <typename T>
//T max(const T& a, const T& b){ return a> b ? a : b; }
//template <typename T>
//T min(const T& a, const T& b){ return a< b ? a : b; }

void getSoverRootN(double& rat, double& ratE, double s, double n, double sE, double nE){
  rat = s; rat /= n > 0 ? sqrt(n) : 1.;
  ratE = ( sE*sE* ( 1. - 0.5*rat)*( 1. - 0.5*rat)  + (nE*nE - sE*sE)*0.25*rat*rat);
  ratE /= n > 0 ? n : 1.;
  ratE = sqrt(ratE);
}

double bOnlyProb(double s, double b, double bE){
  //at some point need to put some protections here or find an appropriate code
  unsigned int lowExp = floor(s+b);
  double pSum = 0;
  for (int i=lowExp; i>=0; --i){
    if (b> 0. && bE/b>0.03 && s>0.1*b && bE<0.5*s){
      pSum+= poisson_smeared_prob(i,b, bE);
    } else if (b>0.){//use regualar Poisson here
      pSum+= TMath::Poisson(i,b);
    }
  }
  return pSum;
}

std::string formatFloat(double x, const char* formatS){
  std::string xS = Form(Form("%s", formatS),x);
  double xB = atof(xS.c_str());
  if (x>0 && xB==0){
    xS = Form(" %6.1g",x);
  }
  return xS;
}

std::string formatFloat(double x, double xRef, int nSignificant, int printStyle, bool latex = false){
  std::string fS;
  bool isInt = false;
  int nDigs = 0;
  if (fabs(xRef) <= x && x>1E-12){
    nDigs =   std::floor(log(fabs(xRef))/log(10.));
    if (fabs(xRef)>1){
      nDigs +=1;
      if (nDigs >= nSignificant){
	fS = std::string("%d");//+std::string(Form("%dd", nDigs));
	isInt = true;
      } else {
	fS = std::string("%")+std::string(Form("%d.%df", nSignificant, nSignificant-nDigs));
      }
    } else {
      nDigs = abs(nDigs);
      fS = std::string("%")+std::string(Form("%d.%df", nDigs+nSignificant, nSignificant+nDigs-1));
    }
  }
  if (x<=1E-12){
    if (nSignificant>1) fS = std::string("%")+std::string(Form("%d.%df", nSignificant, nSignificant-1));
    if (nSignificant<=1) fS = std::string("%d");
  }

  std::string pmSign  = latex ? " \\pm " : " &plusmn; "; 
  std::string mathSep = latex ? "$" : "";
  std::string timesSign = latex ? "\\times" : "X";

  
  std::string xS;
  if (printStyle == 0){
    if (isInt){
      xS = mathSep + std::string(Form(Form("%s", fS.c_str()),(long int)x)) + mathSep;
    } else {
      xS = mathSep + std::string(Form(Form("%s", fS.c_str()),x)) + mathSep;
    }
  }
  if (printStyle == 1){
    if (isInt){
      xS = mathSep + std::string(Form(Form("%s", fS.c_str()),(long int)x)) 
	+ pmSign + std::string(Form(Form("%s", fS.c_str()),(long int)xRef)) + mathSep;
    } else {
      if (fabs(xRef)> 1E-12 && fabs(xRef) < 0.001){
	double xTmp = pow(10,nDigs)*x;
	double xRefTmp = pow(10,nDigs)*xRef;
	fS = std::string("%")+std::string(Form("%d.%df", nSignificant, nSignificant-1));
	xS = mathSep + std::string(Form(Form("(%s", fS.c_str()),xTmp)) 
	  + pmSign + std::string(Form(Form("%s) \\times 10^{-%d}", fS.c_str(), nDigs),xRefTmp)) + mathSep;
      } else {
	xS = mathSep + std::string(Form(Form("%s", fS.c_str()),x))
	  + pmSign + std::string(Form(Form("%s", fS.c_str()),xRef)) + mathSep;
      }
    }
  }
  return xS;
}

void printNJets( bool latex=false, const char* formatS = "%6.1f", const char* signalS= "ttdil"){
  char* suffix[4];
  suffix[0] = "ee";
  suffix[1] = "mm";
  suffix[2] = "em";
  suffix[3] = "all";


  std::string pmSign  = latex ? " \\pm " : " &plusmn; ";
  std::string colSep  = latex ? " & " : " | ";
  std::string beginL  = latex ? ""   : " | ";
  std::string endL    = latex ? " \\\\ " : " | ";
  std::string mathSep = latex ? "$" : "";
  

  if (latex) {
    std::cout << "\\begin{table}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "{\\small" << std::endl;
  }

     
  for (int sample=0; sample<4; sample++) {
    
    hist::stack(Form("st_hnJet_%s", suffix[sample]), Form("^[^_]+_hnJet_%s$", suffix[sample]));
    THStack* thisStack = (THStack*) 
      gROOT->FindObjectAny(Form("st_hnJet_%s", suffix[sample]));
    //    std::cout<<"Found "<<thisStack->GetName()<<std::endl;

    int nHists = thisStack->GetHists()->GetSize();
    if (latex) {
      std::cout << "\\begin{minipage}{0.48\\textwidth}" << std::endl;
      if (sample == 0 || sample == 2){
	std::cout << "\\begin{tabular}{|l|c|c|c|} \\hline" << std::endl;
	std::cout << "\\hline \\hline" << std::endl;
	std::cout << "\\multicolumn{4}{|l|}{" << suffix[sample] << " final state} \\"<<"\\  " << std::endl;
	std::cout << "   & $N_{jets}=0$       & $N_{jets}=1$      & $N_{jets} \\geq 2$ \\"<<"\\ \\hline" <<std::endl;
      } else {
	std::cout << "\\begin{tabular}{||c|c|c|} \\hline" << std::endl;
	std::cout << "\\hline \\hline" << std::endl;
	std::cout << "\\multicolumn{3}{|l|}{" << suffix[sample] << " final state} \\"<<"\\  " << std::endl;
	std::cout << "   $N_{jets}=0$       & $N_{jets}=1$      & $N_{jets} \\geq 2$ \\"<<"\\ \\hline" <<std::endl;
      }
    } else {      
      std::cout<<"===================================================="<<std::endl<<std::endl;
      std::cout<<"Table for "<< suffix[sample]<<std::endl;
      std::cout<<" |   *sample* |       *nJet = 0*       |      *nJet = 1*        |      *nJet >= 2*      |"<<std::endl;
    }
    double n0all =0;
    double n0allE = 0;
    double n1all =0;
    double n1allE = 0;
    double n2all =0;
    double n2allE = 0;

    double n0sig =0;
    double n0sigE = 0;
    double n1sig =0;
    double n1sigE = 0;
    double n2sig =0;
    double n2sigE = 0;
    for(int iH=0; iH< nHists; ++iH){
      TH1F* h1F = (TH1F*)(thisStack->GetHists()->At(iH));
      TObjArray* sampleNs = TString(h1F->GetName()).Tokenize("_");
      bool isSig = std::string(sampleNs->At(0)->GetName()) == std::string(signalS);

      double n0 = h1F->GetBinContent(1); 
      double n0E = h1F->GetBinError(1);
      n0all += n0;
      n0allE = sqrt(n0allE*n0allE + n0E*n0E);
      if (isSig){
	n0sig += n0;
	n0sigE = sqrt(n0sigE*n0sigE + n0E*n0E);
      }

      double n1 = h1F->GetBinContent(2);
      double n1E = h1F->GetBinError(2);
      n1all += n1;
      n1allE = sqrt(n1allE*n1allE + n1E*n1E);
      if (isSig){
	n1sig += n1;
	n1sigE = sqrt(n1sigE*n1sigE + n1E*n1E);
      }

      int nBins = h1F->GetNbinsX();
      double n2 = 0; 
      for (int i=3; i<= nBins+1; ++i) n2+= h1F->GetBinContent(i);
      double n2E = 0; 
      for (int i=3; i<= nBins+1; ++i) n2E += h1F->GetBinError(i)*h1F->GetBinError(i);
      n2E = sqrt(n2E);
      n2all += n2;
      n2allE = sqrt(n2allE*n2allE + n2E*n2E);
      if (isSig){
	n2sig += n2;
	n2sigE = sqrt(n2sigE*n2sigE + n2E*n2E);
      }

      bool showFirstCol = true;
      if (latex && (sample==1 || sample == 3)) showFirstCol=false;
      std::string firstCol= Form("%9s ",sampleNs->At(0)->GetName());
      
      std::cout << beginL;
      if (showFirstCol) std::cout<< firstCol << colSep;
      std::cout << mathSep << formatFloat(n0,formatS) <<pmSign<< formatFloat(n0E,formatS)<<mathSep<<colSep
		<< mathSep << formatFloat(n1,formatS) <<pmSign<< formatFloat(n1E,formatS)<<mathSep<<colSep
		<< mathSep << formatFloat(n2,formatS) <<pmSign<< formatFloat(n2E,formatS)<<mathSep
		<< endL
		<< std::endl;
    }// for(int iH=0; iH< nHists
    // now print the total and related stuff
    if (latex) std::cout<<"\\hline"<<std::endl;
    bool showFirstCol = true;
    if (latex && (sample==1 || sample == 3)) showFirstCol=false;
    
    std::string firstCol = Form("%9s ","Total");
    std::cout << beginL;
    if (showFirstCol) std::cout<< firstCol << colSep;
    std::cout << mathSep << formatFloat(n0all,formatS) <<pmSign<< formatFloat(n0allE,formatS)<<mathSep<<colSep
	      << mathSep << formatFloat(n1all,formatS) <<pmSign<< formatFloat(n1allE,formatS)<<mathSep<<colSep
	      << mathSep << formatFloat(n2all,formatS) <<pmSign<< formatFloat(n2allE,formatS)<<mathSep
	      << endL
	      << std::endl;
    firstCol = Form("%9s ","Total B");
    std::cout << beginL;
    if (showFirstCol) std::cout<< firstCol << colSep;
    std::cout << mathSep << formatFloat(n0all-n0sig,formatS) <<pmSign<< formatFloat(sqrt(n0allE*n0allE-n0sigE*n0sigE),formatS)<<mathSep<<colSep
	      << mathSep << formatFloat(n1all-n1sig,formatS) <<pmSign<< formatFloat(sqrt(n1allE*n1allE-n1sigE*n1sigE),formatS)<<mathSep<<colSep
	      << mathSep << formatFloat(n2all-n2sig,formatS) <<pmSign<< formatFloat(sqrt(n2allE*n2allE-n2sigE*n2sigE),formatS)<<mathSep
	      << endL
	      << std::endl;

    firstCol = Form("%9s ","S/sqrt(S+B)");
    std::cout << beginL;
    if (showFirstCol) std::cout<< firstCol << colSep;
    double rat = 0; double ratE = 0;
    
    getSoverRootN(rat, ratE, n0sig, n0all, n0sigE, n0allE);
    std::cout << mathSep << formatFloat(rat,formatS) <<pmSign<< formatFloat(ratE,formatS)<<mathSep<<colSep;
    getSoverRootN(rat, ratE, n1sig, n1all, n1sigE, n1allE);
    std::cout << mathSep << formatFloat(rat,formatS) <<pmSign<< formatFloat(ratE,formatS)<<mathSep<<colSep;
    getSoverRootN(rat, ratE, n2sig, n2all, n2sigE, n2allE);
    std::cout << mathSep << formatFloat(rat,formatS) <<pmSign<< formatFloat(ratE,formatS)<<mathSep
	      << endL
	      << std::endl;

    firstCol = Form("%9s "," prob 30\\%, s");
    std::cout << beginL;
    if (showFirstCol) std::cout<< firstCol << colSep;
    double prob = 0; double nSigma = 0;
    
    prob = bOnlyProb(n0sig, n0all-n0sig, std::max(0.3*(n0all-n0sig),sqrt(n0allE*n0allE-n0sigE*n0sigE))); //FIXME
    nSigma = TMath::NormQuantile(min(0.5*(1.+prob),1-1e-15));
    std::cout << mathSep << formatFloat((1.-prob),formatS) <<" , "<< formatFloat(nSigma,formatS)<<mathSep<<colSep;

    prob = bOnlyProb(n1sig, n1all-n1sig, std::max(0.3*(n1all-n1sig),sqrt(n1allE*n1allE-n1sigE*n1sigE))); //FIXME
    nSigma = TMath::NormQuantile(min(0.5*(1.+prob),1-1e-15));
    std::cout << mathSep << formatFloat((1.-prob),formatS) <<" , "<< formatFloat(nSigma,formatS)<<mathSep<<colSep;

    prob = bOnlyProb(n2sig, n2all-n2sig, std::max(0.3*(n2all-n2sig),sqrt(n2allE*n2allE-n2sigE*n2sigE))); //FIXME
    nSigma = TMath::NormQuantile(min(0.5*(1.+prob),1-1e-15));
    std::cout << mathSep << formatFloat((1.-prob),formatS) <<" , "<< formatFloat(nSigma,formatS)<<mathSep
	      << endL
	      << std::endl;

    firstCol = Form("%9s "," prob 50\\%, s");
    std::cout << beginL;
    if (showFirstCol) std::cout<< firstCol << colSep;
    prob = 0; nSigma = 0;
    
    prob = bOnlyProb(n0sig, n0all-n0sig, std::max(0.5*(n0all-n0sig),sqrt(n0allE*n0allE-n0sigE*n0sigE))); //FIXME
    nSigma = TMath::NormQuantile(min(0.5*(1.+prob),1-1e-15));
    std::cout << mathSep << formatFloat((1.-prob),formatS) <<" , "<< formatFloat(nSigma,formatS)<<mathSep<<colSep;

    prob = bOnlyProb(n1sig, n1all-n1sig, std::max(0.5*(n1all-n1sig),sqrt(n1allE*n1allE-n1sigE*n1sigE))); //FIXME
    nSigma = TMath::NormQuantile(min(0.5*(1.+prob),1-1e-15));
    std::cout << mathSep << formatFloat((1.-prob),formatS) <<" , "<< formatFloat(nSigma,formatS)<<mathSep<<colSep;

    prob = bOnlyProb(n2sig, n2all-n2sig, std::max(0.5*(n2all-n2sig),sqrt(n2allE*n2allE-n2sigE*n2sigE))); //FIXME
    nSigma = TMath::NormQuantile(min(0.5*(1.+prob),1-1e-15));
    std::cout << mathSep << formatFloat((1.-prob),formatS) <<" , "<< formatFloat(nSigma,formatS)<<mathSep
	      << endL
	      << std::endl;
    
    if (latex){
      std::cout << "\\hline " << std::endl;
      std::cout << "\\end{tabular}" << std::endl;
      std::cout << "\\end{minipage}" << std::endl;
    }else{
      std::cout << " " << std::endl;
    }
  }//over samples
 
  if (latex) {
    std::cout << "}" <<std::endl; //close \small fontsize
    std::cout << "\\caption{\\label{tab:dummyLabel} Put your caption here}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
  }
}

void printN2JetsColumns( const std::vector<std::string>& pfxs, int nPfx, bool latex = false, bool noErrs = false){
  char* suffix[4];
  suffix[0] = "ee";
  suffix[1] = "mm";
  suffix[2] = "em";
  suffix[3] = "all";


  if (latex) {
    std::cout << "\\begin{table}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "{\\small" << std::endl;
  }

     
  for (int sample=0; sample<4; sample++) {
    
    THStack* thisStack[10]; memset(thisStack, 0, sizeof(thisStack));
    for (int iP = 0; iP< nPfx; ++iP){
      hist::stack(Form("st_%s_hnJet_%s", pfxs[iP].c_str(), suffix[sample]), 
		  Form("%s_[a-zA-Z0-9]*_hnJet_%s$", pfxs[iP].c_str(), suffix[sample]));
      thisStack[iP] = (THStack*) 
	gROOT->FindObjectAny(Form("st_%s_hnJet_%s", pfxs[iP].c_str(), suffix[sample]));
    }
    //    std::cout<<"Found "<<thisStack->GetName()<<std::endl;

    int nHists = thisStack[0]->GetHists()->GetSize();
    if (latex) {
      std::cout << "\\begin{minipage}{0.48\\textwidth}" << std::endl;
      std::cout << "\\begin{tabular}{|l|";
      for(int iP=0;iP<nPfx; ++iP) std::cout<<"c|";
      std::cout<<"} \\hline" << std::endl;
      std::cout << "\\hline \\hline" << std::endl;
      std::cout << "\\multicolumn{"<< nPfx+1 <<"}{|l|}{" << suffix[sample] << " final state} \\"<<"\\  " << std::endl;
      std::cout<<" & ";
      for (int iP = 0; iP < nPfx; ++iP) {
	std::cout<< pfxs[iP].c_str();
	if (iP != nPfx - 1) std::cout<<" & ";
      }
      std::cout <<"\\\\ \\hline" <<std::endl;
    } else {
      std::cout<<"===================================================="<<std::endl;
      std::cout<<suffix[sample]<<std::endl;
      std::cout<<" |   sample  |";
      for (int iP = 0; iP < nPfx; ++iP){
	std::cout<< "        "<<pfxs[iP].c_str()<<"        |";
      }
      std::cout<<std::endl;
    }
    double n2all[10]; memset(n2all, 0, sizeof(n2all));
    double n2allE[10]; memset(n2allE, 0, sizeof(n2allE));
    for(int iH=0; iH< nHists; ++iH){
      TH1F* h1F[10];
      for (int iP = 0; iP < nPfx; ++iP){
	h1F[iP] = (TH1F*)(thisStack[iP]->GetHists()->At(iH));
      }
      TObjArray* sampleNs = TString(h1F[0]->GetName()).Tokenize("_");

      //fisrt column first
      if (latex) {
	std::cout<<Form("%9s",sampleNs->At(1)->GetName()) << " & ";
      } else {
	std::cout<<" | "<<Form("%9s",sampleNs->At(1)->GetName())<<" | ";
      }

      int nBins = h1F[0]->GetNbinsX();
      for (int iP = 0; iP < nPfx; ++iP){
	double n2 = 0; 
	for (int i=3; i<= nBins+1; ++i) n2+= h1F[iP]->GetBinContent(i);
	double n2E = 0; 
	for (int i=3; i<= nBins+1; ++i) n2E += h1F[iP]->GetBinError(i)*h1F[iP]->GetBinError(i);
	n2E = sqrt(n2E);
	n2all[iP] += n2;
	n2allE[iP] += n2E*n2E;
	if (latex) {
	  std::cout <<" $"<<Form("%6.1f",n2);
	  if (! noErrs){
	    std::cout <<" \\pm "<< Form("%6.1f",n2E);
	  }
	  std::cout <<" $";
	  if (iP != nPfx - 1){
	    std::cout<< " & ";
	  } else {
	    std::cout <<" \\"<<"\\"<<std::endl;
	  }
	} else {
	  std::cout <<Form("%6.1f",n2) <<" &plusmn; "<<Form("%6.1f",n2E) <<" | ";
	  if (iP == nPfx - 1) std::cout<<std::endl;
	}
      }
    }
    if (latex) {
      std::cout<<"\\hline"<<std::endl;
      std::cout<<Form("%9s","Total")<<" & ";
      for (int iP = 0; iP < nPfx; ++iP){
	std::cout<<" $"<<Form("%6.1f",n2all[iP]);
	if (! noErrs) {
	  std::cout<<" \\pm "<< Form("%6.1f",n2allE[iP]);
	}
	std::cout<<" $";
	if (iP != nPfx - 1){
	  std::cout<< " & ";
	} else {
	  std::cout <<" \\"<<"\\"<<std::endl;
	}
      }
      std::cout << "\\hline " << std::endl;
      std::cout << "\\end{tabular}" << std::endl;
      std::cout << "\\end{minipage}" << std::endl;
    } else {
      std::cout<<" | "<<Form("%9s","Total")<<" | ";
      for (int iP = 0; iP < nPfx; ++iP){
	std::cout<<Form("%6.1f",n2all[iP]) <<" &plusmn; "<<Form("%6.1f",n2allE[iP]) <<" | ";
	if (iP == nPfx - 1) std::cout<<std::endl;
      }
    }
    if(!latex)
      std::cout << " " << std::endl;
  }
 
  if (latex) {
    std::cout << "}" <<std::endl; //close \small fontsize
    std::cout << "\\caption{\\label{tab:dummyLabel} Put your caption here}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
  }
}


void getJESSyst(const char* formatS = "% 4.f"){
  //input the names of the files here
  hist::loadHist("sk_tdilgp032309/myHist_1957888__OS_noDupWt_isoDil08_preDil08noIso_preMet08_outZ08_hltMu9E15.root", 
 		 "mid", "*_hnJet_*");
  hist::loadHist("sk_tdilgp032309/myHist_10346496__OS_noDupWt_isoDil08_preDil08noIso_preMet08_outZ08_hltMu9E15_jets10Up.root", 
 		 "jup", "*_hnJet_*");
  hist::loadHist("sk_tdilgp032309/myHist_18735104__OS_noDupWt_isoDil08_preDil08noIso_preMet08_outZ08_hltMu9E15_jets10Dn.root", 
 		 "jdn", "*_hnJet_*");

  char* suffix[4];
  suffix[0] = "ee";
  suffix[1] = "mm";
  suffix[2] = "em";
  suffix[3] = "all";
  char* suffixLatex[4] = { "$\\eepm$", "$\\mmpm$", "$\\empm$", "All"};

  TH1F* ttdil_mid[4];
  TH1F* ttdil_jup[4];
  TH1F* ttdil_jdn[4];

  TH1F* tw_mid[4];
  TH1F* tw_jup[4];
  TH1F* tw_jdn[4];

  TH1F* DYtautau_mid[4];
  TH1F* DYtautau_jup[4];
  TH1F* DYtautau_jdn[4];

  TH1F* ww_mid[4];
  TH1F* ww_jup[4];
  TH1F* ww_jdn[4];

  TH1F* wz_mid[4];
  TH1F* wz_jup[4];
  TH1F* wz_jdn[4];

  TH1F* zz_mid[4];
  TH1F* zz_jup[4];
  TH1F* zz_jdn[4];

  for (int iSam = 0; iSam < 4; ++iSam){
    ttdil_mid[iSam] = (TH1F*)gDirectory->Get(Form("mid_ttdil_hnJet_%s",suffix[iSam]));
    ttdil_jup[iSam] = (TH1F*)gDirectory->Get(Form("jup_ttdil_hnJet_%s",suffix[iSam]));
    ttdil_jdn[iSam] = (TH1F*)gDirectory->Get(Form("jdn_ttdil_hnJet_%s",suffix[iSam]));

    tw_mid[iSam] = (TH1F*)gDirectory->Get(Form("mid_tW_hnJet_%s",suffix[iSam]));
    tw_jup[iSam] = (TH1F*)gDirectory->Get(Form("jup_tW_hnJet_%s",suffix[iSam]));
    tw_jdn[iSam] = (TH1F*)gDirectory->Get(Form("jdn_tW_hnJet_%s",suffix[iSam]));

    DYtautau_mid[iSam] = (TH1F*)gDirectory->Get(Form("mid_DYtautau_hnJet_%s",suffix[iSam]));
    DYtautau_jup[iSam] = (TH1F*)gDirectory->Get(Form("jup_DYtautau_hnJet_%s",suffix[iSam]));
    DYtautau_jdn[iSam] = (TH1F*)gDirectory->Get(Form("jdn_DYtautau_hnJet_%s",suffix[iSam]));

    ww_mid[iSam] = (TH1F*)gDirectory->Get(Form("mid_ww_hnJet_%s",suffix[iSam]));
    ww_jup[iSam] = (TH1F*)gDirectory->Get(Form("jup_ww_hnJet_%s",suffix[iSam]));
    ww_jdn[iSam] = (TH1F*)gDirectory->Get(Form("jdn_ww_hnJet_%s",suffix[iSam]));

    wz_mid[iSam] = (TH1F*)gDirectory->Get(Form("mid_wz_hnJet_%s",suffix[iSam]));
    wz_jup[iSam] = (TH1F*)gDirectory->Get(Form("jup_wz_hnJet_%s",suffix[iSam]));
    wz_jdn[iSam] = (TH1F*)gDirectory->Get(Form("jdn_wz_hnJet_%s",suffix[iSam]));

    zz_mid[iSam] = (TH1F*)gDirectory->Get(Form("mid_zz_hnJet_%s",suffix[iSam]));
    zz_jup[iSam] = (TH1F*)gDirectory->Get(Form("jup_zz_hnJet_%s",suffix[iSam]));
    zz_jdn[iSam] = (TH1F*)gDirectory->Get(Form("jdn_zz_hnJet_%s",suffix[iSam]));

  }

  //  if(latex)
  // just latex for now
  std::cout<<  "\\begin{table}"                                              <<std::endl;
  std::cout<<  "\\begin{center}"                                             <<std::endl;
  std::cout<<  "\\begin{tabular}{|l|cc|cc|cc|cc|cc|cc|}"                              <<std::endl;
  std::cout<<  "\\hline "                                                    <<std::endl;
  std::cout<<  " & \\multicolumn{4}{|c|}{\\ttbar} & ";
  std::cout<<  "   \\multicolumn{4}{|c|}{single-top} & ";
  std::cout<<  "   \\multicolumn{4}{|c|}{$VV+\\dytt$} \\\\"<<std::endl;

  std::cout<<  " & \\multicolumn{2}{|c|}{$N_\\jets=0,1$} & ";
  std::cout<<  "   \\multicolumn{2}{|c|}{$N_\\jets\\geq 2$} & "<<std::endl;
  std::cout<<  "   \\multicolumn{2}{|c|}{$N_\\jets=0,1$} & ";
  std::cout<<  "   \\multicolumn{2}{|c|}{$N_\\jets\\geq 2$} & "<<std::endl;
  std::cout<<  " \t  \\multicolumn{2}{|c|}{$N_\\jets=0,1$} &";
  std::cout<<  "   \\multicolumn{2}{|c|}{$N_\\jets\\geq 2$} \\\\"<<std::endl;
  std::cout<<  "Channel & $-10\\%$ & $+10\\%$ & $-10\\%$ & $+10\\%$ & ";
  std::cout<<  "  $-10\\%$ & $+10\\%$ & $-10\\%$ & $+10\\%$ & " <<std::endl;
  std::cout<<  "  $-10\\%$ & $+10\\%$ & $-10\\%$ & $+10\\%$ \\\\ " <<std::endl;
  std::cout<<  "\\hline "                                                    <<std::endl;
  for (int iSam=0; iSam< 4; ++iSam){
    if (iSam == 3) std::cout<<"\\hline"<<std::endl;
    std::cout<<suffixLatex[iSam]<<" & ";
    std::cout<< Form("$% 4.2f$ & ", (ttdil_jdn[iSam]->Integral(0,2)/ttdil_mid[iSam]->Integral(0,2) - 1.)*100.);
    std::cout<< Form("$% 4.2f$ & ", (ttdil_jup[iSam]->Integral(0,2)/ttdil_mid[iSam]->Integral(0,2) - 1.)*100.);
    std::cout<< Form("$% 4.2f$ & ", (ttdil_jdn[iSam]->Integral(3,9)/ttdil_mid[iSam]->Integral(3,9) - 1.)*100.);
    std::cout<< Form("$% 4.2f$ & ", (ttdil_jup[iSam]->Integral(3,9)/ttdil_mid[iSam]->Integral(3,9) - 1.)*100.);

    std::cout<< Form("$% 4.2f$ & ", (tw_jdn[iSam]->Integral(0,2)/tw_mid[iSam]->Integral(0,2) - 1.)*100.);
    std::cout<< Form("$% 4.2f$ & ", (tw_jup[iSam]->Integral(0,2)/tw_mid[iSam]->Integral(0,2) - 1.)*100.);
    std::cout<< Form("$% 4.2f$ & ", (tw_jdn[iSam]->Integral(3,9)/tw_mid[iSam]->Integral(3,9) - 1.)*100.);
    std::cout<< Form("$% 4.2f$ & ", (tw_jup[iSam]->Integral(3,9)/tw_mid[iSam]->Integral(3,9) - 1.)*100.);

    //0,1
    double dyttVV_mid = DYtautau_mid[iSam]->Integral(0,2);
    dyttVV_mid += ww_mid[iSam]->Integral(0,2);
    dyttVV_mid += wz_mid[iSam]->Integral(0,2);
    dyttVV_mid += zz_mid[iSam]->Integral(0,2);
    double dyttVV_jup = DYtautau_jup[iSam]->Integral(0,2);
    dyttVV_jup += ww_jup[iSam]->Integral(0,2);
    dyttVV_jup += wz_jup[iSam]->Integral(0,2);
    dyttVV_jup += zz_jup[iSam]->Integral(0,2);
    double dyttVV_jdn = DYtautau_jdn[iSam]->Integral(0,2);
    dyttVV_jdn += ww_jdn[iSam]->Integral(0,2);
    dyttVV_jdn += wz_jdn[iSam]->Integral(0,2);
    dyttVV_jdn += zz_jdn[iSam]->Integral(0,2);
    std::cout<< Form("$% 4.2f$ & ", (dyttVV_jdn/dyttVV_mid - 1.)*100.);
    std::cout<< Form("$% 4.2f$ & ", (dyttVV_jup/dyttVV_mid - 1.)*100.);
    //gt2
    dyttVV_mid = DYtautau_mid[iSam]->Integral(3,9);
    dyttVV_mid += ww_mid[iSam]->Integral(3,9);
    dyttVV_mid += wz_mid[iSam]->Integral(3,9);
    dyttVV_mid += zz_mid[iSam]->Integral(3,9);
    dyttVV_jup = DYtautau_jup[iSam]->Integral(3,9);
    dyttVV_jup += ww_jup[iSam]->Integral(3,9);
    dyttVV_jup += wz_jup[iSam]->Integral(3,9);
    dyttVV_jup += zz_jup[iSam]->Integral(3,9);
    dyttVV_jdn = DYtautau_jdn[iSam]->Integral(3,9);
    dyttVV_jdn += ww_jdn[iSam]->Integral(3,9);
    dyttVV_jdn += wz_jdn[iSam]->Integral(3,9);
    dyttVV_jdn += zz_jdn[iSam]->Integral(3,9);
    std::cout<< Form("$% 4.2f$ & ", (dyttVV_jdn/dyttVV_mid - 1.)*100.);
    std::cout<< Form("$% 4.2f$ \\\\ ", (dyttVV_jup/dyttVV_mid - 1.)*100.);

    std::cout<<std::endl;
  }
  std::cout<<  "\\hline "                                                    <<std::endl;
  std::cout<<  "\\end{tabular}"                                              <<std::endl;
  std::cout<<  "\\caption{\\label{tab:dummyLabel} Relative changes in the number of selected events"<<std::endl;
  std::cout<<  "expressed in $\\%$ for changes in the jet energy scale down and up by $10\\%$.}"    <<std::endl;
  std::cout<<  "\\end{center}"                                               <<std::endl;
  std::cout<<  "\\end{table}"                                                <<std::endl;
}

float histSum(TH1F* h, int minB, int maxB){
  float res = 0;
  for (int iB = minB; iB <= maxB; ++iB){
    res += h->GetBinContent(iB);
  }
  return res;
}
float histSumErr2(TH1F* h, int minB, int maxB){
  float res2 = 0;
  for (int iB = minB; iB <= maxB; ++iB){
    res2 += h->GetBinError(iB)*h->GetBinError(iB);
  }
  return res2;
}
float histSumErr(TH1F* h, int minB, int maxB){
  return sqrt(histSumErr2( h, minB, maxB));
}

void printFlow(std::string fList, const char* formatS = "%6.1f", bool latex = false){
  FILE *confF = fopen(fList.c_str(), "r");
  char inFile[1024]; memset(inFile, 0, sizeof(inFile));
  char conf[20]; memset(conf, 0, sizeof(conf));
  long long bitMask=0;
  
  int nCount = 0;
  std::vector<std::string> inFileSV(20);
  std::vector<std::string> confSV(20);
  std::vector<TFile*> inFileV(20,(TFile*)0);
  while (fscanf(confF, "%s %s", inFile, conf) == 2) {
    //    std::cout<<inFile<<" is "<<conf<<std::endl;
    if (nCount > 19) return;
    inFileSV[nCount] = std::string(inFile);
    confSV[nCount]  = std::string(conf);
    inFileV[nCount] = new TFile(inFileSV[nCount].c_str());
    if (inFileV[nCount] == 0){
      std::cout<<"Failed to read "<<inFile<<std::endl;
      return;
    }
    std::cout<<"Got "<<inFileV[nCount]->GetName()<<" "<<confSV[nCount].c_str()<<std::endl;
    nCount++;
  }

  std::vector<TH1F*> h_topdil_mm(20,(TH1F*)0);
  std::vector<TH1F*> h_topdil_ee(20,(TH1F*)0);
  std::vector<TH1F*> h_topdil_em(20,(TH1F*)0);
  std::vector<TH1F*> h_topdil_all(20,(TH1F*)0);
  for ( int iF = 0; iF < nCount; ++iF){
    h_topdil_mm[iF] = (TH1F*)inFileV[iF]->Get("ttdil_hnJet_mm");
    h_topdil_ee[iF] = (TH1F*)inFileV[iF]->Get("ttdil_hnJet_ee");
    h_topdil_em[iF] = (TH1F*)inFileV[iF]->Get("ttdil_hnJet_em");
    h_topdil_all[iF] = (TH1F*)inFileV[iF]->Get("ttdil_hnJet_all");
  }

  std::string pmSign  = latex ? " \\pm " : " &plusmn; ";
  std::string colSep  = latex ? " & " : " | ";
  std::string hcolSep  = latex ? " & " : "* | *";
  std::string beginL  = latex ? ""   : " | ";
  std::string endL    = latex ? " \\\\ " : " | ";
  std::string mathSep = latex ? "$" : "";
  

  if (latex) {
    std::cout << "\\begin{table}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "{\\small" << std::endl;
  }

  std::cout<<beginL<<"*Mode* ";
  for (int iF = 0; iF< nCount; ++iF) std::cout<<hcolSep<<confSV[iF].c_str();
  std::cout<<endL<<std::endl;

  std::cout<<beginL<<"mm ";
  for (int iF = 0; iF< nCount; ++iF) std::cout<<colSep<<formatFloat(histSum(h_topdil_mm[iF],0,100),formatS);
  std::cout<<endL<<std::endl;
  std::cout<<beginL<<"ee ";
  for (int iF = 0; iF< nCount; ++iF) std::cout<<colSep<<formatFloat(histSum(h_topdil_ee[iF],0,100),formatS);
  std::cout<<endL<<std::endl;
  std::cout<<beginL<<"em ";
  for (int iF = 0; iF< nCount; ++iF) std::cout<<colSep<<formatFloat(histSum(h_topdil_em[iF],0,100),formatS);
  std::cout<<endL<<std::endl;
  std::cout<<beginL<<"all ";
  for (int iF = 0; iF< nCount; ++iF) std::cout<<colSep<<formatFloat(histSum(h_topdil_all[iF],0,100),formatS);
  std::cout<<endL<<std::endl;

}

// read the $hPrefixes_hnJet_$hSuffix histograms from the vector of files inFiles and print out the values
// printStyle: 0 -- value only; 1 -- value +/- error
void printFlowOneLine(const std::vector<TFile*>& inFiles, float minX, float maxX,
		      const std::vector<std::string>& hPrefixes, 
		      const char* modeName, const char* hSuffix, 
		      float scaleF = 1.0,
		      int nSignificant = 2, unsigned int  printStyle = 1,
		      bool latex = false, bool printDebug = false){

  int nCount = inFiles.size();
  int nPfxs = hPrefixes.size();

  TH1F* h_tmp = 0;
  std::vector<double> inValue(nCount, 0);
  std::vector<double> inError(nCount, 0);

  for ( int iF = 0; iF < nCount; ++iF){
    double tVal = 0; double tErr2 = 0;
    for (int iPfx = 0; iPfx  < nPfxs; ++iPfx){
      h_tmp = (TH1F*)inFiles[iF]->Get(Form("%s_hnJet_%s", hPrefixes[iPfx].c_str(), hSuffix));
      if (printDebug) std::cout<<" Got  hName "<<h_tmp->GetName()<<" : Int = "<<h_tmp->Integral(-1,1000)<<std::endl;
      int iMin = h_tmp->FindBin(minX);
      int iMax = h_tmp->FindBin(maxX);
      tVal += histSum(h_tmp, iMin, iMax);
      tErr2 += histSumErr2(h_tmp, iMin, iMax);
    }

    inValue[iF] =tVal*scaleF;
    inError[iF] =sqrt(tErr2)*scaleF;
  }

  std::string pmSign  = latex ? " \\pm " : " &plusmn; ";
  std::string colSep  = latex ? " & " : " | ";
  std::string hcolSep  = latex ? " & " : "* | *";
  std::string beginL  = latex ? ""   : " | ";
  std::string endL    = latex ? " \\\\ " : " | ";
  std::string mathSep = latex ? "$" : "";
  

  std::cout<<beginL<<modeName<<" ";
  for (int iF = 0; iF< nCount; ++iF){
    std::cout<<"\t"<<colSep<<formatFloat(inValue[iF], inError[iF], nSignificant, printStyle, latex);
  }
  std::cout<<endL<<std::endl;
}

void printFlowMCnData(std::string fList, float scaleF = 1., bool latex = false, bool debugPrint = false){
  std::cout<<fList<<std::endl;
  FILE *confF = fopen(fList.c_str(), "r");
  char inFile[1024]; memset(inFile, 0, sizeof(inFile));
  char inFilePP[1024]; memset(inFilePP, 0, sizeof(inFilePP));
  char conf[20]; memset(conf, 0, sizeof(conf));
  long long bitMask=0;
  
  int nCount = 0;
  while (fscanf(confF, "%s %s %s", inFile, inFilePP, conf) == 3) nCount++;
  std::cout<<"Read "<<nCount<<" files"<<std::endl;

  fclose(confF);

  confF = fopen(fList.c_str(), "r");
  std::vector<std::string> inFileSV(nCount);
  std::vector<std::string> inFilePPSV(nCount);
  std::vector<std::string> confSV(nCount);
  std::vector<TFile*> inFileV(nCount, (TFile*)0);
  std::vector<TFile*> inFilePPV(nCount, (TFile*)0);
  nCount = 0;
  while (fscanf(confF, "%s %s %s", inFilePP, inFile, conf) == 3) {
    inFileSV[nCount] = std::string(inFile);
    inFilePPSV[nCount] = std::string(inFilePP);
    confSV[nCount] = std::string(conf);
    inFileV[nCount] = new TFile(inFileSV[nCount].c_str());
    inFilePPV[nCount] = new TFile(inFilePPSV[nCount].c_str());
    if (inFileV[nCount] == 0 || inFilePPV[nCount] == 0){
      std::cout<<"Failed to read "<<inFile<<std::endl;
      return;
    }
    if (debugPrint) std::cout<<"Got "<<inFileV[nCount]->GetName()<<" "<<confSV[nCount].c_str()<<std::endl;
    nCount++;
    memset(inFile, 0, sizeof(inFile));
    memset(conf, 0, sizeof(conf));
  }

  std::string pmSign  = latex ? " \\pm " : " &plusmn; ";
  std::string colSep  = latex ? " & " : " | ";
  std::string hcolSep  = latex ? " & " : "* | *";
  std::string beginL  = latex ? ""   : " | ";
  std::string endL    = latex ? " \\\\ " : " | ";
  std::string mathSep = latex ? "$" : "";
  
  std::vector<std::string> modeSV;
  modeSV.push_back("ee");
  modeSV.push_back("mm");
  modeSV.push_back("em");
  modeSV.push_back("all");

  
  for (int iMode = 0; iMode < modeSV.size(); ++iMode){
    std::string modeS  = modeSV[iMode];
    if (latex) {
      std::cout << "\\begin{table}" << std::endl;
      std::cout << "\\begin{center}" << std::endl;
      std::cout << "{\\small \\caption{Mode "<<modeS<<" }" << std::endl;
      std::cout << "\\begin{tabular}{l";
      for (int i=0; i< nCount; ++i) std::cout<<"c";
      std::cout << "}\\hline"<<std::endl;
    }
    
    if (latex) std::cout<<beginL<<" Channel ";
    else std::cout<<beginL<<" *Channel* ";
    for (int iF = 0; iF< nCount; ++iF) std::cout<<hcolSep<<confSV[iF].c_str();
    std::cout<<endL<<std::endl;
    if (latex) std::cout<<"\\hline"<<std::endl;
    std::vector<std::string> hPrefixes;
    
    hPrefixes.push_back(std::string("ttdil"));
    printFlowOneLine(inFileV, -1, 100, hPrefixes, "Dilepton $t\\bar{t}$", modeS.c_str(), scaleF, 1, 1, true); if(latex) std::cout <<"\\hline"<<std::endl;
    hPrefixes.clear(); hPrefixes.push_back(std::string("VV"));
    printFlowOneLine(inFileV, -1, 100, hPrefixes, "$VV$", modeS.c_str(), scaleF, 1, 1, true);
    hPrefixes.clear(); hPrefixes.push_back(std::string("tw"));
    printFlowOneLine(inFileV, -1, 100, hPrefixes, "Single-top $tW$", modeS.c_str(), scaleF, 1, 1, true);
    hPrefixes.clear(); hPrefixes.push_back(std::string("DYtautau"));
    printFlowOneLine(inFileV, -1, 100, hPrefixes, "Drell-Yan $\\tau\\tau$", modeS.c_str(), scaleF, 1, 1, true);
    hPrefixes.clear(); hPrefixes.push_back(std::string("DYee")); hPrefixes.push_back(std::string("DYmm")); 
    printFlowOneLine(inFileV, -1, 100, hPrefixes, "Drell-Yan $ee/\\mu\\mu$", modeS.c_str(), scaleF, 1, 1, true);
    hPrefixes.clear(); hPrefixes.push_back(std::string("ttotr"));
    printFlowOneLine(inFileV, -1, 100, hPrefixes, "Non-dilepton $t\\bar{t}$", modeS.c_str(), scaleF, 1, 1, true);
    hPrefixes.clear(); hPrefixes.push_back(std::string("wjets"));
    printFlowOneLine(inFileV, -1, 100, hPrefixes, "W+jets", modeS.c_str(), scaleF, 1, 1, true);
    hPrefixes.clear(); hPrefixes.push_back(std::string("qcd15"));  hPrefixes.push_back(std::string("qcd30"));
    printFlowOneLine(inFileV, -1, 100, hPrefixes, "QCD multijets", modeS.c_str(), scaleF, 1, 1, true);  if(latex) std::cout <<"\\hline"<<std::endl;
    hPrefixes.clear(); 
    hPrefixes.push_back(std::string("ttdil")); 
    hPrefixes.push_back(std::string("VV"));
    hPrefixes.push_back(std::string("tw")); 
    hPrefixes.push_back(std::string("DYtautau")); 
    hPrefixes.push_back(std::string("DYee")); hPrefixes.push_back(std::string("DYmm"));
    hPrefixes.push_back(std::string("ttotr"));
    hPrefixes.push_back(std::string("wjets")); 
    hPrefixes.push_back(std::string("qcd15"));  hPrefixes.push_back(std::string("qcd30"));
    printFlowOneLine(inFileV, -1, 100, hPrefixes, "Total simulated", modeS.c_str(), scaleF, 1, 1, true);  if(latex) std::cout <<"\\hline"<<std::endl;
    hPrefixes.clear(); hPrefixes.push_back(std::string("data"));
    printFlowOneLine(inFilePPV, -1, 100, hPrefixes, "Data", modeS.c_str(), 1., 1, 0, true);  if(latex) std::cout <<"\\hline"<<std::endl;
    
    if (latex) {
      std::cout << "\\end{tabular}"<<std::endl;
      std::cout << "}"<<std::endl;
      std::cout << "\\end{center}" << std::endl;
      std::cout << "\\end{table}" << std::endl;
      std::cout << std::endl;
    }
  } // loop over modes

}
