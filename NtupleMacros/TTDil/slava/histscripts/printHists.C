#include "specFun.C"

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
    
    prob = bOnlyProb(n0sig, n0all-n0sig, max(0.3*(n0all-n0sig),sqrt(n0allE*n0allE-n0sigE*n0sigE))); //FIXME
    nSigma = TMath::NormQuantile(min(0.5*(1.+prob),1-1e-15));
    std::cout << mathSep << formatFloat((1.-prob),formatS) <<" , "<< formatFloat(nSigma,formatS)<<mathSep<<colSep;

    prob = bOnlyProb(n1sig, n1all-n1sig, max(0.3*(n1all-n1sig),sqrt(n1allE*n1allE-n1sigE*n1sigE))); //FIXME
    nSigma = TMath::NormQuantile(min(0.5*(1.+prob),1-1e-15));
    std::cout << mathSep << formatFloat((1.-prob),formatS) <<" , "<< formatFloat(nSigma,formatS)<<mathSep<<colSep;

    prob = bOnlyProb(n2sig, n2all-n2sig, max(0.3*(n2all-n2sig),sqrt(n2allE*n2allE-n2sigE*n2sigE))); //FIXME
    nSigma = TMath::NormQuantile(min(0.5*(1.+prob),1-1e-15));
    std::cout << mathSep << formatFloat((1.-prob),formatS) <<" , "<< formatFloat(nSigma,formatS)<<mathSep
	      << endL
	      << std::endl;

    firstCol = Form("%9s "," prob 50\\%, s");
    std::cout << beginL;
    if (showFirstCol) std::cout<< firstCol << colSep;
    prob = 0; nSigma = 0;
    
    prob = bOnlyProb(n0sig, n0all-n0sig, max(0.5*(n0all-n0sig),sqrt(n0allE*n0allE-n0sigE*n0sigE))); //FIXME
    nSigma = TMath::NormQuantile(min(0.5*(1.+prob),1-1e-15));
    std::cout << mathSep << formatFloat((1.-prob),formatS) <<" , "<< formatFloat(nSigma,formatS)<<mathSep<<colSep;

    prob = bOnlyProb(n1sig, n1all-n1sig, max(0.5*(n1all-n1sig),sqrt(n1allE*n1allE-n1sigE*n1sigE))); //FIXME
    nSigma = TMath::NormQuantile(min(0.5*(1.+prob),1-1e-15));
    std::cout << mathSep << formatFloat((1.-prob),formatS) <<" , "<< formatFloat(nSigma,formatS)<<mathSep<<colSep;

    prob = bOnlyProb(n2sig, n2all-n2sig, max(0.5*(n2all-n2sig),sqrt(n2allE*n2allE-n2sigE*n2sigE))); //FIXME
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
