void printNJets( bool latex=false, const char* formatS = "%6.1f"){
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
    
    hist::stack(Form("st_hnJet_%s", suffix[sample]), Form("hnJet_%s$", suffix[sample]));
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
      std::cout<<"===================================================="<<std::endl;
      std::cout<<suffix[sample]<<std::endl;
      std::cout<<" |   sample  |        nJet = 0        |       nJet = 1         |       nJet >= 2       |"<<std::endl;
    }
    double n0all =0;
    double n0allE = 0;
    double n1all =0;
    double n1allE = 0;
    double n2all =0;
    double n2allE = 0;
    for(int iH=0; iH< nHists; ++iH){
      TH1F* h1F = (TH1F*)(thisStack->GetHists()->At(iH));
      TObjArray* sampleNs = TString(h1F->GetName()).Tokenize("_");
      double n0 = h1F->GetBinContent(1); 
      double n0E = h1F->GetBinError(1);
      n0all += n0;
      n0allE += n0E*n0E;

      double n1 = h1F->GetBinContent(2);
      double n1E = h1F->GetBinError(2);
      n1all += n1;
      n1allE += n1E*n1E;

      int nBins = h1F->GetNbinsX();
      double n2 = 0; 
      for (int i=3; i<= nBins+1; ++i) n2+= h1F->GetBinContent(i);
      double n2E = 0; 
      for (int i=3; i<= nBins+1; ++i) n2E += h1F->GetBinError(i)*h1F->GetBinError(i);
      n2E = sqrt(n2E);
      n2all += n2;
      n2allE += n2E*n2E;

      if (latex) {
	if (sample == 0 || sample ==2){
	  std::cout<<Form("%9s & ",sampleNs->At(0)->GetName()) ;
	}
	std::cout<< "  $"<<Form(Form("%s", formatS),n0) <<" \\pm "<< Form(Form("%s", formatS),n0E) 
		 <<"$ & $"<<Form(Form("%s", formatS),n1) <<" \\pm "<< Form(Form("%s", formatS),n1E) 
		 <<"$ & $"<<Form(Form("%s", formatS),n2) <<" \\pm "<< Form(Form("%s", formatS),n2E) 
		 <<"$ \\"<<"\\"  
		 <<std::endl;
      } else {
	std::cout<<" | "<<Form("%9s",sampleNs->At(0)->GetName())
	         <<" | "<<Form(Form("%s", formatS),n0) <<" &plusmn; "<<Form(Form("%s", formatS),n0E)
	         <<" | "<<Form(Form("%s", formatS),n1) <<" &plusmn; "<<Form(Form("%s", formatS),n1E)
	         <<" | "<<Form(Form("%s", formatS),n2) <<" &plusmn; "<<Form(Form("%s", formatS),n2E)
	         <<" | "<<std::endl;
      }
    }
    if (latex) {
      std::cout<<"\\hline"<<std::endl;
      if (sample == 0 || sample == 2){
	std::cout<<Form("%9s & ","Total") ;
      }
      std::cout<< " $"<<Form(Form("%s", formatS),n0all) <<" \\pm "<< Form(Form("%s", formatS),n0allE) 
	       <<"$ & $"<<Form(Form("%s", formatS),n1all) <<" \\pm "<< Form(Form("%s", formatS),n1allE) 
	       <<"$ & $"<<Form(Form("%s", formatS),n2all) <<" \\pm "<< Form(Form("%s", formatS),n2allE) 
	       <<"$ \\"<<"\\"  
	       <<std::endl;
      std::cout << "\\hline " << std::endl;
      std::cout << "\\end{tabular}" << std::endl;
      std::cout << "\\end{minipage}" << std::endl;
    } else {
      std::cout<<" | "<<Form("%9s","Total")
	       <<" | "<<Form(Form("%s", formatS),n0all) <<" &plusmn; "<<Form(Form("%s", formatS),n0allE)
	       <<" | "<<Form(Form("%s", formatS),n1all) <<" &plusmn; "<<Form(Form("%s", formatS),n1allE)
	       <<" | "<<Form(Form("%s", formatS),n2all) <<" &plusmn; "<<Form(Form("%s", formatS),n2allE)
	       <<" | "<<std::endl;
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
